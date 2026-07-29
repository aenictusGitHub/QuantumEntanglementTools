#!/usr/bin/env julia

using Downloads
using SHA
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_MANIFEST = joinpath(REPOSITORY_ROOT, "UpstreamManifest.toml")

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/check_upstream_updates.jl [options]

Compare UpstreamManifest.toml pins with current upstream refs/archives.

Options:
  --manifest PATH       manifest to check (default: UpstreamManifest.toml)
  --source ID           check only one source ID (repeatable)
  --candidate ID=PATH   compare a candidate checkout/archive with its pin
                        (PATH alone is accepted with exactly one --source)
  --provenance PATH     Julia provenance map (default: PROVENANCE.toml)
  --offline             verify ignored local audit copies without networking
  --no-fail-on-update   report updates but exit successfully
  -h, --help            show this help

Exit status:
  0  checks passed and no update was found (or --no-fail-on-update)
  1  at least one upstream update was found
  2  a manifest, network, command, or integrity check failed
""",
    )
end

function parse_options(args)
    manifest = DEFAULT_MANIFEST
    selected_sources = String[]
    candidate_specs = String[]
    provenance = joinpath(REPOSITORY_ROOT, "PROVENANCE.toml")
    offline = false
    fail_on_update = true
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            usage()
            exit(0)
        elseif arg == "--offline"
            offline = true
        elseif arg == "--no-fail-on-update"
            fail_on_update = false
        elseif startswith(arg, "--manifest=")
            manifest = split(arg, "="; limit=2)[2]
        elseif arg == "--manifest"
            i == length(args) && error("--manifest requires a path")
            i += 1
            manifest = args[i]
        elseif startswith(arg, "--candidate=")
            push!(candidate_specs, split(arg, "="; limit=2)[2])
        elseif arg == "--candidate"
            i == length(args) && error("--candidate requires ID=PATH or PATH")
            i += 1
            push!(candidate_specs, args[i])
        elseif startswith(arg, "--provenance=")
            provenance = split(arg, "="; limit=2)[2]
        elseif arg == "--provenance"
            i == length(args) && error("--provenance requires a path")
            i += 1
            provenance = args[i]
        elseif startswith(arg, "--source=")
            push!(selected_sources, split(arg, "="; limit=2)[2])
        elseif arg == "--source"
            i == length(args) && error("--source requires an ID")
            i += 1
            push!(selected_sources, args[i])
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return (
        manifest=abspath(manifest),
        selected_sources=unique(selected_sources),
        candidate_specs,
        provenance=abspath(provenance),
        offline,
        fail_on_update,
    )
end

mutable struct Report
    ok::Vector{String}
    updates::Vector{String}
    warnings::Vector{String}
    blockers::Vector{String}
    errors::Vector{String}
end

Report() = Report(String[], String[], String[], String[], String[])

function record!(report::Report, kind::Symbol, message::AbstractString)
    destination = getfield(report, kind)
    push!(destination, String(message))
    prefix = Dict(
        :ok => "[OK]",
        :updates => "[UPDATE]",
        :warnings => "[WARN]",
        :blockers => "[BLOCKED]",
        :errors => "[ERROR]",
    )[kind]
    return println(prefix, " ", message)
end

file_sha256(path::AbstractString) = bytes2hex(sha256(read(path)))
bytes_sha256(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
normalize_path(path::AbstractString) = replace(String(path), '\\' => '/')

const CP1252_CONTROLS = Dict{UInt8,Char}(
    0x80 => '€',
    0x82 => '‚',
    0x83 => 'ƒ',
    0x84 => '„',
    0x85 => '…',
    0x86 => '†',
    0x87 => '‡',
    0x88 => 'ˆ',
    0x89 => '‰',
    0x8a => 'Š',
    0x8b => '‹',
    0x8c => 'Œ',
    0x8e => 'Ž',
    0x91 => '‘',
    0x92 => '’',
    0x93 => '“',
    0x94 => '”',
    0x95 => '•',
    0x96 => '–',
    0x97 => '—',
    0x98 => '˜',
    0x99 => '™',
    0x9a => 'š',
    0x9b => '›',
    0x9c => 'œ',
    0x9e => 'ž',
    0x9f => 'Ÿ',
)

function decode_source(bytes::Vector{UInt8})
    isvalid(String, bytes) && return String(bytes)
    io = IOBuffer()
    index = 1
    while index <= length(bytes)
        byte = bytes[index]
        if byte <= 0x7f
            write(io, byte)
            index += 1
            continue
        end
        width = if 0xc2 <= byte <= 0xdf
            2
        elseif 0xe0 <= byte <= 0xef
            3
        elseif 0xf0 <= byte <= 0xf4
            4
        else
            0
        end
        valid_sequence =
            width > 0 &&
            index + width - 1 <= length(bytes) &&
            all(
                continuation -> 0x80 <= continuation <= 0xbf,
                bytes[(index + 1):(index + width - 1)],
            )
        if valid_sequence
            write(io, @view bytes[index:(index + width - 1)])
            index += width
        else
            print(io, get(CP1252_CONTROLS, byte, Char(byte)))
            index += 1
        end
    end
    return String(take!(io))
end

function normalize_newlines(text::AbstractString)
    return replace(replace(String(text), "\r\n" => "\n"), "\r" => "\n")
end

function matlab_signature(bytes::Vector{UInt8})
    lines = split(normalize_newlines(decode_source(bytes)), '\n')
    index = findfirst(line -> occursin(r"^\s*function(?:\s|$)"i, line), lines)
    isnothing(index) && return ""
    header = strip(lines[index])
    while occursin(r"\.\.\.\s*(?:%.*)?$", header) && index < length(lines)
        header = replace(header, r"\.\.\.\s*(?:%.*)?$" => " ")
        index += 1
        header *= " " * strip(lines[index])
    end
    return strip(replace(first(split(header, "%"; limit=2)), r"\s+" => " "))
end

function matlab_help_text(bytes::Vector{UInt8})
    lines = split(normalize_newlines(decode_source(bytes)), '\n')
    help_lines = String[]
    saw_function = false
    continuing_signature = false
    for line in lines
        stripped = strip(line)
        if isempty(stripped)
            isempty(help_lines) || push!(help_lines, "")
            continue
        elseif startswith(stripped, "%")
            push!(help_lines, rstrip(stripped))
            continue
        elseif !saw_function && occursin(r"^function(?:\s|$)"i, stripped)
            saw_function = true
            continuing_signature = endswith(stripped, "...")
            continue
        elseif continuing_signature
            continuing_signature = endswith(stripped, "...")
            continue
        end
        break
    end
    while !isempty(help_lines) && isempty(last(help_lines))
        pop!(help_lines)
    end
    return join(help_lines, "\n")
end

function matlab_file_map(root::AbstractString)
    result = Dict{String,Vector{UInt8}}()
    for (directory, subdirectories, files) in walkdir(root)
        filter!(name -> name != ".git", subdirectories)
        for filename in files
            endswith(lowercase(filename), ".m") || continue
            path = joinpath(directory, filename)
            result[normalize_path(relpath(path, root))] = read(path)
        end
    end
    return result
end

function license_file_map(root::AbstractString)
    result = Dict{String,String}()
    pattern = r"(?i)(?:^|[._-])(license|licence|copying|copyright|notice)(?:[._-]|$)"
    for (directory, subdirectories, files) in walkdir(root)
        filter!(name -> name != ".git", subdirectories)
        for filename in files
            occursin(pattern, filename) || continue
            path = joinpath(directory, filename)
            result[normalize_path(relpath(path, root))] = file_sha256(path)
        end
    end
    return result
end

function command_output(argv::Vector{String})
    stdout_buffer = IOBuffer()
    stderr_buffer = IOBuffer()
    process = run(
        pipeline(ignorestatus(Cmd(argv)); stdout=stdout_buffer, stderr=stderr_buffer)
    )
    if !success(process)
        details = strip(String(take!(stderr_buffer)))
        isempty(details) && (details = strip(String(take!(stdout_buffer))))
        error("command failed ($(join(argv, ' '))): $details")
    end
    return String(take!(stdout_buffer))
end

function require_fields(source, fields)
    identifier = get(source, "id", "<missing-id>")
    missing = [field for field in fields if !haskey(source, field)]
    return isempty(missing) ||
           error("source '$identifier' lacks required fields: $(join(missing, ", "))")
end

function local_path(manifest_root::AbstractString, source, field::AbstractString)
    haskey(source, field) || return nothing
    path = source[field]
    path isa AbstractString || error("$field for $(source["id"]) must be a string")
    return isabspath(path) ? normpath(path) : normpath(joinpath(manifest_root, path))
end

function check_hash!(
    report::Report,
    identifier::AbstractString,
    description::AbstractString,
    actual::AbstractString,
    expected::AbstractString,
)
    if actual == expected
        record!(report, :ok, "$identifier $description SHA-256 $actual")
    else
        record!(
            report,
            :errors,
            "$identifier $description SHA-256 mismatch: expected $expected, got $actual",
        )
    end
end

function verify_git_license_files!(report, source, checkout)
    identifier = source["id"]
    license_path = joinpath(checkout, source["license_file_path"])
    if !isfile(license_path)
        record!(report, :errors, "$identifier missing license file $license_path")
    else
        check_hash!(
            report,
            identifier,
            source["license_file_path"],
            file_sha256(license_path),
            source["license_sha256"],
        )
    end
    extra_paths = get(source, "additional_license_file_paths", Any[])
    extra_hashes = get(source, "additional_license_sha256", Any[])
    if length(extra_paths) != length(extra_hashes)
        record!(
            report, :errors, "$identifier additional license path/hash array lengths differ"
        )
        return nothing
    end
    for (path, expected) in zip(extra_paths, extra_hashes)
        full_path = joinpath(checkout, path)
        if !isfile(full_path)
            record!(report, :errors, "$identifier missing additional license $full_path")
        else
            check_hash!(report, identifier, path, file_sha256(full_path), expected)
        end
    end
end

function semantic_version(tag::AbstractString)
    matched = match(r"^v?(\d+)(?:\.(\d+))?(?:\.(\d+))?$", tag)
    isnothing(matched) && return nothing
    return ntuple(index -> parse(Int, something(matched.captures[index], "0")), 3)
end

function latest_remote_tag(url::AbstractString)
    output = command_output(["git", "ls-remote", "--tags", url])
    direct = Dict{String,String}()
    peeled = Dict{String,String}()
    for line in split(chomp(output), '\n')
        isempty(line) && continue
        parts = split(line)
        length(parts) == 2 || continue
        sha, ref = parts
        startswith(ref, "refs/tags/") || continue
        tag = replace(ref, "refs/tags/" => "")
        if endswith(tag, "^{}")
            peeled[tag[1:(end - 3)]] = sha
        else
            direct[tag] = sha
        end
    end
    candidates = [(semantic_version(tag), tag) for tag in keys(direct)]
    filter!(candidate -> !isnothing(candidate[1]), candidates)
    isempty(candidates) && return nothing
    sort!(candidates; by=candidate -> candidate[1])
    tag = last(candidates)[2]
    return (tag=tag, revision=get(peeled, tag, direct[tag]))
end

function latest_local_tag(checkout::AbstractString)
    output = command_output(["git", "-C", checkout, "tag", "--list"])
    tags = filter(!isempty, split(chomp(output), '\n'))
    candidates = [(semantic_version(tag), tag) for tag in tags]
    filter!(candidate -> !isnothing(candidate[1]), candidates)
    isempty(candidates) && return nothing
    sort!(candidates; by=candidate -> candidate[1])
    tag = last(candidates)[2]
    revision = strip(command_output(["git", "-C", checkout, "rev-parse", "$tag^{}"]))
    return (; tag, revision)
end

function compare_release_tag!(report, identifier, source, latest)
    isnothing(latest) &&
        return record!(report, :warnings, "$identifier has no semantic release tags")
    expected_tag = get(source, "latest_release_tag", "")
    expected_revision = get(source, "latest_release_revision", "")
    if latest.tag != expected_tag || latest.revision != expected_revision
        record!(
            report,
            :updates,
            "$identifier latest release is $(latest.tag) at $(latest.revision); " *
            "manifest pins $expected_tag at $expected_revision",
        )
    else
        record!(
            report,
            :ok,
            "$identifier latest release $(latest.tag) resolves to $(latest.revision)",
        )
    end
end

function function_name_from_signature(signature::AbstractString, fallback::AbstractString)
    isempty(signature) && return String(fallback)
    body = strip(replace(signature, r"^function\s*"i => ""))
    equal_at = findfirst(==('='), body)
    right = isnothing(equal_at) ? body : strip(body[nextind(body, equal_at):end])
    matched = match(r"^([A-Za-z][A-Za-z0-9_]*)", right)
    return isnothing(matched) ? String(fallback) : matched.captures[1]
end

function affected_upstream_functions(manifest_root, source, changed_paths, candidate_files)
    seeds = Set{String}()
    inventory_path = joinpath(manifest_root, "porting", "qetlab_inventory.toml")
    records = Any[]
    if source["id"] == "qetlab" && isfile(inventory_path)
        inventory = TOML.parsefile(inventory_path)
        records = get(inventory, "functions", Any[])
        by_path = Dict(
            record["source_path"] => record["function_name"] for record in records
        )
        for path in changed_paths
            haskey(by_path, path) && push!(seeds, by_path[path])
        end
        changed = true
        while changed
            changed = false
            for record in records
                name = record["function_name"]
                name in seeds && continue
                any(dependency -> dependency in seeds, record["upstream_calls"]) || continue
                push!(seeds, name)
                changed = true
            end
        end
    end
    for path in changed_paths
        haskey(candidate_files, path) || continue
        fallback = splitext(basename(path))[1]
        push!(
            seeds,
            function_name_from_signature(matlab_signature(candidate_files[path]), fallback),
        )
    end
    return sort!(collect(seeds); by=lowercase)
end

function provenance_impact(
    provenance_path::AbstractString, source, changed_paths, affected_upstream
)
    isfile(provenance_path) || return nothing
    provenance = TOML.parsefile(provenance_path)
    entries = get(provenance, "functions", Any[])
    entries isa AbstractVector || error("PROVENANCE.toml `functions` must be an array")
    source_keys = Set(
        lowercase.([source["id"], source["name"], replace(source["name"], ".jl" => "")])
    )
    upstream_keys = Set(lowercase.(affected_upstream))
    for path in changed_paths
        push!(upstream_keys, lowercase(path))
        push!(upstream_keys, lowercase(basename(path)))
        push!(upstream_keys, lowercase(splitext(basename(path))[1]))
    end
    julia_functions = String[]
    tests = String[]
    for entry in entries
        entry isa AbstractDict || continue
        project = lowercase(string(get(entry, "source_project", "")))
        !isempty(project) && !(project in source_keys) && continue
        origins = String[]
        for field in (
            "source_function",
            "source_path",
            "source_file",
            "upstream_file",
            "upstream_function",
        )
            haskey(entry, field) || continue
            value = string(entry[field])
            push!(origins, lowercase(value))
            push!(origins, lowercase(basename(value)))
            push!(origins, lowercase(splitext(basename(value))[1]))
        end
        any(origin -> origin in upstream_keys, origins) || continue
        haskey(entry, "julia_name") && push!(julia_functions, string(entry["julia_name"]))
        if haskey(entry, "tests")
            entry_tests = entry["tests"]
            if entry_tests isa AbstractVector
                append!(tests, string.(entry_tests))
            else
                push!(tests, string(entry_tests))
            end
        end
    end
    return (
        julia_functions=sort!(unique(julia_functions); by=lowercase),
        tests=sort!(unique(tests); by=lowercase),
    )
end

function report_path_set!(report, kind, identifier, description, paths)
    for path in paths
        record!(report, kind, "$identifier $description: $path")
    end
end

function compare_candidate_git!(
    report,
    source,
    baseline::AbstractString,
    candidate::AbstractString,
    manifest_root::AbstractString,
    provenance_path::AbstractString,
)
    identifier = source["id"]
    isdir(candidate) || return record!(
        report, :errors, "$identifier candidate is not a directory: $candidate"
    )
    ispath(joinpath(candidate, ".git")) || return record!(
        report, :errors, "$identifier candidate is not a Git checkout: $candidate"
    )
    candidate_revision = strip(
        command_output(["git", "-C", candidate, "rev-parse", "HEAD"])
    )
    record!(report, :ok, "$identifier candidate revision is $candidate_revision")

    baseline_files = matlab_file_map(baseline)
    candidate_files = matlab_file_map(candidate)
    baseline_paths = Set(keys(baseline_files))
    candidate_paths = Set(keys(candidate_files))
    added = sort!(collect(setdiff(candidate_paths, baseline_paths)); by=lowercase)
    removed = sort!(collect(setdiff(baseline_paths, candidate_paths)); by=lowercase)
    common = intersect(baseline_paths, candidate_paths)
    changed = sort!(
        [path for path in common if baseline_files[path] != candidate_files[path]];
        by=lowercase,
    )
    signature_changed = sort!(
        [
            path for path in changed if matlab_signature(baseline_files[path]) !=
            matlab_signature(candidate_files[path])
        ];
        by=lowercase,
    )
    help_changed = sort!(
        [
            path for path in changed if matlab_help_text(baseline_files[path]) !=
            matlab_help_text(candidate_files[path])
        ];
        by=lowercase,
    )

    baseline_licenses = license_file_map(baseline)
    candidate_licenses = license_file_map(candidate)
    baseline_license_paths = Set(keys(baseline_licenses))
    candidate_license_paths = Set(keys(candidate_licenses))
    licenses_added = sort!(
        collect(setdiff(candidate_license_paths, baseline_license_paths)); by=lowercase
    )
    licenses_removed = sort!(
        collect(setdiff(baseline_license_paths, candidate_license_paths)); by=lowercase
    )
    licenses_changed = sort!(
        [
            path for path in intersect(baseline_license_paths, candidate_license_paths) if
            baseline_licenses[path] != candidate_licenses[path]
        ];
        by=lowercase,
    )

    any_matlab_diff = !isempty(added) || !isempty(removed) || !isempty(changed)
    any_license_diff =
        !isempty(licenses_added) || !isempty(licenses_removed) || !isempty(licenses_changed)
    if any_matlab_diff
        record!(
            report,
            :updates,
            "$identifier candidate MATLAB diff: $(length(added)) added, " *
            "$(length(removed)) removed, $(length(changed)) changed; " *
            "$(length(signature_changed)) signatures and $(length(help_changed)) help blocks changed",
        )
        report_path_set!(report, :updates, identifier, "added MATLAB file", added)
        report_path_set!(report, :updates, identifier, "removed MATLAB file", removed)
        report_path_set!(report, :updates, identifier, "changed MATLAB file", changed)
        report_path_set!(
            report, :updates, identifier, "changed MATLAB signature", signature_changed
        )
        report_path_set!(
            report, :updates, identifier, "changed MATLAB help text", help_changed
        )
    else
        record!(
            report,
            :ok,
            "$identifier candidate MATLAB diff is empty (0 added, 0 removed, 0 changed)",
        )
    end
    if any_license_diff
        record!(
            report,
            :updates,
            "$identifier candidate license diff: $(length(licenses_added)) added, " *
            "$(length(licenses_removed)) removed, $(length(licenses_changed)) changed",
        )
        report_path_set!(report, :updates, identifier, "added license file", licenses_added)
        report_path_set!(
            report, :updates, identifier, "removed license file", licenses_removed
        )
        report_path_set!(
            report, :updates, identifier, "changed license file", licenses_changed
        )
    else
        record!(report, :ok, "$identifier candidate license-file diff is empty")
    end

    any_matlab_diff || return nothing
    changed_paths = sort!(unique(vcat(added, removed, changed)); by=lowercase)
    affected_upstream = affected_upstream_functions(
        manifest_root, source, changed_paths, candidate_files
    )
    if isempty(affected_upstream)
        record!(
            report,
            :warnings,
            "$identifier no affected upstream functions could be resolved",
        )
    else
        record!(
            report,
            :updates,
            "$identifier affected upstream functions (including reverse dependencies): " *
            join(affected_upstream, ", "),
        )
    end
    impact = provenance_impact(provenance_path, source, changed_paths, affected_upstream)
    if isnothing(impact)
        record!(
            report,
            :warnings,
            "$identifier PROVENANCE.toml unavailable at $provenance_path; " *
            "affected Julia functions and tests cannot be mapped",
        )
    elseif isempty(impact.julia_functions)
        record!(
            report,
            :warnings,
            "$identifier no Julia functions in PROVENANCE.toml map to the candidate changes",
        )
    else
        record!(
            report,
            :updates,
            "$identifier affected Julia functions: " * join(impact.julia_functions, ", "),
        )
        if isempty(impact.tests)
            record!(
                report,
                :warnings,
                "$identifier affected provenance entries list no tests to rerun",
            )
        else
            record!(
                report, :updates, "$identifier tests to rerun: " * join(impact.tests, ", ")
            )
        end
    end
    return record!(
        report,
        :updates,
        "$identifier regenerate porting/qetlab_inventory.toml, CSV, and dependency graph before adapting code",
    )
end

function check_git_source!(
    report, source, manifest_root, offline, candidate, provenance_path
)
    require_fields(
        source,
        [
            "id",
            "clone_url",
            "default_branch",
            "exact_revision",
            "license_file_path",
            "license_sha256",
        ],
    )
    identifier = source["id"]
    checkout = local_path(manifest_root, source, "local_audit_path")
    if !isnothing(checkout) && isdir(checkout)
        head = strip(command_output(["git", "-C", checkout, "rev-parse", "HEAD"]))
        if head == source["exact_revision"]
            record!(report, :ok, "$identifier local checkout is pinned at $head")
        else
            record!(
                report,
                :errors,
                "$identifier local checkout is $head, expected $(source["exact_revision"])",
            )
        end
        dirty = command_output(["git", "-C", checkout, "status", "--porcelain"])
        if isempty(strip(dirty))
            record!(report, :ok, "$identifier local checkout is clean")
        else
            record!(report, :errors, "$identifier local checkout has uncommitted changes")
        end
        verify_git_license_files!(report, source, checkout)
    elseif offline
        record!(
            report,
            :errors,
            "$identifier offline checkout is unavailable at $(something(checkout, "<unspecified>"))",
        )
        return nothing
    else
        record!(
            report,
            :warnings,
            "$identifier local checkout unavailable; pinned license bytes were not rehashed",
        )
    end

    if !isnothing(candidate)
        if isnothing(checkout) || !isdir(checkout)
            record!(
                report,
                :errors,
                "$identifier candidate comparison requires the pinned local checkout",
            )
        else
            compare_candidate_git!(
                report, source, checkout, candidate, manifest_root, provenance_path
            )
        end
    end

    if offline
        isnothing(checkout) && return nothing
        compare_release_tag!(report, identifier, source, latest_local_tag(checkout))
        return nothing
    end

    remote = strip(
        command_output([
            "git",
            "ls-remote",
            "--exit-code",
            source["clone_url"],
            "refs/heads/$(source["default_branch"])",
        ]),
    )
    remote_revision = first(split(remote))
    if remote_revision == source["exact_revision"]
        record!(
            report,
            :ok,
            "$identifier $(source["default_branch"]) remains at $remote_revision",
        )
    else
        record!(
            report,
            :updates,
            "$identifier $(source["default_branch"]) moved from " *
            "$(source["exact_revision"]) to $remote_revision",
        )
    end
    return compare_release_tag!(
        report, identifier, source, latest_remote_tag(source["clone_url"])
    )
end

function zip_member(archive::AbstractString, member::AbstractString)
    return read(Cmd(["unzip", "-p", archive, member]))
end

function verify_archive_metadata!(
    report,
    source,
    archive,
    manifest_root;
    archive_mismatch_kind::Symbol=:errors,
    emit_known_blocker::Bool=true,
)
    identifier = source["id"]
    actual_archive_hash = file_sha256(archive)
    if actual_archive_hash == source["archive_sha256"]
        record!(report, :ok, "$identifier archive SHA-256 $actual_archive_hash")
    else
        record!(
            report,
            archive_mismatch_kind,
            "$identifier archive SHA-256 changed: expected $(source["archive_sha256"]), " *
            "got $actual_archive_hash",
        )
    end
    if haskey(source, "archive_size_bytes")
        actual_size = filesize(archive)
        expected_size = source["archive_size_bytes"]
        if actual_size == expected_size
            record!(report, :ok, "$identifier archive size is $actual_size bytes")
        else
            record!(
                report,
                :updates,
                "$identifier archive size changed from $expected_size to $actual_size bytes",
            )
        end
    end
    for (path_field, hash_field, label) in [
        ("license_file_path", "license_sha256", "bundled license"),
        ("readme_file_path", "readme_sha256", "README"),
    ]
        haskey(source, path_field) || continue
        try
            actual = bytes_sha256(zip_member(archive, source[path_field]))
            expected = source[hash_field]
            if actual == expected
                record!(report, :ok, "$identifier $label SHA-256 $actual")
            elseif actual_archive_hash == source["archive_sha256"]
                record!(
                    report,
                    :errors,
                    "$identifier pinned archive has unexpected $label hash $actual (expected $expected)",
                )
            else
                record!(
                    report,
                    :updates,
                    "$identifier updated archive has changed $label hash $actual (pinned $expected)",
                )
            end
        catch exception
            record!(
                report,
                :errors,
                "$identifier could not read only permitted $label member: $(sprint(showerror, exception))",
            )
        end
    end
    if haskey(source, "supplied_license_file_path")
        supplied = local_path(manifest_root, source, "supplied_license_file_path")
        if isnothing(supplied) || !isfile(supplied)
            record!(report, :errors, "$identifier supplied license file is missing")
        else
            check_hash!(
                report,
                identifier,
                "supplied license",
                file_sha256(supplied),
                source["supplied_license_sha256"],
            )
        end
    end
    status = get(source, "license_verification_status", "")
    if emit_known_blocker && startswith(status, "MATERIAL_CONFLICT")
        for blocker in get(source, "blockers", Any[])
            record!(report, :blockers, "$identifier: $blocker")
        end
    end
end

function check_archive_source!(report, source, manifest_root, offline, candidate)
    require_fields(
        source,
        ["id", "archive_url", "archive_sha256", "license_file_path", "license_sha256"],
    )
    identifier = source["id"]
    local_archive = local_path(manifest_root, source, "local_audit_path")
    candidate_matches_local =
        !isnothing(candidate) &&
        !isnothing(local_archive) &&
        isfile(candidate) &&
        isfile(local_archive) &&
        samefile(candidate, local_archive)
    if !isnothing(candidate)
        if !isfile(candidate)
            record!(
                report, :errors, "$identifier candidate archive does not exist: $candidate"
            )
        else
            record!(
                report,
                :blockers,
                "$identifier candidate comparison is restricted to archive/license/README metadata; blocked .m content was not read",
            )
            if !candidate_matches_local
                verify_archive_metadata!(
                    report,
                    source,
                    candidate,
                    manifest_root;
                    archive_mismatch_kind=:updates,
                    emit_known_blocker=false,
                )
            end
        end
    end
    if offline
        if isnothing(local_archive) || !isfile(local_archive)
            record!(
                report,
                :errors,
                "$identifier offline archive unavailable at $(something(local_archive, "<unspecified>"))",
            )
            return nothing
        end
        verify_archive_metadata!(report, source, local_archive, manifest_root)
        return nothing
    end

    temporary_path, temporary_io = mktemp()
    close(temporary_io)
    try
        Downloads.download(source["archive_url"], temporary_path)
        actual_hash = file_sha256(temporary_path)
        if actual_hash == source["archive_sha256"]
            record!(report, :ok, "$identifier remote archive remains at $actual_hash")
        else
            record!(
                report,
                :updates,
                "$identifier remote archive changed from $(source["archive_sha256"]) to $actual_hash",
            )
        end
        verify_archive_metadata!(
            report, source, temporary_path, manifest_root; archive_mismatch_kind=:updates
        )
    finally
        isfile(temporary_path) && rm(temporary_path)
    end
end

function resolve_candidates(specs, selected_sources, sources)
    result = Dict{String,String}()
    known = Set(source["id"] for source in sources)
    for spec in specs
        if occursin("=", spec)
            identifier, path = split(spec, "="; limit=2)
        else
            if length(selected_sources) == 1
                identifier = only(selected_sources)
            elseif length(sources) == 1
                identifier = only(sources)["id"]
            else
                error("candidate PATH without ID requires exactly one --source")
            end
            path = spec
        end
        identifier in known || error("candidate names unknown source ID '$identifier'")
        isempty(path) && error("candidate path for '$identifier' is empty")
        haskey(result, identifier) && error("duplicate candidate for '$identifier'")
        result[identifier] = abspath(path)
    end
    return result
end

function main(args)
    options = parse_options(args)
    isfile(options.manifest) || error("manifest does not exist: $(options.manifest)")
    manifest = TOML.parsefile(options.manifest)
    get(manifest, "schema_version", nothing) == 1 ||
        error("unsupported or missing manifest schema_version")
    sources = get(manifest, "sources", nothing)
    sources isa AbstractVector || error("manifest must contain [[sources]] tables")
    known_ids = Set(get(source, "id", "") for source in sources)
    unknown_selected = setdiff(Set(options.selected_sources), known_ids)
    isempty(unknown_selected) ||
        error("unknown --source IDs: $(join(sort!(collect(unknown_selected)), ", "))")
    if !isempty(options.selected_sources)
        selected = Set(options.selected_sources)
        sources = [source for source in sources if source["id"] in selected]
    end
    candidates = resolve_candidates(
        options.candidate_specs, options.selected_sources, sources
    )

    report = Report()
    manifest_root = dirname(options.manifest)
    for source in sources
        identifier = get(source, "id", "<missing-id>")
        try
            kind = get(source, "source_kind", "")
            candidate = get(candidates, identifier, nothing)
            if kind == "git"
                check_git_source!(
                    report,
                    source,
                    manifest_root,
                    options.offline,
                    candidate,
                    options.provenance,
                )
            elseif kind == "archive"
                check_archive_source!(
                    report, source, manifest_root, options.offline, candidate
                )
            else
                record!(report, :errors, "$identifier has unsupported source_kind '$kind'")
            end
        catch exception
            record!(
                report, :errors, "$identifier check failed: $(sprint(showerror, exception))"
            )
        end
    end
    println(
        "summary: ",
        length(report.ok),
        " ok, ",
        length(report.updates),
        " updates, ",
        length(report.warnings),
        " warnings, ",
        length(report.blockers),
        " known blockers, ",
        length(report.errors),
        " errors",
    )
    !isempty(report.errors) && exit(2)
    return options.fail_on_update && !isempty(report.updates) && exit(1)
end

main(ARGS)
