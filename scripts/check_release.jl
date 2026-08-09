#!/usr/bin/env julia

using Dates
using SHA
using Tar
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const REQUIRED_DISTRIBUTION_FILES = (
    "Project.toml",
    "README.md",
    "LICENSE",
    "NOTICE",
    "THIRD_PARTY_LICENSES.md",
    "UpstreamManifest.toml",
    "PROVENANCE.toml",
    "CITATION.cff",
    "CITATION.bib",
    "CHANGELOG.md",
    "SECURITY.md",
    "SUPPORT.md",
    "artifacts/convergence/current_snapshot.toml",
    "docs/BENCHMARK_REPORT.md",
    "docs/BUILD_ENVIRONMENT.md",
    "docs/CONVERGENCE_AUDIT.md",
    "docs/FULL_QETLAB_BASELINE.md",
    "docs/JULIA_CORE_AUDIT.md",
    "docs/LEGAL.md",
    "docs/PORTING_STATUS.md",
    "docs/QETLAB_COMPLETION_PLAN.md",
    "docs/RELEASE_CHECKLIST.md",
    "docs/SESSION_HANDOFF.md",
    "docs/VALIDATION_REPORT.md",
    "docs/adr/0006-general-operator-space-maps.md",
    "ext/QuantumEntanglementToolsEntanglementDetectionExt.jl",
    "ext/entanglement_detection_worker.jl",
    "licenses/QETLAB-LICENSE.txt",
    "licenses/QETLAB-BRUNO-LUONG-HELPERS-LICENSE.txt",
    "porting/qetlab_completion_plan.toml",
    "porting/qetlab_completion_policy.toml",
    "porting/qetlab_completion_queue.toml",
    "porting/qetlab_inventory.toml",
    "porting/qetlab_status.toml",
    "scripts/build_qetlab_completion_plan.jl",
    "scripts/check_qetlab_completion.jl",
    "scripts/check_release.jl",
    "scripts/qetlab_completion_common.jl",
    "scripts/reconcile_project_claims.jl",
    "test/oracle/runtests.jl",
    "test/release_gate_checks.jl",
    "test/runtests.jl",
)
const FORBIDDEN_DISTRIBUTION_PATHS = (
    r"(^|/)\.DS_Store$",
    r"(^|/)(?:Julia)?Manifest(?:-v\d+\.\d+)?\.toml$",
    r"^benchmark/results/",
    r"^coverage/",
    r"^dev/(?:downloads|upstream)/",
    r"^docs/build/",
    r"^test/oracle/generated/",
    r"(?:^|/)[^/]*\.cov$",
    r"\.code-workspace$",
)

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/check_release.jl [options]

Validate unreleased or tagged-release metadata and the exact Git archive.

Options:
  --treeish REF      Validate REF (default: HEAD)
  --release-candidate
                     Require dated release metadata before creating the tag
  --tag TAG          Require an annotated v<version> tag and dated release metadata
  --archive-smoke    Extract REF and load/smoke-test it in a fresh Julia depot
  --registry         Run a partial General-registry preflight
  --allow-dirty      Validate an isolated archive of current worktree bytes
  -h, --help         Show this help

Without `--allow-dirty`, checks read bytes from the exact committed tree and
reject tracked changes and unexpected untracked files. Dirty mode uses an
isolated temporary Git index and object store, includes nonignored untracked
paths, and never changes the repository index. It is local preflight evidence,
not release evidence. Without `--release-candidate` or `--tag`, the changelog
and citation metadata must describe an unreleased development milestone.
Release-candidate mode accepts the dated metadata that must be committed and
tested before an annotated tag exists, but it is not tagged-release evidence.
`--registry` is only a partial preflight: it cannot prove the maintainer's
non-delegable review, remote visibility, or RegistryCI acceptance.
""",
    )
end

function parse_options(args)
    treeish = "HEAD"
    tag = nothing
    registry = false
    archive_smoke = false
    allow_dirty = false
    release_candidate = false
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument in ("-h", "--help")
            usage()
            exit(0)
        elseif argument == "--registry"
            registry = true
        elseif argument == "--archive-smoke"
            archive_smoke = true
        elseif argument == "--allow-dirty"
            allow_dirty = true
        elseif argument == "--release-candidate"
            release_candidate = true
        elseif startswith(argument, "--treeish=")
            treeish = split(argument, "="; limit=2)[2]
        elseif argument == "--treeish"
            index == length(args) && error("--treeish requires a value")
            index += 1
            treeish = args[index]
        elseif startswith(argument, "--tag=")
            tag = split(argument, "="; limit=2)[2]
        elseif argument == "--tag"
            index == length(args) && error("--tag requires a value")
            index += 1
            tag = args[index]
        else
            error("unknown option: $argument")
        end
        index += 1
    end
    if isnothing(tag) &&
        get(ENV, "GITHUB_REF_TYPE", "") == "tag" &&
        !isempty(get(ENV, "GITHUB_REF_NAME", ""))
        tag = ENV["GITHUB_REF_NAME"]
    end
    release_candidate &&
        !isnothing(tag) &&
        error("--release-candidate and --tag are mutually exclusive")
    return (; treeish, tag, registry, archive_smoke, allow_dirty, release_candidate)
end

function git_command(arguments...)
    command = String["git", "-C", REPOSITORY_ROOT]
    append!(command, string.(arguments))
    return Cmd(command)
end

git_bytes(arguments...) = read(git_command(arguments...))
git_string(arguments...) = read(git_command(arguments...), String)

function try_git_string(arguments...)
    command = pipeline(git_command(arguments...); stderr=devnull)
    try
        return strip(read(command, String))
    catch error
        error isa ProcessFailedException || rethrow()
        return nothing
    end
end

function split_lines(output::AbstractString)
    stripped = chomp(output)
    return isempty(stripped) ? String[] : String.(split(stripped, '\n'))
end

git_lines(arguments...) = split_lines(git_string(arguments...))

function is_valid_iso_date(value::AbstractString)
    try
        return string(Date(value, dateformat"yyyy-mm-dd")) == value
    catch
        return false
    end
end

function validate_release_metadata!(
    failures::Vector{String},
    version::VersionNumber,
    citation::AbstractString,
    changelog::AbstractString;
    dated::Bool,
)
    check(condition, message) = condition || push!(failures, message)
    escaped_version = replace(string(version), "." => "\\.")
    release_match = match(
        Regex("(?m)^## \\[$escaped_version\\] - (\\d{4}-\\d{2}-\\d{2})\\s*\$"), changelog
    )
    citation_date_match = match(
        r"(?m)^date-released:\s*[\"']?(\d{4}-\d{2}-\d{2})[\"']?\s*$", citation
    )

    check(
        occursin(r"(?m)^## \[Unreleased\]\s*$", changelog),
        "CHANGELOG.md has no Unreleased heading",
    )
    check(
        occursin(Regex("(?m)^version:\\s*[\"']?$escaped_version[\"']?\\s*\$"), citation),
        "CITATION.cff version does not match Project.toml",
    )
    if dated
        check(
            !isnothing(release_match),
            "CHANGELOG.md has no dated release heading for $version",
        )
        if !isnothing(release_match)
            check(
                is_valid_iso_date(release_match.captures[1]),
                "CHANGELOG.md release date is not a valid ISO calendar date",
            )
        end
        check(!isnothing(citation_date_match), "CITATION.cff has no release date")
        if !isnothing(citation_date_match)
            check(
                is_valid_iso_date(citation_date_match.captures[1]),
                "CITATION.cff release date is not a valid ISO calendar date",
            )
        end
        if !isnothing(release_match) && !isnothing(citation_date_match)
            check(
                release_match.captures[1] == citation_date_match.captures[1],
                "CHANGELOG.md and CITATION.cff release dates disagree",
            )
        end
    else
        check(
            isnothing(release_match),
            "unreleased metadata must not contain a dated $version release heading",
        )
        check(
            isnothing(citation_date_match),
            "unreleased CITATION.cff must not contain date-released",
        )
    end
    check(
        occursin(r"(?m)^repository-code:\s*[\"']?https://github\.com/", citation),
        "CITATION.cff has no canonical repository-code URL",
    )
    return nothing
end

function validate_release_copy!(
    failures::Vector{String}, version::VersionNumber, documents::AbstractDict
)
    check(condition, message) = condition || push!(failures, message)
    required_paths = (
        "AGENTS.md",
        "CONTRIBUTING.md",
        "README.md",
        "SECURITY.md",
        "CITATION.cff",
        "CITATION.bib",
        "docs/CONVERGENCE_AUDIT.md",
        "docs/PORTING_STATUS.md",
        "docs/SESSION_HANDOFF.md",
        "docs/VALIDATION_REPORT.md",
        "docs/src/getting_started.md",
    )
    for path in required_paths
        check(haskey(documents, path), "release-copy check has no $path")
    end
    all(path -> haskey(documents, path), required_paths) || return nothing

    stale_patterns = (
        r"(?i)\bunreleased\b",
        r"(?i)no version has been (?:released|published|tagged)",
        r"(?i)no release, tag, or publication has been made",
    )
    for path in required_paths, pattern in stale_patterns
        checked_content =
            if path in (
                "AGENTS.md",
                "CONTRIBUTING.md",
                "docs/CONVERGENCE_AUDIT.md",
                "docs/PORTING_STATUS.md",
                "docs/SESSION_HANDOFF.md",
                "docs/VALIDATION_REPORT.md",
            )
                # These living maintenance documents retain historical records and
                # instructions that may legitimately mention the Changelog's
                # Unreleased section. Only their release-facing preamble/current
                # status is required to be tag-safe.
                join(Iterators.take(eachline(IOBuffer(documents[path])), 120), '\n')
            else
                documents[path]
            end
        check(
            !occursin(pattern, checked_content),
            "$path retains pre-release wording matching $(repr(pattern))",
        )
    end

    expected_revision = "rev=\"v$version\""
    for path in ("README.md", "docs/src/getting_started.md")
        check(
            occursin(expected_revision, documents[path]),
            "$path must contain the stable installation revision $expected_revision",
        )
    end
    return nothing
end

function split_nul_records(bytes::AbstractVector{UInt8})
    records = Vector{Vector{UInt8}}()
    start = firstindex(bytes)
    for index in eachindex(bytes)
        bytes[index] == 0x00 || continue
        start < index && push!(records, collect(@view bytes[start:(index - 1)]))
        start = index + 1
    end
    start <= lastindex(bytes) && push!(records, collect(@view bytes[start:end]))
    return records
end

function inspect_candidate_bytes!(
    failures::Vector{String},
    path::AbstractString,
    bytes::AbstractVector{UInt8},
    source::AbstractString,
    lfs_pointer_prefix::AbstractVector{UInt8},
)
    length(bytes) >= length(lfs_pointer_prefix) &&
        bytes[1:length(lfs_pointer_prefix)] == lfs_pointer_prefix &&
        push!(
            failures, "Git LFS pointer is distributed instead of content: $path ($source)"
        )
    return nothing
end

function candidate_content(
    path::AbstractString, archive_root::AbstractString, candidate_paths::Set{String}
)
    path in candidate_paths || return ""
    full_path = joinpath(archive_root, path)
    return isfile(full_path) ? read(full_path, String) : ""
end

function dirty_archive_bytes(resolved_commit::AbstractString)
    return mktempdir() do temporary_root
        index_path = joinpath(temporary_root, "index")
        object_path = joinpath(temporary_root, "objects")
        mkpath(object_path)
        git_directory = strip(git_string("rev-parse", "--absolute-git-dir"))
        environment = (
            "GIT_INDEX_FILE" => index_path,
            "GIT_OBJECT_DIRECTORY" => object_path,
            "GIT_ALTERNATE_OBJECT_DIRECTORIES" => joinpath(git_directory, "objects"),
        )
        run(addenv(git_command("read-tree", resolved_commit), environment...))
        run(addenv(git_command("add", "-A", "--", "."), environment...))
        tree = strip(read(addenv(git_command("write-tree"), environment...), String))
        return read(addenv(git_command("archive", "--format=tar", tree), environment...))
    end
end

function run_archive_smoke(archive_bytes::Vector{UInt8})
    return mktempdir() do temporary_root
        source_root = joinpath(temporary_root, "source")
        depot_root = joinpath(temporary_root, "depot")
        mkpath(source_root)
        mkpath(depot_root)
        Tar.extract(IOBuffer(archive_bytes), source_root)
        program = """
        using Pkg
        Pkg.instantiate()
        using QuantumEntanglementTools
        using LinearAlgebra
        psi = bell_state()
        rho_a = partial_trace(psi, (2, 2); trace_out=(2,))
        @assert isapprox(rho_a, Matrix{ComplexF64}(I, 2, 2) / 2)
        """
        command = addenv(
            `$(Base.julia_cmd()) --startup-file=no --project=$source_root -e $program`,
            "JULIA_DEPOT_PATH" => depot_root,
        )
        return run(command)
    end
end

function check_completion_artifacts(candidate_root::AbstractString)
    script = joinpath(candidate_root, "scripts", "build_qetlab_completion_plan.jl")
    isfile(script) || return (
        success=false, output="candidate has no scripts/build_qetlab_completion_plan.jl"
    )
    output = IOBuffer()
    command = `$(Base.julia_cmd()) --startup-file=no --project=$candidate_root $script --check`
    success = try
        run(pipeline(command; stdout=output, stderr=output))
        true
    catch error
        error isa ProcessFailedException || rethrow()
        false
    end
    return (; success, output=strip(String(take!(output))))
end

function check_release(options)
    failures = String[]
    check(condition, message) = condition || push!(failures, message)

    resolved_commit = try_git_string("rev-parse", "--verify", "$(options.treeish)^{commit}")
    if isnothing(resolved_commit)
        println(stderr, "release check failed: invalid treeish $(repr(options.treeish))")
        exit(1)
    end

    archive_bytes = if options.allow_dirty
        dirty_archive_bytes(resolved_commit)
    else
        git_bytes("archive", "--format=tar", resolved_commit)
    end
    archive_digest = bytes2hex(sha256(archive_bytes))
    archive_headers = Tar.list(IOBuffer(archive_bytes))
    archive_root = mktempdir()
    Tar.extract(IOBuffer(archive_bytes), archive_root)
    archive_paths = Set(
        String(header.path) for header in archive_headers if header.type == :file
    )
    candidate_paths = copy(archive_paths)
    content(path) = candidate_content(path, archive_root, candidate_paths)

    status_lines = git_lines("status", "--porcelain=v1", "--untracked-files=all")
    for line in status_lines
        if startswith(line, "?? ")
            path = line[4:end]
            options.allow_dirty || push!(failures, "unexpected untracked path: $path")
        else
            status_code = line[1:2]
            unsafe_dirty_status =
                occursin('D', status_code) ||
                occursin('T', status_code) ||
                status_code in ("AA", "AU", "DD", "DU", "UA", "UD", "UU")
            if unsafe_dirty_status
                push!(
                    failures,
                    "worktree contains a deletion, type change, or unresolved entry: $line",
                )
            elseif !options.allow_dirty
                push!(failures, "tracked worktree changes are present")
            end
        end
    end

    for path in REQUIRED_DISTRIBUTION_FILES
        check(path in candidate_paths, "required file is absent from candidate: $path")
        check(!isempty(content(path)), "required file is absent or empty: $path")
        if !options.allow_dirty
            check(
                path in archive_paths, "required file is excluded from Git archive: $path"
            )
        end
    end

    completion_artifacts = check_completion_artifacts(archive_root)
    if !completion_artifacts.success
        detail = if isempty(completion_artifacts.output)
            "no diagnostic output"
        else
            replace(completion_artifacts.output, '\n' => " | ")
        end
        push!(
            failures, "generated QETLAB completion artifacts are missing or stale: $detail"
        )
    end

    for path in candidate_paths, pattern in FORBIDDEN_DISTRIBUTION_PATHS
        occursin(pattern, path) &&
            push!(failures, "generated or editor-local path is distributed: $path")
    end
    tree_records = split_nul_records(
        git_bytes("ls-tree", "-r", "-z", "--full-tree", resolved_commit)
    )

    project = try
        TOML.parse(content("Project.toml"))
    catch error
        push!(failures, "Project.toml cannot be parsed: $(sprint(showerror, error))")
        Dict{String,Any}()
    end
    name = get(project, "name", "")
    version_text = get(project, "version", "")
    check(name == "QuantumEntanglementTools", "unexpected package name: $(repr(name))")
    check(
        version_text isa AbstractString && occursin(r"^\d+\.\d+\.\d+$", version_text),
        "Project.toml package version must have exactly three numeric components",
    )
    version = try
        VersionNumber(version_text)
    catch
        nothing
    end
    check(version isa VersionNumber, "Project.toml version is not valid SemVer")
    if version isa VersionNumber
        check(
            isempty(version.prerelease), "package version must not contain prerelease data"
        )
        check(isempty(version.build), "package version must not contain build data")
    end
    uuid = get(project, "uuid", "")
    check(
        uuid isa AbstractString && occursin(
            r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", uuid
        ),
        "Project.toml uuid is missing or malformed",
    )
    authors = get(project, "authors", Any[])
    check(authors isa AbstractVector && !isempty(authors), "Project.toml authors is empty")

    compat = get(project, "compat", Dict{String,Any}())
    dependencies = Set{String}()
    for section in ("deps", "weakdeps", "extras")
        union!(dependencies, keys(get(project, section, Dict{String,Any}())))
    end
    for dependency in sort!(collect(dependencies))
        check(haskey(compat, dependency), "missing [compat] entry for $dependency")
    end
    check(haskey(compat, "julia"), "missing [compat] entry for julia")

    extensions = get(project, "extensions", Dict{String,Any}())
    weak_dependencies = get(project, "weakdeps", Dict{String,Any}())
    for (extension, dependency) in extensions
        extension_path = "ext/$extension.jl"
        check(
            extension_path in candidate_paths,
            "declared extension entry point is absent: $extension_path",
        )
        names = dependency isa AbstractVector ? dependency : [dependency]
        for dependency_name in names
            check(
                dependency_name isa AbstractString &&
                    haskey(weak_dependencies, dependency_name),
                "extension $extension names undeclared weak dependency $(repr(dependency_name))",
            )
        end
    end
    entrypoint = "src/$name.jl"
    check(entrypoint in candidate_paths, "package entry point is absent: $entrypoint")

    lfs_pointer_prefix = collect(codeunits("version https://git-lfs.github.com/spec/v1"))
    for path in candidate_paths
        full_path = joinpath(archive_root, path)
        bytes = isfile(full_path) ? read(full_path) : UInt8[]
        inspect_candidate_bytes!(
            failures,
            path,
            bytes,
            options.allow_dirty ? "isolated worktree archive" : "archive",
            lfs_pointer_prefix,
        )
    end
    for header in archive_headers
        header.type == :symlink && push!(
            failures, "symbolic link requires explicit review: $(header.path) (archive)"
        )
    end
    for record in tree_records
        tab_index = findfirst(==(0x09), record)
        if isnothing(tab_index)
            push!(failures, "malformed NUL-delimited git tree entry")
            continue
        end
        metadata_bytes = @view record[1:(tab_index - 1)]
        path_bytes = @view record[(tab_index + 1):end]
        if !isvalid(String, metadata_bytes) || !isvalid(String, path_bytes)
            push!(failures, "git tree contains a non-UTF-8 path or metadata entry")
            continue
        end
        metadata = String(metadata_bytes)
        path = String(path_bytes)
        for pattern in FORBIDDEN_DISTRIBUTION_PATHS
            occursin(pattern, path) && push!(
                failures, "generated or editor-local path is tracked: $path (Git tree)"
            )
        end
        fields = split(metadata, ' ')
        if length(fields) != 3
            push!(failures, "malformed git tree metadata for $(repr(path))")
            continue
        end
        mode, kind, object_id = fields
        mode == "160000" &&
            push!(failures, "git submodules are not permitted in the release tree: $path")
        mode == "120000" &&
            push!(failures, "symbolic link requires explicit review: $path (Git tree)")
        kind == "blob" || continue
        bytes = git_bytes("cat-file", "blob", object_id)
        inspect_candidate_bytes!(failures, path, bytes, "Git tree", lfs_pointer_prefix)
    end

    license = content("LICENSE")
    check(
        occursin("BSD 3-Clause License", license),
        "LICENSE does not identify the declared BSD-3-Clause license",
    )
    notice = content("NOTICE")
    check(
        occursin("QETLAB", notice) && occursin("EntanglementDetection.jl", notice),
        "NOTICE is missing an audited upstream/dependency notice",
    )

    readme = content("README.md")
    check(
        occursin("OpenAI Codex", readme),
        "README must disclose the repository's substantial Codex assistance",
    )
    check(
        occursin(
            r"(?is)(not\s+a\s+claim\s+of|does\s+not\s+claim|no)\s+(complete\s+)?QETLAB\s+parity",
            readme,
        ),
        "README must retain an explicit no-QETLAB-parity statement",
    )

    citation = content("CITATION.cff")
    changelog = content("CHANGELOG.md")
    if version isa VersionNumber
        validate_release_metadata!(
            failures,
            version,
            citation,
            changelog;
            dated=options.release_candidate || !isnothing(options.tag),
        )
        if options.release_candidate || !isnothing(options.tag)
            release_copy_paths = (
                "AGENTS.md",
                "CONTRIBUTING.md",
                "README.md",
                "SECURITY.md",
                "CITATION.cff",
                "CITATION.bib",
                "docs/CONVERGENCE_AUDIT.md",
                "docs/PORTING_STATUS.md",
                "docs/SESSION_HANDOFF.md",
                "docs/VALIDATION_REPORT.md",
                "docs/src/getting_started.md",
            )
            validate_release_copy!(
                failures,
                version,
                Dict(path => content(path) for path in release_copy_paths),
            )
        end
    end

    if !isnothing(options.tag)
        expected_tag = version isa VersionNumber ? "v$version" : ""
        check(options.tag == expected_tag, "tag $(options.tag) must equal $expected_tag")
        tag_commit = try_git_string(
            "rev-parse", "--verify", "refs/tags/$(options.tag)^{commit}"
        )
        check(!isnothing(tag_commit), "tag $(options.tag) does not exist locally")
        tag_object_type = try_git_string("cat-file", "-t", "refs/tags/$(options.tag)")
        check(
            tag_object_type == "tag", "tag $(options.tag) must be an annotated tag object"
        )
        if !isnothing(tag_commit)
            check(
                tag_commit == resolved_commit,
                "tag $(options.tag) does not target candidate commit $resolved_commit",
            )
        end
    end

    if options.registry
        check(
            get(ENV, "QET_HUMAN_REVIEW_CONFIRMED", "") == "true",
            "maintainer human-review attestation is absent; an environment flag is not independent proof",
        )
        expected_url = "https://github.com/aenictusGitHub/$name.jl"
        repository_match = match(
            r"(?m)^repository-code:\s*[\"']?([^\"'\s]+)[\"']?\s*$", citation
        )
        repository_url = isnothing(repository_match) ? "" : repository_match.captures[1]
        check(
            repository_url == expected_url,
            "partial General preflight expects repository-code: $expected_url",
        )
        origin = strip(git_string("remote", "get-url", "origin"))
        normalized_origin = replace(origin, r"\.git$" => "")
        check(
            normalized_origin == expected_url,
            "partial General preflight expects origin $expected_url",
        )
    end

    if isempty(failures) && options.archive_smoke
        try
            run_archive_smoke(archive_bytes)
        catch error
            push!(failures, "fresh-depot archive smoke failed: $(sprint(showerror, error))")
        end
    end

    if isempty(failures)
        if options.allow_dirty
            println(
                if options.release_candidate
                    "working-tree release-candidate preflight passed for "
                else
                    "working-tree preflight passed for "
                end,
                name,
                " v",
                version_text,
                " (not release evidence)",
            )
        else
            if isnothing(options.tag)
                println(
                    if options.release_candidate
                        "release-candidate archive check passed for "
                    else
                        "unreleased archive check passed for "
                    end,
                    name,
                    " v",
                    version_text,
                    " at ",
                    resolved_commit,
                    "; tar SHA-256 ",
                    archive_digest,
                    " (not tagged release evidence)",
                    options.registry ? " (partial General preflight)" : "",
                )
            else
                println(
                    "release archive check passed for ",
                    name,
                    " v",
                    version_text,
                    " at ",
                    resolved_commit,
                    "; tar SHA-256 ",
                    archive_digest,
                    options.registry ? " (partial General preflight)" : "",
                )
            end
        end
        return nothing
    end

    println(stderr, "release check failed:")
    for failure in unique(failures)
        println(stderr, "  - ", failure)
    end
    return exit(1)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    check_release(parse_options(ARGS))
end
