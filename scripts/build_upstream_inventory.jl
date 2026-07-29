#!/usr/bin/env julia

using Dates
using SHA
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_SOURCE = joinpath(REPOSITORY_ROOT, "dev", "upstream", "QETLAB")
const DEFAULT_OUTPUT = joinpath(REPOSITORY_ROOT, "porting")
const DEFAULT_STATUS = joinpath(REPOSITORY_ROOT, "porting", "qetlab_status.toml")
const GENERATOR_VERSION = 1

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/build_upstream_inventory.jl [options] [QETLAB_CHECKOUT]

Build the authoritative QETLAB MATLAB-source inventory and dependency graph.

Options:
  --source PATH       QETLAB checkout (default: dev/upstream/QETLAB)
  --output-dir PATH   output directory (default: porting)
  --status PATH       reviewed status overlay (default: porting/qetlab_status.toml)
  --date YYYY-MM-DD   audit date recorded in generated files (default: today)
  --check             fail instead of writing when generated files are stale
  -h, --help          show this help

Outputs:
  qetlab_inventory.toml
  qetlab_inventory.csv
  qetlab_dependency_graph.dot
""",
    )
end

function parse_options(args)
    source = DEFAULT_SOURCE
    output_dir = DEFAULT_OUTPUT
    status_path = DEFAULT_STATUS
    audit_date = string(Dates.today())
    audit_date_explicit = false
    check = false
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            usage()
            exit(0)
        elseif arg == "--check"
            check = true
        elseif startswith(arg, "--source=")
            source = split(arg, "="; limit=2)[2]
        elseif arg == "--source"
            i == length(args) && error("--source requires a path")
            i += 1
            source = args[i]
        elseif startswith(arg, "--output-dir=")
            output_dir = split(arg, "="; limit=2)[2]
        elseif arg == "--output-dir"
            i == length(args) && error("--output-dir requires a path")
            i += 1
            output_dir = args[i]
        elseif startswith(arg, "--status=")
            status_path = split(arg, "="; limit=2)[2]
        elseif arg == "--status"
            i == length(args) && error("--status requires a path")
            i += 1
            status_path = args[i]
        elseif startswith(arg, "--date=")
            audit_date = split(arg, "="; limit=2)[2]
            audit_date_explicit = true
        elseif arg == "--date"
            i == length(args) && error("--date requires YYYY-MM-DD")
            i += 1
            audit_date = args[i]
            audit_date_explicit = true
        elseif startswith(arg, "-")
            error("unknown option: $arg")
        else
            push!(positional, arg)
        end
        i += 1
    end
    length(positional) <= 1 || error("at most one positional checkout path is accepted")
    isempty(positional) || (source = only(positional))
    output_dir = abspath(output_dir)
    if check && !audit_date_explicit
        inventory_path = joinpath(output_dir, "qetlab_inventory.toml")
        if isfile(inventory_path)
            existing_inventory = TOML.parsefile(inventory_path)
            recorded_date = get(existing_inventory, "audit_date", nothing)
            recorded_date isa AbstractString ||
                error("existing inventory has no string audit_date: $inventory_path")
            audit_date = String(recorded_date)
        end
    end
    try
        Date(audit_date)
    catch
        error("invalid --date value '$audit_date'; expected YYYY-MM-DD")
    end
    return (
        source=abspath(source),
        output_dir,
        status_path=abspath(status_path),
        audit_date,
        check,
    )
end

function normalize_newlines(text::AbstractString)
    return replace(replace(String(text), "\r\n" => "\n"), "\r" => "\n")
end

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

function read_source_text(path::AbstractString)
    bytes = read(path)
    isvalid(String, bytes) && return String(bytes)
    # A few historical MATLAB files contain isolated Windows-1252 punctuation.
    # Preserve valid UTF-8 runs and decode only invalid individual bytes.
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

function readchomp_cmd(argv::Vector{String})
    command = Cmd(argv)
    success(command) || error("command failed: $(join(argv, ' '))")
    return chomp(read(command, String))
end

function git_revision(source::AbstractString)
    git_dir = joinpath(source, ".git")
    ispath(git_dir) || error("QETLAB source is not a Git checkout: $source")
    revision = readchomp_cmd(["git", "-C", source, "rev-parse", "HEAD"])
    occursin(r"^[0-9a-f]{40}$", revision) ||
        error("unexpected Git revision returned for $source: $revision")
    return revision
end

function matlab_files(source::AbstractString)
    files = String[]
    for (directory, subdirectories, names) in walkdir(source)
        filter!(name -> name != ".git", subdirectories)
        for name in names
            endswith(lowercase(name), ".m") || continue
            push!(files, joinpath(directory, name))
        end
    end
    sort!(files; by=path -> lowercase(normalize_path(relpath(path, source))))
    isempty(files) && error("no MATLAB .m files found beneath $source")
    return files
end

function function_header(text::AbstractString)
    lines = split(normalize_newlines(text), '\n')
    index = findfirst(line -> occursin(r"^\s*function(?:\s|$)"i, line), lines)
    isnothing(index) && return ""
    header = strip(lines[index])
    while occursin(r"\.\.\.\s*(?:%.*)?$", header) && index < length(lines)
        header = replace(header, r"\.\.\.\s*(?:%.*)?$" => " ")
        index += 1
        header *= " " * strip(lines[index])
    end
    header = replace(header, r"\s+" => " ")
    return strip(first(split(header, "%"; limit=2)))
end

function split_arguments(arguments::AbstractString)
    stripped = strip(arguments)
    isempty(stripped) && return String[]
    result = String[]
    buffer = IOBuffer()
    nesting = 0
    for char in stripped
        if char in ('(', '[', '{')
            nesting += 1
        elseif char in (')', ']', '}')
            nesting = max(0, nesting - 1)
        end
        if char == ',' && nesting == 0
            push!(result, strip(String(take!(buffer))))
        else
            print(buffer, char)
        end
    end
    push!(result, strip(String(take!(buffer))))
    filter!(!isempty, result)
    return result
end

function parse_signature(header::AbstractString, fallback_name::AbstractString)
    isempty(header) && return (
        name=String(fallback_name),
        arguments=String[],
        outputs=String[],
        variable_outputs=false,
    )
    body = strip(replace(header, r"^function\s*"i => ""))
    equal_at = findfirst(==('='), body)
    if isnothing(equal_at)
        left = ""
        right = body
    else
        left = strip(body[firstindex(body):prevind(body, equal_at)])
        right = strip(body[nextind(body, equal_at):lastindex(body)])
    end
    matched = match(r"^([A-Za-z][A-Za-z0-9_]*)\s*(?:\((.*)\))?$", right)
    if isnothing(matched)
        return (
            name=String(fallback_name),
            arguments=String[],
            outputs=String[],
            variable_outputs=false,
        )
    end
    name = matched.captures[1]
    raw_arguments = something(matched.captures[2], "")
    arguments = split_arguments(raw_arguments)
    outputs = if isempty(left)
        String[]
    elseif startswith(left, "[") && endswith(left, "]")
        filter!(!isempty, split(strip(left[2:(end - 1)]), r"[,\s]+"))
    else
        [left]
    end
    variable_outputs = any(output -> lowercase(output) == "varargout", outputs)
    return (; name, arguments, outputs, variable_outputs)
end

function help_summary(text::AbstractString, function_name::AbstractString)
    lines = split(normalize_newlines(text), '\n')
    for line in lines
        matched = match(r"^\s*%%\s+\S+\s+(.+?)\s*$", line)
        if !isnothing(matched)
            summary = strip(matched.captures[1])
            isempty(summary) || return summary
        end
    end
    for line in lines[1:min(end, 80)]
        matched = match(r"^\s*%\s*(.+?)\s*$", line)
        isnothing(matched) && continue
        candidate = strip(matched.captures[1])
        isempty(candidate) && continue
        occursin(r"^(URL|requires?|authors?|package|last updated)\s*:"i, candidate) &&
            continue
        occursin(r"^[A-Za-z0-9_~,\[\] ]+\s*=", candidate) && continue
        lowercase(candidate) == lowercase(function_name) && continue
        return candidate
    end
    return ""
end

function documentation_url(text::AbstractString)
    for line in split(normalize_newlines(text), '\n')
        matched = match(r"^\s*%\s*URL\s*:\s*(\S+)"i, line)
        isnothing(matched) || return rstrip(matched.captures[1], ('.', ',', ';'))
    end
    return ""
end

function author_annotation(text::AbstractString)
    lines = split(normalize_newlines(text), '\n')
    for (index, line) in pairs(lines)
        matched = match(r"^\s*%\s*authors?\s*:\s*(.*?)\s*$"i, line)
        isnothing(matched) && continue
        pieces = [strip(matched.captures[1])]
        cursor = index + 1
        while cursor <= length(lines)
            continuation = match(r"^\s*%\s{2,}(.+?)\s*$", lines[cursor])
            isnothing(continuation) && break
            candidate = strip(continuation.captures[1])
            occursin(r"^(package|last updated|requires?|URL)\s*:"i, candidate) && break
            push!(pieces, candidate)
            cursor += 1
        end
        return strip(join(filter(!isempty, pieces), " "))
    end
    return ""
end

function requirement_payloads(text::AbstractString)
    lines = split(normalize_newlines(text), '\n')
    payloads = String[]
    index = 1
    while index <= length(lines)
        matched = match(r"^\s*%\s*requires?\s*:\s*(.*?)\s*$"i, lines[index])
        if isnothing(matched)
            index += 1
            continue
        end
        pieces = [strip(matched.captures[1])]
        cursor = index + 1
        while cursor <= length(lines)
            continuation = match(r"^\s*%\s+(.+?)\s*$", lines[cursor])
            isnothing(continuation) && break
            candidate = strip(continuation.captures[1])
            isempty(candidate) && break
            occursin(r"^(authors?|package|last updated|URL|references?)\s*:"i, candidate) &&
                break
            push!(pieces, candidate)
            cursor += 1
        end
        push!(payloads, join(pieces, " "))
        index = max(index + 1, cursor)
    end
    return payloads
end

function unique_sorted(values)
    return sort!(unique!(String[String(value) for value in values]); by=lowercase)
end

function declared_dependencies(text::AbstractString)
    dependencies = String[]
    payloads = requirement_payloads(text)
    for payload in payloads
        for matched in eachmatch(r"\b([A-Za-z][A-Za-z0-9_]*)\.m\b"i, payload)
            push!(dependencies, matched.captures[1])
        end
        lower = lowercase(payload)
        occursin(r"\bcvx\b", lower) && push!(dependencies, "CVX")
        occursin(r"\byalmip\b", lower) && push!(dependencies, "YALMIP")
        occursin(r"\bsedumi\b", lower) && push!(dependencies, "SeDuMi")
        occursin(r"\bmosek\b", lower) && push!(dependencies, "MOSEK")
    end
    return unique_sorted(dependencies)
end

function optional_defaults(text::AbstractString)
    names = String[]
    defaults = String[]
    lines = split(normalize_newlines(text), '\n')
    for line in lines
        matched = match(
            r"^\s*%\s+([A-Za-z][A-Za-z0-9_]*)\s+\(default\s+(.+?)\)\s*(?::.*)?$"i, line
        )
        if !isnothing(matched)
            name = uppercase(strip(matched.captures[1]))
            default = strip(matched.captures[2])
            if !(lowercase(name) in lowercase.(names))
                push!(names, name)
                push!(defaults, "$name=$default")
            end
        end
        inline = match(r"set optional (?:input )?argument defaults?\s*:\s*(.+)$"i, line)
        isnothing(inline) && continue
        for item in
            eachmatch(r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^,;%]+)", inline.captures[1])
            name = uppercase(strip(item.captures[1]))
            default = strip(item.captures[2])
            if !(lowercase(name) in lowercase.(names))
                push!(names, name)
                push!(defaults, "$name=$default")
            end
        end
    end
    return (; names, defaults)
end

function code_without_comments(text::AbstractString)
    lines = split(normalize_newlines(text), '\n')
    return join((first(split(line, "%"; limit=2)) for line in lines), "\n")
end

function called_identifiers(code::AbstractString)
    identifiers = String[]
    for matched in eachmatch(r"\b([A-Za-z][A-Za-z0-9_]*)\s*\(", code)
        offset = matched.offset
        if offset > firstindex(code)
            previous = code[prevind(code, offset)]
            previous == '.' && continue
        end
        push!(identifiers, matched.captures[1])
    end
    return unique_sorted(identifiers)
end

function toolbox_dependencies(code::AbstractString, declared::Vector{String})
    lower = lowercase(code)
    toolboxes = String[]
    (
        "CVX" in declared ||
        occursin(r"\bcvx_(?:begin|end|precision|optval)\b", lower) ||
        occursin(r"\bisa\s*\([^,\n]+,\s*['\"]cvx['\"]", lower)
    ) && push!(toolboxes, "CVX")
    occursin(r"\b(?:fmincon|fminunc|linprog|quadprog|intlinprog|lsqnonlin)\s*\(", lower) &&
        push!(toolboxes, "MATLAB Optimization Toolbox")
    occursin(r"\b(?:syms|sym|vpa|simplify|solve)\s*\(", lower) &&
        push!(toolboxes, "MATLAB Symbolic Math Toolbox")
    occursin(r"\b(?:normrnd|unifrnd|gamrnd|randsample|mvnrnd)\s*\(", lower) &&
        push!(toolboxes, "MATLAB Statistics and Machine Learning Toolbox")
    occursin(r"\bpadarray\s*\(", lower) &&
        push!(toolboxes, "MATLAB Image Processing Toolbox")
    occursin(r"\b(?:de2bi|bi2de)\s*\(", lower) &&
        push!(toolboxes, "MATLAB Communications Toolbox")
    occursin(r"\bparfor\b", lower) && push!(toolboxes, "MATLAB Parallel Computing Toolbox")
    "YALMIP" in declared && push!(toolboxes, "YALMIP")
    "SeDuMi" in declared && push!(toolboxes, "SeDuMi")
    "MOSEK" in declared && push!(toolboxes, "MOSEK")
    return unique_sorted(toolboxes)
end

function random_calls(code::AbstractString)
    calls = String[]
    for matched in
        eachmatch(r"\b(rand|randn|randi|randperm|randg|sprand|sprandn|rng)\s*\("i, code)
        push!(calls, lowercase(matched.captures[1]))
    end
    return unique_sorted(calls)
end

function snake_case(name::AbstractString)
    result = replace(String(name), r"([A-Z]+)([A-Z][a-z])" => s"\1_\2")
    result = replace(result, r"([a-z0-9])([A-Z])" => s"\1_\2")
    result = replace(result, r"[^A-Za-z0-9]+" => "_")
    return lowercase(strip(result, '_'))
end

function classification(relative_path::AbstractString, text::AbstractString)
    filename = lowercase(splitext(basename(relative_path))[1])
    startswith(filename, "example") && return "example"
    if occursin(r"(?im)^\s*%.*\b(?:deprecated|obsolete)\b", text)
        return "deprecated"
    end
    startswith(normalize_path(relative_path), "helpers/") && return "internal"
    return "public"
end

function mathematical_category(
    function_name::AbstractString, class::AbstractString, summary::AbstractString
)
    class == "internal" && return "internal_helper"
    class == "example" && return "example"
    key = lowercase(function_name * " " * summary)
    occursin(r"\brandom|rand", key) && return "random_generation"
    occursin(r"coher|incoh", key) && return "quantum_coherence"
    occursin(r"polynomial|copositive|clique", key) && return "polynomial_optimization"
    occursin(r"nonlocal|bell|game|npa|xor", key) && return "nonlocality_and_games"
    occursin(r"channel|superoperator|choi|kraus|map\b|dephas|depolar", key) &&
        return "quantum_channels_and_maps"
    occursin(
        r"separab|entang|ppt|negativ|concurr|schmidt|block.?positive|upb|witness", key
    ) && return "entanglement"
    occursin(r"partial|permut|swap|realign|tensor|subsystem|kronecker", key) &&
        return "tensor_and_subsystem_operations"
    occursin(r"state|dicke|ghz|werner|isotropic|brauer|gisin|horodecki", key) &&
        return "quantum_states"
    occursin(r"norm|fidelity|entropy|purity|projection|matrix|pauli|gell", key) &&
        return "linear_algebra_and_matrix_analysis"
    return "quantum_information_miscellaneous"
end

function proposed_signature(
    function_name::AbstractString, arguments::Vector{String}, optional_names::Vector{String}
)
    required = [
        snake_case(argument) for argument in arguments if lowercase(argument) != "varargin"
    ]
    optional = [snake_case(argument) * "=…" for argument in optional_names]
    positional = join(required, ", ")
    if isempty(optional)
        return "$(snake_case(function_name))($positional)"
    elseif isempty(positional)
        return "$(snake_case(function_name))(; $(join(optional, ", ")))"
    end
    return "$(snake_case(function_name))($positional; $(join(optional, ", ")))"
end

function source_sha256(path::AbstractString)
    return bytes2hex(sha256(read(path)))
end

function build_base_records(source::AbstractString, revision::AbstractString)
    records = Vector{Dict{String,Any}}()
    for path in matlab_files(source)
        relative_path = normalize_path(relpath(path, source))
        text = normalize_newlines(read_source_text(path))
        header = function_header(text)
        fallback_name = splitext(basename(path))[1]
        signature = parse_signature(header, fallback_name)
        fixed_output_arity = count(
            output -> lowercase(output) != "varargout", signature.outputs
        )
        summary = help_summary(text, signature.name)
        declared = declared_dependencies(text)
        optional = optional_defaults(text)
        code = code_without_comments(text)
        sparse_primitives = unique_sorted([
            lowercase(matched.captures[1]) for matched in
            eachmatch(r"\b(issparse|sparse|spalloc|speye|sprand|sprandn)\s*\("i, code)
        ])
        densifies = occursin(r"\bfull\s*\("i, code)
        random = random_calls(code)
        class = classification(relative_path, text)
        ambiguities = String[]
        any(argument -> lowercase(argument) == "varargin", signature.arguments) &&
            push!(ambiguities, "MATLAB varargin/default semantics require manual review")
        signature.variable_outputs &&
            push!(ambiguities, "output arity varies through MATLAB varargout")
        (occursin(r"\bnargin\b"i, code) || occursin(r"\bnargout\b"i, code)) &&
            push!(ambiguities, "behavior depends on MATLAB input/output arity")
        !isempty(random) && push!(
            ambiguities,
            "implicit MATLAB global RNG; Julia API must accept an explicit rng",
        )
        densifies && push!(
            ambiguities,
            "calls full(...); densification conditions and size limits require review",
        )
        toolboxes = toolbox_dependencies(code, declared)
        "CVX" in toolboxes && push!(
            ambiguities,
            "CVX solver status and failure semantics require an explicit Julia design",
        )
        isempty(documentation_url(text)) &&
            push!(ambiguities, "no upstream documentation URL annotation found")
        license_status = if lowercase(signature.name) in ("spnull", "sporth")
            "separate-bruno-luong-bsd-2-clause-license-verified"
        else
            "qetlab-bsd-2-clause-license-verified"
        end
        record = Dict{String,Any}(
            "source_path" => relative_path,
            "source_sha256" => source_sha256(path),
            "source_revision" => revision,
            "function_name" => signature.name,
            "classification" => class,
            "matlab_signature" => header,
            "matlab_arguments" => signature.arguments,
            "optional_arguments" => optional.names,
            "default_arguments" => optional.defaults,
            "output_arity" => fixed_output_arity,
            "variable_output_arity" => signature.variable_outputs,
            "output_arity_semantics" => if signature.variable_outputs
                "at_least_$(fixed_output_arity)_fixed_outputs_plus_varargout"
            else
                "exactly_$(fixed_output_arity)_outputs"
            end,
            "output_names" => signature.outputs,
            "output_semantics" => summary,
            "documentation_url" => documentation_url(text),
            "author_annotation" => author_annotation(text),
            "declared_dependencies" => declared,
            "unresolved_declared_dependencies" => String[],
            "inferred_dependencies" => String[],
            "upstream_calls" => String[],
            "external_toolboxes" => toolboxes,
            "sparse_support" => isempty(sparse_primitives) ? "not_declared" : "explicit",
            "sparse_primitives" => sparse_primitives,
            "densification" => densifies ? "explicit_full_call" : "none_detected",
            "direct_rng_calls" => random,
            "uses_global_rng_via" => String[],
            "random_behavior" => if isempty(random)
                "deterministic_no_direct_or_transitive_rng_calls_detected"
            else
                "implicit_matlab_global_rng_direct:" * join(random, ",")
            end,
            "mathematical_category" =>
                mathematical_category(signature.name, class, summary),
            "proposed_julia_name" => snake_case(signature.name),
            "proposed_julia_signature" =>
                proposed_signature(signature.name, signature.arguments, optional.names),
            "compatibility_alias" =>
                class in ("public", "deprecated") ? signature.name : "",
            "implementation_status" => "not_started",
            "test_status" => "not_started",
            "documentation_status" => "not_started",
            "benchmark_status" => "not_started",
            "license_provenance_status" => license_status,
            "review_status" => "automated_parse_pending_manual_review",
            "review_notes" => String[],
            "known_numerical_or_semantic_ambiguities" => unique_sorted(ambiguities),
            "line_count" => count(==('\n'), text) + 1,
            "_called_identifiers" => called_identifiers(code),
        )
        push!(records, record)
    end
    return records
end

const STATUS_OVERLAY_FIELDS = Set([
    "implementation_status",
    "test_status",
    "documentation_status",
    "benchmark_status",
    "license_provenance_status",
    "proposed_julia_name",
    "proposed_julia_signature",
    "compatibility_alias",
    "review_status",
    "review_notes",
])

function apply_status_overlay!(records, status_path::AbstractString)
    isfile(status_path) || error(
        "reviewed status overlay does not exist: $status_path\n" *
        "Create it with `schema_version = 1` or pass --status PATH.",
    )
    overlay = TOML.parsefile(status_path)
    top_level_allowed = Set(["schema_version", "functions"])
    unknown_top_level = setdiff(Set(keys(overlay)), top_level_allowed)
    isempty(unknown_top_level) || error(
        "unknown top-level status-overlay keys: $(join(sort!(collect(unknown_top_level)), ", "))",
    )
    get(overlay, "schema_version", nothing) == 1 ||
        error("status overlay schema_version must equal 1: $status_path")
    entries = get(overlay, "functions", Any[])
    entries isa AbstractVector ||
        error("status overlay `functions` must be an array of tables")

    by_name = Dict(record["function_name"] => index for (index, record) in pairs(records))
    by_path = Dict(record["source_path"] => index for (index, record) in pairs(records))
    applied = Set{Int}()
    identity_fields = Set(["function_name", "source_path"])
    for (entry_number, entry) in pairs(entries)
        entry isa AbstractDict || error("status overlay entry $entry_number is not a table")
        unknown = setdiff(Set(keys(entry)), union(identity_fields, STATUS_OVERLAY_FIELDS))
        isempty(unknown) || error(
            "unknown keys in status overlay entry $entry_number: " *
            join(sort!(collect(unknown)), ", "),
        )
        has_name = haskey(entry, "function_name")
        has_path = haskey(entry, "source_path")
        (has_name || has_path) ||
            error("status overlay entry $entry_number needs function_name or source_path")
        index_from_name = if has_name
            name = entry["function_name"]
            name isa AbstractString || error(
                "function_name in status overlay entry $entry_number must be a string",
            )
            haskey(by_name, name) || error(
                "unknown function_name in status overlay entry $entry_number: $name"
            )
            by_name[name]
        else
            nothing
        end
        index_from_path = if has_path
            path = normalize_path(entry["source_path"])
            haskey(by_path, path) || error(
                "unknown source_path in status overlay entry $entry_number: $path"
            )
            by_path[path]
        else
            nothing
        end
        if !isnothing(index_from_name) &&
            !isnothing(index_from_path) &&
            index_from_name != index_from_path
            error(
                "function_name and source_path identify different functions in overlay entry $entry_number",
            )
        end
        index = something(index_from_name, index_from_path)
        index in applied &&
            error("duplicate status overlay entry for $(records[index]["function_name"])")
        push!(applied, index)
        for field in STATUS_OVERLAY_FIELDS
            haskey(entry, field) || continue
            value = entry[field]
            if field == "review_notes"
                value isa AbstractVector && all(note -> note isa AbstractString, value) ||
                    error(
                        "review_notes for $(records[index]["function_name"]) must be an array of strings",
                    )
                records[index][field] = unique_sorted(value)
            else
                value isa AbstractString && !isempty(strip(value)) || error(
                    "$field for $(records[index]["function_name"]) must be a non-empty string",
                )
                records[index][field] = String(value)
            end
        end
    end
    return (entry_count=length(entries), sha256=source_sha256(status_path))
end

function add_dependency_analysis!(records)
    canonical = Dict{String,String}()
    known_external = Set(["cvx", "yalmip", "sedumi", "mosek"])
    for record in records
        canonical[lowercase(record["function_name"])] = record["function_name"]
    end
    for record in records
        self = lowercase(record["function_name"])
        inferred = String[]
        for called in record["_called_identifiers"]
            key = lowercase(called)
            haskey(canonical, key) || continue
            key == self && continue
            push!(inferred, canonical[key])
        end
        declared_upstream = String[]
        unresolved_declared = String[]
        for dependency in record["declared_dependencies"]
            key = lowercase(dependency)
            if !haskey(canonical, key)
                key in known_external || push!(unresolved_declared, dependency)
                continue
            end
            key == self && continue
            push!(declared_upstream, canonical[key])
        end
        unresolved_declared = unique_sorted(unresolved_declared)
        if !isempty(unresolved_declared)
            push!(
                record["known_numerical_or_semantic_ambiguities"],
                "declared dependency absent at pinned revision: " *
                join(unresolved_declared, ", "),
            )
            record["known_numerical_or_semantic_ambiguities"] = unique_sorted(
                record["known_numerical_or_semantic_ambiguities"]
            )
        end
        record["unresolved_declared_dependencies"] = unresolved_declared
        record["inferred_dependencies"] = unique_sorted(inferred)
        record["upstream_calls"] = unique_sorted(vcat(inferred, declared_upstream))
        delete!(record, "_called_identifiers")
    end

    # Propagate implicit global-RNG use through the upstream call graph. The
    # inventory must not describe wrappers such as RandomPOVM as deterministic.
    rng_functions = Set(
        record["function_name"] for
        record in records if !isempty(record["direct_rng_calls"])
    )
    changed = true
    while changed
        changed = false
        for record in records
            record["function_name"] in rng_functions && continue
            any(dependency -> dependency in rng_functions, record["upstream_calls"]) ||
                continue
            push!(rng_functions, record["function_name"])
            changed = true
        end
    end
    for record in records
        gateways = unique_sorted(
            dependency for
            dependency in record["upstream_calls"] if dependency in rng_functions
        )
        record["uses_global_rng_via"] = gateways
        if isempty(record["direct_rng_calls"]) && !isempty(gateways)
            record["random_behavior"] =
                "implicit_matlab_global_rng_transitive_via:" * join(gateways, ",")
            push!(
                record["known_numerical_or_semantic_ambiguities"],
                "transitively uses implicit MATLAB global RNG via " * join(gateways, ", "),
            )
            record["known_numerical_or_semantic_ambiguities"] = unique_sorted(
                record["known_numerical_or_semantic_ambiguities"]
            )
        end
    end
    return records
end

function dependency_graph(records)
    graph = Dict{String,Vector{String}}()
    for record in records
        graph[record["function_name"]] = copy(record["upstream_calls"])
    end
    return graph
end

function strongly_connected_components(graph::Dict{String,Vector{String}})
    visited = Set{String}()
    order = String[]
    function visit(node)
        node in visited && return nothing
        push!(visited, node)
        for neighbour in get(graph, node, String[])
            visit(neighbour)
        end
        return push!(order, node)
    end
    for node in sort!(collect(keys(graph)); by=lowercase)
        visit(node)
    end

    reverse_graph = Dict(node => String[] for node in keys(graph))
    for (node, neighbours) in graph
        for neighbour in neighbours
            haskey(reverse_graph, neighbour) || continue
            push!(reverse_graph[neighbour], node)
        end
    end
    empty!(visited)
    components = Vector{Vector{String}}()
    function reverse_visit(node, component)
        node in visited && return nothing
        push!(visited, node)
        push!(component, node)
        for neighbour in reverse_graph[node]
            reverse_visit(neighbour, component)
        end
    end
    for node in reverse(order)
        node in visited && continue
        component = String[]
        reverse_visit(node, component)
        sort!(component; by=lowercase)
        if length(component) > 1 || node in get(graph, node, String[])
            push!(components, component)
        end
    end
    sort!(components; by=component -> lowercase(join(component, ",")))
    return components
end

function toml_escape(value::AbstractString)
    return replace(
        String(value),
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\b" => "\\b",
        "\t" => "\\t",
        "\n" => "\\n",
        "\f" => "\\f",
        "\r" => "\\r",
    )
end

toml_string(value::AbstractString) = "\"" * toml_escape(value) * "\""
toml_value(value::AbstractString) = toml_string(value)
toml_value(value::Bool) = value ? "true" : "false"
toml_value(value::Integer) = string(value)
toml_value(values::AbstractVector) = "[" * join(toml_value.(values), ", ") * "]"

const TOML_FIELD_ORDER = [
    "source_path",
    "source_sha256",
    "source_revision",
    "function_name",
    "classification",
    "matlab_signature",
    "matlab_arguments",
    "optional_arguments",
    "default_arguments",
    "output_arity",
    "variable_output_arity",
    "output_arity_semantics",
    "output_names",
    "output_semantics",
    "documentation_url",
    "author_annotation",
    "declared_dependencies",
    "unresolved_declared_dependencies",
    "inferred_dependencies",
    "upstream_calls",
    "external_toolboxes",
    "sparse_support",
    "sparse_primitives",
    "densification",
    "direct_rng_calls",
    "uses_global_rng_via",
    "random_behavior",
    "mathematical_category",
    "proposed_julia_name",
    "proposed_julia_signature",
    "compatibility_alias",
    "implementation_status",
    "test_status",
    "documentation_status",
    "benchmark_status",
    "license_provenance_status",
    "review_status",
    "review_notes",
    "known_numerical_or_semantic_ambiguities",
    "line_count",
]

function toml_inventory(
    records,
    revision::AbstractString,
    audit_date::AbstractString,
    graph,
    cycles,
    status_path::AbstractString,
    status_overlay,
)
    io = IOBuffer()
    public_count = count(record -> record["classification"] == "public", records)
    internal_count = count(record -> record["classification"] == "internal", records)
    deprecated_count = count(record -> record["classification"] == "deprecated", records)
    example_count = count(record -> record["classification"] == "example", records)
    edge_count = sum(length, values(graph))
    pending_review_count = count(
        record ->
            get(record, "review_status", "") == "automated_parse_pending_manual_review",
        records,
    )
    inventory_review_status = if iszero(pending_review_count)
        "source_review_complete"
    else
        "automated_parse_pending_manual_review"
    end
    println(io, "# Generated by scripts/build_upstream_inventory.jl; do not edit by hand.")
    println(io, "schema_version = $GENERATOR_VERSION")
    println(io, "source_name = \"QETLAB\"")
    println(io, "review_status = ", toml_string(inventory_review_status))
    println(
        io,
        "status_overlay = ",
        toml_string(normalize_path(relpath(status_path, REPOSITORY_ROOT))),
    )
    println(io, "status_overlay_sha256 = ", toml_string(status_overlay.sha256))
    println(io, "status_overlay_entry_count = ", status_overlay.entry_count)
    println(io, "source_reviewed_count = ", length(records) - pending_review_count)
    println(io, "pending_review_count = ", pending_review_count)
    println(io, "source_revision = ", toml_string(revision))
    println(io, "audit_date = ", toml_string(audit_date))
    println(io, "matlab_file_count = ", length(records))
    println(io, "public_function_count = ", public_count)
    println(io, "internal_helper_count = ", internal_count)
    println(io, "deprecated_function_count = ", deprecated_count)
    println(io, "example_count = ", example_count)
    println(io, "dependency_edge_count = ", edge_count)
    println(io, "dependency_cycle_count = ", length(cycles))
    println(io, "dependency_cycles = ", toml_value(cycles))
    for record in records
        println(io)
        println(io, "[[functions]]")
        for field in TOML_FIELD_ORDER
            println(io, field, " = ", toml_value(record[field]))
        end
    end
    return String(take!(io))
end

function csv_escape(value)
    string_value = if value isa AbstractVector
        join(string.(value), ";")
    else
        string(value)
    end
    return "\"" * replace(string_value, "\"" => "\"\"") * "\""
end

function csv_inventory(records)
    fields = TOML_FIELD_ORDER
    io = IOBuffer()
    println(io, join(csv_escape.(fields), ","))
    for record in records
        println(io, join((csv_escape(record[field]) for field in fields), ","))
    end
    return String(take!(io))
end

function dot_escape(value::AbstractString)
    return replace(String(value), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
end

function dot_dependency_graph(records, revision::AbstractString, cycles)
    io = IOBuffer()
    cycle_nodes = Set(Iterators.flatten(cycles))
    println(io, "// Generated by scripts/build_upstream_inventory.jl; do not edit by hand.")
    println(io, "// QETLAB revision: $revision")
    if isempty(cycles)
        println(io, "// Static dependency analysis found no cycles.")
    else
        for cycle in cycles
            println(io, "// Dependency cycle: ", join(cycle, " -> "))
        end
    end
    println(io, "digraph qetlab_dependencies {")
    println(io, "  graph [rankdir=LR, overlap=false, splines=true];")
    println(io, "  node [fontname=\"Helvetica\", fontsize=9, style=\"filled\"];")
    println(io, "  edge [color=\"#777777\", arrowsize=0.6];")
    by_name = Dict(record["function_name"] => record for record in records)
    for name in sort!(collect(keys(by_name)); by=lowercase)
        record = by_name[name]
        class = record["classification"]
        shape = class == "internal" ? "ellipse" : "box"
        fill = if name in cycle_nodes
            "#ffcccc"
        elseif class == "internal"
            "#eeeeee"
        else
            "#dceeff"
        end
        label = dot_escape(name * "\n" * record["source_path"])
        println(
            io,
            "  \"",
            dot_escape(name),
            "\" [label=\"",
            label,
            "\", shape=",
            shape,
            ", fillcolor=\"",
            fill,
            "\"];",
        )
    end
    for name in sort!(collect(keys(by_name)); by=lowercase)
        for dependency in sort!(copy(by_name[name]["upstream_calls"]); by=lowercase)
            println(io, "  \"", dot_escape(name), "\" -> \"", dot_escape(dependency), "\";")
        end
    end
    println(io, "}")
    return String(take!(io))
end

function write_or_check(path::AbstractString, content::AbstractString, check::Bool)
    if check
        if !isfile(path)
            println(stderr, "missing generated artifact: $path")
            return false
        elseif read(path, String) != content
            println(stderr, "stale generated artifact: $path")
            return false
        end
        return true
    end
    open(path, "w") do io
        return write(io, content)
    end
    println("wrote ", relpath(path, REPOSITORY_ROOT))
    return true
end

function main(args)
    options = parse_options(args)
    isdir(options.source) || error("QETLAB checkout does not exist: $(options.source)")
    revision = git_revision(options.source)
    records = build_base_records(options.source, revision)
    add_dependency_analysis!(records)
    status_overlay = apply_status_overlay!(records, options.status_path)
    graph = dependency_graph(records)
    cycles = strongly_connected_components(graph)
    mkpath(options.output_dir)
    artifacts = [
        (
            joinpath(options.output_dir, "qetlab_inventory.toml"),
            toml_inventory(
                records,
                revision,
                options.audit_date,
                graph,
                cycles,
                options.status_path,
                status_overlay,
            ),
        ),
        (joinpath(options.output_dir, "qetlab_inventory.csv"), csv_inventory(records)),
        (
            joinpath(options.output_dir, "qetlab_dependency_graph.dot"),
            dot_dependency_graph(records, revision, cycles),
        ),
    ]
    success = all(
        write_or_check(path, content, options.check) for (path, content) in artifacts
    )
    println(
        options.check ? "checked" : "inventoried",
        " ",
        length(records),
        " MATLAB files at ",
        revision,
        "; ",
        sum(length, values(graph)),
        " dependency edges; ",
        length(cycles),
        " cycles",
    )
    return success || exit(1)
end

main(ARGS)
