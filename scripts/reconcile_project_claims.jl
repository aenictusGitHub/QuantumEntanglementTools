#!/usr/bin/env julia

# Deterministic reconciliation of repository-owned quantitative claims.
#
# The generated snapshot is deliberately free of wall-clock timestamps and the
# current Git commit: either would make a tracked snapshot self-invalidating.
# Source hashes, runtime exports, ledgers, tests, and maintained prose are the
# inputs that determine the output.

using SHA
using TOML
using Test

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)
using QuantumEntanglementTools: QuantumEntanglementTools
include(joinpath(@__DIR__, "qetlab_completion_common.jl"))
using .QETLABCompletion

const DEFAULT_SNAPSHOT = joinpath(
    REPOSITORY_ROOT, "artifacts", "convergence", "current_snapshot.toml"
)
const INVENTORY_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_inventory.toml")
const STATUS_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_status.toml")
const PROVENANCE_PATH = joinpath(REPOSITORY_ROOT, "PROVENANCE.toml")
const PROJECT_PATH = joinpath(REPOSITORY_ROOT, "Project.toml")
const TEST_RUNNER = joinpath(REPOSITORY_ROOT, "test", "runtests.jl")
const EXTENSION_TEST_RUNNER = joinpath(
    REPOSITORY_ROOT, "test", "extensions", "entanglement_detection", "runtests.jl"
)
const EXTENSION_PROJECT = dirname(EXTENSION_TEST_RUNNER)
const TUTORIAL_RUNNER = joinpath(REPOSITORY_ROOT, "tutorials", "runtests.jl")
const BENCHMARK_RUNNER = joinpath(REPOSITORY_ROOT, "benchmark", "benchmarks.jl")
const CURRENT_CLAIMS_BEGIN = "<!-- qetlab-current-claims: begin -->"
const CURRENT_CLAIMS_END = "<!-- qetlab-current-claims: end -->"

const PUBLIC_MODULES = Dict(
    "QuantumEntanglementTools" => QuantumEntanglementTools,
    "QuantumEntanglementTools.MATLABCompat" => QuantumEntanglementTools.MATLABCompat,
)
const STATUS_FIELDS = (
    "implementation_status",
    "test_status",
    "documentation_status",
    "benchmark_status",
    "license_provenance_status",
    "review_status",
)
const MAINTAINED_CLAIM_FILES = (
    "AGENTS.md",
    "README.md",
    "CHANGELOG.md",
    "CONTRIBUTING.md",
    "CITATION.cff",
    "CITATION.bib",
    "docs/BENCHMARK_REPORT.md",
    "docs/BUILD_ENVIRONMENT.md",
    "docs/CONVERGENCE_AUDIT.md",
    "docs/FULL_QETLAB_BASELINE.md",
    "docs/INVENTORY_REVIEW.md",
    "docs/JULIA_CORE_AUDIT.md",
    "docs/PORTING_STATUS.md",
    "docs/QETLAB_COMPLETION_PLAN.md",
    "docs/RELEASE_CHECKLIST.md",
    "docs/ROADMAP.md",
    "docs/SESSION_HANDOFF.md",
    "docs/VALIDATION_REPORT.md",
    "docs/src/api/index.md",
    "docs/src/entanglement_backends.md",
    "docs/src/getting_started.md",
    "docs/src/index.md",
    "docs/src/migration_from_qetlab.md",
)
const RANDOM_CALLS = Set([
    "rand",
    "rand!",
    "randn",
    "randn!",
    "randexp",
    "randexp!",
    "randperm",
    "shuffle",
    "shuffle!",
    "bitrand",
    "randsubseq",
    "randcycle",
    "randstring",
])
const DENSIFICATION_CALLS = Set(["Matrix", "Array", "collect"])
const MARKER_PATTERN = r"(?i)\b(TODO|FIXME|XXX|not[ _-]?implemented|placeholder)\b"

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/reconcile_project_claims.jl [options]

Compute the QETLAB inventory, public-API, evidence, and source-audit snapshot.

Options:
  --check          Fail unless the tracked snapshot is byte-for-byte current
  --skip-tests     Reuse the tracked dynamic-test section (for release checks
                   that already have a separate full test gate)
  --output PATH    Write/check PATH instead of the default tracked snapshot
  -h, --help       Show this help

Without --check the script runs the core and optional-extension suites, verifies
maintained numeric prose claims, and writes a deterministic TOML snapshot.
""",
    )
end

function parse_options(args)
    check = false
    skip_tests = false
    output = DEFAULT_SNAPSHOT
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument in ("-h", "--help")
            usage()
            exit(0)
        elseif argument == "--check"
            check = true
        elseif argument == "--skip-tests"
            skip_tests = true
        elseif startswith(argument, "--output=")
            output = abspath(split(argument, "="; limit=2)[2])
        elseif argument == "--output"
            index == length(args) && error("--output requires a path")
            index += 1
            output = abspath(args[index])
        else
            error("unknown option: $argument")
        end
        index += 1
    end
    return (; check, skip_tests, output)
end

function public_names(mod::Module)
    return sort!(
        collect(
            String(name) for
            name in names(mod; all=false, imported=false) if name != nameof(mod)
        ),
    )
end

function count_values(entries, field)
    counts = Dict{String,Int}()
    for entry in entries
        value = string(get(entry, field, "missing"))
        counts[value] = get(counts, value, 0) + 1
    end
    return Dict{String,Any}(key => counts[key] for key in sort!(collect(keys(counts))))
end

function implementation_bucket(status::AbstractString)
    status in (
        "compatibility_alias",
        "implemented",
        "superseded_with_documented_mapping",
        "verified",
        "verified_with_documented_upstream_correction",
    ) && return "implemented"
    startswith(status, "partially_implemented") && return "partial"
    status == "deferred" && return "deferred"
    status == "blocked_with_explicit_reason" && return "blocked"
    return "other"
end

function internal_implementation_bucket(status::AbstractString)
    startswith(status, "internal_helper_superseded") && return "superseded"
    startswith(status, "intentionally_not_ported") && return "intentionally_unexposed"
    startswith(status, "partially_superseded") && return "partial"
    startswith(status, "deferred_with_") && return "deferred"
    return "other"
end

function count_buckets(entries, classifier; expected=String[])
    counts = Dict{String,Int}(bucket => 0 for bucket in expected)
    for entry in entries
        bucket = classifier(string(get(entry, "implementation_status", "missing")))
        counts[bucket] = get(counts, bucket, 0) + 1
    end
    return Dict{String,Any}(key => counts[key] for key in sort!(collect(keys(counts))))
end

function test_counts(testset::Test.AbstractTestSet)
    counts = Test.get_test_counts(testset)
    return Dict{String,Any}(
        "passed" => counts.passes + counts.cumulative_passes,
        "failed" => counts.fails + counts.cumulative_fails,
        "errored" => counts.errors + counts.cumulative_errors,
        "broken" => counts.broken + counts.cumulative_broken,
        "assertions" =>
            counts.passes +
            counts.cumulative_passes +
            counts.fails +
            counts.cumulative_fails +
            counts.errors +
            counts.cumulative_errors +
            counts.broken +
            counts.cumulative_broken,
    )
end

function find_testset(testset::Test.AbstractTestSet, description::AbstractString)
    getfield(testset, :description) == description && return testset
    for result in getfield(testset, :results)
        result isa Test.AbstractTestSet || continue
        found = find_testset(result, description)
        isnothing(found) || return found
    end
    return nothing
end

function direct_test_suites(testset::Test.AbstractTestSet)
    suites = Dict{String,Any}[]
    for result in getfield(testset, :results)
        result isa Test.AbstractTestSet || continue
        counts = test_counts(result)
        push!(
            suites,
            Dict{String,Any}(
                "name" => string(getfield(result, :description)),
                "passed" => counts["passed"],
                "failed" => counts["failed"],
                "errored" => counts["errored"],
                "broken" => counts["broken"],
                "assertions" => counts["assertions"],
            ),
        )
    end
    return suites
end

function run_core_tests()
    result = @testset "claim reconciliation package corpus" begin
        Base.include(Main, TEST_RUNNER)
    end
    counts = test_counts(result)
    counts["failed"] == 0 || error("core test measurement observed failures")
    counts["errored"] == 0 || error("core test measurement observed errors")

    tutorial_set = find_testset(result, "Executable tutorials")
    isnothing(tutorial_set) &&
        error("could not locate the Executable tutorials runtime testset")
    tutorial_counts = test_counts(tutorial_set)
    core_assertions = counts["assertions"] - tutorial_counts["assertions"]

    return Dict{String,Any}(
        "measurement" => "runtime Test.get_test_counts",
        "integrated_assertions" => counts["assertions"],
        "integrated_passed" => counts["passed"],
        "integrated_broken" => counts["broken"],
        "core_assertions" => core_assertions,
        "tutorial_assertions" => tutorial_counts["assertions"],
        "suites" => direct_test_suites(result),
    )
end

function extension_probe_program()
    runner = repr(EXTENSION_TEST_RUNNER)
    return """
    using Test
    result = @testset "claim reconciliation extension" begin
        Base.include(Main, $runner)
    end
    counts = Test.get_test_counts(result)
    passed = counts.passes + counts.cumulative_passes
    failed = counts.fails + counts.cumulative_fails
    errored = counts.errors + counts.cumulative_errors
    broken = counts.broken + counts.cumulative_broken
    println(
        "__QET_EXTENSION_COUNTS__",
        passed,
        ",",
        failed,
        ",",
        errored,
        ",",
        broken,
        ",",
        Base.pkgversion(Main.EntanglementDetection),
    )
    """
end

function run_extension_tests()
    executable = String.(Base.julia_cmd().exec)
    append!(
        executable,
        [
            "--startup-file=no",
            "--history-file=no",
            "--project=$(EXTENSION_PROJECT)",
            "-e",
            extension_probe_program(),
        ],
    )
    command = Cmd(Cmd(executable); dir=REPOSITORY_ROOT)
    output = read(pipeline(command; stderr=stderr), String)
    matched = match(
        r"__QET_EXTENSION_COUNTS__(\d+),(\d+),(\d+),(\d+),([0-9A-Za-z.+-]+)", output
    )
    isnothing(matched) && error("extension test probe did not emit its count sentinel")
    passed, failed, errored, broken = parse.(Int, matched.captures[1:4])
    failed == 0 || error("extension test measurement observed $failed failures")
    errored == 0 || error("extension test measurement observed $errored errors")
    return Dict{String,Any}(
        "measurement" => "runtime Test.get_test_counts in isolated environment",
        "assertions" => passed + failed + errored + broken,
        "passed" => passed,
        "failed" => failed,
        "errored" => errored,
        "broken" => broken,
        "backend_version" => matched.captures[5],
    )
end

function run_dynamic_tests()
    return Dict{String,Any}(
        "package" => run_core_tests(),
        "entanglement_detection_extension" => run_extension_tests(),
    )
end

function count_test_macros()
    result = Dict{String,Int}()
    roots = (joinpath(REPOSITORY_ROOT, "test"), joinpath(REPOSITORY_ROOT, "tutorials"))
    for root in roots
        for (directory, _, files) in walkdir(root)
            for file in files
                endswith(file, ".jl") || continue
                path = joinpath(directory, file)
                relative = relpath(path, REPOSITORY_ROOT)
                count = length(
                    collect(
                        eachmatch(
                            r"@test(?:set|_throws|_broken|_skip)?\b", read(path, String)
                        ),
                    ),
                )
                result[relative] = count
            end
        end
    end
    return Dict{String,Any}(
        "measurement" => "static macro occurrences; not an assertion count",
        "total" => sum(values(result)),
        "by_file" =>
            Dict{String,Any}(key => result[key] for key in sort!(collect(keys(result)))),
    )
end

function source_paths()
    paths = String[]
    for root_name in ("src", "ext")
        root = joinpath(REPOSITORY_ROOT, root_name)
        for (directory, _, files) in walkdir(root)
            for file in files
                endswith(file, ".jl") || continue
                push!(paths, joinpath(directory, file))
            end
        end
    end
    return sort!(paths)
end

function call_name(expression)
    expression isa Symbol && return String(expression)
    if expression isa Expr && expression.head == :. && !isempty(expression.args)
        tail = expression.args[end]
        tail isa QuoteNode && tail.value isa Symbol && return String(tail.value)
        tail isa Symbol && return String(tail)
    end
    return ""
end

function contains_rng_identifier(value)
    value isa Symbol && return occursin("rng", lowercase(String(value)))
    value isa QuoteNode && return contains_rng_identifier(value.value)
    value isa Expr || return false
    return any(contains_rng_identifier, value.args)
end

function definition_signature(expression)
    expression isa Expr || return nothing
    if expression.head == :function
        return expression.args[1]
    elseif expression.head == :(=)
        left = expression.args[1]
        left isa Expr && left.head in (:call, :where) && return left
    end
    return nothing
end

function signature_name(signature)
    signature isa Symbol && return String(signature)
    signature isa Expr || return ""
    signature.head == :where && return signature_name(signature.args[1])
    signature.head == :(::) && return signature_name(signature.args[1])
    signature.head == :call || return ""
    callee = signature.args[1]
    callee isa Symbol && return String(callee)
    if callee isa Expr && callee.head == :.
        return call_name(callee)
    end
    return ""
end

function normalized_signature(signature)
    return replace(sprint(show, signature), r"\s+" => " ")
end

function nearby_guard_marker(lines, line)
    start = max(firstindex(lines), line - 8)
    stop = min(lastindex(lines), line + 8)
    context = join(lines[start:stop], "\n")
    return occursin(r"(?i)allow_densify|densif|max_dense|dense budget", context)
end

function scan_expression!(
    expression,
    state,
    relative_path,
    lines;
    line::Int=1,
    exported_names::Set{String}=Set{String}(),
)
    expression isa LineNumberNode && return Int(expression.line)
    expression isa Expr || return line
    current_line = line

    signature = definition_signature(expression)
    if !isnothing(signature)
        name = signature_name(signature)
        if name in exported_names
            key = normalized_signature(signature)
            push!(
                get!(state.definitions, key, Dict{String,Any}[]),
                Dict{String,Any}(
                    "name" => name, "path" => relative_path, "line" => current_line
                ),
            )
        end
    end

    if expression.head == :call && !isempty(expression.args)
        name = call_name(expression.args[1])
        if name in RANDOM_CALLS
            arguments = expression.args[2:end]
            explicit_rng = !isempty(arguments) && contains_rng_identifier(first(arguments))
            if !explicit_rng
                push!(
                    state.global_rng_candidates,
                    Dict{String,Any}(
                        "call" => name, "path" => relative_path, "line" => current_line
                    ),
                )
            end
        elseif name == "seed!"
            push!(
                state.seed_calls,
                Dict{String,Any}("path" => relative_path, "line" => current_line),
            )
        end
        if name in DENSIFICATION_CALLS
            guarded = nearby_guard_marker(lines, current_line)
            push!(
                state.densification_sites,
                Dict{String,Any}(
                    "call" => name,
                    "path" => relative_path,
                    "line" => current_line,
                    "nearby_guard_marker" => guarded,
                ),
            )
        end
    end

    for argument in expression.args
        current_line = scan_expression!(
            argument,
            state,
            relative_path,
            lines;
            line=current_line,
            exported_names=exported_names,
        )
    end
    return current_line
end

function global_refs!(references::Set{Tuple{Module,Symbol}}, value)
    if value isa GlobalRef
        push!(references, (value.mod, value.name))
    elseif value isa Expr
        for argument in value.args
            global_refs!(references, argument)
        end
    elseif value isa Core.CodeInfo
        for statement in value.code
            global_refs!(references, statement)
        end
    elseif value isa AbstractArray
        for item in value
            global_refs!(references, item)
        end
    end
    return references
end

function function_global_refs!(
    references::Set{Tuple{Module,Symbol}}, function_value, seen::Base.IdSet{Any}
)
    function_value isa Function || return references
    function_value in seen && return references
    push!(seen, function_value)
    before = copy(references)
    for method in methods(function_value)
        code = try
            Base.uncompressed_ast(method)
        catch
            nothing
        end
        isnothing(code) || global_refs!(references, code)
    end
    for (owner, name) in setdiff(references, before)
        (
            startswith(String(name), "#") || owner === QuantumEntanglementTools.MATLABCompat
        ) || continue
        isdefined(owner, name) || continue
        nested = getfield(owner, name)
        nested isa Function || continue
        function_global_refs!(references, nested, seen)
    end
    return references
end

function resolve_binding(path::AbstractString)
    parts = split(path, '.')
    isempty(parts) && return nothing
    first(parts) == "QuantumEntanglementTools" || return nothing
    value = QuantumEntanglementTools
    for part in parts[2:end]
        value isa Module || return nothing
        symbol = Symbol(part)
        isdefined(value, symbol) || return nothing
        value = getfield(value, symbol)
    end
    return value
end

function wrapper_bypass_candidates(provenance_entries)
    candidates = Dict{String,Any}[]
    documented_independent = Dict{String,Any}[]
    compat_module = QuantumEntanglementTools.MATLABCompat
    for entry in provenance_entries
        get(entry, "public_module", "") == "QuantumEntanglementTools.MATLABCompat" ||
            continue
        delegates = get(entry, "delegates_to", Any[])
        delegates isa AbstractVector && !isempty(delegates) || continue
        name = string(get(entry, "julia_name", ""))
        symbol = Symbol(name)
        isdefined(compat_module, symbol) || continue
        wrapper = getfield(compat_module, symbol)

        resolved = filter(value -> !isnothing(value), resolve_binding.(String.(delegates)))
        any(delegate -> wrapper === delegate, resolved) && continue

        references = Set{Tuple{Module,Symbol}}()
        if wrapper isa Function
            function_global_refs!(references, wrapper, Base.IdSet{Any}())
        end
        expected = Set{Tuple{Module,Symbol}}()
        for delegate_path in String.(delegates)
            parts = split(delegate_path, '.')
            length(parts) >= 2 || continue
            owner = if length(parts) == 2
                QuantumEntanglementTools
            else
                resolve_binding(join(parts[1:(end - 1)], '.'))
            end
            owner isa Module || continue
            delegate_symbol = Symbol(last(parts))
            push!(expected, (owner, delegate_symbol))
            # `using ..QuantumEntanglementTools: name` is lowered as a
            # GlobalRef owned by MATLABCompat even though the binding resolves
            # to the reviewed native function.
            push!(expected, (compat_module, delegate_symbol))
        end
        isempty(intersect(references, expected)) || continue
        finding = Dict{String,Any}(
            "wrapper" => name,
            "declared_delegates" => sort!(String.(delegates)),
            "reason" => "no declared delegate appears in transitive uncompressed method AST",
        )
        implementation_kind = string(get(entry, "implementation_kind", ""))
        if occursin(
            "raw-parameter-compatibility-wrapper-over-native-formula", implementation_kind
        )
            finding["implementation_kind"] = implementation_kind
            push!(documented_independent, finding)
        else
            push!(candidates, finding)
        end
    end
    return (; candidates, documented_independent)
end

function source_audit(provenance_entries)
    native_names = Set(public_names(QuantumEntanglementTools))
    compat_names = Set(public_names(QuantumEntanglementTools.MATLABCompat))
    state = (
        definitions=Dict{String,Vector{Dict{String,Any}}}(),
        global_rng_candidates=Dict{String,Any}[],
        seed_calls=Dict{String,Any}[],
        densification_sites=Dict{String,Any}[],
    )
    markers = Dict{String,Any}[]

    for path in source_paths()
        relative = relpath(path, REPOSITORY_ROOT)
        content = read(path, String)
        lines = split(content, '\n'; keepempty=true)
        exported_names =
            startswith(relative, joinpath("src", "compat")) ? compat_names : native_names
        expression = Meta.parseall(content)
        scan_expression!(expression, state, relative, lines; exported_names=exported_names)
        for (line_number, text) in pairs(lines)
            occursin(MARKER_PATTERN, text) || continue
            push!(
                markers,
                Dict{String,Any}(
                    "path" => relative, "line" => line_number, "text" => strip(text)
                ),
            )
        end
    end

    duplicate_definitions = Dict{String,Any}[]
    for (signature, sites) in state.definitions
        length(sites) > 1 || continue
        push!(
            duplicate_definitions,
            Dict{String,Any}("signature" => signature, "sites" => sites),
        )
    end
    sort!(duplicate_definitions; by=entry -> entry["signature"])
    sort!(state.global_rng_candidates; by=entry -> (entry["path"], entry["line"]))
    sort!(state.seed_calls; by=entry -> (entry["path"], entry["line"]))
    sort!(state.densification_sites; by=entry -> (entry["path"], entry["line"]))
    sort!(markers; by=entry -> (entry["path"], entry["line"]))

    unguarded_densification = count(
        site -> !site["nearby_guard_marker"], state.densification_sites
    )
    bypass = wrapper_bypass_candidates(provenance_entries)
    sort!(bypass.candidates; by=entry -> entry["wrapper"])
    sort!(bypass.documented_independent; by=entry -> entry["wrapper"])

    return Dict{String,Any}(
        "unimplemented_marker_count" => length(markers),
        "unimplemented_markers" => markers,
        "duplicate_public_definition_candidate_count" => length(duplicate_definitions),
        "duplicate_public_definition_candidates" => duplicate_definitions,
        "compatibility_wrapper_bypass_candidate_count" => length(bypass.candidates),
        "compatibility_wrapper_bypass_candidates" => bypass.candidates,
        "documented_independent_compatibility_implementation_count" =>
            length(bypass.documented_independent),
        "documented_independent_compatibility_implementations" =>
            bypass.documented_independent,
        "global_rng_candidate_count" => length(state.global_rng_candidates),
        "global_rng_candidates" => state.global_rng_candidates,
        "random_seed_call_count" => length(state.seed_calls),
        "random_seed_calls" => state.seed_calls,
        "densification_call_count" => length(state.densification_sites),
        "densification_without_nearby_guard_marker_count" => unguarded_densification,
        "densification_sites" => state.densification_sites,
        "static_scan_limitations" => "Candidates require review: lexical/AST scans do not prove semantic duplication, RNG mutation, or undocumented densification.",
    )
end

function tutorial_summary()
    includes = String[]
    content = read(TUTORIAL_RUNNER, String)
    for matched in eachmatch(r"include\(\"([^\"]+\.jl)\"\)", content)
        push!(includes, matched.captures[1])
    end
    return Dict{String,Any}(
        "count" => length(includes),
        "scripts" => sort!(includes),
        "runner" => relpath(TUTORIAL_RUNNER, REPOSITORY_ROOT),
    )
end

function benchmark_summary()
    content = read(BENCHMARK_RUNNER, String)
    names = String[
        matched.captures[1] for
        matched in eachmatch(r"suite\[\"([^\"]+)\"\]\s*=\s*@benchmarkable", content)
    ]
    return Dict{String,Any}(
        "count" => length(names),
        "cases" => sort!(names),
        "measurement" => "static BenchmarkGroup assignment count; execution is a separate gate",
    )
end

function completion_summary(model, audit)
    terminal_helpers = count(row -> !row["required_for_completion"], model.helper_rows)
    required_helpers = count(row -> row["required_for_completion"], model.helper_rows)
    final_terminal_rows = count(
        row -> row["current_status"] in QETLABCompletion.FINAL_PUBLIC_TERMINAL_STATUSES,
        model.public_rows,
    )
    static_failures = length(audit.failures)
    strict_pass =
        final_terminal_rows == model.metadata.public_row_count &&
        iszero(required_helpers) &&
        iszero(static_failures)
    return Dict{String,Any}(
        "source_revision" => model.metadata.source_revision,
        "public_rows" => model.metadata.public_row_count,
        "implementation_complete_rows" => model.counts["implementation_complete"],
        "final_terminal_rows" => final_terminal_rows,
        "partial_rows" => model.counts["partial"],
        "deferred_rows" => model.counts["deferred"],
        "blocked_rows" => model.counts["blocked"],
        "queued_public_rows" => model.metadata.incomplete_public_count,
        "internal_helpers" => model.metadata.internal_helper_count,
        "terminal_helpers" => terminal_helpers,
        "required_helpers" => required_helpers,
        "static_failure_count" => static_failures,
        "static_failures" => copy(audit.failures),
        "strict_pass" => strict_pass,
    )
end

function provenance_gaps(provenance_entries, exported_keys)
    by_key = Dict{Tuple{String,String},Any}()
    for entry in provenance_entries
        key = (
            string(get(entry, "public_module", "")), string(get(entry, "julia_name", ""))
        )
        by_key[key] = entry
    end
    provenance_keys = Set(keys(by_key))
    missing_provenance = sort!([
        "$module_name.$name" for
        (module_name, name) in setdiff(exported_keys, provenance_keys)
    ])
    nonexported_provenance = sort!([
        "$module_name.$name" for
        (module_name, name) in setdiff(provenance_keys, exported_keys)
    ])
    missing_tests = String[]
    missing_docs = String[]
    missing_test_paths = String[]
    missing_doc_paths = String[]
    for (key, entry) in by_key
        qualified = "$(key[1]).$(key[2])"
        tests = get(entry, "tests", Any[])
        docs = get(entry, "docs", Any[])
        tests isa AbstractVector && !isempty(tests) || push!(missing_tests, qualified)
        docs isa AbstractVector && !isempty(docs) || push!(missing_docs, qualified)
        for path in (tests isa AbstractVector ? tests : Any[])
            ispath(joinpath(REPOSITORY_ROOT, string(path))) ||
                push!(missing_test_paths, "$qualified => $path")
        end
        for path in (docs isa AbstractVector ? docs : Any[])
            startswith(string(path), "http") && continue
            ispath(joinpath(REPOSITORY_ROOT, string(path))) ||
                push!(missing_doc_paths, "$qualified => $path")
        end
    end
    return Dict{String,Any}(
        "exports_without_provenance_count" => length(missing_provenance),
        "exports_without_provenance" => missing_provenance,
        "provenance_without_export_count" => length(nonexported_provenance),
        "provenance_without_export" => nonexported_provenance,
        "exports_without_tests_count" => length(missing_tests),
        "exports_without_tests" => sort!(missing_tests),
        "exports_without_docs_count" => length(missing_docs),
        "exports_without_docs" => sort!(missing_docs),
        "missing_test_path_count" => length(missing_test_paths),
        "missing_test_paths" => sort!(missing_test_paths),
        "missing_doc_path_count" => length(missing_doc_paths),
        "missing_doc_paths" => sort!(missing_doc_paths),
    )
end

function parse_number(value)
    return parse(Int, replace(value, "," => ""))
end

function captured_number(matched::RegexMatch)
    for capture in matched.captures
        isnothing(capture) || return parse_number(capture)
    end
    return error("numeric claim regex matched without a capture")
end

function current_claim_content(relative::AbstractString)
    content = read(joinpath(REPOSITORY_ROOT, relative), String)
    begin_count = count(CURRENT_CLAIMS_BEGIN, content)
    end_count = count(CURRENT_CLAIMS_END, content)
    begin_count == end_count ||
        error("$relative has unbalanced current-claim reconciliation markers")
    iszero(begin_count) && return content

    blocks = String[]
    offset = firstindex(content)
    while true
        opening = findnext(CURRENT_CLAIMS_BEGIN, content, offset)
        opening === nothing && break
        body_start = nextind(content, last(opening))
        closing = findnext(CURRENT_CLAIMS_END, content, body_start)
        closing === nothing &&
            error("$relative has an unclosed current-claim reconciliation block")
        push!(blocks, content[body_start:prevind(content, first(closing))])
        offset = nextind(content, last(closing))
    end
    return join(blocks, "\n")
end

function validate_claims(snapshot)
    completion = snapshot["completion"]
    expected = Dict(
        "inventory_rows" => snapshot["inventory"]["rows"],
        "public_rows" => snapshot["inventory"]["public_rows"],
        "internal_rows" => snapshot["inventory"]["internal_rows"],
        "dependency_edges" => snapshot["inventory"]["dependency_edges"],
        "implemented" => snapshot["inventory"]["public_implementation"]["implemented"],
        "partial" => snapshot["inventory"]["public_implementation"]["partial"],
        "deferred" => snapshot["inventory"]["public_implementation"]["deferred"],
        "blocked" => snapshot["inventory"]["public_implementation"]["blocked"],
        "bindings" => snapshot["api"]["total_exports"],
        "integrated_assertions" =>
            snapshot["dynamic_tests"]["package"]["integrated_assertions"],
        "core_assertions" => snapshot["dynamic_tests"]["package"]["core_assertions"],
        "tutorial_assertions" =>
            snapshot["dynamic_tests"]["package"]["tutorial_assertions"],
        "benchmark_cases" => snapshot["benchmarks"]["count"],
        "extension_assertions" =>
            snapshot["dynamic_tests"]["entanglement_detection_extension"]["assertions"],
        "completion_rows" => completion["final_terminal_rows"],
        "completion_queue" => completion["queued_public_rows"],
        "terminal_helpers" => completion["terminal_helpers"],
        "required_helpers" => completion["required_helpers"],
        "completion_failures" => completion["static_failure_count"],
    )
    patterns = [
        (
            "inventory_rows",
            r"(?i)\b(\d[\d,]*) (?:inventoried QETLAB files|inventory rows|source-reviewed rows)\b",
        ),
        ("public_rows", r"(?i)\b(?:among (?:the |its )|of the )(\d[\d,]*) public rows\b"),
        ("internal_rows", r"(?i)\b(\d[\d,]*) (?:private|internal) helpers\b"),
        ("dependency_edges", r"(?i)\b(\d[\d,]*) dependency edges\b"),
        ("bindings", r"(?i)\b(\d[\d,]*) public bindings\b"),
        (
            "integrated_assertions",
            r"(?i)\bintegrated (\d[\d,]*)-assertion package corpus\b",
        ),
        (
            "integrated_assertions",
            r"(?i)\b(?:a |the )?(\d[\d,]*)-assertion (?:full |local )?package (?:suite|corpus)\b",
        ),
        ("core_assertions", r"(?i)\b(\d[\d,]*) core assertions\b"),
        ("core_assertions", r"(?i)\b(\d[\d,]*) core plus\b"),
        (
            "tutorial_assertions",
            r"(?i)\b(\d[\d,]*) (?:executable[- ]tutorial|tutorial) assertions\b",
        ),
        ("benchmark_cases", r"(?i)\b(\d[\d,]*)(?:-case)? declared quick benchmark cases\b"),
        (
            "extension_assertions",
            r"(?i)\bEntanglementDetection(?:\.jl)?(?: extension)?[^\n]{0,160}?\b(\d[\d,]*)(?:/\d[\d,]*)? (?:focused )?assertions\b",
        ),
        (
            "completion_rows",
            r"(?i)\b(\d[\d,]*)/127 public rows (?:are )?(?:implementation-complete|verified)\b",
        ),
        (
            "completion_queue",
            r"(?i)\bcompletion queue(?: contains| has|:)? (\d[\d,]*) public rows\b",
        ),
        (
            "terminal_helpers",
            r"(?i)\b(\d[\d,]*)/36 internal helpers (?:have|with) terminal dispositions\b",
        ),
        ("required_helpers", r"(?i)\b(\d[\d,]*) required (?:internal )?helpers remain\b"),
        ("completion_failures", r"(?i)\b(\d[\d,]*) static completion failures\b"),
    ]
    status_pattern = r"(?i)\b(\d[\d,]*) (?:mappings? (?:are )?)?implemented,\s*(\d[\d,]*) (?:are )?(?:explicitly )?partial,\s*(\d[\d,]*) (?:are )?deferred,\s*(?:and )?(\d[\d,]*) (?:are )?blocked\b"

    mismatches = Dict{String,Any}[]
    recognized = 0
    by_file = Dict{String,Int}()
    for relative in MAINTAINED_CLAIM_FILES
        content = current_claim_content(relative)
        file_count = 0
        for (metric, pattern) in patterns
            for matched in eachmatch(pattern, content)
                file_count += 1
                recognized += 1
                actual = captured_number(matched)
                actual == expected[metric] && continue
                push!(
                    mismatches,
                    Dict{String,Any}(
                        "path" => relative,
                        "metric" => metric,
                        "claimed" => actual,
                        "computed" => expected[metric],
                        "text" => matched.match,
                    ),
                )
            end
        end
        for matched in eachmatch(status_pattern, content)
            file_count += 1
            recognized += 1
            values = parse_number.(matched.captures)
            metrics = ("implemented", "partial", "deferred", "blocked")
            for (metric, actual) in zip(metrics, values)
                actual == expected[metric] && continue
                push!(
                    mismatches,
                    Dict{String,Any}(
                        "path" => relative,
                        "metric" => metric,
                        "claimed" => actual,
                        "computed" => expected[metric],
                        "text" => matched.match,
                    ),
                )
            end
        end
        by_file[relative] = file_count
    end
    sort!(mismatches; by=entry -> (entry["path"], entry["metric"], entry["text"]))
    return Dict{String,Any}(
        "checked_files" => collect(MAINTAINED_CLAIM_FILES),
        "recognized_claim_count" => recognized,
        "recognized_claims_by_file" =>
            Dict{String,Any}(key => by_file[key] for key in sort!(collect(keys(by_file)))),
        "mismatch_count" => length(mismatches),
        "mismatches" => mismatches,
    )
end

function relevant_input_paths()
    paths = String[]
    fixed = [
        "Project.toml",
        "PROVENANCE.toml",
        "UpstreamManifest.toml",
        "porting/qetlab_inventory.toml",
        "porting/qetlab_status.toml",
        "porting/qetlab_completion_policy.toml",
        "porting/qetlab_completion_plan.toml",
        "porting/qetlab_completion_queue.toml",
        "benchmark/benchmarks.jl",
        "tutorials/runtests.jl",
        "scripts/qetlab_completion_common.jl",
        "scripts/reconcile_project_claims.jl",
    ]
    append!(paths, fixed)
    append!(paths, MAINTAINED_CLAIM_FILES)
    for root_name in ("src", "ext", "test", "tutorials", ".github/workflows")
        root = joinpath(REPOSITORY_ROOT, root_name)
        for (directory, _, files) in walkdir(root)
            for file in files
                path = joinpath(directory, file)
                relative = relpath(path, REPOSITORY_ROOT)
                (
                    endswith(file, ".jl") ||
                    endswith(file, ".yml") ||
                    endswith(file, ".toml")
                ) || continue
                push!(paths, relative)
            end
        end
    end
    return sort!(unique(paths))
end

function input_digest()
    buffer = IOBuffer()
    for relative in relevant_input_paths()
        path = joinpath(REPOSITORY_ROOT, relative)
        isfile(path) || continue
        write(buffer, relative)
        write(buffer, UInt8(0))
        write(buffer, read(path))
        write(buffer, UInt8(0))
    end
    return bytes2hex(sha256(take!(buffer)))
end

function load_preserved_dynamic_tests(output)
    isfile(output) || error(
        "--skip-tests requires an existing snapshot at $(relpath(output, REPOSITORY_ROOT))",
    )
    previous = TOML.parsefile(output)
    haskey(previous, "dynamic_tests") ||
        error("existing snapshot has no dynamic_tests section")
    return previous["dynamic_tests"]
end

function build_snapshot(options)
    inventory = TOML.parsefile(INVENTORY_PATH)
    provenance = TOML.parsefile(PROVENANCE_PATH)
    project = TOML.parsefile(PROJECT_PATH)
    inventory_entries = get(inventory, "functions", Any[])
    provenance_entries = get(provenance, "functions", Any[])
    public_entries = filter(
        entry -> get(entry, "classification", "") == "public", inventory_entries
    )
    internal_entries = filter(
        entry -> get(entry, "classification", "") == "internal", inventory_entries
    )

    native_exports = public_names(QuantumEntanglementTools)
    compat_exports = public_names(QuantumEntanglementTools.MATLABCompat)
    exported_keys = Set{Tuple{String,String}}()
    union!(
        exported_keys,
        (("QuantumEntanglementTools", name) for name in native_exports),
        (("QuantumEntanglementTools.MATLABCompat", name) for name in compat_exports),
    )

    qetlab_export_count = count(
        entry -> get(entry, "source_project", "") == "QETLAB", provenance_entries
    )
    project_native_export_count = count(
        entry -> get(entry, "source_project", "") == "QuantumEntanglementTools",
        provenance_entries,
    )
    status_counts = Dict{String,Any}()
    for field in STATUS_FIELDS
        status_counts[field] = count_values(inventory_entries, field)
    end

    dynamic_tests = if options.skip_tests
        load_preserved_dynamic_tests(options.output)
    else
        run_dynamic_tests()
    end
    completion_model = QETLABCompletion.build_completion_model()
    completion_audit = QETLABCompletion.repository_quality_audit(
        completion_model; runtime_modules=PUBLIC_MODULES
    )
    snapshot = Dict{String,Any}(
        "schema_version" => 1,
        "generation" => Dict{String,Any}(
            "script" => "scripts/reconcile_project_claims.jl",
            "deterministic" => true,
            "input_sha256" => input_digest(),
            "inventory_sha256" => bytes2hex(sha256(read(INVENTORY_PATH))),
            "status_overlay_sha256" => bytes2hex(sha256(read(STATUS_PATH))),
            "provenance_sha256" => bytes2hex(sha256(read(PROVENANCE_PATH))),
        ),
        "project" => Dict{String,Any}(
            "name" => get(project, "name", ""),
            "uuid" => get(project, "uuid", ""),
            "version" => get(project, "version", ""),
            "julia_compat" => get(get(project, "compat", Dict()), "julia", ""),
            "core_dependency_count" => length(get(project, "deps", Dict())),
            "weak_dependency_count" => length(get(project, "weakdeps", Dict())),
            "extension_count" => length(get(project, "extensions", Dict())),
        ),
        "inventory" => Dict{String,Any}(
            "source_revision" => get(inventory, "source_revision", ""),
            "review_status" => get(inventory, "review_status", ""),
            "rows" => length(inventory_entries),
            "public_rows" => length(public_entries),
            "internal_rows" => length(internal_entries),
            "dependency_edges" => get(inventory, "dependency_edge_count", 0),
            "dependency_cycles" => get(inventory, "dependency_cycle_count", 0),
            "source_reviewed_rows" => get(inventory, "source_reviewed_count", 0),
            "pending_review_rows" => get(inventory, "pending_review_count", 0),
            "public_implementation" => count_buckets(
                public_entries,
                implementation_bucket;
                expected=["implemented", "partial", "deferred", "blocked", "other"],
            ),
            "internal_implementation" => count_buckets(
                internal_entries,
                internal_implementation_bucket;
                expected=[
                    "superseded",
                    "intentionally_unexposed",
                    "partial",
                    "deferred",
                    "other",
                ],
            ),
            "status_counts" => status_counts,
        ),
        "api" => Dict{String,Any}(
            "native_exports" => length(native_exports),
            "matlab_compat_exports" => length(compat_exports),
            "total_exports" => length(exported_keys),
            "provenance_entries" => length(provenance_entries),
            "qetlab_provenance_bindings" => qetlab_export_count,
            "project_native_bindings" => project_native_export_count,
            "gaps" => provenance_gaps(provenance_entries, exported_keys),
        ),
        "completion" => completion_summary(completion_model, completion_audit),
        "dynamic_tests" => dynamic_tests,
        "static_test_macros" => count_test_macros(),
        "tutorials" => tutorial_summary(),
        "benchmarks" => benchmark_summary(),
        "source_audit" => source_audit(provenance_entries),
    )
    snapshot["maintained_claims"] = validate_claims(snapshot)
    return snapshot
end

function render_snapshot(snapshot)
    buffer = IOBuffer()
    println(buffer, "# This file is generated by scripts/reconcile_project_claims.jl.")
    println(
        buffer, "# Edit source ledgers/tests/claims, then regenerate; do not edit by hand."
    )
    TOML.print(buffer, snapshot; sorted=true)
    return String(take!(buffer))
end

function write_or_check(options, rendered)
    relative = relpath(options.output, REPOSITORY_ROOT)
    if options.check
        isfile(options.output) || error("generated snapshot is missing: $relative")
        existing = read(options.output, String)
        existing == rendered || error(
            "generated snapshot is stale: $relative\n" *
            "run `julia --startup-file=no --project=. scripts/reconcile_project_claims.jl`",
        )
        println("project-claim snapshot is current: $relative")
    else
        mkpath(dirname(options.output))
        write(options.output, rendered)
        println("wrote project-claim snapshot: $relative")
    end
    return nothing
end

function main(args)
    options = parse_options(args)
    snapshot = build_snapshot(options)
    claim_mismatches = snapshot["maintained_claims"]["mismatch_count"]
    if claim_mismatches != 0
        println(stderr, "maintained quantitative claim mismatches:")
        for mismatch in snapshot["maintained_claims"]["mismatches"]
            println(
                stderr,
                "  - ",
                mismatch["path"],
                ": ",
                mismatch["metric"],
                " claims ",
                mismatch["claimed"],
                ", computed ",
                mismatch["computed"],
                " (",
                repr(mismatch["text"]),
                ")",
            )
        end
        error("$claim_mismatches maintained quantitative claim(s) disagree with evidence")
    end
    write_or_check(options, render_snapshot(snapshot))
    inventory = snapshot["inventory"]
    api = snapshot["api"]
    tests = snapshot["dynamic_tests"]
    println(
        "reconciled ",
        inventory["rows"],
        " inventory rows (",
        inventory["public_rows"],
        " public, ",
        inventory["internal_rows"],
        " internal), ",
        api["total_exports"],
        " exports/provenance entries, ",
        tests["package"]["integrated_assertions"],
        " package assertions, and ",
        tests["entanglement_detection_extension"]["assertions"],
        " extension assertions",
    )
    return nothing
end

main(ARGS)
