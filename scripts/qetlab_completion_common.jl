module QETLABCompletion

using SHA
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const INVENTORY_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_inventory.toml")
const PROVENANCE_PATH = joinpath(REPOSITORY_ROOT, "PROVENANCE.toml")
const POLICY_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_completion_policy.toml")
const PLAN_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_completion_plan.toml")
const QUEUE_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_completion_queue.toml")
const DOCUMENT_PATH = joinpath(REPOSITORY_ROOT, "docs", "QETLAB_COMPLETION_PLAN.md")

const GENERATOR_PATH = "scripts/build_qetlab_completion_plan.jl"
const CHECKER_PATH = "scripts/check_qetlab_completion.jl"

const FINAL_PUBLIC_TERMINAL_STATUSES = Set([
    "superseded_with_documented_mapping",
    "verified",
    "verified_with_documented_upstream_correction",
])

const IMPLEMENTATION_COMPLETE_STATUSES = Set([
    "compatibility_alias", "implemented", FINAL_PUBLIC_TERMINAL_STATUSES...
])

const PARTIAL_PUBLIC_STATUSES = Set([
    "partially_implemented_bipartite_scalar_only",
    "partially_implemented_certificate_first_subset_only",
    "partially_implemented_completely_positive_maps_only",
    "partially_implemented_completely_positive_square_operator_spaces_only",
    "partially_implemented_exact_pure_bipartite_and_mixed_two_qubit_domains_only",
    "partially_implemented_numeric_matrices_only_cvx_symbolic_branch_omitted",
    "partially_implemented_special_cases_via_verified_norm_apis",
    "partially_implemented_square_operator_spaces_and_no_two_sided_kraus_input",
    "partially_implemented_square_operator_spaces_and_one_sided_cp_kraus_only",
    "partially_implemented_square_operator_spaces_and_one_sided_maps_only",
    "partially_implemented_typed_representation_diagnostic_only",
    "partially_implemented_via_schmidt_coefficients_without_dedicated_binding",
    "partially_implemented_von_neumann_alpha_one_only",
    "partially_implemented_without_hermitian_factor_repair_or_matlab_positional_output_parity",
])

const ALLOWED_PUBLIC_STATUSES = union(
    IMPLEMENTATION_COMPLETE_STATUSES,
    PARTIAL_PUBLIC_STATUSES,
    Set(["blocked_with_explicit_reason", "deferred"]),
)

const HELPER_STATUS_CLASS = Dict(
    "deferred_with_bcs_game_scope" => "required_deferred",
    "deferred_with_bell_nonlocal_game_scope" => "required_deferred",
    "deferred_with_commutant_scope" => "required_deferred",
    "deferred_with_distinguishability_upb_and_sk_scope" => "required_deferred",
    "deferred_with_k_incoherence_scope" => "required_deferred",
    "deferred_with_polynomial_optimization_scope" => "required_deferred",
    "deferred_with_sk_norm_and_extension_scope" => "required_deferred",
    "deferred_with_sk_norm_block_positivity_and_separability_scope" => "required_deferred",
    "deferred_with_upb_scope" => "required_deferred",
    "intentionally_not_ported_unreferenced_upstream_private_helper" => "terminal_excluded_private",
    "internal_helper_superseded_by_base_invperm" => "terminal_superseded",
    "internal_helper_superseded_by_caller_specific_factorizations" => "terminal_superseded",
    "internal_helper_superseded_by_checked_qetlab_monomial_index" => "terminal_superseded",
    "internal_helper_superseded_by_commutant_specific_bounded_dense_svd" => "terminal_superseded",
    "internal_helper_superseded_by_corrected_bounded_graph_bandwidth_search" => "terminal_superseded",
    "internal_helper_superseded_by_explicit_rng_bounded_native_sampling" => "terminal_superseded",
    "internal_helper_superseded_by_explicit_rng_bounded_sk_projected_search" => "terminal_superseded",
    "internal_helper_superseded_by_function_specific_state_dispatch" => "terminal_superseded",
    "internal_helper_superseded_by_julia_concatenation_and_fill" => "terminal_superseded",
    "internal_helper_superseded_by_julia_identity_constructors" => "terminal_superseded",
    "internal_helper_superseded_by_julia_methods_defaults_and_keywords" => "terminal_superseded",
    "internal_helper_superseded_by_parent_specific_construction_and_strict_normalization_validation" => "terminal_superseded",
    "internal_helper_superseded_by_parent_specific_iteration" => "terminal_superseded",
    "internal_helper_superseded_by_private_parent_algorithm" => "terminal_superseded",
    "internal_helper_superseded_by_public_typed_bell_behavior_conversion" => "terminal_superseded",
    "internal_helper_superseded_by_public_typed_nonlocal_game_conversion" => "terminal_superseded",
    "internal_helper_superseded_by_public_checked_monomial_exponents" => "terminal_superseded",
    "internal_helper_superseded_by_structured_results_and_standard_io" => "terminal_superseded",
    "internal_helper_superseded_by_typed_parent_specific_canonical_nonlocal_representation" => "terminal_superseded",
    "internal_helper_superseded_by_typed_operator_space_dimensions" => "terminal_superseded",
    "internal_helper_superseded_with_intentional_index_convention_change" => "terminal_superseded",
    "partially_superseded_deferred_general_map_remainder" => "required_partial",
    "partially_superseded_deferred_nonlocal_and_repetition_remainder" => "required_partial",
    "partially_superseded_deferred_polynomial_remainder" => "required_partial",
)

const REQUIRED_PUBLIC_ROW_FIELDS = Set([
    "upstream_function",
    "source_path",
    "current_status",
    "target_status",
    "work_package",
    "native_binding",
    "compatibility_binding",
    "required_backend",
    "prerequisite_rows",
    "specification",
    "known_upstream_defects",
    "implementation_files",
    "test_files",
    "oracle_files",
    "documentation_files",
    "benchmark_files",
    "verification_status",
])

const REQUIRED_PROVENANCE_FIELDS = Set([
    "julia_name",
    "public_module",
    "public_kind",
    "julia_file",
    "source_project",
    "specification",
    "implementation_kind",
    "status",
    "verification_status",
    "tests",
    "docs",
])

const REQUIRED_QETLAB_PROVENANCE_FIELDS = Set([
    "upstream_function", "source_path", "source_revision", "source_sha256"
])

const ALLOWED_PROVENANCE_STATUSES = Set(["implemented", "verified", "compatibility_alias"])
const ALLOWED_PUBLIC_KINDS = Set([
    "module", "type", "function", "constant", "alias", "compatibility_wrapper"
])
const PUBLIC_MODULE_NAMES = (
    "QuantumEntanglementTools", "QuantumEntanglementTools.MATLABCompat"
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
const SOURCE_PLACEHOLDER_PATTERN = r"(?i)\b(TODO|FIXME|XXX|not[ _-]?implemented|placeholder)\b"

# These rows require a solver-aware result contract and an optional extension
# environment before a future implementation-complete disposition is accepted.
# The list is the reviewed CVX/optimization subset of the pinned inventory;
# accepting a row here does not select a Julia optimizer.
const SOLVER_BACKED_PUBLIC_ROWS = Set([
    "AbsPPTConstraints",
    "BCSGameLB",
    "BCSGameValue",
    "BellInequalityMax",
    "BellInequalityMaxQubits",
    "CBNorm",
    "ChannelDistinguishability",
    "CliqueNumber",
    "DiamondNorm",
    "Distinguishability",
    "GenRobustnesskCoherence",
    "IsAbskIncoh",
    "IsAbsPPT",
    "IsBlockPositive",
    "IsCopositive",
    "IsSeparable",
    "IskIncoherent",
    "kpNorm",
    "kpNormDual",
    "LocalDistinguishability",
    "MaximumOutputFidelity",
    "NonlocalGameLB",
    "NPAHierarchy",
    "PolynomialSOS",
    "RobustnessCoherence",
    "SkOperatorNorm",
    "SymmetricExtension",
    "SymmetricInnerExtension",
    "TraceDistanceCoherence",
    "UPBSepDistinguishable",
    "XORGameValue",
])

function normalized_path(path::AbstractString)
    return replace(String(path), '\\' => '/')
end

function file_sha256(path::AbstractString)
    isfile(path) || error("required completion input does not exist: $path")
    return bytes2hex(sha256(read(path)))
end

function unique_sorted(values)
    return sort!(unique!(String[String(value) for value in values]); by=lowercase)
end

function required_string(table, key::AbstractString, context::AbstractString)
    value = get(table, key, nothing)
    value isa AbstractString && !isempty(strip(value)) ||
        error("$context must provide a non-empty string `$key`")
    return String(value)
end

function required_string_vector(table, key::AbstractString, context::AbstractString)
    value = get(table, key, nothing)
    value isa AbstractVector && all(item -> item isa AbstractString, value) ||
        error("$context must provide an array of strings `$key`")
    return String.(value)
end

function public_status_class(status::AbstractString)
    status in IMPLEMENTATION_COMPLETE_STATUSES && return "implementation_complete"
    status in PARTIAL_PUBLIC_STATUSES && return "partial"
    status == "deferred" && return "deferred"
    status == "blocked_with_explicit_reason" && return "blocked"
    return error("unknown public implementation status: $status")
end

function helper_status_class(status::AbstractString)
    haskey(HELPER_STATUS_CLASS, status) ||
        error("unknown internal-helper implementation status: $status")
    return HELPER_STATUS_CLASS[status]
end

function is_required_helper(status_class::AbstractString)
    return startswith(status_class, "required_")
end

function qualified_binding(entry)
    return String(entry["public_module"]) * "." * String(entry["julia_name"])
end

function provenance_paths(entries, key::AbstractString)
    paths = String[]
    for entry in entries
        values = get(entry, key, Any[])
        values isa AbstractVector || error("provenance field `$key` must be an array")
        append!(paths, String.(values))
    end
    return unique_sorted(paths)
end

function provenance_strings(entries, key::AbstractString)
    values = String[]
    for entry in entries
        haskey(entry, key) || continue
        value = entry[key]
        value isa AbstractString || error("provenance field `$key` must be a string")
        push!(values, String(value))
    end
    return unique_sorted(values)
end

function dependency_closure(name::AbstractString, graph)
    seen = Set{String}()
    pending = copy(get(graph, String(name), String[]))
    while !isempty(pending)
        dependency = pop!(pending)
        dependency in seen && continue
        push!(seen, dependency)
        append!(pending, get(graph, dependency, String[]))
    end
    return seen
end

function call_graph_ranks(incomplete_names, inventory_by_name)
    incomplete = Set(String.(incomplete_names))
    ranks = Dict{String,Int}()
    visiting = Set{String}()

    function rank(name::String)
        haskey(ranks, name) && return ranks[name]
        name in visiting && error("cycle in incomplete-public source-call graph at $name")
        push!(visiting, name)
        dependencies = [
            dependency for dependency in inventory_by_name[name]["upstream_calls"] if
            dependency in incomplete
        ]
        value = isempty(dependencies) ? 0 : 1 + maximum(rank.(dependencies))
        delete!(visiting, name)
        ranks[name] = value
        return value
    end

    for name in incomplete
        rank(name)
    end
    return ranks
end

function validate_policy(policy, incomplete_names, required_helper_names)
    get(policy, "schema_version", nothing) == 1 ||
        error("completion policy schema_version must equal 1")
    allowed_top_level = Set(["schema_version", "baseline_capabilities", "work_packages"])
    unknown_top_level = setdiff(Set(keys(policy)), allowed_top_level)
    isempty(unknown_top_level) || error(
        "unknown completion-policy top-level fields: " *
        join(sort!(collect(unknown_top_level)), ", "),
    )

    baseline = required_string_vector(policy, "baseline_capabilities", "completion policy")
    length(baseline) == length(unique(baseline)) ||
        error("completion policy repeats a baseline capability")
    packages_raw = get(policy, "work_packages", nothing)
    packages_raw isa AbstractVector ||
        error("completion policy must provide `work_packages`")
    length(packages_raw) == 12 ||
        error("completion policy must define exactly 12 work packages")

    allowed_package_fields = Set([
        "id",
        "order",
        "title",
        "summary",
        "architecture_only",
        "prerequisite_capabilities",
        "provides_capabilities",
        "functions",
        "internal_helpers",
    ])
    packages = NamedTuple[]
    for (index, package) in pairs(packages_raw)
        package isa AbstractDict || error("completion work package $index is not a table")
        unknown = setdiff(Set(keys(package)), allowed_package_fields)
        isempty(unknown) || error(
            "unknown fields in completion work package $index: " *
            join(sort!(collect(unknown)), ", "),
        )
        id = required_string(package, "id", "completion work package $index")
        title = required_string(package, "title", "completion work package $id")
        summary = required_string(package, "summary", "completion work package $id")
        order = get(package, "order", nothing)
        order isa Integer ||
            error("completion work package $id must provide integer `order`")
        architecture_only = get(package, "architecture_only", nothing)
        architecture_only isa Bool ||
            error("completion work package $id must provide Boolean `architecture_only`")
        prerequisites = required_string_vector(
            package, "prerequisite_capabilities", "completion work package $id"
        )
        provides = required_string_vector(
            package, "provides_capabilities", "completion work package $id"
        )
        functions = required_string_vector(
            package, "functions", "completion work package $id"
        )
        internal_helpers = required_string_vector(
            package, "internal_helpers", "completion work package $id"
        )
        push!(
            packages,
            (;
                id,
                order=Int(order),
                title,
                summary,
                architecture_only,
                prerequisite_capabilities=unique_sorted(prerequisites),
                provides_capabilities=unique_sorted(provides),
                functions=String.(functions),
                internal_helpers=String.(internal_helpers),
            ),
        )
    end

    ids = [package.id for package in packages]
    length(ids) == length(unique(ids)) ||
        error("completion policy repeats a work-package id")
    orders = [package.order for package in packages]
    sort(orders) == collect(1:length(packages)) ||
        error("completion work-package orders must be exactly 1:$(length(packages))")
    sort!(packages; by=package -> package.order)

    available = Set(baseline)
    provided = Set{String}()
    for package in packages
        missing = setdiff(Set(package.prerequisite_capabilities), available)
        isempty(missing) || error(
            "completion work package $(package.id) has unavailable prerequisite " *
            "capabilities: $(join(sort!(collect(missing)), ", "))",
        )
        duplicates = intersect(Set(package.provides_capabilities), provided)
        isempty(duplicates) || error(
            "completion work package $(package.id) repeats provided capabilities: " *
            join(sort!(collect(duplicates)), ", "),
        )
        union!(provided, package.provides_capabilities)
        union!(available, package.provides_capabilities)
    end

    assignments = Dict{String,String}()
    for package in packages, function_name in package.functions
        haskey(assignments, function_name) &&
            error("completion policy assigns $function_name more than once")
        assignments[function_name] = package.id
    end
    expected = Set(String.(incomplete_names))
    assigned = Set(keys(assignments))
    missing = sort!(collect(setdiff(expected, assigned)); by=lowercase)
    unexpected = sort!(collect(setdiff(assigned, expected)); by=lowercase)
    isempty(missing) ||
        error("completion policy omits incomplete public rows: " * join(missing, ", "))
    isempty(unexpected) || error(
        "completion policy assigns rows that are not incomplete: " * join(unexpected, ", "),
    )

    helper_assignments = Dict{String,String}()
    for package in packages, helper_name in package.internal_helpers
        haskey(helper_assignments, helper_name) &&
            error("completion policy assigns internal helper $helper_name more than once")
        helper_assignments[helper_name] = package.id
    end
    expected_helpers = Set(String.(required_helper_names))
    assigned_helpers = Set(keys(helper_assignments))
    missing_helpers = sort!(
        collect(setdiff(expected_helpers, assigned_helpers)); by=lowercase
    )
    unexpected_helpers = sort!(
        collect(setdiff(assigned_helpers, expected_helpers)); by=lowercase
    )
    isempty(missing_helpers) || error(
        "completion policy omits required internal helpers: " * join(missing_helpers, ", "),
    )
    isempty(unexpected_helpers) || error(
        "completion policy assigns terminal or unknown internal helpers: " *
        join(unexpected_helpers, ", "),
    )
    return (;
        baseline_capabilities=unique_sorted(baseline),
        packages,
        assignments,
        helper_assignments,
    )
end

function build_completion_model()
    inventory = TOML.parsefile(INVENTORY_PATH)
    provenance = TOML.parsefile(PROVENANCE_PATH)
    policy_raw = TOML.parsefile(POLICY_PATH)

    get(inventory, "schema_version", nothing) == 1 ||
        error("qetlab_inventory.toml schema_version must equal 1")
    get(provenance, "schema_version", nothing) == 1 ||
        error("PROVENANCE.toml schema_version must equal 1")
    inventory_entries = get(inventory, "functions", nothing)
    inventory_entries isa AbstractVector ||
        error("qetlab_inventory.toml must provide `functions`")
    provenance_entries = get(provenance, "functions", nothing)
    provenance_entries isa AbstractVector ||
        error("PROVENANCE.toml must provide `functions`")

    inventory_by_name = Dict{String,Any}()
    for entry in inventory_entries
        name = required_string(entry, "function_name", "inventory row")
        haskey(inventory_by_name, name) && error("duplicate inventory function_name: $name")
        inventory_by_name[name] = entry
    end

    public_inventory = [
        entry for entry in inventory_entries if entry["classification"] == "public"
    ]
    internal_inventory = [
        entry for entry in inventory_entries if entry["classification"] == "internal"
    ]
    length(public_inventory) == get(inventory, "public_function_count", nothing) ||
        error("inventory public_function_count is stale")
    length(internal_inventory) == get(inventory, "internal_helper_count", nothing) ||
        error("inventory internal_helper_count is stale")

    for entry in public_inventory
        status = required_string(entry, "implementation_status", "public inventory row")
        status in ALLOWED_PUBLIC_STATUSES || error(
            "unknown public implementation status for $(entry["function_name"]): $status",
        )
    end
    for entry in internal_inventory
        status = required_string(entry, "implementation_status", "internal inventory row")
        helper_status_class(status)
    end

    incomplete_inventory = [
        entry for entry in public_inventory if
        !(entry["implementation_status"] in IMPLEMENTATION_COMPLETE_STATUSES)
    ]
    incomplete_names = String[entry["function_name"] for entry in incomplete_inventory]
    required_internal_names = Set(
        String(entry["function_name"]) for entry in internal_inventory if
        is_required_helper(helper_status_class(entry["implementation_status"]))
    )
    policy = validate_policy(policy_raw, incomplete_names, required_internal_names)
    package_by_id = Dict(package.id => package for package in policy.packages)

    provenance_by_upstream = Dict{String,Vector{Any}}()
    for entry in provenance_entries
        get(entry, "source_project", "") == "QETLAB" || continue
        upstream_function = required_string(
            entry, "upstream_function", "QETLAB provenance row"
        )
        haskey(inventory_by_name, upstream_function) ||
            error("QETLAB provenance names an unknown inventory row: $upstream_function")
        push!(get!(provenance_by_upstream, upstream_function, Any[]), entry)
    end

    public_names = Set(String[entry["function_name"] for entry in public_inventory])
    internal_names = Set(String[entry["function_name"] for entry in internal_inventory])
    graph = Dict(
        String(entry["function_name"]) => String.(entry["upstream_calls"]) for
        entry in inventory_entries
    )
    ranks = call_graph_ranks(incomplete_names, inventory_by_name)

    public_rows = Dict{String,Any}[]
    for entry in sort!(public_inventory; by=row -> lowercase(row["function_name"]))
        function_name = String(entry["function_name"])
        mapped = get(provenance_by_upstream, function_name, Any[])
        native_entries = [
            provenance_entry for provenance_entry in mapped if
            provenance_entry["public_module"] == "QuantumEntanglementTools"
        ]
        compatibility_entries = [
            provenance_entry for provenance_entry in mapped if
            provenance_entry["public_module"] == "QuantumEntanglementTools.MATLABCompat"
        ]
        all_test_paths = provenance_paths(mapped, "tests")
        oracle_files = [path for path in all_test_paths if startswith(path, "test/oracle/")]
        test_files = [path for path in all_test_paths if !startswith(path, "test/oracle/")]
        specification = provenance_strings(mapped, "specification")
        documentation_url = String(get(entry, "documentation_url", ""))
        isempty(documentation_url) || push!(specification, documentation_url)
        specification = unique_sorted(specification)
        benchmark_files = if occursin("recorded", entry["benchmark_status"])
            ["benchmark/benchmarks.jl"]
        else
            String[]
        end
        status = String(entry["implementation_status"])
        work_package = if status in IMPLEMENTATION_COMPLETE_STATUSES
            "completed_scope"
        else
            policy.assignments[function_name]
        end
        prerequisite_rows = unique_sorted(
            dependency for
            dependency in entry["upstream_calls"] if dependency in public_names
        )
        row = Dict{String,Any}(
            "upstream_function" => function_name,
            "source_path" => String(entry["source_path"]),
            "current_status" => status,
            "target_status" => "verified",
            "work_package" => work_package,
            "native_binding" => unique_sorted(qualified_binding.(native_entries)),
            "compatibility_binding" =>
                unique_sorted(qualified_binding.(compatibility_entries)),
            "required_backend" => unique_sorted(entry["external_toolboxes"]),
            "prerequisite_rows" => prerequisite_rows,
            "specification" => specification,
            "known_upstream_defects" =>
                unique_sorted(entry["known_numerical_or_semantic_ambiguities"]),
            "implementation_files" => provenance_strings(mapped, "julia_file"),
            "test_files" => unique_sorted(test_files),
            "oracle_files" => unique_sorted(oracle_files),
            "documentation_files" => provenance_paths(mapped, "docs"),
            "benchmark_files" => benchmark_files,
            "verification_status" => String(entry["test_status"]),
            "provenance_verification_statuses" =>
                provenance_strings(mapped, "verification_status"),
            "mathematical_category" => String(entry["mathematical_category"]),
            "documentation_status" => String(entry["documentation_status"]),
            "benchmark_status" => String(entry["benchmark_status"]),
            "review_notes" => unique_sorted(entry["review_notes"]),
            "missing_declared_dependencies" =>
                unique_sorted(entry["unresolved_declared_dependencies"]),
        )
        push!(public_rows, row)
    end

    incomplete_set = Set(incomplete_names)
    helper_rows = Dict{String,Any}[]
    for entry in sort!(internal_inventory; by=row -> lowercase(row["function_name"]))
        helper_name = String(entry["function_name"])
        status = String(entry["implementation_status"])
        disposition = helper_status_class(status)
        direct_public_consumers = unique_sorted(
            row["function_name"] for
            row in public_inventory if helper_name in row["upstream_calls"]
        )
        completion_consumers = unique_sorted(
            name for
            name in incomplete_names if helper_name in dependency_closure(name, graph)
        )
        packages = String[
            policy.assignments[name] for
            name in completion_consumers if haskey(policy.assignments, name)
        ]
        if is_required_helper(disposition)
            push!(packages, policy.helper_assignments[helper_name])
        end
        packages = unique_sorted(packages)
        push!(
            helper_rows,
            Dict{String,Any}(
                "function_name" => helper_name,
                "source_path" => String(entry["source_path"]),
                "current_status" => status,
                "disposition_class" => disposition,
                "required_for_completion" => is_required_helper(disposition),
                "direct_public_consumers" => direct_public_consumers,
                "completion_consumers" => completion_consumers,
                "work_packages" => packages,
                "review_notes" => unique_sorted(entry["review_notes"]),
            ),
        )
    end

    queue_rows = Dict{String,Any}[]
    for row in public_rows
        name = row["upstream_function"]
        name in incomplete_set || continue
        inventory_entry = inventory_by_name[name]
        incomplete_dependencies = unique_sorted(
            dependency for
            dependency in inventory_entry["upstream_calls"] if dependency in incomplete_set
        )
        internal_dependencies = unique_sorted(
            dependency for dependency in inventory_entry["upstream_calls"] if
            dependency in required_internal_names
        )
        package = package_by_id[row["work_package"]]
        push!(
            queue_rows,
            Dict{String,Any}(
                "upstream_function" => name,
                "source_path" => row["source_path"],
                "current_status" => row["current_status"],
                "target_status" => row["target_status"],
                "work_package" => row["work_package"],
                "work_package_order" => package.order,
                "call_graph_rank" => ranks[name],
                "capability_prerequisites" => package.prerequisite_capabilities,
                "incomplete_public_call_dependencies" => incomplete_dependencies,
                "unresolved_internal_call_dependencies" => internal_dependencies,
                "missing_declared_dependencies" => row["missing_declared_dependencies"],
                "required_backend" => row["required_backend"],
                "reason" => row["review_notes"],
            ),
        )
    end
    sort!(
        queue_rows;
        by=row -> (
            row["work_package_order"],
            row["call_graph_rank"],
            lowercase(row["upstream_function"]),
        ),
    )
    for (index, row) in pairs(queue_rows)
        row["queue_index"] = index
    end

    counts = Dict(
        class =>
            count(row -> public_status_class(row["current_status"]) == class, public_rows)
        for class in ("implementation_complete", "partial", "deferred", "blocked")
    )
    helper_counts = Dict(
        class => count(row -> row["disposition_class"] == class, helper_rows) for
        class in unique(values(HELPER_STATUS_CLASS))
    )
    hashes = (
        inventory=file_sha256(INVENTORY_PATH),
        status=String(inventory["status_overlay_sha256"]),
        provenance=file_sha256(PROVENANCE_PATH),
        policy=file_sha256(POLICY_PATH),
    )
    metadata = (
        source_revision=String(inventory["source_revision"]),
        audit_date=String(inventory["audit_date"]),
        dependency_edge_count=Int(inventory["dependency_edge_count"]),
        dependency_cycle_count=Int(inventory["dependency_cycle_count"]),
        public_row_count=length(public_rows),
        internal_helper_count=length(helper_rows),
        incomplete_public_count=length(queue_rows),
    )
    model = (;
        metadata,
        hashes,
        policy,
        public_rows,
        helper_rows,
        queue_rows,
        counts,
        helper_counts,
    )
    validate_model(model)
    return model
end

function validate_model(model)
    length(model.public_rows) == model.metadata.public_row_count ||
        error("completion model public-row count mismatch")
    names = String[row["upstream_function"] for row in model.public_rows]
    length(names) == length(unique(names)) || error("completion model repeats a public row")
    for row in model.public_rows
        missing = setdiff(REQUIRED_PUBLIC_ROW_FIELDS, Set(keys(row)))
        isempty(missing) || error(
            "completion row $(row["upstream_function"]) is missing required fields: " *
            join(sort!(collect(missing)), ", "),
        )
    end
    queue_names = Set(String[row["upstream_function"] for row in model.queue_rows])
    incomplete_names = Set(
        String[
            row["upstream_function"] for row in model.public_rows if
            !(row["current_status"] in IMPLEMENTATION_COMPLETE_STATUSES)
        ],
    )
    queue_names == incomplete_names ||
        error("completion queue does not exactly match incomplete public rows")
    length(model.helper_rows) == model.metadata.internal_helper_count ||
        error("completion model internal-helper count mismatch")
    return model
end

function audit_failure!(failures, category::AbstractString, message::AbstractString)
    push!(failures, "[$category] $message")
    return failures
end

function public_names(mod::Module)
    return Set(
        String(name) for
        name in names(mod; all=false, imported=false) if name != nameof(mod)
    )
end

function runtime_public_modules(repository_root::AbstractString=REPOSITORY_ROOT)
    isdefined(Main, :QuantumEntanglementTools) || error(
        "load QuantumEntanglementTools before running the static release-evidence audit"
    )
    package = getfield(Main, :QuantumEntanglementTools)
    package isa Module || error("QuantumEntanglementTools did not load as a module")
    isdefined(package, :MATLABCompat) ||
        error("QuantumEntanglementTools.MATLABCompat is not defined")
    compatibility = getfield(package, :MATLABCompat)
    compatibility isa Module ||
        error("QuantumEntanglementTools.MATLABCompat is not a module")
    return Dict(
        "QuantumEntanglementTools" => package,
        "QuantumEntanglementTools.MATLABCompat" => compatibility,
    )
end

function resolve_qualified_binding(path::AbstractString, modules)
    parts = split(String(path), '.')
    isempty(parts) && return nothing
    first(parts) == "QuantumEntanglementTools" || return nothing
    value = get(modules, "QuantumEntanglementTools", nothing)
    isnothing(value) && return nothing
    for part in parts[2:end]
        value isa Module || return nothing
        symbol = Symbol(part)
        isdefined(value, symbol) || return nothing
        value = getfield(value, symbol)
    end
    return value
end

function binding_kind(value)
    value isa Module && return "module"
    (value isa DataType || value isa UnionAll) && return "type"
    value isa Function && return "function"
    value isa Enum && return "constant"
    return "other"
end

function kind_matches(declared::AbstractString, actual::AbstractString)
    declared == actual && return true
    declared in ("alias", "compatibility_wrapper") && actual == "function" && return true
    return false
end

function existing_repository_reference(
    repository_root::AbstractString, reference::AbstractString
)
    startswith(reference, "http://") && return true
    startswith(reference, "https://") && return true
    isabspath(reference) && return false
    return isfile(joinpath(repository_root, reference))
end

function method_accepts_explicit_rng(value)
    value isa Function || return false
    return any(method -> occursin("AbstractRNG", sprint(show, method.sig)), methods(value))
end

function source_paths(repository_root::AbstractString)
    paths = String[]
    for root_name in ("src", "ext")
        root = joinpath(repository_root, root_name)
        isdir(root) || continue
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
    return call_name(signature.args[1])
end

function normalized_signature(signature)
    return replace(sprint(show, signature), r"\s+" => " ")
end

function scan_source_expression!(
    expression, state, relative_path; line::Int=1, exported_names::Set{String}=Set{String}()
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
    end

    for argument in expression.args
        current_line = scan_source_expression!(
            argument, state, relative_path; line=current_line, exported_names=exported_names
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
    references::Set{Tuple{Module,Symbol}},
    function_value,
    seen::Base.IdSet{Any},
    package::Module,
    compatibility::Module,
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
        (owner === package || owner === compatibility || startswith(String(name), "#")) ||
            continue
        isdefined(owner, name) || continue
        nested = getfield(owner, name)
        nested isa Function || continue
        function_global_refs!(references, nested, seen, package, compatibility)
    end
    return references
end

function compatibility_wrapper_bypass_audit(provenance_entries, modules)
    candidates = Dict{String,Any}[]
    documented_independent = Dict{String,Any}[]
    package = modules["QuantumEntanglementTools"]
    compatibility = modules["QuantumEntanglementTools.MATLABCompat"]
    for entry in provenance_entries
        get(entry, "public_module", "") == "QuantumEntanglementTools.MATLABCompat" ||
            continue
        get(entry, "public_kind", "") in ("alias", "compatibility_wrapper") || continue
        delegates = get(entry, "delegates_to", Any[])
        delegates isa AbstractVector && !isempty(delegates) || continue
        name = string(get(entry, "julia_name", ""))
        symbol = Symbol(name)
        isdefined(compatibility, symbol) || continue
        wrapper = getfield(compatibility, symbol)

        resolved = filter(
            value -> !isnothing(value),
            resolve_qualified_binding.(String.(delegates), Ref(modules)),
        )
        any(delegate -> wrapper === delegate, resolved) && continue

        references = Set{Tuple{Module,Symbol}}()
        if wrapper isa Function
            function_global_refs!(
                references, wrapper, Base.IdSet{Any}(), package, compatibility
            )
        end
        expected = Set{Tuple{Module,Symbol}}()
        for delegate_path in String.(delegates)
            parts = split(delegate_path, '.')
            length(parts) >= 2 || continue
            owner = if length(parts) == 2
                package
            else
                resolve_qualified_binding(join(parts[1:(end - 1)], '.'), modules)
            end
            owner isa Module || continue
            delegate_symbol = Symbol(last(parts))
            push!(expected, (owner, delegate_symbol))
            push!(expected, (compatibility, delegate_symbol))
        end
        isempty(intersect(references, expected)) || continue
        finding = Dict{String,Any}(
            "wrapper" => name,
            "declared_delegates" => sort!(String.(delegates)),
            "reason" => "no declared delegate appears in the transitive method AST",
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
    sort!(candidates; by=entry -> entry["wrapper"])
    sort!(documented_independent; by=entry -> entry["wrapper"])
    return (; candidates, documented_independent)
end

function scan_repository_sources(
    repository_root::AbstractString, exported_names_by_module=Dict{String,Set{String}}()
)
    state = (
        definitions=Dict{String,Vector{Dict{String,Any}}}(),
        global_rng_candidates=Dict{String,Any}[],
        seed_calls=Dict{String,Any}[],
    )
    markers = Dict{String,Any}[]
    paths = source_paths(repository_root)
    native_names = get(exported_names_by_module, "QuantumEntanglementTools", Set{String}())
    compatibility_names = get(
        exported_names_by_module, "QuantumEntanglementTools.MATLABCompat", Set{String}()
    )

    for path in paths
        relative = normalized_path(relpath(path, repository_root))
        content = read(path, String)
        lines = split(content, '\n'; keepempty=true)
        exported_names = if startswith(relative, "src/compat/")
            compatibility_names
        else
            native_names
        end
        expression = Meta.parseall(content)
        scan_source_expression!(expression, state, relative; exported_names=exported_names)
        for (line_number, text) in pairs(lines)
            occursin(SOURCE_PLACEHOLDER_PATTERN, text) || continue
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
    sort!(markers; by=entry -> (entry["path"], entry["line"]))
    return (;
        source_file_count=length(paths),
        markers,
        duplicate_definitions,
        global_rng_candidates=state.global_rng_candidates,
        seed_calls=state.seed_calls,
    )
end

function combined_file_content(repository_root::AbstractString, paths)
    io = IOBuffer()
    for reference in unique_sorted(paths)
        startswith(reference, "http://") && continue
        startswith(reference, "https://") && continue
        isabspath(reference) && continue
        path = joinpath(repository_root, reference)
        isfile(path) || continue
        write(io, read(path))
        write(io, '\n')
    end
    return String(take!(io))
end

function extension_environment_for_test(
    repository_root::AbstractString, test_path::AbstractString
)
    parts = split(normalized_path(test_path), '/')
    length(parts) >= 4 || return nothing
    parts[1:2] == ["test", "extensions"] || return nothing
    environment = joinpath(repository_root, parts[1], parts[2], parts[3])
    isfile(joinpath(environment, "Project.toml")) || return nothing
    return normalized_path(relpath(environment, repository_root))
end

function optimization_evidence_failures(
    row, mapped_entries, repository_root::AbstractString
)
    failures = String[]
    name = String(row["upstream_function"])
    name in SOLVER_BACKED_PUBLIC_ROWS || return failures
    row["current_status"] in IMPLEMENTATION_COMPLETE_STATUSES || return failures

    test_paths = provenance_paths(mapped_entries, "tests")
    extension_tests = filter(path -> startswith(path, "test/extensions/"), test_paths)
    isempty(extension_tests) && push!(
        failures,
        "implemented solver-backed row $name has no dedicated optional-extension test",
    )
    environments = filter(
        environment -> !isnothing(environment),
        extension_environment_for_test.(Ref(repository_root), extension_tests),
    )
    isempty(environments) && push!(
        failures,
        "implemented solver-backed row $name has no test extension environment with Project.toml",
    )

    bindings = String[]
    for entry in mapped_entries
        push!(
            bindings,
            string(get(entry, "public_module", ""), ".", get(entry, "julia_name", "")),
        )
    end
    simple_names = unique_sorted(last(split(binding, '.')) for binding in bindings)
    extension_match = false
    extension_paths = filter(
        path -> occursin("/ext/", "/" * normalized_path(path)),
        source_paths(repository_root),
    )
    for path in extension_paths
        content = read(path, String)
        any(name -> occursin(name, content), simple_names) || continue
        extension_match = true
        break
    end
    !extension_match && push!(
        failures, "implemented solver-backed row $name has no matching method in ext/"
    )

    tests = lowercase(combined_file_content(repository_root, extension_tests))
    status_tokens = (
        "termination_status", "primal_status", "dual_status", "criterionunknown", ":unknown"
    )
    any(token -> occursin(token, tests), status_tokens) || push!(
        failures,
        "implemented solver-backed row $name lacks extension tests for structured solver status",
    )
    failure_tokens = (
        "@test_throws", "infeasible", "time_limit", "failure", "failed", "unknown"
    )
    any(token -> occursin(token, tests), failure_tokens) || push!(
        failures,
        "implemented solver-backed row $name lacks extension tests for solver failure/inconclusive behavior",
    )
    architecture_paths = (
        joinpath(repository_root, "docs", "OPTIMIZATION_ARCHITECTURE.md"),
        joinpath(repository_root, "docs", "src", "optimization_architecture.md"),
    )
    any(isfile, architecture_paths) || push!(
        failures,
        "implemented solver-backed row $name requires a documented optimization architecture",
    )
    return failures
end

function repository_quality_audit(
    model=build_completion_model();
    repository_root::AbstractString=REPOSITORY_ROOT,
    provenance_entries=nothing,
    runtime_modules=nothing,
    scan_sources::Bool=true,
)
    failures = String[]
    limitations = [
        "Static AST and token checks identify evidence and bypass candidates; they do not prove semantic correctness or MATLAB parity.",
        "Test files and verification labels are checked for required evidence, but this lightweight checker does not execute package, solver, oracle, extension, or documentation suites.",
        "Optimization checks require package-owned status/failure evidence and a dedicated optional-extension environment; they do not validate a solver certificate.",
        "Absence of a placeholder or global-RNG candidate is not a proof that every runtime branch is reachable or deterministic.",
    ]

    provenance = if isnothing(provenance_entries)
        path = joinpath(repository_root, "PROVENANCE.toml")
        isfile(path) || error("required provenance ledger does not exist: $path")
        parsed = TOML.parsefile(path)
        get(parsed, "schema_version", nothing) == 1 || audit_failure!(
            failures, "provenance", "PROVENANCE.toml schema_version is not 1"
        )
        get(parsed, "functions", Any[])
    else
        provenance_entries
    end
    provenance isa AbstractVector ||
        error("PROVENANCE.toml `functions` must be an array of tables")

    modules = if isnothing(runtime_modules)
        try
            runtime_public_modules(repository_root)
        catch exception
            audit_failure!(
                failures,
                "exports",
                "could not load the public modules: $(sprint(showerror, exception))",
            )
            Dict{String,Module}()
        end
    else
        runtime_modules
    end
    runtime_available = all(
        module_name -> haskey(modules, module_name), PUBLIC_MODULE_NAMES
    )
    exported_names_by_module = Dict{String,Set{String}}()
    exported_keys = Set{Tuple{String,String}}()
    for module_name in PUBLIC_MODULE_NAMES
        mod = get(modules, module_name, nothing)
        mod isa Module || continue
        module_names = public_names(mod)
        exported_names_by_module[module_name] = module_names
        union!(exported_keys, ((module_name, name) for name in module_names))
    end

    inventory_path = joinpath(repository_root, "porting", "qetlab_inventory.toml")
    inventory_entries = if isfile(inventory_path)
        get(TOML.parsefile(inventory_path), "functions", Any[])
    else
        Any[]
    end
    inventory_by_name = Dict{String,Any}()
    for entry in inventory_entries
        name = string(get(entry, "function_name", ""))
        isempty(name) || (inventory_by_name[name] = entry)
    end

    provenance_by_key = Dict{Tuple{String,String},Any}()
    qetlab_by_upstream = Dict{String,Vector{Any}}()
    for (entry_number, entry) in pairs(provenance)
        entry isa AbstractDict || begin
            audit_failure!(failures, "provenance", "entry $entry_number is not a table")
            continue
        end
        missing = setdiff(REQUIRED_PROVENANCE_FIELDS, Set(keys(entry)))
        isempty(missing) || audit_failure!(
            failures,
            "provenance",
            "entry $entry_number is missing fields: $(join(sort!(collect(missing)), ", "))",
        )
        name = get(entry, "julia_name", nothing)
        module_name = get(entry, "public_module", nothing)
        name isa AbstractString && module_name isa AbstractString || begin
            audit_failure!(
                failures,
                "provenance",
                "entry $entry_number must have string julia_name and public_module",
            )
            continue
        end
        qualified = "$module_name.$name"
        key = (String(module_name), String(name))
        if haskey(provenance_by_key, key)
            audit_failure!(failures, "provenance", "duplicate binding entry for $qualified")
        else
            provenance_by_key[key] = entry
        end

        module_name in PUBLIC_MODULE_NAMES ||
            audit_failure!(failures, "exports", "unknown public module for $qualified")
        mod = get(modules, String(module_name), nothing)
        if mod isa Module
            symbol = Symbol(name)
            if !isdefined(mod, symbol)
                audit_failure!(
                    failures, "exports", "provenance binding is undefined: $qualified"
                )
            else
                declared_kind = string(get(entry, "public_kind", ""))
                actual_kind = binding_kind(getfield(mod, symbol))
                declared_kind in ALLOWED_PUBLIC_KINDS || audit_failure!(
                    failures,
                    "provenance",
                    "$qualified has unsupported public_kind '$declared_kind'",
                )
                kind_matches(declared_kind, actual_kind) || audit_failure!(
                    failures,
                    "exports",
                    "$qualified declares '$declared_kind' but runtime kind is '$actual_kind'",
                )
            end
        end

        status = string(get(entry, "status", ""))
        status in ALLOWED_PROVENANCE_STATUSES || audit_failure!(
            failures, "provenance", "$qualified has unsupported status '$status'"
        )
        verification = lowercase(string(get(entry, "verification_status", "")))
        occursin("pass", verification) || audit_failure!(
            failures,
            "implemented-claim",
            "$qualified has no passing verification_status evidence",
        )
        for field in ("julia_file", "specification")
            value = get(entry, field, nothing)
            value isa AbstractString && !isempty(strip(value)) || begin
                audit_failure!(
                    failures, "files", "$qualified has no non-empty $field reference"
                )
                continue
            end
            existing_repository_reference(repository_root, value) || audit_failure!(
                failures, "files", "$qualified references missing $field path: $value"
            )
        end
        for field in ("tests", "docs")
            values = get(entry, field, nothing)
            values isa AbstractVector &&
            !isempty(values) &&
            all(value -> value isa AbstractString, values) || begin
                audit_failure!(
                    failures,
                    "files",
                    "$qualified must provide a non-empty array of $field paths",
                )
                continue
            end
            for value in values
                existing_repository_reference(repository_root, value) || audit_failure!(
                    failures,
                    "files",
                    "$qualified references missing $field path: $value",
                )
            end
        end

        declared_kind = string(get(entry, "public_kind", ""))
        if declared_kind in ("alias", "compatibility_wrapper") && runtime_available
            delegates = get(entry, "delegates_to", nothing)
            delegates isa AbstractVector &&
            !isempty(delegates) &&
            all(delegate -> delegate isa AbstractString, delegates) || begin
                audit_failure!(
                    failures,
                    "compatibility",
                    "$qualified must declare a non-empty delegates_to array",
                )
                delegates = String[]
            end
            for delegate in delegates
                isnothing(resolve_qualified_binding(delegate, modules)) && audit_failure!(
                    failures,
                    "compatibility",
                    "$qualified declares unresolved delegate $delegate",
                )
            end
        end

        source_project = string(get(entry, "source_project", ""))
        if source_project == "QETLAB"
            missing_qetlab = setdiff(REQUIRED_QETLAB_PROVENANCE_FIELDS, Set(keys(entry)))
            isempty(missing_qetlab) || audit_failure!(
                failures,
                "provenance",
                "$qualified QETLAB evidence is missing: $(join(sort!(collect(missing_qetlab)), ", "))",
            )
            upstream = string(get(entry, "upstream_function", ""))
            isempty(upstream) || push!(get!(qetlab_by_upstream, upstream, Any[]), entry)
            inventory_entry = get(inventory_by_name, upstream, nothing)
            if isnothing(inventory_entry)
                audit_failure!(
                    failures,
                    "provenance",
                    "$qualified names absent QETLAB inventory row $upstream",
                )
            else
                for field in ("source_path", "source_revision", "source_sha256")
                    get(entry, field, nothing) == get(inventory_entry, field, nothing) ||
                        audit_failure!(
                            failures,
                            "provenance",
                            "$qualified $field disagrees with inventory row $upstream",
                        )
                end
            end
        elseif source_project == "QuantumEntanglementTools"
            origin = get(entry, "project_origin", nothing)
            origin isa AbstractString && !isempty(strip(origin)) || audit_failure!(
                failures,
                "provenance",
                "$qualified project-native evidence has no project_origin",
            )
        else
            audit_failure!(
                failures,
                "provenance",
                "$qualified has unsupported source_project '$source_project'",
            )
        end
    end

    provenance_keys = Set(keys(provenance_by_key))
    if runtime_available
        for (module_name, name) in sort!(collect(setdiff(exported_keys, provenance_keys)))
            audit_failure!(
                failures,
                "exports",
                "exported binding has no provenance: $module_name.$name",
            )
        end
        for (module_name, name) in sort!(collect(setdiff(provenance_keys, exported_keys)))
            audit_failure!(
                failures, "exports", "provenance entry is not exported: $module_name.$name"
            )
        end
    end

    complete_rows = [
        row for row in model.public_rows if
        row["current_status"] in IMPLEMENTATION_COMPLETE_STATUSES
    ]
    randomized_complete_rows = Dict{String,Any}[]
    solver_complete_rows = Dict{String,Any}[]
    for row in complete_rows
        upstream = String(row["upstream_function"])
        mapped = get(qetlab_by_upstream, upstream, Any[])
        isempty(mapped) && audit_failure!(
            failures,
            "implemented-claim",
            "implementation-complete row $upstream has no direct QETLAB provenance",
        )
        native = [
            entry for
            entry in mapped if get(entry, "public_module", "") == "QuantumEntanglementTools"
        ]
        compatibility = [
            entry for entry in mapped if
            get(entry, "public_module", "") == "QuantumEntanglementTools.MATLABCompat"
        ]
        isempty(native) && audit_failure!(
            failures,
            "implemented-claim",
            "implementation-complete row $upstream has no native binding evidence",
        )
        isempty(compatibility) && audit_failure!(
            failures,
            "implemented-claim",
            "implementation-complete row $upstream has no MATLABCompat binding evidence",
        )
        for field in (
            "native_binding",
            "compatibility_binding",
            "implementation_files",
            "test_files",
            "documentation_files",
            "specification",
        )
            isempty(row[field]) && audit_failure!(
                failures,
                "implemented-claim",
                "implementation-complete row $upstream has no $field evidence",
            )
        end
        occursin("pass", lowercase(String(row["verification_status"]))) || audit_failure!(
            failures,
            "implemented-claim",
            "implementation-complete row $upstream lacks a passing inventory test disposition",
        )
        for entry in mapped
            verification = lowercase(string(get(entry, "verification_status", "")))
            occursin("pass", verification) || audit_failure!(
                failures,
                "implemented-claim",
                "$(qualified_binding(entry)) lacks a passing verification status for $upstream",
            )
        end

        inventory_entry = get(inventory_by_name, upstream, nothing)
        if !isnothing(inventory_entry)
            proposed = string(get(inventory_entry, "proposed_julia_name", ""))
            compatibility_alias = string(get(inventory_entry, "compatibility_alias", ""))
            native_key = ("QuantumEntanglementTools", proposed)
            compatibility_key = (
                "QuantumEntanglementTools.MATLABCompat", compatibility_alias
            )
            isempty(proposed) ||
                native_key in provenance_keys ||
                audit_failure!(
                    failures,
                    "implemented-claim",
                    "$upstream proposed native binding lacks provenance: $(join(native_key, '.'))",
                )
            isempty(compatibility_alias) ||
                compatibility_key in provenance_keys ||
                audit_failure!(
                    failures,
                    "implemented-claim",
                    "$upstream compatibility binding lacks provenance: $(join(compatibility_key, '.'))",
                )

            random_behavior = string(get(inventory_entry, "random_behavior", ""))
            if !startswith(random_behavior, "deterministic_")
                push!(randomized_complete_rows, row)
                binding_values = filter(
                    value -> !isnothing(value),
                    resolve_qualified_binding.(
                        vcat(row["native_binding"], row["compatibility_binding"]),
                        Ref(modules),
                    ),
                )
                if runtime_available
                    any(method_accepts_explicit_rng, binding_values) || audit_failure!(
                        failures,
                        "rng",
                        "implemented randomized row $upstream exposes no AbstractRNG method",
                    )
                end
                test_paths = provenance_paths(mapped, "tests")
                test_content = combined_file_content(repository_root, test_paths)
                binding_names = unique_sorted(
                    string(get(entry, "julia_name", "")) for entry in mapped
                )
                explicit_rng_tokens = (
                    "Xoshiro(", "MersenneTwister(", "AbstractRNG", "default_rng()"
                )
                (
                    any(token -> occursin(token, test_content), explicit_rng_tokens) &&
                    any(name -> occursin(name, test_content), binding_names)
                ) || audit_failure!(
                    failures,
                    "rng",
                    "implemented randomized row $upstream lacks an explicit-RNG binding test",
                )
                global_isolation_tokens = ("Random.seed!", "copy(Random.default_rng())")
                any(token -> occursin(token, test_content), global_isolation_tokens) ||
                    audit_failure!(
                        failures,
                        "rng",
                        "implemented randomized row $upstream lacks a global-RNG isolation test",
                    )
            end
        end

        if upstream in SOLVER_BACKED_PUBLIC_ROWS
            push!(solver_complete_rows, row)
            for message in optimization_evidence_failures(row, mapped, repository_root)
                audit_failure!(failures, "optimization", message)
            end
        end
    end

    source_scan = if scan_sources
        scan_repository_sources(repository_root, exported_names_by_module)
    else
        (;
            source_file_count=0,
            markers=Dict{String,Any}[],
            duplicate_definitions=Dict{String,Any}[],
            global_rng_candidates=Dict{String,Any}[],
            seed_calls=Dict{String,Any}[],
        )
    end
    for marker in source_scan.markers
        audit_failure!(
            failures,
            "placeholder",
            "$(marker["path"]):$(marker["line"]) contains $(repr(marker["text"]))",
        )
    end
    for finding in source_scan.duplicate_definitions
        locations = join(
            ("$(site["path"]):$(site["line"])" for site in finding["sites"]), ", "
        )
        audit_failure!(
            failures,
            "duplicate-definition",
            "$(finding["signature"]) is defined at $locations",
        )
    end
    for finding in source_scan.global_rng_candidates
        audit_failure!(
            failures,
            "rng",
            "$(finding["path"]):$(finding["line"]) calls $(finding["call"]) without a lexical RNG argument",
        )
    end
    for finding in source_scan.seed_calls
        audit_failure!(
            failures,
            "rng",
            "$(finding["path"]):$(finding["line"]) calls Random.seed!/seed! in library code",
        )
    end

    bypass = if runtime_available
        compatibility_wrapper_bypass_audit(provenance, modules)
    else
        (; candidates=Dict{String,Any}[], documented_independent=Dict{String,Any}[])
    end
    for finding in bypass.candidates
        audit_failure!(
            failures,
            "compatibility",
            "MATLABCompat.$(finding["wrapper"]) bypasses declared delegates $(join(finding["declared_delegates"], ", "))",
        )
    end

    failures = sort!(unique!(failures); by=lowercase)
    metrics = Dict{String,Int}(
        "runtime_export_count" => length(exported_keys),
        "provenance_entry_count" => length(provenance),
        "implementation_complete_rows_checked" => length(complete_rows),
        "randomized_complete_rows_checked" => length(randomized_complete_rows),
        "solver_backed_complete_rows_checked" => length(solver_complete_rows),
        "source_files_scanned" => source_scan.source_file_count,
        "placeholder_markers" => length(source_scan.markers),
        "duplicate_public_definition_candidates" =>
            length(source_scan.duplicate_definitions),
        "global_rng_candidates" => length(source_scan.global_rng_candidates),
        "library_seed_calls" => length(source_scan.seed_calls),
        "compatibility_wrapper_bypass_candidates" => length(bypass.candidates),
        "documented_independent_compatibility_implementations" =>
            length(bypass.documented_independent),
        "failure_count" => length(failures),
    )
    return (; failures, metrics, limitations, source_scan, compatibility_bypass=bypass)
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

toml_value(value::AbstractString) = "\"" * toml_escape(value) * "\""
toml_value(value::Bool) = value ? "true" : "false"
toml_value(value::Integer) = string(value)
function toml_value(values::AbstractVector)
    return "[" * join((toml_value(value) for value in values), ", ") * "]"
end

function write_toml_field(io::IO, key::AbstractString, value)
    return println(io, key, " = ", toml_value(value))
end

function write_source_metadata(io::IO, model)
    write_toml_field(io, "source_revision", model.metadata.source_revision)
    write_toml_field(io, "audit_date", model.metadata.audit_date)
    write_toml_field(io, "inventory_sha256", model.hashes.inventory)
    write_toml_field(io, "status_overlay_sha256", model.hashes.status)
    write_toml_field(io, "provenance_sha256", model.hashes.provenance)
    write_toml_field(io, "completion_policy_sha256", model.hashes.policy)
    write_toml_field(io, "dependency_edge_count", model.metadata.dependency_edge_count)
    write_toml_field(io, "dependency_cycle_count", model.metadata.dependency_cycle_count)
    return write_toml_field(io, "source_call_graph_is_implementation_dag", false)
end

function render_completion_plan(model)
    io = IOBuffer()
    println(io, "# Generated by $GENERATOR_PATH; do not edit by hand.")
    println(io, "schema_version = 1")
    write_toml_field(io, "generated_by", GENERATOR_PATH)
    write_source_metadata(io, model)
    write_toml_field(io, "public_row_count", model.metadata.public_row_count)
    write_toml_field(
        io, "implementation_complete_count", model.counts["implementation_complete"]
    )
    write_toml_field(io, "partial_count", model.counts["partial"])
    write_toml_field(io, "deferred_count", model.counts["deferred"])
    write_toml_field(io, "blocked_count", model.counts["blocked"])
    write_toml_field(io, "incomplete_public_count", model.metadata.incomplete_public_count)
    write_toml_field(io, "internal_helper_count", model.metadata.internal_helper_count)

    for package in model.policy.packages
        println(io)
        println(io, "[[work_packages]]")
        write_toml_field(io, "id", package.id)
        write_toml_field(io, "order", package.order)
        write_toml_field(io, "title", package.title)
        write_toml_field(io, "summary", package.summary)
        write_toml_field(io, "architecture_only", package.architecture_only)
        write_toml_field(io, "prerequisite_capabilities", package.prerequisite_capabilities)
        write_toml_field(io, "provides_capabilities", package.provides_capabilities)
        write_toml_field(io, "functions", package.functions)
        write_toml_field(io, "internal_helpers", package.internal_helpers)
    end

    public_field_order = [
        "upstream_function",
        "source_path",
        "current_status",
        "target_status",
        "work_package",
        "native_binding",
        "compatibility_binding",
        "required_backend",
        "prerequisite_rows",
        "specification",
        "known_upstream_defects",
        "implementation_files",
        "test_files",
        "oracle_files",
        "documentation_files",
        "benchmark_files",
        "verification_status",
        "provenance_verification_statuses",
        "mathematical_category",
        "documentation_status",
        "benchmark_status",
        "review_notes",
        "missing_declared_dependencies",
    ]
    for row in model.public_rows
        println(io)
        println(io, "[[public_rows]]")
        for field in public_field_order
            write_toml_field(io, field, row[field])
        end
    end

    helper_field_order = [
        "function_name",
        "source_path",
        "current_status",
        "disposition_class",
        "required_for_completion",
        "direct_public_consumers",
        "completion_consumers",
        "work_packages",
        "review_notes",
    ]
    for row in model.helper_rows
        println(io)
        println(io, "[[internal_helpers]]")
        for field in helper_field_order
            write_toml_field(io, field, row[field])
        end
    end
    return String(take!(io))
end

function render_completion_queue(model)
    io = IOBuffer()
    println(io, "# Generated by $GENERATOR_PATH; do not edit by hand.")
    println(io, "schema_version = 1")
    write_toml_field(io, "generated_by", GENERATOR_PATH)
    write_source_metadata(io, model)
    write_toml_field(io, "incomplete_public_count", model.metadata.incomplete_public_count)
    write_toml_field(io, "work_package_count", length(model.policy.packages))
    write_toml_field(
        io,
        "ordering",
        "capability-led work-package order, then source-call rank, then function name",
    )

    field_order = [
        "queue_index",
        "upstream_function",
        "source_path",
        "current_status",
        "target_status",
        "work_package",
        "work_package_order",
        "call_graph_rank",
        "capability_prerequisites",
        "incomplete_public_call_dependencies",
        "unresolved_internal_call_dependencies",
        "missing_declared_dependencies",
        "required_backend",
        "reason",
    ]
    for row in model.queue_rows
        println(io)
        println(io, "[[tasks]]")
        for field in field_order
            write_toml_field(io, field, row[field])
        end
    end
    return String(take!(io))
end

function markdown_escape(value::AbstractString)
    return replace(replace(String(value), "|" => "\\|"), "\n" => " ")
end

function markdown_code_list(values)
    isempty(values) && return "—"
    return join(("`$(markdown_escape(value))`" for value in values), ", ")
end

function render_completion_document(model)
    io = IOBuffer()
    println(io, "<!-- Generated by $GENERATOR_PATH; do not edit by hand. -->")
    println(io, "# QETLAB completion plan")
    println(io)
    println(
        io,
        "This plan is generated from the reviewed status overlay, the generated ",
        "QETLAB inventory, public provenance, and the small capability policy in ",
        "`porting/qetlab_completion_policy.toml`.",
    )
    println(io)
    println(io, "## Current snapshot")
    println(io)
    println(io, "| Measure | Rows |")
    println(io, "|---|---:|")
    println(io, "| Public inventory | $(model.metadata.public_row_count) |")
    println(io, "| Implementation-complete | $(model.counts["implementation_complete"]) |")
    println(io, "| Partial | $(model.counts["partial"]) |")
    println(io, "| Deferred | $(model.counts["deferred"]) |")
    println(io, "| Blocked with an explicit reason | $(model.counts["blocked"]) |")
    println(io, "| Incomplete public queue | $(model.metadata.incomplete_public_count) |")
    println(io, "| Internal helpers | $(model.metadata.internal_helper_count) |")
    println(io)
    println(io, "- QETLAB revision: `$(model.metadata.source_revision)`")
    println(io, "- Inventory SHA-256: `$(model.hashes.inventory)`")
    println(io, "- Status-overlay SHA-256: `$(model.hashes.status)`")
    println(io, "- Provenance SHA-256: `$(model.hashes.provenance)`")
    println(io, "- Completion-policy SHA-256: `$(model.hashes.policy)`")
    println(io)
    println(
        io,
        "The all-public-row machine ledger is ",
        "[`porting/qetlab_completion_plan.toml`](../porting/qetlab_completion_plan.toml). ",
        "The executable queue is ",
        "[`porting/qetlab_completion_queue.toml`](../porting/qetlab_completion_queue.toml).",
    )
    println(io)
    println(io, "## Ordering semantics")
    println(io)
    println(
        io,
        "The queue is ordered by reviewed capability package, then by rank in the ",
        "static upstream source-call graph, then by function name. The 503-edge ",
        "graph is evidence about calls in the pinned MATLAB source. It is not an ",
        "implementation DAG: it cannot express package-owned architecture, and a ",
        "call to a partially implemented row may use a branch that is already ",
        "available. Capability prerequisites therefore control package order; ",
        "the source-call dependencies remain visible on every task for review.",
    )

    queue_by_package = Dict(
        package.id =>
            [row for row in model.queue_rows if row["work_package"] == package.id] for
        package in model.policy.packages
    )
    for package in model.policy.packages
        println(io)
        println(io, "## $(package.order). $(package.title)")
        println(io)
        println(io, package.summary)
        println(io)
        println(io, "**Requires:** ", markdown_code_list(package.prerequisite_capabilities))
        println(io)
        println(io, "**Provides:** ", markdown_code_list(package.provides_capabilities))
        tasks = queue_by_package[package.id]
        if isempty(tasks)
            println(io)
            if !isempty(package.internal_helpers)
                println(
                    io,
                    "No public row is queued for this work package. Required private ",
                    "helper dispositions remain: ",
                    markdown_code_list(package.internal_helpers),
                    ".",
                )
            elseif package.architecture_only
                println(
                    io,
                    "This is an architecture gate. It intentionally owns no upstream ",
                    "row; later solver-backed packages depend on its result contracts.",
                )
            else
                println(
                    io,
                    "This work package is complete; no public row or required private ",
                    "helper remains queued.",
                )
            end
            continue
        end
        println(io)
        println(
            io, "| # | Upstream row | Current status | Source-call prerequisites | Reason |"
        )
        println(io, "|---:|---|---|---|---|")
        for row in tasks
            reason = if isempty(row["reason"])
                "No review note recorded."
            else
                join(row["reason"], " ")
            end
            dependencies = unique_sorted(
                vcat(
                    row["incomplete_public_call_dependencies"],
                    row["unresolved_internal_call_dependencies"],
                    row["missing_declared_dependencies"],
                ),
            )
            println(
                io,
                "| $(row["queue_index"]) | `$(row["upstream_function"])` | ",
                "`$(row["current_status"])` | ",
                "$(markdown_code_list(dependencies)) | ",
                "$(markdown_escape(reason)) |",
            )
        end
    end

    terminal_helpers = [row for row in model.helper_rows if !row["required_for_completion"]]
    required_helpers = [row for row in model.helper_rows if row["required_for_completion"]]
    println(io)
    println(io, "## Internal-helper dispositions")
    println(io)
    if isempty(required_helpers)
        println(
            io,
            "The 36 private helpers are not public parity promises. All ",
            "$(length(terminal_helpers)) have terminal superseded/excluded ",
            "dispositions; none remains required by an incomplete public ",
            "capability. Direct and transitive source-call consumers are recorded ",
            "where detected.",
        )
    else
        println(
            io,
            "The 36 private helpers are not public parity promises. ",
            "$(length(terminal_helpers)) have terminal superseded/excluded ",
            "dispositions; $(length(required_helpers)) retain required ",
            "partial/deferred dispositions and are assigned to completion work ",
            "packages. Direct and transitive source-call consumers are recorded ",
            "where detected.",
        )
    end
    println(io)
    println(
        io,
        "- Terminal: ",
        markdown_code_list(sort!([row["function_name"] for row in terminal_helpers])),
    )
    println(
        io,
        "- Required: ",
        markdown_code_list(sort!([row["function_name"] for row in required_helpers])),
    )

    println(io)
    println(io, "## Checker scope and limitations")
    println(io)
    println(io)
    if iszero(model.metadata.incomplete_public_count) && isempty(required_helpers)
        println(
            io,
            "- `julia --startup-file=no --project=. $CHECKER_PATH` validates the ",
            "completed static ledger and its repository evidence. That result is ",
            "an implementation-completeness claim, not a MATLAB-parity or release ",
            "approval claim.",
        )
    else
        println(
            io,
            "- `julia --startup-file=no --project=. $CHECKER_PATH` validates the ",
            "ledger plus static release evidence and reports the current state ",
            "without claiming completion.",
        )
    end
    if iszero(model.metadata.incomplete_public_count) && isempty(required_helpers)
        println(
            io,
            "- `julia --startup-file=no --project=. $CHECKER_PATH --strict` is an ",
            "implementation-completeness gate and passes for the generated ledger: ",
            "all $(model.metadata.public_row_count) public rows and all ",
            "$(model.metadata.internal_helper_count) private-helper dispositions are ",
            "terminal.",
        )
    else
        println(
            io,
            "- `julia --startup-file=no --project=. $CHECKER_PATH --strict` is an ",
            "implementation-completeness gate and currently fails on ",
            "$(model.metadata.incomplete_public_count) queued public rows and ",
            "$(length(required_helpers)) required private helpers.",
        )
    end
    println(
        io,
        "- The static evidence gate checks runtime exports against unique ",
        "provenance entries, referenced implementation/specification/test/doc ",
        "files, support for implementation-complete dispositions, duplicate ",
        "public definitions, compatibility delegation, placeholder markers, ",
        "library RNG calls, and explicit-RNG tests for completed randomized rows.",
    )
    println(
        io,
        "- Every implementation-complete solver-backed row must have a matching ",
        "package extension, a dedicated test environment with its own ",
        "`Project.toml`, documented architecture, and status plus ",
        "failure/inconclusive test evidence.",
    )
    println(
        io,
        "- `julia --startup-file=no --project=. test/qetlab_completion_checker.jl` ",
        "runs focused mutation fixtures proving that missing files/provenance, ",
        "RNG evidence, solver-extension evidence, duplicate definitions, ",
        "placeholders, and global-RNG calls are rejected.",
    )
    if iszero(model.metadata.incomplete_public_count) && isempty(required_helpers)
        println(
            io,
            "- This strict pass does not by itself establish MATLAB parity, semantic ",
            "correctness, solver certificate validity, supported-platform coverage, ",
            "or the maintainer's non-delegable release review.",
        )
    else
        println(
            io,
            "- A future strict pass would not by itself establish MATLAB parity, ",
            "semantic correctness, solver certificate validity, supported-platform ",
            "coverage, or the maintainer's non-delegable release review.",
        )
    end
    println(
        io,
        "- `required_backend` preserves upstream toolbox facts from the inventory; ",
        "it does not select or endorse a Julia optimizer.",
    )
    println(
        io,
        "- `known_upstream_defects` carries the inventory's recorded numerical or ",
        "semantic ambiguity strings. An empty list is not proof that upstream has ",
        "no defect.",
    )
    println(
        io,
        "- Binding and file lists come only from direct QETLAB provenance rows. ",
        "An empty list can coexist with documented project-native component ",
        "evidence and is not filled by name guessing.",
    )
    println(
        io,
        "- Per-row benchmark files point to the shared benchmark driver only when ",
        "the inventory says a benchmark was recorded; the present provenance ",
        "schema has no case-level benchmark mapping.",
    )
    return String(take!(io))
end

function rendered_outputs(model)
    plan = render_completion_plan(model)
    queue = render_completion_queue(model)
    document = render_completion_document(model)

    parsed_plan = TOML.parse(plan)
    parsed_queue = TOML.parse(queue)
    length(get(parsed_plan, "public_rows", Any[])) == model.metadata.public_row_count ||
        error("rendered completion plan lost public rows")
    length(get(parsed_plan, "internal_helpers", Any[])) ==
    model.metadata.internal_helper_count ||
        error("rendered completion plan lost internal helpers")
    length(get(parsed_queue, "tasks", Any[])) == model.metadata.incomplete_public_count ||
        error("rendered completion queue lost tasks")
    for row in parsed_plan["public_rows"]
        missing = setdiff(REQUIRED_PUBLIC_ROW_FIELDS, Set(keys(row)))
        isempty(missing) || error(
            "rendered public row $(row["upstream_function"]) lacks fields: " *
            join(sort!(collect(missing)), ", "),
        )
    end
    return Dict(PLAN_PATH => plan, QUEUE_PATH => queue, DOCUMENT_PATH => document)
end

function write_or_check_outputs(; check::Bool)
    model = build_completion_model()
    outputs = rendered_outputs(model)
    success = true
    for path in (PLAN_PATH, QUEUE_PATH, DOCUMENT_PATH)
        expected = outputs[path]
        if check
            if !isfile(path)
                println(stderr, "missing generated completion artifact: ", path)
                success = false
            elseif read(path, String) != expected
                println(stderr, "stale generated completion artifact: ", path)
                success = false
            end
        else
            mkpath(dirname(path))
            open(path, "w") do io
                return write(io, expected)
            end
            println("wrote ", normalized_path(relpath(path, REPOSITORY_ROOT)))
        end
    end
    success || return false
    check && println(
        "checked ",
        model.metadata.public_row_count,
        " public completion rows, ",
        model.metadata.incomplete_public_count,
        " queued rows, and ",
        model.metadata.internal_helper_count,
        " internal-helper dispositions",
    )
    return true
end

function completion_report(
    model=build_completion_model(); audit=repository_quality_audit(model), io::IO=stdout
)
    println(
        io,
        "QETLAB completion ledger at ",
        model.metadata.source_revision,
        ": ",
        model.metadata.public_row_count,
        " public rows (",
        model.counts["implementation_complete"],
        " implementation-complete, ",
        model.counts["partial"],
        " partial, ",
        model.counts["deferred"],
        " deferred, ",
        model.counts["blocked"],
        " blocked); ",
        model.metadata.incomplete_public_count,
        " queued.",
    )
    terminal_count = count(row -> !row["required_for_completion"], model.helper_rows)
    required_count = count(row -> row["required_for_completion"], model.helper_rows)
    println(
        io,
        "Internal helpers: ",
        terminal_count,
        " terminal superseded/excluded and ",
        required_count,
        " assigned to reviewed completion work packages; source-call consumers ",
        "are recorded where detected.",
    )
    println(
        io,
        "Static release evidence: ",
        audit.metrics["runtime_export_count"],
        " runtime exports, ",
        audit.metrics["provenance_entry_count"],
        " provenance entries, ",
        audit.metrics["implementation_complete_rows_checked"],
        " implementation-complete rows, ",
        audit.metrics["randomized_complete_rows_checked"],
        " completed randomized rows, and ",
        audit.metrics["solver_backed_complete_rows_checked"],
        " completed solver-backed rows checked; ",
        audit.metrics["failure_count"],
        " failure(s).",
    )
    if !isempty(audit.failures)
        println(io, "Static release-evidence failures:")
        for failure in audit.failures
            println(io, "  - ", failure)
        end
    end
    println(io, "Static checker limitations:")
    for limitation in audit.limitations
        println(io, "  - ", limitation)
    end
    return model
end

function strict_completion_check(
    model=build_completion_model(); audit=repository_quality_audit(model), io::IO=stdout
)
    incomplete = [
        row for row in model.public_rows if
        !(row["current_status"] in FINAL_PUBLIC_TERMINAL_STATUSES)
    ]
    required_helpers = [row for row in model.helper_rows if row["required_for_completion"]]
    if isempty(incomplete) && isempty(required_helpers) && isempty(audit.failures)
        println(io, "strict QETLAB implementation-completeness gate passed")
        return true
    end
    println(io, "strict QETLAB implementation-completeness gate failed:")
    println(io, "  incomplete public rows: ", length(incomplete))
    nonterminal_complete = sort!(
        [
            row["upstream_function"] for
            row in incomplete if row["current_status"] in IMPLEMENTATION_COMPLETE_STATUSES
        ];
        by=lowercase,
    )
    isempty(nonterminal_complete) || println(
        io,
        "  implementation-complete but not final-terminal (",
        length(nonterminal_complete),
        "): ",
        join(nonterminal_complete, ", "),
    )
    for class in ("partial", "deferred", "blocked")
        names = sort!(
            [
                row["upstream_function"] for
                row in incomplete if public_status_class(row["current_status"]) == class
            ];
            by=lowercase,
        )
        isempty(names) ||
            println(io, "  ", class, " (", length(names), "): ", join(names, ", "))
    end
    helper_suffix = if isempty(required_helpers)
        ""
    else
        names = sort!([row["function_name"] for row in required_helpers]; by=lowercase)
        " (" * join(names, ", ") * ")"
    end
    println(io, "  required internal helpers: ", length(required_helpers), helper_suffix)
    println(io, "  static release-evidence failures: ", length(audit.failures))
    for failure in audit.failures
        println(io, "    - ", failure)
    end
    println(
        io,
        "This gate checks dispositions and static release evidence. It does not execute ",
        "tests or establish MATLAB parity, solver-certificate validity, or human release review.",
    )
    return false
end

export build_completion_model
export completion_report
export optimization_evidence_failures
export repository_quality_audit
export scan_repository_sources
export strict_completion_check
export write_or_check_outputs

end
