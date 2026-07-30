#!/usr/bin/env julia

# Read-only consistency check for exported bindings, public provenance, and the
# reviewed QETLAB inventory. This script deliberately does not rewrite ledgers.

using SHA
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)
using QuantumEntanglementTools: QuantumEntanglementTools

const PROVENANCE_PATH = joinpath(REPOSITORY_ROOT, "PROVENANCE.toml")
const INVENTORY_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_inventory.toml")
const STATUS_PATH = joinpath(REPOSITORY_ROOT, "porting", "qetlab_status.toml")
const UPSTREAM_MANIFEST_PATH = joinpath(REPOSITORY_ROOT, "UpstreamManifest.toml")

const PUBLIC_MODULES = Dict(
    "QuantumEntanglementTools" => QuantumEntanglementTools,
    "QuantumEntanglementTools.MATLABCompat" => QuantumEntanglementTools.MATLABCompat,
)

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

const QETLAB_PROVENANCE_FIELDS = Set([
    "upstream_function", "source_path", "source_revision", "source_sha256"
])

const ALLOWED_PUBLIC_KINDS = Set([
    "module", "type", "function", "constant", "alias", "compatibility_wrapper"
])
const ALLOWED_API_STATUSES = Set(["implemented", "verified", "compatibility_alias"])
const IMPLEMENTED_INVENTORY_STATUSES = Set([
    "implemented",
    "verified",
    "verified_with_documented_upstream_correction",
    "compatibility_alias",
    "superseded_with_documented_mapping",
])
const STATUS_IDENTITY_FIELDS = Set(["function_name", "source_path"])
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

function public_names(mod::Module)
    return Set(
        String(name) for
        name in names(mod; all=false, imported=false) if name != nameof(mod)
    )
end

function qualified_key(module_name::AbstractString, julia_name::AbstractString)
    return (String(module_name), String(julia_name))
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

function resolve_qualified_binding(path::AbstractString)
    parts = split(String(path), '.')
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

function existing_local_path(path::AbstractString)
    startswith(path, "http://") && return true
    startswith(path, "https://") && return true
    return ispath(joinpath(REPOSITORY_ROOT, path))
end

function normalized_overlay_value(field::AbstractString, value)
    if field == "review_notes" && value isa AbstractVector
        return sort!(String.(value); by=lowercase)
    end
    return value
end

function qetlab_source(manifest)
    matches = filter(
        source -> get(source, "id", nothing) == "qetlab", get(manifest, "sources", Any[])
    )
    length(matches) == 1 ||
        error("UpstreamManifest.toml must contain exactly one source with id = \"qetlab\"")
    return only(matches)
end

function main(args)
    isempty(args) || error("scripts/check_public_api.jl does not accept arguments")
    required_files = [PROVENANCE_PATH, INVENTORY_PATH, STATUS_PATH, UPSTREAM_MANIFEST_PATH]
    for path in required_files
        isfile(path) || error("required ledger does not exist: $path")
    end

    failures = String[]
    fail(message) = push!(failures, String(message))

    provenance = TOML.parsefile(PROVENANCE_PATH)
    inventory = TOML.parsefile(INVENTORY_PATH)
    status_overlay = TOML.parsefile(STATUS_PATH)
    upstream_manifest = TOML.parsefile(UPSTREAM_MANIFEST_PATH)

    get(provenance, "schema_version", nothing) == 1 ||
        fail("PROVENANCE.toml schema_version must equal 1")
    get(inventory, "schema_version", nothing) == 1 ||
        fail("qetlab_inventory.toml schema_version must equal 1")
    get(status_overlay, "schema_version", nothing) == 1 ||
        fail("qetlab_status.toml schema_version must equal 1")

    provenance_entries = get(provenance, "functions", Any[])
    provenance_entries isa AbstractVector ||
        error("PROVENANCE.toml `functions` must be an array of tables")
    inventory_entries = get(inventory, "functions", Any[])
    inventory_entries isa AbstractVector ||
        error("qetlab_inventory.toml `functions` must be an array of tables")
    status_entries = get(status_overlay, "functions", Any[])
    status_entries isa AbstractVector ||
        error("qetlab_status.toml `functions` must be an array of tables")

    exported = Set{Tuple{String,String}}()
    for (module_name, mod) in PUBLIC_MODULES
        union!(exported, (qualified_key(module_name, name) for name in public_names(mod)))
    end

    inventory_by_name = Dict{String,Any}()
    inventory_by_path = Dict{String,Any}()
    for (entry_number, entry) in pairs(inventory_entries)
        entry isa AbstractDict || error("inventory entry $entry_number is not a table")
        name = get(entry, "function_name", nothing)
        path = get(entry, "source_path", nothing)
        name isa AbstractString ||
            error("inventory entry $entry_number has no string function_name")
        path isa AbstractString ||
            error("inventory entry $entry_number has no string source_path")
        haskey(inventory_by_name, name) && fail("duplicate inventory function_name: $name")
        haskey(inventory_by_path, path) && fail("duplicate inventory source_path: $path")
        inventory_by_name[name] = entry
        inventory_by_path[path] = entry
    end

    provenance_by_key = Dict{Tuple{String,String},Any}()
    qetlab_provenance_by_function = Dict{String,Vector{Any}}()
    for (entry_number, entry) in pairs(provenance_entries)
        entry isa AbstractDict || error("provenance entry $entry_number is not a table")
        missing = setdiff(REQUIRED_PROVENANCE_FIELDS, Set(keys(entry)))
        isempty(missing) || fail(
            "provenance entry $entry_number is missing fields: " *
            join(sort!(collect(missing)), ", "),
        )
        haskey(entry, "julia_name") && haskey(entry, "public_module") || continue
        name = entry["julia_name"]
        module_name = entry["public_module"]
        name isa AbstractString ||
            (fail("provenance entry $entry_number julia_name must be a string"); continue)
        module_name isa AbstractString || (
            fail("provenance entry $entry_number public_module must be a string");
            continue
        )
        key = qualified_key(module_name, name)
        if haskey(provenance_by_key, key)
            fail("duplicate provenance entry for $(module_name).$(name)")
            continue
        end
        provenance_by_key[key] = entry

        haskey(PUBLIC_MODULES, module_name) ||
            fail("unknown public_module in provenance: $module_name")
        if haskey(PUBLIC_MODULES, module_name)
            mod = PUBLIC_MODULES[module_name]
            symbol = Symbol(name)
            if !isdefined(mod, symbol)
                fail("provenance binding is undefined: $(module_name).$(name)")
            else
                actual_kind = binding_kind(getfield(mod, symbol))
                declared_kind = get(entry, "public_kind", "")
                declared_kind in ALLOWED_PUBLIC_KINDS ||
                    fail("invalid public_kind '$declared_kind' for $(module_name).$(name)")
                kind_matches(declared_kind, actual_kind) || fail(
                    "public_kind '$declared_kind' does not match runtime kind " *
                    "'$actual_kind' for $(module_name).$(name)",
                )
            end
        end

        status = get(entry, "status", "")
        status in ALLOWED_API_STATUSES ||
            fail("unsupported status '$status' for $(module_name).$(name)")
        for field in ("julia_file", "specification")
            value = get(entry, field, "")
            value isa AbstractString && !isempty(value) ||
                fail("$field must be a non-empty string for $(module_name).$(name)")
            value isa AbstractString && existing_local_path(value) ||
                fail("$field path does not exist for $(module_name).$(name): $value")
        end
        for field in ("tests", "docs")
            values = get(entry, field, nothing)
            values isa AbstractVector && all(value -> value isa AbstractString, values) ||
                begin
                    fail("$field must be an array of strings for $(module_name).$(name)")
                    continue
                end
            for value in values
                existing_local_path(value) ||
                    fail("$field path does not exist for $(module_name).$(name): $value")
            end
        end

        declared_kind = get(entry, "public_kind", "")
        if declared_kind in ("alias", "compatibility_wrapper")
            delegates = get(entry, "delegates_to", nothing)
            delegates isa AbstractVector &&
            !isempty(delegates) &&
            all(delegate -> delegate isa AbstractString, delegates) || begin
                fail("$(module_name).$(name) must declare non-empty delegates_to")
                delegates = String[]
            end
            for delegate in delegates
                isnothing(resolve_qualified_binding(delegate)) &&
                    fail("unresolved delegate for $(module_name).$(name): $delegate")
            end
        end

        source_project = get(entry, "source_project", "")
        if source_project == "QETLAB"
            missing_qetlab = setdiff(QETLAB_PROVENANCE_FIELDS, Set(keys(entry)))
            isempty(missing_qetlab) || fail(
                "QETLAB provenance for $(module_name).$(name) is missing: " *
                join(sort!(collect(missing_qetlab)), ", "),
            )
            haskey(entry, "upstream_function") || continue
            upstream_function = entry["upstream_function"]
            source_path = get(entry, "source_path", "")
            if !haskey(inventory_by_name, upstream_function)
                fail(
                    "QETLAB provenance for $(module_name).$(name) names absent " *
                    "inventory function $upstream_function",
                )
                continue
            end
            inventory_record = inventory_by_name[upstream_function]
            get(inventory_record, "source_path", nothing) == source_path || fail(
                "source_path mismatch for $(module_name).$(name): " *
                "$source_path versus inventory $(inventory_record["source_path"])",
            )
            get(inventory_record, "source_revision", nothing) ==
            get(entry, "source_revision", nothing) ||
                fail("source_revision mismatch for $(module_name).$(name)")
            get(inventory_record, "source_sha256", nothing) ==
            get(entry, "source_sha256", nothing) ||
                fail("source_sha256 mismatch for $(module_name).$(name)")
            push!(get!(qetlab_provenance_by_function, upstream_function, Any[]), entry)
        elseif source_project == "QuantumEntanglementTools"
            origin = get(entry, "project_origin", nothing)
            origin isa AbstractString && !isempty(strip(origin)) || fail(
                "project-native provenance for $(module_name).$(name) " *
                "must declare project_origin",
            )
        else
            fail("unsupported source_project '$source_project' for $(module_name).$(name)")
        end
    end

    provenance_keys = Set(keys(provenance_by_key))
    for (module_name, name) in sort!(collect(setdiff(exported, provenance_keys)))
        fail("exported binding has no provenance: $(module_name).$(name)")
    end
    for (module_name, name) in sort!(collect(setdiff(provenance_keys, exported)))
        fail("provenance entry is not exported: $(module_name).$(name)")
    end

    qetlab = qetlab_source(upstream_manifest)
    pinned_revision = get(qetlab, "exact_revision", nothing)
    get(inventory, "source_revision", nothing) == pinned_revision ||
        fail("inventory source_revision does not match UpstreamManifest QETLAB pin")
    get(provenance, "qetlab_revision", nothing) == pinned_revision ||
        fail("PROVENANCE.toml qetlab_revision does not match UpstreamManifest QETLAB pin")
    for entries in values(qetlab_provenance_by_function), entry in entries
        get(entry, "source_revision", nothing) == pinned_revision || fail(
            "provenance source_revision for $(entry["public_module"])." *
            "$(entry["julia_name"]) does not match UpstreamManifest QETLAB pin",
        )
    end

    status_keys = Set{String}()
    for (entry_number, overlay_entry) in pairs(status_entries)
        overlay_entry isa AbstractDict || error("status entry $entry_number is not a table")
        unknown = setdiff(
            Set(keys(overlay_entry)), union(STATUS_IDENTITY_FIELDS, STATUS_OVERLAY_FIELDS)
        )
        isempty(unknown) || fail(
            "status entry $entry_number has unsupported fields: " *
            join(sort!(collect(unknown)), ", "),
        )
        by_name = if haskey(overlay_entry, "function_name")
            get(inventory_by_name, overlay_entry["function_name"], nothing)
        else
            nothing
        end
        by_path = if haskey(overlay_entry, "source_path")
            get(inventory_by_path, overlay_entry["source_path"], nothing)
        else
            nothing
        end
        if isnothing(by_name) && isnothing(by_path)
            fail("status entry $entry_number does not identify an inventory row")
            continue
        elseif !isnothing(by_name) && !isnothing(by_path) && by_name !== by_path
            fail("status entry $entry_number function_name/source_path disagree")
            continue
        end
        inventory_record = something(by_name, by_path)
        function_name = inventory_record["function_name"]
        function_name in status_keys && fail("duplicate status entry for $function_name")
        push!(status_keys, function_name)
        for field in STATUS_OVERLAY_FIELDS
            haskey(overlay_entry, field) || continue
            expected = normalized_overlay_value(field, overlay_entry[field])
            actual = normalized_overlay_value(field, get(inventory_record, field, nothing))
            expected == actual ||
                fail("generated inventory is stale for $function_name field $field")
        end
    end

    expected_status_hash = bytes2hex(sha256(read(STATUS_PATH)))
    get(inventory, "status_overlay_sha256", nothing) == expected_status_hash ||
        fail("generated inventory status_overlay_sha256 is stale")
    get(inventory, "status_overlay_entry_count", nothing) == length(status_entries) ||
        fail("generated inventory status_overlay_entry_count is stale")

    missing_status_rows = sort!(collect(setdiff(Set(keys(inventory_by_name)), status_keys)))
    isempty(missing_status_rows) || fail(
        "inventory rows without reviewed status entries: " *
        join(missing_status_rows, ", "),
    )
    pending_review_rows = sort!([
        record["function_name"] for record in inventory_entries if
        get(record, "review_status", "") == "automated_parse_pending_manual_review"
    ])
    isempty(pending_review_rows) || fail(
        "inventory rows still pending source review: " * join(pending_review_rows, ", ")
    )
    get(inventory, "source_reviewed_count", nothing) ==
    length(inventory_entries) - length(pending_review_rows) ||
        fail("generated inventory source_reviewed_count is stale")
    get(inventory, "pending_review_count", nothing) == length(pending_review_rows) ||
        fail("generated inventory pending_review_count is stale")
    get(inventory, "review_status", nothing) == "source_review_complete" ||
        fail("generated inventory top-level review_status is not source_review_complete")

    for inventory_record in inventory_entries
        get(inventory_record, "implementation_status", "not_started") in
        IMPLEMENTED_INVENTORY_STATUSES || continue
        function_name = inventory_record["function_name"]
        provenance_entries_for_function = get(
            qetlab_provenance_by_function, function_name, Any[]
        )
        isempty(provenance_entries_for_function) &&
            fail("implemented inventory row has no public provenance: $function_name")

        native_name = get(inventory_record, "proposed_julia_name", "")
        native_key = qualified_key("QuantumEntanglementTools", native_name)
        native_key in exported || fail(
            "implemented inventory row maps to non-exported native name: $function_name -> $native_name",
        )
        native_key in provenance_keys || fail(
            "implemented inventory native mapping lacks provenance: $function_name -> $native_name",
        )

        compatibility_name = get(inventory_record, "compatibility_alias", "")
        compatibility_key = qualified_key(
            "QuantumEntanglementTools.MATLABCompat", compatibility_name
        )
        compatibility_key in exported || fail(
            "implemented inventory row maps to non-exported compatibility name: " *
            "$function_name -> $compatibility_name",
        )
        compatibility_key in provenance_keys || fail(
            "implemented inventory compatibility mapping lacks provenance: " *
            "$function_name -> $compatibility_name",
        )
    end

    if isempty(failures)
        println(
            "checked ",
            length(exported),
            " public bindings (",
            length(public_names(QuantumEntanglementTools)),
            " native/module, ",
            length(public_names(QuantumEntanglementTools.MATLABCompat)),
            " MATLABCompat), ",
            length(provenance_entries),
            " provenance entries, and ",
            length(status_entries),
            " reviewed inventory rows",
        )
        return nothing
    end

    println(stderr, "public API/provenance consistency check failed:")
    for message in failures
        println(stderr, "  - ", message)
    end
    return exit(1)
end

main(ARGS)
