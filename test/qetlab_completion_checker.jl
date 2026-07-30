#!/usr/bin/env julia

# Focused mutation tests for the static completion checker. Run directly with:
# julia --startup-file=no --project=. test/qetlab_completion_checker.jl

using Test
using TOML
using QuantumEntanglementTools

include(joinpath(@__DIR__, "..", "scripts", "qetlab_completion_common.jl"))
using .QETLABCompletion

const TEST_RUNTIME_MODULES = Dict(
    "QuantumEntanglementTools" => QuantumEntanglementTools,
    "QuantumEntanglementTools.MATLABCompat" => QuantumEntanglementTools.MATLABCompat,
)

function mutated_model(model, function_name, status)
    rows = deepcopy(model.public_rows)
    row = only(filter(row -> row["upstream_function"] == function_name, rows))
    row["current_status"] = status
    return merge(model, (; public_rows=rows))
end

function terminal_disposition_model(model)
    public_rows = deepcopy(model.public_rows)
    for row in public_rows
        row["current_status"] = "verified"
        row["work_package"] = "completed_scope"
    end

    helper_rows = deepcopy(model.helper_rows)
    for row in helper_rows
        row["current_status"] = "internal_helper_superseded_by_parent_specific_iteration"
        row["disposition_class"] = "terminal_superseded"
        row["required_for_completion"] = false
        row["completion_consumers"] = String[]
        row["work_packages"] = String[]
    end

    packages = [
        merge(package, (; functions=String[], internal_helpers=String[])) for
        package in model.policy.packages
    ]
    policy = merge(
        model.policy,
        (;
            packages,
            assignments=Dict{String,String}(),
            helper_assignments=Dict{String,String}(),
        ),
    )
    counts = Dict(
        "implementation_complete" => length(public_rows),
        "partial" => 0,
        "deferred" => 0,
        "blocked" => 0,
    )
    helper_counts = Dict(class => 0 for class in keys(model.helper_counts))
    helper_counts["terminal_superseded"] = length(helper_rows)
    metadata = merge(model.metadata, (; incomplete_public_count=0))
    return merge(
        model,
        (;
            metadata,
            policy,
            public_rows,
            helper_rows,
            queue_rows=Dict{String,Any}[],
            counts,
            helper_counts,
        ),
    )
end

@testset "QETLAB completion checker" begin
    model = QETLABCompletion.build_completion_model()
    provenance = TOML.parsefile(QETLABCompletion.PROVENANCE_PATH)["functions"]

    baseline = QETLABCompletion.repository_quality_audit(
        model; provenance_entries=provenance, runtime_modules=TEST_RUNTIME_MODULES
    )
    @test isempty(baseline.failures)
    @test baseline.metrics["runtime_export_count"] == length(provenance)
    expected_complete = count(
        row ->
            QETLABCompletion.public_status_class(row["current_status"]) ==
            "implementation_complete",
        model.public_rows,
    )
    @test baseline.metrics["implementation_complete_rows_checked"] == expected_complete
    @test baseline.metrics["randomized_complete_rows_checked"] >= 7
    @test QETLABCompletion.public_status_class(
        "verified_with_documented_upstream_correction"
    ) == "implementation_complete"
    @test QETLABCompletion.strict_completion_check(
        mutated_model(model, "Tensor", "verified"); audit=baseline, io=IOBuffer()
    )
    @test !QETLABCompletion.strict_completion_check(
        mutated_model(model, "Tensor", "implemented"); audit=baseline, io=IOBuffer()
    )
    for status in (
        "internal_helper_superseded_by_parent_specific_construction_and_strict_normalization_validation",
        "internal_helper_superseded_by_public_typed_bell_behavior_conversion",
        "internal_helper_superseded_by_public_typed_nonlocal_game_conversion",
        "internal_helper_superseded_by_typed_parent_specific_canonical_nonlocal_representation",
    )
        @test QETLABCompletion.helper_status_class(status) == "terminal_superseded"
    end

    terminal_model = terminal_disposition_model(model)
    @test QETLABCompletion.validate_model(terminal_model) === terminal_model
    @test QETLABCompletion.strict_completion_check(
        terminal_model; audit=baseline, io=IOBuffer()
    )
    terminal_queue = QETLABCompletion.render_completion_queue(terminal_model)
    @test !occursin("[[tasks]]", terminal_queue)
    terminal_document = QETLABCompletion.render_completion_document(terminal_model)
    @test occursin(
        "implementation-completeness gate and passes for the generated ledger",
        terminal_document,
    )
    @test !occursin("currently fails", terminal_document)
    architecture_package_count = count(
        package -> package.architecture_only, terminal_model.policy.packages
    )
    @test length(collect(eachmatch(r"This is an architecture gate", terminal_document))) ==
        architecture_package_count
    @test occursin(
        "This work package is complete; no public row or required private helper remains queued.",
        terminal_document,
    )
    @test !occursin("Required private helper dispositions remain", terminal_document)
    @test occursin(
        "All 36 have terminal superseded/excluded dispositions", terminal_document
    )
    @test occursin(
        "Every implementation-complete solver-backed row must", terminal_document
    )
    @test occursin(
        "This strict pass does not by itself establish MATLAB parity", terminal_document
    )
    @test occursin(
        "validates the completed static ledger and its repository evidence",
        terminal_document,
    )

    helper_pending_rows = deepcopy(terminal_model.helper_rows)
    helper_pending = first(helper_pending_rows)
    helper_pending["current_status"] = "deferred_with_bcs_game_scope"
    helper_pending["disposition_class"] = "required_deferred"
    helper_pending["required_for_completion"] = true
    helper_pending["work_packages"] = ["algebraic_iterative_foundations"]
    helper_pending_packages = [
        if package.id == "algebraic_iterative_foundations"
            merge(package, (; internal_helpers=[helper_pending["function_name"]]))
        else
            package
        end for package in terminal_model.policy.packages
    ]
    helper_pending_policy = merge(
        terminal_model.policy,
        (;
            packages=helper_pending_packages,
            helper_assignments=Dict(
                helper_pending["function_name"] => "algebraic_iterative_foundations"
            ),
        ),
    )
    helper_pending_model = merge(
        terminal_model, (; helper_rows=helper_pending_rows, policy=helper_pending_policy)
    )
    helper_pending_document = QETLABCompletion.render_completion_document(
        helper_pending_model
    )
    @test occursin(
        "No public row is queued for this work package. Required private helper dispositions remain",
        helper_pending_document,
    )
    @test occursin(
        "implementation-completeness gate and currently fails on 0 queued public rows and 1 required private helpers",
        helper_pending_document,
    )

    missing_tests = deepcopy(provenance)
    tensor_entry = only(
        filter(
            entry ->
                get(entry, "source_project", "") == "QETLAB" &&
                get(entry, "upstream_function", "") == "Tensor" &&
                get(entry, "public_module", "") == "QuantumEntanglementTools" &&
                get(entry, "julia_name", "") == "tensor_product",
            missing_tests,
        ),
    )
    tensor_entry["tests"] = String[]
    tensor_entry["docs"] = ["docs/missing_completion_checker_fixture.md"]
    tensor_entry["julia_file"] = "src/missing_completion_checker_fixture.jl"
    tensor_entry["verification_status"] = "not_run"
    test_audit = QETLABCompletion.repository_quality_audit(
        model;
        provenance_entries=missing_tests,
        runtime_modules=TEST_RUNTIME_MODULES,
        scan_sources=false,
    )
    @test any(
        failure ->
            occursin("[files]", failure) &&
            occursin("tensor_product", failure) &&
            occursin("tests", failure),
        test_audit.failures,
    )
    @test any(
        failure ->
            occursin("[files]", failure) &&
            occursin("tensor_product", failure) &&
            occursin("missing docs path", failure),
        test_audit.failures,
    )
    @test any(
        failure ->
            occursin("[files]", failure) &&
            occursin("tensor_product", failure) &&
            occursin("missing julia_file path", failure),
        test_audit.failures,
    )
    @test any(
        failure ->
            occursin("[implemented-claim]", failure) &&
            occursin("tensor_product", failure) &&
            occursin("verification_status", failure),
        test_audit.failures,
    )
    all_complete_rows = deepcopy(model.public_rows)
    for row in all_complete_rows
        row["current_status"] = "verified"
    end
    terminal_helpers = deepcopy(model.helper_rows)
    for row in terminal_helpers
        row["required_for_completion"] = false
    end
    disposition_complete_model = merge(
        model, (; public_rows=all_complete_rows, helper_rows=terminal_helpers)
    )
    @test !QETLABCompletion.strict_completion_check(
        disposition_complete_model; audit=test_audit, io=IOBuffer()
    )

    missing_provenance = filter(
        entry -> !(
            get(entry, "source_project", "") == "QETLAB" &&
            get(entry, "upstream_function", "") == "Tensor"
        ),
        deepcopy(provenance),
    )
    provenance_audit = QETLABCompletion.repository_quality_audit(
        model;
        provenance_entries=missing_provenance,
        runtime_modules=TEST_RUNTIME_MODULES,
        scan_sources=false,
    )
    @test any(
        failure ->
            occursin("[implemented-claim]", failure) &&
            occursin("Tensor", failure) &&
            occursin("no direct QETLAB provenance", failure),
        provenance_audit.failures,
    )

    missing_rng_tests = deepcopy(provenance)
    for entry in missing_rng_tests
        get(entry, "source_project", "") == "QETLAB" || continue
        get(entry, "upstream_function", "") == "RandomGraph" || continue
        entry["tests"] = String[]
    end
    rng_audit = QETLABCompletion.repository_quality_audit(
        model;
        provenance_entries=missing_rng_tests,
        runtime_modules=TEST_RUNTIME_MODULES,
        scan_sources=false,
    )
    @test any(
        failure ->
            occursin("[rng]", failure) &&
            occursin("RandomGraph", failure) &&
            occursin("explicit-RNG binding test", failure),
        rng_audit.failures,
    )

    missing_solver_tests = deepcopy(provenance)
    for entry in missing_solver_tests
        get(entry, "source_project", "") == "QETLAB" || continue
        get(entry, "upstream_function", "") == "DiamondNorm" || continue
        entry["tests"] = filter(
            path -> !startswith(path, "test/extensions/"), get(entry, "tests", String[])
        )
    end
    solver_audit = QETLABCompletion.repository_quality_audit(
        model;
        provenance_entries=missing_solver_tests,
        runtime_modules=TEST_RUNTIME_MODULES,
        scan_sources=false,
    )
    @test any(
        failure ->
            occursin("[optimization]", failure) &&
            occursin("DiamondNorm", failure) &&
            occursin("optional-extension test", failure),
        solver_audit.failures,
    )

    mktempdir() do root
        source_dir = joinpath(root, "src")
        mkpath(source_dir)
        open(joinpath(source_dir, "first.jl"), "w") do io
            write(
                io,
                """
                # TODO replace this placeholder
                duplicated(x) = rand()
                Random.seed!(1)
                """,
            )
        end
        open(joinpath(source_dir, "second.jl"), "w") do io
            write(io, "duplicated(x) = x\n")
        end
        scan = QETLABCompletion.scan_repository_sources(
            root, Dict("QuantumEntanglementTools" => Set(["duplicated"]))
        )
        @test length(scan.markers) == 1
        @test length(scan.global_rng_candidates) == 1
        @test length(scan.seed_calls) == 1
        @test length(scan.duplicate_definitions) == 1
    end
end
