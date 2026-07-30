using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const SymExtOptQET = QuantumEntanglementTools

if !isdefined(SymExtOptQET, :SymmetricExtensionStatus)
    Base.include(
        SymExtOptQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "entanglement", "symmetric_extensions.jl"
        ),
    )
end

function symext_hypatia_backend(; options=NamedTuple(), kwargs...)
    return SymExtOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_options=options,
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=5.0e-7,
        rtol=5.0e-7,
        kwargs...,
    )
end

function symext_scs_backend(; kwargs...)
    return SymExtOptQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=2.0e-6, eps_rel=2.0e-6),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=1.0e-5,
        rtol=1.0e-5,
        kwargs...,
    )
end

function symext_bell()
    vector = [inv(sqrt(2.0)), 0, 0, inv(sqrt(2.0))]
    return vector * adjoint(vector)
end

@testset "WP3 symmetric-extension optimization paths" begin
    mixed = Matrix{Float64}(I, 4, 4) / 4
    hypatia = symext_hypatia_backend()

    @testset "outer primal extensions and residuals" begin
        general = SymExtOptQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            order=2,
            ppt=false,
            bosonic=false,
            backend=hypatia,
            prefer_analytic=false,
            allow_densify=true,
        )
        @test general.status === SymExtOptQET.SymmetricExtensionSolverPresent
        @test general.verdict === true
        @test general.certificate_kind === :validated_primal_extension
        @test general.optimization_result.status === SymExtOptQET.OptimizationOptimal
        @test general.optimization_result.termination_status === :optimal
        @test general.optimization_result.primal_status === :feasible_point
        @test general.optimization_result.dual_status === :feasible_point
        @test general.optimization_result.primal !== nothing
        @test general.optimization_result.dual !== nothing
        @test general.optimization_result.objective_value ≈ 0 atol = 1e-12
        @test general.optimization_result.primal_residual <= general.tolerance
        @test general.optimization_result.dual_residual <= general.tolerance
        @test general.residuals.valid
        @test general.residuals.minimum_eigenvalue >= -general.tolerance
        @test general.residuals.trace_residual <= general.tolerance
        @test general.residuals.marginal_residual <= general.tolerance
        @test general.residuals.permutation_residual <= general.tolerance
        @test size(general.extension) == (8, 8)

        bosonic_ppt = SymExtOptQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            order=2,
            ppt=true,
            bosonic=true,
            backend=hypatia,
            prefer_analytic=false,
            allow_densify=true,
        )
        @test bosonic_ppt.status === SymExtOptQET.SymmetricExtensionSolverPresent
        @test bosonic_ppt.verdict === true
        @test bosonic_ppt.residuals.valid
        @test bosonic_ppt.residuals.bosonic_support_residual <= bosonic_ppt.tolerance
        @test bosonic_ppt.residuals.permutation_residual <= bosonic_ppt.tolerance
        @test length(bosonic_ppt.residuals.ppt_minimum_eigenvalues) == 2
        @test bosonic_ppt.residuals.ppt_violation <= bosonic_ppt.tolerance

        scs = SymExtOptQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            order=2,
            bosonic=true,
            backend=symext_scs_backend(),
            prefer_analytic=false,
            allow_densify=true,
        )
        @test scs.status in (
            SymExtOptQET.SymmetricExtensionSolverPresent,
            SymExtOptQET.SymmetricExtensionNumericalBoundary,
        )
        @test scs.optimization_result.optimizer.reported_optimizer_name == "SCS"
        @test scs.extension !== nothing
        @test scs.residuals.marginal_residual <= 5e-5
    end

    @testset "outer dual separator" begin
        bell = symext_bell()
        result = SymExtOptQET.symmetric_extension(
            bell;
            dims=(2, 2),
            order=2,
            ppt=false,
            bosonic=false,
            backend=hypatia,
            prefer_analytic=false,
            allow_densify=true,
        )
        @test result.optimization_result.status === SymExtOptQET.OptimizationInfeasible
        @test result.status === SymExtOptQET.SymmetricExtensionSolverAbsent
        @test result.verdict === false
        @test result.certificate_kind === :validated_symmetric_extension_separator
        @test result.extension === nothing
        @test result.witness !== nothing
        @test result.witness.separator_validated
        @test result.witness.entanglement_witness
        @test result.witness.dual_stationarity_residual <= result.tolerance
        @test real(dot(result.witness.operator, bell)) ≈ -1 atol = 2e-12
        @test result.witness.normalization_residual <= 2e-12
    end

    @testset "inner hierarchy primal and dual semantics" begin
        for ppt in (false, true)
            result = SymExtOptQET.symmetric_inner_extension(
                mixed; dims=(2, 2), order=2, ppt=ppt, backend=hypatia, allow_densify=true
            )
            @test result.status === SymExtOptQET.SymmetricExtensionSolverPresent
            @test result.verdict === true
            @test result.hierarchy === :inner
            @test result.bosonic
            @test result.residuals.valid
            @test result.residuals.marginal_residual <= result.tolerance
            @test result.residuals.permutation_residual <= result.tolerance
            @test result.residuals.bosonic_support_residual <= result.tolerance
            @test result.optimization_result.status === SymExtOptQET.OptimizationOptimal
            if ppt
                @test result.problem.mixing_parameter ≈ 1 - inv(sqrt(3)) atol = 2e-14
                @test length(result.residuals.ppt_minimum_eigenvalues) == 1
                @test result.residuals.ppt_violation <= result.tolerance
            end
        end

        bell = symext_bell()
        outside = SymExtOptQET.symmetric_inner_extension(
            bell; dims=(2, 2), order=2, ppt=false, backend=hypatia, allow_densify=true
        )
        @test outside.optimization_result.status === SymExtOptQET.OptimizationInfeasible
        @test outside.status === SymExtOptQET.SymmetricExtensionSolverAbsent
        @test outside.verdict === false
        @test outside.certificate_kind === :validated_inner_cone_separator
        @test outside.witness.separator_validated
        @test !outside.witness.entanglement_witness
        @test occursin("not automatically", outside.witness.warning)
        @test occursin("does not prove entanglement", outside.message)
        @test real(dot(outside.witness.operator, bell)) ≈ -1 atol = 2e-12
    end

    @testset "MATLAB compatibility solver paths" begin
        bell = symext_bell()
        outer = SymExtOptQET.MATLABCompat.SymmetricExtension(
            bell,
            2,
            [2, 2],
            0,
            0,
            nothing;
            return_witness=true,
            backend=hypatia,
            allow_densify=true,
        )
        @test outer.ex == 0
        @test real(dot(outer.wit, bell)) ≈ -1 atol = 2e-12

        inner = SymExtOptQET.MATLABCompat.SymmetricInnerExtension(
            bell,
            2,
            [2, 2],
            0,
            nothing;
            return_witness=true,
            backend=hypatia,
            allow_densify=true,
        )
        @test inner.ex == 0
        @test real(dot(inner.wit, bell)) ≈ -1 atol = 2e-12

        structured_inner = SymExtOptQET.MATLABCompat.SymmetricInnerExtension(
            bell, 2, [2, 2]; structured=true, backend=hypatia, allow_densify=true
        )
        @test structured_inner.status === SymExtOptQET.SymmetricExtensionSolverAbsent
        @test !structured_inner.witness.entanglement_witness
        @test occursin("does not prove entanglement", structured_inner.message)
    end

    @testset "backend limit and malformed factory" begin
        limited = SymExtOptQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            backend=symext_hypatia_backend(options=(iter_limit=0,)),
            prefer_analytic=false,
            allow_densify=true,
        )
        @test limited.optimization_result.status === SymExtOptQET.OptimizationLimit
        @test limited.status === SymExtOptQET.SymmetricExtensionResourceLimit
        @test limited.verdict === nothing

        malformed_backend = SymExtOptQET.JuMPBackend(
            () -> 42; optimizer_name="malformed test factory", allow_densify=true
        )
        malformed = SymExtOptQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            backend=malformed_backend,
            prefer_analytic=false,
            allow_densify=true,
        )
        @test malformed.optimization_result.status ===
            SymExtOptQET.OptimizationMalformedBackend
        @test malformed.status === SymExtOptQET.SymmetricExtensionBackendFailure
        @test malformed.verdict === nothing
        @test malformed.extension === nothing
        @test occursin("failed", malformed.message)
    end
end
