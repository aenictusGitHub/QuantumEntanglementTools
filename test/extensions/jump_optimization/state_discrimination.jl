using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const StateDiscriminationQET = QuantumEntanglementTools

if !isdefined(StateDiscriminationQET, :StateDiscriminationResult)
    Base.include(
        StateDiscriminationQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "state_discrimination.jl"
        ),
    )
end

function discrimination_hypatia_backend(; kwargs...)
    return StateDiscriminationQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=2.0e-7,
        rtol=2.0e-7,
        kwargs...,
    )
end

function discrimination_scs_backend(; kwargs...)
    return StateDiscriminationQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=8.0e-6,
        rtol=8.0e-6,
        kwargs...,
    )
end

function trine_states()
    omega = cis(2pi / 3)
    return hcat(
        ComplexF64[1, 1] / sqrt(2),
        ComplexF64[1, omega] / sqrt(2),
        ComplexF64[1, omega ^ 2] / sqrt(2),
    )
end

@testset "WP4 state-discrimination optional optimization" begin
    @testset "Hypatia trine optimum and POVM" begin
        result = StateDiscriminationQET.state_distinguishability(
            trine_states(); backend=discrimination_hypatia_backend()
        )
        positional = StateDiscriminationQET.state_distinguishability(
            trine_states(), discrimination_hypatia_backend()
        )
        @test positional.status === result.status
        @test positional.success_probability ≈ result.success_probability
        @test result.status === StateDiscriminationQET.StateDiscriminationSolverOptimal
        @test result.optimization_result.status ===
            StateDiscriminationQET.OptimizationOptimal
        @test result.optimization_result.termination_status === :optimal
        @test result.optimization_result.primal_status === :feasible_point
        @test result.optimization_result.dual_status === :feasible_point
        @test result.success_probability ≈ 2 / 3 atol = 2e-7
        @test result.lower_bound ≈ 2 / 3 atol = 2e-7
        @test result.upper_bound ≈ 2 / 3 atol = 2e-7
        @test result.measurement !== nothing
        @test length(result.measurement) == 3
        @test result.residuals.valid
        @test result.residuals.completeness_residual <= 2e-7
        @test result.residuals.positivity_violation <= 2e-7
        @test result.residuals.solver_objective_residual <= 2e-7
        @test result.optimization_result.primal_residual <= 2e-7
        @test result.optimization_result.dual_residual <= 2e-7
        @test result.optimization_result.absolute_gap <= 2e-7
        @test result.optimization_result.optimizer.configured_optimizer_name == "Hypatia"
        @test !result.certified
        @test result.certificate_kind === nothing
    end

    @testset "SCS cross-check and identical states" begin
        scs = StateDiscriminationQET.state_distinguishability(
            trine_states(); backend=discrimination_scs_backend()
        )
        @test scs.status in (
            StateDiscriminationQET.StateDiscriminationSolverOptimal,
            StateDiscriminationQET.StateDiscriminationSolverFeasible,
        )
        @test scs.lower_bound ≈ 2 / 3 atol = 3e-5
        @test scs.measurement !== nothing
        @test scs.residuals.completeness_residual <= 3e-5
        @test scs.residuals.positivity_violation <= 3e-5
        @test scs.optimization_result.optimizer.reported_optimizer_name == "SCS"

        rho = ComplexF64[0.7 0.1; 0.1 0.3]
        identical = StateDiscriminationQET.state_distinguishability(
            (rho, rho, rho); backend=discrimination_hypatia_backend()
        )
        @test identical.status === StateDiscriminationQET.StateDiscriminationSolverOptimal
        @test identical.success_probability ≈ 1 / 3 atol = 2e-7
        @test identical.lower_bound ≈ 1 / 3 atol = 2e-7
        @test identical.upper_bound ≈ 1 / 3 atol = 2e-7
        @test identical.residuals.valid
    end

    @testset "backend limit and malformed factory" begin
        limited = StateDiscriminationQET.state_distinguishability(
            trine_states();
            backend=discrimination_hypatia_backend(optimizer_options=(iter_limit=0,)),
        )
        @test limited.status === StateDiscriminationQET.StateDiscriminationResourceLimit
        @test limited.optimization_result.status ===
            StateDiscriminationQET.OptimizationLimit
        @test limited.success_probability === nothing

        malformed_backend = StateDiscriminationQET.JuMPBackend(
            () -> error("intentional state-discrimination factory failure");
            optimizer_name="intentional malformed factory",
            allow_densify=true,
        )
        malformed = StateDiscriminationQET.state_distinguishability(
            trine_states(); backend=malformed_backend
        )
        @test malformed.status === StateDiscriminationQET.StateDiscriminationBackendFailure
        @test malformed.optimization_result.status ===
            StateDiscriminationQET.OptimizationMalformedBackend
        @test malformed.measurement === nothing
        @test malformed.success_probability === nothing
    end
end
