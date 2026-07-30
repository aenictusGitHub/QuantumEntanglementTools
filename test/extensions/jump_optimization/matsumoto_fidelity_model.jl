using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const MatsumotoOptQET = QuantumEntanglementTools

if !isdefined(MatsumotoOptQET, :MatsumotoFidelityModel)
    Base.include(
        MatsumotoOptQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "matsumoto_fidelity_model.jl"
        ),
    )
end

function matsumoto_hypatia_backend(; kwargs...)
    return MatsumotoOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=2.0e-7,
        rtol=2.0e-7,
        kwargs...,
    )
end

function matsumoto_scs_backend(; kwargs...)
    return MatsumotoOptQET.JuMPBackend(
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

@testset "WP3 Matsumoto-fidelity optional optimization" begin
    rho = Diagonal([0.7, 0.3])
    sigma = Diagonal([0.4, 0.6])
    expected = MatsumotoOptQET.matsumoto_fidelity(rho, sigma)
    problem = MatsumotoOptQET.matsumoto_fidelity_problem(rho, sigma)

    hypatia = MatsumotoOptQET.solve_optimization(problem, matsumoto_hypatia_backend())
    @test hypatia.status === MatsumotoOptQET.OptimizationOptimal
    @test hypatia.termination_status === :optimal
    @test hypatia.primal_status === :feasible_point
    @test hypatia.dual_status === :feasible_point
    @test hypatia.objective_value ≈ expected atol = 2e-7
    @test hypatia.objective_bound ≈ expected atol = 2e-7
    @test hypatia.absolute_gap <= 2e-7
    @test hypatia.primal_residual <= 2e-7
    @test hypatia.dual_residual <= 2e-7
    @test hypatia.primal !== nothing
    @test hypatia.dual !== nothing
    @test hasproperty(hypatia.primal.views, :matsumoto_fidelity_coupling)
    coupling = hypatia.primal.views.matsumoto_fidelity_coupling
    @test ishermitian(coupling)
    @test real(tr(coupling)) ≈ expected atol = 2e-7
    @test hypatia.optimizer.configured_optimizer_name == "Hypatia"
    @test !hypatia.certified

    scs = MatsumotoOptQET.solve_optimization(problem, matsumoto_scs_backend())
    @test scs.status in
        (MatsumotoOptQET.OptimizationOptimal, MatsumotoOptQET.OptimizationFeasible)
    @test scs.termination_status in (:optimal, :almost_optimal)
    @test scs.objective_value ≈ expected atol = 3e-5
    @test scs.primal_residual <= 3e-5
    @test scs.optimizer.reported_optimizer_name == "SCS"

    affine_rho = MatsumotoOptQET.hermitian_variable(:rho, 2; coefficient_type=Float64)
    affine_sigma = MatsumotoOptQET.HermitianAffineMatrix(
        :sigma, Matrix(sigma), Int[], Matrix{Float64}[], affine_rho.variable_count
    )
    fixed_rho = Matrix(rho)
    affine_problem = MatsumotoOptQET.matsumoto_fidelity_problem(
        affine_rho,
        affine_sigma;
        equalities=MatsumotoOptQET.hermitian_equalities(
            affine_rho; target=fixed_rho, name_prefix=:fixed_rho
        ),
        psd_constraints=[affine_rho],
        primal_views=[affine_rho],
        known_feasible_point=[0.7, 0.3, 0.0, 0.0],
        name=:affine_matsumoto,
    )
    affine_result = MatsumotoOptQET.solve_optimization(
        affine_problem, matsumoto_hypatia_backend()
    )
    @test affine_result.status === MatsumotoOptQET.OptimizationOptimal
    @test affine_result.objective_value ≈ expected atol = 2e-7
    @test affine_result.primal_residual <= 2e-7
    @test affine_result.primal.views.rho ≈ Matrix(rho) atol = 2e-7
    @test hasproperty(affine_result.primal.views, :affine_matsumoto_coupling)

    limited = MatsumotoOptQET.solve_optimization(
        problem, matsumoto_hypatia_backend(optimizer_options=(iter_limit=0,))
    )
    @test limited.status === MatsumotoOptQET.OptimizationLimit
    @test limited.termination_status === :iteration_limit
    @test !limited.certified

    malformed_backend = MatsumotoOptQET.JuMPBackend(
        () -> error("intentional Matsumoto-fidelity optimizer failure");
        optimizer_name="intentional malformed factory",
        allow_densify=true,
    )
    malformed = MatsumotoOptQET.solve_optimization(problem, malformed_backend)
    @test malformed.status === MatsumotoOptQET.OptimizationMalformedBackend
    @test malformed.termination_status === :backend_exception
    @test malformed.primal === nothing
    @test malformed.dual === nothing
    @test occursin("failed", malformed.message)
end
