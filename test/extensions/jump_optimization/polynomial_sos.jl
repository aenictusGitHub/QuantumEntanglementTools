using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SCS
using Test

const PolynomialSOSOptQET = QuantumEntanglementTools

if !isdefined(PolynomialSOSOptQET, :PolynomialSOSResult)
    Base.include(
        PolynomialSOSOptQET,
        joinpath(@__DIR__, "..", "..", "..", "src", "optimization", "polynomial_sos.jl"),
    )
end

const PolynomialSOSExtension = Base.get_extension(
    PolynomialSOSOptQET, :QuantumEntanglementToolsJuMPExt
)
isnothing(PolynomialSOSExtension) && error("QuantumEntanglementToolsJuMPExt is not loaded")
if !isdefined(PolynomialSOSExtension, :_POLYNOMIAL_SOS_EXTENSION_LOADED)
    Base.include(
        PolynomialSOSExtension,
        joinpath(@__DIR__, "..", "..", "..", "ext", "polynomial_sos.jl"),
    )
end

function polynomial_sos_hypatia_backend(; kwargs...)
    return PolynomialSOSOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=2.0e-7,
        rtol=2.0e-7,
        kwargs...,
    )
end

function polynomial_sos_scs_backend(; kwargs...)
    return PolynomialSOSOptQET.JuMPBackend(
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

@testset "WP5 polynomial SOS optional optimization" begin
    polynomial = PolynomialSOSOptQET.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
    exact_maximum = 2 + sqrt(2)
    exact_minimum = 2 - sqrt(2)

    maximum = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(1),
        polynomial,
        polynomial_sos_hypatia_backend();
        sense=:max,
        inner_samples=8,
    )
    @test maximum.status === PolynomialSOSOptQET.OptimizationOptimal
    @test maximum.optimization_result.termination_status === :optimal
    @test maximum.optimization_result.primal_status === :feasible_point
    @test maximum.optimization_result.dual_status === :feasible_point
    @test maximum.outer_bound ≈ exact_maximum atol = 3e-7
    @test maximum.inner_bound <= exact_maximum + 2e-14
    @test maximum.moment_matrix !== nothing
    @test ishermitian(maximum.moment_matrix)
    @test real(tr(maximum.moment_matrix)) ≈ 1 atol = 2e-7
    @test minimum(eigvals(Hermitian(maximum.moment_matrix))) >= -2e-7
    @test maximum.optimization_result.primal_residual <= 2e-7
    @test maximum.optimization_result.dual_residual <= 2e-7
    @test maximum.optimization_result.absolute_gap <= 3e-7
    @test maximum.optimization_result.optimizer.configured_optimizer_name == "Hypatia"
    @test !maximum.certified_outer

    minimum_result = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(2), polynomial; backend=polynomial_sos_hypatia_backend(), sense=:min
    )
    @test minimum_result.status === PolynomialSOSOptQET.OptimizationOptimal
    @test minimum_result.outer_bound ≈ exact_minimum atol = 3e-7
    @test minimum_result.inner_bound === nothing

    level_one = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(3), polynomial; backend=polynomial_sos_hypatia_backend(), level=1
    )
    @test level_one.status === PolynomialSOSOptQET.OptimizationOptimal
    @test level_one.outer_bound ≈ exact_maximum atol = 4e-7
    @test level_one.problem.metadata.hierarchy_level == 1
    @test level_one.problem.metadata.moment_dimension == 3

    scs = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(4), polynomial; backend=polynomial_sos_scs_backend()
    )
    @test scs.status in (
        PolynomialSOSOptQET.OptimizationOptimal, PolynomialSOSOptQET.OptimizationFeasible
    )
    @test scs.outer_bound ≈ exact_maximum atol = 3e-5
    @test scs.optimization_result.primal_residual <= 3e-5
    @test scs.optimization_result.optimizer.reported_optimizer_name == "SCS"

    target_rng = MersenneTwister(99)
    untouched_rng = copy(target_rng)
    target = PolynomialSOSOptQET.polynomial_sos_bounds(
        target_rng,
        polynomial;
        backend=polynomial_sos_hypatia_backend(),
        target=4.0,
        inner_samples=10,
    )
    @test target.target_status === :outer_proves_at_most
    @test target.inner_kind === :skipped_by_outer_target
    @test target.samples_evaluated == 0
    @test rand(target_rng) == rand(untouched_rng)

    limited = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(5),
        polynomial;
        backend=polynomial_sos_hypatia_backend(optimizer_options=(iter_limit=0,)),
    )
    @test limited.status === PolynomialSOSOptQET.OptimizationLimit
    @test limited.optimization_result.termination_status === :iteration_limit
    @test limited.outer_bound === nothing
    @test limited.moment_matrix !== nothing

    malformed_backend = PolynomialSOSOptQET.JuMPBackend(
        () -> error("intentional polynomial-SOS optimizer failure");
        optimizer_name="intentional malformed factory",
        allow_densify=true,
    )
    malformed = PolynomialSOSOptQET.polynomial_sos_bounds(
        MersenneTwister(6), polynomial; backend=malformed_backend
    )
    @test malformed.status === PolynomialSOSOptQET.OptimizationMalformedBackend
    @test malformed.optimization_result.termination_status === :backend_exception
    @test malformed.outer_bound === nothing
    @test malformed.moment_matrix === nothing
    @test occursin("did not produce", malformed.message)
end
