using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SCS
using Test

const CopCliqueOptQET = QuantumEntanglementTools

if !isdefined(CopCliqueOptQET, :CopositivityResult)
    Base.include(
        CopCliqueOptQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "copositivity_clique.jl"
        ),
    )
end

const CopCliqueExtension = Base.get_extension(
    CopCliqueOptQET, :QuantumEntanglementToolsJuMPExt
)
isnothing(CopCliqueExtension) && error("QuantumEntanglementToolsJuMPExt is not loaded")
if !isdefined(CopCliqueExtension, :_COPOSITIVITY_CLIQUE_EXTENSION_LOADED)
    Base.include(
        CopCliqueExtension,
        joinpath(@__DIR__, "..", "..", "..", "ext", "copositivity_clique.jl"),
    )
end

function copclique_hypatia_backend(; kwargs...)
    return CopCliqueOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=2.0e-7,
        rtol=2.0e-7,
        kwargs...,
    )
end

function copclique_scs_backend(; kwargs...)
    return CopCliqueOptQET.JuMPBackend(
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

function copclique_cycle_adjacency(vertices::Int)
    adjacency = zeros(Int, vertices, vertices)
    for vertex in 1:vertices
        neighbor = mod1(vertex + 1, vertices)
        adjacency[vertex, neighbor] = 1
        adjacency[neighbor, vertex] = 1
    end
    return adjacency
end

@testset "WP5 copositivity and clique-number optional optimization" begin
    numerical_psd = [1.0 -0.25; -0.25 1.0]
    numerical_psd_copy = copy(numerical_psd)
    exact_minimum = 3 / 8

    hypatia = CopCliqueOptQET.copositivity_criterion(
        MersenneTwister(1), numerical_psd, copclique_hypatia_backend()
    )
    @test hypatia.status === CopCliqueOptQET.CopositivityHierarchyUnknown
    @test hypatia.verdict === nothing
    @test !hypatia.certified
    @test hypatia.lower_bound ≈ exact_minimum atol = 4e-7
    @test hypatia.upper_bound === nothing
    @test hypatia.lower_kind === :sos_relaxation_numerical
    @test hypatia.hierarchy_result.status === CopCliqueOptQET.OptimizationOptimal
    @test hypatia.hierarchy_result.optimization_result.termination_status === :optimal
    @test hypatia.hierarchy_result.optimization_result.primal_status === :feasible_point
    @test hypatia.hierarchy_result.optimization_result.dual_status === :feasible_point
    @test hypatia.hierarchy_result.optimization_result.primal_residual <= 2e-7
    @test hypatia.hierarchy_result.optimization_result.dual_residual <= 2e-7
    @test hypatia.hierarchy_result.optimization_result.optimizer.configured_optimizer_name ==
        "Hypatia"
    @test numerical_psd == numerical_psd_copy

    scs = CopCliqueOptQET.copositivity_criterion(
        MersenneTwister(2), numerical_psd, copclique_scs_backend()
    )
    @test scs.status === CopCliqueOptQET.CopositivityHierarchyUnknown
    @test scs.verdict === nothing
    @test scs.lower_bound ≈ exact_minimum atol = 3e-5
    @test scs.hierarchy_result.status in
        (CopCliqueOptQET.OptimizationOptimal, CopCliqueOptQET.OptimizationFeasible)
    @test scs.hierarchy_result.optimization_result.primal_residual <= 3e-5
    @test scs.hierarchy_result.optimization_result.optimizer.reported_optimizer_name ==
        "SCS"

    sampled = CopCliqueOptQET.copositivity_criterion(
        MersenneTwister(3),
        numerical_psd,
        copclique_hypatia_backend();
        allow_densify=true,
        inner_samples=8,
    )
    @test sampled.status === CopCliqueOptQET.CopositivityHierarchyUnknown
    @test sampled.samples_requested == 8
    @test sampled.samples_evaluated == 8
    @test sampled.upper_bound !== nothing
    @test sampled.upper_bound >= exact_minimum - 2e-14
    @test sampled.upper_kind === :sampled_feasible
    @test sampled.hierarchy_result.best_point !== nothing

    cycle_five = copclique_cycle_adjacency(5)
    hypatia_clique = CopCliqueOptQET.clique_number_bounds(
        MersenneTwister(4), cycle_five, copclique_hypatia_backend()
    )
    @test hypatia_clique.status === CopCliqueOptQET.CliqueNumberNumericalHierarchy
    @test hypatia_clique.lower_bound == 2
    @test hypatia_clique.upper_bound == 3
    @test hypatia_clique.bounds_certified
    @test !hypatia_clique.exact
    @test hypatia_clique.continuous_upper_bound ≈ 1 - inv(sqrt(5)) atol = 5e-7
    @test hypatia_clique.uncertified_upper_candidate == 2
    @test hypatia_clique.hierarchy_result.status === CopCliqueOptQET.OptimizationOptimal
    @test hypatia_clique.hierarchy_result.optimization_result.primal_residual <= 2e-7
    @test occursin("not promoted", hypatia_clique.message)

    scs_clique = CopCliqueOptQET.clique_number_bounds(
        MersenneTwister(5), cycle_five, copclique_scs_backend()
    )
    @test scs_clique.status === CopCliqueOptQET.CliqueNumberNumericalHierarchy
    @test scs_clique.lower_bound == 2
    @test scs_clique.upper_bound == 3
    @test scs_clique.continuous_upper_bound ≈ 1 - inv(sqrt(5)) atol = 3e-5
    @test scs_clique.uncertified_upper_candidate == 2
    @test scs_clique.hierarchy_result.status in
        (CopCliqueOptQET.OptimizationOptimal, CopCliqueOptQET.OptimizationFeasible)
    @test scs_clique.hierarchy_result.optimization_result.primal_residual <= 3e-5

    limited = CopCliqueOptQET.copositivity_criterion(
        MersenneTwister(6),
        numerical_psd,
        copclique_hypatia_backend(optimizer_options=(iter_limit=0,)),
    )
    @test limited.status === CopCliqueOptQET.CopositivityResourceLimit
    @test limited.verdict === nothing
    @test limited.lower_bound === nothing
    @test limited.hierarchy_result.status === CopCliqueOptQET.OptimizationLimit
    @test limited.hierarchy_result.optimization_result.termination_status ===
        :iteration_limit

    malformed_backend = CopCliqueOptQET.JuMPBackend(
        () -> error("intentional copositivity optimizer failure");
        optimizer_name="intentional malformed factory",
        allow_densify=true,
    )
    malformed = CopCliqueOptQET.copositivity_criterion(
        MersenneTwister(7), numerical_psd, malformed_backend
    )
    @test malformed.status === CopCliqueOptQET.CopositivityBackendFailure
    @test malformed.verdict === nothing
    @test malformed.lower_bound === nothing
    @test malformed.hierarchy_result.status === CopCliqueOptQET.OptimizationMalformedBackend
    @test malformed.hierarchy_result.optimization_result.termination_status ===
        :backend_exception

    malformed_clique = CopCliqueOptQET.clique_number_bounds(
        MersenneTwister(8), cycle_five, malformed_backend
    )
    @test malformed_clique.status === CopCliqueOptQET.CliqueNumberBackendFailure
    @test malformed_clique.lower_bound == 2
    @test malformed_clique.upper_bound == 3
    @test malformed_clique.bounds_certified
    @test malformed_clique.continuous_upper_bound === nothing
    @test malformed_clique.uncertified_upper_candidate === nothing
end
