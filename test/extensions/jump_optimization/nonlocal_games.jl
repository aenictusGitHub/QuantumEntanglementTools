using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SCS
using Test

const NonlocalOptQET = QuantumEntanglementTools

for source_file in ("scenarios.jl", "npa.jl", "game_values.jl", "nonlocal_optimization.jl")
    symbol = if source_file == "scenarios.jl"
        :BellScenario
    elseif source_file == "npa.jl"
        :NPAWord
    elseif source_file == "game_values.jl"
        :NonlocalValueResult
    else
        :NonlocalLowerBoundResult
    end
    if !isdefined(NonlocalOptQET, symbol)
        Base.include(
            NonlocalOptQET,
            joinpath(@__DIR__, "..", "..", "..", "src", "nonlocal_games", source_file),
        )
    end
end

const NonlocalExtension = Base.get_extension(
    NonlocalOptQET, :QuantumEntanglementToolsJuMPExt
)
isnothing(NonlocalExtension) && error("QuantumEntanglementToolsJuMPExt is not loaded")
if !isdefined(NonlocalExtension, :_NONLOCAL_GAMES_EXTENSION_LOADED)
    Base.include(
        NonlocalExtension, joinpath(@__DIR__, "..", "..", "..", "ext", "nonlocal_games.jl")
    )
end

function nonlocal_hypatia_backend(; kwargs...)
    return NonlocalOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=2.0e-7,
        rtol=2.0e-7,
        kwargs...,
    )
end

function nonlocal_scs_backend(; kwargs...)
    return NonlocalOptQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=2.0e-7, eps_rel=2.0e-7),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=2.0e-5,
        rtol=2.0e-5,
        kwargs...,
    )
end

function _nonlocal_pr_box()
    probabilities = zeros(Float64, 2, 2, 2, 2)
    for x in 0:1, y in 0:1, a in 0:1, b in 0:1
        xor(a, b) == x * y || continue
        probabilities[a + 1, b + 1, x + 1, y + 1] = 0.5
    end
    return NonlocalOptQET.FullProbabilityBehavior(
        probabilities, NonlocalOptQET.BellScenario(2, 2, 2, 2); atol=0, rtol=0
    )
end

@testset "WP6 nonlocal-game optional JuMP optimization" begin
    hypatia = nonlocal_hypatia_backend()
    scs = nonlocal_scs_backend()
    probabilities = fill(0.25, 2, 2)
    parity = [0 0; 0 1]
    expected_xor = (2 + sqrt(2)) / 4

    @testset "XOR SDP and explicit extension dispatch" begin
        hypatia_xor = NonlocalOptQET.xor_game_value(
            probabilities, parity, hypatia; regime=:quantum
        )
        @test hypatia_xor.status === NonlocalOptQET.NonlocalValueNumericalUpperBound
        @test hypatia_xor.value === nothing
        @test hypatia_xor.upper_bound ≈ expected_xor atol = 5e-7
        @test !hypatia_xor.exact
        @test !hypatia_xor.certified_upper
        @test hypatia_xor.witness isa AbstractMatrix
        @test size(hypatia_xor.witness) == (4, 4)
        @test ishermitian(hypatia_xor.witness)
        @test minimum(eigvals(Hermitian(hypatia_xor.witness))) >= -3e-7
        @test diag(hypatia_xor.witness) ≈ ones(4) atol = 3e-7
        @test_throws Base.CanonicalIndexError setindex!(hypatia_xor.witness, 0, 1, 1)
        @test hypatia_xor.optimization_result.status === NonlocalOptQET.OptimizationOptimal
        @test hypatia_xor.optimization_result.termination_status === :optimal
        @test hypatia_xor.optimization_result.primal_status === :feasible_point
        @test hypatia_xor.optimization_result.dual_status === :feasible_point
        @test hypatia_xor.optimization_result.primal_residual <= 2e-7
        @test hypatia_xor.optimization_result.optimizer.configured_optimizer_name ==
            "Hypatia"

        scs_xor = NonlocalOptQET.xor_game_value(probabilities, parity, scs; regime=:quantum)
        @test scs_xor.status === NonlocalOptQET.NonlocalValueNumericalUpperBound
        @test scs_xor.upper_bound ≈ expected_xor atol = 3e-5
        @test scs_xor.optimization_result.status in
            (NonlocalOptQET.OptimizationOptimal, NonlocalOptQET.OptimizationFeasible)
        @test scs_xor.optimization_result.primal_residual <= 3e-5
        @test scs_xor.optimization_result.optimizer.reported_optimizer_name == "SCS"
    end

    @testset "Bell no-signalling and NPA upper bounds" begin
        chsh_fc = [
            0.0 0.0 0.0
            0.0 1.0 1.0
            0.0 1.0 -1.0
        ]
        scenario = NonlocalOptQET.BellScenario(2, 2, 2, 2)
        no_signalling = NonlocalOptQET.bell_inequality_bound(
            chsh_fc, scenario, hypatia; notation=:fc, regime=:no_signalling
        )
        @test no_signalling.status === NonlocalOptQET.NonlocalValueNumericalUpperBound
        @test no_signalling.upper_bound ≈ 4.0 atol = 5e-7
        @test no_signalling.witness isa AbstractArray{<:Real,4}
        @test minimum(no_signalling.witness) >= -3e-7
        @test maximum(abs.(sum(no_signalling.witness; dims=(1, 2)) .- ones(1, 1, 2, 2))) <=
            3e-7
        @test no_signalling.optimization_result.primal_residual <= 2e-7

        npa_bound = NonlocalOptQET.bell_inequality_bound(
            chsh_fc, scenario, hypatia; notation=:fc, regime=:quantum, level=1
        )
        @test npa_bound.status === NonlocalOptQET.NonlocalValueNumericalUpperBound
        @test npa_bound.upper_bound ≈ 2sqrt(2) atol = 1e-6
        @test npa_bound.witness isa AbstractMatrix
        @test size(npa_bound.witness) == (5, 5)
        @test npa_bound.problem.program.metadata.level == "1"
        @test npa_bound.optimization_result.primal_residual <= 2e-7
        @test !npa_bound.certified_upper

        scs_npa = NonlocalOptQET.bell_inequality_bound(
            chsh_fc, scenario, scs; notation=:fc, regime=:quantum, level=1
        )
        @test scs_npa.upper_bound ≈ 2sqrt(2) atol = 5e-5
        @test scs_npa.optimization_result.status in
            (NonlocalOptQET.OptimizationOptimal, NonlocalOptQET.OptimizationFeasible)
        @test scs_npa.optimization_result.primal_residual <= 3e-5
    end

    @testset "NPA membership, exclusion, and failures" begin
        deterministic = zeros(Float64, 2, 2, 2, 2)
        deterministic[1, 1, :, :] .= 1
        local_behavior = NonlocalOptQET.FullProbabilityBehavior(
            deterministic, NonlocalOptQET.BellScenario(2, 2, 2, 2); atol=0, rtol=0
        )
        feasible = NonlocalOptQET.npa_membership(local_behavior, hypatia; level=1)
        @test feasible.status === NonlocalOptQET.NPANumericallyFeasible
        @test feasible.verdict === nothing
        @test !feasible.certified
        @test feasible.moment_matrix isa AbstractMatrix
        @test size(feasible.moment_matrix) == (5, 5)
        @test feasible.optimization_result.status === NonlocalOptQET.OptimizationOptimal
        @test feasible.optimization_result.primal_residual <= 2e-7

        pr_result = NonlocalOptQET.npa_membership(_nonlocal_pr_box(), hypatia; level=1)
        @test pr_result.status in
            (NonlocalOptQET.NPANumericallyInfeasible, NonlocalOptQET.NPABackendFailure)
        @test pr_result.verdict === nothing
        @test !pr_result.certified
        @test pr_result.optimization_result.status in (
            NonlocalOptQET.OptimizationInfeasible,
            NonlocalOptQET.OptimizationUnknown,
            NonlocalOptQET.OptimizationNumericalFailure,
        )

        limited_backend = nonlocal_hypatia_backend(optimizer_options=(iter_limit=0,))
        limited = NonlocalOptQET.npa_membership(local_behavior, limited_backend; level=1)
        @test limited.status === NonlocalOptQET.NPAResourceLimit
        @test limited.optimization_result.status === NonlocalOptQET.OptimizationLimit
        @test limited.optimization_result.termination_status === :iteration_limit

        malformed_backend = NonlocalOptQET.JuMPBackend(
            () -> error("intentional nonlocal optimizer failure");
            optimizer_name="intentional malformed factory",
            allow_densify=true,
        )
        malformed = NonlocalOptQET.npa_membership(
            local_behavior, malformed_backend; level=1
        )
        @test malformed.status === NonlocalOptQET.NPABackendFailure
        @test malformed.optimization_result.status ===
            NonlocalOptQET.OptimizationMalformedBackend
        @test malformed.optimization_result.termination_status === :backend_exception
    end

    @testset "bounded see-saw lower bounds and RNG isolation" begin
        payoff = zeros(Float64, 2, 2, 2, 2)
        for a in 1:2, b in 1:2, x in 1:2, y in 1:2
            payoff[a, b, x, y] = a == b
        end
        game = NonlocalOptQET.NonlocalGame(probabilities, payoff)
        first = NonlocalOptQET.nonlocal_game_lower_bound(
            MersenneTwister(20), 2, game, hypatia; max_iterations=4
        )
        second = NonlocalOptQET.nonlocal_game_lower_bound(
            MersenneTwister(20), 2, game, hypatia; max_iterations=4
        )
        @test first.status === NonlocalOptQET.NonlocalValueNumericalLowerBound
        @test first.lower_bound ≈ 1.0 atol = 2e-6
        @test !first.certified_lower
        @test first.converged
        @test first.iterations <= 4
        @test first.objective_history == second.objective_history
        @test first.lower_bound == second.lower_bound
        @test size(first.alice_operators) == (2, 2)
        @test size(first.bob_measurements) == (2, 2)
        @test_throws Base.CanonicalIndexError setindex!(
            first.alice_operators[1, 1], 0, 1, 1
        )
        @test_throws Base.CanonicalIndexError setindex!(first.objective_history, 0, 1)
        @test all(
            matrix -> minimum(eigvals(Hermitian(matrix))) >= -3e-7, first.alice_operators
        )
        @test all(
            matrix -> minimum(eigvals(Hermitian(matrix))) >= -3e-7, first.bob_measurements
        )
        @test all(y -> norm(sum(first.bob_measurements[:, y]) - I, Inf) <= 3e-7, 1:2)
        @test length(first.optimization_history) == 2 * first.iterations
        @test all(
            result ->
                result.termination_status === :optimal &&
                result.primal_status === :feasible_point &&
                result.dual_status === :feasible_point,
            first.optimization_history,
        )

        Random.seed!(20260730)
        expected_global = rand()
        Random.seed!(20260730)
        _ = NonlocalOptQET.nonlocal_game_lower_bound(
            MersenneTwister(21), 2, game, hypatia; max_iterations=2
        )
        @test rand() == expected_global
        @test_throws ErrorException NonlocalOptQET.MATLABCompat.NonlocalGameLB(
            MersenneTwister(24),
            2,
            probabilities,
            payoff,
            0;
            backend=hypatia,
            max_iterations=4,
            structured=false,
        )

        limited = NonlocalOptQET.nonlocal_game_lower_bound(
            MersenneTwister(22), 2, game, hypatia; max_iterations=1, atol=0, rtol=0
        )
        @test limited.status in (
            NonlocalOptQET.NonlocalValueResourceLimit,
            NonlocalOptQET.NonlocalValueNumericalLowerBound,
        )
        @test limited.lower_bound !== nothing

        malformed_backend = NonlocalOptQET.JuMPBackend(
            () -> error("intentional see-saw failure");
            optimizer_name="intentional malformed factory",
            allow_densify=true,
        )
        failed = NonlocalOptQET.nonlocal_game_lower_bound(
            MersenneTwister(23), 2, game, malformed_backend; max_iterations=2
        )
        @test failed.status === NonlocalOptQET.NonlocalValueBackendFailure
        @test failed.lower_bound === nothing
        @test failed.diagnostics.failed_step === :alice_state
        @test failed.optimization_history[1].termination_status === :backend_exception
    end

    @testset "BCS lower bound and corrected rectangular qubit SDP" begin
        constraints = [[1 0; 0 1], [0 1; 1 0]]
        bcs = NonlocalOptQET.bcs_game_lower_bound(
            MersenneTwister(30),
            2,
            constraints,
            hypatia;
            max_iterations=4,
            coefficient_type=Float64,
        )
        @test bcs.lower_bound !== nothing
        @test bcs.lower_bound >= 0.75 - 2e-6
        @test bcs.lower_bound <= 1 + 2e-6
        @test !bcs.certified_lower
        @test eltype(bcs.game.probabilities) === Float64
        @test bcs.game.tolerance > 0

        rectangular = NonlocalOptQET.bell_inequality_qubit_bound(
            ones(1, 2), zeros(1), zeros(2), [-1.0, 1.0], [-1.0, 1.0], hypatia
        )
        @test rectangular.status === NonlocalOptQET.NonlocalValueNumericalUpperBound
        @test rectangular.upper_bound ≈ 2.0 atol = 2e-6
        @test !rectangular.certified_upper
        @test rectangular.relaxation_state isa AbstractMatrix
        @test size(rectangular.relaxation_state) == (8, 8)
        @test tr(rectangular.relaxation_state) ≈ 1 atol = 3e-7
        @test minimum(eigvals(Hermitian(rectangular.relaxation_state))) >= -3e-7
        @test rectangular.problem.grouped_dimensions == (4, 2)
        @test length(rectangular.problem.ppt_partitions) == 1
        @test rectangular.diagnostics.alice_settings == 1
        @test rectangular.diagnostics.bob_settings == 2
        @test rectangular.diagnostics.upstream_rectangular_loop_corrected
        @test rectangular.optimization_result.primal_residual <= 2e-7

        scs_rectangular = NonlocalOptQET.bell_inequality_qubit_bound(
            ones(1, 2), zeros(1), zeros(2), [-1.0, 1.0], [-1.0, 1.0], scs
        )
        @test scs_rectangular.upper_bound ≈ 2.0 atol = 5e-5
        @test scs_rectangular.optimization_result.status in
            (NonlocalOptQET.OptimizationOptimal, NonlocalOptQET.OptimizationFeasible)
        @test scs_rectangular.optimization_result.primal_residual <= 3e-5
    end
end
