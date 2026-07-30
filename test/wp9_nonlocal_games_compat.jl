using QuantumEntanglementTools
using Random
using Test

const NonlocalCompatQET = QuantumEntanglementTools

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
    if !isdefined(NonlocalCompatQET, symbol)
        Base.include(
            NonlocalCompatQET,
            joinpath(@__DIR__, "..", "src", "nonlocal_games", source_file),
        )
    end
end

module _NonlocalCompatHarness
using QuantumEntanglementTools:
    BellScenario,
    CollinsGisinBehavior,
    NoOptimizationBackend,
    bcs_game_lower_bound,
    bcs_game_value,
    bell_inequality_bound,
    bell_inequality_qubit_bound,
    nonlocal_game_lower_bound,
    npa_membership,
    xor_game_value

include(joinpath(@__DIR__, "..", "src", "compat", "nonlocal_games.jl"))
end

const NonlocalCompat = _NonlocalCompatHarness

@testset "WP9 nonlocal-game MATLAB compatibility contracts" begin
    probabilities = fill(0.25, 2, 2)
    parity = [0 0; 0 1]
    xor_structured = NonlocalCompat.XORGameValue(probabilities, parity)
    @test xor_structured.status === NonlocalCompatQET.NonlocalValueExact
    @test xor_structured.value == 0.75
    @test NonlocalCompat.XORGameValue(
        probabilities, parity, "CLASSICAL"; structured=false
    ) == 0.75
    xor_quantum = NonlocalCompat.XORGameValue(probabilities, parity, "quantum")
    @test xor_quantum.status === NonlocalCompatQET.NonlocalValueBackendUnavailable
    @test_throws ErrorException NonlocalCompat.XORGameValue(
        probabilities, parity, "quantum"; structured=false
    )
    @test_throws ArgumentError NonlocalCompat.XORGameValue(
        probabilities, parity, "no-signal"
    )

    deterministic = zeros(Float64, 2, 2, 2, 2)
    deterministic[1, 1, :, :] .= 1
    behavior = NonlocalCompatQET.FullProbabilityBehavior(
        deterministic, NonlocalCompatQET.BellScenario(2, 2, 2, 2); atol=0, rtol=0
    )
    cg = NonlocalCompatQET.collins_gisin_behavior(behavior; atol=0, rtol=0)
    npa_zero = NonlocalCompat.NPAHierarchy(cg.coefficients, [2, 2, 2, 2], 0)
    @test npa_zero.status === NonlocalCompatQET.NPABasicConditionsSatisfied
    @test NonlocalCompat.NPAHierarchy(cg.coefficients, [2, 2, 2, 2], 0; structured=false) ==
        1
    npa_one = NonlocalCompat.NPAHierarchy(cg.coefficients, [2, 2, 2, 2], 1)
    @test npa_one.status === NonlocalCompatQET.NPABackendUnavailable
    @test_throws ErrorException NonlocalCompat.NPAHierarchy(
        cg.coefficients, [2, 2, 2, 2], 1; structured=false
    )

    chsh_fc = [
        0.0 0.0 0.0
        0.0 1.0 1.0
        0.0 1.0 -1.0
    ]
    bell = NonlocalCompat.BellInequalityMax(chsh_fc, [2, 2, 2, 2], "FC")
    @test bell.status === NonlocalCompatQET.NonlocalValueExact
    @test bell.value ≈ 2.0
    @test NonlocalCompat.BellInequalityMax(
        chsh_fc, [2, 2, 2, 2], :fc, :classical, 1; structured=false
    ) ≈ 2.0
    @test_throws ArgumentError NonlocalCompat.BellInequalityMax(chsh_fc, [2, 2, 2, 2], :bad)

    rectangular = NonlocalCompat.BellInequalityMaxQubits(
        ones(1, 2), zeros(1), zeros(2), [-1.0, 1.0], [-1.0, 1.0]
    )
    @test rectangular.status === NonlocalCompatQET.NonlocalValueBackendUnavailable
    @test rectangular.problem.program.metadata.alice_settings == 1
    @test rectangular.problem.program.metadata.bob_settings == 2
    @test_throws ErrorException NonlocalCompat.BellInequalityMaxQubits(
        ones(1, 2), zeros(1), zeros(2), [-1.0, 1.0], [-1.0, 1.0]; structured=false
    )

    constraints = [[1 0; 0 1], [0 1; 1 0]]
    bcs = NonlocalCompat.BCSGameValue(constraints)
    @test bcs.status === NonlocalCompatQET.NonlocalValueExact
    @test bcs.value == 0.75
    @test NonlocalCompat.BCSGameValue(constraints, "classical", 1; structured=false) == 0.75
    bcs_quantum = NonlocalCompat.BCSGameValue(constraints, "quantum")
    @test bcs_quantum.status === NonlocalCompatQET.NonlocalValueBackendUnavailable
    @test_throws ErrorException NonlocalCompat.BCSGameValue(
        constraints, "quantum"; structured=false
    )

    payoff = zeros(Float64, 2, 2, 2, 2)
    for a in 1:2, b in 1:2, x in 1:2, y in 1:2
        payoff[a, b, x, y] = a == b
    end
    rng = MersenneTwister(5)
    control = copy(rng)
    missing = NonlocalCompat.NonlocalGameLB(
        rng, 2, probabilities, payoff, 0; max_iterations=1
    )
    @test missing.status === NonlocalCompatQET.NonlocalValueBackendUnavailable
    @test rand(rng) == rand(control)
    @test_throws ErrorException NonlocalCompat.NonlocalGameLB(
        MersenneTwister(6), 2, probabilities, payoff; structured=false, max_iterations=1
    )
    @test_throws ArgumentError NonlocalCompat.NonlocalGameLB(
        MersenneTwister(7), 2, probabilities, payoff, 2
    )

    bcs_rng = MersenneTwister(8)
    bcs_control = copy(bcs_rng)
    bcs_missing = NonlocalCompat.BCSGameLB(bcs_rng, 2, constraints, false; max_iterations=1)
    @test bcs_missing.status === NonlocalCompatQET.NonlocalValueBackendUnavailable
    @test rand(bcs_rng) == rand(bcs_control)
    @test_throws ErrorException NonlocalCompat.BCSGameLB(
        MersenneTwister(9), 2, constraints; structured=false, max_iterations=1
    )
end
