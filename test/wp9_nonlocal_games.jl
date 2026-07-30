using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const NonlocalQET = QuantumEntanglementTools

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
    if !isdefined(NonlocalQET, symbol)
        Base.include(
            NonlocalQET, joinpath(@__DIR__, "..", "src", "nonlocal_games", source_file)
        )
    end
end

struct _NonlocalZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end
Base.size(matrix::_NonlocalZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_NonlocalZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_NonlocalZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_NonlocalZeroBasedMatrix, i::Int, j::Int)
    return matrix.storage[i + 1, j + 1]
end

function _nonlocal_chsh_behavior(::Type{T}=Float64) where {T}
    scenario = NonlocalQET.BellScenario(2, 2, 2, 2)
    probabilities = zeros(T, 2, 2, 2, 2)
    for x in 1:2, y in 1:2
        probabilities[1, 1, x, y] = one(T)
    end
    return NonlocalQET.FullProbabilityBehavior(
        probabilities, scenario; atol=zero(T), rtol=zero(T)
    )
end

@testset "WP9 Bell, nonlocal-game, XOR, BCS, and NPA core" begin
    @testset "scenario, behavior, and notation ownership" begin
        scenario = NonlocalQET.BellScenario(2, 3, 2, 1)
        @test Tuple(scenario) == (2, 3, 2, 1)
        @test NonlocalQET.BellScenario([2, 3, 2, 1]) == scenario
        @test_throws ArgumentError NonlocalQET.BellScenario(true, 2, 2, 2)
        @test_throws ArgumentError NonlocalQET.BellScenario(0, 2, 2, 2)
        @test_throws DimensionMismatch NonlocalQET.BellScenario([2, 2, 2])

        deterministic = zeros(Rational{Int}, 2, 3, 2, 1)
        deterministic[1, 2, :, 1] .= 1
        original = copy(deterministic)
        behavior = NonlocalQET.FullProbabilityBehavior(
            deterministic, scenario; atol=0, rtol=0
        )
        deterministic[1, 2, 1, 1] = 0
        @test behavior.probabilities == original
        @test behavior.minimum_probability == 0
        @test behavior.normalization_residual == 0
        @test behavior.no_signalling_residual == 0
        @test !behavior.boundary
        @test_throws Base.CanonicalIndexError behavior.probabilities[1, 2, 1, 1] = 0
        @test_throws MethodError NonlocalQET.FullProbabilityBehavior(
            scenario, original, 0, 0, 0, 0, false
        )

        cg = NonlocalQET.collins_gisin_behavior(behavior; atol=0, rtol=0)
        round_trip = NonlocalQET.full_probability_behavior(cg; atol=0, rtol=0)
        @test round_trip.probabilities == original
        @test cg.coefficients[1, 1] == 1
        @test size(cg.coefficients) == (3, 3)

        bad_normalization = copy(original)
        bad_normalization[1, 2, 1, 1] = 2
        @test_throws ArgumentError NonlocalQET.FullProbabilityBehavior(
            bad_normalization, scenario; atol=0, rtol=0
        )
        negative = copy(original)
        negative[1, 2, 1, 1] = -1
        negative[2, 2, 1, 1] = 2
        @test_throws DomainError NonlocalQET.FullProbabilityBehavior(
            negative, scenario; atol=0, rtol=0
        )
        signalling_scenario = NonlocalQET.BellScenario(2, 2, 1, 2)
        signalling = zeros(Int, 2, 2, 1, 2)
        signalling[1, 1, 1, 1] = 1
        signalling[2, 1, 1, 2] = 1
        @test_throws ArgumentError NonlocalQET.FullProbabilityBehavior(
            signalling, signalling_scenario; atol=0, rtol=0
        )

        coefficients = reshape(collect(1:12), 2, 3, 2, 1)
        coefficients_copy = copy(coefficients)
        functional = NonlocalQET.BellFunctional(coefficients, scenario; notation=:fp)
        @test functional.coefficients == coefficients_copy
        coefficients .= 0
        @test functional.coefficients == coefficients_copy
        @test NonlocalQET.evaluate_bell_functional(functional, behavior) ==
            sum(coefficients_copy .* original)
        @test NonlocalQET.evaluate_bell_functional(functional, cg) ==
            NonlocalQET.evaluate_bell_functional(functional, behavior)

        binary_scenario = NonlocalQET.BellScenario(2, 2, 2, 3)
        fc = [
            0 0 0 0
            0 1 1 -1
            0 1 -1 1
        ]
        fc_functional = NonlocalQET.BellFunctional(fc, binary_scenario; notation=:fc)
        @test size(fc_functional.coefficients) == (2, 2, 2, 3)
        affine_fc = [
            3 2 -1
            -2 1 1
            1 1 -1
        ]
        affine_functional = NonlocalQET.BellFunctional(
            affine_fc, NonlocalQET.BellScenario(2, 2, 2, 2); notation=:fc
        )
        @test eltype(affine_functional.coefficients) === Rational{BigInt}
        @test NonlocalQET.bell_inequality_bound(
            affine_fc, NonlocalQET.BellScenario(2, 2, 2, 2); notation=:fc, regime=:classical
        ).value == 11

        cg_functional = NonlocalQET.BellFunctional(
            affine_fc, NonlocalQET.BellScenario(2, 2, 2, 2); notation=:cg
        )
        @test eltype(cg_functional.coefficients) === Rational{BigInt}
        @test_throws ArgumentError NonlocalQET.BellFunctional(
            zeros(3, 2), scenario; notation=:fc
        )
        @test_throws ArgumentError NonlocalQET.BellFunctional(
            zeros(2, 3, 2, 1), scenario; notation=:unknown
        )

        zero_based = _NonlocalZeroBasedMatrix(fill(0.25, 2, 2))
        @test_throws ArgumentError NonlocalQET.NonlocalGame(zero_based, ones(2, 2, 2, 2))
    end

    @testset "exact classical values and strategy guards" begin
        probabilities = fill(1 // 4, 2, 2)
        parity = [0 0; 0 1]
        classical = NonlocalQET.xor_game_value(
            probabilities, parity; regime=:classical, atol=0, rtol=0
        )
        @test classical isa NonlocalQET.NonlocalValueResult
        @test classical.status === NonlocalQET.NonlocalValueExact
        @test classical.value == 3 // 4
        @test classical.lower_bound == classical.upper_bound == 3 // 4
        @test classical.exact
        @test classical.certified_lower
        @test classical.certified_upper
        @test classical.certificate_kind === :exhaustive_xor_strategy_reduction
        @test classical.strategy isa NonlocalQET.DeterministicStrategy
        @test classical.diagnostics.bias == 1 // 2
        @test classical.diagnostics.strategies_evaluated == 4
        @test classical.diagnostics.joint_strategy_count == 16
        @test occursin("status=", sprint(show, classical))
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities, parity; max_strategies=3, atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities, parity; max_work=15, atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities, parity; max_entries=3, atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            fill(1 // 3, 2, 2), parity; atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities, [0 2; 0 1]; atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities, parity; regime=:nosignal, atol=0, rtol=0
        )

        chsh_fc = [
            0.0 0.0 0.0
            0.0 1.0 1.0
            0.0 1.0 -1.0
        ]
        scenario = NonlocalQET.BellScenario(2, 2, 2, 2)
        bell = NonlocalQET.bell_inequality_bound(
            chsh_fc, scenario; notation=:fc, regime=:classical
        )
        @test bell.status === NonlocalQET.NonlocalValueExact
        @test bell.value ≈ 2.0
        @test bell.strategy isa NonlocalQET.DeterministicStrategy
        @test bell.diagnostics.strategies_evaluated == 4
        @test bell.diagnostics.joint_strategy_count == 16
        @test_throws ArgumentError NonlocalQET.bell_inequality_bound(
            chsh_fc, scenario; notation=:fc, regime=:classical, max_strategies=3
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_bound(
            chsh_fc, scenario; notation=:fc, regime=:classical, max_work=31
        )

        payoff = zeros(Int, 2, 2, 2, 2)
        for a in 1:2, b in 1:2, x in 1:2, y in 1:2
            payoff[a, b, x, y] = a == b
        end
        game = NonlocalQET.NonlocalGame(probabilities, payoff; atol=0, rtol=0)
        game_value = NonlocalQET.nonlocal_game_value(game)
        @test game_value.value == 1
        @test game_value.exact
        @test game_value.strategy.alice_outputs == [1, 1]
        @test game_value.strategy.bob_outputs == [1, 1]
        @test game.probabilities isa AbstractMatrix{Rational{Int}}
        @test game.payoff isa AbstractArray{Int,4}
        @test_throws Base.CanonicalIndexError game.probabilities[1, 1] = 0
        @test_throws Base.CanonicalIndexError game.payoff[1, 1, 1, 1] = 0
        @test_throws Base.CanonicalIndexError game_value.strategy.alice_outputs[1] = 2

        sparse_game = NonlocalQET.NonlocalGame(
            sparse(probabilities), payoff; atol=0, rtol=0
        )
        @test issparse(sparse_game.probabilities)
        @test NonlocalQET.nonlocal_game_value(sparse_game).value == 1
        @test_throws DomainError NonlocalQET.NonlocalGame(
            [1.0 -0.1; 0.0 0.1], payoff; atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.NonlocalGame(
            fill(0.2, 2, 2), payoff; atol=0, rtol=0
        )
    end

    @testset "BCS conversion and missing-dispatch reconstruction" begin
        equality = [1 0; 0 1]
        inequality = [0 1; 1 0]
        constraints = [equality, inequality]
        bcs = NonlocalQET.BCSGame(constraints)
        @test bcs.variable_count == 2
        game = NonlocalQET.nonlocal_game(bcs)
        @test Tuple(game.scenario) == (4, 2, 2, 2)
        @test sum(game.probabilities) == 1
        @test size(game.payoff) == (4, 2, 2, 2)
        @test game.probabilities[1, 1] == 1 // 4
        @test game.probabilities[2, 2] == 1 // 4
        # The pinned helper labels Alice's assignments with the first active
        # variable as the most-significant bit (`00`, `01`, `10`, `11`).
        @test game.payoff[2, 1, 2, 1] == 1
        @test game.payoff[2, 2, 2, 2] == 1
        @test game.payoff[3, 2, 2, 1] == 1
        @test game.payoff[3, 1, 2, 2] == 1
        partial = NonlocalQET.nonlocal_game(NonlocalQET.BCSGame([[1 1; 0 0], [1 0; 1 0]]))
        @test partial.probabilities == [1 // 2 0; 0 1 // 2]
        @test vec(partial.payoff) == [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
        classical = NonlocalQET.bcs_game_value(bcs; regime=:classical)
        @test classical.status === NonlocalQET.NonlocalValueExact
        @test classical.value == 3 // 4
        @test classical.certificate_kind === :exhaustive_deterministic_strategies
        constraints[1][1, 1] = 0
        @test bcs.constraints[1][1, 1] == 1
        @test_throws Base.CanonicalIndexError bcs.constraints[1][1, 1] = 0

        @test_throws ArgumentError NonlocalQET.BCSGame(Any[])
        @test_throws DimensionMismatch NonlocalQET.BCSGame([
            ones(Int, 2, 2), ones(Int, 2, 2, 2)
        ])
        @test_throws ArgumentError NonlocalQET.BCSGame([[1 2; 0 1]])
        inactive = NonlocalQET.BCSGame([ones(Int, 2, 2)])
        @test_throws ArgumentError NonlocalQET.nonlocal_game(inactive)
        @test_throws ArgumentError NonlocalQET.nonlocal_game(bcs; max_entries=10)
    end

    @testset "canonical NPA words and bounded models" begin
        A10 = NonlocalQET.NPALetter(:alice, 1, 1)
        A20 = NonlocalQET.NPALetter(:alice, 2, 1)
        B10 = NonlocalQET.NPALetter(:bob, 1, 1)
        @test_throws ArgumentError NonlocalQET.NPALetter(:charlie, 1, 1)
        @test_throws MethodError NonlocalQET.NPALetter(0xff, -1, -1)
        @test NonlocalQET.npa_word((A10, A10)) == NonlocalQET.NPAWord((A10,))
        @test NonlocalQET.npa_word((B10, A10)) == NonlocalQET.NPAWord((A10, B10))
        @test NonlocalQET.npa_word((A10, NonlocalQET.NPALetter(:alice, 1, 2))) === nothing
        @test NonlocalQET.npa_word((A10, A20, A10)) == NonlocalQET.NPAWord((A10, A20, A10))

        integer_level = NonlocalQET.NPALevel(2)
        @test integer_level.base_level == 2
        @test integer_level.maximum_length == 2
        intermediate = NonlocalQET.NPALevel("1+ab+aab")
        @test intermediate.base_level == 1
        @test intermediate.maximum_length == 3
        @test intermediate.components == ((1, 1), (2, 1))
        @test_throws ArgumentError NonlocalQET.NPALevel("1+ac")
        @test_throws ArgumentError NonlocalQET.NPALevel(-1)
        @test_throws MethodError NonlocalQET.NPALevel("forged", -1, -1, ())

        behavior = _nonlocal_chsh_behavior()
        level_zero = NonlocalQET.npa_membership(behavior; level=0)
        @test level_zero.status === NonlocalQET.NPABasicConditionsSatisfied
        @test level_zero.verdict === true
        @test level_zero.certified
        @test level_zero.certificate_kind ===
            :exact_probability_and_no_signalling_conditions
        @test level_zero.problem === nothing

        missing = NonlocalQET.npa_membership(behavior; level=1)
        @test missing.status === NonlocalQET.NPABackendUnavailable
        @test missing.verdict === nothing
        @test !missing.certified
        @test missing.problem isa NonlocalQET.NPAProblem
        @test length(missing.problem.words) == 5
        @test missing.problem.program.variable_count == 25
        @test missing.problem.program.sense === :feasibility
        @test missing.optimization_result.status ===
            NonlocalQET.OptimizationBackendUnavailable
        @test occursin("not an exact", missing.message) == false
        @test occursin("no optimization backend", missing.optimization_result.message)

        cg = NonlocalQET.collins_gisin_behavior(behavior; atol=0, rtol=0)
        direct = NonlocalQET.npa_membership(cg.coefficients, [2, 2, 2, 2]; level=1)
        @test direct.status === NonlocalQET.NPABackendUnavailable
        @test direct.problem.behavior.coefficients == cg.coefficients
        @test_throws ArgumentError NonlocalQET.npa_membership(
            behavior; level=2, max_word_generation_work=10
        )
        @test_throws ArgumentError NonlocalQET.npa_membership(
            behavior; level=2, max_words=5
        )
        @test_throws ArgumentError NonlocalQET.npa_membership(
            behavior; level=1, limits=NonlocalQET.OptimizationLimits(max_variables=24)
        )
        @test_throws MethodError NonlocalQET.NPAResult(
            level_zero.status,
            level_zero.verdict,
            level_zero.certified,
            level_zero.certificate_kind,
            level_zero.level,
            level_zero.problem,
            level_zero.moment_matrix,
            level_zero.optimization_result,
            level_zero.residuals,
            level_zero.tolerance,
            level_zero.message,
        )
    end

    @testset "solver-neutral Bell and XOR problems" begin
        @test NonlocalQET._nonlocal_solver_upper((
            objective_bound=nothing, dual_objective_value=nothing
        )) === nothing
        @test NonlocalQET._nonlocal_solver_upper((
            objective_bound=2.0, dual_objective_value=3.0
        )) == 2.0
        @test NonlocalQET._nonlocal_solver_upper((
            objective_bound=nothing, dual_objective_value=3.0
        )) == 3.0
        probabilities = fill(0.25, 2, 2)
        parity = [0 0; 0 1]
        quantum = NonlocalQET.xor_game_value(probabilities, parity; regime=:quantum)
        @test quantum.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test quantum.value === nothing
        @test quantum.upper_bound === nothing
        @test quantum.problem.program.name === :xor_quantum_value
        @test quantum.problem.program.variable_count == 16
        @test length(quantum.problem.program.equalities) == 4
        @test length(quantum.problem.program.psd_constraints) == 1
        @test quantum.optimization_result.termination_status === :backend_unavailable
        @test_throws ArgumentError NonlocalQET.xor_game_value(
            probabilities,
            parity;
            regime=:quantum,
            limits=NonlocalQET.OptimizationLimits(max_variables=15),
        )

        scenario = NonlocalQET.BellScenario(2, 2, 2, 2)
        coefficients = zeros(Float64, 2, 2, 2, 2)
        coefficients[1, 1, :, :] .= 1
        no_signalling = NonlocalQET.bell_inequality_bound(
            coefficients, scenario; notation=:fp, regime=:no_signalling
        )
        @test no_signalling.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test no_signalling.problem.program.sense === :maximize
        @test no_signalling.problem.program.variable_count == 16
        @test length(no_signalling.problem.program.intervals) == 16
        @test length(no_signalling.problem.program.equalities) == 12
        @test isempty(no_signalling.problem.program.psd_constraints)
        @test NonlocalQET.primal_residual(
            no_signalling.problem.program,
            no_signalling.problem.program.known_feasible_point;
            allow_densify=true,
        ) == 0

        npa = NonlocalQET.bell_inequality_bound(
            coefficients, scenario; notation=:fp, regime=:quantum, level=1
        )
        @test npa.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test npa.problem isa NonlocalQET.NPAProblem
        @test npa.problem.objective_kind === :bell_upper_bound
        @test npa.problem.program.sense === :maximize
        @test npa.problem.program.metadata.level == "1"

        @test_throws ArgumentError NonlocalQET.bell_inequality_bound(
            coefficients, scenario; notation=:fp, regime=:quantum, level=0
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_bound(
            coefficients, scenario; notation=:fp, regime=:other
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_bound(
            Rational{Int}.(coefficients), scenario; notation=:fp, regime=:no_signalling
        )
    end

    @testset "explicit RNG and bounded lower-bound paths" begin
        probabilities = fill(0.25, 2, 2)
        payoff = zeros(Float64, 2, 2, 2, 2)
        for a in 1:2, b in 1:2, x in 1:2, y in 1:2
            payoff[a, b, x, y] = a == b
        end
        game = NonlocalQET.NonlocalGame(probabilities, payoff)

        zero_rng = MersenneTwister(40)
        zero_control = copy(zero_rng)
        zero_work = NonlocalQET.nonlocal_game_lower_bound(
            zero_rng, 2, game; max_iterations=0
        )
        @test zero_work.status === NonlocalQET.NonlocalValueResourceLimit
        @test zero_work.lower_bound === nothing
        @test zero_work.iterations == 0
        @test rand(zero_rng) == rand(zero_control)

        missing_rng = MersenneTwister(41)
        missing_control = copy(missing_rng)
        missing = NonlocalQET.nonlocal_game_lower_bound(
            missing_rng, 2, game; max_iterations=3
        )
        @test missing.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test missing.lower_bound === nothing
        @test rand(missing_rng) == rand(missing_control)
        @test_throws ArgumentError NonlocalQET.nonlocal_game_lower_bound(
            MersenneTwister(42), 0, game
        )
        @test_throws ArgumentError NonlocalQET.nonlocal_game_lower_bound(
            MersenneTwister(42), 2, game; max_iterations=-1
        )
        resource_rng = MersenneTwister(45)
        resource_control = copy(resource_rng)
        @test_throws ArgumentError NonlocalQET.nonlocal_game_lower_bound(
            resource_rng,
            2,
            game;
            backend=NonlocalQET.JuMPBackend(
                () -> nothing; optimizer_name="resource sentinel"
            ),
            max_dense_entries=1,
        )
        @test rand(resource_rng) == rand(resource_control)
        model_rng = MersenneTwister(46)
        model_control = copy(model_rng)
        @test_throws ArgumentError NonlocalQET.nonlocal_game_lower_bound(
            model_rng,
            2,
            game;
            backend=NonlocalQET.JuMPBackend(
                () -> nothing; optimizer_name="model-resource sentinel"
            ),
            limits=NonlocalQET.OptimizationLimits(max_model_entries=1),
        )
        @test rand(model_rng) == rand(model_control)

        exact_game = NonlocalQET.NonlocalGame(
            fill(1 // 4, 2, 2), Int.(payoff); atol=0, rtol=0
        )
        @test_throws ArgumentError NonlocalQET.nonlocal_game_lower_bound(
            MersenneTwister(43),
            2,
            exact_game;
            max_iterations=1,
            backend=NonlocalQET.JuMPBackend(
                () -> nothing; optimizer_name="type-check sentinel"
            ),
        )

        bcs_constraints = [[1 0; 0 1], [0 1; 1 0]]
        bcs_rng = MersenneTwister(44)
        bcs_control = copy(bcs_rng)
        bcs_missing = NonlocalQET.bcs_game_lower_bound(
            bcs_rng, 2, bcs_constraints; max_iterations=1
        )
        @test bcs_missing.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test rand(bcs_rng) == rand(bcs_control)
    end

    @testset "corrected rectangular qubit Bell model" begin
        joint = [1.0 -1.0 0.5; 0.25 1.5 -0.75]
        alice = [0.5, -0.25]
        bob = [0.1, 0.2, -0.3]
        result = NonlocalQET.bell_inequality_qubit_bound(
            joint, alice, bob, [-1.0, 1.0], [-1.0, 1.0]
        )
        @test result.status === NonlocalQET.NonlocalValueBackendUnavailable
        @test result.upper_bound === nothing
        @test result.problem.grouped_dimensions == (4, 2, 2, 2)
        @test result.problem.program.metadata.alice_settings == 2
        @test result.problem.program.metadata.bob_settings == 3
        @test result.problem.program.metadata.upstream_rectangular_loop_corrected
        @test result.diagnostics.upstream_rectangular_loop_corrected
        @test result.diagnostics.alice_settings == 2
        @test result.diagnostics.bob_settings == 3
        @test length(result.problem.ppt_partitions) == 7
        @test length(result.problem.program.psd_constraints) == 8
        @test result.problem.program.variable_count == 32^2
        @test_throws Base.CanonicalIndexError result.problem.alice_values[1] = 0
        @test NonlocalQET.primal_residual(
            result.problem.program,
            result.problem.program.known_feasible_point;
            allow_densify=true,
        ) <= 1e-14
        @test_throws DimensionMismatch NonlocalQET.bell_inequality_qubit_bound(
            joint, [1.0], bob, [-1.0, 1.0], [-1.0, 1.0]
        )
        @test_throws DimensionMismatch NonlocalQET.bell_inequality_qubit_bound(
            joint, alice, bob, [-1.0, 0.0, 1.0], [-1.0, 1.0]
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_qubit_bound(
            joint, alice, bob, [-1.0, 1.0], [-1.0, 1.0]; max_dimension=31
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_qubit_bound(
            joint, alice, bob, [-1.0, 1.0], [-1.0, 1.0]; max_ppt_constraints=6
        )
        @test_throws ArgumentError NonlocalQET.bell_inequality_qubit_bound(
            joint,
            alice,
            bob,
            [-1.0, 1.0],
            [-1.0, 1.0];
            limits=NonlocalQET.OptimizationLimits(max_variables=1023),
        )
    end
end
