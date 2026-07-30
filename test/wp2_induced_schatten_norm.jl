using LinearAlgebra
using Random
using SparseArrays

const QETInducedSchatten = QuantumEntanglementTools
const CompatInducedSchatten = QuantumEntanglementTools.MATLABCompat

@testset "WP2 induced Schatten exact Frobenius branch" begin
    identity_map = QETInducedSchatten.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
    result = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(1), identity_map, 2
    )
    @test result isa QETInducedSchatten.InducedSchattenNormResult
    @test result.value ≈ 1
    @test result.p == result.q == 2
    @test result.bound_kind === :exact
    @test result.status === :exact
    @test result.exact
    @test result.converged
    @test result.iterations == 0
    @test result.iteration_residual === nothing
    @test result.objective_history == [result.value]
    @test result.witness_input_norm ≈ 1
    @test result.witness_output_norm ≈ result.value
    @test result.normalization_residual <= 8eps(Float64)
    @test result.value_residual <= 8eps(Float64)
    @test occursin("bound_kind=exact", sprint(show, result))

    scale = 2.5
    scaled_map = QETInducedSchatten.KrausRepresentation([
        sqrt(scale) * Matrix{Float64}(I, 2, 2)
    ])
    scaled = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(2), scaled_map, 2
    )
    @test scaled.value ≈ scale
    @test QETInducedSchatten.apply_channel(scaled.witness, scaled_map) |>
        matrix -> norm(matrix) ≈ scale

    rectangular_kraus = [
        1.0 0.0
        0.0 2.0
        0.5 -1.0
    ]
    rectangular_map = QETInducedSchatten.KrausRepresentation([rectangular_kraus])
    rectangular = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(3), rectangular_map, 2
    )
    transfer = QETInducedSchatten.superoperator_matrix(rectangular_map)
    @test rectangular.value ≈ first(svdvals(transfer))
    @test size(rectangular.witness) == (2, 2)
    @test norm(QETInducedSchatten.apply_channel(rectangular.witness, rectangular_map)) ≈
        rectangular.value

    choi_map = QETInducedSchatten.choi_representation(rectangular_map)
    transfer_map = QETInducedSchatten.superoperator_representation(rectangular_map)
    @test QETInducedSchatten.induced_schatten_lower_bound(MersenneTwister(4), choi_map, 2).value ≈
        rectangular.value
    @test QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(5), transfer_map, 2
    ).value ≈ rectangular.value

    float32_map = QETInducedSchatten.KrausRepresentation([Matrix{Float32}(I, 2, 2)])
    float32_result = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(6), float32_map, Float32(2)
    )
    @test float32_result.value isa Float32
    @test eltype(float32_result.witness) == Float32

    exact_rng = MersenneTwister(0x5348415454454e)
    exact_control = copy(exact_rng)
    QETInducedSchatten.induced_schatten_lower_bound(exact_rng, identity_map, 2)
    @test rand(exact_rng, UInt64) == rand(exact_control, UInt64)

    sparse_map = QETInducedSchatten.KrausRepresentation([sparse(Matrix{Float64}(I, 2, 2))])
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(7), sparse_map, 2
    )
    sparse_exact = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(8), sparse_map, 2; allow_densify=true
    )
    @test sparse_exact.value ≈ 1
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(9), identity_map, 2; max_dense_entries=15
    )
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(10), identity_map, 2; max_work=1
    )
end

@testset "WP2 induced Schatten witnessed lower bounds" begin
    identity_map = QETInducedSchatten.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
    initial = [
        1.0 2.0
        -3.0 4.0
    ]

    three_to_two = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(11),
        identity_map,
        3;
        q=2,
        initial_matrix=initial,
        tolerance=1e-13,
        max_iterations=100,
    )
    @test three_to_two.bound_kind === :lower_bound
    @test !three_to_two.exact
    @test three_to_two.status === :converged_lower_bound
    @test three_to_two.converged
    @test three_to_two.value ≈ 2.0^(1 / 6) rtol = 1e-10
    @test three_to_two.witness_input_norm ≈ 1
    @test three_to_two.witness_output_norm ≈ three_to_two.value
    @test three_to_two.iteration_residual <= three_to_two.tolerance
    @test length(three_to_two.objective_history) == three_to_two.iterations + 1
    @test all(diff(three_to_two.objective_history) .>= -64eps(Float64))

    infinity_to_one = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(12),
        identity_map,
        Inf;
        q=1,
        initial_matrix=initial,
        tolerance=0,
        max_iterations=10,
    )
    @test infinity_to_one.value ≈ 2
    @test infinity_to_one.witness_input_norm ≈ 1
    @test infinity_to_one.witness_output_norm ≈ 2

    one_to_infinity = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(13),
        identity_map,
        1;
        q=Inf,
        initial_matrix=initial,
        max_iterations=10,
    )
    @test one_to_infinity.value ≈ 1
    @test one_to_infinity.witness_input_norm ≈ 1

    limited = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(14), identity_map, 3; q=2, initial_matrix=initial, max_iterations=0
    )
    @test limited.status === :iteration_limit
    @test limited.iterations == 0
    @test limited.iteration_residual === nothing
    @test length(limited.objective_history) == 1

    initialization_work = BigInt(24)
    work_limited = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(15),
        identity_map,
        3;
        q=2,
        initial_matrix=initial,
        max_work=initialization_work,
    )
    @test work_limited.status === :work_limit
    @test work_limited.work_used == initialization_work
    @test work_limited.iterations == 0
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(16),
        identity_map,
        3;
        q=2,
        initial_matrix=initial,
        max_work=initialization_work - 1,
    )

    projection = [
        1.0 0.0
        0.0 0.0
    ]
    projection_map = QETInducedSchatten.KrausRepresentation([projection])
    stationary = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(17), projection_map, 1; q=1, initial_matrix=[0.0 0.0; 0.0 1.0]
    )
    @test stationary.status === :stationary_zero
    @test stationary.value == 0
    @test !stationary.exact
    @test !stationary.converged
    @test occursin("not proof", stationary.message)

    first_random = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(0x524e4753), identity_map, 3; q=2, max_iterations=5
    )
    second_random = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(0x524e4753), identity_map, 3; q=2, max_iterations=5
    )
    @test first_random.value == second_random.value
    @test first_random.witness == second_random.witness
    @test first_random.objective_history == second_random.objective_history

    global_control = copy(Random.default_rng())
    QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(0x474c4f42414c), identity_map, 3; q=2, max_iterations=2
    )
    @test rand(Random.default_rng(), UInt64) == rand(global_control, UInt64)

    supplied_rng = MersenneTwister(0x535550504c494544)
    supplied_control = copy(supplied_rng)
    QETInducedSchatten.induced_schatten_lower_bound(
        supplied_rng, identity_map, 3; q=2, initial_matrix=initial, max_iterations=2
    )
    @test rand(supplied_rng, UInt64) == rand(supplied_control, UInt64)

    original = copy(initial)
    QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(18), identity_map, 3; q=2, initial_matrix=initial, max_iterations=2
    )
    @test initial == original
end

@testset "WP2 InducedSchattenNorm compatibility and validation" begin
    identity_matrix = Matrix{Float64}(I, 2, 2)
    kraus = [identity_matrix]
    identity_map = QETInducedSchatten.KrausRepresentation(kraus)
    initial = [
        1.0 2.0
        -3.0 4.0
    ]

    compatibility_exact = CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(20), kraus, "FRO"
    )
    @test compatibility_exact.status === :exact
    @test compatibility_exact.value ≈ 1

    compatibility_lower = CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(21), kraus, 3, 2, [2, 2], 1e-12, initial; max_iterations=100
    )
    native_lower = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(99),
        identity_map,
        3;
        q=2,
        tolerance=1e-12,
        initial_matrix=initial,
        max_iterations=100,
    )
    @test compatibility_lower.value == native_lower.value
    @test compatibility_lower.witness == native_lower.witness

    sentinel = CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(22), kraus, 3, 2, [2, 2], 1e-8, -1; max_iterations=2
    )
    omitted = CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(22), kraus, 3, 2, [2, 2], 1e-8; max_iterations=2
    )
    @test sentinel.value == omitted.value
    @test sentinel.witness == omitted.witness

    raw_choi = QETInducedSchatten.choi_matrix(identity_map)
    choi_result = CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(23), raw_choi, 2, 2, [2, 2]
    )
    @test choi_result.value ≈ 1
    @test_throws MethodError CompatInducedSchatten.InducedSchattenNorm(kraus, 2)
    @test_throws DimensionMismatch CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(24), kraus, 2, 2, [3, 3]
    )
    @test_throws DimensionMismatch CompatInducedSchatten.InducedSchattenNorm(
        MersenneTwister(25), kraus, 3, 2, [2, 2], 1e-8, ones(3, 3)
    )

    integer_map = QETInducedSchatten.KrausRepresentation([Matrix{Int}(I, 2, 2)])
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(30), integer_map, 2
    )
    nonsquare_map = QETInducedSchatten.OperatorSumRepresentation([ones(2, 2)], [ones(3, 2)])
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(31), nonsquare_map, 2
    )

    for invalid_order in (0, -1, NaN, -Inf, true)
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(32), identity_map, invalid_order
        )
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(33), identity_map, 2; q=invalid_order
        )
    end
    for invalid_tolerance in (-1, Inf, NaN, true)
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(34), identity_map, 3; q=2, tolerance=invalid_tolerance
        )
    end
    for invalid_iterations in (-1, true, 1.5)
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(35), identity_map, 3; q=2, max_iterations=invalid_iterations
        )
    end
    for invalid_work in (0, -1, true, 1.5)
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(36), identity_map, 3; q=2, max_work=invalid_work
        )
    end
    for invalid_entries in (0, -1, true, 1.5)
        @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
            MersenneTwister(37), identity_map, 3; q=2, max_dense_entries=invalid_entries
        )
    end
    @test_throws DomainError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(38), identity_map, 3; q=2, initial_matrix=zeros(2, 2)
    )
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(39), identity_map, 3; q=2, initial_matrix=[1.0 NaN; 0 1]
    )
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(40), identity_map, 3; q=2, initial_matrix=Float32[1 0; 0 1]
    )
    @test_throws ArgumentError QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(41), identity_map, 3; q=2, initial_matrix=sparse(initial)
    )
    sparse_initial = QETInducedSchatten.induced_schatten_lower_bound(
        MersenneTwister(42),
        identity_map,
        3;
        q=2,
        initial_matrix=sparse(initial),
        allow_densify=true,
        max_iterations=2,
    )
    @test sparse_initial.value > 0
end
