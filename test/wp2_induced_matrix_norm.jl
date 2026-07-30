using LinearAlgebra
using Random
using SparseArrays

const QETInduced = QuantumEntanglementTools
const CompatInduced = QuantumEntanglementTools.MATLABCompat

@testset "WP2 induced matrix norm exact branches" begin
    matrix = [
        1.0 -2.0
        3.0 4.0
        -5.0 1.0
    ]

    one_result = QETInduced.induced_matrix_norm(MersenneTwister(1), matrix, 1)
    @test one_result isa QETInduced.InducedMatrixNormResult
    @test one_result.value == 9
    @test one_result.p == one_result.q == 1
    @test one_result.bound_kind === :exact
    @test one_result.status === :exact
    @test one_result.exact
    @test one_result.converged
    @test one_result.witness == [1.0, 0.0]
    @test one_result.witness_input_norm == 1
    @test one_result.witness_output_norm == one_result.value
    @test one_result.normalization_residual == 0
    @test one_result.value_residual == 0
    @test one_result.iteration_residual === nothing
    @test one_result.iterations == 0
    @test one_result.work_used > 0
    @test occursin("bound_kind=exact", sprint(show, one_result))

    two_result = QETInduced.induced_matrix_norm(MersenneTwister(2), matrix, 2)
    @test two_result.value ≈ opnorm(matrix, 2)
    @test two_result.bound_kind === :exact
    @test two_result.status === :exact
    @test two_result.exact
    @test norm(two_result.witness) ≈ 1
    @test norm(matrix * two_result.witness) ≈ two_result.value
    @test two_result.value_residual <= 16eps(Float64) * two_result.value

    infinity_result = QETInduced.induced_matrix_norm(MersenneTwister(3), matrix, Inf)
    @test infinity_result.value == 7
    @test infinity_result.bound_kind === :exact
    @test infinity_result.witness == [1.0, 1.0]
    @test norm(infinity_result.witness, Inf) == 1
    @test norm(matrix * infinity_result.witness, Inf) == 7

    complex_matrix = ComplexF64[
        1+im -2im
        -3+4im 1-im
    ]
    complex_infinity = QETInduced.induced_matrix_norm(
        MersenneTwister(4), complex_matrix, Inf
    )
    @test complex_infinity.value ≈ 5 + sqrt(2)
    @test norm(complex_infinity.witness, Inf) == 1
    @test norm(complex_matrix * complex_infinity.witness, Inf) ≈ complex_infinity.value
    @test complex_infinity.value_residual <= 16eps(Float64) * complex_infinity.value

    float32_matrix = Float32[1 -2; 3 4; -5 1]
    float32_one = QETInduced.induced_matrix_norm(
        MersenneTwister(5), float32_matrix, Float32(1)
    )
    float32_two = QETInduced.induced_matrix_norm(
        MersenneTwister(6), float32_matrix, Float32(2)
    )
    @test float32_one.value isa Float32
    @test eltype(float32_one.witness) == Float32
    @test float32_two.value isa Float32
    @test eltype(float32_two.witness) == Float32

    exact_rng = MersenneTwister(0x4558414354)
    exact_control = copy(exact_rng)
    QETInduced.induced_matrix_norm(exact_rng, matrix, 1)
    @test rand(exact_rng, UInt64) == rand(exact_control, UInt64)

    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(7), matrix, 1; max_work=1
    )
end

@testset "WP2 induced matrix norm witnessed lower bounds" begin
    diagonal = Diagonal([3.0, 2.0])
    exact_three_to_two = (3.0^6 + 2.0^6)^(1 / 6)
    lower = QETInduced.induced_matrix_norm(
        MersenneTwister(10),
        diagonal,
        3;
        q=2,
        initial_vector=[1.0, 1.0],
        tolerance=1e-13,
        max_iterations=200,
    )
    @test lower.bound_kind === :lower_bound
    @test !lower.exact
    @test lower.status === :converged_lower_bound
    @test lower.converged
    @test 0 < lower.iterations <= 200
    @test lower.iteration_residual <= lower.tolerance
    @test norm(lower.witness, 3) ≈ 1
    @test lower.witness_output_norm ≈ norm(diagonal * lower.witness, 2)
    @test lower.value == lower.witness_output_norm
    @test lower.value <= exact_three_to_two * (1 + 16eps(Float64))
    @test lower.value ≈ exact_three_to_two rtol = 1e-10

    identity_result = QETInduced.induced_matrix_norm(
        MersenneTwister(11),
        Matrix{Float64}(I, 3, 3),
        3;
        initial_vector=[1.0, 0.0, 0.0],
        tolerance=0,
        max_iterations=3,
    )
    @test identity_result.status === :converged_lower_bound
    @test identity_result.bound_kind === :lower_bound
    @test !identity_result.exact
    @test identity_result.iterations == 1
    @test identity_result.iteration_residual == 0
    @test identity_result.value == 1

    rectangular = [
        1.0 2.0
        -3.0 0.5
        0.25 -1.0
    ]
    limited = QETInduced.induced_matrix_norm(
        MersenneTwister(12),
        rectangular,
        3;
        q=2,
        initial_vector=[1.0, -2.0],
        max_iterations=0,
    )
    @test limited.status === :iteration_limit
    @test limited.bound_kind === :lower_bound
    @test !limited.converged
    @test limited.iterations == 0
    @test limited.iteration_residual === nothing
    @test limited.value == norm(rectangular * limited.witness, 2)

    initialization_work =
        BigInt(size(rectangular, 1) * size(rectangular, 2)) + sum(size(rectangular))
    work_limited = QETInduced.induced_matrix_norm(
        MersenneTwister(13),
        rectangular,
        3;
        q=2,
        initial_vector=[1.0, -2.0],
        max_iterations=100,
        max_work=initialization_work,
    )
    @test work_limited.status === :work_limit
    @test work_limited.iterations == 0
    @test work_limited.work_used == initialization_work
    @test work_limited.max_work == initialization_work
    @test work_limited.value == norm(rectangular * work_limited.witness, 2)
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(13), rectangular, 3; q=2, max_work=initialization_work - 1
    )

    zero_result = QETInduced.induced_matrix_norm(
        MersenneTwister(14), zeros(Float32, 3, 2), Float32(3); q=Float32(4)
    )
    @test zero_result.status === :exact_zero
    @test zero_result.bound_kind === :exact
    @test zero_result.exact
    @test zero_result.value === 0.0f0
    @test zero_result.witness == Float32[1, 0]
    @test zero_result.iteration_residual === nothing

    complex_matrix = ComplexF64[1 im -1; 2-im -0.5im 3]
    complex_result = QETInduced.induced_matrix_norm(
        MersenneTwister(15), complex_matrix, Inf; q=1, max_iterations=30
    )
    @test complex_result.bound_kind === :lower_bound
    @test norm(complex_result.witness, Inf) ≈ 1
    @test complex_result.value ≈ norm(complex_matrix * complex_result.witness, 1)

    one_to_infinity = QETInduced.induced_matrix_norm(
        MersenneTwister(16),
        rectangular,
        1;
        q=Inf,
        initial_vector=[1.0, 0.0],
        max_iterations=10,
    )
    @test one_to_infinity.bound_kind === :lower_bound
    @test length(one_to_infinity.witness) == size(rectangular, 2)
    @test norm(one_to_infinity.witness, 1) == 1
    @test one_to_infinity.value <= maximum(abs, rectangular)
end

@testset "WP2 induced matrix norm RNG and sparse contracts" begin
    matrix = [
        1.0 2.0
        -3.0 0.5
        0.25 -1.0
    ]
    first_result = QETInduced.induced_matrix_norm(
        MersenneTwister(0x524e4741), matrix, 3; q=2, max_iterations=20
    )
    second_result = QETInduced.induced_matrix_norm(
        MersenneTwister(0x524e4741), matrix, 3; q=2, max_iterations=20
    )
    @test first_result.value == second_result.value
    @test first_result.witness == second_result.witness
    @test first_result.status === second_result.status
    @test first_result.iterations == second_result.iterations

    supplied_rng = MersenneTwister(0x524e4742)
    supplied_control = copy(supplied_rng)
    supplied_result = QETInduced.induced_matrix_norm(
        supplied_rng, matrix, 3; q=2, initial_vector=[1.0, -2.0], max_iterations=20
    )
    @test rand(supplied_rng, UInt64) == rand(supplied_control, UInt64)
    @test supplied_result.value > 0

    Random.seed!(0x474c4f42414c)
    expected_global = rand(UInt64, 4)
    Random.seed!(0x474c4f42414c)
    QETInduced.induced_matrix_norm(
        MersenneTwister(0x4c4f43414c), matrix, 3; q=2, max_iterations=5
    )
    @test rand(UInt64, 4) == expected_global

    sparse_matrix = sparse(matrix)
    sparse_one = QETInduced.induced_matrix_norm(MersenneTwister(20), sparse_matrix, 1)
    sparse_infinity = QETInduced.induced_matrix_norm(
        MersenneTwister(21), sparse_matrix, Inf
    )
    sparse_lower = QETInduced.induced_matrix_norm(
        MersenneTwister(22),
        sparse_matrix,
        3;
        q=2,
        initial_vector=[1.0, -2.0],
        max_iterations=20,
    )
    @test sparse_one.value == maximum(vec(sum(abs, matrix; dims=1)))
    @test sparse_infinity.value == maximum(vec(sum(abs, matrix; dims=2)))
    @test sparse_lower.value ≈ supplied_result.value
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(23), sparse_matrix, 2
    )
    sparse_two = QETInduced.induced_matrix_norm(
        MersenneTwister(24), sparse_matrix, 2; allow_densify=true
    )
    @test sparse_two.value ≈ opnorm(matrix, 2)
end

@testset "WP2 InducedMatrixNorm compatibility and validation" begin
    matrix = [
        1.0 2.0
        -3.0 0.5
        0.25 -1.0
    ]
    compatibility_exact = CompatInduced.InducedMatrixNorm(
        MersenneTwister(30), matrix, "FRO"
    )
    @test compatibility_exact.bound_kind === :exact
    @test compatibility_exact.value ≈ opnorm(matrix, 2)

    row_start = reshape([1.0, -2.0], 1, :)
    compatibility_lower = CompatInduced.InducedMatrixNorm(
        MersenneTwister(31), matrix, 3, 2, 1e-10, row_start; max_iterations=30
    )
    native_lower = QETInduced.induced_matrix_norm(
        MersenneTwister(99),
        matrix,
        3;
        q=2,
        tolerance=1e-10,
        initial_vector=vec(row_start),
        max_iterations=30,
    )
    @test compatibility_lower.value == native_lower.value
    @test compatibility_lower.witness == native_lower.witness
    @test compatibility_lower.bound_kind === :lower_bound

    sentinel = CompatInduced.InducedMatrixNorm(
        MersenneTwister(32), matrix, 3, 2, 1e-8, -1; max_iterations=10
    )
    omitted = CompatInduced.InducedMatrixNorm(
        MersenneTwister(32), matrix, 3, 2, 1e-8; max_iterations=10
    )
    @test sentinel.value == omitted.value
    @test sentinel.witness == omitted.witness
    @test_throws MethodError CompatInduced.InducedMatrixNorm(matrix, 1)
    @test_throws DimensionMismatch CompatInduced.InducedMatrixNorm(
        MersenneTwister(33), matrix, 3, 2, 1e-8, ones(2, 2)
    )
    @test_throws ArgumentError CompatInduced.InducedMatrixNorm(
        MersenneTwister(34), matrix, "nuclear"
    )

    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(40), zeros(Float64, 0, 2), 2
    )
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(41), [1 2; 3 4], 2
    )
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(42), Matrix{BigFloat}(I, 2, 2), 2
    )
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(43), [1.0 NaN; 0 1], 2
    )
    for invalid_order in (0, -1, NaN, -Inf, true)
        @test_throws ArgumentError QETInduced.induced_matrix_norm(
            MersenneTwister(44), matrix, invalid_order
        )
        @test_throws ArgumentError QETInduced.induced_matrix_norm(
            MersenneTwister(45), matrix, 2; q=invalid_order
        )
    end
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(46), Float32.(matrix), 1.0e100
    )
    for invalid_tolerance in (-1, Inf, NaN, true)
        @test_throws ArgumentError QETInduced.induced_matrix_norm(
            MersenneTwister(47), matrix, 3; q=2, tolerance=invalid_tolerance
        )
    end
    for invalid_iterations in (-1, true, 1.5)
        @test_throws ArgumentError QETInduced.induced_matrix_norm(
            MersenneTwister(48), matrix, 3; q=2, max_iterations=invalid_iterations
        )
    end
    for invalid_work in (0, -1, true, 1.5)
        @test_throws ArgumentError QETInduced.induced_matrix_norm(
            MersenneTwister(49), matrix, 3; q=2, max_work=invalid_work
        )
    end
    @test_throws DimensionMismatch QETInduced.induced_matrix_norm(
        MersenneTwister(50), matrix, 3; q=2, initial_vector=ones(3)
    )
    @test_throws DomainError QETInduced.induced_matrix_norm(
        MersenneTwister(51), matrix, 3; q=2, initial_vector=zeros(2)
    )
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(52), matrix, 3; q=2, initial_vector=[1.0, NaN]
    )
    @test_throws ArgumentError QETInduced.induced_matrix_norm(
        MersenneTwister(53), matrix, 3; q=2, initial_vector=Float32[1, 2]
    )

    original_matrix = copy(matrix)
    original_start = ComplexF64[1 + im, -2im]
    start_copy = copy(original_start)
    complex_start = QETInduced.induced_matrix_norm(
        MersenneTwister(54), matrix, 3; q=2, initial_vector=original_start, max_iterations=3
    )
    @test eltype(complex_start.witness) == ComplexF64
    @test matrix == original_matrix
    @test original_start == start_copy
end
