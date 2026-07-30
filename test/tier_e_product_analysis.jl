using LinearAlgebra
using Random
using SparseArrays

const QETProduct = QuantumEntanglementTools

struct _TierEProductZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_TierEProductZeroBasedVector) = size(vector.storage)
Base.axes(vector::_TierEProductZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_TierEProductZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_TierEProductZeroBasedVector, index::Int) = vector.storage[index + 1]

struct _TierEProductZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_TierEProductZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_TierEProductZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_TierEProductZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierEProductZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

@testset "Tier E product-analysis array axes validation" begin
    vector = _TierEProductZeroBasedVector([1.0, 0.0, 0.0, 0.0])
    operator = _TierEProductZeroBasedMatrix(Matrix{Float64}(I, 4, 4))
    @test_throws ArgumentError QETProduct.is_product_vector(vector, (2, 2))
    @test_throws ArgumentError QETProduct.operator_schmidt_decomposition(operator, (2, 2))
    @test_throws ArgumentError QETProduct.is_product_operator(operator, (2, 2))
    @test_throws ArgumentError QETProduct.entanglement_of_formation(vector, (2, 2))
end

@testset "Tier E operator Schmidt analysis" begin
    left = ComplexF64[
        1 2im -1
        3 + im 0 2
    ]
    right = ComplexF64[
        0 1 - im
        2 3im
    ]
    rectangular_product = QETProduct.tensor_product(left, right)
    decomposition = QETProduct.operator_schmidt_decomposition(
        rectangular_product, (2, 2), (3, 2)
    )
    @test decomposition isa QETProduct.OperatorSchmidtDecompositionResult
    @test decomposition.row_dims == (2, 2)
    @test decomposition.column_dims == (3, 2)
    @test length(decomposition.coefficients) == 4
    @test all(size(factor) == (2, 3) for factor in decomposition.left_factors)
    @test all(size(factor) == (2, 2) for factor in decomposition.right_factors)
    @test QETProduct.tensor_sum(
        decomposition.left_factors,
        decomposition.right_factors;
        weights=decomposition.coefficients,
    ) ≈ rectangular_product
    @test QETProduct.operator_schmidt_coefficients(rectangular_product, (2, 2), (3, 2)) ≈
        decomposition.coefficients
    @test QETProduct.operator_schmidt_rank(rectangular_product, (2, 2), (3, 2)) == 1
    @test norm(decomposition.left_factors[1]) ≈ 1
    @test norm(decomposition.right_factors[1]) ≈ 1

    first_unit = [1.0 0.0; 0.0 0.0]
    second_unit = [0.0 0.0; 0.0 1.0]
    rank_two =
        QETProduct.tensor_product(first_unit, first_unit) +
        2QETProduct.tensor_product(second_unit, second_unit)
    @test QETProduct.operator_schmidt_coefficients(rank_two, (2, 2)) ≈ [2.0, 1.0, 0.0, 0.0]
    @test QETProduct.operator_schmidt_rank(rank_two, (2, 2)) == 2

    near_product =
        QETProduct.tensor_product(first_unit, first_unit) +
        1e-7QETProduct.tensor_product(second_unit, second_unit)
    @test QETProduct.operator_schmidt_rank(near_product, (2, 2); atol=1e-6, rtol=0) == 1
    @test QETProduct.operator_schmidt_rank(near_product, (2, 2); atol=1e-8, rtol=0) == 2

    sparse_product = sparse(rectangular_product)
    @test_throws ArgumentError QETProduct.operator_schmidt_decomposition(
        sparse_product, (2, 2), (3, 2)
    )
    @test QETProduct.operator_schmidt_rank(
        sparse_product, (2, 2), (3, 2); allow_densify=true
    ) == 1

    @test_throws ArgumentError QETProduct.operator_schmidt_decomposition(
        rectangular_product, (4,), (6,)
    )
    @test_throws DimensionMismatch QETProduct.operator_schmidt_decomposition(
        rectangular_product, (2, 2), (3, 2, 1)
    )
    @test_throws DimensionMismatch QETProduct.operator_schmidt_decomposition(
        rectangular_product, (2, 3), (3, 2)
    )
    @test_throws DimensionMismatch QETProduct.operator_schmidt_decomposition(
        rectangular_product, (2, 2)
    )
    @test_throws ArgumentError QETProduct.operator_schmidt_rank(
        rectangular_product, (2, 2), (3, 2); atol=-1
    )
    @test_throws ArgumentError QETProduct.operator_schmidt_rank(
        rectangular_product, (2, 2), (3, 2); rtol=Inf
    )
    invalid = copy(rectangular_product)
    invalid[1, 1] = NaN
    @test_throws ArgumentError QETProduct.operator_schmidt_coefficients(
        invalid, (2, 2), (3, 2)
    )
    @test_throws ArgumentError QETProduct.operator_schmidt_coefficients(
        Complex{BigFloat}.(rectangular_product), (2, 2), (3, 2)
    )
    @test occursin("OperatorSchmidtDecompositionResult", sprint(show, decomposition))
end

@testset "Tier E product-vector analysis" begin
    first = ComplexF64[1, 2im]
    second = ComplexF64[3, -1]
    third = ComplexF64[2, im]
    vector = QETProduct.tensor_product(first, second, third)
    result = QETProduct.is_product_vector(vector, (2, 2, 2))
    @test result isa QETProduct.ProductAnalysisResult
    @test result.status === :within_tolerance
    @test length(result.factors) == 3
    @test map(length, result.factors) == (2, 2, 2)
    @test QETProduct.tensor_product(result.factors...) ≈ vector
    @test result.reconstruction_residual ≤ 1e-12
    @test result.approximation_threshold > 0
    @test length(result.cut_residuals) == 2
    @test length(result.cut_thresholds) == 2
    @test all(result.cut_residuals .<= result.cut_thresholds)

    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    bell_result = QETProduct.is_product_vector(bell, (2, 2); atol=0, rtol=0)
    @test bell_result.status === :outside_tolerance
    @test bell_result.cut_residuals[1] ≈ inv(sqrt(2))
    @test bell_result.reconstruction_residual ≈ inv(sqrt(2))

    near_product = Float64[1, 0, 0, 1e-7]
    @test QETProduct.is_product_vector(near_product, (2, 2); atol=1e-6, rtol=0).status ===
        :within_tolerance
    @test QETProduct.is_product_vector(near_product, (2, 2); atol=1e-8, rtol=0).status ===
        :outside_tolerance
    boundary = QETProduct.is_product_vector(near_product, (2, 2); atol=1e-7, rtol=0)
    @test boundary.status === :boundary
    @test boundary.cut_residuals == boundary.cut_thresholds

    scaled_entangled_remainder = 1e-12QETProduct.tensor_product(ComplexF64[1, 0], bell)
    scaled_result = QETProduct.is_product_vector(
        scaled_entangled_remainder, (2, 2, 2); atol=1e-10, rtol=0
    )
    @test scaled_result.status === :within_tolerance
    @test scaled_result.reconstruction_residual ≤ scaled_result.approximation_threshold
    @test scaled_result.cut_residuals[2] < 1e-12

    tail_spectrum = vcat(1.0, fill(0.009, 15))
    many_small_tails = vec(Matrix(Diagonal(tail_spectrum)))
    tail_result = QETProduct.is_product_vector(
        many_small_tails, (16, 16); atol=0.01, rtol=0
    )
    @test tail_result.status === :outside_tolerance
    @test tail_result.cut_residuals[1] ≈ norm(tail_spectrum[2:end])
    @test tail_result.cut_residuals[1] > tail_result.cut_thresholds[1]

    sparse_vector = sparsevec(findall(!iszero, vector), vector[findall(!iszero, vector)], 8)
    @test_throws ArgumentError QETProduct.is_product_vector(sparse_vector, (2, 2, 2))
    @test QETProduct.is_product_vector(sparse_vector, (2, 2, 2); allow_densify=true).status ===
        :within_tolerance

    @test_throws DomainError QETProduct.is_product_vector(zeros(ComplexF64, 4), (2, 2))
    @test_throws ArgumentError QETProduct.is_product_vector(ComplexF64[1, 0, 0, 0], (4,))
    @test_throws DimensionMismatch QETProduct.is_product_vector(
        ComplexF64[1, 0, 0, 0], (2, 3)
    )
    @test_throws ArgumentError QETProduct.is_product_vector(
        ComplexF64[1, NaN, 0, 0], (2, 2)
    )
    @test_throws ArgumentError QETProduct.is_product_vector(BigFloat[1, 0, 0, 0], (2, 2))
    @test_throws ArgumentError QETProduct.is_product_vector(
        ComplexF64[1, 0, 0, 0], (2, 2); atol=true
    )
    @test_throws ArgumentError QETProduct.is_product_vector(
        ComplexF64[1, 0, 0, 0], (2, 2); rtol=-1
    )
    @test occursin("ProductAnalysisResult", sprint(show, result))
end

@testset "Tier E product-operator analysis" begin
    first = ComplexF64[
        1 2
        3im 4
    ]
    second = ComplexF64[
        0 1 2
        3 4im 5
    ]
    third = reshape(ComplexF64[1, 2], 2, 1)
    operator = QETProduct.tensor_product(first, second, third)
    result = QETProduct.is_product_operator(operator, (2, 2, 2), (2, 3, 1))
    @test result.status === :within_tolerance
    @test map(size, result.factors) == ((2, 2), (2, 3), (2, 1))
    @test QETProduct.tensor_product(result.factors...) ≈ operator
    @test result.reconstruction_residual ≤ 1e-11
    @test all(result.cut_residuals .<= result.cut_thresholds)

    first_unit = [1.0 0.0; 0.0 0.0]
    second_unit = [0.0 0.0; 0.0 1.0]
    nonproduct =
        QETProduct.tensor_product(first_unit, first_unit) +
        QETProduct.tensor_product(second_unit, second_unit)
    @test QETProduct.is_product_operator(nonproduct, (2, 2); atol=0, rtol=0).status ===
        :outside_tolerance

    near_product =
        QETProduct.tensor_product(first_unit, first_unit) +
        1e-7QETProduct.tensor_product(second_unit, second_unit)
    @test QETProduct.is_product_operator(near_product, (2, 2); atol=1e-6, rtol=0).status ===
        :within_tolerance
    @test QETProduct.is_product_operator(near_product, (2, 2); atol=1e-8, rtol=0).status ===
        :outside_tolerance
    @test QETProduct.is_product_operator(near_product, (2, 2); atol=1e-7, rtol=0).status ===
        :boundary

    scaled_nonproduct = 1e-12QETProduct.tensor_product(first_unit, nonproduct)
    scaled_result = QETProduct.is_product_operator(
        scaled_nonproduct, (2, 2, 2); atol=1e-10, rtol=0
    )
    @test scaled_result.status === :within_tolerance
    @test scaled_result.reconstruction_residual ≤ scaled_result.approximation_threshold
    @test scaled_result.cut_residuals[2] < 1e-12

    sparse_operator = sparse(operator)
    @test_throws ArgumentError QETProduct.is_product_operator(
        sparse_operator, (2, 2, 2), (2, 3, 1)
    )
    @test QETProduct.is_product_operator(
        sparse_operator, (2, 2, 2), (2, 3, 1); allow_densify=true
    ).status === :within_tolerance

    @test_throws DomainError QETProduct.is_product_operator(zeros(ComplexF64, 4, 4), (2, 2))
    @test_throws ArgumentError QETProduct.is_product_operator(
        Matrix{Float64}(I, 4, 4), (4,)
    )
    @test_throws DimensionMismatch QETProduct.is_product_operator(
        operator, (2, 2, 2), (2, 3)
    )
    @test_throws DimensionMismatch QETProduct.is_product_operator(
        operator, (2, 2, 3), (2, 3, 1)
    )
    invalid = copy(operator)
    invalid[1, 1] = Inf
    @test_throws ArgumentError QETProduct.is_product_operator(invalid, (2, 2, 2), (2, 3, 1))
    @test_throws ArgumentError QETProduct.is_product_operator(
        Complex{BigFloat}.(operator), (2, 2, 2), (2, 3, 1)
    )
    @test_throws ArgumentError QETProduct.is_product_operator(
        operator, (2, 2, 2), (2, 3, 1); rtol=NaN
    )
end

@testset "Tier E entanglement of formation" begin
    product = ComplexF64[1, 0, 0, 0]
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    @test QETProduct.entanglement_of_formation(product, (2, 2)) == 0
    @test QETProduct.entanglement_of_formation(bell, (2, 2)) ≈ 1
    @test QETProduct.entanglement_of_formation(bell, (2, 2); base=exp(1)) ≈ log(2)
    bell32 = ComplexF32[1, 0, 0, 1] / sqrt(Float32(2))
    bell32_eof = QETProduct.entanglement_of_formation(bell32, (2, 2))
    @test bell32_eof isa Float32
    @test bell32_eof ≈ 1.0f0 atol = 8eps(Float32)
    @test 0.0f0 ≤ bell32_eof ≤ 1.0f0 + sqrt(eps(Float32))

    angle = 0.31
    rectangular_pure = ComplexF64[cos(angle), 0, 0, 0, sin(angle), 0]
    expected = -cos(angle)^2 * log2(cos(angle)^2) - sin(angle)^2 * log2(sin(angle)^2)
    @test QETProduct.entanglement_of_formation(rectangular_pure, (2, 3)) ≈ expected
    exact_rectangular_density = zeros(6, 6)
    exact_rectangular_density[1, 1] = 0.5
    exact_rectangular_density[1, 5] = 0.5
    exact_rectangular_density[5, 1] = 0.5
    exact_rectangular_density[5, 5] = 0.5
    @test QETProduct.entanglement_of_formation(exact_rectangular_density, (2, 3)) ≈ 1
    @test QETProduct.entanglement_of_formation(
        exact_rectangular_density, (2, 3); base=exp(1)
    ) ≈ log(2)
    exact_product_density = Diagonal([1.0, 0, 0, 0, 0, 0])
    @test QETProduct.entanglement_of_formation(exact_product_density, (2, 3)) == 0

    numerical_rectangular_density = rectangular_pure * adjoint(rectangular_pure)
    @test_throws DomainError QETProduct.entanglement_of_formation(
        numerical_rectangular_density, (2, 3)
    )
    @test QETProduct.entanglement_of_formation(
        numerical_rectangular_density,
        (2, 3);
        psd_boundary_policy=:project,
        rank_boundary_policy=:project,
    ) ≈ expected

    bell_density = bell * adjoint(bell)
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    @test QETProduct.entanglement_of_formation(bell_density, (2, 2)) ≈ 1
    @test QETProduct.entanglement_of_formation(maximally_mixed, (2, 2)) == 0
    @test QETProduct.entanglement_of_formation(Diagonal([1.0, 0, 0, 0]), (2, 2)) == 0
    @test QETProduct.entanglement_of_formation(Diagonal([0.5, 0, 0, 0.5]), (2, 2)) == 0
    bell_mixture = 0.7bell_density + 0.3maximally_mixed
    mixture_concurrence = 0.55
    mixture_probability = (1 + sqrt(1 - mixture_concurrence^2)) / 2
    expected_mixture =
        -mixture_probability * log2(mixture_probability) -
        (1 - mixture_probability) * log2(1 - mixture_probability)
    @test QETProduct.entanglement_of_formation(bell_mixture, (2, 2)) ≈ expected_mixture

    rng = MersenneTwister(0x514554)
    for _ in 1:12
        random_pure = randn(rng, ComplexF64, 4)
        random_pure ./= norm(random_pure)
        random_density = random_pure * adjoint(random_pure)
        matrix_eof = QETProduct.entanglement_of_formation(
            random_density,
            (2, 2);
            psd_boundary_policy=:project,
            range_boundary_policy=:project,
        )
        @test isfinite(matrix_eof)
        @test matrix_eof ≈ QETProduct.entanglement_of_formation(random_pure, (2, 2)) atol =
            2e-7
    end
    first_rank_two = randn(rng, ComplexF64, 4)
    first_rank_two ./= norm(first_rank_two)
    second_rank_two = randn(rng, ComplexF64, 4)
    second_rank_two ./= norm(second_rank_two)
    rank_two_density =
        0.37(first_rank_two * adjoint(first_rank_two)) +
        0.63(second_rank_two * adjoint(second_rank_two))
    rank_two_eof = QETProduct.entanglement_of_formation(
        rank_two_density, (2, 2); psd_boundary_policy=:project
    )
    @test isfinite(rank_two_eof)
    @test 0 ≤ rank_two_eof ≤ 1 + 1e-12

    for _ in 1:8
        local_a = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        local_b = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
        rotated_bell = kron(local_a, local_b) * bell
        rotated_density = rotated_bell * adjoint(rotated_bell)
        @test QETProduct.entanglement_of_formation(
            rotated_density,
            (2, 2);
            psd_boundary_policy=:project,
            range_boundary_policy=:project,
        ) ≈ 1 atol = 2e-7
    end

    sparse_bell = sparsevec([1, 4], fill(inv(sqrt(2)), 2), 4)
    @test_throws ArgumentError QETProduct.entanglement_of_formation(sparse_bell, (2, 2))
    @test QETProduct.entanglement_of_formation(sparse_bell, (2, 2); allow_densify=true) ≈ 1
    sparse_mixed = sparse(maximally_mixed)
    @test_throws ArgumentError QETProduct.entanglement_of_formation(sparse_mixed, (2, 2))
    @test QETProduct.entanglement_of_formation(sparse_mixed, (2, 2); allow_densify=true) ==
        0
    sparse_rectangular_density = sparse(exact_rectangular_density)
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        sparse_rectangular_density, (2, 3)
    )
    @test QETProduct.entanglement_of_formation(
        sparse_rectangular_density, (2, 3); allow_densify=true
    ) ≈ 1

    rank_boundary = Diagonal([1 - 1e-10, 1e-10, 0.0, 0.0, 0.0, 0.0])
    @test_throws DomainError QETProduct.entanglement_of_formation(
        rank_boundary, (2, 3); atol=1e-9, rtol=0
    )
    @test QETProduct.entanglement_of_formation(
        rank_boundary, (2, 3); atol=1e-9, rtol=0, rank_boundary_policy=:project
    ) == 0
    psd_boundary = Diagonal([1 + 1e-10, -1e-10, 0.0, 0.0, 0.0, 0.0])
    @test_throws DomainError QETProduct.entanglement_of_formation(
        psd_boundary, (2, 3); atol=1e-9, rtol=0
    )
    @test QETProduct.entanglement_of_formation(
        psd_boundary,
        (2, 3);
        atol=1e-9,
        rtol=0,
        psd_boundary_policy=:project,
        rank_boundary_policy=:project,
    ) == 0
    high_dimensional_mixed = Diagonal([0.5, 0.5, 0.0, 0.0, 0.0, 0.0])
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        high_dimensional_mixed, (2, 3)
    )
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        high_dimensional_mixed,
        (2, 3);
        psd_boundary_policy=:project,
        rank_boundary_policy=:project,
    )

    two_qubit_psd_boundary = Diagonal([-1e-10, 0.25, 0.25, 0.5000000001])
    @test_throws DomainError QETProduct.entanglement_of_formation(
        two_qubit_psd_boundary, (2, 2); atol=1e-9, rtol=0
    )
    @test QETProduct.entanglement_of_formation(
        two_qubit_psd_boundary, (2, 2); atol=1e-9, rtol=0, psd_boundary_policy=:project
    ) == 0

    @test_throws ArgumentError QETProduct.entanglement_of_formation(0.9bell, (2, 2))
    @test_throws ArgumentError QETProduct.entanglement_of_formation(bell, (2, 2); base=1)
    @test_throws ArgumentError QETProduct.entanglement_of_formation(bell, (2, 2); base=true)
    @test_throws ArgumentError QETProduct.entanglement_of_formation(bell, (4,))
    @test_throws DimensionMismatch QETProduct.entanglement_of_formation(bell, (2, 3))
    @test_throws DimensionMismatch QETProduct.entanglement_of_formation(ones(4, 3), (2, 2))
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        2maximally_mixed, (2, 2)
    )
    @test_throws DomainError QETProduct.entanglement_of_formation(
        Diagonal([-0.1, 0.4, 0.3, 0.4]), (2, 2)
    )
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        Complex{BigFloat}.(bell), (2, 2)
    )
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        bell_density, (2, 2); psd_boundary_policy=:unsupported
    )
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        bell_density, (2, 2); rank_boundary_policy=:unsupported
    )
    @test_throws ArgumentError QETProduct.entanglement_of_formation(
        bell_density, (2, 2); range_boundary_policy=:unsupported
    )
end

@testset "Tier E sufficient separable ball" begin
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    mixed_result = QETProduct.in_separable_ball(maximally_mixed, (2, 2))
    @test mixed_result isa QETProduct.SeparableBallResult
    @test mixed_result.status === :separable_certified
    @test mixed_result.purity == 0.25
    @test mixed_result.boundary ≈ 1 / 3
    @test occursin("certifies", mixed_result.message)
    @test QETProduct.in_separable_ball(Diagonal(fill(0.25, 4)), (2, 2)).status ===
        :separable_certified
    @test QETProduct.in_separable_ball(Diagonal(fill(BigFloat(1) / 4, 4)), (2, 2)).status ===
        :separable_certified

    near_hermitian = ComplexF64.(maximally_mixed)
    near_hermitian[1, 2] = 1e-9
    @test QETProduct.in_separable_ball(near_hermitian, (2, 2)).status === :unknown
    small_imaginary_trace = ComplexF64.(maximally_mixed)
    small_imaginary_trace[1, 1] += 1e-9im
    @test QETProduct.in_separable_ball(small_imaginary_trace, (2, 2)).status === :unknown
    @test QETProduct.in_separable_ball(
        Diagonal(ComplexF64[0.25 + 1e-9im, 0.25, 0.25, 0.25]), (2, 2)
    ).status === :unknown
    @test QETProduct.in_separable_ball([0.250000001, 0.25, 0.25, 0.25], (2, 2)).status ===
        :unknown

    pure_result = QETProduct.in_separable_ball([1.0, 0, 0, 0], (2, 2))
    @test pure_result.status === :outside_ball
    @test !occursin("entangled", lowercase(pure_result.message))
    @test occursin("no entanglement conclusion", lowercase(pure_result.message))

    boundary_values = [1 / 3, 1 / 3, 1 / 3, 0.0]
    boundary_result = QETProduct.in_separable_ball(boundary_values, (2, 2))
    @test boundary_result.status === :unknown
    @test QETProduct.in_separable_ball(boundary_values, (2, 2); atol=0, rtol=0).status ===
        :separable_certified

    sparse_mixed = sparse(maximally_mixed)
    @test_throws ArgumentError QETProduct.in_separable_ball(sparse_mixed, (2, 2))
    @test QETProduct.in_separable_ball(sparse_mixed, (2, 2); allow_densify=true).status ===
        :separable_certified
    sparse_values = sparsevec([1, 2, 3], fill(1 / 3, 3), 4)
    @test QETProduct.in_separable_ball(sparse_values, (2, 2)).status === :unknown
    @test QETProduct.in_separable_ball(fill(BigFloat(1) / 4, 4), (2, 2)).status ===
        :separable_certified

    boundary_invalid = Diagonal([-1e-10, 0.5, 0.3, 0.2000000001])
    @test QETProduct.in_separable_ball(boundary_invalid, (2, 2); atol=1e-9, rtol=0).status ===
        :unknown
    @test QETProduct.in_separable_ball(
        diag(boundary_invalid), (2, 2); atol=1e-9, rtol=0
    ).status === :unknown

    @test_throws ArgumentError QETProduct.in_separable_ball(2maximally_mixed, (2, 2))
    @test_throws ArgumentError QETProduct.in_separable_ball(fill(0.5, 4), (2, 2))
    @test_throws DomainError QETProduct.in_separable_ball([-0.1, 0.4, 0.3, 0.4], (2, 2))
    @test_throws ArgumentError QETProduct.in_separable_ball(
        ComplexF64[0.25, 0.25, 0.25, 0.25], (2, 2)
    )
    @test_throws ArgumentError QETProduct.in_separable_ball(
        Rational{Int}[1 // 4, 1 // 4, 1 // 4, 1 // 4], (2, 2)
    )
    @test_throws ArgumentError QETProduct.in_separable_ball(
        [0.25, 0.25, 0.25, 0.25], (2, 2); allow_densify=true
    )
    @test_throws ArgumentError QETProduct.in_separable_ball([0.25, 0.25, 0.25, Inf], (2, 2))
    @test_throws ArgumentError QETProduct.in_separable_ball(maximally_mixed, (4,))
    @test_throws DimensionMismatch QETProduct.in_separable_ball(maximally_mixed, (2, 3))
    @test_throws ArgumentError QETProduct.in_separable_ball(reshape([1.0], 1, 1), (1, 1))
    @test_throws ArgumentError QETProduct.in_separable_ball([1.0], (1, 1))
    @test_throws ArgumentError QETProduct.in_separable_ball(
        [0.25, 0.25, 0.25, 0.25], (2, 2); atol=-1
    )
    @test_throws ArgumentError QETProduct.in_separable_ball(
        [0.25, 0.25, 0.25, 0.25], (2, 2); rtol=true
    )
    @test occursin("SeparableBallResult", sprint(show, mixed_result))
end
