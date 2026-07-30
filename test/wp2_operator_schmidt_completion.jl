using LinearAlgebra
using Random
using SparseArrays

const QETOperatorSchmidt = QuantumEntanglementTools
const CompatOperatorSchmidt = QuantumEntanglementTools.MATLABCompat

function _wp2_operator_schmidt_reconstruction(result)
    return tensor_sum(
        result.left_factors, result.right_factors; weights=result.coefficients
    )
end

function _wp2_operator_schmidt_gram(factors)
    vectors = hcat(vec.(factors)...)
    return adjoint(vectors) * vectors
end

function _wp2_operator_schmidt_projector(factors, count)
    vectors = hcat(vec.(factors[1:count])...)
    return vectors * adjoint(vectors)
end

@testset "WP2 OperatorSchmidtDecomposition completion" begin
    @testset "default, scalar, and positional DIM contracts" begin
        left = [1.0 0.0; 0.0 2.0]
        right = [0.0 1.0; 1.0 0.0]
        operator = tensor_product(left, right)

        inferred = CompatOperatorSchmidt.OperatorSchmidtDecomposition(operator)
        scalar = CompatOperatorSchmidt.OperatorSchmidtDecomposition(operator, 2)
        singleton = CompatOperatorSchmidt.OperatorSchmidtDecomposition(operator, (2,))
        vector = CompatOperatorSchmidt.OperatorSchmidtDecomposition(operator, [2, 2])
        @test keys(inferred) == (:coefficients, :left_factors, :right_factors)
        @test length(inferred) == 3
        @test inferred.coefficients ≈ scalar.coefficients
        @test inferred.coefficients ≈ singleton.coefficients
        @test inferred.coefficients ≈ vector.coefficients
        @test _wp2_operator_schmidt_reconstruction(inferred) ≈ operator
        @test all(ishermitian, inferred.left_factors)
        @test all(ishermitian, inferred.right_factors)

        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            Matrix{Float64}(I, 6, 6)
        )
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            operator, 3
        )
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            ones(4, 6), 2
        )
    end

    @testset "Hermitian convention with unequal local dimensions" begin
        rng = MersenneTwister(0x4f50455241544f52)
        raw = randn(rng, ComplexF64, 6, 6)
        operator = raw + adjoint(raw)

        general = QETOperatorSchmidt.operator_schmidt_decomposition(operator, (2, 3))
        hermitian = QETOperatorSchmidt.operator_schmidt_decomposition(
            operator, (2, 3); hermitian_factors=true
        )
        @test general.factor_convention === :general
        @test general.coordinate_imaginary_residual === nothing
        @test general.coordinate_imaginary_tolerance === nothing
        @test hermitian.factor_convention === :hermitian
        @test hermitian.row_dims == (2, 3)
        @test hermitian.column_dims == (2, 3)
        @test hermitian.coordinate_imaginary_residual <=
            hermitian.coordinate_imaginary_tolerance
        @test hermitian.coordinate_imaginary_tolerance > 0
        @test length(hermitian.coefficients) == 4
        @test hermitian.coefficients ≈ general.coefficients atol = 2e-12
        @test _wp2_operator_schmidt_reconstruction(hermitian) ≈ operator atol = 2e-12
        @test all(ishermitian, hermitian.left_factors)
        @test all(ishermitian, hermitian.right_factors)
        @test _wp2_operator_schmidt_gram(hermitian.left_factors) ≈
            Matrix{ComplexF64}(I, 4, 4) atol = 2e-12
        @test _wp2_operator_schmidt_gram(hermitian.right_factors) ≈
            Matrix{ComplexF64}(I, 4, 4) atol = 2e-12

        compatibility = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            operator, (2, 3), -1
        )
        coefficients, left_factors, right_factors = compatibility
        @test coefficients ≈ hermitian.coefficients atol = 2e-12
        @test length(left_factors) == length(coefficients)
        @test length(right_factors) == length(coefficients)
        @test all(ishermitian, left_factors)
        @test all(ishermitian, right_factors)
        @test _wp2_operator_schmidt_reconstruction(compatibility) ≈ operator atol = 2e-12

        compatibility_general = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            operator, (2, 3), -1; hermitian_factors=false
        )
        @test _wp2_operator_schmidt_reconstruction(compatibility_general) ≈ operator atol =
            2e-12
        @test occursin("factor_convention=hermitian", sprint(show, hermitian))
    end

    @testset "degenerate gauges and K selection" begin
        left_identity = Matrix{Float64}(I, 2, 2) / sqrt(2)
        left_z = [1.0 0.0; 0.0 -1.0] / sqrt(2)
        right_one = Diagonal([1.0, -1.0, 0.0]) / sqrt(2)
        right_two = Diagonal([1.0, 1.0, -2.0]) / sqrt(6)
        rank_two =
            3tensor_product(left_identity, right_one) + 2tensor_product(left_z, right_two)

        nonzero = CompatOperatorSchmidt.OperatorSchmidtDecomposition(rank_two, (2, 3))
        full = CompatOperatorSchmidt.OperatorSchmidtDecomposition(rank_two, (2, 3), -1)
        leading = CompatOperatorSchmidt.OperatorSchmidtDecomposition(rank_two, (2, 3), 1)
        @test nonzero.coefficients ≈ [3.0, 2.0]
        @test length(full.coefficients) == 4
        @test full.coefficients[1:2] ≈ [3.0, 2.0]
        @test full.coefficients[3:4] ≈ zeros(2) atol = 1e-14
        @test leading.coefficients ≈ [3.0]
        @test _wp2_operator_schmidt_reconstruction(nonzero) ≈ rank_two atol = 1e-12
        @test _wp2_operator_schmidt_reconstruction(full) ≈ rank_two atol = 1e-12
        @test all(ishermitian, full.left_factors)
        @test all(ishermitian, full.right_factors)

        threshold_operator =
            tensor_product(left_identity, right_one) +
            eps(Float64) * tensor_product(left_z, right_two)
        threshold_nonzero = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            threshold_operator, (2, 3), 0
        )
        threshold_full = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            threshold_operator, (2, 3), -1
        )
        @test length(threshold_nonzero.coefficients) == 1
        @test length(threshold_full.coefficients) == 4
        @test threshold_full.coefficients[2] <= 9eps(first(threshold_full.coefficients))

        degenerate =
            2tensor_product(left_identity, right_one) + 2tensor_product(left_z, right_two)
        degenerate_result = QETOperatorSchmidt.operator_schmidt_decomposition(
            degenerate, (2, 3); hermitian_factors=true
        )
        expected_left = hcat(vec(left_identity), vec(left_z))
        expected_right = hcat(vec(right_one), vec(right_two))
        @test degenerate_result.coefficients[1:2] ≈ [2.0, 2.0]
        @test _wp2_operator_schmidt_projector(degenerate_result.left_factors, 2) ≈
            expected_left * adjoint(expected_left) atol = 1e-12
        @test _wp2_operator_schmidt_projector(degenerate_result.right_factors, 2) ≈
            expected_right * adjoint(expected_right) atol = 1e-12
        @test _wp2_operator_schmidt_reconstruction(degenerate_result) ≈ degenerate atol =
            1e-12

        zero_nonzero = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            zeros(4, 4), (2, 2), 0
        )
        zero_full = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            zeros(4, 4), (2, 2), -1
        )
        @test isempty(zero_nonzero.coefficients)
        @test isempty(zero_nonzero.left_factors)
        @test isempty(zero_nonzero.right_factors)
        @test length(zero_full.coefficients) == 4
        @test all(iszero, zero_full.coefficients)
        @test _wp2_operator_schmidt_reconstruction(zero_full) == zeros(4, 4)

        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            rank_two, (2, 3), -2
        )
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            rank_two, (2, 3), 5
        )
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            rank_two, (2, 3), true
        )
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            rank_two, (2, 3), big(typemax(Int)) + 1
        )
    end

    @testset "rectangular layouts and exact Hermiticity boundary" begin
        left = Float32[1 2 0; 0 1 3]
        right = Float32[0 1; 2 3; 1 -1]
        locally_rectangular = tensor_product(left, right)
        general = QETOperatorSchmidt.operator_schmidt_decomposition(
            locally_rectangular, (2, 3), (3, 2)
        )
        @test general.factor_convention === :general
        @test eltype(general.coefficients) == Float32
        @test all(factor -> eltype(factor) == Float32, general.left_factors)
        @test all(factor -> eltype(factor) == Float32, general.right_factors)
        @test all(size(factor) == (2, 3) for factor in general.left_factors)
        @test all(size(factor) == (3, 2) for factor in general.right_factors)
        @test _wp2_operator_schmidt_reconstruction(general) ≈ locally_rectangular atol =
            2.0f-5
        @test _wp2_operator_schmidt_gram(general.left_factors) ≈ Matrix{Float32}(I, 6, 6) atol =
            2.0f-5
        @test _wp2_operator_schmidt_gram(general.right_factors) ≈ Matrix{Float32}(I, 6, 6) atol =
            2.0f-5
        @test_throws ArgumentError QETOperatorSchmidt.operator_schmidt_decomposition(
            locally_rectangular, (2, 3), (3, 2); hermitian_factors=true
        )

        globally_hermitian = Matrix{Float64}(I, 6, 6)
        cross_layout = QETOperatorSchmidt.operator_schmidt_decomposition(
            globally_hermitian, (2, 3), (3, 2)
        )
        @test _wp2_operator_schmidt_reconstruction(cross_layout) ≈ globally_hermitian
        cross_compatibility = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            globally_hermitian, [2 3; 3 2], -1
        )
        @test all(size(factor) == (2, 3) for factor in cross_compatibility.left_factors)
        @test all(size(factor) == (3, 2) for factor in cross_compatibility.right_factors)
        @test _wp2_operator_schmidt_reconstruction(cross_compatibility) ≈ globally_hermitian
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            globally_hermitian, [2 3; 3 2], -1; hermitian_factors=true
        )

        almost_hermitian = copy(globally_hermitian)
        almost_hermitian[1, 2] = 1e-12
        automatic = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            almost_hermitian, (2, 3), -1
        )
        @test _wp2_operator_schmidt_reconstruction(automatic) ≈ almost_hermitian atol =
            1e-14
        @test_throws ArgumentError QETOperatorSchmidt.operator_schmidt_decomposition(
            almost_hermitian, (2, 3); hermitian_factors=true
        )
    end

    @testset "sparse, precision, and invalid boundaries" begin
        hermitian_float32 = Float32[
            1 2 0 0
            2 3 1 0
            0 1 4 2
            0 0 2 5
        ]
        repaired_float32 = QETOperatorSchmidt.operator_schmidt_decomposition(
            hermitian_float32, (2, 2); hermitian_factors=true
        )
        @test eltype(repaired_float32.coefficients) == Float32
        @test all(factor -> eltype(factor) == ComplexF32, repaired_float32.left_factors)
        @test all(factor -> eltype(factor) == ComplexF32, repaired_float32.right_factors)
        @test _wp2_operator_schmidt_reconstruction(repaired_float32) ≈ hermitian_float32 atol =
            2.0f-5
        @test all(ishermitian, repaired_float32.left_factors)
        @test all(ishermitian, repaired_float32.right_factors)

        sparse_operator = sparse(hermitian_float32)
        @test_throws ArgumentError QETOperatorSchmidt.operator_schmidt_decomposition(
            sparse_operator, (2, 2); hermitian_factors=true
        )
        sparse_repaired = QETOperatorSchmidt.operator_schmidt_decomposition(
            sparse_operator, (2, 2); allow_densify=true, hermitian_factors=true
        )
        @test _wp2_operator_schmidt_reconstruction(sparse_repaired) ≈ hermitian_float32 atol =
            2.0f-5
        @test all(ishermitian, sparse_repaired.left_factors)
        @test_throws ArgumentError CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            sparse_operator, (2, 2), -1
        )
        sparse_compatibility = CompatOperatorSchmidt.OperatorSchmidtDecomposition(
            sparse_operator, (2, 2), -1; allow_densify=true
        )
        @test all(ishermitian, sparse_compatibility.left_factors)
        @test all(ishermitian, sparse_compatibility.right_factors)

        nonfinite = Matrix{ComplexF64}(I, 4, 4)
        nonfinite[1, 1] = Inf
        @test_throws ArgumentError QETOperatorSchmidt.operator_schmidt_decomposition(
            nonfinite, (2, 2); hermitian_factors=true
        )
        @test_throws ArgumentError QETOperatorSchmidt.operator_schmidt_decomposition(
            Matrix{BigFloat}(I, 4, 4), (2, 2); hermitian_factors=true
        )

        legacy = QETOperatorSchmidt.OperatorSchmidtDecompositionResult(
            [1.0], [ones(1, 1)], [ones(1, 1)], (1, 1), (1, 1)
        )
        @test legacy.factor_convention === :general
        @test legacy.coordinate_imaginary_residual === nothing
        @test legacy.coordinate_imaginary_tolerance === nothing
    end
end
