using LinearAlgebra
using SparseArrays

const CompatProduct = QuantumEntanglementTools.MATLABCompat

@testset "Tier E product-analysis compatibility wrappers" begin
    @testset "operator Schmidt wrappers" begin
        e11 = [1.0 0.0; 0.0 0.0]
        e22 = [0.0 0.0; 0.0 1.0]
        basis_12 = zeros(3, 3)
        basis_12[1, 2] = 1
        basis_23 = zeros(3, 3)
        basis_23[2, 3] = 1
        operator = 3tensor_product(e11, basis_12) + 2tensor_product(e22, basis_23)

        nonzero = CompatProduct.OperatorSchmidtDecomposition(operator, (2, 3))
        full = CompatProduct.OperatorSchmidtDecomposition(operator, (2, 3), -1)
        leading = CompatProduct.OperatorSchmidtDecomposition(operator, (2, 3), 1)
        @test nonzero.coefficients ≈ [3.0, 2.0]
        @test full.coefficients ≈ [3.0, 2.0, 0.0, 0.0]
        @test leading.coefficients ≈ [3.0]
        @test tensor_sum(full.left_factors, full.right_factors; weights=full.coefficients) ≈
            operator
        @test CompatProduct.OperatorSchmidtRank(operator, (2, 3)) == 2

        rectangular_left = [1.0 2.0 0.0; 0.0 1.0 3.0]
        rectangular_right = [0.0 1.0; 2.0 3.0]
        rectangular = tensor_product(rectangular_left, rectangular_right)
        rectangular_dim = [2 2; 3 2]
        rectangular_result = CompatProduct.OperatorSchmidtDecomposition(
            rectangular, rectangular_dim, -1
        )
        @test tensor_sum(
            rectangular_result.left_factors,
            rectangular_result.right_factors;
            weights=rectangular_result.coefficients,
        ) ≈ rectangular
        @test CompatProduct.OperatorSchmidtRank(rectangular, rectangular_dim) == 1

        hermitian = Matrix(Diagonal(collect(1.0:6.0)))
        hermitian_result = CompatProduct.OperatorSchmidtDecomposition(hermitian, (2, 3), -1)
        @test tensor_sum(
            hermitian_result.left_factors,
            hermitian_result.right_factors;
            weights=hermitian_result.coefficients,
        ) ≈ hermitian

        sparse_operator = sparse(operator)
        @test_throws ArgumentError CompatProduct.OperatorSchmidtRank(
            sparse_operator, (2, 3)
        )
        @test CompatProduct.OperatorSchmidtRank(
            sparse_operator, (2, 3); allow_densify=true
        ) == 2
        @test_throws ArgumentError CompatProduct.OperatorSchmidtDecomposition(
            operator, (2, 3), -2
        )
        @test_throws ArgumentError CompatProduct.OperatorSchmidtDecomposition(
            operator, (2, 3), 5
        )
        @test_throws ArgumentError CompatProduct.OperatorSchmidtDecomposition(
            rectangular, (2, 2)
        )
        @test_throws DimensionMismatch CompatProduct.OperatorSchmidtRank(
            operator, (1, 1, 6)
        )
    end

    @testset "structured product wrappers" begin
        product_vector = tensor_product(
            ComplexF64[1, 2im], ComplexF64[3, -1], ComplexF64[2, im]
        )
        vector_result = CompatProduct.IsProductVector(product_vector, (2, 2, 2))
        @test vector_result.status === :within_tolerance
        @test tensor_product(vector_result.factors...) ≈ product_vector
        @test CompatProduct.IsProductVector(transpose(product_vector), (2, 2, 2)).status ===
            :within_tolerance

        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_result = CompatProduct.IsProductVector(bell, 2)
        @test bell_result.status === :outside_tolerance
        @test bell_result.reconstruction_residual > bell_result.approximation_threshold

        near_product = Float64[1, 0, 0, 1e-7]
        @test CompatProduct.IsProductVector(near_product, (2, 2); atol=1e-7, rtol=0).status ===
            :boundary

        sparse_vector = sparse(product_vector)
        @test_throws ArgumentError CompatProduct.IsProductVector(sparse_vector, (2, 2, 2))
        @test CompatProduct.IsProductVector(
            sparse_vector, (2, 2, 2); allow_densify=true
        ).status === :within_tolerance
        @test_throws DomainError CompatProduct.IsProductVector(zeros(4), (2, 2))
        @test_throws DimensionMismatch CompatProduct.IsProductVector(ones(4), (2, 3))

        first = [1.0 2.0; 3.0 4.0]
        second = [0.0 1.0 2.0; 3.0 4.0 5.0]
        third = reshape([1.0, 2.0], 2, 1)
        product_operator = tensor_product(first, second, third)
        product_dim = [2 2 2; 2 3 1]
        operator_result = CompatProduct.IsProductOperator(product_operator, product_dim)
        @test operator_result.status === :within_tolerance
        @test tensor_product(operator_result.factors...) ≈ product_operator

        e11 = [1.0 0.0; 0.0 0.0]
        e22 = [0.0 0.0; 0.0 1.0]
        rank_two = tensor_product(e11, e11) + tensor_product(e22, e22)
        @test CompatProduct.IsProductOperator(rank_two, (2, 2)).status ===
            :outside_tolerance
        near_operator = tensor_product(e11, e11) + 1e-7tensor_product(e22, e22)
        @test CompatProduct.IsProductOperator(near_operator, (2, 2); atol=1e-7, rtol=0).status ===
            :boundary

        @test_throws ArgumentError CompatProduct.IsProductOperator(
            sparse(product_operator), product_dim
        )
        @test CompatProduct.IsProductOperator(
            sparse(product_operator), product_dim; allow_densify=true
        ).status === :within_tolerance
        @test_throws DomainError CompatProduct.IsProductOperator(zeros(4, 4), (2, 2))
        @test_throws DimensionMismatch CompatProduct.IsProductOperator(
            product_operator, [2 2; 2 3]
        )
    end

    @testset "entanglement of formation wrapper" begin
        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        @test CompatProduct.EntFormation(bell, 2) ≈ 1
        @test CompatProduct.EntFormation(transpose(bell), (2, 2)) ≈ 1
        bell_density = bell * adjoint(bell)
        mixed_bell = 0.7bell_density + 0.3Matrix{ComplexF64}(I, 4, 4) / 4
        @test CompatProduct.EntFormation(mixed_bell, (2, 2)) ≈
            entanglement_of_formation(mixed_bell, (2, 2))
        @test CompatProduct.EntFormation(Matrix{Float64}(I, 4, 4) / 4, (2, 2)) == 0

        angle = 0.31
        rectangular_pure = [cos(angle), 0.0, 0.0, 0.0, sin(angle), 0.0]
        @test CompatProduct.EntFormation(rectangular_pure, (2, 3)) ≈
            entanglement_of_formation(rectangular_pure, (2, 3))
        @test CompatProduct.EntFormation(rectangular_pure) ≈
            entanglement_of_formation(rectangular_pure, (2, 3))
        @test CompatProduct.EntFormation(reshape(rectangular_pure, 1, :)) ≈
            entanglement_of_formation(rectangular_pure, (2, 3))

        exact_rectangular_density = zeros(6, 6)
        exact_rectangular_density[1, 1] = 0.5
        exact_rectangular_density[1, 5] = 0.5
        exact_rectangular_density[5, 1] = 0.5
        exact_rectangular_density[5, 5] = 0.5
        @test CompatProduct.EntFormation(exact_rectangular_density) ≈ 1
        numerical_rectangular_density = rectangular_pure * adjoint(rectangular_pure)
        @test_throws DomainError CompatProduct.EntFormation(
            numerical_rectangular_density, (2, 3)
        )
        @test CompatProduct.EntFormation(
            numerical_rectangular_density,
            (2, 3);
            psd_boundary_policy=:project,
            rank_boundary_policy=:project,
        ) ≈ entanglement_of_formation(rectangular_pure, (2, 3))
        @test_throws ArgumentError CompatProduct.EntFormation(
            Diagonal([0.5, 0.5, 0.0, 0.0, 0.0, 0.0]), (2, 3)
        )
        @test_throws ArgumentError CompatProduct.EntFormation(0.9bell, (2, 2))
        sparse_bell = sparsevec([1, 4], fill(inv(sqrt(2)), 2), 4)
        @test_throws ArgumentError CompatProduct.EntFormation(sparse_bell, (2, 2))
        @test CompatProduct.EntFormation(sparse_bell, (2, 2); allow_densify=true) ≈ 1
        sparse_rectangular_density = sparse(exact_rectangular_density)
        @test_throws ArgumentError CompatProduct.EntFormation(
            sparse_rectangular_density, (2, 3)
        )
        @test CompatProduct.EntFormation(
            sparse_rectangular_density, (2, 3); allow_densify=true
        ) ≈ 1
        @test_throws ArgumentError CompatProduct.EntFormation(
            bell; psd_boundary_policy=:project
        )
        @test_throws ArgumentError CompatProduct.EntFormation(
            bell_density; rank_boundary_policy=:unsupported
        )
    end

    @testset "structured separable-ball wrapper" begin
        unnormalized_identity = 2Matrix{Float64}(I, 4, 4)
        certified = CompatProduct.InSeparableBall(unnormalized_identity)
        @test certified.status === :separable_certified
        @test certified.purity == 0.25

        outside = CompatProduct.InSeparableBall([2.0, 0.0, 0.0, 0.0])
        @test outside.status === :outside_ball
        @test occursin("no entanglement conclusion", lowercase(outside.message))
        @test !occursin("entangled", lowercase(outside.message))

        boundary = CompatProduct.InSeparableBall([1 / 3, 1 / 3, 1 / 3, 0.0], (2, 2))
        @test boundary.status === :unknown
        @test_throws DomainError CompatProduct.InSeparableBall(zeros(4))
        @test_throws DomainError CompatProduct.InSeparableBall([-0.1, 0.4, 0.3, 0.4],)
        @test_throws DimensionMismatch CompatProduct.InSeparableBall(ones(2, 3))

        sparse_identity = sparse(unnormalized_identity)
        @test_throws ArgumentError CompatProduct.InSeparableBall(sparse_identity)
        @test CompatProduct.InSeparableBall(sparse_identity; allow_densify=true).status ===
            :separable_certified
    end
end
