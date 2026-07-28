using LinearAlgebra
using SparseArrays

const CompatMatrixPredicates = QuantumEntanglementTools.MATLABCompat

@testset "Tier E matrix-predicate compatibility" begin
    @testset "IsPSD preserves structured boundary semantics" begin
        positive = CompatMatrixPredicates.IsPSD([2.0 1.0; 1.0 2.0])
        @test positive isa MatrixPredicateResult
        @test positive.status === MatrixPredicateSatisfied

        boundary = CompatMatrixPredicates.IsPSD(zeros(2, 2))
        @test boundary.status === MatrixPredicateUnknown
        @test boundary.witness !== nothing

        exact_boundary = CompatMatrixPredicates.IsPSD(Rational{Int}[1 1; 1 1])
        @test exact_boundary.status === MatrixPredicateSatisfied
        @test exact_boundary.tolerance == 0

        indefinite = CompatMatrixPredicates.IsPSD([1.0 2.0; 2.0 1.0])
        @test indefinite.status === MatrixPredicateViolated

        nonhermitian = [1.0 1.0e-8; 0.0 1.0]
        saved = copy(nonhermitian)
        uncertain = CompatMatrixPredicates.IsPSD(nonhermitian, 2.0e-8)
        @test uncertain.status === MatrixPredicateUnknown
        @test uncertain.witness.kind === :hermiticity_boundary
        @test nonhermitian == saved
        @test CompatMatrixPredicates.IsPSD(nonhermitian, 1.0e-9).status ===
            MatrixPredicateViolated

        sparse_identity = sparse(Matrix{Float64}(I, 2, 2))
        @test_throws ArgumentError CompatMatrixPredicates.IsPSD(sparse_identity)
        @test CompatMatrixPredicates.IsPSD(sparse_identity; allow_densify=true).status ===
            MatrixPredicateSatisfied
        @test_throws ArgumentError CompatMatrixPredicates.IsPSD(
            Matrix{Float64}(I, 2, 2), -1
        )
    end

    @testset "IsLocallyPSD preserves indices and guard" begin
        local_positive = CompatMatrixPredicates.IsLocallyPSD(
            Rational{Int}[1 0 0; 0 2 0; 0 0 3], 2
        )
        @test local_positive.status === MatrixPredicateSatisfied
        @test local_positive.planned == 3

        local_violation = CompatMatrixPredicates.IsLocallyPSD(
            Rational{Int}[1 0 0; 0 2 0; 0 0 -1], 2
        )
        @test local_violation.status === MatrixPredicateViolated
        @test local_violation.witness.indices == (1, 3)

        local_boundary = CompatMatrixPredicates.IsLocallyPSD(zeros(3, 3), 1; atol=0)
        @test local_boundary.status === MatrixPredicateUnknown
        @test local_boundary.checked == 3

        @test_throws ArgumentError CompatMatrixPredicates.IsLocallyPSD(
            Matrix{Float64}(I, 20, 20), 10; max_submatrices=100
        )
        @test_throws ArgumentError CompatMatrixPredicates.IsLocallyPSD(
            sparse(Matrix{Float64}(I, 2, 2)), 1
        )
    end

    @testset "all-minor positional compatibility" begin
        pascal = Rational{Int}[1 1 1; 1 2 3; 1 3 6]
        positive = CompatMatrixPredicates.IsTotallyPositive(pascal, [1, 2, 3])
        @test positive.status === MatrixPredicateSatisfied
        @test positive.checked == 19

        order_two = CompatMatrixPredicates.IsTotallyPositive(pascal, 2)
        @test order_two.status === MatrixPredicateSatisfied
        @test order_two.planned == 9

        strict_zero = CompatMatrixPredicates.IsTotallyPositive(
            reshape([0.0], 1, 1), nothing, 1.0e-8
        )
        @test strict_zero.status === MatrixPredicateViolated

        positive_boundary = CompatMatrixPredicates.IsTotallyPositive(
            reshape([1.0e-9], 1, 1), nothing, 1.0e-8
        )
        @test positive_boundary.status === MatrixPredicateUnknown

        nonsingular = CompatMatrixPredicates.IsTotallyNonsingular(Rational{Int}[1 2; 3 5])
        @test nonsingular.status === MatrixPredicateSatisfied
        @test nonsingular.planned == 5

        singular = CompatMatrixPredicates.IsTotallyNonsingular(
            [1.0 2.0; 2.0 4.0], 2, 1.0e-8
        )
        @test singular.status === MatrixPredicateViolated
        @test singular.witness.order == 2

        near_singular = CompatMatrixPredicates.IsTotallyNonsingular(
            [1.0 1.0; 1.0 1.0 + 1.0e-12], 2, 1.0e-8
        )
        @test near_singular.status === MatrixPredicateUnknown

        sparse_positive = sparse([1.0 1.0; 1.0 2.0])
        @test_throws ArgumentError CompatMatrixPredicates.IsTotallyPositive(sparse_positive)
        @test CompatMatrixPredicates.IsTotallyPositive(
            sparse_positive; allow_densify=true
        ).status === MatrixPredicateSatisfied

        @test_throws ArgumentError CompatMatrixPredicates.IsTotallyNonsingular(
            ones(10, 10); max_minors=100
        )
        @test_throws ArgumentError CompatMatrixPredicates.IsTotallyPositive(
            ComplexF64[1 2; 3 4]
        )
    end
end
