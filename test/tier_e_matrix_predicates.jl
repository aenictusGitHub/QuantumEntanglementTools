using LinearAlgebra
using SparseArrays

struct _TierEMatrixPredicateZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_TierEMatrixPredicateZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_TierEMatrixPredicateZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_TierEMatrixPredicateZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierEMatrixPredicateZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

@testset "Tier E native matrix predicates" begin
    @testset "array axes validation" begin
        matrix = _TierEMatrixPredicateZeroBasedMatrix([1.0 0.0; 0.0 2.0])
        @test_throws ArgumentError is_positive_semidefinite(matrix)
        @test_throws ArgumentError is_locally_positive_semidefinite(matrix, 1)
        @test_throws ArgumentError is_totally_positive(matrix)
        @test_throws ArgumentError is_totally_nonsingular(matrix)
    end

    @testset "structured result vocabulary" begin
        result = is_positive_semidefinite([1.0 0.0; 0.0 2.0])
        @test result isa MatrixPredicateResult
        @test result.predicate === :positive_semidefinite
        @test result.status === MatrixPredicateSatisfied
        @test result.value ≈ 1.0
        @test result.tolerance > 0
        @test result.witness === nothing
        @test result.checked == 1
        @test result.planned == 1
        @test occursin("MatrixPredicateResult", sprint(show, result))
        @test MatrixPredicateSatisfied != MatrixPredicateViolated
        @test MatrixPredicateUnknown != MatrixPredicateSatisfied
    end

    @testset "positive semidefinite floating arithmetic" begin
        for value_type in (Float32, Float64)
            identity_matrix = Matrix{value_type}(I, 3, 3)
            satisfied = is_positive_semidefinite(identity_matrix)
            @test satisfied.status === MatrixPredicateSatisfied
            @test satisfied.value == one(value_type)

            indefinite = value_type[1 0 0; 0 -1 0; 0 0 2]
            violated = is_positive_semidefinite(indefinite)
            @test violated.status === MatrixPredicateViolated
            @test violated.value < 0
            @test violated.witness !== nothing
            @test real(dot(violated.witness, indefinite * violated.witness)) < 0

            tolerance = value_type(1e-3)
            boundary = is_positive_semidefinite(
                Diagonal(value_type[tolerance, 1]); atol=tolerance, rtol=zero(value_type)
            )
            @test boundary.status === MatrixPredicateUnknown
            @test boundary.value == tolerance
            @test boundary.witness.index == 1

            positive_margin = is_positive_semidefinite(
                Diagonal(value_type[2tolerance, 1]); atol=tolerance, rtol=zero(value_type)
            )
            @test positive_margin.status === MatrixPredicateSatisfied

            negative_margin = is_positive_semidefinite(
                Diagonal(value_type[-2tolerance, 1]); atol=tolerance, rtol=zero(value_type)
            )
            @test negative_margin.status === MatrixPredicateViolated
        end

        hermitian_complex = ComplexF64[2 im; -im 2]
        @test is_positive_semidefinite(hermitian_complex).status ===
            MatrixPredicateSatisfied
        complex_indefinite = ComplexF64[1 2im; -2im 1]
        @test is_positive_semidefinite(complex_indefinite).status ===
            MatrixPredicateViolated

        near_hermitian = [1.0 1.0e-8; 0.0 1.0]
        saved = copy(near_hermitian)
        uncertain = is_positive_semidefinite(near_hermitian; atol=2.0e-8, rtol=0.0)
        @test uncertain.status === MatrixPredicateUnknown
        @test uncertain.witness.kind === :hermiticity_boundary
        @test near_hermitian == saved

        rejected = is_positive_semidefinite(near_hermitian; atol=1.0e-9, rtol=0.0)
        @test rejected.status === MatrixPredicateViolated
        @test rejected.witness.kind === :nonhermitian
    end

    @testset "positive semidefinite exact and generic arithmetic" begin
        rational_semidefinite = Rational{Int}[1 1; 1 1]
        exact_result = is_positive_semidefinite(rational_semidefinite)
        @test exact_result.status === MatrixPredicateSatisfied
        @test exact_result.tolerance == 0

        zero_pivot_indefinite = Rational{Int}[0 1; 1 0]
        zero_pivot_result = is_positive_semidefinite(zero_pivot_indefinite)
        @test zero_pivot_result.status === MatrixPredicateViolated
        @test zero_pivot_result.value == -1
        @test zero_pivot_result.witness.kind === :zero_pivot_nonzero_schur_entry

        exact_nonhermitian = Rational{Int}[1 1; 0 1]
        nonhermitian_result = is_positive_semidefinite(exact_nonhermitian)
        @test nonhermitian_result.status === MatrixPredicateViolated
        @test nonhermitian_result.witness.kind === :nonhermitian

        large = typemax(Int)
        large_exact = [large 0; 0 large - 1]
        @test is_positive_semidefinite(large_exact).status === MatrixPredicateSatisfied

        big_positive = BigFloat[2 1; 1 2]
        big_result = is_positive_semidefinite(big_positive)
        @test big_result.status === MatrixPredicateSatisfied
        @test big_result.value == BigFloat("1.5")

        big_indefinite = BigFloat[1 2; 2 1]
        @test is_positive_semidefinite(big_indefinite).status === MatrixPredicateViolated

        big_boundary = is_positive_semidefinite(BigFloat[1 1; 1 1])
        @test big_boundary.status === MatrixPredicateUnknown
        @test iszero(big_boundary.value)

        diagonal_storage = Diagonal(BigFloat[0, 1, 2])
        diagonal_result = is_positive_semidefinite(diagonal_storage)
        @test diagonal_result.status === MatrixPredicateUnknown
        @test diagonal_result.witness.index == 1
        @test diagonal_storage.diag == BigFloat[0, 1, 2]
    end

    @testset "positive semidefinite validation and sparse policy" begin
        @test_throws DimensionMismatch is_positive_semidefinite(zeros(2, 3))
        @test_throws ArgumentError is_positive_semidefinite(zeros(0, 0))
        @test_throws ArgumentError is_positive_semidefinite([1.0 NaN; NaN 1.0])
        @test_throws ArgumentError is_positive_semidefinite([1.0 Inf; Inf 1.0])
        @test_throws ArgumentError is_positive_semidefinite(
            Matrix{Float64}(I, 2, 2); atol=-1
        )
        @test_throws ArgumentError is_positive_semidefinite(
            Matrix{Float64}(I, 2, 2); rtol=Inf
        )
        @test_throws ArgumentError is_positive_semidefinite(
            Matrix{Float64}(I, 2, 2); atol=true
        )
        @test_throws ArgumentError is_positive_semidefinite(
            Rational{Int}[1 0; 0 1]; atol=1.0e-12
        )

        sparse_identity = sparse(Matrix{Float64}(I, 3, 3))
        sparse_copy = copy(sparse_identity)
        @test_throws ArgumentError is_positive_semidefinite(sparse_identity)
        sparse_result = is_positive_semidefinite(sparse_identity; allow_densify=true)
        @test sparse_result.status === MatrixPredicateSatisfied
        @test sparse_identity == sparse_copy
        @test SparseArrays.issparse(sparse_identity)
    end

    @testset "local positive semidefiniteness" begin
        positive_diagonal = Rational{Int}[
            1 0 0
            0 2 0
            0 0 3
        ]
        local_result = is_locally_positive_semidefinite(positive_diagonal, 2)
        @test local_result.status === MatrixPredicateSatisfied
        @test local_result.checked == 3
        @test local_result.planned == 3

        negative_diagonal = Rational{Int}[
            1 0 0
            0 2 0
            0 0 -1
        ]
        local_violation = is_locally_positive_semidefinite(negative_diagonal, 2)
        @test local_violation.status === MatrixPredicateViolated
        @test local_violation.witness.indices == (1, 3)
        @test local_violation.witness.subresult.status === MatrixPredicateViolated
        @test local_violation.checked == 2
        @test local_violation.planned == 3

        nonsymmetric = Rational{Int}[1 2; 3 1]
        @test is_locally_positive_semidefinite(nonsymmetric, 1).status ===
            MatrixPredicateSatisfied
        @test is_locally_positive_semidefinite(nonsymmetric, 2).status ===
            MatrixPredicateViolated

        boundary_then_violation = Diagonal([0.0, -1.0])
        precedence = is_locally_positive_semidefinite(
            boundary_then_violation, 1; atol=0.0, rtol=0.0
        )
        @test precedence.status === MatrixPredicateViolated
        @test precedence.checked == 2
        @test precedence.witness.indices == (2,)

        all_boundary = is_locally_positive_semidefinite(zeros(3, 3), 1; atol=0.0, rtol=0.0)
        @test all_boundary.status === MatrixPredicateUnknown
        @test all_boundary.checked == 3
        @test all_boundary.witness.indices == (1,)

        full_order = is_locally_positive_semidefinite(Matrix{Float64}(I, 4, 4), 4)
        @test full_order.status === MatrixPredicateSatisfied
        @test full_order.planned == 1

        @test_throws DimensionMismatch is_locally_positive_semidefinite(zeros(2, 3), 1)
        @test_throws ArgumentError is_locally_positive_semidefinite(zeros(0, 0), 1)
        @test_throws ArgumentError is_locally_positive_semidefinite(
            Matrix{Float64}(I, 3, 3), 0
        )
        @test_throws ArgumentError is_locally_positive_semidefinite(
            Matrix{Float64}(I, 3, 3), 4
        )
        @test_throws ArgumentError is_locally_positive_semidefinite(
            Matrix{Float64}(I, 3, 3), true
        )
        @test_throws ArgumentError is_locally_positive_semidefinite(
            Matrix{Float64}(I, 20, 20), 10; max_submatrices=100
        )

        sparse_local = sparse(Matrix{Float64}(I, 3, 3))
        @test_throws ArgumentError is_locally_positive_semidefinite(sparse_local, 2)
        @test is_locally_positive_semidefinite(sparse_local, 2; allow_densify=true).status ===
            MatrixPredicateSatisfied
    end

    @testset "total positivity" begin
        pascal = Rational{Int}[
            1 1 1
            1 2 3
            1 3 6
        ]
        exact_positive = is_totally_positive(pascal)
        @test exact_positive.status === MatrixPredicateSatisfied
        @test exact_positive.checked == 19
        @test exact_positive.planned == 19

        row_scaling = Diagonal(Rational{Int}[2, 3, 5])
        column_scaling = Diagonal(Rational{Int}[7, 11, 13])
        scaled_pascal = Matrix(row_scaling * pascal * column_scaling)
        @test is_totally_positive(scaled_pascal).status === MatrixPredicateSatisfied

        rectangular_cauchy = Rational{Int}[
            1//2 1//3 1//4
            1//3 1//4 1//5
        ]
        rectangular_result = is_totally_positive(rectangular_cauchy)
        @test rectangular_result.status === MatrixPredicateSatisfied
        @test rectangular_result.planned == 9

        zero_minor = Rational{Int}[1 1; 1 1]
        zero_result = is_totally_positive(zero_minor)
        @test zero_result.status === MatrixPredicateViolated
        @test zero_result.value == 0
        @test zero_result.witness.order == 2
        @test zero_result.witness.rows == (1, 2)
        @test zero_result.witness.columns == (1, 2)

        negative_entry = Rational{Int}[1 -1; 1 2]
        negative_result = is_totally_positive(negative_entry)
        @test negative_result.status === MatrixPredicateViolated
        @test negative_result.witness.order == 1
        @test negative_result.value == -1

        large = typemax(Int)
        large_positive = [large large - 1; large - 1 large]
        @test is_totally_positive(large_positive).status === MatrixPredicateSatisfied

        second_order_only = is_totally_positive(pascal; orders=2)
        @test second_order_only.status === MatrixPredicateSatisfied
        @test second_order_only.checked == 9
        @test second_order_only.planned == 9

        for value_type in (Float32, Float64, BigFloat)
            floating_pascal = value_type.(pascal)
            @test is_totally_positive(floating_pascal).status === MatrixPredicateSatisfied
        end

        positive_boundary = is_totally_positive(
            reshape([1.0e-9], 1, 1); atol=1.0e-8, rtol=0.0
        )
        @test positive_boundary.status === MatrixPredicateUnknown
        @test positive_boundary.value == 1.0e-9

        negative_boundary = is_totally_positive(
            reshape([-1.0e-9], 1, 1); atol=1.0e-8, rtol=0.0
        )
        @test negative_boundary.status === MatrixPredicateUnknown

        exact_zero_float = is_totally_positive(reshape([0.0], 1, 1); atol=1.0e-8, rtol=0.0)
        @test exact_zero_float.status === MatrixPredicateViolated

        robust_negative = is_totally_positive(
            reshape([-2.0e-8], 1, 1); atol=1.0e-8, rtol=0.0
        )
        @test robust_negative.status === MatrixPredicateViolated
    end

    @testset "total positivity validation and complexity guard" begin
        @test_throws ArgumentError is_totally_positive(zeros(0, 2))
        @test_throws ArgumentError is_totally_positive(ComplexF64[1 2; 3 4])
        @test_throws ArgumentError is_totally_positive([1.0 NaN; 2.0 3.0])
        @test_throws ArgumentError is_totally_positive([1.0 2.0; 3.0 4.0]; orders=0)
        @test_throws ArgumentError is_totally_positive([1.0 2.0; 3.0 4.0]; orders=3)
        @test_throws ArgumentError is_totally_positive(
            [1.0 2.0; 3.0 4.0]; orders=Bool[true]
        )
        @test_throws ArgumentError is_totally_positive([1.0 2.0; 3.0 4.0]; orders=Int[])
        @test_throws ArgumentError is_totally_positive(ones(10, 10); max_minors=100)
        @test_throws ArgumentError is_totally_positive(ones(2, 2); max_minors=0)
        @test_throws ArgumentError is_totally_positive(
            Rational{Int}[1 1; 1 2]; rtol=1.0e-12
        )

        sparse_positive = sparse([1.0 1.0; 1.0 2.0])
        sparse_copy = copy(sparse_positive)
        @test_throws ArgumentError is_totally_positive(sparse_positive)
        @test is_totally_positive(sparse_positive; allow_densify=true).status ===
            MatrixPredicateSatisfied
        @test sparse_positive == sparse_copy

        diagonal_structured = Diagonal(Rational{Int}[1, 2])
        structured_result = is_totally_positive(diagonal_structured)
        @test structured_result.status === MatrixPredicateViolated
        @test structured_result.witness.order == 1
    end

    @testset "total nonsingularity" begin
        rectangular_cauchy = Rational{Int}[
            1//2 1//3 1//4
            1//3 1//4 1//5
        ]
        exact_result = is_totally_nonsingular(rectangular_cauchy)
        @test exact_result.status === MatrixPredicateSatisfied
        @test exact_result.checked == 9
        @test exact_result.planned == 9

        row_scaling = Diagonal(Rational{Int}[2, 3])
        column_scaling = Diagonal(Rational{Int}[5, 7, 11])
        scaled = Matrix(row_scaling * rectangular_cauchy * column_scaling)
        @test is_totally_nonsingular(scaled).status === MatrixPredicateSatisfied

        singular_exact = Rational{Int}[1 2; 2 4]
        singular_result = is_totally_nonsingular(singular_exact; orders=2)
        @test singular_result.status === MatrixPredicateViolated
        @test singular_result.value == 0
        @test singular_result.witness.order == 2

        full_order_only = is_totally_nonsingular(Rational{Int}[0 1; 1 0]; orders=2)
        @test full_order_only.status === MatrixPredicateSatisfied
        @test full_order_only.planned == 1
        @test is_totally_nonsingular(Rational{Int}[0 1; 1 0]).status ===
            MatrixPredicateViolated

        large = typemax(Int)
        large_nonsingular = [large large - 1; large - 1 large - 2]
        @test is_totally_nonsingular(large_nonsingular).status === MatrixPredicateSatisfied

        for value_type in (Float32, Float64, BigFloat)
            floating = value_type.(Rational{Int}[1 2; 3 5])
            @test is_totally_nonsingular(floating).status === MatrixPredicateSatisfied
        end

        complex_result = is_totally_nonsingular(ComplexF64[1 im; 2 3])
        @test complex_result.status === MatrixPredicateSatisfied

        near_singular = is_totally_nonsingular(
            [1.0 1.0; 1.0 1.0 + 1.0e-12]; orders=2, atol=1.0e-8, rtol=0.0
        )
        @test near_singular.status === MatrixPredicateUnknown
        @test near_singular.witness.order == 2

        represented_singular = is_totally_nonsingular(
            [1.0 2.0; 2.0 4.0]; orders=2, atol=1.0e-8, rtol=0.0
        )
        @test represented_singular.status === MatrixPredicateViolated

        zero_entry = is_totally_nonsingular(reshape([0.0], 1, 1); atol=1.0e-8, rtol=0.0)
        @test zero_entry.status === MatrixPredicateViolated

        nonzero_boundary = is_totally_nonsingular(
            reshape([1.0e-9], 1, 1); atol=1.0e-8, rtol=0.0
        )
        @test nonzero_boundary.status === MatrixPredicateUnknown
    end

    @testset "total nonsingularity validation and sparse policy" begin
        @test_throws ArgumentError is_totally_nonsingular(zeros(2, 0))
        @test_throws ArgumentError is_totally_nonsingular([1.0 Inf; 2.0 3.0])
        @test_throws ArgumentError is_totally_nonsingular(ones(12, 12); max_minors=100)
        @test_throws ArgumentError is_totally_nonsingular(ones(2, 2); orders=(1, 3))
        @test_throws ArgumentError is_totally_nonsingular(ones(2, 2); atol=-1)

        sparse_matrix = sparse([1.0 2.0; 3.0 5.0])
        sparse_copy = copy(sparse_matrix)
        @test_throws ArgumentError is_totally_nonsingular(sparse_matrix)
        sparse_result = is_totally_nonsingular(sparse_matrix; allow_densify=true)
        @test sparse_result.status === MatrixPredicateSatisfied
        @test sparse_matrix == sparse_copy
        @test SparseArrays.issparse(sparse_matrix)

        guarded_override = is_totally_nonsingular(
            Rational{Int}[1 2; 3 5]; max_minors=nothing
        )
        @test guarded_override.status === MatrixPredicateSatisfied
        @test guarded_override.planned == 5
    end
end
