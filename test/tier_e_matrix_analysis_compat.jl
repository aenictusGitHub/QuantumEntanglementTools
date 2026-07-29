using LinearAlgebra
using SparseArrays
using Test

const CompatMatrixAnalysis = QuantumEntanglementTools.MATLABCompat

struct _TierEMatrixCompatZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_TierEMatrixCompatZeroBasedVector) = size(vector.storage)
Base.axes(vector::_TierEMatrixCompatZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_TierEMatrixCompatZeroBasedVector}) = IndexLinear()
function Base.getindex(vector::_TierEMatrixCompatZeroBasedVector, index::Int)
    return vector.storage[index + 1]
end

struct _TierEMatrixCompatZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_TierEMatrixCompatZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_TierEMatrixCompatZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_TierEMatrixCompatZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierEMatrixCompatZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

@testset "Tier E matrix-analysis compatibility" begin
    @testset "array axes validation" begin
        vector = _TierEMatrixCompatZeroBasedVector([2.0, 1.0])
        row_matrix = _TierEMatrixCompatZeroBasedMatrix(reshape([2.0, 1.0], 1, :))
        @test_throws ArgumentError CompatMatrixAnalysis.Majorizes(vector, [1.5, 1.5])
        @test_throws ArgumentError CompatMatrixAnalysis.ElemSymPoly(row_matrix, 1)
    end

    @testset "Majorizes preserves pinned weak semantics" begin
        @test CompatMatrixAnalysis.Majorizes([4, 1, 1], [3, 2, 1]; rtol=0)
        @test CompatMatrixAnalysis.Majorizes([2, 0], [1, 0]; rtol=0)

        # Pinned QETLAB sorts before appending zeros and treats row matrices
        # as vectors. The native strong-majorization results are opposite.
        @test !CompatMatrixAnalysis.Majorizes([1, -1], [1, 0, -1]; rtol=0)
        @test !CompatMatrixAnalysis.Majorizes([1.0 1.0], [sqrt(2.0)]; atol=1e-14, rtol=0)
        @test majorizes([1, -1], [1, 0, -1]; rtol=0)
        @test majorizes([1.0 1.0], [sqrt(2.0)]; atol=1e-14, rtol=0)

        delta = 1e-8
        first = [0.6 - delta, 0.3 + delta, 0.1]
        second = [0.6, 0.3, 0.1]
        @test CompatMatrixAnalysis.Majorizes(first, second; atol=2delta, rtol=0)
        @test CompatMatrixAnalysis.Majorizes(second, first; atol=2delta, rtol=0)
        @test_throws ArgumentError CompatMatrixAnalysis.Majorizes(
            sparsevec([1], [1.0], 2), [1.0, 0.0]
        )
        @test CompatMatrixAnalysis.Majorizes(
            sparsevec([1], [1.0], 2), [1.0, 0.0]; rtol=0, allow_densify=true
        )
    end

    @testset "ElemSymPoly delegates with exact arithmetic" begin
        @test CompatMatrixAnalysis.ElemSymPoly([1, 2, 3, 4], 2) == 35
        @test CompatMatrixAnalysis.ElemSymPoly(reshape([1, 2, 3, 4], 1, :), 3) == 50
        @test CompatMatrixAnalysis.ElemSymPoly(reshape([1, 2, 3, 4], :, 1), 4) == 24

        rationals = Rational{Int}[1 // 2, 2 // 3, 3 // 4]
        result = CompatMatrixAnalysis.ElemSymPoly(rationals, 2)
        @test result == 29 // 24
        @test result isa Rational{Int}
        @test CompatMatrixAnalysis.ElemSymPoly(
            sparsevec([1, 4], Rational{Int}[2 // 3, -3 // 2], 5), 2
        ) == -1
        @test_throws DimensionMismatch CompatMatrixAnalysis.ElemSymPoly([1 2; 3 4], 1)
    end

    @testset "CompoundMatrix preserves pinned high-order shape" begin
        matrix = [1 2 3; 4 5 6; 7 8 10]
        @test CompatMatrixAnalysis.CompoundMatrix(matrix, 2) == compound_matrix(matrix, 2)
        @test CompatMatrixAnalysis.CompoundMatrix(matrix, 0) == reshape([1], 1, 1)
        @test size(CompatMatrixAnalysis.CompoundMatrix(zeros(Int, 2, 3), 3)) == (0, 0)
        @test size(compound_matrix(zeros(Int, 2, 3), 3)) == (0, 1)

        rational_matrix = Rational{Int}[1//2 2//3; -3//4 5//7]
        rational_result = CompatMatrixAnalysis.CompoundMatrix(rational_matrix, 2)
        @test rational_result[1, 1] == det(rational_matrix)
        @test eltype(rational_result) == Rational{Int}
        @test issparse(
            CompatMatrixAnalysis.CompoundMatrix(sparse(matrix), 2; sparse_output=true)
        )
    end

    @testset "AdditiveCompoundMatrix preserves pinned order-zero error" begin
        matrix = [1 2 3; 4 5 6; 7 8 10]
        @test CompatMatrixAnalysis.AdditiveCompoundMatrix(matrix, 2) ==
            additive_compound_matrix(matrix, 2)
        @test_throws ArgumentError CompatMatrixAnalysis.AdditiveCompoundMatrix(matrix, 0)
        @test additive_compound_matrix(matrix, 0) == reshape([0], 1, 1)
        @test size(CompatMatrixAnalysis.AdditiveCompoundMatrix(matrix, 4)) == (0, 0)

        rational_matrix = Rational{Int}[1//2 2//3; -3//4 5//7]
        rational_result = CompatMatrixAnalysis.AdditiveCompoundMatrix(rational_matrix, 2)
        @test rational_result == reshape([tr(rational_matrix)], 1, 1)
        @test eltype(rational_result) == Rational{Int}
        @test issparse(
            CompatMatrixAnalysis.AdditiveCompoundMatrix(
                sparse(matrix), 2; sparse_output=true
            ),
        )
        @test_throws DimensionMismatch CompatMatrixAnalysis.AdditiveCompoundMatrix(
            zeros(2, 3), 1
        )
    end
end
