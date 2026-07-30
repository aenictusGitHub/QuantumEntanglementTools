using LinearAlgebra
using SparseArrays

const QETEntangledSubspace = QuantumEntanglementTools
const CompatEntangledSubspace = QuantumEntanglementTools.MATLABCompat

struct _WP2ZeroBasedDimensionVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end
Base.size(vector::_WP2ZeroBasedDimensionVector) = size(vector.storage)
Base.axes(vector::_WP2ZeroBasedDimensionVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_WP2ZeroBasedDimensionVector}) = IndexLinear()
Base.getindex(vector::_WP2ZeroBasedDimensionVector, index::Int) = vector.storage[index + 1]

function _wp2_exact_rank(matrix)
    reduced = Rational{BigInt}.(matrix)
    rows, columns = size(reduced)
    pivot_row = 1
    for column in 1:columns
        pivot = findfirst(row -> !iszero(reduced[row, column]), pivot_row:rows)
        pivot === nothing && continue
        selected = pivot_row + pivot - 1
        if selected != pivot_row
            temporary = copy(reduced[pivot_row, :])
            reduced[pivot_row, :] .= reduced[selected, :]
            reduced[selected, :] .= temporary
        end
        pivot_value = reduced[pivot_row, column]
        reduced[pivot_row, :] ./= pivot_value
        for row in 1:rows
            row == pivot_row && continue
            factor = reduced[row, column]
            iszero(factor) && continue
            reduced[row, :] .-= factor .* reduced[pivot_row, :]
        end
        pivot_row += 1
        pivot_row > rows && break
    end
    return pivot_row - 1
end

function _wp2_minimum_schmidt_rank_on_grid(basis, dims, coefficient_values)
    minimum_rank = min(dims...)
    coefficient_ranges = ntuple(_ -> coefficient_values, size(basis, 2))
    for coefficients in Iterators.product(coefficient_ranges...)
        all(iszero, coefficients) && continue
        vector = basis * collect(coefficients)
        coefficient_matrix = reshape(vector, dims[2], dims[1])
        minimum_rank = min(minimum_rank, _wp2_exact_rank(coefficient_matrix))
    end
    return minimum_rank
end

@testset "WP2 entangled subspaces" begin
    @testset "sharp dimensions and exact Schmidt-rank property" begin
        basis_2x3 = QETEntangledSubspace.entangled_subspace(
            2, (2, 3); r=1, coefficient_type=BigInt
        )
        @test size(basis_2x3) == (6, 2)
        @test issparse(basis_2x3)
        @test eltype(basis_2x3) == BigInt
        @test _wp2_exact_rank(basis_2x3) == 2
        @test _wp2_minimum_schmidt_rank_on_grid(basis_2x3, (2, 3), -2:2) == 2

        basis_3x3 = QETEntangledSubspace.entangled_subspace(
            4, (3, 3); r=1, coefficient_type=Rational{BigInt}
        )
        @test size(basis_3x3) == (9, 4)
        @test _wp2_exact_rank(basis_3x3) == 4
        @test _wp2_minimum_schmidt_rank_on_grid(basis_3x3, (3, 3), -1:1) >= 2

        basis_rank_three = QETEntangledSubspace.entangled_subspace(
            2, (3, 4); r=2, coefficient_type=BigInt
        )
        @test size(basis_rank_three) == (12, 2)
        @test _wp2_minimum_schmidt_rank_on_grid(basis_rank_three, (3, 4), -2:2) == 3
    end

    @testset "construction order, prefixes, and scalar dimensions" begin
        full = QETEntangledSubspace.entangled_subspace(
            6, (2, 3); r=0, coefficient_type=BigInt
        )
        @test size(full) == (6, 6)
        @test _wp2_exact_rank(full) == 6

        maximum = QETEntangledSubspace.entangled_subspace(
            3, (2, 4); r=1, coefficient_type=BigInt
        )
        prefix = QETEntangledSubspace.entangled_subspace(
            2, (2, 4); r=1, coefficient_type=BigInt
        )
        @test prefix == maximum[:, 1:2]
        @test maximum[:, 1] == sparsevec([3, 8], BigInt[1, 1], 8)

        square = QETEntangledSubspace.entangled_subspace(4, 3; coefficient_type=BigInt)
        @test size(square) == (9, 4)
        @test _wp2_minimum_schmidt_rank_on_grid(square, (3, 3), -1:1) >= 2
    end

    @testset "types, wrapper, and explicit budgets" begin
        float32_basis = QETEntangledSubspace.entangled_subspace(
            2, (2, 3); coefficient_type=Float32
        )
        @test eltype(float32_basis) == Float32
        big_basis = QETEntangledSubspace.entangled_subspace(
            2, (2, 3); coefficient_type=BigFloat
        )
        @test eltype(big_basis) == BigFloat
        complex_basis = QETEntangledSubspace.entangled_subspace(
            2, (2, 3); coefficient_type=ComplexF64
        )
        @test complex_basis == ComplexF64.(float32_basis)

        compatibility = CompatEntangledSubspace.EntangledSubspace(2, [2, 3], 1)
        @test compatibility == Float64.(float32_basis)
        @test issparse(compatibility)

        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            2, (2, 3); max_nonzeros=3
        )
        @test nnz(QETEntangledSubspace.entangled_subspace(2, (2, 3); max_nonzeros=4)) == 4
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            2, (2, 3); max_work=3
        )
        @test size(
            QETEntangledSubspace.entangled_subspace(
                2, (2, 3); max_nonzeros=nothing, max_work=nothing
            ),
        ) == (6, 2)
    end

    @testset "validation and overflow" begin
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(0, (2, 2))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(-1, (2, 2))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(true, (2, 2))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2.0, (2, 2))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2, (2,))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2, (2, 3, 4))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2, (2.0, 3))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2, true)
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(2, (0, 3))
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            2, _WP2ZeroBasedDimensionVector([2, 3])
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(3, (2, 3); r=1)
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(1, (2, 3); r=2)
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(1, (2, 3); r=-1)
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (2, 3); r=true
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (2, 3); coefficient_type=Number
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (2, 3); coefficient_type=Bool
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (2, 3); max_nonzeros=0
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (2, 3); max_work=true
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            1, (typemax(Int), 2); r=1, max_nonzeros=nothing, max_work=nothing
        )
        @test_throws ArgumentError QETEntangledSubspace.entangled_subspace(
            121,
            (12, 12);
            r=1,
            coefficient_type=Int8,
            max_nonzeros=nothing,
            max_work=nothing,
        )
    end
end
