using LinearAlgebra
using Random
using SparseArrays

struct _WP2CommutantZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_WP2CommutantZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_WP2CommutantZeroBasedMatrix)
    return ntuple(index -> 0:(size(matrix, index) - 1), 2)
end
Base.IndexStyle(::Type{<:_WP2CommutantZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_WP2CommutantZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function _wp2_commutant_vectors(basis)
    isempty(basis) && return zeros(0, 0)
    return reduce(hcat, vec.(basis))
end

function _wp2_commutant_projector(basis)
    vectors = _wp2_commutant_vectors(basis)
    return vectors * adjoint(vectors)
end

function _wp2_check_commutant_basis(generators, basis; atol)
    dimension = size(first(generators), 1)
    @test all(size(matrix) == (dimension, dimension) for matrix in basis)
    vectors = _wp2_commutant_vectors(basis)
    @test adjoint(vectors) * vectors ≈
        Matrix{eltype(vectors)}(I, length(basis), length(basis)) atol = atol rtol = atol
    for generator in generators, matrix in basis
        scale = max(norm(generator) * norm(matrix), one(atol))
        @test norm(generator * matrix - matrix * generator) <= atol * scale
    end
    return nothing
end

@testset "WP2 solver-free commutant" begin
    @testset "dimensions, residuals, and algebra structure" begin
        diagonal_generator = Diagonal([1.0, 2.0, 4.0])
        diagonal_basis = commutant(diagonal_generator)
        @test length(diagonal_basis) == 3
        _wp2_check_commutant_basis([diagonal_generator], diagonal_basis; atol=1e-12)

        full_basis = commutant(Matrix{Float64}(I, 3, 3))
        @test length(full_basis) == 9
        _wp2_check_commutant_basis([Matrix{Float64}(I, 3, 3)], full_basis; atol=1e-12)

        block_generator = Diagonal([1.0, 1.0, 2.0])
        block_basis = commutant(block_generator)
        @test length(block_basis) == 5
        _wp2_check_commutant_basis([block_generator], block_basis; atol=1e-12)

        matrix_units = [
            [1.0 0.0; 0.0 0.0], [0.0 1.0; 0.0 0.0], [0.0 0.0; 1.0 0.0], [0.0 0.0; 0.0 1.0]
        ]
        scalar_basis = commutant(matrix_units)
        @test length(scalar_basis) == 1
        _wp2_check_commutant_basis(matrix_units, scalar_basis; atol=1e-12)
        identity_vector = vec(Matrix{Float64}(I, 2, 2)) / sqrt(2)
        @test _wp2_commutant_projector(scalar_basis) ≈
            identity_vector * adjoint(identity_vector) atol = 1e-12

        pauli_x = ComplexF64[0 1; 1 0]
        pauli_z = ComplexF64[1 0; 0 -1]
        irreducible_basis = commutant((pauli_x, pauli_z))
        @test length(irreducible_basis) == 1
        _wp2_check_commutant_basis([pauli_x, pauli_z], irreducible_basis; atol=1e-12)

        distinct_diagonal = Diagonal([1.0, 2.0, 3.0])
        cyclic_shift = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
        multiple_basis = commutant([distinct_diagonal, cyclic_shift])
        @test length(multiple_basis) == 1
        _wp2_check_commutant_basis(
            [distinct_diagonal, cyclic_shift], multiple_basis; atol=1e-12
        )
    end

    @testset "basis-change and wrapper subspace invariance" begin
        rng = MersenneTwister(0x434f4d4d5554414e)
        generators = [ComplexF64[1 0 0; 0 1 0; 0 0 2], ComplexF64[0 1 0; 1 0 0; 0 0 3]]
        basis = commutant(generators)
        unitary = Matrix(qr(randn(rng, ComplexF64, 3, 3)).Q)
        rotated_generators = [unitary * matrix * adjoint(unitary) for matrix in generators]
        rotated_basis = commutant(rotated_generators)
        vectorized_conjugation = kron(conj(unitary), unitary)
        expected_projector =
            vectorized_conjugation *
            _wp2_commutant_projector(basis) *
            adjoint(vectorized_conjugation)
        @test _wp2_commutant_projector(rotated_basis) ≈ expected_projector atol = 2e-12
        _wp2_check_commutant_basis(rotated_generators, rotated_basis; atol=2e-12)

        compatibility_basis = MATLABCompat.Commutant(generators)
        @test length(compatibility_basis) == length(basis)
        @test _wp2_commutant_projector(compatibility_basis) ≈
            _wp2_commutant_projector(basis) atol = 2e-12
    end

    @testset "numeric types and explicit precision boundary" begin
        for value_type in (Float32, Float64, ComplexF32, ComplexF64)
            generator = value_type[1 2; 3 4]
            basis = commutant(generator)
            @test eltype(first(basis)) == value_type
            tolerance = value_type <: Union{Float32,ComplexF32} ? 2.0f-5 : 1e-12
            _wp2_check_commutant_basis([generator], basis; atol=tolerance)
        end

        integer_basis = commutant([1 2; 3 4])
        @test eltype(first(integer_basis)) == Float64
        @test_throws ArgumentError commutant(BigFloat[1 0; 0 2])
        @test_throws ArgumentError commutant(Complex{BigFloat}[1 0; 0 2])
    end

    @testset "sparse opt-in and allocation budgets" begin
        sparse_generator = spdiagm(0 => [1.0, 2.0, 3.0])
        @test_throws ArgumentError commutant(sparse_generator)
        sparse_basis = commutant(sparse_generator; allow_densify=true)
        @test length(sparse_basis) == 3
        @test all(matrix -> matrix isa Matrix{Float64}, sparse_basis)
        _wp2_check_commutant_basis([sparse_generator], sparse_basis; atol=1e-12)

        compatibility_basis = MATLABCompat.Commutant(sparse_generator; allow_densify=true)
        @test all(issparse, compatibility_basis)
        @test _wp2_commutant_projector(compatibility_basis) ≈
            _wp2_commutant_projector(sparse_basis) atol = 1e-12

        mixed_generators = [sparse_generator, Matrix(sparse_generator)]
        @test_throws ArgumentError commutant(mixed_generators)
        @test length(commutant(mixed_generators; allow_densify=true)) == 3

        @test_throws ArgumentError commutant(ones(2, 2); max_entries=15)
        @test_throws ArgumentError commutant(ones(2, 2); max_work=63)
        @test length(commutant(ones(2, 2); max_entries=nothing, max_work=nothing)) == 2
        for invalid_limit in (0, -1, true, 1.5)
            @test_throws ArgumentError commutant(ones(2, 2); max_entries=invalid_limit)
            @test_throws ArgumentError commutant(ones(2, 2); max_work=invalid_limit)
        end
    end

    @testset "invalid inputs and tolerances" begin
        @test_throws ArgumentError commutant(Matrix{Float64}[])
        @test_throws ArgumentError commutant(())
        @test_throws ArgumentError commutant(zeros(0, 0))
        @test_throws DimensionMismatch commutant(zeros(2, 3))
        @test_throws DimensionMismatch commutant([zeros(2, 2), zeros(3, 3)])
        @test_throws ArgumentError commutant([1.0, 2.0])
        @test_throws ArgumentError commutant("not a matrix")
        @test_throws ArgumentError commutant(Any[1.0 "x"; 0.0 1.0])
        @test_throws ArgumentError commutant([1.0 NaN; 0.0 1.0])
        @test_throws ArgumentError commutant(ComplexF64[1 Inf * im; 0 1])
        @test_throws ArgumentError commutant(
            _WP2CommutantZeroBasedMatrix([1.0 0.0; 0.0 2.0])
        )

        for invalid_tolerance in (-1, Inf, NaN, true, 1 + im)
            @test_throws ArgumentError commutant([1.0 0.0; 0.0 2.0]; atol=invalid_tolerance)
            @test_throws ArgumentError commutant([1.0 0.0; 0.0 2.0]; rtol=invalid_tolerance)
        end
    end
end
