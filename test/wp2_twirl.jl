using LinearAlgebra
using Random
using SparseArrays
using Test

const TwirlCompat = QuantumEntanglementTools.MATLABCompat

struct _WP2TwirlZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_WP2TwirlZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_WP2TwirlZeroBasedMatrix)
    return ntuple(index -> 0:(size(matrix, index) - 1), 2)
end
Base.IndexStyle(::Type{<:_WP2TwirlZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_WP2TwirlZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function _wp2_twirl_random_positive(rng, dimension)
    factor = randn(rng, ComplexF64, dimension, dimension)
    return factor * adjoint(factor)
end

function _wp2_twirl_tensor_power(matrix, copies)
    return foldl(kron, ntuple(_ -> matrix, copies))
end

function _wp2_direct_project(operator, basis)
    gram = [tr(adjoint(left) * right) for left in basis, right in basis]
    overlaps = [tr(adjoint(matrix) * operator) for matrix in basis]
    coefficients = gram \ overlaps
    return sum(coefficient * matrix for (coefficient, matrix) in zip(coefficients, basis))
end

@testset "WP2 guarded group twirls" begin
    rng = MersenneTwister(0x545749524c)

    @testset "bipartite projector formulas" begin
        dimension = 3
        total = dimension^2
        operator = randn(rng, ComplexF64, total, total)
        symmetric = Matrix(symmetric_projector(dimension, 2))
        antisymmetric = Matrix(antisymmetric_projector(dimension, 2))
        expected_werner =
            (tr(symmetric * operator) / tr(symmetric)) * symmetric +
            (tr(antisymmetric * operator) / tr(antisymmetric)) * antisymmetric
        actual_werner = twirl(operator; kind=:werner)
        @test actual_werner ≈ expected_werner atol = 2e-12 rtol = 2e-12

        entangled = maximally_entangled(dimension)
        entangled_projector = entangled * adjoint(entangled)
        overlap = dot(entangled, operator * entangled)
        expected_isotropic =
            overlap * entangled_projector +
            (tr(operator) - overlap) *
            (Matrix{ComplexF64}(I, total, total) - entangled_projector) / (total - 1)
        actual_isotropic = twirl(operator; kind=:isotropic)
        @test actual_isotropic ≈ expected_isotropic atol = 2e-12 rtol = 2e-12

        identity_matrix = Matrix{ComplexF64}(I, total, total)
        swap = Matrix(swap_operator((dimension, dimension); T=ComplexF64))
        omega = vec(Matrix{ComplexF64}(I, dimension, dimension))
        omega_projector = omega * adjoint(omega)
        expected_real = _wp2_direct_project(
            operator, [identity_matrix, swap, omega_projector]
        )
        actual_real = twirl(operator; kind=:real)
        @test actual_real ≈ expected_real atol = 3e-12 rtol = 3e-12
    end

    @testset "Pauli finite-group agreement" begin
        operator = randn(rng, ComplexF64, 4, 4)
        direct = zeros(ComplexF64, 4, 4)
        for label in 0:3
            pauli_operator = pauli(label)
            bilateral = kron(pauli_operator, pauli_operator)
            direct += bilateral * operator * adjoint(bilateral)
        end
        direct /= 4
        actual = twirl(operator; kind=:pauli)
        @test actual ≈ direct atol = 1e-13 rtol = 1e-13
        @test TwirlCompat.Twirl(operator, "PaUlI", 2) ≈ direct atol = 1e-13 rtol = 1e-13
        @test TwirlCompat.Twirl(operator, :PAULI, 2) ≈ direct atol = 1e-13 rtol = 1e-13
    end

    @testset "idempotence, trace, positivity, and group fixed points" begin
        for kind in (:werner, :isotropic, :real, :pauli)
            positive = _wp2_twirl_random_positive(rng, 4)
            projected = twirl(positive; kind)
            @test twirl(projected; kind) ≈ projected atol = 4e-12 rtol = 4e-12
            @test tr(projected) ≈ tr(positive) atol = 4e-12 rtol = 4e-12
            @test ishermitian(projected)
            @test minimum(eigvals(Hermitian(projected))) >= -4e-12
        end

        for copies in (2, 3)
            local_dimension = 2
            total = local_dimension^copies
            positive = _wp2_twirl_random_positive(rng, total)

            werner = twirl(positive; kind=:werner, copies)
            unitary = Matrix(qr(randn(rng, ComplexF64, 2, 2)).Q)
            unitary_power = _wp2_twirl_tensor_power(unitary, copies)
            @test unitary_power * werner * adjoint(unitary_power) ≈ werner atol = 8e-12 rtol =
                8e-12
            @test tr(werner) ≈ tr(positive) atol = 8e-12 rtol = 8e-12
            @test minimum(eigvals(Hermitian(werner))) >= -8e-12

            real_twirl = twirl(positive; kind=:real, copies)
            orthogonal = Matrix(qr(randn(rng, 2, 2)).Q)
            orthogonal_power = _wp2_twirl_tensor_power(orthogonal, copies)
            @test orthogonal_power * real_twirl * transpose(orthogonal_power) ≈ real_twirl atol =
                8e-12 rtol = 8e-12
            @test tr(real_twirl) ≈ tr(positive) atol = 8e-12 rtol = 8e-12
            @test minimum(eigvals(Hermitian(real_twirl))) >= -8e-12
        end

        permutation = permutation_operator((2, 2, 2), (2, 3, 1); T=ComplexF64)
        @test twirl(permutation; kind=:werner, copies=3) ≈ permutation atol = 1e-13

        identity4 = Matrix{Float64}(I, 4, 4)
        swap4 = Matrix(swap_operator((2, 2)))
        omega = vec(Matrix{Float64}(I, 2, 2))
        for fixed_operator in (identity4, swap4, omega * transpose(omega))
            @test twirl(fixed_operator; kind=:real) ≈ fixed_operator atol = 1e-13
        end
    end

    @testset "exact, generic-precision, and sparse arithmetic" begin
        exact_operator = Rational{Int}[
            1 2 0 1
            3 5 1 0
            0 2 4 1
            1 0 3 6
        ]
        for kind in (:werner, :isotropic, :real, :pauli)
            projected = twirl(exact_operator; kind)
            @test eltype(projected) == Rational{BigInt}
            @test twirl(projected; kind) == projected
            @test tr(projected) == tr(exact_operator)
        end

        for value_type in (Float32, Float64, ComplexF32, ComplexF64)
            typed_operator = value_type.(exact_operator)
            for kind in (:werner, :isotropic, :real, :pauli)
                projected = twirl(typed_operator; kind)
                @test eltype(projected) == value_type
                tolerance = value_type <: Union{Float32,ComplexF32} ? 3.0f-5 : 1e-12
                @test norm(twirl(projected; kind) - projected) <= tolerance
            end
        end

        big_operator = BigFloat.(exact_operator)
        big_projected = twirl(big_operator; kind=:real)
        @test eltype(big_projected) == BigFloat
        @test norm(twirl(big_projected; kind=:real) - big_projected) <= big"1e-60"

        exact_complex = complex.(exact_operator, reverse(exact_operator; dims=2))
        complex_projected = twirl(exact_complex; kind=:pauli)
        @test eltype(complex_projected) == Complex{Rational{BigInt}}
        @test twirl(complex_projected; kind=:pauli) == complex_projected

        sparse_operator = sparse(exact_operator)
        sparse_copy = copy(sparse_operator)
        sparse_projected = twirl(sparse_operator; kind=:werner)
        @test issparse(sparse_projected)
        @test eltype(sparse_projected) == Rational{BigInt}
        @test sparse_operator == sparse_copy
        @test Matrix(sparse_projected) == twirl(exact_operator; kind=:werner)
        dense_projected = twirl(
            sparse_operator; kind=:werner, sparse_output=false, allow_densify=true
        )
        @test dense_projected isa Matrix{Rational{BigInt}}
        @test dense_projected == Matrix(sparse_projected)
    end

    @testset "one-copy and one-dimensional boundaries" begin
        operator = Rational{Int}[1 2 3; 4 5 6; 7 8 10]
        expected = Rational{BigInt}(tr(operator) / 3) * Matrix{Rational{BigInt}}(I, 3, 3)
        @test twirl(operator; kind=:werner, copies=1) == expected
        @test twirl(operator; kind=:real, copies=1) == expected

        scalar = reshape(Rational{Int}[7 // 3], 1, 1)
        expected_scalar = reshape(Rational{BigInt}[7 // 3], 1, 1)
        @test twirl(scalar; kind=:werner, copies=100) == expected_scalar
        @test twirl(scalar; kind=:real, copies=100) == expected_scalar
        @test twirl(scalar; kind=:isotropic) == expected_scalar
        @test twirl(scalar; kind=:pauli) == expected_scalar
    end

    @testset "compatibility defaults and keyword forwarding" begin
        operator = reshape(Float64.(1:16), 4, 4)
        @test TwirlCompat.Twirl(operator) == twirl(operator)
        @test issparse(TwirlCompat.Twirl(operator))
        @test TwirlCompat.Twirl(operator, "werner") == twirl(operator)
        @test TwirlCompat.Twirl(operator, "real", 2) == twirl(operator; kind=:real)
        @test TwirlCompat.Twirl(
            sparse(operator), "isotropic", 2; sparse_output=false, allow_densify=true
        ) == twirl(
            sparse(operator); kind=:isotropic, sparse_output=false, allow_densify=true
        )
        @test_throws ArgumentError TwirlCompat.Twirl(operator, "unknown", 2)
        @test_throws ArgumentError TwirlCompat.Twirl(operator, 1, 2)
        @test_throws ArgumentError TwirlCompat.Twirl(operator, "werner", 2.0)
    end

    @testset "strict validation and resource guards" begin
        operator = ones(4, 4)
        @test_throws DimensionMismatch twirl(ones(2, 3))
        @test_throws ArgumentError twirl(zeros(0, 0))
        @test_throws ArgumentError twirl([1.0 NaN; 0.0 1.0]; copies=1)
        @test_throws ArgumentError twirl(Any[1.0 0.0; 0.0 1.0]; copies=1)
        @test_throws ArgumentError twirl(Bool[1 0; 0 1]; copies=1)
        @test_throws ArgumentError twirl(operator; kind="werner")
        @test_throws ArgumentError twirl(operator; kind=:Werner)
        @test_throws ArgumentError twirl(operator; kind=:unknown)
        @test_throws ArgumentError twirl(operator; copies=true)
        @test_throws ArgumentError twirl(operator; copies=0)
        @test_throws ArgumentError twirl(operator; kind=:isotropic, copies=1)
        @test_throws ArgumentError twirl(ones(8, 8); kind=:isotropic, copies=3)
        @test_throws ArgumentError twirl(operator; kind=:pauli, copies=1)
        @test_throws DimensionMismatch twirl(ones(9, 9); kind=:pauli)
        @test_throws DimensionMismatch twirl(ones(8, 8); kind=:werner, copies=2)
        @test_throws ArgumentError twirl(
            _WP2TwirlZeroBasedMatrix(Matrix{Float64}(I, 2, 2)); copies=1
        )

        @test_throws ArgumentError twirl(
            ones(16, 16); kind=:werner, copies=4, max_basis_size=23
        )
        @test_throws ArgumentError twirl(
            ones(8, 8); kind=:real, copies=3, max_basis_size=14
        )
        @test_throws ArgumentError twirl(ones(16, 16); kind=:pauli, max_basis_size=15)
        @test_throws ArgumentError twirl(ones(4, 4); kind=:isotropic, max_basis_size=1)
        @test_throws ArgumentError twirl(
            ones(8, 8); kind=:werner, copies=3, max_nonzeros=47
        )
        @test_throws ArgumentError twirl(ones(4, 4); kind=:isotropic, max_nonzeros=7)
        @test_throws ArgumentError twirl(
            ones(8, 8); kind=:werner, copies=3, max_dense_entries=107
        )
        @test_throws ArgumentError twirl(sparse(ones(4, 4)); sparse_output=false)
        @test_throws ArgumentError twirl(
            sparse(ones(4, 4));
            sparse_output=false,
            allow_densify=true,
            max_dense_entries=15,
        )
        @test_throws ArgumentError twirl(ones(4, 4); kind=:isotropic, max_work=31)

        for keyword in (:max_basis_size, :max_nonzeros, :max_dense_entries, :max_work)
            for invalid in (0, -1, true, 1.5)
                keywords = (; keyword => invalid)
                @test_throws ArgumentError twirl(operator; keywords...)
            end
        end

        unrestricted = twirl(
            sparse(operator);
            kind=:real,
            max_basis_size=nothing,
            max_nonzeros=nothing,
            max_dense_entries=nothing,
            max_work=nothing,
        )
        @test unrestricted == twirl(sparse(operator); kind=:real)
    end
end
