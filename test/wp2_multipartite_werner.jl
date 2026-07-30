using LinearAlgebra
using SparseArrays
using Test

const WernerCompat = QuantumEntanglementTools.MATLABCompat

function _reviewed_three_party_parameters(::Type{T}=Float64) where {T}
    return T[0.05, 0.04, 0.03, 0.03, 0.02]
end

function _direct_multipartite_werner(dim, alpha)
    permutations = [(1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]
    local_dims = (dim, dim, dim)
    total = dim^3
    coefficient_type = eltype(alpha)
    raw = spdiagm(0 => fill(one(coefficient_type), total))
    for (coefficient, permutation) in zip(alpha, permutations)
        raw -=
            coefficient * permutation_operator(local_dims, permutation; T=coefficient_type)
    end
    return raw / real(tr(raw))
end

@testset "WP2 multipartite Werner states" begin
    @testset "documented permutation sum" begin
        parameters = _reviewed_three_party_parameters()
        parameters_copy = copy(parameters)
        state = werner_state(3, parameters)
        @test issparse(state)
        @test size(state) == (27, 27)
        @test eltype(state) == Float64
        @test ishermitian(state)
        @test tr(state) ≈ 1
        @test minimum(eigvals(Hermitian(Matrix(state)))) >= -1e-14
        @test state == _direct_multipartite_werner(3, parameters)
        @test parameters == parameters_copy

        rational_parameters = Rational{Int}[1 // 20, 1 // 25, 1 // 30, 1 // 30, 1 // 50]
        rational_state = werner_state(2, rational_parameters)
        @test eltype(rational_state) == Rational{Int}
        @test tr(rational_state) == 1
        @test rational_state == _direct_multipartite_werner(2, rational_parameters)

        float32_state = werner_state(2, _reviewed_three_party_parameters(Float32))
        @test eltype(float32_state) == Float32
        @test tr(float32_state) ≈ 1.0f0

        dense_state = werner_state(2, parameters; sparse_output=false)
        @test dense_state isa Matrix
        @test dense_state == Matrix(werner_state(2, parameters))
    end

    @testset "complex inverse-pair coefficients" begin
        parameters = ComplexF64[0.04, 0.03, 0.02im, -0.02im, 0.01]
        state = werner_state(2, parameters)
        @test eltype(state) == ComplexF64
        @test ishermitian(state)
        @test tr(state) ≈ 1
        @test minimum(eigvals(Hermitian(Matrix(state)))) >= -1e-14
        broken = copy(parameters)
        broken[4] = 0.02im
        @test_throws ArgumentError werner_state(2, broken)
    end

    @testset "compatibility and bipartite vector" begin
        parameters = _reviewed_three_party_parameters()
        @test WernerCompat.WernerState(2, parameters) == werner_state(2, parameters)
        @test WernerCompat.WernerState(3, [0.2]) == werner_state(3, 0.2)
        @test WernerCompat.WernerState(3, 0.2) == werner_state(3, 0.2)
    end

    @testset "physicality, shape, and resource guards" begin
        parameters = _reviewed_three_party_parameters()
        @test_throws ArgumentError werner_state(2, Float64[])
        @test_throws ArgumentError werner_state(2, [0.1, 0.2])
        @test_throws ArgumentError werner_state(2, [0.1, 0.2, 0.3])
        @test_throws ArgumentError werner_state(1, parameters)
        @test_throws ArgumentError werner_state(2, Any[true, 0, 0, 0, 0])
        @test_throws ArgumentError werner_state(2, [0.1, 0.1, NaN, NaN, 0.1])
        @test_throws ArgumentError werner_state(2, parameters; max_permutations=5)
        @test_throws ArgumentError werner_state(2, parameters; max_nonzeros=47)
        @test_throws ArgumentError werner_state(2, parameters; max_work=47)
        @test_throws ArgumentError werner_state(
            2, parameters; sparse_output=false, max_dense_entries=63
        )
        @test_throws DomainError werner_state(3, [2.0, 0, 0, 0, 0])
        @test_throws DomainError werner_state(3, [3.0, 0, 0, 0, 0])
        @test werner_state(
            2,
            parameters;
            max_permutations=nothing,
            max_nonzeros=nothing,
            max_dense_entries=nothing,
            max_work=nothing,
        ) == werner_state(2, parameters)
    end
end
