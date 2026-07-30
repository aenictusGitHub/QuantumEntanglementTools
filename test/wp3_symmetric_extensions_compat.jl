using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const SymExtCompatQET = QuantumEntanglementTools
const SymExtCompat = SymExtCompatQET.MATLABCompat

@testset "WP3 symmetric-extension MATLAB compatibility" begin
    mixed = Matrix{Float64}(I, 4, 4) / 4
    bell_vector = [inv(sqrt(2.0)), 0, 0, inv(sqrt(2.0))]
    bell = bell_vector * transpose(bell_vector)

    @test SymExtCompat.SymmetricExtension(mixed, 2, [2, 2]) == 1
    @test SymExtCompat.SymmetricExtension(bell, 2, [2, 2]) == 0

    one_copy = SymExtCompat.SymmetricExtension(
        mixed, 1, [2, 2], 0, 0, nothing; return_witness=true
    )
    @test one_copy.ex == 1
    @test one_copy.wit == mixed
    ex, witness = one_copy
    @test ex == 1
    @test witness == mixed

    default_dimensions = SymExtCompat.SymmetricExtension(
        Matrix{Float64}(I, 6, 6) / 6, 1; structured=true
    )
    @test default_dimensions.dimensions == (2, 3)

    unavailable = SymExtCompat.SymmetricExtension(
        Matrix{Float64}(I, 6, 6) / 6, 2, [2, 3]; structured=true
    )
    @test unavailable.status === SymExtCompatQET.SymmetricExtensionBackendUnavailable
    @test unavailable.verdict === nothing
    @test_throws DomainError SymExtCompat.SymmetricExtension(
        Matrix{Float64}(I, 6, 6) / 6, 2, [2, 3]
    )

    zero_inner = SymExtCompat.SymmetricInnerExtension(
        zeros(4, 4), 2, [2, 2], 1, nothing; return_witness=true
    )
    @test zero_inner.ex == 1
    @test iszero(norm(zero_inner.wit))
    unavailable_inner = SymExtCompat.SymmetricInnerExtension(
        mixed, 2, [2, 2]; structured=true, allow_densify=true
    )
    @test unavailable_inner.status === SymExtCompatQET.SymmetricExtensionBackendUnavailable
    @test_throws DomainError SymExtCompat.SymmetricInnerExtension(
        mixed, 2, [2, 2]; allow_densify=true
    )

    @test_throws ArgumentError SymExtCompat.SymmetricExtension(mixed, 2, [2, 2], 2)
    @test_throws ArgumentError SymExtCompat.SymmetricExtension(mixed, 2, [2, 2], 0, 0, -1)
    @test_throws ArgumentError SymExtCompat.SymmetricExtension(
        mixed, 2, [2, 2]; return_witness=true, structured=true
    )
    @test_throws ArgumentError SymExtCompat.SymmetricInnerExtension(mixed, 1, [2, 2])

    first = SymExtCompat.RandomPPTState(Xoshiro(8201), [2, 3], [3, 4])
    second = SymExtCompat.RandomPPTState(Xoshiro(8201), [2, 3], [3, 4])
    @test first == second
    @test ishermitian(first)
    @test tr(first) ≈ 1 atol = 1e-12
    @test eigmin(Hermitian(first)) >= -1e-12
    @test eigmin(
        Hermitian(Matrix(SymExtCompatQET.partial_transpose(first, (2, 3); systems=(2,))))
    ) >= -1e-12

    structured_random = SymExtCompat.RandomPPTState(
        Xoshiro(8202), 2, 2, 1e-10, Inf; structured=true, real=true
    )
    @test structured_random.status === SymExtCompatQET.RandomPPTConstructed
    @test structured_random.verified
    @test structured_random.max_iterations == 1
    @test eltype(structured_random.state) === Float64

    limited = SymExtCompat.RandomPPTState(
        Xoshiro(8203), 2, nothing, 1e-12, 0; structured=true
    )
    @test limited.status === SymExtCompatQET.RandomPPTIterationLimit
    @test_throws DomainError SymExtCompat.RandomPPTState(
        Xoshiro(8203), 2, nothing, 1e-12, 0
    )
    @test_throws ArgumentError SymExtCompat.RandomPPTState(Xoshiro(1), 2, nothing, -1)
    @test_throws ArgumentError SymExtCompat.RandomPPTState(
        Xoshiro(1), 2, nothing, 1e-12, 1.5
    )
end
