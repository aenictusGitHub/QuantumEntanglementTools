using LinearAlgebra
using Random
using Test

using QuantumEntanglementTools
using QuantumEntanglementTools.MATLABCompat

@testset "WP4 optimization compatibility mappings" begin
    ket0 = ComplexF64[1, 0]
    ketplus = ComplexF64[1, 1] / sqrt(2)
    states = hcat(ket0, ketplus)

    structured = Distinguishability(states)
    @test structured isa StateDiscriminationResult
    @test structured.status === QuantumEntanglementTools.StateDiscriminationHelstromOptimal
    @test structured.success_probability ≈ (1 + inv(sqrt(2))) / 2
    @test structured.measurement !== nothing

    legacy = Distinguishability(states; structured=false)
    @test keys(legacy) == (:dist, :meas)
    @test legacy.dist == structured.success_probability
    @test legacy.meas == structured.measurement

    trine = hcat(
        ComplexF64[1, 1] / sqrt(2),
        ComplexF64[1, cis(2pi / 3)] / sqrt(2),
        ComplexF64[1, cis(4pi / 3)] / sqrt(2),
    )
    unavailable = Distinguishability(trine)
    @test unavailable.status ===
        QuantumEntanglementTools.StateDiscriminationBackendUnavailable
    @test_throws DomainError Distinguishability(trine; structured=false)
    @test_throws ArgumentError Distinguishability(2states)
    @test_throws DomainError Distinguishability(states, [1.1, -0.1])

    affine = complex_affine_variable(:X, 2, 3)
    atom = kpNorm(affine, 2, 3)
    @test atom isa TopKPNormEpigraph
    @test atom.k == 2
    @test atom.p == 3

    rho = hermitian_variable(:rho, 2)
    sigma = HermitianAffineMatrix(
        :sigma, Matrix{Float64}(I, 2, 2) / 2, Int[], Matrix{Float64}[], 4
    )
    psd = IsPSD(rho)
    @test psd isa HermitianAffineMatrix
    @test psd !== rho
    @test psd.constant == rho.constant
    @test_throws ArgumentError IsPSD(rho, 1.0e-8)

    fidelity_atom = MatsumotoFidelity(rho, sigma)
    @test fidelity_atom isa MatsumotoFidelityModel
    @test fidelity_atom.input_variable_count == 4
    @test fidelity_atom.dimension == 2

    polynomial = [1.0, 2.0, 3.0]
    sos = PolynomialSOS(MersenneTwister(1), polynomial, 2, 1, 0; inner_samples=3)
    @test sos isa PolynomialSOSResult
    @test sos.status === QuantumEntanglementTools.OptimizationBackendUnavailable
    @test sos.inner_bound !== nothing
    @test sos.samples_evaluated == 3

    @test_throws ArgumentError PolynomialSOS(
        MersenneTwister(1), polynomial, 2, 1, 0, "invalid"
    )
    @test_throws DomainError PolynomialSOS(
        MersenneTwister(1), polynomial, 2, 1, 0; structured=false
    )
end
