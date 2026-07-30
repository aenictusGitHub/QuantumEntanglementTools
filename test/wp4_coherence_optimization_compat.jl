using LinearAlgebra
using Test

using QuantumEntanglementTools
using QuantumEntanglementTools.MATLABCompat

@testset "coherence optimization compatibility mappings" begin
    diagonal = Diagonal(ComplexF64[0.5, 0.3, 0.2]) |> Matrix
    diagonal_result = IskIncoherent(diagonal, 1)
    @test diagonal_result isa CoherenceCriterionResult
    @test diagonal_result.verdict === true
    @test IskIncoherent(diagonal, 1; structured=false) == 1

    coherent = fill(inv(sqrt(3)) + 0im, 3)
    coherent_rho = coherent * coherent'
    @test IskIncoherent(coherent_rho, 1; structured=false) == 0
    missing = IskIncoherent(coherent_rho, 2; strategy=:sdp)
    @test missing.verdict === nothing
    @test_throws DomainError IskIncoherent(coherent_rho, 2; strategy=:sdp, structured=false)

    maximally_mixed = Matrix{ComplexF64}(I, 4, 4) / 4
    absolute = IsAbskIncoh(maximally_mixed, 1)
    @test absolute.verdict === true
    @test IsAbskIncoh(maximally_mixed, 1; structured=false) == 1
    undecided = Diagonal(ComplexF64[0.55, 0.18, 0.12, 0.09, 0.06]) |> Matrix
    @test IsAbskIncoh(undecided, 3).verdict === nothing
    @test_throws DomainError IsAbskIncoh(undecided, 3; structured=false)

    plus = fill(0.5 + 0im, 4)
    robustness = RobustnessCoherence(reshape(plus, 1, :))
    @test robustness isa CoherenceOptimizationResult
    @test robustness.value ≈ 3
    @test RobustnessCoherence(plus; structured=false) ≈ 3

    distance = TraceDistanceCoherence(plus; structured=false)
    @test keys(distance) == (:tdc, :diagonal)
    @test distance.tdc ≈ 1.5
    @test distance.diagonal ≈ fill(0.25, 4)

    generalized = GenRobustnesskCoherence(plus, 2; structured=false)
    @test keys(generalized) == (:robk, :sig)
    @test generalized.robk ≈ 1
    @test generalized.sig !== nothing
    @test real(tr(generalized.sig)) ≈ 1

    free = GenRobustnesskCoherence(diagonal, 3; structured=false)
    @test free.robk == 0
    @test free.sig === nothing

    mixed = ComplexF64[
        0.40 0.05 0.02im
        0.05 0.35 0.03
        -0.02im 0.03 0.25
    ]
    @test_throws DomainError RobustnessCoherence(mixed; strategy=:sdp, structured=false)
    @test_throws DomainError TraceDistanceCoherence(mixed; strategy=:sdp, structured=false)
    @test_throws DomainError GenRobustnesskCoherence(
        mixed, 2; strategy=:sdp, structured=false
    )
end
