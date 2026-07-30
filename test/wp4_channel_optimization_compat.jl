using LinearAlgebra
using Test

using QuantumEntanglementTools
const ChannelOptimizationCompat = QuantumEntanglementTools.MATLABCompat

@testset "WP4 channel optimization compatibility mappings" begin
    identity_kraus = [Matrix{ComplexF64}(I, 2, 2)]
    dephasing_kraus = [ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]]
    reset_zero = [ComplexF64[1 0; 0 0], ComplexF64[0 1; 0 0]]
    phase_flip = [0.8Matrix{ComplexF64}(I, 2, 2), 0.6ComplexF64[1 0; 0 -1]]

    diamond = ChannelOptimizationCompat.DiamondNorm(identity_kraus)
    @test diamond isa ChannelNormResult
    @test diamond.status === QuantumEntanglementTools.ChannelOptimizationAnalyticOptimal
    @test diamond.value == 1
    @test ChannelOptimizationCompat.DiamondNorm(identity_kraus, [2, 2]; structured=false) ==
        1

    completely_bounded = ChannelOptimizationCompat.CBNorm(reset_zero)
    @test completely_bounded.quantity === :completely_bounded
    @test completely_bounded.value == 2
    @test ChannelOptimizationCompat.CBNorm(reset_zero; structured=false) == 2

    discrimination = ChannelOptimizationCompat.ChannelDistinguishability(
        identity_kraus, identity_kraus, [0.8, 0.2]
    )
    @test discrimination isa ChannelDistinguishabilityResult
    @test discrimination.success_probability == 0.8
    @test discrimination.certificate_kind === :identical_channels
    @test ChannelOptimizationCompat.ChannelDistinguishability(
        identity_kraus, identity_kraus, [0.8, 0.2], [2, 2]; structured=false
    ) == 0.8
    @test_throws DomainError ChannelOptimizationCompat.ChannelDistinguishability(
        identity_kraus, dephasing_kraus, [1.1, -0.1]
    )

    fidelity = ChannelOptimizationCompat.MaximumOutputFidelity(identity_kraus, phase_flip)
    @test fidelity isa MaximumOutputFidelityResult
    @test fidelity.value == 1
    @test fidelity.certificate_kind === :common_basis_output
    @test ChannelOptimizationCompat.MaximumOutputFidelity(
        identity_kraus, phase_flip; structured=false
    ) == 1

    transpose_choi = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    unresolved = ChannelOptimizationCompat.DiamondNorm(transpose_choi, [2, 2])
    @test unresolved.status ===
        QuantumEntanglementTools.ChannelOptimizationBackendUnavailable
    @test unresolved.value === nothing
    @test_throws DomainError ChannelOptimizationCompat.DiamondNorm(
        transpose_choi, [2, 2]; structured=false
    )

    @test_throws DimensionMismatch ChannelOptimizationCompat.DiamondNorm(
        identity_kraus, [3, 3]
    )
    @test_throws ArgumentError ChannelOptimizationCompat.DiamondNorm(
        transpose_choi, [2, 2]; max_variables=1
    )
end
