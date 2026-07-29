using Test

include("subsystem_reductions.jl")
include("local_channel_noise.jl")
include("entanglement_certificates.jl")
include("separability_examples.jl")

@testset "Executable tutorials" begin
    @testset "subsystem reductions" begin
        output = IOBuffer()
        result = TutorialSubsystemReductions.run(; io=output)

        @test result.state_dimension == 8
        @test result.reduced_pair_size == (4, 4)
        @test isapprox(result.reduced_pair_trace, 1; atol=1e-12, rtol=0)
        @test isapprox(result.reduced_pair_purity, 1; atol=1e-12, rtol=0)
        @test isapprox(result.reduced_qubit_purity, 0.5; atol=1e-12, rtol=0)
        @test isapprox(
            result.minimum_partial_transpose_eigenvalue, -0.5; atol=1e-12, rtol=0
        )
    end

    @testset "local channel noise" begin
        output = IOBuffer()
        result = TutorialLocalChannelNoise.run(; io=output)

        @test result.completely_positive
        @test result.trace_preserving
        @test result.unital
        @test result.representation_error <= 1e-12
        @test result.initial_status === :entangled
        @test result.initial_certificate === :negative_partial_transpose_witness
        @test result.noisy_status === :separable
        @test result.noisy_certificate === :ppt_low_dimension_theorem
        @test isapprox(result.noisy_purity, 0.25; atol=1e-12, rtol=0)
    end

    @testset "entanglement certificates" begin
        output = IOBuffer()
        result = TutorialEntanglementCertificates.run(; io=output)

        @test result.bell_status === :entangled
        @test result.bell_certificate === :pure_state_schmidt_rank
        @test result.product_status === :separable
        @test result.product_certificate === :pure_product_decomposition
        @test result.horodecki_status === :entangled
        @test result.horodecki_certificate === :realignment_cross_norm_violation
        @test result.horodecki_attempts == (:ppt => :unknown, :realignment => :entangled)
        @test result.mixed_status === :unknown
        @test !result.mixed_certified
        @test result.mixed_attempts ==
            (:ppt => :unknown, :realignment => :unknown, :reduction => :unknown)
    end

    @testset "separability examples" begin
        output = IOBuffer()
        result = TutorialSeparabilityExamples.run(; io=output)

        @test result.pure_status === :separable
        @test result.pure_certified
        @test result.pure_certificate === :pure_product_decomposition
        @test result.mixed_status === :separable
        @test result.mixed_certificate === :ppt_low_dimension_theorem
        @test result.mixed_ball_status === :separable_certified
        @test result.higher_pipeline_status === :unknown
        @test !result.higher_pipeline_certified
        @test result.higher_attempts ==
            (:ppt => :unknown, :realignment => :unknown, :reduction => :unknown)
        @test result.higher_ball_status === :separable_certified
        @test result.product_ball_status === :outside_ball
    end
end
