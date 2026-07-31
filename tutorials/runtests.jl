using Test

include("subsystem_reductions.jl")
include("local_channel_noise.jl")
include("entanglement_certificates.jl")
include("separability_examples.jl")
include("symmetric_sappt_witnesses.jl")
include("qetlab_intro_schmidt.jl")
include("qetlab_intro_tiles.jl")

@testset "Executable tutorials" begin
    @testset "QETLAB introduction: Schmidt decomposition" begin
        output = IOBuffer()
        result = TutorialQETLABIntroSchmidt.run(; io=output)

        @test result.local_dimensions == (3, 3)
        @test result.state_dimension == 9
        @test result.term_count == 3
        @test result.numerical_schmidt_rank == 3
        @test result.coefficients_nonnegative
        @test result.coefficients_descending
        @test result.replay_exact
        @test result.normalization_error <= 1e-12
        @test result.coefficient_normalization_error <= 1e-12
        @test result.left_orthogonality_error <= 1e-12
        @test result.right_orthogonality_error <= 1e-12
        @test result.manual_reconstruction_error <= 1e-12
        @test result.tensor_sum_reconstruction_error <= 1e-12
        @test result.reconstruction_agreement_error <= 1e-12
    end

    @testset "QETLAB introduction: Tiles bound entanglement" begin
        output = IOBuffer()
        result = TutorialQETLABIntroTiles.run(; io=output)

        @test result.catalog_family === :tiles
        @test result.local_dimensions == (3, 3)
        @test result.product_vector_count == 5
        @test result.catalog_projector_error <= 1e-14
        @test result.upb_status === :upb
        @test result.upb_reason === :unextendible
        @test result.upb_certificate === :exact
        @test result.partitions_examined == 20
        @test result.complement_rank == 4
        @test result.exact_trace == 1
        @test result.exact_projector
        @test result.exact_density_operator
        @test result.exact_ppt
        @test result.exact_range_entanglement
        @test result.bound_entangled
        @test result.numeric_ppt_status === TutorialQETLABIntroTiles.CriterionUnknown
        @test result.entanglement_status === :entangled
        @test result.entanglement_certified
        @test result.entanglement_method === :realignment
        @test result.entanglement_certificate === :realignment_cross_norm_violation
        @test result.attempts == (:ppt => :unknown, :realignment => :entangled)
        @test result.realignment_margin > result.realignment_tolerance
    end

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

    @testset "symmetric SAPPT witnesses" begin
        output = IOBuffer()
        result = TutorialSymmetricSAPPTWitnesses.run(; io=output)

        @test result.p_sappt_min == 30 // 31
        @test result.same_spectrum_error <= 1e-14
        @test result.separable_decomposition_terms == 19
        @test isapprox(result.separable_decomposition_weight_sum, 1; atol=1e-14, rtol=0)
        @test result.separable_decomposition_error <= 1e-14
        @test result.p_sappt_min < result.p_demo < result.witness_detection_limit
        @test min(result.restricted_pt_minimum_1_4, result.restricted_pt_minimum_2_3) > 0
        @test result.witness_expectation_demo < 0
        @test isapprox(
            result.witness_detection_limit, 0.9686241592915386; atol=1e-14, rtol=0
        )
        @test isapprox(result.product_witness_minimum, 0.0027637875; atol=1e-13, rtol=0)
        @test result.npt_minimum < 0 && isapprox(
            result.npt_witness_expectation, result.npt_minimum; atol=1e-14, rtol=0
        )
        @test result.unmatched_minus_witness_expectation > 0 && isapprox(
            result.minus_phase_witness_expectation,
            result.witness_expectation_at_pmin;
            atol=1e-13,
            rtol=0,
        )
    end
end
