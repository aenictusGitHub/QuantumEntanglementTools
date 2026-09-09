using Test

include("subsystem_reductions.jl")
include("local_channel_noise.jl")
include("entanglement_certificates.jl")
include("separability_examples.jl")
include("dephasing_channel_threshold.jl")
include("depolarizing_threshold_entanglement.jl")
include("pauli_channel_twirling.jl")
include("entanglement_monotone_noise_sweep.jl")
include("multipartite_entanglement_motifs.jl")
include("multipartite_entanglement_monogamy.jl")
include("entangled_state_families.jl")
include("bound_entangled_states.jl")
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
        @test result.numeric_ppt_status === TutorialQETLABIntroTiles.CriterionUnknown &&
            result.entanglement_input === :full_rank_depolarized_neighbor &&
            result.depolarizing_weight == 1 // 1024 &&
            result.neighbor_minimum_eigenvalue > 1e-5 &&
            result.neighbor_ppt_status === TutorialQETLABIntroTiles.CriterionSatisfied
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

    @testset "dephasing channel threshold" begin
        output = IOBuffer()
        result = TutorialDephasingNoiseThreshold.run(; io=output)

        @test result.is_cptp
        @test result.is_unital
        @test result.complete_positive
        @test result.fully_dephased_status === :separable ||
            result.fully_dephased_status === :unknown
        @test result.fully_dephased_ppt === TutorialDephasingNoiseThreshold.CriterionSatisfied ||
            result.fully_dephased_ppt === TutorialDephasingNoiseThreshold.CriterionUnknown
        @test result.partial_dephased_status === :entangled
        @test result.partial_dephased_ppt ===
            TutorialDephasingNoiseThreshold.CriterionEntanglementDetected
        @test result.partial_dephased_certified
        @test isapprox(result.partial_offdiagonal, 0.25; atol=1e-12, rtol=0)
        @test result.partial_negativity > 0
    end

    @testset "depolarizing threshold entanglement" begin
        output = IOBuffer()
        result = TutorialDepolarizingThresholdEntanglement.run(; io=output)

        @test result.is_cptp_0
        @test result.is_cptp_2_3
        @test result.separable_ppt ===
            TutorialDepolarizingThresholdEntanglement.CriterionSatisfied
        @test result.entangled_ppt ===
            TutorialDepolarizingThresholdEntanglement.CriterionEntanglementDetected
        @test result.separable_status === :separable
        @test result.entangled_status === :entangled
        @test isapprox(result.separable_concurrence, 0; atol=1e-12, rtol=0)
        @test result.entangled_concurrence > 0.3
        @test result.separable_overlap <= 0.25
        @test result.entangled_overlap > 0.66
    end

    @testset "pauli channel twirling" begin
        output = IOBuffer()
        result = TutorialPauliChannelTwirl.run(; io=output)

        @test result.identity_status === :entangled
        @test result.twirled_status === :entangled
        @test result.fully_twirled_status === :separable
        @test result.identity_is_cptp
        @test result.twirled_is_cptp
        @test result.fully_twirled_is_cptp
        @test result.twirled_ppt === TutorialPauliChannelTwirl.CriterionEntanglementDetected
        @test result.fully_twirled_ppt === TutorialPauliChannelTwirl.CriterionSatisfied
        @test isapprox(result.identity_concurrence, 1; atol=1e-12, rtol=0)
        @test result.twirled_concurrence >= 0.4 - 1e-12
        @test result.fully_twirled_concurrence <= 1e-12
        @test isapprox(result.fully_mixed_overlap, 0.25; atol=1e-12, rtol=0)
        @test result.twirled_overlap > 0.6
    end

    @testset "entanglement monotone noise sweep" begin
        output = IOBuffer()
        result = TutorialEntanglementMonotoneNoiseSweep.run(; io=output)

        @test result.local_dimensions == (2, 2)
        @test result.channel_points == [1, 1 / 2, 1 / 4, 0]
        @test result.concurrence_values[1] > result.concurrence_values[2]
        @test result.concurrence_values[2] > result.concurrence_values[3]
        @test result.concurrence_values[3] > result.concurrence_values[4]
        @test isapprox(result.concurrence_values[1], 1; atol=1e-12, rtol=0)
        @test result.concurrence_values[4] <= 1e-12
        @test result.negativity_values[1] > result.negativity_values[2]
        @test result.negativity_values[2] > result.negativity_values[3]
        @test result.negativity_values[3] > result.negativity_values[4]
        @test result.negativity_values[4] <= 1e-12
        @test result.log_neg_values[1] > result.log_neg_values[2]
        @test result.log_neg_values[2] > result.log_neg_values[3]
        @test result.log_neg_values[3] > result.log_neg_values[4]
        @test result.log_neg_values[4] <= 1e-12
        @test result.statuses[1:3] == [:entangled, :entangled, :entangled]
        @test result.statuses[4] === :separable || result.statuses[4] === :unknown
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

    @testset "multipartite entanglement motifs" begin
        output = IOBuffer()
        result = TutorialMultipartiteEntanglementMotifs.run(; io=output)

        @test result.local_dimensions == (2, 2, 2)
        @test isapprox(result.ghz_single_purity, 1 / 2; atol=1e-12, rtol=0)
        @test isapprox(result.w_single_purity, 5 / 9; atol=1e-12, rtol=0)
        @test result.ghz_pair_status === :separable || result.ghz_pair_status === :unknown
        @test result.w_pair_status === :entangled
        @test !result.ghz_pair_certified || result.ghz_pair_status === :entangled
        @test result.w_pair_certified
        @test result.ghz_pair_negativity <= 1e-12
        @test result.w_pair_negativity > 0
        @test result.ghz_single_entropy > 0
        @test result.w_single_entropy > 0
    end

    @testset "multipartite entanglement monogamy" begin
        output = IOBuffer()
        result = TutorialMultipartiteEntanglementMonogamy.run(; io=output)

        @test result.local_dimensions == (2, 2, 2)
        @test isapprox(result.ghz_pair_concurrence_ab, 0; atol=1e-12, rtol=0)
        @test isapprox(result.ghz_pair_concurrence_ac, 0; atol=1e-12, rtol=0)
        @test result.ghz_one_vs_rest_concurrence > 0.99
        @test isapprox(result.w_pair_concurrence_ab, result.w_pair_concurrence_ac; atol=1e-12, rtol=0)
        @test result.w_pair_concurrence_ab > 0.6
        @test result.ghz_tau > result.w_tau
        @test result.w_tau <= 1e-12
        @test result.w_tau >= -1e-12
        @test isapprox(result.ghz_entropy, 1; atol=1e-12, rtol=0)
        @test result.w_entropy > 0.9
        @test result.w_entropy < 1
    end

    @testset "entangled state families" begin
        output = IOBuffer()
        result = TutorialEntangledStateFamilies.run(; io=output)

        @test result.local_dimensions == (2, 2)
        @test result.iso_sep_ppt === TutorialEntangledStateFamilies.CriterionSatisfied
        @test result.iso_ent_ppt === TutorialEntangledStateFamilies.CriterionEntanglementDetected
        @test result.werner_sep_ppt === TutorialEntangledStateFamilies.CriterionSatisfied
        @test result.werner_ent_ppt === TutorialEntangledStateFamilies.CriterionEntanglementDetected
        @test result.iso_sep_status === :separable
        @test result.iso_ent_status === :entangled
        @test result.werner_sep_status === :separable
        @test result.werner_ent_status === :entangled
        @test result.iso_ent_negativity > 0.49
        @test result.werner_ent_negativity > 0.49
    end

    @testset "bound entangled states" begin
        output = IOBuffer()
        result = TutorialBoundEntangledStates.run(; io=output)

        @test result.family_parameter == 0.3
        @test result.local_dimensions == (3, 3)
        @test result.minimum_partial_transpose_eigenvalue >= -1e-12
        @test result.ppt_status === TutorialBoundEntangledStates.CriterionSatisfied ||
            result.ppt_status === TutorialBoundEntangledStates.CriterionUnknown
        @test result.entanglement_status === :entangled
        @test result.entanglement_certified
        @test result.entanglement_method === :realignment
        @test result.entanglement_certificate === :realignment_cross_norm_violation
        @test result.attempts[1][1] === :ppt
        @test result.attempts[1][2] === :satisfied || result.attempts[1][2] === :unknown
        @test result.attempts[2] === (:realignment => :entangled)
        @test result.realignment_value > result.realignment_threshold
        @test result.realignment_value > result.realignment_tolerance
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
