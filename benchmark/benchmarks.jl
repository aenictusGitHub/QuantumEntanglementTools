#!/usr/bin/env julia

using BenchmarkTools
using Dates
using InteractiveUtils
using LinearAlgebra
using Random
using SparseArrays
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)

using QuantumEntanglementTools

function parse_options(args)
    quick = false
    save_results = true
    output_prefix = nothing
    for arg in args
        if arg == "--quick"
            quick = true
        elseif arg == "--no-save"
            save_results = false
        elseif startswith(arg, "--output=")
            output_prefix = split(arg, "="; limit=2)[2]
        elseif arg in ("-h", "--help")
            println(
                """
                Usage: julia --project=benchmark benchmark/benchmarks.jl [options]

                  --quick          short smoke benchmark
                  --no-save        print results without writing raw artifacts
                  --output=PREFIX  artifact prefix (default: benchmark/results/local/<timestamp>)
                """,
            )
            exit(0)
        else
            error("unknown option: $arg")
        end
    end
    return (; quick, save_results, output_prefix)
end

function benchmark_suite()
    rng = MersenneTwister(0x5145544c41424a55)
    suite = BenchmarkGroup()

    dense_dims = (2, 3, 4, 2)
    dense_dimension = prod(dense_dims)
    dense_state = randn(rng, ComplexF64, dense_dimension)
    dense_operator = randn(rng, ComplexF64, dense_dimension, dense_dimension)
    trace_plan = PartialTracePlan(dense_dims, (2, 4))
    transpose_plan = PartialTransposePlan(dense_dims, (2, 4))
    permutation_plan = SubsystemPermutationPlan(dense_dims, (4, 2, 1, 3))

    suite["partial_trace/dense_pure_plan_reuse"] = @benchmarkable partial_trace(
        $dense_state, $trace_plan
    )
    suite["partial_trace/dense_matrix_plan_reuse"] = @benchmarkable partial_trace(
        $dense_operator, $trace_plan
    )
    suite["partial_trace/plan_construction"] = @benchmarkable PartialTracePlan(
        $dense_dims, (2, 4)
    )
    suite["partial_transpose/dense_plan_reuse"] = @benchmarkable partial_transpose(
        $dense_operator, $transpose_plan
    )
    suite["permutation/dense_matrix_plan_reuse"] = @benchmarkable permute_subsystems(
        $dense_operator, $permutation_plan
    )

    sparse_dims = (4, 4, 4, 4)
    sparse_dimension = prod(sparse_dims)
    sparse_operator = sprand(rng, ComplexF64, sparse_dimension, sparse_dimension, 0.005)
    sparse_trace_plan = PartialTracePlan(sparse_dims, (2, 4))
    sparse_transpose_plan = PartialTransposePlan(sparse_dims, (2, 4))
    suite["partial_trace/sparse_matrix_plan_reuse"] = @benchmarkable partial_trace(
        $sparse_operator, $sparse_trace_plan
    )
    suite["partial_transpose/sparse_plan_reuse"] = @benchmarkable partial_transpose(
        $sparse_operator, $sparse_transpose_plan
    )

    bipartite_dims = (8, 8)
    bipartite_operator = randn(rng, ComplexF64, 64, 64)
    realignment_plan = RealignmentPlan(bipartite_dims)
    suite["realignment/dense_plan_reuse"] = @benchmarkable realign(
        $bipartite_operator, $realignment_plan
    )

    factor_a = randn(rng, ComplexF64, 8, 8)
    factor_b = randn(rng, ComplexF64, 6, 6)
    suite["tensor_product/dense_two_factor"] = @benchmarkable tensor_product(
        $factor_a, $factor_b
    )
    repeated_game = rand(rng, Float64, 2, 2, 3, 2)
    suite["nonlocal/parallel_repetition_2x2x3x2_copies2"] = @benchmarkable parallel_repetition(
        $repeated_game, 2
    )
    benchmark_xor_probabilities = fill(0.25, 2, 2)
    benchmark_xor_parity = [0 0; 0 1]
    benchmark_bell_scenario = BellScenario(2, 2, 2, 2)
    benchmark_chsh_fc = [
        0.0 0.0 0.0
        0.0 1.0 1.0
        0.0 1.0 -1.0
    ]
    benchmark_behavior_probabilities = zeros(Float64, 2, 2, 2, 2)
    benchmark_behavior_probabilities[1, 1, :, :] .= 1
    benchmark_behavior = FullProbabilityBehavior(
        benchmark_behavior_probabilities, benchmark_bell_scenario; atol=0, rtol=0
    )
    benchmark_game = NonlocalGame(
        benchmark_xor_probabilities,
        reshape(
            Float64[
                iseven((a - 1) + (b - 1)) == iszero(benchmark_xor_parity[x, y]) for
                a in 1:2, b in 1:2, x in 1:2, y in 1:2
            ],
            2,
            2,
            2,
            2,
        );
        atol=0,
        rtol=0,
    )
    benchmark_bcs_constraints = [
        reshape(Int[1, 0, 0, 1], 2, 2), reshape(Int[0, 1, 1, 0], 2, 2)
    ]
    suite["nonlocal/xor_classical_chsh"] = @benchmarkable xor_game_value(
        $benchmark_xor_probabilities,
        $benchmark_xor_parity;
        regime=:classical,
        atol=0,
        rtol=0,
    )
    suite["nonlocal/bell_classical_chsh"] = @benchmarkable bell_inequality_bound(
        $benchmark_chsh_fc, $benchmark_bell_scenario; notation=:fc, regime=:classical
    )
    suite["nonlocal/npa_problem_level1"] = @benchmarkable npa_problem(
        $benchmark_behavior; level=1
    )
    suite["nonlocal/game_lower_bound_no_backend"] = @benchmarkable nonlocal_game_lower_bound(
        MersenneTwister(0x4e4c47), 2, $benchmark_game
    )
    suite["nonlocal/bcs_lower_bound_no_backend"] = @benchmarkable bcs_game_lower_bound(
        MersenneTwister(0x424353), 2, $benchmark_bcs_constraints
    )
    suite["nonlocal/bcs_classical_value"] = @benchmarkable bcs_game_value(
        $benchmark_bcs_constraints; regime=:classical
    )
    suite["nonlocal/bell_qubit_rectangular_problem"] = @benchmarkable bell_inequality_qubit_bound(
        ones(1, 2), zeros(1), zeros(2), [-1.0, 1.0], [-1.0, 1.0]
    )
    suite["projector/symmetric_d4_p4_sparse"] = @benchmarkable symmetric_projector(4, 4)
    benchmark_psd_variable = hermitian_variable(:benchmark_psd, 4)
    suite["optimization/psd_constraint_affine_d4"] = @benchmarkable positive_semidefinite_constraint(
        $benchmark_psd_variable
    )
    benchmark_kp_variable = complex_affine_variable(:benchmark_kp, 3, 2)
    suite["optimization/top_k_p_epigraph_3x2_k2_p3"] = @benchmarkable top_k_p_norm_epigraph(
        $benchmark_kp_variable, 2, 3
    )
    suite["optimization/top_k_p_dual_epigraph_3x2_k2_p3"] = @benchmarkable top_k_p_norm_dual_epigraph(
        $benchmark_kp_variable, 2, 3
    )
    benchmark_sk_diagonal = Diagonal([4.0, 3.0, 2.0, 1.0])
    benchmark_block_violation = Diagonal([-1.0, 2.0, 2.0, 2.0])
    suite["optimization/sk_norm_exact_d2_k2"] = @benchmarkable sk_operator_norm(
        MersenneTwister(0x534b01), $benchmark_sk_diagonal; k=2, dims=(2, 2), strength=0
    )
    suite["optimization/sk_norm_bounds_d2_k1"] = @benchmarkable sk_operator_norm(
        MersenneTwister(0x534b02), $benchmark_sk_diagonal; k=1, dims=(2, 2), strength=0
    )
    suite["entanglement/block_positive_identity_d2"] = @benchmarkable is_block_positive(
        MersenneTwister(0x534b03), Matrix{Float64}(I, 4, 4); dims=(2, 2), strength=0
    )
    suite["entanglement/block_positive_violation_d2"] = @benchmarkable is_block_positive(
        MersenneTwister(0x534b04), $benchmark_block_violation; dims=(2, 2), strength=0
    )
    benchmark_separable_diagonal = Diagonal([0.4, 0.1, 0.2, 0.3])
    benchmark_local_states = hcat(ComplexF64[1, 0, 0, 0], ComplexF64[0, 0, 0, 1])
    benchmark_tiles = upb(:tiles).local_factors
    suite["entanglement/separability_diagonal_d2"] = @benchmarkable is_separable(
        $benchmark_separable_diagonal, (2, 2)
    )
    suite["optimization/local_discrimination_problem_d2_order2"] = @benchmarkable local_distinguishability_problem(
        $benchmark_local_states, (2, 2); order=2
    )
    suite["entanglement/upb_separable_discrimination_tiles"] = @benchmarkable upb_sep_distinguishable(
        $benchmark_tiles
    )

    suite["operators/fourier_d32"] = @benchmarkable fourier_matrix(32)
    suite["states/werner_d8_sparse"] = @benchmarkable werner_state(8, 0.2)
    werner_three_party_parameters = [0.05, 0.04, 0.03, 0.03, 0.02]
    suite["states/werner_three_party_d3_sparse"] = @benchmarkable werner_state(
        3, $werner_three_party_parameters
    )
    suite["states/entangled_subspace_8x8_r2_dim36"] = @benchmarkable entangled_subspace(
        36, (8, 8); r=2
    )
    cnot_gate = [
        1.0 0 0 0
        0 1 0 0
        0 0 0 1
        0 0 1 0
    ]
    suite["entanglement/is_entangling_gate_cnot"] = @benchmarkable is_entangling_gate(
        $cnot_gate, (2, 2)
    )
    absolute_ppt_spectrum9 = collect(9.0:-1:1)
    absolute_ppt_spectrum4 = [0.45, 0.35, 0.1, 0.1]
    suite["entanglement/abs_ppt_constraints_p3"] = @benchmarkable abs_ppt_constraints(
        $absolute_ppt_spectrum9; dims=(3, 3)
    )
    suite["entanglement/is_abs_ppt_p2_exhaustive"] = @benchmarkable is_abs_ppt(
        $absolute_ppt_spectrum4; dims=(2, 2)
    )
    symmetric_extension_mixed = Matrix{Float64}(I, 4, 4) / 4
    discrimination_ket0 = ComplexF64[1, 0]
    discrimination_ket_plus = ComplexF64[1, 1] / sqrt(2)
    discrimination_states = hcat(discrimination_ket0, discrimination_ket_plus)
    suite["optimization/state_discrimination_helstrom_qubit"] = @benchmarkable state_distinguishability(
        $discrimination_states
    )
    suite["entanglement/symmetric_extension_two_qubit_analytic"] = @benchmarkable symmetric_extension(
        $symmetric_extension_mixed; dims=(2, 2), order=2
    )
    suite["entanglement/symmetric_extension_problem_bosonic_ppt_k2"] = @benchmarkable symmetric_extension_problem(
        $symmetric_extension_mixed;
        dims=(2, 2),
        order=2,
        ppt=true,
        bosonic=true,
        allow_densify=true,
    )
    is_upb_tiles = (
        [
            1 1 0 0 1
            0 -1 0 1 1
            0 0 1 -1 1
        ],
        [
            1 0 0 1 1
            -1 0 1 0 1
            0 1 -1 0 1
        ],
    )
    suite["entanglement/is_upb_tiles_exact"] = @benchmarkable is_upb(
        $is_upb_tiles; normalization=:allow
    )
    suite["entanglement/minimum_upb_size_twelve_qubits"] = @benchmarkable minimum_upb_size(
        fill(2, 12)
    )
    suite["entanglement/upb_tiles"] = @benchmarkable upb(:tiles)
    suite["entanglement/upb_generalized_tiles_1_d8"] = @benchmarkable upb(
        :generalized_tiles_1, 8
    )
    suite["entanglement/upb_johnston_twelve_qubits"] = @benchmarkable upb(
        :johnston_2_power_4k, 12
    )
    alon_lovasz_seed = UInt64(0x5145544c55504231)
    chen_johnston_seed = UInt64(0x5145544c55504232)
    suite["entanglement/upb_alon_lovasz_3x4"] = @benchmarkable upb(
        MersenneTwister($alon_lovasz_seed), (3, 4)
    )
    suite["entanglement/upb_chen_johnston_6x6"] = @benchmarkable upb(
        MersenneTwister($chen_johnston_seed), (6, 6)
    )

    density_rng = MersenneTwister(0x5145544c41424431)
    unitary_rng = MersenneTwister(0x5145544c41424432)
    povm_rng = MersenneTwister(0x5145544c41424433)
    suite["random/density_d32_rank16"] = @benchmarkable random_density_matrix(
        $density_rng, 32; rank=16
    )
    suite["random/unitary_d32"] = @benchmarkable random_unitary($unitary_rng, 32)
    suite["random/povm_d8_outcomes4"] = @benchmarkable random_povm($povm_rng, 8, 4)
    suite["random/ppt_state_2x3_full"] = @benchmarkable random_ppt_state(
        MersenneTwister(0x51505054), (2, 3)
    )

    channel_input = randn(rng, ComplexF64, 16, 16)
    channel = dephasing_channel(16, 0.25)
    dephasing_correlation = fill(0.25, 16, 16)
    dephasing_correlation[diagind(dephasing_correlation)] .= 1.0
    dephasing_factor = Matrix(cholesky(Hermitian(dephasing_correlation)).L)
    channel_kraus = KrausRepresentation([
        Diagonal(dephasing_factor[:, index]) for index in axes(dephasing_factor, 2)
    ])
    induced_schatten_rng = MersenneTwister(0x5145544c53484154)
    induced_schatten_start = randn(rng, ComplexF64, 16, 16)
    twirl_three_copy = randn(rng, ComplexF64, 8, 8)
    twirl_sparse_pauli = sprand(rng, ComplexF64, 16, 16, 0.2)
    full_rank_channel = depolarizing_channel(8, 0.25)
    channel_optimization_identity = KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
    channel_optimization_reset = KrausRepresentation([
        ComplexF64[1 0; 0 0], ComplexF64[0 1; 0 0]
    ])
    channel_optimization_phase_flip = KrausRepresentation([
        0.8Matrix{ComplexF64}(I, 2, 2), 0.6ComplexF64[1 0; 0 -1]
    ])
    channel_optimization_transpose = ChoiRepresentation(
        ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1], 2, 2
    )
    suite["channels/apply_dephasing_d16"] = @benchmarkable apply_channel(
        $channel_input, $channel
    )
    suite["channels/kraus_to_choi_dephasing_d16"] = @benchmarkable choi_representation(
        $channel_kraus
    )
    suite["channels/kraus_to_superoperator_dephasing_d16"] = @benchmarkable superoperator_representation(
        $channel_kraus
    )
    suite["channels/induced_schatten_3to2_dephasing_d16"] = @benchmarkable induced_schatten_lower_bound(
        $induced_schatten_rng,
        $channel_kraus,
        3;
        q=2,
        initial_matrix=($induced_schatten_start),
        max_iterations=10,
    )
    suite["channels/random_superoperator_tp_4to4_rank16"] = @benchmarkable random_superoperator(
        Xoshiro(0x51535031), 4; kraus_rank=16
    )
    suite["channels/random_superoperator_bistochastic_d3_rank3"] = @benchmarkable random_superoperator(
        Xoshiro(0x51535032), 3; unital=true, kraus_rank=3
    )
    suite["channels/random_superoperator_proportional_2to3_rank2"] = @benchmarkable random_superoperator(
        Xoshiro(0x51535033),
        (2, 3);
        trace_preserving=true,
        proportional_unital=true,
        kraus_rank=2,
    )
    suite["channels/twirl_werner_p3_d2_dense"] = @benchmarkable twirl(
        $twirl_three_copy; kind=:werner, copies=3
    )
    suite["channels/twirl_pauli_d4_sparse"] = @benchmarkable twirl(
        $twirl_sparse_pauli; kind=:pauli
    )
    suite["channels/diamond_norm_identity_d2_analytic"] = @benchmarkable diamond_norm(
        $channel_optimization_identity
    )
    suite["channels/diamond_norm_transpose_d2_model"] = @benchmarkable diamond_norm_problem(
        $channel_optimization_transpose
    )
    suite["channels/cb_norm_reset_d2_analytic"] = @benchmarkable cb_norm(
        $channel_optimization_reset
    )
    suite["channels/distinguishability_identical_d2_analytic"] = @benchmarkable channel_distinguishability(
        $channel_optimization_identity, $channel_optimization_identity; priors=(0.8, 0.2)
    )
    suite["channels/maximum_output_fidelity_d2_analytic"] = @benchmarkable maximum_output_fidelity(
        $channel_optimization_identity, $channel_optimization_phase_flip
    )
    operator_sum_input = randn(rng, ComplexF64, 16, 12)
    operator_sum_left = [randn(rng, ComplexF64, 20, 16) for _ in 1:4]
    operator_sum_right = [randn(rng, ComplexF64, 10, 12) for _ in 1:4]
    operator_sum_map = OperatorSumRepresentation(operator_sum_left, operator_sum_right)
    suite["channels/operator_sum_direct_16x12_to_20x10_terms4"] = @benchmarkable apply_channel(
        $operator_sum_input, $operator_sum_map
    )
    suite["channels/operator_sum_to_choi_16x12_to_20x10_terms4"] = @benchmarkable choi_representation(
        $operator_sum_map
    )
    suite["channels/operator_sum_dual_16x12_to_20x10_terms4"] = @benchmarkable dual_channel(
        $operator_sum_map
    )
    partial_prefix = randn(rng, ComplexF64, 2, 3)
    partial_operator_sum_input = kron(partial_prefix, operator_sum_input)
    suite["channels/partial_operator_sum_2x16_by_3x12_terms4"] = @benchmarkable partial_map(
        $partial_operator_sum_input, $operator_sum_map, 2, (2, 16), (3, 12)
    )
    suite["channels/is_cp_depolarizing_d8"] = @benchmarkable is_completely_positive(
        $full_rank_channel; allow_densify=true
    )
    canonical_left = [randn(rng, ComplexF64, 5, 4) for _ in 1:4]
    canonical_right = [randn(rng, ComplexF64, 6, 3) for _ in 1:4]
    canonical_general_map = OperatorSumRepresentation(canonical_left, canonical_right)
    suite["channels/canonical_general_4x3_to_5x6_terms4"] = @benchmarkable canonical_map_decomposition(
        $canonical_general_map
    )
    complementary_left = [randn(rng, ComplexF64, 8, 6) for _ in 1:4]
    complementary_right = [randn(rng, ComplexF64, 8, 5) for _ in 1:4]
    complementary_map = OperatorSumRepresentation(complementary_left, complementary_right)
    suite["channels/complement_operator_sum_6x5_to_8x8_terms4"] = @benchmarkable complementary_channel(
        $complementary_map
    )

    tier_d_operator = randn(rng, ComplexF64, 64, 48)
    tier_d_factor = randn(rng, ComplexF64, 64, 64)
    tier_d_density = tier_d_factor * adjoint(tier_d_factor)
    tier_d_density ./= real(tr(tier_d_density))
    tier_d_other_factor = randn(rng, ComplexF64, 64, 64)
    tier_d_other_density = tier_d_other_factor * adjoint(tier_d_other_factor)
    tier_d_other_density ./= real(tr(tier_d_other_density))
    diagonal_probabilities = collect(1.0:512.0)
    diagonal_probabilities ./= sum(diagonal_probabilities)
    diagonal_density = Diagonal(diagonal_probabilities)
    reverse_diagonal_density = Diagonal(reverse(diagonal_probabilities))
    induced_norm_rng = MersenneTwister(0x5145544c4142494d)
    induced_norm_start = randn(rng, ComplexF64, size(tier_d_operator, 2))
    fidelity_factor = randn(rng, ComplexF64, 32, 32)
    fidelity_density = fidelity_factor * adjoint(fidelity_factor)
    fidelity_density ./= real(tr(fidelity_density))
    fidelity_other_factor = randn(rng, ComplexF64, 32, 32)
    fidelity_other_density = fidelity_other_factor * adjoint(fidelity_other_factor)
    fidelity_other_density ./= real(tr(fidelity_other_density))
    schmidt_vector = randn(rng, ComplexF64, 32 * 32)
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    concurrence_density =
        0.7 * (bell * adjoint(bell)) + 0.3 * Matrix{ComplexF64}(I, 4, 4) / 4

    suite["measures/trace_norm_64x48"] = @benchmarkable trace_norm($tier_d_operator)
    suite["measures/schatten_p3_64x48"] = @benchmarkable schatten_norm($tier_d_operator, 3)
    suite["measures/ky_fan_k16_64x48"] = @benchmarkable ky_fan_norm($tier_d_operator, 16)
    suite["measures/induced_matrix_norm_3to2_64x48"] = @benchmarkable induced_matrix_norm(
        $induced_norm_rng,
        $tier_d_operator,
        3;
        q=2,
        initial_vector=($induced_norm_start),
        max_iterations=25,
    )
    suite["measures/purity_density_d64"] = @benchmarkable purity($tier_d_density)
    suite["measures/purity_diagonal_d512"] = @benchmarkable purity($diagonal_density)
    suite["measures/entropy_density_d64"] = @benchmarkable von_neumann_entropy(
        $tier_d_density; base=2
    )
    suite["measures/fidelity_density_d32"] = @benchmarkable fidelity(
        $fidelity_density, $fidelity_other_density
    )
    suite["measures/fidelity_diagonal_d512"] = @benchmarkable fidelity(
        $diagonal_density, $reverse_diagonal_density
    )
    suite["measures/matsumoto_fidelity_density_d32"] = @benchmarkable matsumoto_fidelity(
        $fidelity_density, $fidelity_other_density
    )
    matsumoto_model_rho = Diagonal([0.7, 0.3])
    matsumoto_model_sigma = Diagonal([0.4, 0.6])
    suite["optimization/matsumoto_fidelity_problem_d2"] = @benchmarkable matsumoto_fidelity_problem(
        $matsumoto_model_rho, $matsumoto_model_sigma
    )
    suite["measures/trace_distance_density_d64"] = @benchmarkable trace_distance(
        $tier_d_density, $tier_d_other_density
    )
    suite["measures/trace_distance_diagonal_d512"] = @benchmarkable trace_distance(
        $diagonal_density, $reverse_diagonal_density
    )
    suite["measures/concurrence_mixed_two_qubit"] = @benchmarkable concurrence(
        $concurrence_density
    )
    suite["entanglement/negativity_8x8"] = @benchmarkable negativity(
        $tier_d_density, (8, 8)
    )
    suite["entanglement/logarithmic_negativity_8x8"] = @benchmarkable logarithmic_negativity(
        $tier_d_density, (8, 8); base=2
    )
    suite["entanglement/schmidt_decomposition_32x32"] = @benchmarkable schmidt_decomposition(
        $schmidt_vector, (32, 32)
    )
    suite["entanglement/schmidt_coefficients_32x32"] = @benchmarkable schmidt_coefficients(
        $schmidt_vector, (32, 32)
    )
    suite["entanglement/schmidt_rank_32x32"] = @benchmarkable schmidt_rank(
        $schmidt_vector, (32, 32)
    )
    suite["criteria/ppt_8x8"] = @benchmarkable ppt_criterion($tier_d_density, (8, 8))
    suite["criteria/realignment_8x8"] = @benchmarkable realignment_criterion(
        $tier_d_density, (8, 8)
    )
    suite["criteria/reduction_8x8"] = @benchmarkable reduction_criterion(
        $tier_d_density, (8, 8)
    )

    coherence_plus = fill(ComplexF64(inv(sqrt(4096))), 4096)
    coherence_basis = fourier_matrix(64)
    coherence_state = randn(rng, ComplexF64, 64)
    coherence_state ./= norm(coherence_state)
    suite["coherence/l1_pure_d4096"] = @benchmarkable l1_coherence($coherence_plus)
    suite["coherence/l1_diagonal_d512"] = @benchmarkable l1_coherence($diagonal_density)
    suite["coherence/relative_entropy_pure_d4096"] = @benchmarkable relative_entropy_coherence(
        $coherence_plus; base=2
    )
    suite["coherence/relative_entropy_diagonal_d512"] = @benchmarkable relative_entropy_coherence(
        $diagonal_density; base=2
    )
    suite["coherence/rank_basis_transform_d64"] = @benchmarkable coherence_rank(
        $coherence_state; basis=($coherence_basis)
    )
    suite["coherence/pure_k_robustness_d4096_k64"] = @benchmarkable pure_k_coherence_robustness(
        $coherence_plus, 64
    )
    coherence_optimization_plus = fill(ComplexF64(inv(sqrt(64))), 64)
    coherence_optimization_diagonal = Matrix(Diagonal(fill(ComplexF64(1 / 16), 16)))
    suite["coherence/is_k_incoherent_diagonal_d16_k1"] = @benchmarkable is_k_incoherent(
        $coherence_optimization_diagonal, 1
    )
    suite["coherence/is_absolutely_k_incoherent_mixed_d16_k1"] = @benchmarkable is_absolutely_k_incoherent(
        $coherence_optimization_diagonal, 1
    )
    suite["coherence/robustness_pure_d64"] = @benchmarkable robustness_coherence(
        $coherence_optimization_plus
    )
    suite["coherence/trace_distance_pure_d64"] = @benchmarkable trace_distance_coherence(
        $coherence_optimization_plus
    )
    suite["coherence/generalized_robustness_pure_d64_k8"] = @benchmarkable generalized_robustness_k_coherence(
        $coherence_optimization_plus, 8
    )

    product_operator = tensor_product(
        randn(rng, ComplexF64, 4, 4),
        randn(rng, ComplexF64, 8, 8),
        randn(rng, ComplexF64, 2, 2),
    )
    rank_one_rectangular_density = zeros(6, 6)
    rank_one_rectangular_density[[1, 5], [1, 5]] .= 0.5
    sinkhorn_density = 0.9fidelity_density + 0.1Matrix{ComplexF64}(I, 32, 32) / 32
    suite["product/operator_schmidt_decomposition_8x8"] = @benchmarkable operator_schmidt_decomposition(
        $tier_d_density, (8, 8)
    )
    suite["product/operator_schmidt_hermitian_factors_4x8"] = @benchmarkable operator_schmidt_decomposition(
        $fidelity_density, (4, 8); hermitian_factors=true
    )
    suite["product/operator_sinkhorn_4x8"] = @benchmarkable operator_sinkhorn(
        $sinkhorn_density, (4, 8)
    )
    suite["product/filter_normal_form_4x8"] = @benchmarkable filter_normal_form(
        $sinkhorn_density, (4, 8)
    )
    suite["product/operator_analysis_4x8x2"] = @benchmarkable is_product_operator(
        $product_operator, (4, 8, 2)
    )
    suite["product/formation_mixed_two_qubit"] = @benchmarkable entanglement_of_formation(
        $concurrence_density, (2, 2)
    )
    suite["product/formation_rank_one_2x3"] = @benchmarkable entanglement_of_formation(
        $rank_one_rectangular_density, (2, 3)
    )

    majorization_first = randn(rng, 32, 32)
    majorization_second = randn(rng, 32, 32)
    compound_input = randn(rng, 8, 8)
    additive_compound_input = randn(rng, 16, 16)
    commutant_generators = [randn(rng, ComplexF64, 8, 8), randn(rng, ComplexF64, 8, 8)]
    copositive_diagonal = Diagonal(collect(1.0:8.0))
    polynomial_quartic = HomogeneousPolynomial([1.0, 0, 2, 0, 3], 2, 4)
    polynomial_quadratic = HomogeneousPolynomial([1.0, 2, 3], 2, 2)
    polynomial_rng = MersenneTwister(0x5145544c504f4c59)
    benchmark_copositive = Rational{Int}[1 -1; -1 1]
    benchmark_cycle_five = [
        0 1 0 0 1
        1 0 1 0 0
        0 1 0 1 0
        0 0 1 0 1
        1 0 0 1 0
    ]
    suite["matrix_analysis/majorizes_singular_values_32x32"] = @benchmarkable majorizes(
        $majorization_first, $majorization_second
    )
    suite["matrix_analysis/compound_dense_8_order3"] = @benchmarkable compound_matrix(
        $compound_input, 3
    )
    suite["matrix_analysis/additive_compound_dense_16_order2"] = @benchmarkable additive_compound_matrix(
        $additive_compound_input, 2
    )
    suite["matrix_analysis/commutant_dense_two_generators_d8"] = @benchmarkable commutant(
        $commutant_generators
    )
    suite["polynomials/copositive_quartic_d8"] = @benchmarkable copositive_polynomial(
        $copositive_diagonal
    )
    suite["polynomials/symmetric_matrix_n2_degree4_level2"] = @benchmarkable polynomial_as_matrix(
        $polynomial_quartic; level=2
    )
    suite["polynomials/bounds_n2_degree2_level1"] = @benchmarkable polynomial_bounds(
        $polynomial_rng, $polynomial_quadratic; level=1, inner_samples=0, allow_densify=true
    )
    suite["polynomials/sos_problem_n2_degree2_level0"] = @benchmarkable polynomial_sos_problem(
        $polynomial_quadratic
    )
    suite["polynomials/copositivity_exact_psd_d2"] = @benchmarkable copositivity_criterion(
        MersenneTwister(0xc001), $benchmark_copositive
    )
    suite["polynomials/clique_cycle5_certified_bounds"] = @benchmarkable clique_number_bounds(
        MersenneTwister(0xc002), $benchmark_cycle_five
    )

    return suite
end

function command_output(args...)
    try
        stdout = IOBuffer()
        process = run(pipeline(ignorestatus(Cmd(collect(args))); stdout, stderr=devnull))
        success(process) || return "unavailable"
        return chomp(String(take!(stdout)))
    catch
        return "unavailable"
    end
end

function result_rows(results::BenchmarkGroup)
    rows = Vector{Dict{String,Any}}()
    for name in sort!(collect(keys(results)); by=string)
        estimate = minimum(results[name])
        push!(
            rows,
            Dict(
                "name" => string(name),
                "minimum_time_ns" => estimate.time,
                "minimum_gctime_ns" => estimate.gctime,
                "memory_bytes" => estimate.memory,
                "allocations" => estimate.allocs,
                "samples" => length(results[name].times),
            ),
        )
    end
    return rows
end

function print_results(rows)
    println("case\tminimum time (ns)\tmemory (bytes)\tallocations\tsamples")
    for row in rows
        println(
            row["name"],
            '\t',
            row["minimum_time_ns"],
            '\t',
            row["memory_bytes"],
            '\t',
            row["allocations"],
            '\t',
            row["samples"],
        )
    end
end

function metadata(rows, options)
    return Dict(
        "schema_version" => 1,
        "recorded_at_utc" => string(now(UTC)),
        "quick" => options.quick,
        "repository_commit" =>
            command_output("git", "-C", REPOSITORY_ROOT, "rev-parse", "HEAD"),
        "repository_status" =>
            command_output("git", "-C", REPOSITORY_ROOT, "status", "--short"),
        "julia_version" => string(VERSION),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "cpu_name" => Sys.CPU_NAME,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_configuration" => sprint(show, BLAS.get_config()),
        "benchmarktools_version" => string(Base.pkgversion(BenchmarkTools)),
        "rng" => "MersenneTwister",
        "rng_seed_hex" => "0x5145544c41424a55",
        "cases" => rows,
    )
end

function main(args)
    options = parse_options(args)
    seconds = options.quick ? 0.2 : 2.0
    samples = options.quick ? 20 : 10_000
    suite = benchmark_suite()
    results = run(suite; verbose=(!options.quick), seconds, samples, evals=1)
    rows = result_rows(results)
    print_results(rows)

    if options.save_results
        default_stamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS")
        prefix = something(
            options.output_prefix, joinpath(@__DIR__, "results", "local", default_stamp)
        )
        mkpath(dirname(prefix))
        raw_path = prefix * ".json"
        metadata_path = prefix * ".toml"
        BenchmarkTools.save(raw_path, results)
        open(metadata_path, "w") do io
            return TOML.print(io, metadata(rows, options); sorted=true)
        end
        println("raw results: ", relpath(raw_path, REPOSITORY_ROOT))
        println("metadata: ", relpath(metadata_path, REPOSITORY_ROOT))
    end
end

main(ARGS)
