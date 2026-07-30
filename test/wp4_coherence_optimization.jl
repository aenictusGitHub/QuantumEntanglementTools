using LinearAlgebra
using SparseArrays
using Test
using QuantumEntanglementTools

const QETCoherenceOptimization = QuantumEntanglementTools
if !isdefined(QETCoherenceOptimization, :is_k_incoherent)
    Base.include(
        QETCoherenceOptimization,
        joinpath(@__DIR__, "..", "src", "coherence", "coherence_optimization.jl"),
    )
end

function _coherence_reconstruct(decomposition, dimension)
    blocks = if hasproperty(decomposition, :blocks)
        decomposition.blocks
    else
        decomposition.unnormalized_blocks
    end
    matrix = zeros(ComplexF64, dimension, dimension)
    for (support, block) in zip(decomposition.supports, blocks)
        indices = collect(support)
        matrix[indices, indices] .+= block
    end
    return matrix
end

@testset "coherence optimization convergence slice" begin
    @testset "strict density and parameter validation" begin
        state = ComplexF64[1, im] / sqrt(2)
        rho = state * adjoint(state)
        @test_throws DimensionMismatch QETCoherenceOptimization.is_k_incoherent(
            ones(2, 3), 1
        )
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(
            ComplexF64[1 0.1; 0 0], 1
        )
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(2rho, 1)
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(
            ComplexF64[1.1 0; 0 -0.1], 1
        )
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(rho, true)
        @test_throws DomainError QETCoherenceOptimization.is_k_incoherent(rho, 0)
        @test_throws DomainError QETCoherenceOptimization.is_k_incoherent(rho, 3)
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(
            rho, 1; max_subsets=big(typemax(Int)) + 1
        )
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(
            rho, 1; max_band_search_nodes=0
        )
        @test_throws ArgumentError QETCoherenceOptimization.robustness_coherence(0.9state)
        @test_throws ArgumentError QETCoherenceOptimization.trace_distance_coherence(
            state; strategy=:invented
        )
        @test_throws ArgumentError QETCoherenceOptimization.is_k_incoherent(sparse(rho), 1)
        sparse_result = QETCoherenceOptimization.is_k_incoherent(
            sparse(rho), 2; allow_densify=true
        )
        @test sparse_result.verdict === true
        @test sparse_result.diagnostics.densified

        original = copy(rho)
        QETCoherenceOptimization.robustness_coherence(rho)
        @test rho == original
    end

    @testset "k-incoherence theorem branches and corrected band search" begin
        diagonal = Diagonal(ComplexF64[0.4, 0.35, 0.25])
        diagonal_result = QETCoherenceOptimization.is_k_incoherent(Matrix(diagonal), 1)
        @test diagonal_result.verdict === true
        @test diagonal_result.method === :diagonal
        @test diagonal_result.exact

        coherent = fill(inv(sqrt(3)) + 0im, 3)
        coherent_rho = coherent * adjoint(coherent)
        level_one = QETCoherenceOptimization.is_k_incoherent(coherent_rho, 1)
        @test level_one.verdict === false
        @test level_one.method === :offdiagonal_entry
        @test QETCoherenceOptimization.is_k_incoherent(coherent_rho, 3).verdict === true

        comparison_state = ComplexF64[
            0.40 0.10 0.00
            0.10 0.35 0.05
            0.00 0.05 0.25
        ]
        comparison_result = QETCoherenceOptimization.is_k_incoherent(comparison_state, 2)
        @test comparison_result.verdict === true
        @test comparison_result.method === :comparison_matrix
        comparison_violation = QETCoherenceOptimization.is_k_incoherent(coherent_rho, 2)
        @test comparison_violation.verdict === false
        @test comparison_violation.certificate_kind ===
            :comparison_matrix_2_incoherence_characterization

        band_vector = ComplexF64[1, 1, 1, 0] / sqrt(3)
        band_state =
            0.9 .* (band_vector * adjoint(band_vector)) +
            0.1 .* Diagonal(ComplexF64[0, 0, 0, 1])
        band_result = QETCoherenceOptimization.is_k_incoherent(band_state, 3)
        @test band_result.verdict === true
        @test band_result.method === :band_ordering
        @test band_result.band_ordering == collect(1:4)

        plus = fill(0.5 + 0im, 4)
        dephasing_state =
            0.5 .* Matrix{ComplexF64}(I, 4, 4) / 4 + 0.5 .* (plus * adjoint(plus))
        dephasing_result = QETCoherenceOptimization.is_k_incoherent(dephasing_state, 3)
        @test dephasing_result.verdict === true
        @test dephasing_result.method === :dephasing_criterion

        path = ComplexF64[
            1 1 0 0
            1 1 1 0
            0 1 1 1
            0 0 1 1
        ]
        found = QETCoherenceOptimization._cohopt_band_ordering(
            path, 2; max_search_nodes=10_000
        )
        @test found.status === :found
        @test QETCoherenceOptimization._cohopt_bandwidth(
            .!iszero.(path - Diagonal(diag(path))), found.ordering
        ) <= 2

        # Pinned has_band_k_ordering returns false for this labelled path
        # because the recursive reversal-pruning call swaps `candidate` and
        # `num_placed`. The valid ordering is [2, 1, 3] (or its reverse).
        pinned_defect_graph = ComplexF64[0 1 1; 1 0 0; 1 0 0]
        corrected = QETCoherenceOptimization._cohopt_band_ordering(
            pinned_defect_graph, 2; max_search_nodes=100
        )
        @test corrected.status === :found
        @test QETCoherenceOptimization._cohopt_bandwidth(
            .!iszero.(pinned_defect_graph), corrected.ordering
        ) == 2

        cycle = ones(ComplexF64, 5, 5) - Matrix{ComplexF64}(I, 5, 5)
        limited = QETCoherenceOptimization._cohopt_band_ordering(
            cycle, 2; max_search_nodes=1
        )
        @test limited.status === :limit
        @test !limited.exact

        # Exhaustive small-graph comparison: the independent recognizer must
        # agree with brute-force permutation enumeration.
        for mask in 0:(2 ^ 6 - 1), k in 1:4
            adjacency = falses(4, 4)
            bit = 0
            for column in 2:4, row in 1:(column - 1)
                edge = !iszero(mask & (1 << bit))
                adjacency[row, column] = edge
                adjacency[column, row] = edge
                bit += 1
            end
            matrix = Float64.(adjacency)
            exact = QETCoherenceOptimization._cohopt_band_ordering(
                matrix, k; max_search_nodes=100_000
            )
            orders = (
                [a, b, c, d] for a in 1:4 for b in 1:4 for c in 1:4 for
                d in 1:4 if length(unique((a, b, c, d))) == 4
            )
            brute = any(
                QETCoherenceOptimization._cohopt_bandwidth(adjacency, order) <= k for
                order in orders
            )
            @test exact.status === (brute ? :found : :not_found)
        end

        missing = QETCoherenceOptimization.is_k_incoherent(coherent_rho, 2; strategy=:sdp)
        @test missing.status === :backend_unavailable
        @test missing.verdict === nothing
        @test missing.problem isa QETCoherenceOptimization.SemidefiniteProgram
        limited_model = QETCoherenceOptimization.is_k_incoherent(
            coherent_rho,
            2;
            strategy=:sdp,
            limits=QETCoherenceOptimization.OptimizationLimits(max_variables=1),
        )
        @test limited_model.status === :resource_limit
        @test limited_model.problem === nothing
        psd_block_limit = QETCoherenceOptimization.is_k_incoherent(
            coherent_rho,
            2;
            strategy=:sdp,
            limits=QETCoherenceOptimization.OptimizationLimits(max_psd_blocks=1),
        )
        @test psd_block_limit.status === :resource_limit
        @test occursin("max_psd_blocks", psd_block_limit.message)
    end

    @testset "absolute k-incoherence theorem and SDP branches" begin
        maximally_mixed = Matrix{ComplexF64}(I, 4, 4) / 4
        @test QETCoherenceOptimization.is_absolutely_k_incoherent(maximally_mixed, 1).verdict ===
            true
        @test QETCoherenceOptimization.is_absolutely_k_incoherent(maximally_mixed, 4).verdict ===
            true

        pure = ComplexF64[1, 0, 0, 0] * ComplexF64[1, 0, 0, 0]'
        rank_failure = QETCoherenceOptimization.is_absolutely_k_incoherent(pure, 2)
        @test rank_failure.verdict === false
        @test rank_failure.method === :rank

        tight_rank = Diagonal(ComplexF64[1 / 3, 1 / 3, 1 / 3, 0])
        tight_result = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Matrix(tight_rank), 2
        )
        @test tight_result.verdict === true
        @test tight_result.method === :equal_nonzero_spectrum

        maximum_bound = Diagonal(ComplexF64[0.30, 0.25, 0.25, 0.20])
        maximum_result = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Matrix(maximum_bound), 2
        )
        @test maximum_result.verdict === true
        @test maximum_result.method === :maximum_eigenvalue

        purity_failure = Diagonal(ComplexF64[0.7, 0.2, 0.1])
        purity_result = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Matrix(purity_failure), 2
        )
        @test purity_result.verdict === false
        @test purity_result.method === :purity

        near_mixed = Diagonal(
            ComplexF64[0.25 + 2eps(Float64), 0.25 - 2eps(Float64), 0.25, 0.25]
        )
        boundary = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Matrix(near_mixed), 1; atol=16eps(Float64), rtol=0
        )
        @test boundary.status === :boundary
        @test boundary.verdict === nothing

        undecided = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Diagonal(ComplexF64[0.55, 0.18, 0.12, 0.09, 0.06]) |> Matrix, 3
        )
        @test undecided.status === :unknown
        @test undecided.verdict === nothing

        theorem_eight = QETCoherenceOptimization.is_absolutely_k_incoherent(
            Diagonal(ComplexF64[0.4, 0.3, 0.2, 0.1]) |> Matrix, 3; strategy=:sdp
        )
        @test theorem_eight.status === :backend_unavailable
        @test theorem_eight.problem.metadata.formulation ===
            :johnston_moein_pereira_plosker_theorem_8
    end

    @testset "closed-form robustness and trace-distance certificates" begin
        maximally_coherent = fill(0.5 + 0im, 4)
        robustness = QETCoherenceOptimization.robustness_coherence(maximally_coherent)
        @test robustness.status === :analytic_exact
        @test robustness.value ≈ 3
        @test robustness.noise_state !== nothing
        @test real(tr(robustness.free_state)) ≈ 1
        @test minimum(eigvals(Hermitian(robustness.unnormalized_noise))) ≥ -1e-12
        @test real(tr(robustness.unnormalized_noise)) ≈ robustness.value
        reconstructed = _coherence_reconstruct(robustness.decomposition, 4)
        @test reconstructed ≈ robustness.free_state * (1 + robustness.value)

        qubit = ComplexF64[0.6 0.2im; -0.2im 0.4]
        qubit_robustness = QETCoherenceOptimization.robustness_coherence(qubit)
        @test qubit_robustness.value ≈ 0.4
        @test qubit_robustness.method === :qubit
        @test isdiag(qubit_robustness.free_state)

        diagonal = Diagonal(ComplexF64[0.5, 0.3, 0.2]) |> Matrix
        zero_robustness = QETCoherenceOptimization.robustness_coherence(diagonal)
        @test zero_robustness.value == 0
        @test zero_robustness.noise_state === nothing
        @test iszero(zero_robustness.unnormalized_noise)

        trace_result = QETCoherenceOptimization.trace_distance_coherence(maximally_coherent)
        @test trace_result.value ≈ 1.5
        @test trace_result.free_state ≈ Matrix{Float64}(I, 4, 4) / 4
        @test trace_result.diagnostics.branch_index == 4

        basis = ComplexF64[1, 0, 0, 0]
        basis_trace = QETCoherenceOptimization.trace_distance_coherence(basis)
        @test basis_trace.value == 0
        @test basis_trace.free_state == Diagonal([1.0, 0, 0, 0])

        qubit_trace = QETCoherenceOptimization.trace_distance_coherence(qubit)
        @test qubit_trace.value ≈ 0.4
        @test qubit_trace.free_state == Matrix(Diagonal([0.6, 0.4]))

        mixed = ComplexF64[
            0.40 0.05 0.02im
            0.05 0.35 0.03
            -0.02im 0.03 0.25
        ]
        robustness_missing = QETCoherenceOptimization.robustness_coherence(
            mixed; strategy=:sdp
        )
        trace_missing = QETCoherenceOptimization.trace_distance_coherence(
            mixed; strategy=:sdp
        )
        @test robustness_missing.status === :backend_unavailable
        @test trace_missing.status === :backend_unavailable
        @test robustness_missing.value === nothing
        @test trace_missing.value === nothing
        @test length(robustness_missing.problem.psd_constraints) == 1
        @test length(trace_missing.problem.psd_constraints) == 2

        trace_limit = QETCoherenceOptimization.trace_distance_coherence(
            mixed;
            strategy=:sdp,
            limits=QETCoherenceOptimization.OptimizationLimits(max_variables=1),
        )
        @test trace_limit.status === :resource_limit
        equality_limit = QETCoherenceOptimization.trace_distance_coherence(
            mixed;
            strategy=:sdp,
            limits=QETCoherenceOptimization.OptimizationLimits(max_equalities=1),
        )
        @test equality_limit.status === :resource_limit
        @test occursin("max_equalities", equality_limit.message)
    end

    @testset "generalized robustness and reconstructed free cone" begin
        maximally_coherent = fill(0.5 + 0im, 4)
        expected = (3.0, 1.0, 1 / 3, 0.0)
        for k in 1:4
            result = QETCoherenceOptimization.generalized_robustness_k_coherence(
                maximally_coherent, k
            )
            @test result.value ≈ expected[k]
            @test result.status === :analytic_exact
            @test real(tr(result.free_state)) ≈ 1
            @test minimum(eigvals(Hermitian(result.unnormalized_noise))) ≥ -1e-12
            if iszero(result.value)
                @test result.noise_state === nothing
            else
                @test real(tr(result.noise_state)) ≈ 1
                reconstructed = _coherence_reconstruct(result.decomposition, 4)
                @test reconstructed ≈ result.free_state * (1 + result.value)
            end
        end

        state = ComplexF64[0.13 + 0.19im, -0.47im, 0.31, -0.28 + 0.08im]
        state ./= norm(state)
        for k in 2:3
            generalized = QETCoherenceOptimization.generalized_robustness_k_coherence(
                state, k
            )
            theorem = QETCoherenceOptimization.pure_k_coherence_robustness(state, k)
            @test generalized.value ≈ theorem.value atol = 2e-14
        end

        free_mixed = ComplexF64[
            0.40 0.10 0.00
            0.10 0.35 0.05
            0.00 0.05 0.25
        ]
        zero = QETCoherenceOptimization.generalized_robustness_k_coherence(free_mixed, 2)
        @test zero.value == 0
        @test zero.noise_state === nothing
        @test zero.diagnostics.free_test_method === :comparison_matrix

        mixed = ComplexF64[
            0.40 0.05 0.02im
            0.05 0.35 0.03
            -0.02im 0.03 0.25
        ]
        missing = QETCoherenceOptimization.generalized_robustness_k_coherence(
            mixed, 2; strategy=:sdp
        )
        @test missing.status === :backend_unavailable
        @test missing.problem.metadata.formulation === :generalized_factor_width_robustness

        subset_limit = QETCoherenceOptimization.generalized_robustness_k_coherence(
            Matrix{ComplexF64}(I, 8, 8) / 8, 4; strategy=:sdp, max_subsets=10
        )
        @test subset_limit.status === :resource_limit
        @test occursin("max_subsets", subset_limit.message)
    end

    @testset "floating precision and result surfaces" begin
        state32 = ComplexF32[1, im, -1] / sqrt(3.0f0)
        robust32 = QETCoherenceOptimization.robustness_coherence(state32)
        trace32 = QETCoherenceOptimization.trace_distance_coherence(state32)
        generalized32 = QETCoherenceOptimization.generalized_robustness_k_coherence(
            state32, 2
        )
        @test robust32.value isa Float32
        @test trace32.value isa Float32
        @test generalized32.value isa Float32
        @test occursin("quantity=robustness_coherence", sprint(show, robust32))

        criterion = QETCoherenceOptimization.is_k_incoherent(state32 * adjoint(state32), 2)
        @test criterion.margin isa Union{Nothing,Float32}
        @test occursin("property=k_incoherent", sprint(show, criterion))

        setprecision(BigFloat, 128) do
            state_big = Complex{BigFloat}[1, im, -1] / sqrt(big"3")
            robust_big = QETCoherenceOptimization.robustness_coherence(
                state_big; atol=big"1e-35", rtol=big"1e-30"
            )
            trace_big = QETCoherenceOptimization.trace_distance_coherence(
                state_big; atol=big"1e-35", rtol=big"1e-30"
            )
            generalized_big = QETCoherenceOptimization.generalized_robustness_k_coherence(
                state_big, 2; atol=big"1e-35", rtol=big"1e-30"
            )
            @test robust_big.value isa BigFloat
            @test trace_big.value isa BigFloat
            @test generalized_big.value isa BigFloat
            @test generalized_big.diagnostics.certificate_minimum_eigenvalue isa
                Union{Nothing,BigFloat}
        end
    end
end
