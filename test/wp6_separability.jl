using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const SeparabilityQET = QuantumEntanglementTools

if !isdefined(SeparabilityQET, :StateDiscriminationResult)
    Base.include(
        SeparabilityQET,
        joinpath(@__DIR__, "..", "src", "optimization", "state_discrimination.jl"),
    )
end
if !isdefined(SeparabilityQET, :LocalDistinguishabilityResult)
    Base.include(
        SeparabilityQET,
        joinpath(@__DIR__, "..", "src", "entanglement", "separability_optimization.jl"),
    )
end

function wp6_bell_state(dimension=2)
    state = zeros(Float64, dimension^2, dimension^2)
    basis_indices = [1 + (index - 1) * (dimension + 1) for index in 1:dimension]
    state[basis_indices, basis_indices] .= inv(Float64(dimension))
    return state
end

function wp6_isotropic_state(weight)
    bell = wp6_bell_state(3)
    state = (1 - weight) * Matrix{Float64}(I, 9, 9) / 9 + weight * bell
    state[end, end] +=
        one(eltype(state)) - real(SeparabilityQET._sepopt_represented_trace(state))
    return state
end

struct WP6ReadCountingMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::M
    reads::Base.RefValue{Int}
end

Base.size(matrix::WP6ReadCountingMatrix) = size(matrix.data)
Base.IndexStyle(::Type{<:WP6ReadCountingMatrix}) = IndexCartesian()

function Base.getindex(matrix::WP6ReadCountingMatrix, index::CartesianIndex{2})
    matrix.reads[] += 1
    return matrix.data[index]
end

function Base.getindex(matrix::WP6ReadCountingMatrix, row::Int, column::Int)
    matrix.reads[] += 1
    return matrix.data[row, column]
end

@testset "WP6 separability and local discrimination core" begin
    @testset "separability strategy discovery" begin
        strategies = SeparabilityQET.available_separability_strategies()
        @test strategies == (
            :ppt,
            :low_rank_ppt,
            :realignment,
            :centered_realignment,
            :reduction,
            :qubit_qudit,
            :rank4_chow,
            :separable_ball,
            :rank_one_identity,
            :operator_schmidt_rank,
            :positive_maps,
            :filter_covariance,
            :symmetric_extension,
            :symmetric_inner_extension,
            :randomized_subtraction,
        )
        @test length(unique(strategies)) == length(strategies)
        @test SeparabilityQET._separability_strategies(:full) == strategies
        @test SeparabilityQET._separability_strategies(:qetlab_deterministic) ==
            filter(!=(:randomized_subtraction), strategies)

        ppt = SeparabilityQET.describe_strategy(:ppt)
        @test ppt.name === :ppt
        @test ppt.deterministic
        @test !ppt.rng_required
        @test ppt.cost === :spectral_cubic
        @test ppt.certificate_directions == (:entangled, :separable)
        @test !ppt.optional_dependency
        @test ppt.dependency === nothing

        randomized = SeparabilityQET.describe_strategy(:randomized_subtraction)
        @test !randomized.deterministic
        @test randomized.rng_required
        @test randomized.certificate_directions == (:separable,)
        @test occursin("Explicit-RNG", randomized.description)

        outer = SeparabilityQET.describe_strategy(:symmetric_extension)
        @test outer.optional_dependency
        @test outer.dependency === :optimization_backend
        @test outer.certificate_directions == (:entangled,)
        inner = SeparabilityQET.describe_strategy(:symmetric_inner_extension)
        @test inner.certificate_directions == (:separable,)

        conservative = SeparabilityQET.describe_strategy(:rank_one_identity)
        @test isempty(conservative.certificate_directions)
        @test occursin("never promoted", conservative.description)

        @test_throws ArgumentError SeparabilityQET.describe_strategy("ppt")
        @test_throws ArgumentError SeparabilityQET.describe_strategy(:full)
        @test_throws ArgumentError SeparabilityQET.describe_strategy(:not_a_strategy)
    end

    @testset "local-discrimination validation and theorem branch" begin
        rho = Matrix{Float64}(I, 4, 4) / 4
        trivial = SeparabilityQET.local_distinguishability((rho,), (2, 2))
        @test trivial.status === SeparabilityQET.LocalDistinguishabilityTrivialOptimal
        @test trivial.relaxation_value == 1
        @test trivial.separable_lower_bound == 1
        @test trivial.separable_upper_bound == 1
        @test trivial.certified
        @test trivial.certificate_kind === :single_guess_separable_povm
        @test length(trivial.measurement) == 1
        @test trivial.measurement[1] == Matrix{ComplexF64}(I, 4, 4)
        @test trivial.dual_solution.source === :analytic_single_guess
        @test trivial.dual_solution.completeness_operator == rho
        @test trivial.residuals.valid

        deterministic = SeparabilityQET.local_distinguishability(
            (rho, rho), (2, 2); priors=[1.0, 0.0]
        )
        @test deterministic.certified
        @test deterministic.relaxation_value == 1
        @test deterministic.measurement[1] == Matrix{ComplexF64}(I, 4, 4)
        @test iszero(deterministic.measurement[2])
        @test deterministic.dual_solution.completeness_operator == rho
        @test_throws Exception (deterministic.measurement[1][1, 1] = 0)
        @test_throws Exception (deterministic.priors[1] = 0)

        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            (rho, rho), (2, 2); priors=[1.0, 1.0e-10]
        )
        tolerance_normalized_column = reshape([sqrt(1 + 1.0e-10), 0.0, 0.0, 0.0], 4, 1)
        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            tolerance_normalized_column, (2, 2)
        )

        local_fields = ntuple(
            index -> getfield(deterministic, index), fieldcount(typeof(deterministic))
        )
        @test_throws MethodError SeparabilityQET.LocalDistinguishabilityResult(
            local_fields...
        )

        @test_throws ArgumentError SeparabilityQET.local_distinguishability((2rho,), (2, 2))
        @test_throws DomainError SeparabilityQET.local_distinguishability(
            (rho, rho), (2, 2); priors=[1.1, -0.1]
        )
        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            (rho, rho), (2, 2); priors=[0.4, 0.4]
        )
        @test_throws DimensionMismatch SeparabilityQET.local_distinguishability(
            (rho, rho), (2, 3)
        )
        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            sparse(rho), (2, 2)
        )
        sparse_opt_in = SeparabilityQET.local_distinguishability(
            (sparse(rho),), (2, 2); allow_densify=true
        )
        @test sparse_opt_in.status === SeparabilityQET.LocalDistinguishabilityTrivialOptimal
        @test sparse_opt_in.certified
        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            (sparse(rho),), (2, 2); allow_densify=true, max_dense_entries=15
        )
        @test_throws ArgumentError SeparabilityQET.local_distinguishability(
            hcat([1.0, 0, 0, 0], [0.0, 0, 0, 2]), (2, 2)
        )
    end

    @testset "solver-neutral local hierarchy" begin
        columns = hcat([1.0, 0, 0, 0], [0.0, 0, 0, 1])
        problem = SeparabilityQET.local_distinguishability_problem(
            columns, (2, 2); order=2, ppt=true, bosonic=true
        )
        @test problem.dimensions == (2, 2)
        @test problem.order == 2
        @test problem.ppt
        @test problem.bosonic
        @test problem.input_kind === :pure_columns
        @test problem.program.sense === :maximize
        @test problem.program.variable_count == 72
        @test length(problem.measurement_views) == 2
        @test length(problem.extension_views) == 2
        @test length(problem.variable_views) == 2
        @test length(problem.program.psd_constraints) == 6
        @test SeparabilityQET.primal_residual(
            problem.program, problem.program.known_feasible_point; allow_densify=true
        ) <= 2e-15

        unavailable = SeparabilityQET.local_distinguishability(columns, (2, 2); order=2)
        @test unavailable.status ===
            SeparabilityQET.LocalDistinguishabilityBackendUnavailable
        @test unavailable.relaxation_value === nothing
        @test unavailable.measurement === nothing
        @test unavailable.problem isa SeparabilityQET.LocalDistinguishabilityProblem
        @test unavailable.optimization_result.status ===
            SeparabilityQET.OptimizationBackendUnavailable
        @test unavailable.separable_lower_bound == 1 / 2
        @test unavailable.separable_upper_bound == 1

        limited = SeparabilityQET.local_distinguishability(
            columns, (2, 2); order=4, max_order=3
        )
        @test limited.status === SeparabilityQET.LocalDistinguishabilityResourceLimit
        @test limited.problem === nothing
        @test occursin("max_order", limited.message)

        model_limited = SeparabilityQET.local_distinguishability(
            columns, (2, 2); limits=SeparabilityQET.OptimizationLimits(max_variables=10)
        )
        @test model_limited.status === SeparabilityQET.LocalDistinguishabilityResourceLimit
        @test occursin("variables", model_limited.message)

        general = SeparabilityQET.local_distinguishability_problem(
            columns, (2, 2); order=2, ppt=false, bosonic=false
        )
        @test general.program.variable_count == 128
        @test length(general.program.psd_constraints) == 2
        @test SeparabilityQET.primal_residual(
            general.program, general.program.known_feasible_point; allow_densify=true
        ) <= 2e-15

        float32_problem = SeparabilityQET.local_distinguishability_problem(
            Float32.(columns), (2, 2)
        )
        @test eltype(float32_problem.priors) === Float32
        @test eltype(float32_problem.program.objective.coefficients) === Float32
        @test eltype(float32_problem.program.known_feasible_point) === Float32
        @test float32_problem.tolerance isa Float32
    end

    @testset "certificate-first separability conclusions" begin
        diagonal = Diagonal([0.4, 0.1, 0.2, 0.3])
        product_report = SeparabilityQET.is_separable(diagonal, (2, 2); strategies=(:ppt,))
        @test product_report.status === :separable
        @test product_report.certified
        @test product_report.certificate_kind === :computational_basis_product_decomposition
        @test product_report.method === :product_decomposition
        @test product_report.evidence.reconstruction == diagonal
        @test sum(component.weight for component in product_report.evidence.components) == 1
        @test product_report.attempts isa Tuple
        @test_throws MethodError empty!(product_report.attempts)
        @test_throws Exception (product_report.evidence.reconstruction[1, 1] = 0)
        @test product_report.evidence.reconstruction == diagonal

        @test_throws ArgumentError SeparabilityQET.EntanglementAttempt(
            :forged,
            :user,
            :separable,
            true,
            :forged_certificate,
            (reconstruction=Matrix(diagonal),),
            "forged",
        )
        @test_throws ArgumentError SeparabilityQET.EntanglementReport(
            :separable,
            true,
            :forged_certificate,
            :forged,
            :user,
            (reconstruction=Matrix(diagonal),),
            collect(product_report.attempts),
            "forged",
        )

        @test_throws ArgumentError SeparabilityQET.is_separable(
            Diagonal([0.4, 0.1, 0.2, 0.3000000001]), (2, 2)
        )
        @test_throws ArgumentError SeparabilityQET.is_separable(
            Diagonal([0.4, 0.1, 0.2, 0.2999999999]), (2, 2)
        )
        trace_error = try
            SeparabilityQET.is_separable(Diagonal([0.4, 0.1, 0.2, 0.3000000001]), (2, 2))
            nothing
        catch error
            error
        end
        @test trace_error isa ArgumentError
        trace_message = sprint(showerror, trace_error)
        @test occursin("trace(rho)=", trace_message)
        @test occursin("unit-trace residual", trace_message)
        @test occursin("validate_density_matrix", trace_message)

        sparse_product = SeparabilityQET.is_separable(
            sparse(diagonal), (2, 2); allow_densify=true
        )
        @test sparse_product.status === :separable
        @test sparse_product.certificate_kind === :computational_basis_product_decomposition
        @test_throws ArgumentError SeparabilityQET.is_separable(
            sparse(diagonal), (2, 2); allow_densify=true, max_dense_entries=15
        )

        bell = wp6_bell_state()
        entangled = SeparabilityQET.is_separable(bell, (2, 2); strategies=(:ppt,))
        @test entangled.status === :entangled
        @test entangled.certified
        @test entangled.certificate_kind === :negative_partial_transpose_witness
        @test entangled.method === :ppt
        @test entangled.evidence.status === SeparabilityQET.CriterionEntanglementDetected
        @test real(
            dot(
                entangled.evidence.witness,
                SeparabilityQET.partial_transpose(bell, (2, 2); systems=(2,)),
                entangled.evidence.witness,
            ),
        ) < 0

        x = [0.0 1.0; 1.0 0.0]
        inside_ball = Matrix{Float64}(I, 4, 4) / 4 + 0.01 * kron(x, x) / 4
        ball_report = SeparabilityQET.is_separable(
            inside_ball, (2, 2); strategies=(:separable_ball,)
        )
        @test ball_report.status === :separable
        @test ball_report.certificate_kind === :gurvits_barnum_separable_ball
        @test ball_report.evidence.status === :separable_certified

        ppt_unknown = SeparabilityQET.is_separable(
            wp6_isotropic_state(0.1), (3, 3); strategies=(:ppt,)
        )
        @test ppt_unknown.status === :unknown
        @test !ppt_unknown.certified
        @test length(ppt_unknown.attempts) == 1
        @test ppt_unknown.attempts[1].status === :unknown
        @test ppt_unknown.attempts[1].raw_result.status ===
            SeparabilityQET.CriterionSatisfied

        centered = SeparabilityQET.is_separable(
            bell, (2, 2); strategies=(:centered_realignment,)
        )
        @test centered.status === :entangled
        @test centered.certificate_kind === :centered_realignment_covariance_violation
        @test centered.evidence.margin > 0

        rank_one_lookalike = 0.0625 * Matrix{Float64}(I, 9, 9)
        rank_one_lookalike[1, 1] += 0.4375
        rank_one_lookalike[2, 4] = 1.0e-30
        rank_one_lookalike[4, 2] = 1.0e-30
        exact_centered = Rational{BigInt}.(
            rank_one_lookalike - 0.0625 * Matrix{Float64}(I, 9, 9)
        )
        @test !iszero(det(exact_centered[[1, 2, 4], [1, 2, 4]]))
        rank_one_report = SeparabilityQET.is_separable(
            rank_one_lookalike, (3, 3); strategies=(:rank_one_identity,)
        )
        @test rank_one_report.status === :unknown
        @test !rank_one_report.certified
        rank_one_attempt = only(rank_one_report.attempts)
        @test rank_one_attempt.method === :rank_one_identity
        @test !rank_one_attempt.raw_result.input_level_exact_proof
        @test occursin("input-level exact", rank_one_attempt.message)

        phi32 = vec(Matrix{Float32}(I, 3, 3)) / sqrt(Float32(3))
        float32_isotropic =
            Float32(0.9) * Matrix{Float32}(I, 9, 9) / Float32(9) +
            Float32(0.1) * (phi32 * adjoint(phi32))
        float32_isotropic[end, end] +=
            Float32(1) - real(SeparabilityQET._sepopt_represented_trace(float32_isotropic))
        positive_maps = SeparabilityQET.is_separable(
            float32_isotropic, (3, 3); strategies=(:positive_maps,)
        )
        @test positive_maps.status === :unknown
        @test length(positive_maps.attempts) == 1
        @test length(positive_maps.attempts[1].raw_result.attempts) == 19
        @test all(
            attempt -> eltype(attempt.mapped) === Float32,
            positive_maps.attempts[1].raw_result.attempts,
        )

        ppt_boundary = SeparabilityQET.is_separable(
            wp6_isotropic_state(0.25), (3, 3); strategies=(:rank4_chow,)
        )
        @test ppt_boundary.status === :unknown
        @test length(ppt_boundary.attempts) == 1
        @test ppt_boundary.attempts[1].method === :rank4_chow
        @test occursin("PPT prerequisite", ppt_boundary.attempts[1].message)

        outer_unavailable = SeparabilityQET.is_separable(
            wp6_isotropic_state(0.1), (3, 3); strategies=(:symmetric_extension,)
        )
        @test outer_unavailable.status === :unknown
        @test length(outer_unavailable.attempts) == 1
        @test outer_unavailable.attempts[1].backend === :optimization
        @test outer_unavailable.attempts[1].raw_result.status ===
            SeparabilityQET.SymmetricExtensionBackendUnavailable

        @test_throws ArgumentError SeparabilityQET.is_separable(
            wp6_isotropic_state(0.1), (3, 3); strategies=(:randomized_subtraction,)
        )
        rng = MersenneTwister(77)
        untouched = copy(rng)
        limited_random = SeparabilityQET.is_separable(
            rng,
            wp6_isotropic_state(0.1),
            (3, 3);
            strategies=(:randomized_subtraction,),
            max_work=0,
        )
        @test limited_random.status === :unknown
        @test limited_random.attempts[1].method === :randomized_subtraction
        @test limited_random.attempts[1].raw_result.work_used == 0
        @test rand(rng) == rand(untouched)

        global_rng_before = copy(Random.default_rng())
        SeparabilityQET.is_separable(
            MersenneTwister(770),
            wp6_isotropic_state(0.1),
            (3, 3);
            strategies=(:randomized_subtraction,),
            max_work=0,
        )
        @test copy(Random.default_rng()) == global_rng_before

        spectral_budget = SeparabilityQET.is_separable(
            MersenneTwister(78),
            float32_isotropic,
            (3, 3);
            strategies=(:randomized_subtraction,),
            max_subtractions=1,
            max_product_restarts=1,
            max_product_iterations=0,
            max_work=12,
        )
        @test spectral_budget.status === :unknown
        @test spectral_budget.attempts[1].raw_result.work_used == 12
        @test occursin("before validating", spectral_budget.attempts[1].raw_result.message)

        unused_rng = MersenneTwister(93)
        unused_copy = copy(unused_rng)
        SeparabilityQET.is_separable(
            unused_rng, wp6_isotropic_state(0.1), (3, 3); strategies=(:ppt,)
        )
        @test rand(unused_rng) == rand(unused_copy)

        @test_throws ArgumentError SeparabilityQET.is_separable(
            2 * Matrix{Float64}(I, 4, 4) / 4, (2, 2)
        )
        nonhermitian = ComplexF64[0.5 0.1im; 0 0.5]
        @test_throws ArgumentError SeparabilityQET.is_separable(nonhermitian, (1, 2))
        @test_throws ArgumentError SeparabilityQET.is_separable(
            sparse(wp6_isotropic_state(0.1)), (3, 3)
        )
        @test_throws ArgumentError SeparabilityQET.is_separable(
            wp6_isotropic_state(0.1), (3, 3); strategies=(:not_a_strategy,)
        )
    end

    @testset "UPB replacement cone and resource semantics" begin
        tiles = SeparabilityQET.upb(:tiles)
        unavailable = SeparabilityQET.upb_sep_distinguishable(tiles.local_factors)
        @test unavailable.status === :backend_unavailable
        @test unavailable.feasibility === :unknown
        @test unavailable.separably_distinguishable === nothing
        @test !unavailable.certified
        @test unavailable.certificate_kind === nothing
        @test unavailable.candidates_generated == 30
        @test unavailable.partitions_examined == 30
        @test unavailable.work_used == 8060
        @test length(unavailable.replacement_vectors) == 30
        @test all(
            replacement -> norm(replacement.global_vector) ≈ 1,
            unavailable.replacement_vectors,
        )
        @test all(
            replacement -> replacement.orthogonality_residual <= 8 * unavailable.tolerance,
            unavailable.replacement_vectors,
        )
        @test all(
            replacement ->
                norm(
                    replacement.projection -
                    replacement.global_vector * adjoint(replacement.global_vector),
                ) <= 1e-14,
            unavailable.replacement_vectors,
        )
        first_replacement = first(unavailable.replacement_vectors)
        @test_throws Exception (first_replacement.global_vector[1] = 0)
        @test_throws Exception (first_replacement.projection[1, 1] = 0)
        @test all(
            iszero,
            first_replacement.projection -
            first_replacement.global_vector * adjoint(first_replacement.global_vector),
        )

        replacement_fields = ntuple(
            index -> getfield(first_replacement, index),
            fieldcount(typeof(first_replacement)),
        )
        @test_throws MethodError SeparabilityQET.UPBReplacementVector(replacement_fields...)
        upb_result_fields = ntuple(
            index -> getfield(unavailable, index), fieldcount(typeof(unavailable))
        )
        @test_throws MethodError SeparabilityQET.UPBSeparableDiscriminationResult(
            upb_result_fields...
        )

        problem = SeparabilityQET.upb_sep_distinguishability_problem(tiles.local_factors)
        @test problem.dimensions == (3, 3)
        @test problem.state_count == 5
        @test problem.candidate_count == 30
        @test problem.program.sense === :feasibility
        @test problem.program.variable_count == 30
        @test length(problem.program.equalities) == 81
        @test length(problem.program.intervals) == 30
        @test isempty(problem.program.psd_constraints)
        @test length(problem.program.primal_views) == 1

        phases = Diagonal(ComplexF64[1, im, cis(0.3)])
        complex_factors = (
            phases * tiles.local_factors[1], adjoint(phases) * tiles.local_factors[2]
        )
        complex_result = SeparabilityQET.upb_sep_distinguishable(complex_factors)
        @test complex_result.status === :backend_unavailable
        @test complex_result.candidates_generated == 30
        @test all(
            replacement ->
                replacement.orthogonality_residual <= 8 * complex_result.tolerance,
            complex_result.replacement_vectors,
        )

        limited = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_partitions=0
        )
        @test limited.status === :resource_limit
        @test limited.partitions_examined == 0
        @test limited.candidates_generated == 0
        @test limited.optimization_result === nothing

        work_limited_before_root = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_work=0
        )
        @test work_limited_before_root.status === :resource_limit
        @test work_limited_before_root.partitions_examined == 0
        @test work_limited_before_root.candidates_generated == 0
        @test work_limited_before_root.work_used == 0

        work_limited_at_root = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_work=1
        )
        @test work_limited_at_root.status === :resource_limit
        @test work_limited_at_root.partitions_examined == 0
        @test work_limited_at_root.candidates_generated == 0
        @test work_limited_at_root.work_used == 1

        work_limited_before_first = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_work=139
        )
        @test work_limited_before_first.status === :resource_limit
        @test work_limited_before_first.partitions_examined == 0
        @test work_limited_before_first.candidates_generated == 0
        @test work_limited_before_first.work_used == 5

        work_limited_after_first = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_work=140
        )
        @test work_limited_after_first.status === :resource_limit
        @test work_limited_after_first.partitions_examined == 1
        @test work_limited_after_first.candidates_generated == 1
        @test work_limited_after_first.work_used == 140

        state_limited = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_states=4
        )
        @test state_limited.status === :resource_limit
        @test state_limited.feasibility === :unknown
        @test state_limited.upb_analysis.reason === :state_limit
        @test state_limited.partitions_examined == 0
        @test state_limited.candidates_generated == 0
        @test state_limited.work_used == 0
        @test occursin("max_states=4", state_limited.message)

        read_counts = ntuple(_ -> Ref(0), length(tiles.local_factors))
        counted_factors = ntuple(length(tiles.local_factors)) do party
            WP6ReadCountingMatrix(tiles.local_factors[party], read_counts[party])
        end
        counted_state_limited = SeparabilityQET.upb_sep_distinguishable(
            counted_factors; max_states=4
        )
        @test counted_state_limited.status === :resource_limit
        @test all(
            read_counts[party][] == length(tiles.local_factors[party]) for
            party in eachindex(read_counts)
        )
        @test all(
            counted_factors[party].data == tiles.local_factors[party] for
            party in eachindex(counted_factors)
        )

        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_states=0
        )
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_states=4.0
        )
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_states=4, max_work=-1
        )
        @test_throws DimensionMismatch SeparabilityQET.upb_sep_distinguishable(
            (tiles.local_factors[1], tiles.local_factors[2][:, 1:4]); max_states=4
        )
        nonfinite_factors = (copy(tiles.local_factors[1]), copy(tiles.local_factors[2]))
        nonfinite_factors[1][1, 1] = Inf
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            nonfinite_factors; max_states=4
        )

        candidate_limited = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_candidates=1
        )
        @test candidate_limited.status === :resource_limit
        @test candidate_limited.candidates_generated == 1

        model_limited = SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; limits=SeparabilityQET.OptimizationLimits(max_variables=1)
        )
        @test model_limited.status === :resource_limit
        @test model_limited.problem === nothing
        @test occursin("variables", model_limited.message)

        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            tiles.local_factors; max_dense_entries=1
        )
        overflowing_factors = ntuple(_ -> reshape([1.0, 0.0], 2, 1), 64)
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(
            overflowing_factors
        )

        nonminimal = SeparabilityQET.upb_sep_distinguishable(
            SeparabilityQET.upb(:minimum_4x4).local_factors
        )
        @test nonminimal.status === :unsupported_nonminimal
        @test nonminimal.feasibility === :unknown
        @test nonminimal.separably_distinguishable === nothing
        @test nonminimal.candidates_generated == 0
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishability_problem(
            SeparabilityQET.upb(:minimum_4x4).local_factors
        )

        bad = (Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 3, 3))
        invalid = SeparabilityQET.upb_sep_distinguishable(bad)
        @test invalid.status === :invalid_input
        @test invalid.separably_distinguishable === nothing

        float32_tiles = SeparabilityQET.upb(:tiles; real_type=Float32)
        float32_result = SeparabilityQET.upb_sep_distinguishable(
            float32_tiles.local_factors
        )
        @test float32_result.status === :backend_unavailable
        @test float32_result.tolerance isa Float32
        @test float32_result.candidates_generated == 30

        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable((
            BigFloat.(tiles.local_factors[1]), BigFloat.(tiles.local_factors[2])
        ))
        rational_tiles = ntuple(length(tiles.local_factors)) do party
            rationalize.(tiles.local_factors[party])
        end
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable(rational_tiles)
        @test_throws ArgumentError SeparabilityQET.upb_sep_distinguishable((
            sparse(tiles.local_factors[1]), tiles.local_factors[2]
        ))
    end
end
