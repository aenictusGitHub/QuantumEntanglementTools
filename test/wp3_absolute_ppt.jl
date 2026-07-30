using LinearAlgebra
using Random
using SparseArrays
using Test

const AbsPPTQET = QuantumEntanglementTools

if !isdefined(AbsPPTQET, :abs_ppt_constraints)
    Base.include(
        AbsPPTQET, joinpath(@__DIR__, "..", "src", "entanglement", "absolute_ppt.jl")
    )
end

struct _AbsPPTZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_AbsPPTZeroBasedVector) = size(vector.storage)
Base.axes(vector::_AbsPPTZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_AbsPPTZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_AbsPPTZeroBasedVector, index::Int) = vector.storage[index + 1]

function _abs_ppt_direct_criss_cross_valid(ordering)
    p = maximum(last, ordering.positive_pairs)
    ranks = zeros(Int, p, p)
    for (rank, (row, column)) in enumerate(ordering.positive_pairs)
        ranks[row, column] = rank
        ranks[column, row] = rank
    end
    for i in 1:p, j in 1:p, k in 1:p, l in 1:p, m in 1:p, n in 1:p
        if ranks[i, j] > ranks[k, l] &&
            ranks[l, n] > ranks[j, m] &&
            ranks[i, n] < ranks[k, m]
            return false
        end
    end
    return true
end

function _abs_ppt_random_unitary(rng, dimension)
    factorization = qr(randn(rng, ComplexF64, dimension, dimension))
    return Matrix(factorization.Q)
end

@testset "WP3 absolute-PPT finite LMI family" begin
    @testset "Hildebrand p=2 and p=3 formulas" begin
        spectrum4 = Rational{Int}[4, 3, 2, 1]
        family2 = AbsPPTQET.abs_ppt_constraints(
            spectrum4;
            dims=(2, 2),
            max_constraints=nothing,
            max_work=nothing,
            max_entries=nothing,
        )
        @test family2.status === AbsPPTQET.AbsPPTEnumerationExhaustive
        @test family2.exhaustive
        @test length(family2.constraints) == 1
        @test family2.known_criss_cross_count == 1
        @test family2.monotone_ordering_count == 1
        @test family2.constraints[1].matrix == Rational{Int}[2 -2; -2 6]
        @test family2.constraints[1].ordering.positive_pairs == [(1, 1), (1, 2), (2, 2)]
        @test family2.constraints[1].ordering.negative_pairs == [(1, 2)]

        spectrum9 = Rational{Int}[9, 8, 7, 6, 5, 4, 3, 2, 1]
        family3 = AbsPPTQET.abs_ppt_constraints(
            spectrum9;
            dims=(3, 3),
            max_constraints=nothing,
            max_work=nothing,
            max_entries=nothing,
        )
        @test family3.exhaustive
        @test length(family3.constraints) == 2
        @test family3.candidate_orderings_checked == 2
        expected = Set([
            Rational{Int}[2 -7 -4; -7 6 -2; -4 -2 12],
            Rational{Int}[2 -7 -5; -7 8 -2; -5 -2 12],
        ])
        @test Set(constraint.matrix for constraint in family3.constraints) == expected

        for constraint in family3.constraints
            ordering = constraint.ordering
            @test _abs_ppt_direct_criss_cross_valid(ordering)
            @test ordering.negative_pairs ==
                filter(pair -> first(pair) != last(pair), ordering.positive_pairs)
            ranks = Dict(
                pair => rank for (rank, pair) in enumerate(ordering.positive_pairs)
            )
            @test all(
                ranks[(row, column)] < ranks[(row, column + 1)] for row in 1:3 for
                column in row:2
            )
            @test all(
                ranks[(row, column)] < ranks[(row + 1, column)] for column in 2:3 for
                row in 1:(column - 1)
            )
            @test AbsPPTQET.abs_ppt_lmi_matrix(ordering, spectrum9) == constraint.matrix
        end
    end

    @testset "combinatorial counts and preflight metadata" begin
        expected_constraints = (0, 1, 2, 10, 114)
        expected_monotone = (1, 1, 2, 12, 286)
        for p in 1:5
            dimension = p^2
            family = AbsPPTQET.abs_ppt_constraints(
                collect(Float64, dimension:-1:1);
                dims=(p, p),
                max_constraints=nothing,
                max_work=nothing,
                max_entries=nothing,
            )
            @test family.exhaustive
            @test length(family.constraints) == expected_constraints[p]
            @test family.monotone_ordering_count == expected_monotone[p]
            @test family.known_criss_cross_count == (p == 1 ? 1 : expected_constraints[p])
            @test family.entries_stored == BigInt(expected_constraints[p] * p^2)
            @test family.work_used >= family.candidate_orderings_checked
        end

        # QETLAB's criss-cross family has 2612 p=6 orderings (four are known
        # redundant rather than realizable). This is intentionally the only
        # heavier combinatorial case in the focused suite.
        family6 = AbsPPTQET.abs_ppt_constraints(
            collect(Float32, 36:-1:1);
            dims=(6, 6),
            max_constraints=2612,
            max_work=500_000_000,
            max_entries=1_000_000,
        )
        @test family6.exhaustive
        @test length(family6.constraints) == 2612
        @test family6.candidate_orderings_checked == 33_592
        @test family6.known_criss_cross_count == 2612

        constraint_capped = AbsPPTQET.abs_ppt_constraints(
            collect(9.0:-1:1); dims=(3, 3), max_constraints=1
        )
        @test constraint_capped.status === AbsPPTQET.AbsPPTEnumerationConstraintLimit
        @test !constraint_capped.exhaustive
        @test length(constraint_capped.constraints) == 1

        work_capped = AbsPPTQET.abs_ppt_constraints(
            collect(9.0:-1:1); dims=(3, 3), max_work=10
        )
        @test work_capped.status === AbsPPTQET.AbsPPTEnumerationWorkLimit
        @test work_capped.work_used == 10
        @test isempty(work_capped.constraints)

        entry_capped = AbsPPTQET.abs_ppt_constraints(
            collect(9.0:-1:1); dims=(3, 3), max_entries=8
        )
        @test entry_capped.status === AbsPPTQET.AbsPPTEnumerationEntryLimit
        @test entry_capped.entries_stored == 0
        @test isempty(entry_capped.constraints)
    end

    @testset "type, affine, sparse, and ownership policy" begin
        rational_spectrum = Rational{Int}[4, 3, 2, 1]
        rational_family = AbsPPTQET.abs_ppt_constraints(
            @view(rational_spectrum[:]); dims=(2, 2)
        )
        @test eltype(rational_family.constraints[1].matrix) == Rational{Int}

        float32_family = AbsPPTQET.abs_ppt_constraints(
            Float32[0.4, 0.3, 0.2, 0.1]; dims=(2, 2)
        )
        @test eltype(float32_family.constraints[1].matrix) == Float32

        big_family = AbsPPTQET.abs_ppt_constraints(
            BigFloat[big"0.4", big"0.3", big"0.2", big"0.1"]; dims=(2, 2)
        )
        @test eltype(big_family.constraints[1].matrix) == BigFloat

        sparse_spectrum = sparsevec([1, 4], [0.75, 0.25], 4)
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            sparse_spectrum; dims=(2, 2)
        )
        sparse_family = AbsPPTQET.abs_ppt_constraints(
            sparse_spectrum; dims=(2, 2), allow_densify=true, sparse_output=true
        )
        @test issparse(sparse_family.constraints[1].matrix)

        sparse_diagonal = spdiagm(0 => [0.4, 0.3, 0.2, 0.1])
        diagonal_family = AbsPPTQET.abs_ppt_constraints(sparse_diagonal; dims=(2, 2))
        @test diagonal_family.constraints[1].matrix ==
            AbsPPTQET.abs_ppt_constraints([0.4, 0.3, 0.2, 0.1]; dims=(2, 2)).constraints[1].matrix
        sparse_nondiagonal = sparse([0.4 0.01 0 0; 0.01 0.3 0 0; 0 0 0.2 0; 0 0 0 0.1])
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            sparse_nondiagonal; dims=(2, 2)
        )
        @test AbsPPTQET.abs_ppt_constraints(
            sparse_nondiagonal; dims=(2, 2), allow_densify=true
        ).exhaustive

        variable_count = 4
        affine_spectrum = [
            AbsPPTQET.AffineScalar(
                0.0, [index == variable ? 1.0 : 0.0 for index in 1:variable_count]
            ) for variable in 1:variable_count
        ]
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            affine_spectrum; dims=(2, 2)
        )
        affine_family = AbsPPTQET.abs_ppt_constraints(
            affine_spectrum; dims=(2, 2), assume_ordered=true
        )
        @test affine_family.symbolic
        @test affine_family.input_kind === :affine_spectrum
        @test affine_family.constraints[1].matrix isa AbsPPTQET.HermitianAffineMatrix
        point = [0.4, 0.3, 0.2, 0.1]
        @test Matrix(
            AbsPPTQET.evaluate_affine(affine_family.constraints[1].matrix, point)
        ) ≈ rational_family.constraints[1].matrix ./ 10

        ordering = rational_family.constraints[1].ordering
        exposed_copy = ordering.positive_pairs.storage
        exposed_copy[1] = (2, 2)
        @test ordering.positive_pairs[1] == (1, 1)
        @test_throws Base.CanonicalIndexError setindex!(ordering.positive_pairs, (2, 2), 1)
        family_copy = rational_family.constraints.storage
        empty!(family_copy)
        @test length(rational_family.constraints) == 1
    end

    @testset "validation and dimension conventions" begin
        @test AbsPPTQET.abs_ppt_constraints(ones(6); dims=2, max_constraints=1).dimensions ==
            (2, 3)
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(ones(6))
        @test_throws DimensionMismatch AbsPPTQET.abs_ppt_constraints(ones(6); dims=(2, 2))
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(ones(4); dims=(2,))
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(ones(4); dims=true)
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(Float64[])
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            [1.0, 0.0, NaN, 0.0]; dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            ComplexF64[1, 0, 0, 0]; dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            Bool[1, 0, 0, 0]; dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            _AbsPPTZeroBasedVector([1.0, 0.0, 0.0, 0.0]); dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_lmi_matrix(
            AbsPPTQET.abs_ppt_constraints([0.4, 0.3, 0.2, 0.1]; dims=(2, 2)).constraints[1].ordering,
            [0.1, 0.2, 0.3, 0.4],
        )
        @test_throws DimensionMismatch AbsPPTQET.abs_ppt_constraints(
            ones(2, 3); dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            [0.5 0.1; 0.0 0.5]; dims=(1, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            Matrix{Number}(I, 4, 4); dims=(2, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            BigFloat[0.4 0.01; 0.01 0.6]; dims=(1, 2)
        )
        @test_throws ArgumentError AbsPPTQET.abs_ppt_constraints(
            [0.4 0.01 0 0; 0.01 0.3 0 0; 0 0 0.2 0; 0 0 0 0.1]; dims=(2, 2), max_entries=15
        )
    end

    @testset "certificate-aware IsAbsPPT branches" begin
        maximally_mixed = fill(0.25, 4)
        sufficient = AbsPPTQET.is_abs_ppt(maximally_mixed; dims=(2, 2))
        @test sufficient.status === AbsPPTQET.AbsolutePPTSufficientTestPassed
        @test sufficient.verdict === true
        @test sufficient.certificate_kind === :gurvits_barnum_separable_ball
        @test sufficient.family === nothing

        gershgorin_only = [
            0.22249259167616695,
            0.10984449916591696,
            0.10666491354947588,
            0.1049086338262903,
            0.10176586116805718,
            0.09440888233881005,
            0.0912571004163801,
            0.08723008803210949,
            0.08142742982679317,
        ]
        gershgorin = AbsPPTQET.is_abs_ppt(gershgorin_only; dims=(3, 3))
        @test gershgorin.status === AbsPPTQET.AbsolutePPTSufficientTestPassed
        @test gershgorin.certificate_kind === :gershgorin_hildebrand_lmis
        @test gershgorin.margin > 0

        exhaustive = AbsPPTQET.is_abs_ppt([0.45, 0.35, 0.1, 0.1]; dims=(2, 2))
        @test exhaustive.status === AbsPPTQET.AbsolutePPTExhaustiveCertified
        @test exhaustive.verdict === true
        @test exhaustive.certificate_kind === :exhaustive_hildebrand_lmis
        @test exhaustive.family.exhaustive

        certified_not = AbsPPTQET.is_abs_ppt([1.0, 0.0, 0.0, 0.0]; dims=(2, 2))
        @test certified_not.status === AbsPPTQET.AbsolutePPTCertifiedNot
        @test certified_not.verdict === false
        @test certified_not.violating_constraint.diagnostic.value < 0
        @test certified_not.ordering_certificate isa AbsPPTQET.AbsPPTOrderingCertificate
        @test certified_not.ordering_certificate.minimum_gap > 0
        @test certified_not.ordering_certificate.minimum_product_margin > 0
        @test certified_not.backend_result === nothing

        boundary = AbsPPTQET.is_abs_ppt([0.5, 1 / 6, 1 / 6, 1 / 6]; dims=(2, 2))
        @test boundary.status === AbsPPTQET.AbsolutePPTNumericalBoundary
        @test boundary.verdict === nothing
        @test boundary.family.first_boundary !== nothing

        capped = AbsPPTQET.is_abs_ppt(
            fill(1 / 9, 9); dims=(3, 3), use_sufficient_tests=false, max_constraints=1
        )
        @test capped.status === AbsPPTQET.AbsolutePPTCappedUnknown
        @test capped.verdict === nothing
        @test capped.family.status === AbsPPTQET.AbsPPTEnumerationConstraintLimit
        @test capped.checked_constraints == 1
        @test capped.planned_orderings == 2

        backend_missing = AbsPPTQET.is_abs_ppt(
            [1.0, 0.0, 0.0, 0.0]; dims=(2, 2), max_realization_iterations=0
        )
        @test backend_missing.status === AbsPPTQET.AbsolutePPTBackendUnavailable
        @test backend_missing.verdict === nothing
        @test backend_missing.backend_result.status ===
            AbsPPTQET.OptimizationBackendUnavailable
        @test backend_missing.violating_constraint !== nothing

        model_capped = AbsPPTQET.is_abs_ppt(
            [1.0, 0.0, 0.0, 0.0];
            dims=(2, 2),
            max_realization_iterations=0,
            optimization_limits=AbsPPTQET.OptimizationLimits(max_intervals=1),
        )
        @test model_capped.status === AbsPPTQET.AbsolutePPTCappedUnknown
        @test model_capped.verdict === nothing
        @test model_capped.backend_result === nothing
        @test occursin("pre-allocation limits", model_capped.message)

        @test_throws ArgumentError AbsPPTQET.is_abs_ppt(
            maximally_mixed; dims=(2, 2), max_constraints=0
        )
        @test_throws ArgumentError AbsPPTQET.is_abs_ppt(
            maximally_mixed; dims=(2, 2), max_work=0
        )
        @test_throws ArgumentError AbsPPTQET.is_abs_ppt(
            maximally_mixed; dims=(2, 2), max_entries=0
        )
        @test_throws ArgumentError AbsPPTQET.is_abs_ppt(
            maximally_mixed; dims=(2, 2), max_realization_iterations=-1
        )

        exact_boundary = AbsPPTQET.is_abs_ppt(
            Rational{Int}[1 // 2, 1 // 6, 1 // 6, 1 // 6]; dims=(2, 2)
        )
        @test exact_boundary.verdict === true
        @test exact_boundary.status === AbsPPTQET.AbsolutePPTSufficientTestPassed

        zero_operator = AbsPPTQET.is_abs_ppt(zeros(4); dims=(2, 2))
        @test zero_operator.status === AbsPPTQET.AbsolutePPTAnalyticCertified
        @test zero_operator.verdict === true
        @test zero_operator.certificate_kind === :zero_operator

        one_factor = AbsPPTQET.is_abs_ppt([0.7, 0.3]; dims=(1, 2))
        @test one_factor.status === AbsPPTQET.AbsolutePPTAnalyticCertified
        @test one_factor.verdict === true
        @test one_factor.certificate_kind === :one_dimensional_factor
    end

    @testset "matrix spectra, permutation invariance, and no repair" begin
        rng = MersenneTwister(0x41505054)
        spectrum = [0.45, 0.35, 0.1, 0.1]
        diagonal_state = Diagonal(spectrum)
        unitary = _abs_ppt_random_unitary(rng, 4)
        rotated_state = unitary * diagonal_state * adjoint(unitary)

        diagonal_result = AbsPPTQET.is_abs_ppt(diagonal_state; dims=(2, 2))
        rotated_result = AbsPPTQET.is_abs_ppt(
            Hermitian(rotated_state); dims=(2, 2), allow_densify=true
        )
        @test diagonal_result.status === rotated_result.status
        @test diagonal_result.verdict === rotated_result.verdict === true
        @test rotated_result.spectrum ≈ sort(spectrum; rev=true)
        @test diagonal_result.family.constraints[1].matrix ≈
            rotated_result.family.constraints[1].matrix atol = 2e-14
        @test AbsPPTQET.is_abs_ppt(spectrum; dims=(4, 1)).verdict === true

        scaled = AbsPPTQET.is_abs_ppt(7 .* spectrum; dims=(2, 2))
        @test scaled.verdict === true
        @test scaled.status === AbsPPTQET.AbsolutePPTExhaustiveCertified

        pure_bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_state = pure_bell * adjoint(pure_bell)
        partial = AbsPPTQET.partial_transpose(bell_state, (2, 2); systems=(2,))
        @test minimum(eigvals(Hermitian(partial))) ≈ -0.5 atol = 2e-15
        @test AbsPPTQET.is_abs_ppt(bell_state; dims=(2, 2)).verdict === false

        for _ in 1:20
            trial_unitary = _abs_ppt_random_unitary(rng, 4)
            trial = trial_unitary * diagonal_state * adjoint(trial_unitary)
            transposed = AbsPPTQET.partial_transpose(trial, (2, 2); systems=(2,))
            @test minimum(eigvals(Hermitian(transposed))) >= -4e-14
        end

        near_hermitian = [0.5 eps(Float64); 0.0 0.5]
        matrix_boundary = AbsPPTQET.is_abs_ppt(near_hermitian; dims=(1, 2))
        @test matrix_boundary.status === AbsPPTQET.AbsolutePPTNumericalBoundary
        @test matrix_boundary.verdict === nothing
        @test near_hermitian == [0.5 eps(Float64); 0.0 0.5]

        @test_throws ArgumentError AbsPPTQET.is_abs_ppt([0.5 0.1; 0.0 0.5]; dims=(1, 2))
        @test_throws DomainError AbsPPTQET.is_abs_ppt([1.1, -0.1, 0.0, 0.0]; dims=(2, 2))
        negative_boundary = AbsPPTQET.is_abs_ppt(
            [0.5, 0.5, 0.0, -eps(Float64)]; dims=(2, 2)
        )
        @test negative_boundary.status === AbsPPTQET.AbsolutePPTNumericalBoundary
        @test negative_boundary.verdict === nothing
        @test negative_boundary.spectrum[end] == -eps(Float64)
    end

    @testset "ordering feasibility model and exact verification" begin
        family = AbsPPTQET.abs_ppt_constraints(collect(9.0:-1:1); dims=(3, 3))
        for constraint in family.constraints
            problem = AbsPPTQET.abs_ppt_ordering_program(
                constraint.ordering;
                limits=AbsPPTQET.OptimizationLimits(
                    max_variables=2, max_intervals=10, max_model_entries=100
                ),
            )
            @test problem.variable_count == 2
            @test problem.sense === :feasibility
            @test problem.metadata.formulation === :hildebrand_log_product_ordering
            @test length(problem.intervals) == 7
            certificate = AbsPPTQET._abs_ppt_deterministic_ordering_certificate(
                constraint.ordering; max_iterations=10_000
            )
            @test certificate !== nothing
            @test certificate.minimum_gap > 0
            @test certificate.minimum_product_margin > 0
            scale_factor =
                1 /
                Float64(min(certificate.minimum_gap, certificate.minimum_product_margin))
            @test AbsPPTQET.primal_residual(
                problem, scale_factor .* Float64.(certificate.log_gaps); allow_densify=true
            ) == 0
        end
        @test_throws ArgumentError AbsPPTQET.abs_ppt_ordering_program(
            AbsPPTQET.AbsPPTOrdering(Tuple{Int,Int}[], Tuple{Int,Int}[])
        )
    end
end
