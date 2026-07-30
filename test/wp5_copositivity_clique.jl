using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const CopCliqueQET = QuantumEntanglementTools

if !isdefined(CopCliqueQET, :CopositivityResult)
    Base.include(
        CopCliqueQET,
        joinpath(@__DIR__, "..", "src", "optimization", "copositivity_clique.jl"),
    )
end

struct _CopCliqueZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_CopCliqueZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_CopCliqueZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_CopCliqueZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_CopCliqueZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function cycle_adjacency(vertices::Int)
    adjacency = zeros(Int, vertices, vertices)
    for vertex in 1:vertices
        neighbor = mod1(vertex + 1, vertices)
        adjacency[vertex, neighbor] = 1
        adjacency[neighbor, vertex] = 1
    end
    return adjacency
end

@testset "WP5 copositivity and clique-number contracts" begin
    @testset "exact copositivity certificates and boundaries" begin
        entrywise = [2 -0; 0 3]
        entrywise_copy = copy(entrywise)
        entrywise_rng = MersenneTwister(10)
        untouched_rng = copy(entrywise_rng)
        entrywise_result = CopCliqueQET.copositivity_criterion(
            entrywise_rng, entrywise; inner_samples=100
        )
        @test entrywise_result isa CopCliqueQET.CopositivityResult
        @test entrywise_result.status === CopCliqueQET.CopositivityCertifiedTrue
        @test entrywise_result.verdict === true
        @test entrywise_result.certified
        @test entrywise_result.certificate_kind === :entrywise_nonnegative
        @test entrywise_result.lower_bound == 0
        @test entrywise_result.upper_bound == 2
        @test entrywise_result.polynomial === nothing
        @test entrywise_result.hierarchy_result === nothing
        @test entrywise_result.samples_requested == 100
        @test entrywise_result.samples_evaluated == 0
        @test entrywise == entrywise_copy
        @test rand(entrywise_rng) == rand(untouched_rng)
        @test occursin("verdict=true", sprint(show, entrywise_result))

        exact_psd = Rational{Int}[1 -1; -1 1]
        psd_result = CopCliqueQET.copositivity_criterion(MersenneTwister(11), exact_psd)
        @test psd_result.status === CopCliqueQET.CopositivityCertifiedTrue
        @test psd_result.verdict === true
        @test psd_result.certified
        @test psd_result.certificate_kind === :exact_positive_semidefinite
        @test psd_result.psd_diagnostic.status === CopCliqueQET.MatrixPredicateSatisfied
        @test psd_result.polynomial === nothing

        negative_diagonal = Rational{Int}[-1 0; 0 2]
        diagonal_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(12), negative_diagonal
        )
        @test diagonal_result.status === CopCliqueQET.CopositivityCertifiedFalse
        @test diagonal_result.verdict === false
        @test diagonal_result.certified
        @test diagonal_result.witness.source === :coordinate_ray
        @test diagonal_result.witness.simplex_vector == [1, 0]
        @test diagonal_result.witness.value == -1
        @test diagonal_result.upper_bound == -1
        @test diagonal_result.polynomial === nothing
        @test diagonal_result.samples_evaluated == 0

        negative_pair = [1 -2; -2 1]
        pair_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(13), negative_pair
        )
        @test pair_result.status === CopCliqueQET.CopositivityCertifiedFalse
        @test pair_result.verdict === false
        @test pair_result.certificate_kind === :exact_two_coordinate_ray_witness
        @test pair_result.witness.simplex_vector == [1 // 2, 1 // 2]
        @test pair_result.witness.value == -1 // 2
        @test pair_result.polynomial === nothing
        pair_polynomial = CopCliqueQET.copositive_polynomial(negative_pair)
        @test CopCliqueQET.evaluate_polynomial(
            pair_polynomial, sqrt.(Float64.(pair_result.witness.simplex_vector))
        ) ≈ -0.5

        # Cheap exact branches run before polynomial construction. In
        # particular, they are not subject to hierarchy term limits and
        # arithmetic stays exact even at the edge of the input integer type.
        huge_entrywise = fill(typemax(Int), 2, 2)
        huge_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(131), huge_entrywise
        )
        @test huge_result.verdict === true
        @test huge_result.certificate_kind === :entrywise_nonnegative
        @test huge_result.polynomial === nothing

        large_identity = Matrix{Int}(I, 38, 38)
        large_identity_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(132), large_identity
        )
        @test large_identity_result.verdict === true
        @test large_identity_result.polynomial === nothing

        large_negative = copy(large_identity)
        large_negative[1, 1] = -1
        large_negative_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(133), large_negative
        )
        @test large_negative_result.verdict === false
        @test large_negative_result.witness.source === :coordinate_ray
        @test large_negative_result.polynomial === nothing

        # The pinned routine returns true here because it compares its lower
        # bound with a fixed -1e-9 threshold. The native contract preserves
        # the tolerance boundary as unknown.
        boundary_matrix = [-5.0e-10 0.0; 0.0 1.0]
        boundary_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(14), boundary_matrix
        )
        @test boundary_result.status === CopCliqueQET.CopositivityNumericalBoundary
        @test boundary_result.verdict === nothing
        @test !boundary_result.certified
        @test boundary_result.tolerance > abs(boundary_matrix[1, 1])
        @test occursin("tolerance band", boundary_result.message)

        robust_float = [-1.0e-4 0.0; 0.0 1.0]
        robust_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(15), robust_float
        )
        @test robust_result.status === CopCliqueQET.CopositivityCertifiedFalse
        @test robust_result.witness.value ==
            rationalize(BigInt, robust_float[1, 1]; tol=0.0)
    end

    @testset "hierarchy evidence, RNG ownership, and limits" begin
        horn = [
            1 -1 1 1 -1
            -1 1 -1 1 1
            1 -1 1 -1 1
            1 1 -1 1 -1
            -1 1 1 -1 1
        ]
        unavailable = CopCliqueQET.copositivity_criterion(MersenneTwister(20), horn)
        @test unavailable.status === CopCliqueQET.CopositivityBackendUnavailable
        @test unavailable.verdict === nothing
        @test unavailable.lower_bound === nothing
        @test unavailable.upper_bound === nothing
        @test unavailable.hierarchy_result.status ===
            CopCliqueQET.OptimizationBackendUnavailable
        @test unavailable.psd_diagnostic.status === CopCliqueQET.MatrixPredicateViolated

        first = CopCliqueQET.copositivity_criterion(
            MersenneTwister(21),
            horn;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=12,
        )
        second = CopCliqueQET.copositivity_criterion(
            MersenneTwister(21),
            horn;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=12,
        )
        @test first.status === CopCliqueQET.CopositivityHierarchyUnknown
        @test first.verdict === nothing
        @test first.lower_kind === :hierarchy_outer_bound
        @test first.upper_kind === :sampled_feasible
        @test first.lower_bound == second.lower_bound
        @test first.upper_bound == second.upper_bound
        @test first.hierarchy_result.best_point == second.hierarchy_result.best_point
        @test first.samples_requested == 12
        @test first.samples_evaluated == 12

        # No coordinate, equal-pair, or uniform ray detects this matrix.
        # A sample is rationalized and re-evaluated exactly before it is
        # allowed to become a negative certificate.
        sampled_negative_matrix = [1.0 -2.2; -2.2 4.0]
        sampled_negative = CopCliqueQET.copositivity_criterion(
            MersenneTwister(2),
            sampled_negative_matrix;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=8,
        )
        @test sampled_negative.status === CopCliqueQET.CopositivityCertifiedFalse
        @test sampled_negative.verdict === false
        @test sampled_negative.certified
        @test sampled_negative.certificate_kind === :exact_rationalized_sample_witness
        @test sampled_negative.witness.source === :rationalized_sample
        @test sum(sampled_negative.witness.simplex_vector) == 1
        @test sampled_negative.witness.value < 0
        @test sampled_negative.samples_evaluated == 8

        Random.seed!(20260730)
        expected_global_draw = rand()
        Random.seed!(20260730)
        CopCliqueQET.copositivity_criterion(
            MersenneTwister(22), horn; hierarchy=:nosdp, allow_densify=true, inner_samples=2
        )
        @test rand() == expected_global_draw

        numeric_psd = [1.0 -0.25; -0.25 1.0]
        numeric_psd_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(23), numeric_psd
        )
        @test numeric_psd_result.psd_diagnostic.status ===
            CopCliqueQET.MatrixPredicateSatisfied
        @test numeric_psd_result.verdict === nothing
        @test numeric_psd_result.status === CopCliqueQET.CopositivityBackendUnavailable

        sparse_psd = sparse(Rational{Int}[1 -1; -1 1])
        sparse_unknown = CopCliqueQET.copositivity_criterion(
            MersenneTwister(24), sparse_psd
        )
        @test sparse_unknown.verdict === nothing
        @test sparse_unknown.psd_diagnostic === nothing
        sparse_true = CopCliqueQET.copositivity_criterion(
            MersenneTwister(24), sparse_psd; allow_densify=true
        )
        @test sparse_true.verdict === true
        @test sparse_true.certificate_kind === :exact_positive_semidefinite

        guarded_rng = MersenneTwister(25)
        untouched_guarded_rng = copy(guarded_rng)
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            guarded_rng, horn; max_matrix_dimension=4
        )
        @test rand(guarded_rng) == rand(untouched_guarded_rng)
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(26), horn; inner_samples=3, max_samples=2
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(27), horn; hierarchy=:nosdp
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(28), horn; inner_samples=1
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(29),
            horn;
            hierarchy=:nosdp,
            backend=CopCliqueQET.JuMPBackend(() -> nothing),
            allow_densify=true,
        )
    end

    @testset "certified graph bounds and Motzkin--Straus evidence" begin
        isolated = zeros(Int, 3, 3)
        isolated_rng = MersenneTwister(30)
        untouched_isolated_rng = copy(isolated_rng)
        isolated_result = CopCliqueQET.clique_number_bounds(
            isolated_rng, isolated; inner_samples=50
        )
        @test isolated_result.status === CopCliqueQET.CliqueNumberExact
        @test isolated_result.lower_bound == 1
        @test isolated_result.upper_bound == 1
        @test isolated_result.exact
        @test isolated_result.bounds_certified
        @test isolated_result.edge_count == 0
        @test isolated_result.best_clique == [1]
        @test isolated_result.polynomial === nothing
        @test isolated_result.hierarchy_result === nothing
        @test isolated_result.samples_evaluated == 0
        @test rand(isolated_rng) == rand(untouched_isolated_rng)

        complete = ones(Int, 4, 4) - Matrix{Int}(I, 4, 4)
        complete_result = CopCliqueQET.clique_number_bounds(MersenneTwister(31), complete)
        @test complete_result.status === CopCliqueQET.CliqueNumberExact
        @test complete_result.lower_bound == 4
        @test complete_result.upper_bound == 4
        @test complete_result.best_clique == collect(1:4)
        @test complete_result.edge_count == 6
        @test complete_result.polynomial === nothing

        path = [
            0 1 0 0
            1 0 1 0
            0 1 0 1
            0 0 1 0
        ]
        path_result = CopCliqueQET.clique_number_bounds(MersenneTwister(32), path)
        @test path_result.exact
        @test path_result.lower_bound == path_result.upper_bound == 2
        @test path_result.polynomial === nothing
        @test path_result.upper_certificate.edge_count_bound == 3
        @test path_result.upper_certificate.maximum_degree == 2
        @test path_result.upper_certificate.maximum_degree_bound == 3
        @test path_result.upper_certificate.greedy_coloring_bound == 2
        @test path_result.upper_certificate_kinds == (:greedy_coloring,)
        path_colors = path_result.upper_certificate.coloring
        @test all(
            path[left, right] == 0 || path_colors[left] != path_colors[right] for
            left in axes(path, 1), right in axes(path, 2)
        )

        large_edgeless = zeros(Int, 38, 38)
        large_edgeless_result = CopCliqueQET.clique_number_bounds(
            MersenneTwister(134), large_edgeless
        )
        @test (large_edgeless_result.lower_bound, large_edgeless_result.upper_bound) ==
            (1, 1)
        @test large_edgeless_result.polynomial === nothing

        large_complete = ones(Int, 38, 38) - Matrix{Int}(I, 38, 38)
        large_complete_result = CopCliqueQET.clique_number_bounds(
            MersenneTwister(135), large_complete
        )
        @test (large_complete_result.lower_bound, large_complete_result.upper_bound) ==
            (38, 38)
        @test large_complete_result.polynomial === nothing

        cycle_five = cycle_adjacency(5)
        unavailable = CopCliqueQET.clique_number_bounds(MersenneTwister(33), cycle_five)
        @test unavailable.status === CopCliqueQET.CliqueNumberBackendUnavailable
        @test unavailable.lower_bound == 2
        @test unavailable.upper_bound == 3
        @test !unavailable.exact
        @test unavailable.bounds_certified
        @test unavailable.lower_certificate_kind === :greedy_clique
        @test unavailable.upper_certificate_kinds ==
            (:edge_count, :maximum_degree, :greedy_coloring)
        @test unavailable.hierarchy_result.status ===
            CopCliqueQET.OptimizationBackendUnavailable
        @test occursin("bounds=2:3", sprint(show, unavailable))

        sampled_one = CopCliqueQET.clique_number_bounds(
            MersenneTwister(34),
            cycle_five;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=16,
        )
        sampled_two = CopCliqueQET.clique_number_bounds(
            MersenneTwister(34),
            cycle_five;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=16,
        )
        @test sampled_one.status === CopCliqueQET.CliqueNumberNumericalHierarchy
        @test sampled_one.lower_bound == 2
        @test sampled_one.upper_bound == 3
        @test sampled_one.continuous_upper_bound == sampled_two.continuous_upper_bound
        @test sampled_one.continuous_lower_bound == sampled_two.continuous_lower_bound
        @test sampled_one.motzkin_straus_witness.simplex_vector ==
            sampled_two.motzkin_straus_witness.simplex_vector
        @test sum(sampled_one.motzkin_straus_witness.simplex_vector) == 1
        @test sampled_one.motzkin_straus_witness.implied_lower_bound == 2
        @test sampled_one.uncertified_upper_candidate !== nothing
        @test sampled_one.samples_requested == 16
        @test sampled_one.samples_evaluated == 16

        # Certificate payloads are internally validated and use owned,
        # read-only arrays. Public constructor calls cannot forge a result, and
        # modifying a storage copy cannot invalidate the retained certificate.
        sealed_cop_result = CopCliqueQET.copositivity_criterion(
            MersenneTwister(136), [1 -2; -2 1]
        )
        copositivity_witness_fields = Tuple(
            getfield(sealed_cop_result.witness, name) for
            name in fieldnames(typeof(sealed_cop_result.witness))
        )
        @test_throws MethodError CopCliqueQET.CopositivityWitness(
            copositivity_witness_fields...
        )
        forged_token = CopCliqueQET._ValidatedConstructorToken()
        @test_throws ArgumentError CopCliqueQET.CopositivityWitness(
            forged_token, copositivity_witness_fields...
        )
        @test_throws Base.CanonicalIndexError setindex!(
            sealed_cop_result.witness.simplex_vector, 0 // 1, 1
        )
        witness_storage_copy = sealed_cop_result.witness.simplex_vector.storage
        witness_storage_copy[1] = 0
        @test sealed_cop_result.witness.simplex_vector == [1 // 2, 1 // 2]
        @test sealed_cop_result.witness.value == -1 // 2

        copositivity_result_fields = Tuple(
            getfield(sealed_cop_result, name) for
            name in fieldnames(typeof(sealed_cop_result))
        )
        @test_throws MethodError CopCliqueQET.CopositivityResult(
            copositivity_result_fields...
        )
        @test_throws ArgumentError CopCliqueQET.CopositivityResult(
            forged_token, copositivity_result_fields...
        )

        motzkin_witness = sampled_one.motzkin_straus_witness
        motzkin_witness_fields = Tuple(
            getfield(motzkin_witness, name) for name in fieldnames(typeof(motzkin_witness))
        )
        @test_throws MethodError CopCliqueQET.MotzkinStrausWitness(
            motzkin_witness_fields...
        )
        @test_throws ArgumentError CopCliqueQET.MotzkinStrausWitness(
            forged_token, motzkin_witness_fields...
        )
        @test_throws Base.CanonicalIndexError setindex!(
            motzkin_witness.simplex_vector, 0 // 1, 1
        )
        motzkin_storage_copy = motzkin_witness.simplex_vector.storage
        motzkin_storage_copy[1] = 0
        @test sum(motzkin_witness.simplex_vector) == 1

        clique_result_fields = Tuple(
            getfield(path_result, name) for name in fieldnames(typeof(path_result))
        )
        @test_throws MethodError CopCliqueQET.CliqueNumberResult(clique_result_fields...)
        @test_throws ArgumentError CopCliqueQET.CliqueNumberResult(
            forged_token, clique_result_fields...
        )
        @test_throws Base.CanonicalIndexError setindex!(path_result.best_clique, 1, 2)
        certified_clique = copy(path_result.best_clique)
        clique_storage_copy = path_result.best_clique.storage
        clique_storage_copy[2] = clique_storage_copy[1]
        @test path_result.best_clique == certified_clique
        @test_throws Base.CanonicalIndexError setindex!(
            path_result.upper_certificate.coloring, 1, 1
        )
        coloring_storage_copy = path_result.upper_certificate.coloring.storage
        coloring_storage_copy[1] = coloring_storage_copy[2]
        @test all(
            path[left, right] == 0 ||
                path_result.upper_certificate.coloring[left] !=
                path_result.upper_certificate.coloring[right] for
            left in axes(path, 1), right in axes(path, 2)
        )

        float_cycle = Float64.(cycle_five)
        float_result = CopCliqueQET.clique_number_bounds(MersenneTwister(35), float_cycle)
        @test float_result.lower_bound == 2
        @test float_result.upper_bound == 3
        @test float_result.tolerance > 0

        sparse_cycle = sparse(cycle_five)
        sparse_result = CopCliqueQET.clique_number_bounds(MersenneTwister(36), sparse_cycle)
        @test sparse_result.lower_bound == 2
        @test sparse_result.upper_bound == 3
        @test sparse_result.polynomial.coefficients isa SparseVector
    end

    @testset "strict graph and matrix validation" begin
        @test_throws DimensionMismatch CopCliqueQET.copositivity_criterion(
            MersenneTwister(40), ones(2, 3)
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(41), zeros(0, 0)
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(42), [1.0 2.0; 3.0 1.0]
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(43), [1.0 NaN; NaN 1.0]
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(44), Bool[1 0; 0 1]
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(44), fill(π, 1, 1)
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(45), Rational{Int}[1 -1; -1 1]; atol=1 // 10
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(46), [1 0; 0 1]; hierarchy=:invalid
        )
        @test_throws ArgumentError CopCliqueQET.copositivity_criterion(
            MersenneTwister(47), _CopCliqueZeroBasedMatrix([1 0; 0 1])
        )

        @test_throws DimensionMismatch CopCliqueQET.clique_number_bounds(
            MersenneTwister(50), zeros(Int, 2, 3)
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(51), zeros(Int, 0, 0)
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(52), [0 1; 0 0]
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(53), [1 0; 0 0]
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(54), [0.0 0.5; 0.5 0.0]
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(55), [0.0 Inf; Inf 0.0]
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(56), ComplexF64[0 1; 1 0]
        )
        bool_graph = CopCliqueQET.clique_number_bounds(MersenneTwister(56), Bool[0 1; 1 0])
        @test (bool_graph.lower_bound, bool_graph.upper_bound) == (2, 2)
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(57), _CopCliqueZeroBasedMatrix([0 1; 1 0])
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(58), cycle_adjacency(5); max_matrix_entries=24
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(59), cycle_adjacency(5); max_greedy_work=1
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(60), cycle_adjacency(5); hierarchy=:nosdp, allow_densify=false
        )
        @test_throws ArgumentError CopCliqueQET.clique_number_bounds(
            MersenneTwister(61), cycle_adjacency(5); inner_samples=3, max_samples=2
        )
    end
end
