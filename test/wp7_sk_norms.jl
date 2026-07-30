using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const SKNormQET = QuantumEntanglementTools

if !isdefined(SKNormQET, :TopKPNormDualEpigraph)
    Base.include(
        SKNormQET,
        joinpath(@__DIR__, "..", "src", "optimization", "top_k_p_norm_dual_epigraph.jl"),
    )
end
if !isdefined(SKNormQET, :SKOperatorNormResult)
    Base.include(
        SKNormQET, joinpath(@__DIR__, "..", "src", "optimization", "sk_operator_norm.jl")
    )
end
if !isdefined(SKNormQET, :BlockPositivityResult)
    Base.include(
        SKNormQET, joinpath(@__DIR__, "..", "src", "entanglement", "block_positivity.jl")
    )
end

@testset "WP7 top-k p-norm dual solver-neutral epigraph" begin
    variable = SKNormQET.complex_affine_variable(:X, 2, 3; coefficient_type=Float32)
    atom = SKNormQET.top_k_p_norm_dual_epigraph(variable, 20, 3)
    @test atom isa SKNormQET.TopKPNormDualEpigraph{Float32,Int,Float64}
    @test atom.matrix === variable
    @test atom.k == 2
    @test atom.p == 3
    @test atom.conjugate_order == 3 / 2
    @test atom.singular_value_count == 2
    @test atom.dilation.dimension == 5
    @test atom.psd_block_count == 4
    @test atom.auxiliary_variable_count == 59
    @test atom.scalar_constraint_count == 9

    endpoint_one = SKNormQET.top_k_p_norm_dual_epigraph(variable, 1, 1)
    endpoint_inf = SKNormQET.top_k_p_norm_dual_epigraph(variable, 1, Inf)
    @test isinf(endpoint_one.conjugate_order)
    @test endpoint_inf.conjugate_order == 1
    @test endpoint_one.scalar_constraint_count == 6
    @test endpoint_inf.scalar_constraint_count == 4

    vector = SKNormQET.complex_affine_variable(:v, 1, 4)
    vector_atom = SKNormQET.top_k_p_norm_dual_epigraph(vector, 2, 2)
    @test vector_atom.singular_value_count == 4
    @test vector_atom.dilation.dimension == 8
    coordinates = [3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    value = vec(Matrix(SKNormQET.evaluate_affine(vector, coordinates)))
    @test SKNormQET.top_k_p_norm_dual(value, 2, 2) ≈ sqrt(18)

    huge_k = big(typemax(Int)) + 10
    @test SKNormQET.top_k_p_norm_dual_epigraph(vector, huge_k, 2).k == 4
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(vector, 0, 2)
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(vector, 1, 0.5)
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(vector, 1, NaN)
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(
        variable, 1, 2; limits=SKNormQET.OptimizationLimits(max_variables=60)
    )
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(
        variable, 1, 2; limits=SKNormQET.OptimizationLimits(max_psd_blocks=3)
    )
    @test_throws ArgumentError SKNormQET.top_k_p_norm_dual_epigraph(
        variable, 1, 2; limits=SKNormQET.OptimizationLimits(max_intervals=7)
    )
    @test_throws ArgumentError SKNormQET.add_top_k_p_norm_dual_epigraph!(
        nothing, atom, zeros(variable.variable_count); allow_densify=true
    )
end

@testset "WP7 S(k) exact branches and witnesses" begin
    zero_rng = MersenneTwister(701)
    zero_control = copy(zero_rng)
    zero_result = SKNormQET.sk_operator_norm(zero_rng, zeros(4, 4); dims=(2, 2))
    @test zero_result.status === SKNormQET.SKOperatorNormExact
    @test zero_result.exact
    @test zero_result.lower_bound == 0
    @test zero_result.upper_bound == 0
    @test zero_result.lower_witness.validated
    @test rand(zero_rng) == rand(zero_control)

    zero_block_rng = MersenneTwister(7001)
    zero_block_control = copy(zero_block_rng)
    zero_block = SKNormQET.is_block_positive(zero_block_rng, zeros(4, 4); k=2, dims=(2, 2))
    @test zero_block.status === SKNormQET.BlockPositivityCertified
    @test zero_block.verdict === true
    @test zero_block.exact
    @test zero_block.certificate_kind === :zero_operator_exact
    @test rand(zero_block_rng) == rand(zero_block_control)

    diagonal = Diagonal([4.0, 3.0, 2.0, 1.0])
    full_rng = MersenneTwister(702)
    full_control = copy(full_rng)
    full = SKNormQET.sk_operator_norm(full_rng, diagonal; k=99, dims=(2, 2))
    @test full.status === SKNormQET.SKOperatorNormExact
    @test full.lower_bound == 4
    @test full.upper_bound == 4
    @test full.k == 2
    @test full.lower_witness.left_schmidt_rank <= 2
    @test full.lower_witness.right_schmidt_rank <= 2
    @test rand(full_rng) == rand(full_control)

    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    product = ComplexF64[1, 0, 0, 0]
    rank_one = 2bell * product'
    rank_one_result = SKNormQET.sk_operator_norm(
        MersenneTwister(703), rank_one; k=1, dims=(2, 2)
    )
    @test rank_one_result.status === SKNormQET.SKOperatorNormExact
    @test rank_one_result.exact
    @test rank_one_result.lower_bound ≈ sqrt(2) atol = 2e-14
    @test rank_one_result.upper_bound ≈ sqrt(2) atol = 2e-14
    @test rank_one_result.lower_witness.validated
    @test rank_one_result.lower_witness.left_schmidt_rank == 1
    @test rank_one_result.lower_witness.right_schmidt_rank == 1

    target_rng = MersenneTwister(704)
    target_control = copy(target_rng)
    above = SKNormQET.sk_operator_norm(target_rng, diagonal; k=1, dims=(2, 2), target=3.5)
    @test above.status === SKNormQET.SKOperatorNormTargetAbove
    @test above.target_decision === :above
    @test above.random_restarts == 0
    @test rand(target_rng) == rand(target_control)

    below_rng = MersenneTwister(705)
    below_control = copy(below_rng)
    below = SKNormQET.sk_operator_norm(
        below_rng, Matrix{Float64}(I, 4, 4); k=1, dims=(2, 2), target=2
    )
    @test below.status === SKNormQET.SKOperatorNormTargetBelow
    @test below.target_decision === :below
    @test rand(below_rng) == rand(below_control)
end

@testset "WP7 explicit RNG leaves the global stream untouched" begin
    Random.seed!(0x534b4e4f524d)
    expected_after_sk = rand(UInt64)
    expected_after_block = rand(UInt64)

    Random.seed!(0x534b4e4f524d)
    SKNormQET.sk_operator_norm(
        MersenneTwister(0x7011), Matrix{Float64}(I, 4, 4); k=1, dims=(2, 2)
    )
    @test rand(UInt64) == expected_after_sk

    SKNormQET.is_block_positive(
        MersenneTwister(0x7012), Matrix{Float64}(I, 4, 4); k=1, dims=(2, 2)
    )
    @test rand(UInt64) == expected_after_block
end

@testset "WP7 S(k) bounded search and resource semantics" begin
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    positive = bell * bell' + 0.2Matrix{ComplexF64}(I, 4, 4)
    rng = MersenneTwister(706)
    control = copy(rng)
    bounded = SKNormQET.sk_operator_norm(
        rng, positive; k=1, dims=(2, 2), strength=1, max_restarts=3, max_iterations=15
    )
    @test bounded.status in
        (SKNormQET.SKOperatorNormBounds, SKNormQET.SKOperatorNormResourceLimited)
    @test bounded.lower_bound <= 0.7 + bounded.tolerance
    @test bounded.upper_bound >= 0.7 - bounded.tolerance
    @test bounded.lower_witness.validated
    @test bounded.lower_witness.quadratic
    @test bounded.lower_witness.left_schmidt_rank <= 1
    @test bounded.witness_lower_bound ≈
        abs(dot(bounded.lower_witness.left, positive * bounded.lower_witness.right))
    @test rand(rng) != rand(control)
    @test bounded.optimization_result.status === SKNormQET.OptimizationBackendUnavailable

    no_random_rng = MersenneTwister(707)
    no_random_control = copy(no_random_rng)
    no_random = SKNormQET.sk_operator_norm(
        no_random_rng,
        positive;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=3,
        max_iterations=15,
        max_work=100,
    )
    @test no_random.status === SKNormQET.SKOperatorNormResourceLimited
    @test no_random.random_restarts == 0
    @test any(occursin("skipped", warning) for warning in no_random.warnings)
    @test rand(no_random_rng) == rand(no_random_control)

    zero_iterations_rng = MersenneTwister(7071)
    zero_iterations_control = copy(zero_iterations_rng)
    zero_iterations = SKNormQET.sk_operator_norm(
        zero_iterations_rng, positive; k=1, dims=(2, 2), strength=2, max_iterations=0
    )
    @test zero_iterations.random_restarts == 0
    @test zero_iterations.random_iterations == 0
    @test rand(zero_iterations_rng) == rand(zero_iterations_control)

    nonnormal = ComplexF64[0 2 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
    general = SKNormQET.sk_operator_norm(
        MersenneTwister(708),
        nonnormal;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=2,
        max_iterations=8,
    )
    @test general.lower_bound <= general.upper_bound + general.tolerance
    @test general.lower_witness.validated
    @test !general.lower_witness.quadratic

    sparse_positive = sparse(positive)
    @test_throws ArgumentError SKNormQET.sk_operator_norm(
        MersenneTwister(709), sparse_positive; dims=(2, 2)
    )
    sparse_result = SKNormQET.sk_operator_norm(
        MersenneTwister(709), sparse_positive; dims=(2, 2), strength=0, allow_densify=true
    )
    @test sparse_result.lower_bound <= sparse_result.upper_bound
    @test_throws ArgumentError SKNormQET.sk_operator_norm(
        MersenneTwister(710), positive; dims=(2, 2), max_dense_entries=15
    )
    @test_throws DimensionMismatch SKNormQET.sk_operator_norm(
        MersenneTwister(710), positive; dims=(3, 2)
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm(
        MersenneTwister(710), positive; dims=(2, 2), strength=-1
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm(
        MersenneTwister(710), fill(ComplexF64(Inf), 4, 4); dims=(2, 2)
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm(
        MersenneTwister(710), BigFloat.(Matrix{Float64}(I, 4, 4)); dims=(2, 2)
    )
end

@testset "WP7 package-owned S(k) SDP models" begin
    positive = Diagonal([4.0, 3.0, 2.0, 1.0])
    ppt = SKNormQET.sk_operator_norm_problem(positive; k=1, dims=(2, 2), allow_densify=true)
    @test ppt isa SKNormQET.SKOperatorNormProblem
    @test ppt.formulation === :ppt_relaxation
    @test ppt.program.sense === :maximize
    @test ppt.program.variable_count == 16
    @test length(ppt.program.psd_constraints) == 2
    @test length(ppt.program.intervals) == 1
    @test ppt.metadata.exact_low_dimension_theorem
    @test ppt.program.known_feasible_point == zeros(16)

    reduction = SKNormQET.sk_operator_norm_problem(
        Diagonal([6.0, 5, 4, 3, 2, 1]); k=2, dims=(2, 3), allow_densify=true
    )
    @test reduction.formulation === :schmidt_number_reduction_relaxation
    @test reduction.k == 2

    hierarchy = SKNormQET.sk_operator_norm_problem(
        positive; k=1, dims=(2, 2), hierarchy_level=2, allow_densify=true
    )
    @test hierarchy.formulation === :bosonic_ppt_symmetric_extension
    @test hierarchy.ambient_dimension == 8
    @test hierarchy.variable_matrix_dimension == 6
    @test hierarchy.program.variable_count == 36
    @test size(hierarchy.bosonic_basis) == (8, 6)
    @test hierarchy.state_view.dimension == 4
    @test hierarchy.extension_view.dimension == 8

    @test_throws ArgumentError SKNormQET.sk_operator_norm_problem(
        positive; k=2, dims=(2, 2), hierarchy_level=2, allow_densify=true
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm_problem(
        Diagonal([1.0, 1, 1, -0.1]); dims=(2, 2), allow_densify=true
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm_problem(
        ComplexF64[1 1; 0 1]; dims=(1, 2), allow_densify=true
    )
    @test_throws ArgumentError SKNormQET.sk_operator_norm_problem(
        positive;
        dims=(2, 2),
        allow_densify=true,
        limits=SKNormQET.OptimizationLimits(max_variables=15),
    )
end

@testset "WP7 block positivity certificates and unknown boundaries" begin
    positive_rng = MersenneTwister(711)
    positive_control = copy(positive_rng)
    positive = SKNormQET.is_block_positive(
        positive_rng, Matrix{Float64}(I, 4, 4); dims=(2, 2)
    )
    @test positive.status === SKNormQET.BlockPositivityCertified
    @test positive.verdict === true
    @test positive.certified
    @test positive.certificate_kind === :positive_definite_spectral_certificate
    @test rand(positive_rng) == rand(positive_control)

    negative_rng = MersenneTwister(712)
    negative_control = copy(negative_rng)
    negative = SKNormQET.is_block_positive(
        negative_rng, Diagonal([-1.0, 2, 2, 2]); dims=(2, 2)
    )
    @test negative.status === SKNormQET.BlockPositivityViolated
    @test negative.verdict === false
    @test negative.witness.validated
    @test negative.witness.schmidt_rank == 1
    @test negative.witness.expectation == -1
    @test rand(negative_rng) == rand(negative_control)

    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    boundary_witness = 0.5Matrix{ComplexF64}(I, 4, 4) - bell * bell'
    block_positive = SKNormQET.is_block_positive(
        MersenneTwister(713), boundary_witness; dims=(2, 2), strength=0
    )
    @test block_positive.status in
        (SKNormQET.BlockPositivityUnknown, SKNormQET.BlockPositivityNumericalBoundary)
    @test block_positive.verdict === nothing
    @test !block_positive.certified

    boundary = SKNormQET.is_block_positive(
        MersenneTwister(714), Diagonal([0.0, 1, 1, 1]); k=2, dims=(2, 2)
    )
    @test boundary.status === SKNormQET.BlockPositivityNumericalBoundary
    @test boundary.verdict === nothing
    @test !boundary.certified

    entangled_negative = 0.4Matrix{ComplexF64}(I, 4, 4) - bell * bell'
    found = SKNormQET.is_block_positive(
        MersenneTwister(715),
        entangled_negative;
        dims=(2, 2),
        strength=1,
        max_restarts=8,
        max_iterations=30,
    )
    @test found.status in (
        SKNormQET.BlockPositivityViolated,
        SKNormQET.BlockPositivityUnknown,
        SKNormQET.BlockPositivityBackendFailure,
    )
    found.verdict === false && @test found.witness.validated
    @test found.verdict !== true

    @test_throws ArgumentError SKNormQET.is_block_positive(
        MersenneTwister(716), ComplexF64[1 1; 0 1]; dims=(1, 2)
    )
    @test_throws ArgumentError SKNormQET.is_block_positive(
        MersenneTwister(716), sparse(Matrix{Float64}(I, 4, 4)); dims=(2, 2)
    )
end
