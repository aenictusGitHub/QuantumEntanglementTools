using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const SymExtQET = QuantumEntanglementTools

if !isdefined(SymExtQET, :SymmetricExtensionStatus)
    Base.include(
        SymExtQET,
        joinpath(@__DIR__, "..", "src", "entanglement", "symmetric_extensions.jl"),
    )
end

function symext_bell_state(::Type{T}=Float64) where {T<:AbstractFloat}
    vector = T[inv(sqrt(T(2))), 0, 0, inv(sqrt(T(2)))]
    return vector * adjoint(vector)
end

function assert_ppt_state(result)
    @test result.status === SymExtQET.RandomPPTConstructed
    @test result.verified
    @test result.state !== nothing
    @test ishermitian(result.state)
    @test tr(result.state) ≈ 1 atol = 8result.tolerance
    @test eigmin(Hermitian(result.state)) >= -result.tolerance
    transposed = SymExtQET.partial_transpose(result.state, result.dimensions; systems=(2,))
    @test eigmin(Hermitian(Matrix(transposed))) >= -result.tolerance
    @test result.numerical_ranks[1] <= result.requested_ranks[1]
    @test result.numerical_ranks[2] <= result.requested_ranks[2]
    return nothing
end

@testset "WP3 symmetric-extension core contracts" begin
    @testset "private Jacobi recurrence and inner coefficient" begin
        @test SymExtQET._jacobi_polynomial_coefficients(0 // 1, 0 // 1, 0) == [1 // 1]
        @test SymExtQET._jacobi_polynomial_coefficients(0 // 1, 0 // 1, 1) ==
            [1 // 1, 0 // 1]
        @test SymExtQET._jacobi_polynomial_coefficients(0 // 1, 0 // 1, 2) ==
            [3 // 2, 0 // 1, -1 // 2]
        @test SymExtQET._jacobi_polynomial_coefficients(1 // 1, 0 // 1, 1) ==
            [3 // 2, 1 // 2]
        roots = SymExtQET._jacobi_roots(0.0, 0.0, 3)
        @test roots.roots ≈ [-sqrt(3 / 5), 0, sqrt(3 / 5)] atol = 2e-14
        @test roots.residual <= 2e-14
        mixing = SymExtQET._symext_inner_mixing_parameter(2, 2, Float64)
        @test mixing.value ≈ 1 - inv(sqrt(3)) atol = 2e-14
        @test all(-1 .< mixing.roots .< 1)
        @test_throws ArgumentError SymExtQET._jacobi_polynomial_coefficients(0, 0, -1)
    end

    @testset "outer analytic and exact branches" begin
        mixed = Matrix{Float64}(I, 4, 4) / 4
        one_copy = SymExtQET.symmetric_extension(mixed; order=1, dims=(2, 2))
        @test one_copy.status === SymExtQET.SymmetricExtensionExactPresent
        @test one_copy.verdict === true
        @test one_copy.certificate_kind === :one_copy_identity_extension
        @test one_copy.extension == mixed
        @test one_copy.residuals.valid

        analytic_present = SymExtQET.symmetric_extension(mixed; dims=(2, 2))
        @test analytic_present.status === SymExtQET.SymmetricExtensionAnalyticPresent
        @test analytic_present.verdict === true
        @test analytic_present.certificate_kind === :two_qubit_symmetric_extension_theorem
        @test analytic_present.residuals.margin > analytic_present.tolerance

        bell = symext_bell_state()
        analytic_absent = SymExtQET.symmetric_extension(bell; dims=(2, 2))
        @test analytic_absent.status === SymExtQET.SymmetricExtensionAnalyticAbsent
        @test analytic_absent.verdict === false
        @test analytic_absent.residuals.margin < -analytic_absent.tolerance

        ppt_present = SymExtQET.symmetric_extension(
            mixed; dims=(2, 2), order=4, ppt=true, bosonic=true
        )
        @test ppt_present.status === SymExtQET.SymmetricExtensionAnalyticPresent
        @test ppt_present.verdict === true
        @test ppt_present.certificate_kind === :low_dimensional_ppt_separability_theorem

        ppt_absent = SymExtQET.symmetric_extension(bell; dims=(2, 2), order=3, ppt=true)
        @test ppt_absent.status === SymExtQET.SymmetricExtensionAnalyticAbsent
        @test ppt_absent.verdict === false
        @test ppt_absent.witness.separator_validated
        @test ppt_absent.witness.entanglement_witness
        @test real(dot(ppt_absent.witness.operator, bell)) ≈ -1 atol = 2e-14
        @test ppt_absent.witness.normalization_residual <= 2e-14

        high_dimensional_bell_vector = zeros(9)
        high_dimensional_bell_vector[[1, 5]] .= inv(sqrt(2.0))
        high_dimensional_bell =
            high_dimensional_bell_vector * transpose(high_dimensional_bell_vector)
        high_dimensional_npt = SymExtQET.symmetric_extension(
            high_dimensional_bell; dims=(3, 3), order=2, ppt=true
        )
        @test high_dimensional_npt.status === SymExtQET.SymmetricExtensionAnalyticAbsent
        @test high_dimensional_npt.verdict === false
        @test high_dimensional_npt.certificate_kind === :negative_partial_transpose

        robust_ppt = Matrix{Float64}(I, 8, 8) / 8
        one_copy_ppt = SymExtQET.symmetric_extension(
            robust_ppt; order=1, dims=(2, 4), ppt=true
        )
        @test one_copy_ppt.status === SymExtQET.SymmetricExtensionAnalyticPresent
        @test one_copy_ppt.verdict === true
        @test one_copy_ppt.certificate_kind === :one_copy_ppt_extension
        @test one_copy_ppt.extension == robust_ppt
        @test one_copy_ppt.residuals.valid

        missing = SymExtQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            prefer_analytic=false,
            allow_densify=true,
            backend=SymExtQET.NoOptimizationBackend(),
        )
        @test missing.status === SymExtQET.SymmetricExtensionBackendUnavailable
        @test missing.verdict === nothing
        @test missing.optimization_result.status ===
            SymExtQET.OptimizationBackendUnavailable
        @test missing.problem isa SymExtQET.SymmetricExtensionProblem

        limited = SymExtQET.symmetric_extension(mixed; dims=(2, 2), order=5, max_order=4)
        @test limited.status === SymExtQET.SymmetricExtensionResourceLimit
        @test limited.verdict === nothing
        @test occursin("max_order", limited.message)

        model_limited = SymExtQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            prefer_analytic=false,
            limits=SymExtQET.OptimizationLimits(max_variables=10),
        )
        @test model_limited.status === SymExtQET.SymmetricExtensionResourceLimit
        @test occursin("max_variables", model_limited.message)

        entry_limited = SymExtQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            ppt=true,
            prefer_analytic=false,
            limits=SymExtQET.OptimizationLimits(max_model_entries=300),
        )
        @test entry_limited.status === SymExtQET.SymmetricExtensionResourceLimit
        @test occursin("real-block PSD", entry_limited.message)

        stored_entry_limited = SymExtQET.symmetric_extension(
            mixed;
            dims=(2, 2),
            prefer_analytic=false,
            limits=SymExtQET.OptimizationLimits(max_model_entries=200),
        )
        @test stored_entry_limited.status === SymExtQET.SymmetricExtensionResourceLimit
        @test occursin("model stores", stored_entry_limited.message)
    end

    @testset "solver-neutral problem structure" begin
        mixed = Matrix{Float64}(I, 4, 4) / 4
        general = SymExtQET.symmetric_extension_problem(
            mixed; dims=(2, 2), order=2, allow_densify=true, bosonic=false, ppt=true
        )
        @test general.dimensions == (2, 2)
        @test general.ambient_dimension == 8
        @test general.variable_matrix_dimension == 8
        @test general.program.variable_count == 64
        @test general.marginal_equality_count == 16
        @test general.symmetry_equality_count == 64
        @test general.ppt_block_count == 2
        @test length(general.program.psd_constraints) == 3
        @test general.program.metadata.subsystem_order === :A_then_B_copies

        bosonic = SymExtQET.symmetric_extension_problem(
            mixed; dims=(2, 2), order=2, allow_densify=true, bosonic=true, ppt=true
        )
        @test bosonic.ambient_dimension == 8
        @test bosonic.variable_matrix_dimension == 6
        @test bosonic.program.variable_count == 36
        @test bosonic.symmetry_equality_count == 0
        @test size(bosonic.bosonic_basis) == (8, 6)
        @test Matrix(adjoint(bosonic.bosonic_basis) * bosonic.bosonic_basis) ≈
            Matrix{Float64}(I, 6, 6)

        inner = SymExtQET.symmetric_inner_extension_problem(
            mixed; dims=(2, 2), order=2, ppt=true, allow_densify=true
        )
        @test inner isa SymExtQET.SymmetricInnerExtensionProblem
        @test inner.mixing_parameter ≈ 1 - inv(sqrt(3)) atol = 2e-14
        @test inner.ambient_dimension == 8
        @test inner.variable_matrix_dimension == 6
        @test inner.marginal_equality_count == 16
        @test inner.ppt_block_count == 1
        @test length(inner.program.psd_constraints) == 2
        @test inner.program.metadata.formulation === :navascues_owari_plenio_inner_extension
        @test inner.metadata.preflight.psd_blocks == 2
        @test inner.metadata.preflight.ppt_blocks == 1

        mixed32 = Matrix{Float32}(I, 4, 4) / 4
        float32_problem = SymExtQET.symmetric_extension_problem(
            mixed32; dims=(2, 2), order=2, bosonic=true, allow_densify=true
        )
        @test eltype(float32_problem.state) === Float32
        @test float32_problem.tolerance isa Float32
        @test eltype(float32_problem.program.objective.coefficients) === Float32
        @test eltype(first(float32_problem.program.psd_constraints).constant) === ComplexF32

        missing_inner = SymExtQET.symmetric_inner_extension(
            mixed; dims=(2, 2), allow_densify=true
        )
        @test missing_inner.status === SymExtQET.SymmetricExtensionBackendUnavailable
        @test missing_inner.verdict === nothing
        @test missing_inner.hierarchy === :inner

        zero_inner = SymExtQET.symmetric_inner_extension(
            zeros(ComplexF64, 4, 4); dims=(2, 2), order=3, ppt=true
        )
        @test zero_inner.status === SymExtQET.SymmetricExtensionExactPresent
        @test zero_inner.verdict === true
        @test size(zero_inner.extension) == (16, 16)
        @test iszero(norm(zero_inner.extension))
        @test zero_inner.residuals.valid

        zero_limited = SymExtQET.symmetric_inner_extension(
            zeros(4, 4); dims=(2, 2), order=5, max_dense_entries=100
        )
        @test zero_limited.status === SymExtQET.SymmetricExtensionResourceLimit
        @test zero_limited.extension === nothing
        @test occursin("explicit zero extension", zero_limited.message)
    end

    @testset "validation and malformed inputs" begin
        mixed = Matrix{Float64}(I, 4, 4) / 4
        @test_throws DimensionMismatch SymExtQET.symmetric_extension(
            ones(2, 3); dims=(2, 2)
        )
        @test_throws DimensionMismatch SymExtQET.symmetric_extension(mixed; dims=(2, 3))
        @test_throws ArgumentError SymExtQET.symmetric_extension(
            ComplexF64[1 im; im 1]; dims=(1, 2)
        )
        @test_throws ArgumentError SymExtQET.symmetric_extension(
            Diagonal([1.0, -0.1]); dims=(1, 2)
        )
        @test_throws ArgumentError SymExtQET.symmetric_extension(
            mixed; order=0, dims=(2, 2)
        )
        @test_throws ArgumentError SymExtQET.symmetric_inner_extension(
            mixed; order=1, dims=(2, 2)
        )
        @test_throws ArgumentError SymExtQET.symmetric_extension(sparse(mixed); dims=(2, 2))
        sparse_allowed = SymExtQET.symmetric_extension(
            sparse(mixed); dims=(2, 2), allow_densify=true
        )
        @test sparse_allowed.verdict === true
        @test_throws ArgumentError SymExtQET.symmetric_extension(
            fill(NaN, 4, 4); dims=(2, 2)
        )
        @test_throws ArgumentError SymExtQET.symmetric_extension(
            Matrix{Int}(I, 4, 4); dims=(2, 2)
        )

        negative_boundary = Diagonal([1.0, 1.0, 1.0, -1.0e-10])
        outer_boundary = SymExtQET.symmetric_extension(
            negative_boundary; dims=(2, 2), order=1
        )
        @test outer_boundary.status === SymExtQET.SymmetricExtensionNumericalBoundary
        @test outer_boundary.verdict === nothing
        @test outer_boundary.residuals.psd_violation == 1.0e-10
        inner_boundary = SymExtQET.symmetric_inner_extension(
            negative_boundary; dims=(2, 2), order=2
        )
        @test inner_boundary.status === SymExtQET.SymmetricExtensionNumericalBoundary
        @test inner_boundary.verdict === nothing
        @test_throws ArgumentError SymExtQET.symmetric_extension_problem(
            negative_boundary; dims=(2, 2), order=2
        )
        @test_throws ArgumentError SymExtQET.symmetric_inner_extension_problem(
            negative_boundary; dims=(2, 2), order=2
        )
    end

    @testset "bounded explicit-RNG PPT construction" begin
        full = SymExtQET.random_ppt_state(Xoshiro(7101), (2, 3))
        assert_ppt_state(full)
        @test full.construction === :shifted_induced
        @test full.requested_ranks == (6, 6)
        @test full.normalized_by_construction
        @test length(full.convergence_history) == 1

        low_rank = SymExtQET.random_ppt_state(Xoshiro(7102), (2, 3); ranks=(3, 4))
        assert_ppt_state(low_rank)
        @test low_rank.construction === :separable_mixture
        @test low_rank.numerical_ranks[1] <= 3
        @test low_rank.numerical_ranks[2] <= 4

        scalar = SymExtQET.random_ppt_state(Xoshiro(7103), 2; ranks=2, real=true, T=Float32)
        assert_ppt_state(scalar)
        @test scalar.dimensions == (2, 2)
        @test eltype(scalar.state) === Float32
        @test scalar.real_output

        first_draw = SymExtQET.random_ppt_state(Xoshiro(7104), (2, 2); ranks=3)
        second_draw = SymExtQET.random_ppt_state(Xoshiro(7104), (2, 2); ranks=3)
        @test first_draw.state == second_draw.state

        global_before = copy(Random.default_rng())
        SymExtQET.random_ppt_state(Xoshiro(7105), (2, 2); ranks=3)
        @test rand(global_before, UInt64) == rand(Random.default_rng(), UInt64)

        limited_rng = Xoshiro(7106)
        limited_reference = copy(limited_rng)
        limited = SymExtQET.random_ppt_state(limited_rng, (3, 3); max_dense_entries=10)
        @test limited.status === SymExtQET.RandomPPTResourceLimit
        @test limited.state === nothing
        @test limited.candidate === nothing
        @test rand(limited_rng, UInt64) == rand(limited_reference, UInt64)

        iteration_rng = Xoshiro(7107)
        iteration_reference = copy(iteration_rng)
        iteration_limited = SymExtQET.random_ppt_state(
            iteration_rng, (2, 2); max_iterations=0
        )
        @test iteration_limited.status === SymExtQET.RandomPPTIterationLimit
        @test !iteration_limited.verified
        @test rand(iteration_rng, UInt64) == rand(iteration_reference, UInt64)

        work_limited = SymExtQET.random_ppt_state(Xoshiro(7108), (2, 2); max_work=1)
        @test work_limited.status === SymExtQET.RandomPPTResourceLimit
        @test occursin("max_work", work_limited.message)

        @test_throws ArgumentError SymExtQET.random_ppt_state(Xoshiro(1), (2, 3); ranks=7)
        @test_throws ArgumentError SymExtQET.random_ppt_state(
            Xoshiro(1), (2, 3); ranks=(2, 3, 4)
        )
        @test_throws ArgumentError SymExtQET.random_ppt_state(
            Xoshiro(1), (2, 3); ranks=3, construction=:shifted_induced
        )
        @test_throws ArgumentError SymExtQET.random_ppt_state(
            Xoshiro(1), (2, 3); construction=:unknown
        )
        @test_throws ArgumentError SymExtQET.random_ppt_state(
            Xoshiro(1), (2, 3); T=BigFloat
        )
        @test_throws ArgumentError SymExtQET.random_ppt_state(Xoshiro(1), (2, 0))
    end
end
