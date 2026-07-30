using LinearAlgebra
using Random
using SparseArrays

const QETFilterNormalForm = QuantumEntanglementTools
const CompatFilterNormalForm = QuantumEntanglementTools.MATLABCompat

function _wp2_filter_reconstruction(result)
    identity_matrix = Matrix{eltype(result.filtered_operator)}(
        I, result.input_dimension, result.input_dimension
    )
    tensor_component = tensor_sum(
        result.left_operators, result.right_operators; weights=result.coefficients
    )
    return (result.normal_form_trace * identity_matrix + tensor_component) /
           result.input_dimension
end

function _wp2_filter_action(result, rho)
    full_filter = tensor_product(result.local_filters...)
    return full_filter * rho * adjoint(full_filter)
end

function _wp2_filter_marginal(operator, dims, kept)
    traced = Tuple(index for index in eachindex(dims) if index != kept)
    return partial_trace(operator, dims; trace_out=traced)
end

@testset "WP2 FilterNormalForm bounded decomposition" begin
    @testset "successful full-rank form and both identities" begin
        bell = [1.0, 0, 0, 1] / sqrt(2)
        rho = 0.8 * (bell * transpose(bell)) + 0.2 * Matrix{Float64}(I, 4, 4) / 4
        original = copy(rho)
        result = QETFilterNormalForm.filter_normal_form(rho, (2, 2))

        @test result isa QETFilterNormalForm.FilterNormalFormResult
        @test result.status === :converged
        @test result.converged
        @test result.sinkhorn_result.status === :converged
        @test result.sinkhorn_result.iterations == 0
        @test result.input_dimension == 4
        @test result.input_numerical_rank == 4
        @test result.input_numerically_full_rank
        @test result.input_minimum_eigenvalue ≈ 0.05
        @test result.input_maximum_eigenvalue ≈ 0.85
        @test result.normal_form_trace ≈ 1
        @test length(result.coefficients) == 4
        @test result.coefficient_numerical_rank == 3
        @test all(result.coefficients[1:3] .≈ 1.6)
        @test result.coefficients[4] <= result.coefficient_threshold
        @test result.left_filter === result.local_filters[1]
        @test result.right_filter === result.local_filters[2]
        @test _wp2_filter_action(result, rho) ≈ result.filtered_operator atol = 2e-15
        @test _wp2_filter_reconstruction(result) ≈ result.filtered_operator atol = 2e-15
        @test result.filter_identity_residual <= result.filter_identity_tolerance
        @test result.normal_form_reconstruction_residual <=
            result.normal_form_reconstruction_tolerance
        @test result.decomposition_hermiticity_residual <=
            result.decomposition_hermiticity_tolerance
        @test result.work_used > result.sinkhorn_result.work_used
        @test result.required_next_stage_work == 0
        @test result.failed_subsystem === nothing
        @test rho == original
        @test occursin("status=converged", sprint(show, result))

        for factors in (result.left_operators, result.right_operators)
            @test all(ishermitian, factors)
            gram = [
                tr(adjoint(factors[row]) * factors[column]) for
                row in eachindex(factors), column in eachindex(factors)
            ]
            @test gram ≈ Matrix{Float64}(I, length(factors), length(factors)) atol = 3e-15
        end
    end

    @testset "rectangular complex balancing and trace gauge" begin
        rng = MersenneTwister(0x464e46)
        factor = randn(rng, ComplexF64, 6, 6)
        rho =
            2 * (factor * adjoint(factor) + Matrix{ComplexF64}(I, 6, 6)) /
            real(tr(factor * adjoint(factor) + Matrix{ComplexF64}(I, 6, 6)))
        original = copy(rho)
        result = QETFilterNormalForm.filter_normal_form(rho, (2, 3); max_iterations=300)

        @test result.status === :converged
        @test 0 < result.sinkhorn_result.iterations <= 300
        @test result.normal_form_trace ≈ 2 atol = 2e-14
        @test length(result.coefficients) == 4
        @test length(result.left_operators) == 4
        @test length(result.right_operators) == 4
        @test all(factor -> size(factor) == (2, 2), result.left_operators)
        @test all(factor -> size(factor) == (3, 3), result.right_operators)
        @test all(ishermitian, result.left_operators)
        @test all(ishermitian, result.right_operators)
        @test _wp2_filter_action(result, rho) ≈ result.filtered_operator atol = 8e-13
        @test _wp2_filter_reconstruction(result) ≈ result.filtered_operator atol = 8e-13
        @test result.filter_identity_residual <= result.filter_identity_tolerance
        @test result.normal_form_reconstruction_residual <=
            result.normal_form_reconstruction_tolerance
        @test result.decomposition_hermiticity_residual > 0
        @test result.centered_operator !==
            (result.centered_operator + adjoint(result.centered_operator)) / 2
        for (subsystem, dimension) in enumerate((2, 3))
            target = result.normal_form_trace * Matrix{ComplexF64}(I, dimension, dimension)
            target /= dimension
            @test _wp2_filter_marginal(result.filtered_operator, (2, 3), subsystem) ≈ target atol =
                2e-7
        end
        @test rho == original
    end

    @testset "rank deficiency is diagnostic when a form exists" begin
        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        rho = bell * adjoint(bell)
        result = QETFilterNormalForm.filter_normal_form(rho, (2, 2))

        @test result.status === :converged
        @test result.input_numerical_rank == 1
        @test !result.input_numerically_full_rank
        @test result.input_minimum_eigenvalue == 0
        @test result.sinkhorn_result.iterations == 0
        @test all(>(0), result.sinkhorn_result.minimum_marginal_eigenvalues)
        @test _wp2_filter_action(result, rho) ≈ result.filtered_operator atol = 2e-15
        @test _wp2_filter_reconstruction(result) ≈ result.filtered_operator atol = 2e-15

        strict_rank = QETFilterNormalForm.filter_normal_form(
            rho, (2, 2); rank_atol=0.6, rank_rtol=0
        )
        @test strict_rank.status === :converged
        @test strict_rank.input_numerical_rank == 1
        @test strict_rank.input_rank_threshold == 0.6
    end

    @testset "all bounded failure statuses" begin
        singular = zeros(4, 4)
        singular[1, 1] = 1
        singular_result = QETFilterNormalForm.filter_normal_form(singular, (2, 2))
        @test singular_result.status === :singular_marginal
        @test !singular_result.converged
        @test singular_result.input_numerical_rank == 1
        @test singular_result.failed_subsystem == 1
        @test singular_result.coefficients === nothing
        @test singular_result.sinkhorn_result.status === :singular_marginal
        @test _wp2_filter_action(singular_result, singular) ≈
            singular_result.filtered_operator

        ill_conditioned = tensor_product(
            Diagonal([1 - 1e-12, 1e-12]), Matrix{Float64}(I, 2, 2) / 2
        )
        ill_result = QETFilterNormalForm.filter_normal_form(
            ill_conditioned, (2, 2); max_condition_number=100
        )
        @test ill_result.status === :ill_conditioned
        @test !ill_result.converged
        @test ill_result.failed_subsystem == 1
        @test ill_result.sinkhorn_result.maximum_filter_condition >
            ill_result.sinkhorn_result.condition_limit

        unbalanced = Diagonal([0.4, 0.3, 0.2, 0.1])
        bounded = QETFilterNormalForm.filter_normal_form(
            unbalanced, (2, 2); max_iterations=0
        )
        @test bounded.status === :max_iterations
        @test bounded.sinkhorn_result.iterations == 0
        @test last(bounded.sinkhorn_result.residual_history) >
            bounded.sinkhorn_result.convergence_threshold

        sweep_limited = QETFilterNormalForm.filter_normal_form(
            unbalanced, (2, 2); max_work=176
        )
        @test sweep_limited.status === :work_limit
        @test sweep_limited.sinkhorn_result.status === :work_limit
        @test sweep_limited.work_used == 176
        @test sweep_limited.required_next_stage_work > 0

        preflight_limited = QETFilterNormalForm.filter_normal_form(
            Matrix{Float64}(I, 4, 4) / 4, (2, 2); max_work=175
        )
        @test preflight_limited.status === :work_limit
        @test preflight_limited.sinkhorn_result === nothing
        @test preflight_limited.work_used == 0
        @test preflight_limited.required_next_stage_work == 176

        postprocess_limited = QETFilterNormalForm.filter_normal_form(
            Matrix{Float64}(I, 4, 4) / 4, (2, 2); max_work=463
        )
        @test postprocess_limited.status === :work_limit
        @test postprocess_limited.sinkhorn_result.status === :converged
        @test postprocess_limited.filtered_operator !== nothing
        @test postprocess_limited.coefficients === nothing
        @test postprocess_limited.required_next_stage_work ==
            postprocess_limited.estimated_postprocessing_work

        bell = [1.0, 0, 0, 1] / sqrt(2)
        noisy_bell = 0.8 * (bell * transpose(bell)) + 0.2 * Matrix{Float64}(I, 4, 4) / 4
        numerical = QETFilterNormalForm.filter_normal_form(
            noisy_bell, (2, 2); verification_atol=0, verification_rtol=0
        )
        @test numerical.status === :numerical_failure
        @test !numerical.converged
        @test numerical.coefficients !== nothing
        @test numerical.normal_form_reconstruction_residual >
            numerical.normal_form_reconstruction_tolerance
        @test occursin("reconstruction", numerical.message)
    end

    @testset "compatibility dimensions, tolerance correction, and five outputs" begin
        bell = [1.0, 0, 0, 1] / sqrt(2)
        rho = 0.8 * (bell * transpose(bell)) + 0.2 * Matrix{Float64}(I, 4, 4) / 4
        inferred = CompatFilterNormalForm.FilterNormalForm(rho)
        scalar = CompatFilterNormalForm.FilterNormalForm(rho, 2)
        singleton = CompatFilterNormalForm.FilterNormalForm(rho, (2,))
        vector = CompatFilterNormalForm.FilterNormalForm(rho, [2, 2])

        @test keys(inferred) == (:xi, :GA, :GB, :FA, :FB)
        @test length(inferred) == 5
        xi, GA, GB, FA, FB = inferred
        @test xi === inferred.xi
        @test GA === inferred.GA
        @test GB === inferred.GB
        @test FA === inferred.FA
        @test FB === inferred.FB
        @test inferred.xi ≈ scalar.xi
        @test inferred.xi ≈ singleton.xi
        @test inferred.xi ≈ vector.xi
        @test length(xi) == 3
        filtered = tensor_product(FA, FB) * rho * adjoint(tensor_product(FA, FB))
        @test filtered ≈ (Matrix{Float64}(I, 4, 4) + tensor_sum(GA, GB; weights=xi)) / 4 atol =
            2e-15

        nearly_balanced = Diagonal([0.26, 0.24, 0.25, 0.25])
        loose = CompatFilterNormalForm.FilterNormalForm(nearly_balanced, (2, 2), 0.1)
        tight = CompatFilterNormalForm.FilterNormalForm(nearly_balanced, (2, 2), 1e-12)
        @test loose.FA == Matrix{Float64}(I, 2, 2)
        @test loose.FB == Matrix{Float64}(I, 2, 2)
        @test norm(tight.FA - I) > 0
        @test norm(tight.FB - I) > 0

        singular = zeros(4, 4)
        singular[1, 1] = 1
        error = try
            CompatFilterNormalForm.FilterNormalForm(singular, (2, 2))
            nothing
        catch caught
            caught
        end
        @test error isa DomainError
        @test error.val isa QETFilterNormalForm.FilterNormalFormResult
        @test error.val.status === :singular_marginal
        @test_throws DomainError CompatFilterNormalForm.FilterNormalForm(
            Diagonal([0.4, 0.3, 0.2, 0.1]), (2, 2); max_iterations=0
        )
        @test_throws ArgumentError CompatFilterNormalForm.FilterNormalForm(rho, (2, 2), 0)
        @test_throws DimensionMismatch CompatFilterNormalForm.FilterNormalForm(
            rho, [2, 2, 1]
        )
        @test_throws ArgumentError CompatFilterNormalForm.FilterNormalForm(rho, [2 2; 2 2])
    end

    @testset "types, sparse policy, guards, and invalid input" begin
        bell32 = ComplexF32[1, 0, 0, 1] / sqrt(Float32(2))
        rho32 =
            Float32(0.7) * (bell32 * adjoint(bell32)) +
            Float32(0.3) * Matrix{ComplexF32}(I, 4, 4) / 4
        result32 = QETFilterNormalForm.filter_normal_form(rho32, (2, 2))
        @test result32.status === :converged
        @test eltype(result32.filtered_operator) == ComplexF32
        @test eltype(result32.centered_operator) == ComplexF32
        @test eltype(result32.coefficients) == Float32
        @test all(factor -> eltype(factor) == ComplexF32, result32.left_operators)
        @test typeof(result32.input_rank_threshold) == Float32
        @test typeof(result32.coefficient_threshold) == Float32

        sparse_rho = sparse(Matrix{Float64}(I, 4, 4) / 4)
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            sparse_rho, (2, 2)
        )
        sparse_result = QETFilterNormalForm.filter_normal_form(
            sparse_rho, (2, 2); allow_densify=true
        )
        @test sparse_result.status === :converged
        @test sparse_result.filtered_operator isa Matrix{Float64}
        sparse_compat = CompatFilterNormalForm.FilterNormalForm(
            sparse_rho, (2, 2); allow_densify=true
        )
        @test keys(sparse_compat) == (:xi, :GA, :GB, :FA, :FB)

        balanced = Matrix{Float64}(I, 4, 4) / 4
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_entries=95
        )
        @test QETFilterNormalForm.filter_normal_form(balanced, (2, 2); max_entries=96).status ===
            :converged
        @test QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_entries=nothing, max_work=nothing
        ).status === :converged

        @test_throws DimensionMismatch QETFilterNormalForm.filter_normal_form(
            zeros(2, 3), (2, 2)
        )
        @test_throws DimensionMismatch QETFilterNormalForm.filter_normal_form(
            balanced, (2, 3)
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(balanced, (4,))
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            Matrix{Float64}(I, 8, 8) / 8, (2, 2, 2)
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced + ComplexF64[0 im 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], (2, 2)
        )
        indefinite = Diagonal([0.4, 0.3, 0.2, -0.1])
        @test_throws DomainError QETFilterNormalForm.filter_normal_form(indefinite, (2, 2))
        @test_throws DomainError QETFilterNormalForm.filter_normal_form(zeros(4, 4), (2, 2))
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            fill(NaN, 4, 4), (2, 2)
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            Matrix{BigFloat}(I, 4, 4) / 4, (2, 2)
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            Matrix{Int}(I, 4, 4), (2, 2)
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_iterations=-1
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_condition_number=0
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); balance_rtol=-1
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); rank_atol=Inf
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); coefficient_rtol=-1
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); verification_atol=NaN
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_entries=0
        )
        @test_throws ArgumentError QETFilterNormalForm.filter_normal_form(
            balanced, (2, 2); max_work=0
        )
    end
end
