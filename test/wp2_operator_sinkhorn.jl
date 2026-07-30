using LinearAlgebra
using Random
using SparseArrays

const QETOperatorSinkhorn = QuantumEntanglementTools
const CompatOperatorSinkhorn = QuantumEntanglementTools.MATLABCompat

function _wp2_sinkhorn_reconstruction(result, input)
    filter = tensor_product(result.local_filters...)
    return filter * input * adjoint(filter)
end

function _wp2_sinkhorn_marginal(operator, dims, kept_subsystem)
    traced = Tuple(
        subsystem for subsystem in eachindex(dims) if subsystem != kept_subsystem
    )
    return partial_trace(operator, dims; trace_out=traced)
end

@testset "WP2 OperatorSinkhorn bounded normalization" begin
    @testset "balanced and low-rank fixed points" begin
        balanced = 2Matrix{Float64}(I, 6, 6) / 6
        original = copy(balanced)
        result = QETOperatorSinkhorn.operator_sinkhorn(balanced, (2, 3))

        @test result isa QETOperatorSinkhorn.OperatorSinkhornResult
        @test result.status === :converged
        @test result.converged
        @test result.iterations == 0
        @test length(result.residual_history) == 1
        @test only(result.residual_history) <= result.convergence_threshold
        @test result.final_marginal_residuals ≈ zeros(2) atol = 2e-15
        @test result.scaled_operator ≈ balanced atol = 2e-15
        @test tr(result.scaled_operator) ≈ 2
        @test _wp2_sinkhorn_reconstruction(result, balanced) ≈ result.scaled_operator atol =
            2e-15
        @test result.local_filters[1] == Matrix{Float64}(I, 2, 2)
        @test result.local_filters[2] == Matrix{Float64}(I, 3, 3)
        @test result.left_filter === result.local_filters[1]
        @test result.right_filter === result.local_filters[2]
        @test result.filter_condition_numbers == ones(2)
        @test result.maximum_filter_condition == 1
        @test result.failed_subsystem === nothing
        @test result.work_used > 0
        @test balanced == original
        @test occursin("status=converged", sprint(show, result))
        absolute_only = QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 3); atol=1e-6, rtol=0
        )
        @test absolute_only.convergence_threshold == 1e-6

        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * adjoint(bell)
        low_rank = QETOperatorSinkhorn.operator_sinkhorn(bell_density, (2, 2))
        @test low_rank.status === :converged
        @test low_rank.iterations == 0
        @test low_rank.scaled_operator ≈ bell_density atol = 2e-15
        @test all(value -> value > 0, low_rank.minimum_marginal_eigenvalues)
    end

    @testset "bipartite convergence, filters, and complex input" begin
        left_density = Diagonal([0.8, 0.2])
        right_density = Diagonal([0.6, 0.3, 0.1])
        product_density = tensor_product(left_density, right_density)
        original = copy(product_density)
        product_result = QETOperatorSinkhorn.operator_sinkhorn(product_density, (2, 3))

        @test product_result.status === :converged
        @test product_result.iterations == 1
        @test length(product_result.residual_history) == product_result.iterations + 1
        @test first(product_result.residual_history) > product_result.convergence_threshold
        @test last(product_result.residual_history) <= product_result.convergence_threshold
        @test product_result.scaled_operator ≈ Matrix{Float64}(I, 6, 6) / 6 atol = 2e-14
        @test _wp2_sinkhorn_reconstruction(product_result, product_density) ≈
            product_result.scaled_operator atol = 2e-14
        @test _wp2_sinkhorn_marginal(product_result.scaled_operator, (2, 3), 1) ≈
            Matrix{Float64}(I, 2, 2) / 2 atol = 2e-14
        @test _wp2_sinkhorn_marginal(product_result.scaled_operator, (2, 3), 2) ≈
            Matrix{Float64}(I, 3, 3) / 3 atol = 2e-14
        @test product_result.filter_condition_numbers ≈ [2.0, sqrt(6)]
        @test product_result.maximum_filter_condition ≈ sqrt(6)
        @test all(isfinite, product_result.minimum_marginal_eigenvalues)
        @test product_result.maximum_marginal_hermiticity_residual <= 64eps(Float64)
        @test product_density == original

        rng = MersenneTwister(0x53494e4b484f524e)
        factor = randn(rng, ComplexF64, 6, 6)
        complex_density = factor * adjoint(factor) + Matrix{ComplexF64}(I, 6, 6)
        complex_density ./= real(tr(complex_density))
        complex_original = copy(complex_density)
        complex_result = QETOperatorSinkhorn.operator_sinkhorn(
            complex_density, (2, 3); max_iterations=200
        )
        @test complex_result.status === :converged
        @test 0 < complex_result.iterations <= 200
        @test last(complex_result.residual_history) <= complex_result.convergence_threshold
        @test _wp2_sinkhorn_reconstruction(complex_result, complex_density) ≈
            complex_result.scaled_operator atol = 3e-13
        @test all(factor -> eltype(factor) == ComplexF64, complex_result.local_filters)
        for subsystem in 1:2
            target = Matrix{ComplexF64}(I, (2, 3)[subsystem], (2, 3)[subsystem])
            target ./= (2, 3)[subsystem]
            @test _wp2_sinkhorn_marginal(
                complex_result.scaled_operator, (2, 3), subsystem
            ) ≈ target atol = 5e-8
        end
        @test complex_density == complex_original
    end

    @testset "multipartite branch and preserved input trace" begin
        factors = (Diagonal([0.7, 0.3]), Diagonal([0.6, 0.4]), Diagonal([0.8, 0.2]))
        input = 3tensor_product(factors...)
        result = QETOperatorSinkhorn.operator_sinkhorn(input, (2, 2, 2))
        @test result.status === :converged
        @test result.iterations == 1
        @test length(result.local_filters) == 3
        @test result.left_filter === result.local_filters[1]
        @test result.right_filter === result.local_filters[2]
        @test result.input_trace ≈ 3
        @test tr(result.scaled_operator) ≈ 3 atol = 2e-14
        @test result.scaled_operator ≈ 3Matrix{Float64}(I, 8, 8) / 8 atol = 3e-14
        @test _wp2_sinkhorn_reconstruction(result, input) ≈ result.scaled_operator atol =
            3e-14
        for subsystem in 1:3
            @test _wp2_sinkhorn_marginal(
                result.scaled_operator / result.input_trace, (2, 2, 2), subsystem
            ) ≈ Matrix{Float64}(I, 2, 2) / 2 atol = 3e-14
        end
    end

    @testset "bounded and controlled failure statuses" begin
        singular = zeros(4, 4)
        singular[1, 1] = 1
        singular_result = QETOperatorSinkhorn.operator_sinkhorn(singular, (2, 2))
        @test singular_result.status === :singular_marginal
        @test !singular_result.converged
        @test singular_result.iterations == 0
        @test singular_result.failed_subsystem == 1
        @test singular_result.minimum_marginal_eigenvalues == zeros(2)
        @test _wp2_sinkhorn_reconstruction(singular_result, singular) ≈
            singular_result.scaled_operator

        ill_conditioned = tensor_product(
            Diagonal([1 - 1e-12, 1e-12]), Matrix{Float64}(I, 2, 2) / 2
        )
        ill_result = QETOperatorSinkhorn.operator_sinkhorn(
            ill_conditioned, (2, 2); max_condition_number=100
        )
        @test ill_result.status === :ill_conditioned
        @test !ill_result.converged
        @test ill_result.failed_subsystem == 1
        @test ill_result.maximum_filter_condition > ill_result.condition_limit
        @test ill_result.filter_condition_numbers == ones(2)

        unbalanced = Diagonal([0.4, 0.3, 0.2, 0.1])
        bounded = QETOperatorSinkhorn.operator_sinkhorn(
            unbalanced, (2, 2); max_iterations=0
        )
        @test bounded.status === :max_iterations
        @test bounded.iterations == 0
        @test length(bounded.residual_history) == 1
        @test last(bounded.residual_history) > bounded.convergence_threshold

        work_limited = QETOperatorSinkhorn.operator_sinkhorn(
            unbalanced, (2, 2); max_work=112
        )
        @test work_limited.status === :work_limit
        @test work_limited.iterations == 0
        @test work_limited.work_used == 112
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            unbalanced, (2, 2); max_work=111
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            unbalanced, (2, 2); max_entries=15
        )
    end

    @testset "MATLABCompat successful positional outputs only" begin
        balanced = Matrix{Float64}(I, 4, 4) / 4
        inferred = CompatOperatorSinkhorn.OperatorSinkhorn(balanced)
        scalar = CompatOperatorSinkhorn.OperatorSinkhorn(balanced, 2)
        singleton = CompatOperatorSinkhorn.OperatorSinkhorn(balanced, (2,))
        vector = CompatOperatorSinkhorn.OperatorSinkhorn(balanced, [2, 2])
        @test keys(inferred) == (:sigma, :filters)
        @test length(inferred) == 2
        sigma, filters = inferred
        @test sigma == inferred.sigma
        @test filters === inferred.filters
        @test inferred.sigma ≈ scalar.sigma
        @test inferred.sigma ≈ singleton.sigma
        @test inferred.sigma ≈ vector.sigma
        @test length(filters) == 2
        explicit_tolerance = CompatOperatorSinkhorn.OperatorSinkhorn(
            tensor_product(Diagonal([0.8, 0.2]), Diagonal([0.6, 0.3, 0.1])), (2, 3), 1e-10
        )
        @test explicit_tolerance.sigma ≈ Matrix{Float64}(I, 6, 6) / 6 atol = 2e-14
        explicit_filter = tensor_product(explicit_tolerance.filters...)
        @test explicit_tolerance.sigma ≈
            explicit_filter *
              tensor_product(Diagonal([0.8, 0.2]), Diagonal([0.6, 0.3, 0.1])) *
              adjoint(explicit_filter) atol = 2e-14

        multipartite_density = tensor_product(
            Diagonal([0.7, 0.3]), Diagonal([0.6, 0.4]), Diagonal([0.8, 0.2])
        )
        multipartite = CompatOperatorSinkhorn.OperatorSinkhorn(
            multipartite_density, (2, 2, 2)
        )
        @test length(multipartite.filters) == 3
        @test multipartite.sigma ≈ Matrix{Float64}(I, 8, 8) / 8 atol = 3e-14

        singular = zeros(4, 4)
        singular[1, 1] = 1
        @test_throws DomainError CompatOperatorSinkhorn.OperatorSinkhorn(singular, (2, 2))
        @test_throws DomainError CompatOperatorSinkhorn.OperatorSinkhorn(
            Diagonal([0.4, 0.3, 0.2, 0.1]), (2, 2); max_iterations=0
        )
        @test_throws ArgumentError CompatOperatorSinkhorn.OperatorSinkhorn(
            balanced, (2, 2), 0
        )
        @test_throws ArgumentError CompatOperatorSinkhorn.OperatorSinkhorn(
            balanced, (2, 2), Inf
        )
        @test_throws ArgumentError CompatOperatorSinkhorn.OperatorSinkhorn(
            balanced, [2 2; 2 2]
        )
    end

    @testset "sparse, precision, and invalid input policies" begin
        float32_input = tensor_product(
            Diagonal(Float32[0.75, 0.25]), Diagonal(Float32[0.6, 0.4])
        )
        float32_result = QETOperatorSinkhorn.operator_sinkhorn(float32_input, (2, 2))
        @test float32_result.status === :converged
        @test eltype(float32_result.scaled_operator) == Float32
        @test all(factor -> eltype(factor) == Float32, float32_result.local_filters)
        @test eltype(float32_result.residual_history) == Float32

        complex32_factor = ComplexF32[
            1 0 0 0
            1im 2 0 0
            0 1 2 0
            1 0 1im 3
        ]
        complex32_input =
            complex32_factor * adjoint(complex32_factor) + Matrix{ComplexF32}(I, 4, 4)
        complex32_input ./= real(tr(complex32_input))
        complex32_result = QETOperatorSinkhorn.operator_sinkhorn(complex32_input, (2, 2))
        @test complex32_result.status === :converged
        @test eltype(complex32_result.scaled_operator) == ComplexF32
        @test all(factor -> eltype(factor) == ComplexF32, complex32_result.local_filters)

        sparse_input = sparse(float32_input)
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            sparse_input, (2, 2)
        )
        sparse_result = QETOperatorSinkhorn.operator_sinkhorn(
            sparse_input, (2, 2); allow_densify=true
        )
        @test sparse_result.status === :converged
        @test sparse_result.scaled_operator isa Matrix{Float32}
        @test_throws ArgumentError CompatOperatorSinkhorn.OperatorSinkhorn(
            sparse_input, (2, 2)
        )
        sparse_compatibility = CompatOperatorSinkhorn.OperatorSinkhorn(
            sparse_input, (2, 2); allow_densify=true
        )
        @test sparse_compatibility.sigma isa Matrix{Float32}

        balanced = Matrix{Float64}(I, 4, 4) / 4
        nonhermitian = copy(balanced)
        nonhermitian[1, 2] = 1e-14
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            nonhermitian, (2, 2)
        )
        @test_throws DomainError QETOperatorSinkhorn.operator_sinkhorn(
            Diagonal([0.4, 0.3, 0.4, -0.1]), (2, 2)
        )
        @test_throws DomainError QETOperatorSinkhorn.operator_sinkhorn(zeros(4, 4), (2, 2))
        @test_throws DimensionMismatch QETOperatorSinkhorn.operator_sinkhorn(
            ones(4, 3), (2, 2)
        )
        @test_throws DimensionMismatch QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 3)
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(balanced, (4,))
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            Matrix{BigFloat}(balanced), (2, 2)
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            Matrix{Int}(I, 4, 4), (2, 2)
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); atol=-1
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); rtol=Inf
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); max_iterations=-1
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); max_iterations=true
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); max_condition_number=0
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); max_entries=true
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            balanced, (2, 2); max_work=0
        )
        @test_throws ArgumentError QETOperatorSinkhorn.operator_sinkhorn(
            Diagonal(fill(floatmax(Float64), 4)), (2, 2)
        )
    end
end
