using LinearAlgebra
using Random
using Test

const QETRandomSuperoperator = QuantumEntanglementTools

if !isdefined(QETRandomSuperoperator, :RandomSuperoperatorResult)
    Base.include(
        QETRandomSuperoperator,
        joinpath(@__DIR__, "..", "src", "channels", "random_superoperator.jl"),
    )
end

struct _RandomSuperoperatorZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_RandomSuperoperatorZeroBasedVector) = size(vector.storage)
function Base.axes(vector::_RandomSuperoperatorZeroBasedVector)
    return (0:(length(vector.storage) - 1),)
end
Base.IndexStyle(::Type{<:_RandomSuperoperatorZeroBasedVector}) = IndexLinear()
function Base.getindex(vector::_RandomSuperoperatorZeroBasedVector, index::Int)
    return vector.storage[index + 1]
end

function _random_superoperator_marginals(result)
    factors = QETRandomSuperoperator.kraus_operators(result.kraus_representation)
    input_dimension = result.input_dimension
    output_dimension = result.output_dimension
    input_marginal = zeros(eltype(first(factors)), input_dimension, input_dimension)
    output_marginal = zeros(eltype(first(factors)), output_dimension, output_dimension)
    for factor in factors
        input_marginal += adjoint(factor) * factor
        output_marginal += factor * adjoint(factor)
    end
    return input_marginal, output_marginal
end

function _random_superoperator_choi(map)
    return QETRandomSuperoperator.choi_matrix(
        QETRandomSuperoperator.choi_representation(map)
    )
end

@testset "WP2 bounded random superoperators" begin
    @testset "default trace-preserving construction" begin
        result = QETRandomSuperoperator.random_superoperator(Xoshiro(0x5153), 3)
        @test result isa QETRandomSuperoperator.RandomSuperoperatorResult
        @test result.status === :success
        @test result.last_attempt_status === :success
        @test result.succeeded
        @test result.attempts == 1
        @test result.input_dimension == 3
        @test result.output_dimension == 3
        @test result.requested_kraus_rank == 9
        @test result.numerical_kraus_rank == 9
        @test result.representation isa QETRandomSuperoperator.KrausRepresentation
        @test result.complete_positivity_guaranteed
        @test result.trace_preservation_guaranteed
        @test !result.unitality_guaranteed
        @test !result.proportional_unitality_guaranteed
        @test result.trace_preservation_residual <= result.trace_preservation_tolerance
        @test result.choi_trace ≈ 3 atol = 2e-13
        @test QETRandomSuperoperator.is_trace_preserving(result.representation)
        @test QETRandomSuperoperator.is_completely_positive(result.representation).status ===
            QETRandomSuperoperator.MatrixPredicateSatisfied
        @test occursin("status=success", sprint(show, result))

        input_marginal, _ = _random_superoperator_marginals(result)
        @test input_marginal ≈ Matrix{ComplexF64}(I, 3, 3) atol = 3e-13
        @test length(result.factor_singular_values) == 9
        @test all(>(result.rank_threshold), result.factor_singular_values)
        @test result.work_used > 0
    end

    @testset "all four pinned marginal branches" begin
        unconstrained = QETRandomSuperoperator.random_superoperator(
            Xoshiro(1), (2, 3); trace_preserving=false, unital=false, kraus_rank=3
        )
        @test unconstrained.status === :success
        @test unconstrained.choi_trace ≈ 1 atol = 2e-14
        @test !unconstrained.trace_preservation_guaranteed
        @test !unconstrained.unitality_guaranteed

        trace_preserving = QETRandomSuperoperator.random_superoperator(
            Xoshiro(2), (3, 2); kraus_rank=2
        )
        input_marginal, output_marginal = _random_superoperator_marginals(trace_preserving)
        @test input_marginal ≈ Matrix{ComplexF64}(I, 3, 3) atol = 3e-13
        @test norm(output_marginal - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test trace_preserving.choi_trace ≈ 3 atol = 3e-13

        unital = QETRandomSuperoperator.random_superoperator(
            Xoshiro(3), (2, 3); trace_preserving=false, unital=true, kraus_rank=2
        )
        input_marginal, output_marginal = _random_superoperator_marginals(unital)
        @test output_marginal ≈ Matrix{ComplexF64}(I, 3, 3) atol = 4e-13
        @test norm(input_marginal - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test unital.choi_trace ≈ 3 atol = 4e-13
        @test QETRandomSuperoperator.is_unital(unital.representation)

        bistochastic = QETRandomSuperoperator.random_superoperator(
            Xoshiro(4), 3; unital=true, kraus_rank=1
        )
        input_marginal, output_marginal = _random_superoperator_marginals(bistochastic)
        @test input_marginal ≈ Matrix{ComplexF64}(I, 3, 3) atol = 5e-13
        @test output_marginal ≈ Matrix{ComplexF64}(I, 3, 3) atol = 5e-13
        @test bistochastic.trace_preservation_guaranteed
        @test bistochastic.unitality_guaranteed
        @test bistochastic.balancing_iterations >= 1
        @test !isempty(bistochastic.balancing_residual_history)
        @test QETRandomSuperoperator.is_trace_preserving(bistochastic.representation)
        @test QETRandomSuperoperator.is_unital(bistochastic.representation)
    end

    @testset "unequal-dimensional upstream correction" begin
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(5), (2, 3); trace_preserving=true, unital=true
        )

        proportional = QETRandomSuperoperator.random_superoperator(
            Xoshiro(5),
            (2, 3);
            trace_preserving=true,
            unital=false,
            proportional_unital=true,
            kraus_rank=2,
        )
        input_marginal, output_marginal = _random_superoperator_marginals(proportional)
        @test proportional.status === :success
        @test proportional.trace_preservation_guaranteed
        @test !proportional.unitality_guaranteed
        @test proportional.proportional_unitality_guaranteed
        @test proportional.proportional_unital_factor == 2 / 3
        @test input_marginal ≈ Matrix{ComplexF64}(I, 2, 2) atol = 4e-8
        @test output_marginal ≈ (2 / 3) * Matrix{ComplexF64}(I, 3, 3) atol = 4e-8
        @test norm(output_marginal - Matrix{ComplexF64}(I, 3, 3)) > 0.1
        @test proportional.trace_preservation_residual <=
            proportional.trace_preservation_tolerance
        @test proportional.proportional_unitality_residual <=
            proportional.unitality_tolerance
    end

    @testset "representation selection and construction certificate" begin
        for representation in (:kraus, :choi, :superoperator)
            result = QETRandomSuperoperator.random_superoperator(
                Xoshiro(
                    0x600 +
                    Int(representation === :choi) +
                    2Int(representation === :superoperator),
                ),
                (2, 3);
                trace_preserving=false,
                kraus_rank=4,
                representation,
            )
            @test result.status === :success
            @test QETRandomSuperoperator.input_dimension(result.representation) == 2
            @test QETRandomSuperoperator.output_dimension(result.representation) == 3
            @test _random_superoperator_choi(result.representation) ≈
                _random_superoperator_choi(result.kraus_representation) atol = 3e-14
            @test result.representation isa (
                if representation === :kraus
                    QETRandomSuperoperator.KrausRepresentation
                elseif representation === :choi
                    QETRandomSuperoperator.ChoiRepresentation
                else
                    QETRandomSuperoperator.SuperoperatorRepresentation
                end
            )
        end
    end

    @testset "real and Float32 paths" begin
        result = QETRandomSuperoperator.random_superoperator(
            Xoshiro(7),
            2;
            unital=true,
            real=true,
            kraus_rank=2,
            representation=:choi,
            T=Float32,
        )
        @test result.status === :success
        @test result.real_output
        @test eltype(result.representation) == Float32
        @test eltype(result.kraus_representation) == Float32
        @test result.choi_trace isa Float32
        @test eltype(result.factor_singular_values) == Float32
        @test eltype(result.balancing_residual_history) == Float32
        @test isreal(_random_superoperator_choi(result.representation))
        @test result.trace_preservation_residual <= result.trace_preservation_tolerance
        @test result.unitality_residual <= result.unitality_tolerance
    end

    @testset "reproducibility and RNG isolation" begin
        first = QETRandomSuperoperator.random_superoperator(
            Xoshiro(0x715), (2, 3); kraus_rank=3, representation=:choi
        )
        second = QETRandomSuperoperator.random_superoperator(
            Xoshiro(0x715), (2, 3); kraus_rank=3, representation=:choi
        )
        @test first.status == second.status
        @test _random_superoperator_choi(first.representation) ==
            _random_superoperator_choi(second.representation)
        @test first.factor_singular_values == second.factor_singular_values

        explicit_rng = Xoshiro(0x716)
        control_rng = copy(explicit_rng)
        QETRandomSuperoperator.random_superoperator(explicit_rng, 2; trace_preserving=false)
        @test rand(explicit_rng, UInt64) != rand(control_rng, UInt64)

        Random.seed!(0x717)
        reference_first = rand()
        reference_second = rand()
        Random.seed!(0x717)
        observed_first = rand()
        QETRandomSuperoperator.random_superoperator(Xoshiro(0x718), 2; unital=true)
        observed_second = rand()
        @test observed_first == reference_first
        @test observed_second == reference_second
        @test !hasmethod(QETRandomSuperoperator.random_superoperator, Tuple{Int})
    end

    @testset "bounded failure is explicit data" begin
        nonconverged = QETRandomSuperoperator.random_superoperator(
            Xoshiro(8), 3; unital=true, kraus_rank=2, max_attempts=1, max_iterations=0
        )
        @test nonconverged.status === :max_attempts
        @test nonconverged.last_attempt_status === :max_iterations
        @test !nonconverged.succeeded
        @test nonconverged.attempts == 1
        @test nonconverged.complete_positivity_guaranteed
        @test !nonconverged.trace_preservation_guaranteed
        @test !nonconverged.unitality_guaranteed
        @test occursin("not established", nonconverged.message)
        @test QETRandomSuperoperator.is_completely_positive(
            nonconverged.kraus_representation
        ).status === QETRandomSuperoperator.MatrixPredicateSatisfied

        work_limited = QETRandomSuperoperator.random_superoperator(
            Xoshiro(9), 2; unital=true, max_attempts=1, max_iterations=100, max_work=400
        )
        @test work_limited.status === :work_limit
        @test work_limited.last_attempt_status === :work_limit
        @test !work_limited.succeeded
        @test work_limited.work_used <= 400
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(9), 2; unital=true, max_work=300
        )
    end

    @testset "strict validation and allocation preflight" begin
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), (3, 1); kraus_rank=2
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), (1, 3); trace_preserving=false, unital=true, kraus_rank=2
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; kraus_rank=5
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), (2, 3); trace_preserving=false, proportional_unital=true
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; unital=true, proportional_unital=true
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; proportional_unital=true
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; trace_preserving=1
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; unital=0
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; real=1
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; representation="choi"
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; representation=:matrix
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; T=BigFloat
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 0
        )
        @test_throws DimensionMismatch QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), (2, 3, 4)
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), _RandomSuperoperatorZeroBasedVector([2, 2])
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 3; max_dimension=2
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; max_entries=1
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; max_attempts=0
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; max_iterations=-1
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; max_condition_number=0
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; atol=-1
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; rtol=Inf
        )
        @test_throws ArgumentError QETRandomSuperoperator.random_superoperator(
            Xoshiro(10), 2; max_work=true
        )
    end

    @testset "MATLAB compatibility wrapper" begin
        raw = QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5253), (2, 3), 1, 0, 1, 2
        )
        diagnostic = QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5253), (2, 3), 1, 0, 1, 2; diagnostics=true
        )
        @test raw == QETRandomSuperoperator.choi_matrix(diagnostic.representation)
        @test size(raw) == (6, 6)
        @test eltype(raw) == Float64
        @test diagnostic.status === :success
        @test diagnostic.trace_preservation_guaranteed

        proportional = QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5254),
            (2, 3),
            1,
            1,
            0,
            2;
            diagnostics=true,
            allow_proportional_unital=true,
        )
        @test proportional.status === :success
        @test proportional.proportional_unitality_guaranteed
        @test !proportional.unitality_guaranteed
        @test proportional.proportional_unitality_residual <=
            proportional.unitality_tolerance

        @test_throws ArgumentError QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5255), (2, 3), 1, 1
        )
        @test_throws ArgumentError QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5255), 2, 2
        )
        @test_throws ArgumentError QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5255), 2, 1, 0, 0, 5
        )
        @test_throws ArgumentError QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5255), 2; diagnostics=1
        )
        @test_throws DomainError QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5256), 3, 1, 1, 0, 2; max_attempts=1, max_iterations=0
        )
        failed = QETRandomSuperoperator.MATLABCompat.RandomSuperoperator(
            Xoshiro(0x5256),
            3,
            1,
            1,
            0,
            2;
            diagnostics=true,
            max_attempts=1,
            max_iterations=0,
        )
        @test failed.status === :max_attempts
        @test !failed.succeeded
        @test !hasmethod(
            QETRandomSuperoperator.MATLABCompat.RandomSuperoperator, Tuple{Int}
        )
    end
end
