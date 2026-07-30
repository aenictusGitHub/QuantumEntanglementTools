using LinearAlgebra
using SparseArrays

const QETWP1 = QuantumEntanglementTools

struct _WP1ZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_WP1ZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_WP1ZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_WP1ZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_WP1ZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

@testset "WP1 general operator-space maps" begin
    @testset "validated operator spaces and compatibility accessors" begin
        rectangular = QETWP1.OperatorSpace(2, 3, 3, 4)
        @test isimmutable(rectangular)
        @test QETWP1.input_size(rectangular) == (2, 3)
        @test QETWP1.output_size(rectangular) == (3, 4)
        @test_throws ArgumentError QETWP1.input_dimension(rectangular)
        @test_throws ArgumentError QETWP1.output_dimension(rectangular)

        square = QETWP1.OperatorSpace((2, 2), (3, 3))
        @test square == QETWP1.OperatorSpace([2, 2], [3, 3])
        @test QETWP1.input_dimension(square) == 2
        @test QETWP1.output_dimension(square) == 3

        @test_throws ArgumentError QETWP1.OperatorSpace(0, 2, 3, 4)
        @test_throws ArgumentError QETWP1.OperatorSpace(true, 2, 3, 4)
        @test_throws ArgumentError QETWP1.OperatorSpace(2.0, 2, 3, 4)
        @test_throws ArgumentError QETWP1.OperatorSpace((2, 3, 4), (3, 4))
        @test_throws ArgumentError QETWP1.OperatorSpace(typemax(Int), 1, 2, 1)

        forged_token = QETWP1._ValidatedConstructorToken()
        @test_throws MethodError QETWP1.OperatorSpace(Val(:validated), 2, 3, 3, 4)
        @test_throws ArgumentError QETWP1.OperatorSpace(forged_token, 2, 3, 3, 4)

        legacy_choi = QETWP1.ChoiRepresentation(zeros(6, 6), 2, 3)
        legacy_super = QETWP1.SuperoperatorRepresentation(zeros(9, 4), 2, 3)
        for map in (legacy_choi, legacy_super)
            @test QETWP1.input_size(map) == (2, 2)
            @test QETWP1.output_size(map) == (3, 3)
            @test QETWP1.input_dimension(map) == 2
            @test QETWP1.output_dimension(map) == 3
            @test map.input_dim == 2
            @test map.output_dim == 3
        end
    end

    @testset "genuinely rectangular analytic operator sum" begin
        left1 = zeros(ComplexF64, 3, 2)
        left1[1, 1] = 1
        right1 = zeros(ComplexF64, 4, 3)
        right1[4, 3] = 1
        left2 = zeros(ComplexF64, 3, 2)
        left2[3, 2] = 1
        right2 = zeros(ComplexF64, 4, 3)
        right2[2, 1] = 2 - im

        map = QETWP1.OperatorSumRepresentation([left1, left2], [right1, right2])
        space = QETWP1.OperatorSpace((2, 3), (3, 4))
        @test QETWP1.operator_space(map) == space
        @test QETWP1.input_size(map) == (2, 3)
        @test QETWP1.output_size(map) == (3, 4)
        @test !hasproperty(map, :input_dim)
        @test !hasproperty(map, :output_dim)
        @test_throws ArgumentError QETWP1.input_dimension(map)
        @test_throws ArgumentError QETWP1.output_dimension(map)

        input = ComplexF64[1 2 3; 4 5 6]
        expected_output = zeros(ComplexF64, 3, 4)
        expected_output[1, 4] = 3
        expected_output[3, 2] = 8 + 4im

        expected_choi = zeros(ComplexF64, 6, 12)
        expected_choi[1, 12] = 1
        expected_choi[6, 2] = 2 + im
        expected_super = zeros(ComplexF64, 12, 6)
        expected_super[10, 5] = 1
        expected_super[6, 2] = 2 + im

        choi = QETWP1.choi_representation(map)
        superoperator = QETWP1.superoperator_representation(map)
        @test size(QETWP1.choi_matrix(map)) == (6, 12)
        @test size(QETWP1.superoperator_matrix(map)) == (12, 6)
        @test QETWP1.choi_matrix(map) == expected_choi
        @test QETWP1.superoperator_matrix(map) == expected_super
        @test QETWP1.operator_space(choi) == space
        @test QETWP1.operator_space(superoperator) == space
        @test_throws ArgumentError choi.input_dim
        @test_throws ArgumentError superoperator.output_dim

        @test QETWP1.apply_channel(input, map) == expected_output
        @test QETWP1.apply_channel(input, choi) == expected_output
        @test QETWP1.apply_channel(input, superoperator) == expected_output
        @test vec(expected_output) == expected_super * vec(input)
        @test_throws DimensionMismatch QETWP1.apply_channel(zeros(3, 2), map)

        @test QETWP1.superoperator_matrix(
            QETWP1.ChoiRepresentation(expected_choi, space)
        ) == expected_super
        @test QETWP1.choi_matrix(
            QETWP1.SuperoperatorRepresentation(expected_super, space)
        ) == expected_choi
        @test QETWP1.choi_matrix(
            QETWP1.choi_representation(QETWP1.superoperator_representation(choi))
        ) == expected_choi

        for input_column in 1:3, input_row in 1:2
            matrix_unit = zeros(ComplexF64, 2, 3)
            matrix_unit[input_row, input_column] = 1
            row_range = ((input_row - 1) * 3 + 1):(input_row * 3)
            column_range = ((input_column - 1) * 4 + 1):(input_column * 4)
            @test expected_choi[row_range, column_range] ==
                QETWP1.apply_channel(matrix_unit, map)
        end

        dual = QETWP1.dual_channel(map)
        output_test = reshape(ComplexF64.(1:12), 3, 4)
        @test QETWP1.input_size(dual) == (3, 4)
        @test QETWP1.output_size(dual) == (2, 3)
        @test dot(vec(output_test), vec(QETWP1.apply_channel(input, map))) ≈
            dot(vec(QETWP1.apply_channel(output_test, dual)), vec(input))
        for representation in (choi, superoperator)
            representation_dual = QETWP1.dual_channel(representation)
            @test QETWP1.input_size(representation_dual) == (3, 4)
            @test QETWP1.output_size(representation_dual) == (2, 3)
            forward_inner_product = dot(
                vec(output_test), vec(QETWP1.apply_channel(input, representation))
            )
            dual_inner_product = dot(
                vec(QETWP1.apply_channel(output_test, representation_dual)), vec(input)
            )
            @test forward_inner_product ≈ dual_inner_product
        end

        for diagnostic in (
            QETWP1.is_completely_positive,
            QETWP1.is_trace_preserving,
            QETWP1.is_unital,
            QETWP1.kraus_representation,
        )
            @test_throws ArgumentError diagnostic(map)
        end
        @test_throws DimensionMismatch QETWP1.partial_map(zeros(6, 6), map, 1, (2, 3))
    end

    @testset "ownership and validation" begin
        source_left = [1.0 0.0; 0.0 1.0]
        source_right = [1.0 2.0; 3.0 4.0]
        map = QETWP1.OperatorSumRepresentation([source_left], [source_right])
        expected_choi = QETWP1.choi_matrix(map)
        source_left[1, 1] = 99
        source_right[1, 1] = 88
        @test QETWP1.choi_matrix(map) == expected_choi
        @test_throws Base.CanonicalIndexError setindex!(map.left_operators, ones(2, 2), 1)
        @test_throws Base.CanonicalIndexError setindex!(map.right_operators, ones(2, 2), 1)

        factors = QETWP1.operator_sum_factors(map)
        @test propertynames(factors) == (:left, :right)
        factors.left[1][1, 1] = -7
        factors.right[1][1, 1] = -8
        @test QETWP1.choi_matrix(map) == expected_choi

        @test_throws ArgumentError QETWP1.OperatorSumRepresentation([], [])
        @test_throws DimensionMismatch QETWP1.OperatorSumRepresentation([ones(2, 2)], [])
        @test_throws DimensionMismatch QETWP1.OperatorSumRepresentation(
            [ones(2, 2), ones(3, 2)], [ones(2, 2), ones(2, 2)]
        )
        @test_throws DimensionMismatch QETWP1.OperatorSumRepresentation(
            [ones(2, 2), ones(2, 2)], [ones(2, 2), ones(3, 2)]
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation([1], [ones(1, 1)])
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation(
            [fill(NaN, 2, 2)], [ones(2, 2)]
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation(
            [ones(2, 2)], [fill(Inf, 2, 2)]
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation(
            [Matrix{Number}([1 0; 0 1])], [ones(2, 2)]
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation(
            [_WP1ZeroBasedMatrix(ones(2, 2))], [ones(2, 2)]
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation(
            [zeros(0, 2)], [zeros(2, 2)]
        )

        promoted = QETWP1.OperatorSumRepresentation(
            [ones(Float32, 2, 2)], [ones(Float64, 2, 2)]
        )
        @test eltype(promoted) == Float64
        mixed_left = AbstractMatrix{Float64}[ones(2, 2), sparse(ones(2, 2))]
        mixed_right = AbstractMatrix{Float64}[sparse(ones(2, 2)), ones(2, 2)]
        mixed = QETWP1.OperatorSumRepresentation(mixed_left, mixed_right)
        @test !issparse(mixed.left_operators[1])
        @test issparse(mixed.left_operators[2])
        @test issparse(mixed.right_operators[1])
        @test !issparse(mixed.right_operators[2])

        forged_token = QETWP1._ValidatedConstructorToken()
        left = [ones(2, 2)]
        right = [ones(2, 2)]
        @test_throws MethodError QETWP1.OperatorSumRepresentation{
            Float64,typeof(left),typeof(right)
        }(
            Val(:validated), left, right, QETWP1.OperatorSpace(2, 2, 2, 2)
        )
        @test_throws ArgumentError QETWP1.OperatorSumRepresentation{
            Float64,typeof(left),typeof(right)
        }(
            forged_token, left, right, QETWP1.OperatorSpace(2, 2, 2, 2)
        )

        corrupted = QETWP1.OperatorSumRepresentation([ones(2, 2)], [ones(2, 2)])
        corrupted.left_operators[1][1, 1] = NaN
        @test_throws ArgumentError QETWP1.apply_channel(ones(2, 2), corrupted)

        source_choi = reshape(ComplexF64.(1:72), 6, 12)
        stored_choi = QETWP1.ChoiRepresentation(
            source_choi, QETWP1.OperatorSpace((2, 3), (3, 4))
        )
        source_choi[1, 1] = -999
        @test QETWP1.choi_matrix(stored_choi)[1, 1] != -999
        source_super = reshape(ComplexF64.(1:72), 12, 6)
        stored_super = QETWP1.SuperoperatorRepresentation(
            source_super, QETWP1.OperatorSpace((2, 3), (3, 4))
        )
        source_super[1, 1] = -999
        @test QETWP1.superoperator_matrix(stored_super)[1, 1] != -999
        @test_throws DimensionMismatch QETWP1.ChoiRepresentation(
            zeros(6, 11), QETWP1.OperatorSpace((2, 3), (3, 4))
        )
        @test_throws DimensionMismatch QETWP1.SuperoperatorRepresentation(
            zeros(11, 6), QETWP1.OperatorSpace((2, 3), (3, 4))
        )
    end

    @testset "exact, sparse, and arbitrary-precision paths" begin
        left = Rational{Int}[1 0; 0 1; 1 -1]
        right = Rational{Int}[1 0 0; 0 1 1]
        input = Rational{Int}[1 2 3; 4 5 6]
        exact_map = QETWP1.OperatorSumRepresentation([left], [right])
        expected = left * input * right'
        exact_choi = QETWP1.choi_representation(exact_map)
        exact_super = QETWP1.superoperator_representation(exact_map)
        @test QETWP1.apply_channel(input, exact_map) == expected
        @test QETWP1.apply_channel(input, exact_choi) == expected
        @test QETWP1.apply_channel(input, exact_super) == expected
        @test QETWP1.choi_matrix(QETWP1.choi_representation(exact_super)) ==
            QETWP1.choi_matrix(exact_map)
        @test eltype(QETWP1.choi_matrix(exact_map)) == Rational{Int}
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(exact_choi)

        sparse_map = QETWP1.OperatorSumRepresentation([sparse(left)], [sparse(right)])
        sparse_input = sparse(input)
        sparse_choi = QETWP1.choi_representation(sparse_map)
        sparse_super = QETWP1.superoperator_representation(sparse_map)
        @test issparse(sparse_map.left_operators[1])
        @test issparse(sparse_map.right_operators[1])
        @test issparse(sparse_choi.matrix)
        @test issparse(sparse_super.matrix)
        @test issparse(QETWP1.apply_channel(sparse_input, sparse_map))
        @test QETWP1.apply_channel(sparse_input, sparse_map) == sparse(expected)
        @test QETWP1.choi_matrix(QETWP1.choi_representation(sparse_super)) ==
            QETWP1.choi_matrix(sparse_map)

        big_left = Complex{BigFloat}.(left)
        big_right = Complex{BigFloat}.(right)
        big_right[1, 2] = Complex{BigFloat}(1, 2)
        big_input = Complex{BigFloat}.(input)
        big_map = QETWP1.OperatorSumRepresentation([big_left], [big_right])
        big_expected = big_left * big_input * big_right'
        big_choi = QETWP1.choi_representation(big_map)
        big_super = QETWP1.superoperator_representation(big_map)
        @test eltype(QETWP1.choi_matrix(big_map)) == Complex{BigFloat}
        @test eltype(QETWP1.superoperator_matrix(big_map)) == Complex{BigFloat}
        @test QETWP1.apply_channel(big_input, big_map) == big_expected
        @test QETWP1.apply_channel(big_input, big_choi) == big_expected
        @test QETWP1.apply_channel(big_input, big_super) == big_expected
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(big_choi)
    end

    @testset "SVD decomposition policy and diagnostics" begin
        for scalar_type in (Float32, Float64, ComplexF32, ComplexF64)
            left1 = scalar_type[1 0; 0 1; 1 -1]
            right1 = scalar_type[1 0 0; 0 1 1]
            left2 = scalar_type[0 1; 1 0; 0 1]
            right2 = scalar_type[0 1 0; 1 0 -1]
            source = QETWP1.OperatorSumRepresentation([left1, left2], [right1, right2])
            choi = QETWP1.choi_representation(source)
            result = QETWP1.operator_sum_decomposition(choi)
            @test result.retained_rank == rank(QETWP1.choi_matrix(source))
            @test result.discarded_frobenius_norm <=
                result.threshold * sqrt(length(result.singular_values))
            @test QETWP1.choi_matrix(result.representation) ≈ QETWP1.choi_matrix(source)
            @test QETWP1.choi_matrix(
                QETWP1.operator_sum_representation(
                    QETWP1.superoperator_representation(source)
                ),
            ) ≈ QETWP1.choi_matrix(source)
        end

        transpose_super = zeros(Float64, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QETWP1.SuperoperatorRepresentation(transpose_super, 2, 2)
        @test QETWP1.is_completely_positive(transpose_map; atol=0, rtol=0).status ===
            QETWP1.MatrixPredicateViolated
        transpose_result = QETWP1.operator_sum_decomposition(transpose_map; atol=0, rtol=0)
        @test QETWP1.superoperator_matrix(transpose_result.representation) ≈ transpose_super

        nonhermitian_matrix = zeros(ComplexF64, 4, 4)
        nonhermitian_matrix[1, 2] = 1 + 2im
        nonhermitian_map = QETWP1.ChoiRepresentation(nonhermitian_matrix, 2, 2)
        @test QETWP1.is_completely_positive(nonhermitian_map).status ===
            QETWP1.MatrixPredicateViolated
        nonhermitian_result = QETWP1.operator_sum_decomposition(
            nonhermitian_map; atol=0, rtol=0
        )
        @test QETWP1.choi_matrix(nonhermitian_result.representation) ≈ nonhermitian_matrix

        rectangular_space = QETWP1.OperatorSpace((2, 3), (3, 4))
        zero_map = QETWP1.ChoiRepresentation(zeros(Float64, 6, 12), rectangular_space)
        zero_result = QETWP1.operator_sum_decomposition(zero_map)
        @test zero_result.retained_rank == 0
        @test zero_result.discarded_frobenius_norm == 0
        @test length(zero_result.representation.left_operators) == 1
        @test iszero(only(zero_result.representation.left_operators))
        @test iszero(only(zero_result.representation.right_operators))
        @test iszero(QETWP1.choi_matrix(zero_result.representation))

        small_singular_value = 1.0e-9
        truncated = QETWP1.ChoiRepresentation(
            Matrix(Diagonal([1.0, small_singular_value, 0.0, 0.0])), 2, 2
        )
        truncated_result = QETWP1.operator_sum_decomposition(truncated; atol=1.0e-8, rtol=0)
        @test truncated_result.threshold == 1.0e-8
        @test truncated_result.retained_rank == 1
        @test truncated_result.discarded_frobenius_norm ≈ small_singular_value
        @test norm(
            QETWP1.choi_matrix(truncated) -
            QETWP1.choi_matrix(truncated_result.representation),
        ) ≈ truncated_result.discarded_frobenius_norm

        sparse_source = QETWP1.ChoiRepresentation(sparse(Matrix{Float64}(I, 4, 4)), 2, 2)
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(sparse_source)
        sparse_result = QETWP1.operator_sum_decomposition(
            sparse_source; allow_densify=true, atol=0, rtol=0
        )
        @test QETWP1.choi_matrix(sparse_result.representation) ≈ Matrix{Float64}(I, 4, 4)
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(truncated; atol=-1)
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(truncated; rtol=Inf)
        @test_throws ArgumentError QETWP1.operator_sum_decomposition(
            truncated; atol=floatmax(Float64), rtol=floatmax(Float64)
        )
    end

    @testset "unequal square algebras retain legacy semantics" begin
        isometry = ComplexF64[1 0; 0 1; 0 0]
        general = QETWP1.OperatorSumRepresentation([isometry], [isometry])
        @test QETWP1.input_size(general) == (2, 2)
        @test QETWP1.output_size(general) == (3, 3)
        @test QETWP1.input_dimension(general) == 2
        @test QETWP1.output_dimension(general) == 3
        @test QETWP1.is_completely_positive(general).status ===
            QETWP1.MatrixPredicateUnknown
        @test QETWP1.is_trace_preserving(general)
        @test !QETWP1.is_unital(general)

        kraus = QETWP1.KrausRepresentation([isometry])
        paired = QETWP1.operator_sum_representation(kraus)
        @test paired isa QETWP1.OperatorSumRepresentation
        @test QETWP1.choi_matrix(paired) == QETWP1.choi_matrix(kraus)
        @test QETWP1.superoperator_matrix(paired) == QETWP1.superoperator_matrix(kraus)
    end
end
