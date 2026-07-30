using LinearAlgebra
using SparseArrays

const QETWP1CP = QuantumEntanglementTools
const QETWP1CPCompat = QuantumEntanglementTools.MATLABCompat

function _wp1_operator_sum_from_cells(cells)
    size(cells, 2) == 2 || error("test fixture must have two factor columns")
    return QETWP1CP.OperatorSumRepresentation(
        [cells[index, 1] for index in axes(cells, 1)],
        [cells[index, 2] for index in axes(cells, 1)],
    )
end

function _wp1_factors_are_orthogonal(factors; atol=1.0e-10)
    for first_index in eachindex(factors), second_index in eachindex(factors)
        first_index == second_index && continue
        abs(dot(vec(factors[first_index]), vec(factors[second_index]))) <= atol ||
            return false
    end
    return true
end

@testset "WP1 pinned IsCP, KrausOperators, and ComplementaryMap" begin
    @testset "status-rich complete positivity" begin
        identity_kraus = QETWP1CP.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
        guaranteed = QETWP1CP.is_completely_positive(identity_kraus)
        @test guaranteed isa QETWP1CP.MatrixPredicateResult
        @test guaranteed.status === QETWP1CP.MatrixPredicateSatisfied
        @test guaranteed.witness === nothing

        positive_choi = QETWP1CP.ChoiRepresentation(
            Matrix(Diagonal([4.0, 3.0, 2.0, 1.0])), 2, 2
        )
        positive = QETWP1CP.is_completely_positive(positive_choi)
        @test positive.status === QETWP1CP.MatrixPredicateSatisfied
        @test positive.value > positive.tolerance

        transpose_super = zeros(Float64, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QETWP1CP.SuperoperatorRepresentation(transpose_super, 2, 2)
        transpose_result = QETWP1CP.is_completely_positive(transpose_map; atol=0, rtol=0)
        @test transpose_result.status === QETWP1CP.MatrixPredicateViolated
        @test transpose_result.witness.kind === :choi_positive_semidefiniteness

        nonhermitian_choi = ComplexF64[1 1.0e-4; 0 1]
        nonhermitian_map = QETWP1CP.ChoiRepresentation(
            nonhermitian_choi, QETWP1CP.OperatorSpace((1, 1), (2, 2))
        )
        nonhermitian = QETWP1CP.is_completely_positive(
            nonhermitian_map; atol=1.0e-8, rtol=0
        )
        @test nonhermitian.status === QETWP1CP.MatrixPredicateViolated
        @test nonhermitian.witness.kind === :choi_hermiticity

        hermiticity_boundary = QETWP1CP.is_completely_positive(
            nonhermitian_map; atol=1.0e-3, rtol=0
        )
        @test hermiticity_boundary.status === QETWP1CP.MatrixPredicateUnknown
        @test hermiticity_boundary.witness.result.status === QETWP1CP.MatrixPredicateUnknown

        negative_boundary_matrix = Matrix(Diagonal([1.0, 1.0, 1.0, -1.0e-10]))
        negative_boundary_map = QETWP1CP.ChoiRepresentation(negative_boundary_matrix, 2, 2)
        negative_boundary = QETWP1CP.is_completely_positive(
            negative_boundary_map; atol=1.0e-8, rtol=0
        )
        @test negative_boundary.status === QETWP1CP.MatrixPredicateUnknown
        @test QETWP1CP.choi_matrix(negative_boundary_map) == negative_boundary_matrix
        @test QETWP1CP.is_completely_positive(
            negative_boundary_map; atol=1.0e-12, rtol=0
        ).status === QETWP1CP.MatrixPredicateViolated

        exact_positive = QETWP1CP.ChoiRepresentation(
            Rational{Int}[1 0; 0 0], QETWP1CP.OperatorSpace((1, 1), (2, 2))
        )
        exact_negative = QETWP1CP.ChoiRepresentation(
            Rational{Int}[1 0; 0 -1], QETWP1CP.OperatorSpace((1, 1), (2, 2))
        )
        @test QETWP1CP.is_completely_positive(exact_positive).status ===
            QETWP1CP.MatrixPredicateSatisfied
        @test QETWP1CP.is_completely_positive(exact_negative).status ===
            QETWP1CP.MatrixPredicateViolated

        big_boundary_map = QETWP1CP.ChoiRepresentation(
            BigFloat[1 0; 0 -big"1e-50"], QETWP1CP.OperatorSpace((1, 1), (2, 2))
        )
        @test QETWP1CP.is_completely_positive(big_boundary_map).status ===
            QETWP1CP.MatrixPredicateUnknown

        sparse_positive = QETWP1CP.ChoiRepresentation(
            sparse(Matrix{Float64}(I, 4, 4)), 2, 2
        )
        @test_throws ArgumentError QETWP1CP.is_completely_positive(sparse_positive)
        @test QETWP1CP.is_completely_positive(sparse_positive; allow_densify=true).status ===
            QETWP1CP.MatrixPredicateSatisfied

        rectangular = QETWP1CP.OperatorSumRepresentation([ones(3, 2)], [ones(4, 3)])
        @test_throws ArgumentError QETWP1CP.is_completely_positive(rectangular)
        @test_throws ArgumentError QETWP1CP.is_completely_positive(
            identity_kraus; atol=true
        )
        @test_throws ArgumentError QETWP1CP.is_completely_positive(identity_kraus; rtol=Inf)
        @test_throws ArgumentError QETWP1CP.is_completely_positive(exact_positive; atol=1)
    end

    @testset "IsCP compatibility branches" begin
        identity = Matrix{Float64}(I, 2, 2)
        cp_column = _wp1_cp_column([identity])
        @test QETWP1CPCompat.IsCP(cp_column).status === QETWP1CP.MatrixPredicateSatisfied
        @test QETWP1CPCompat.IsCP([identity]).status === QETWP1CP.MatrixPredicateSatisfied

        full_rank = Matrix(Diagonal([4.0, 3.0, 2.0, 1.0]))
        @test QETWP1CPCompat.IsCP(full_rank).status === QETWP1CP.MatrixPredicateSatisfied

        transpose_map = QETWP1CP.SuperoperatorRepresentation(
            [
                1.0 0 0 0
                0 0 1 0
                0 1 0 0
                0 0 0 1
            ],
            2,
            2,
        )
        transpose_choi = QETWP1CP.choi_matrix(transpose_map)
        @test QETWP1CPCompat.IsCP(transpose_choi, 0).status ===
            QETWP1CP.MatrixPredicateViolated

        paired_cp = _wp1_factor_cells([identity], [identity])
        @test QETWP1CPCompat.IsCP(paired_cp).status === QETWP1CP.MatrixPredicateUnknown

        nonhermitian = ComplexF64[1 1.0e-10; 0 1]
        @test QETWP1CPCompat.IsCP(nonhermitian, 1.0e-8).status ===
            QETWP1CP.MatrixPredicateUnknown
        @test QETWP1CPCompat.IsCP(ones(2, 3)).status === QETWP1CP.MatrixPredicateViolated
        @test QETWP1CPCompat.IsCP(Rational{Int}[1 0; 0 0]).status ===
            QETWP1CP.MatrixPredicateSatisfied
        @test QETWP1CPCompat.IsCP(BigFloat[1 0; 0 -big"1e-50"]).status ===
            QETWP1CP.MatrixPredicateUnknown

        sparse_full_rank = sparse(full_rank)
        @test_throws ArgumentError QETWP1CPCompat.IsCP(sparse_full_rank)
        @test QETWP1CPCompat.IsCP(sparse_full_rank; allow_densify=true).status ===
            QETWP1CP.MatrixPredicateSatisfied
        @test_throws ArgumentError QETWP1CPCompat.IsCP(full_rank, -1)
        @test_throws ArgumentError QETWP1CPCompat.IsCP(full_rank, true)
        @test_throws DimensionMismatch QETWP1CPCompat.IsCP(
            _wp1_cp_row([identity, identity, identity])
        )
    end

    @testset "canonical CP, signed HP, and general SVD factors" begin
        full_rank = QETWP1CP.ChoiRepresentation(
            Matrix(Diagonal([4.0, 3.0, 2.0, 1.0])), 2, 2
        )
        cp_result = QETWP1CP.canonical_map_decomposition(full_rank)
        @test cp_result isa QETWP1CP.CanonicalMapDecompositionResult
        @test cp_result.classification === :completely_positive
        @test cp_result.complete_positivity.status === QETWP1CP.MatrixPredicateSatisfied
        @test cp_result.retained_rank == 4
        @test cp_result.discarded_frobenius_norm == 0
        @test cp_result.reconstruction_frobenius_norm <= 1.0e-12
        cp_factors = QETWP1CP.operator_sum_factors(cp_result.representation)
        @test cp_factors.left == cp_factors.right
        @test _wp1_factors_are_orthogonal(cp_factors.left)
        @test QETWP1CP.choi_matrix(cp_result.representation) ≈
            QETWP1CP.choi_matrix(full_rank)

        amplitude = 0.3
        supplied_kraus = QETWP1CP.KrausRepresentation([
            [1.0 0; 0 sqrt(1 - amplitude)], [0 sqrt(amplitude); 0 0]
        ])
        supplied_result = QETWP1CP.canonical_map_decomposition(supplied_kraus)
        @test supplied_result.classification === :completely_positive
        @test supplied_result.complete_positivity.status ===
            QETWP1CP.MatrixPredicateSatisfied
        @test QETWP1CP.choi_matrix(supplied_result.representation) ≈
            QETWP1CP.choi_matrix(supplied_kraus)

        raw_rank_deficient = QETWP1CP.choi_representation(supplied_kraus)
        boundary_result = QETWP1CP.canonical_map_decomposition(
            raw_rank_deficient; atol=0, rtol=0
        )
        @test boundary_result.classification === :hermiticity_preserving
        @test boundary_result.complete_positivity.status in
            (QETWP1CP.MatrixPredicateUnknown, QETWP1CP.MatrixPredicateViolated)
        @test QETWP1CP.choi_matrix(boundary_result.representation) ≈
            QETWP1CP.choi_matrix(raw_rank_deficient)

        transpose_super = zeros(Float64, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QETWP1CP.SuperoperatorRepresentation(transpose_super, 2, 2)
        hp_result = QETWP1CP.canonical_map_decomposition(transpose_map)
        @test hp_result.classification === :hermiticity_preserving
        @test hp_result.hermiticity_preservation.status ===
            QETWP1CP.MatrixPredicateSatisfied
        @test hp_result.complete_positivity.status === QETWP1CP.MatrixPredicateViolated
        hp_factors = QETWP1CP.operator_sum_factors(hp_result.representation)
        @test _wp1_factors_are_orthogonal(hp_factors.left)
        @test _wp1_factors_are_orthogonal(hp_factors.right)
        pair_signs = [
            if isapprox(hp_factors.left[index], hp_factors.right[index])
                1
            elseif isapprox(hp_factors.left[index], -hp_factors.right[index])
                -1
            else
                0
            end for index in eachindex(hp_factors.left)
        ]
        @test all(!iszero, pair_signs)
        first_negative = findfirst(==(-1), pair_signs)
        @test first_negative !== nothing
        @test all(==(1), pair_signs[1:(first_negative - 1)])
        @test all(==(-1), pair_signs[first_negative:end])
        @test QETWP1CP.superoperator_matrix(hp_result.representation) ≈ transpose_super

        rectangular_space = QETWP1CP.OperatorSpace((2, 3), (3, 4))
        general_choi = zeros(ComplexF64, 6, 12)
        general_choi[1, 2] = 2
        general_choi[4, 9] = 3im
        general_choi[6, 12] = -1
        general_map = QETWP1CP.ChoiRepresentation(general_choi, rectangular_space)
        general_result = QETWP1CP.canonical_map_decomposition(general_map)
        @test general_result.classification === :general
        @test general_result.complete_positivity === nothing
        @test general_result.hermiticity_preservation === nothing
        @test general_result.retained_rank == 3
        @test general_result.reconstruction_frobenius_norm <= 1.0e-12
        general_factors = QETWP1CP.operator_sum_factors(general_result.representation)
        @test all(==((3, 2)), size.(general_factors.left))
        @test all(==((4, 3)), size.(general_factors.right))
        @test _wp1_factors_are_orthogonal(general_factors.left)
        @test _wp1_factors_are_orthogonal(general_factors.right)
        @test QETWP1CP.choi_matrix(general_result.representation) ≈ general_choi

        zero_map = QETWP1CP.ChoiRepresentation(zeros(4, 4), 2, 2)
        zero_result = QETWP1CP.canonical_map_decomposition(zero_map)
        @test zero_result.classification === :hermiticity_preserving
        @test zero_result.retained_rank == 0
        zero_factors = QETWP1CP.operator_sum_factors(zero_result.representation)
        @test length(zero_factors.left) == 1
        @test iszero(only(zero_factors.left))
        @test iszero(only(zero_factors.right))

        sparse_map = QETWP1CP.ChoiRepresentation(
            sparse(Matrix(Diagonal([4.0, 3.0, 2.0, 1.0]))), 2, 2
        )
        @test_throws ArgumentError QETWP1CP.canonical_map_decomposition(sparse_map)
        @test QETWP1CP.canonical_map_decomposition(sparse_map; allow_densify=true).classification ===
            :completely_positive
        @test_throws ArgumentError QETWP1CP.canonical_map_decomposition(
            QETWP1CP.ChoiRepresentation(
                Rational{Int}[1 0; 0 1], QETWP1CP.OperatorSpace((1, 1), (2, 2))
            ),
        )
        @test_throws ArgumentError QETWP1CP.canonical_map_decomposition(
            QETWP1CP.ChoiRepresentation(
                BigFloat[1 0; 0 1], QETWP1CP.OperatorSpace((1, 1), (2, 2))
            ),
        )
        @test_throws ArgumentError QETWP1CP.canonical_map_decomposition(full_rank; atol=-1)
    end

    @testset "KrausOperators raw shapes and reconstruction" begin
        amplitude = 0.3
        supplied = [[1.0 0; 0 sqrt(1 - amplitude)], [0 sqrt(amplitude); 0 0]]
        canonical_cp = QETWP1CPCompat.KrausOperators(supplied)
        @test canonical_cp isa Vector
        @test QETWP1CP.choi_matrix(QETWP1CP.KrausRepresentation(canonical_cp)) ≈
            QETWP1CP.choi_matrix(QETWP1CP.KrausRepresentation(supplied))

        raw_rank_deficient = QETWP1CP.choi_matrix(QETWP1CP.KrausRepresentation(supplied))
        canonical_boundary = QETWP1CPCompat.KrausOperators(raw_rank_deficient)
        @test canonical_boundary isa Matrix
        @test size(canonical_boundary, 2) == 2
        @test QETWP1CP.choi_matrix(_wp1_operator_sum_from_cells(canonical_boundary)) ≈
            raw_rank_deficient
        diagnostics = QETWP1CPCompat.KrausOperators(raw_rank_deficient; diagnostics=true)
        @test diagnostics.classification === :hermiticity_preserving
        @test diagnostics.reconstruction_frobenius_norm <= 1.0e-12

        full_rank = Matrix(Diagonal([4.0, 3.0, 2.0, 1.0]))
        canonical_full_rank = QETWP1CPCompat.KrausOperators(full_rank)
        @test canonical_full_rank isa Vector
        @test QETWP1CP.choi_matrix(QETWP1CP.KrausRepresentation(canonical_full_rank)) ≈
            full_rank

        transpose_super = [
            1.0 0 0 0
            0 0 1 0
            0 1 0 0
            0 0 0 1
        ]
        transpose_choi = QETWP1CP.choi_matrix(
            QETWP1CP.SuperoperatorRepresentation(transpose_super, 2, 2)
        )
        signed = QETWP1CPCompat.KrausOperators(transpose_choi)
        @test signed isa Matrix
        @test size(signed) == (4, 2)
        @test QETWP1CP.choi_matrix(_wp1_operator_sum_from_cells(signed)) ≈ transpose_choi

        rectangular_choi = zeros(ComplexF64, 6, 12)
        rectangular_choi[1, 2] = 2
        rectangular_choi[4, 9] = 3im
        rectangular_choi[6, 12] = -1
        general = QETWP1CPCompat.KrausOperators(rectangular_choi, [2 3; 3 4])
        @test size(general) == (3, 2)
        general_map = _wp1_operator_sum_from_cells(general)
        @test QETWP1CP.input_size(general_map) == (2, 3)
        @test QETWP1CP.output_size(general_map) == (3, 4)
        @test QETWP1CP.choi_matrix(general_map) ≈ rectangular_choi

        zero_factors = QETWP1CPCompat.KrausOperators(zeros(4, 4))
        @test size(zero_factors) == (1, 2)
        @test iszero(zero_factors[1, 1])
        @test iszero(zero_factors[1, 2])

        sparse_full_rank = sparse(full_rank)
        @test_throws ArgumentError QETWP1CPCompat.KrausOperators(sparse_full_rank)
        @test QETWP1CPCompat.KrausOperators(sparse_full_rank; allow_densify=true) isa Vector
        @test_throws ArgumentError QETWP1CPCompat.KrausOperators(
            Rational{Int}[1 0; 0 1], (1, 2)
        )
        @test_throws DimensionMismatch QETWP1CPCompat.KrausOperators(
            rectangular_choi, (2, 3)
        )
        @test_throws DimensionMismatch QETWP1CPCompat.KrausOperators(
            _wp1_cp_row([ones(2, 2) for _ in 1:3])
        )
    end

    @testset "complementary dilation semantics" begin
        identity = Matrix{Rational{Int}}(I, 2, 2)
        minimal = QETWP1CP.KrausRepresentation([identity])
        redundant = QETWP1CP.KrausRepresentation([
            identity / 2, identity / 2, identity / 2, identity / 2
        ])
        minimal_complement = QETWP1CP.complementary_channel(minimal)
        redundant_complement = QETWP1CP.complementary_channel(redundant)
        @test QETWP1CP.output_dimension(minimal_complement) == 1
        @test QETWP1CP.output_dimension(redundant_complement) == 4
        input = Rational{Int}[1 2; 3 4]
        @test QETWP1CP.apply_channel(input, minimal_complement) ==
            reshape([tr(input)], 1, 1)
        @test QETWP1CP.is_completely_positive(minimal_complement).status ===
            QETWP1CP.MatrixPredicateSatisfied

        left = [Rational{Int}[1 0; 0 1], Rational{Int}[0 1; 1 0]]
        right = [Rational{Int}[1 0 0; 0 1 0], Rational{Int}[0 1 0; 0 0 1]]
        paired = QETWP1CP.OperatorSumRepresentation(left, right)
        paired_complement = QETWP1CP.complementary_channel(paired)
        @test paired_complement isa QETWP1CP.OperatorSumRepresentation
        @test QETWP1CP.input_size(paired_complement) == (2, 3)
        @test QETWP1CP.output_size(paired_complement) == (2, 2)
        rectangular_input = Rational{Int}[1 2 3; 4 5 6]
        expected = Matrix{Rational{Int}}(undef, 2, 2)
        for column in 1:2, row in 1:2
            expected[row, column] = tr(
                left[row] * rectangular_input * adjoint(right[column])
            )
        end
        @test QETWP1CP.apply_channel(rectangular_input, paired_complement) == expected

        paired_factors = QETWP1CP.operator_sum_factors(paired_complement)
        @test length(paired_factors.left) == 2
        @test all(==((2, 2)), size.(paired_factors.left))
        @test all(==((2, 3)), size.(paired_factors.right))

        sparse_paired = QETWP1CP.OperatorSumRepresentation(sparse.(left), sparse.(right))
        sparse_complement = QETWP1CP.complementary_channel(sparse_paired)
        sparse_factors = QETWP1CP.operator_sum_factors(sparse_complement)
        @test all(issparse, sparse_factors.left)
        @test all(issparse, sparse_factors.right)
        @test issparse(QETWP1CP.apply_channel(sparse(rectangular_input), sparse_complement))

        big_paired = QETWP1CP.OperatorSumRepresentation(
            Base.map(matrix -> Complex{BigFloat}.(matrix), left),
            Base.map(matrix -> Complex{BigFloat}.(matrix), right),
        )
        big_complement = QETWP1CP.complementary_channel(big_paired)
        @test eltype(big_complement) == Complex{BigFloat}
        @test QETWP1CP.apply_channel(
            Complex{BigFloat}.(rectangular_input), big_complement
        ) == Complex{BigFloat}.(expected)

        nonsquare_output = QETWP1CP.OperatorSumRepresentation([ones(2, 2)], [ones(3, 3)])
        @test_throws ArgumentError QETWP1CP.complementary_channel(nonsquare_output)

        paired_cells = _wp1_factor_cells(left, right)
        raw_paired_complement = QETWP1CPCompat.ComplementaryMap(paired_cells)
        @test size(raw_paired_complement) == (2, 2)
        @test QETWP1CP.choi_matrix(_wp1_operator_sum_from_cells(raw_paired_complement)) ==
            QETWP1CP.choi_matrix(paired_complement)

        cp_column = _wp1_cp_column([Float64.(identity)])
        raw_column_complement = QETWP1CPCompat.ComplementaryMap(cp_column)
        @test raw_column_complement isa Matrix
        @test size(raw_column_complement) == (2, 1)
        raw_vector_complement = QETWP1CPCompat.ComplementaryMap([Float64.(identity)])
        @test raw_vector_complement isa Vector
        @test length(raw_vector_complement) == 2

        paired_choi = QETWP1CP.choi_representation(
            QETWP1CP.OperatorSumRepresentation(
                Base.map(matrix -> Float64.(matrix), left),
                Base.map(matrix -> Float64.(matrix), right),
            ),
        )
        choi_complement = QETWP1CP.complementary_channel(paired_choi)
        @test choi_complement isa QETWP1CP.ChoiRepresentation
        @test QETWP1CP.input_size(choi_complement) == (2, 3)
        raw_choi_complement = QETWP1CPCompat.ComplementaryMap(
            QETWP1CP.choi_matrix(paired_choi), [2 2; 3 2]
        )
        @test raw_choi_complement isa AbstractMatrix

        sparse_choi = QETWP1CP.ChoiRepresentation(
            sparse(Matrix(Diagonal([4.0, 3.0, 2.0, 1.0]))), 2, 2
        )
        @test_throws ArgumentError QETWP1CP.complementary_channel(sparse_choi)
        @test QETWP1CP.complementary_channel(sparse_choi; allow_densify=true) isa
            QETWP1CP.ChoiRepresentation
        @test_throws DimensionMismatch QETWP1CPCompat.ComplementaryMap(
            _wp1_cp_row([ones(2, 2) for _ in 1:3])
        )
        @test_throws ArgumentError QETWP1CPCompat.ComplementaryMap(
            _wp1_factor_cells([ones(2, 2)], [ones(3, 3)])
        )
    end
end
