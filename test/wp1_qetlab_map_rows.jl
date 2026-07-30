using LinearAlgebra
using Random
using SparseArrays

const QETWP1Rows = QuantumEntanglementTools
const QETWP1Compat = QuantumEntanglementTools.MATLABCompat

function _wp1_factor_cells(left, right)
    length(left) == length(right) || error("test fixture factor counts differ")
    cells = Matrix{AbstractMatrix}(undef, length(left), 2)
    for index in eachindex(left, right)
        cells[index, 1] = left[index]
        cells[index, 2] = right[index]
    end
    return cells
end

function _wp1_cp_column(operators)
    cells = Matrix{AbstractMatrix}(undef, length(operators), 1)
    for index in eachindex(operators)
        cells[index, 1] = operators[index]
    end
    return cells
end

function _wp1_cp_row(operators)
    cells = Matrix{AbstractMatrix}(undef, 1, length(operators))
    for index in eachindex(operators)
        cells[1, index] = operators[index]
    end
    return cells
end

@testset "WP1 pinned ApplyMap, ChoiMatrix, DualMap, PartialMap, IsHermPreserving" begin
    @testset "direct two-sided application and ApplyMap branches" begin
        rng = MersenneTwister(0x4150504c594d4150)
        left = [randn(rng, ComplexF64, 3, 2) for _ in 1:4]
        right = [randn(rng, ComplexF64, 4, 3) for _ in 1:4]
        input = randn(rng, ComplexF64, 2, 3)
        map = QETWP1Rows.OperatorSumRepresentation(left, right)
        expected = left[1] * input * adjoint(right[1])
        for index in 2:length(left)
            expected += left[index] * input * adjoint(right[index])
        end

        compatible, scalar_type = QETWP1Rows._can_use_strided_operator_sum_kernel(
            input, left, right
        )
        @test compatible
        @test scalar_type === ComplexF64
        @test QETWP1Rows.apply_channel(input, map) ≈ expected
        @test QETWP1Rows.apply_channel(input, QETWP1Rows.choi_representation(map)) ≈
            expected
        @test QETWP1Rows.apply_channel(
            input, QETWP1Rows.superoperator_representation(map)
        ) ≈ expected

        cells = _wp1_factor_cells(left, right)
        @test QETWP1Compat.ApplyMap(input, cells) ≈ expected
        @test QETWP1Compat.ApplyMap(input, QETWP1Rows.choi_matrix(map)) ≈ expected

        kraus = [randn(rng, ComplexF64, 3, 2) for _ in 1:3]
        square_input = randn(rng, ComplexF64, 2, 2)
        cp_expected = sum(operator * square_input * adjoint(operator) for operator in kraus)
        @test QETWP1Compat.ApplyMap(square_input, kraus) ≈ cp_expected
        @test QETWP1Compat.ApplyMap(square_input, _wp1_cp_column(kraus)) ≈ cp_expected
        @test QETWP1Compat.ApplyMap(square_input, _wp1_cp_row(kraus)) ≈ cp_expected

        rational_left = [Rational{Int}[1 0; 0 1; 1 -1], Rational{Int}[0 1; 1 0; 1 1]]
        rational_right = [Rational{Int}[1 0 0; 0 1 1], Rational{Int}[0 1 0; 1 0 -1]]
        rational_input = Rational{Int}[1 2 3; 4 5 6]
        rational_map = QETWP1Rows.OperatorSumRepresentation(rational_left, rational_right)
        rational_expected =
            rational_left[1] * rational_input * adjoint(rational_right[1]) +
            rational_left[2] * rational_input * adjoint(rational_right[2])
        @test QETWP1Rows.apply_channel(rational_input, rational_map) == rational_expected
        @test QETWP1Compat.ApplyMap(
            rational_input, _wp1_factor_cells(rational_left, rational_right)
        ) == rational_expected

        big_left = Base.map(matrix -> Complex{BigFloat}.(matrix), rational_left)
        big_right = Base.map(matrix -> Complex{BigFloat}.(matrix), rational_right)
        big_input = Complex{BigFloat}.(rational_input)
        big_map = QETWP1Rows.OperatorSumRepresentation(big_left, big_right)
        big_expected =
            big_left[1] * big_input * adjoint(big_right[1]) +
            big_left[2] * big_input * adjoint(big_right[2])
        @test QETWP1Rows.apply_channel(big_input, big_map) == big_expected

        sparse_left = sparse.(rational_left)
        sparse_right = sparse.(rational_right)
        sparse_input = sparse(rational_input)
        sparse_map = QETWP1Rows.OperatorSumRepresentation(sparse_left, sparse_right)
        sparse_output = QETWP1Rows.apply_channel(sparse_input, sparse_map)
        @test issparse(sparse_output)
        @test sparse_output == sparse(rational_expected)
        @test issparse(
            QETWP1Compat.ApplyMap(
                sparse_input, _wp1_factor_cells(sparse_left, sparse_right)
            ),
        )
        @test !issparse(
            QETWP1Compat.ApplyMap(
                Matrix(rational_input), _wp1_factor_cells(sparse_left, sparse_right)
            ),
        )

        mixed_map = QETWP1Rows.OperatorSumRepresentation(
            [Float32.(rational_left[1]), Float64.(rational_left[2])],
            [Float64.(rational_right[1]), Float32.(rational_right[2])],
        )
        mixed_output = QETWP1Rows.apply_channel(Float64.(rational_input), mixed_map)
        @test mixed_output ≈
            Float32.(rational_left[1]) *
              Float64.(rational_input) *
              adjoint(Float64.(rational_right[1])) +
              Float64.(rational_left[2]) *
              Float64.(rational_input) *
              adjoint(Float32.(rational_right[2]))

        malformed = Matrix{AbstractMatrix}(undef, 2, 3)
        fill!(malformed, ones(2, 2))
        @test_throws DimensionMismatch QETWP1Compat.ApplyMap(ones(2, 2), malformed)
        @test_throws DimensionMismatch QETWP1Compat.ApplyMap(ones(3, 2), cells)
        @test_throws ArgumentError QETWP1Compat.ApplyMap(
            _WP1ZeroBasedMatrix(ones(2, 3)), cells
        )
    end

    @testset "ChoiMatrix conventions and raw early return" begin
        left = [ComplexF64[1 0; 0 1; 1 -1], ComplexF64[0 1; 1 0; 1 1]]
        right = [ComplexF64[1 0 0; 0 1 1], ComplexF64[0 1 0; 1 0 -1]]
        cells = _wp1_factor_cells(left, right)
        map = QETWP1Rows.OperatorSumRepresentation(left, right)
        standard = QETWP1Rows.choi_matrix(map)
        @test QETWP1Compat.ChoiMatrix(cells) == standard
        @test QETWP1Compat.ChoiMatrix(cells, 2) == standard

        row_plan = QETWP1Rows.SubsystemPermutationPlan((2, 3), (2, 1))
        column_plan = QETWP1Rows.SubsystemPermutationPlan((3, 2), (2, 1))
        swapped = QETWP1Rows.permute_subsystems(standard, row_plan, column_plan)
        @test QETWP1Compat.ChoiMatrix(cells, 1) == swapped
        @test size(swapped) == (6, 6)

        raw = reshape(ComplexF64.(1:72), 6, 12)
        copied = QETWP1Compat.ChoiMatrix(raw, 0)
        @test copied == raw
        @test copied !== raw

        sparse_cells = _wp1_factor_cells(sparse.(left), sparse.(right))
        @test issparse(QETWP1Compat.ChoiMatrix(sparse_cells))
        @test !issparse(QETWP1Compat.ChoiMatrix(cells))
        @test_throws ArgumentError QETWP1Compat.ChoiMatrix(cells, 0)
        @test_throws DimensionMismatch QETWP1Compat.ChoiMatrix(
            _wp1_cp_row([ones(2, 2) for _ in 1:3])
        )
    end

    @testset "DualMap raw representations and Hilbert--Schmidt identity" begin
        rng = MersenneTwister(0x4455414c4d415000)
        left = [randn(rng, ComplexF64, 3, 2) for _ in 1:2]
        right = [randn(rng, ComplexF64, 4, 3) for _ in 1:2]
        cells = _wp1_factor_cells(left, right)
        raw_dual = QETWP1Compat.DualMap(cells)
        @test size(raw_dual) == size(cells)
        for index in eachindex(cells)
            @test raw_dual[index] == adjoint(cells[index])
        end

        arbitrary_cells = Matrix{AbstractMatrix}(undef, 2, 3)
        for index in eachindex(arbitrary_cells)
            arbitrary_cells[index] = fill(ComplexF64(index), 2, 2)
        end
        arbitrary_dual = QETWP1Compat.DualMap(arbitrary_cells, [99, 101])
        @test size(arbitrary_dual) == (2, 3)
        @test all(
            arbitrary_dual[index] == adjoint(arbitrary_cells[index]) for
            index in eachindex(arbitrary_cells)
        )

        map = QETWP1Rows.OperatorSumRepresentation(left, right)
        choi = QETWP1Rows.choi_matrix(map)
        dimensions = [2 3; 3 4]
        dual_choi = QETWP1Compat.DualMap(choi, dimensions)
        expected_dual = QETWP1Rows.choi_matrix(QETWP1Rows.dual_channel(map))
        @test dual_choi ≈ expected_dual

        typed_dual = QETWP1Compat.DualMap(map)
        @test typed_dual isa QETWP1Rows.OperatorSumRepresentation
        @test QETWP1Rows.choi_matrix(typed_dual) ≈ expected_dual

        input = randn(rng, ComplexF64, 2, 3)
        output = randn(rng, ComplexF64, 3, 4)
        @test dot(vec(output), vec(QETWP1Rows.apply_channel(input, map))) ≈
            dot(vec(QETWP1Rows.apply_channel(output, typed_dual)), vec(input))

        cp = [randn(rng, ComplexF64, 3, 2) for _ in 1:2]
        raw_cp_dual = QETWP1Compat.DualMap(cp)
        @test raw_cp_dual == [copy(adjoint(operator)) for operator in cp]
        @test_throws ArgumentError QETWP1Compat.DualMap(choi)
        @test_throws DimensionMismatch QETWP1Compat.DualMap(choi, [2 2; 3 3])
    end

    @testset "rectangular PartialMap layouts and subsystem ordering" begin
        row_dimensions = (2, 2, 3)
        column_dimensions = (3, 3, 2)
        factor1 = reshape(Rational{Int}.(1:6), 2, 3)
        factor2 = reshape(Rational{Int}.(2:7), 2, 3)
        factor3 = reshape(Rational{Int}.(1:6), 3, 2)
        input = kron(factor1, factor2, factor3)

        left2 = [Rational{Int}[1 0; 0 1; 1 -1], Rational{Int}[0 1; 1 0; 1 1]]
        right2 = [
            Rational{Int}[1 0 0; 0 1 1; 1 0 -1; 0 1 0],
            Rational{Int}[0 1 0; 1 0 -1; 0 0 1; 1 1 0],
        ]
        map2 = QETWP1Rows.OperatorSumRepresentation(left2, right2)
        mapped2 =
            left2[1] * factor2 * adjoint(right2[1]) +
            left2[2] * factor2 * adjoint(right2[2])
        expected2 = kron(factor1, mapped2, factor3)
        result2 = QETWP1Rows.partial_map(input, map2, 2, row_dimensions, column_dimensions)
        @test result2 == expected2
        @test size(result2) == (18, 24)
        @test QETWP1Rows.partial_map(
            input,
            QETWP1Rows.choi_representation(map2),
            2,
            row_dimensions,
            column_dimensions,
        ) == expected2
        @test QETWP1Rows.partial_map(
            input,
            QETWP1Rows.superoperator_representation(map2),
            2,
            row_dimensions,
            column_dimensions,
        ) == expected2

        left1 = Rational{Int}[1 0; 0 1; 1 1; 1 -1]
        right1 = Rational{Int}[1 0 0; 0 1 1]
        map1 = QETWP1Rows.OperatorSumRepresentation([left1], [right1])
        expected1 = kron(left1 * factor1 * adjoint(right1), factor2, factor3)
        @test QETWP1Rows.partial_map(input, map1, 1, row_dimensions, column_dimensions) ==
            expected1

        left3 = Rational{Int}[1 0 0; 0 1 1]
        right3 = Rational{Int}[1 0; 0 1; 1 -1; 0 1; 1 1]
        map3 = QETWP1Rows.OperatorSumRepresentation([left3], [right3])
        expected3 = kron(factor1, factor2, left3 * factor3 * adjoint(right3))
        @test QETWP1Rows.partial_map(input, map3, 3, row_dimensions, column_dimensions) ==
            expected3

        dimensions = [2 2 3; 3 3 2]
        cells2 = _wp1_factor_cells(left2, right2)
        @test QETWP1Compat.PartialMap(input, cells2, 2, dimensions) == expected2
        @test QETWP1Compat.PartialMap(input, QETWP1Rows.choi_matrix(map2), 2, dimensions) ==
            expected2

        two_factor_input = kron(factor1, factor2)
        two_factor_expected = kron(factor1, mapped2)
        @test QETWP1Compat.PartialMap(two_factor_input, cells2) == two_factor_expected

        square_factor1 = Rational{Int}[1 2; 3 4]
        square_factor2 = Rational{Int}[2 0; 1 3]
        square_input = kron(square_factor1, square_factor2)
        cp_operators = [
            Rational{Int}[1 0; 0 1], Rational{Int}[0 1; 1 0], Rational{Int}[1 0; 0 -1]
        ]
        mapped_square_factor2 = sum(
            operator * square_factor2 * adjoint(operator) for operator in cp_operators
        )
        square_expected = kron(square_factor1, mapped_square_factor2)
        @test QETWP1Compat.PartialMap(square_input, _wp1_cp_column(cp_operators), 2, 2) ==
            square_expected
        @test QETWP1Compat.PartialMap(square_input, _wp1_cp_row(cp_operators), 2, (2, 2)) ==
            square_expected

        sparse_input = sparse(input)
        sparse_map = QETWP1Rows.OperatorSumRepresentation(sparse.(left2), sparse.(right2))
        sparse_result = QETWP1Rows.partial_map(
            sparse_input, sparse_map, 2, row_dimensions, column_dimensions
        )
        @test issparse(sparse_result)
        @test sparse_result == sparse(expected2)
        @test issparse(
            QETWP1Compat.PartialMap(
                sparse_input,
                _wp1_factor_cells(sparse.(left2), sparse.(right2)),
                2,
                dimensions,
            ),
        )

        big_result = QETWP1Rows.partial_map(
            Complex{BigFloat}.(input),
            QETWP1Rows.OperatorSumRepresentation(
                Base.map(matrix -> Complex{BigFloat}.(matrix), left2),
                Base.map(matrix -> Complex{BigFloat}.(matrix), right2),
            ),
            2,
            row_dimensions,
            column_dimensions,
        )
        @test big_result == Complex{BigFloat}.(expected2)

        @test_throws DimensionMismatch QETWP1Rows.partial_map(
            input, map2, 2, (2, 2), column_dimensions
        )
        @test_throws DimensionMismatch QETWP1Rows.partial_map(
            input, map2, 2, (2, 3, 2), column_dimensions
        )
        @test_throws DimensionMismatch QETWP1Rows.partial_map(
            input, map2, 3, row_dimensions, column_dimensions
        )
        @test_throws ArgumentError QETWP1Rows.partial_map(
            input, map2, (1, 2), row_dimensions, column_dimensions
        )
        @test_throws ArgumentError QETWP1Compat.PartialMap(input, cells2, 0, dimensions)
        @test_throws DimensionMismatch QETWP1Compat.PartialMap(
            input, cells2, 2, [2 2 3; 3 3 3]
        )
        @test_throws DimensionMismatch QETWP1Compat.PartialMap(input, cells2, 2, (2, 2, 3))
        @test_throws DimensionMismatch QETWP1Compat.PartialMap(input, cells2, 2, 2)
        @test_throws ArgumentError QETWP1Compat.PartialMap(
            input, cells2, 2, _WP1ZeroBasedMatrix(dimensions)
        )
    end

    @testset "Hermiticity-preservation diagnostics" begin
        identity_map = QETWP1Rows.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
        identity_result = QETWP1Rows.is_hermiticity_preserving(identity_map)
        @test identity_result isa QETWP1Rows.MatrixPredicateResult
        @test identity_result.status === QETWP1Rows.MatrixPredicateSatisfied
        @test identity_result.value == 0
        @test identity_result.witness === nothing

        transpose_super = zeros(Rational{Int}, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QETWP1Rows.SuperoperatorRepresentation(transpose_super, 2, 2)
        transpose_result = QETWP1Rows.is_hermiticity_preserving(transpose_map)
        @test transpose_result.status === QETWP1Rows.MatrixPredicateSatisfied
        @test QETWP1Rows.is_completely_positive(transpose_map; atol=0, rtol=0).status ===
            QETWP1Rows.MatrixPredicateViolated

        nonpreserving = QETWP1Rows.OperatorSumRepresentation(
            [ComplexF64[1 0; 0 1]], [ComplexF64[1 0; 0 im]]
        )
        violated = QETWP1Rows.is_hermiticity_preserving(nonpreserving; atol=0, rtol=0)
        @test violated.status === QETWP1Rows.MatrixPredicateViolated
        @test violated.witness.kind === :nonhermitian
        @test violated.value > 0

        defect = 1.0e-10
        boundary_choi = ComplexF64[1 defect; 0 1]
        boundary_map = QETWP1Rows.ChoiRepresentation(
            boundary_choi, QETWP1Rows.OperatorSpace((1, 1), (2, 2))
        )
        boundary = QETWP1Rows.is_hermiticity_preserving(boundary_map; atol=defect, rtol=0)
        @test boundary.status === QETWP1Rows.MatrixPredicateUnknown
        @test boundary.witness.kind === :hermiticity_boundary
        @test QETWP1Rows.choi_matrix(boundary_map) == boundary_choi
        outside = QETWP1Rows.is_hermiticity_preserving(
            boundary_map; atol=defect / 2, rtol=0
        )
        @test outside.status === QETWP1Rows.MatrixPredicateViolated

        exact_hermitian = Rational{Int}[1 1//2; 1//2 2]
        exact_nonhermitian = Rational{Int}[1 1//2; 0 2]
        @test QETWP1Compat.IsHermPreserving(exact_hermitian).status ===
            QETWP1Rows.MatrixPredicateSatisfied
        @test QETWP1Compat.IsHermPreserving(exact_nonhermitian, 100).status ===
            QETWP1Rows.MatrixPredicateViolated

        sparse_result = QETWP1Compat.IsHermPreserving(sparse(exact_hermitian))
        @test sparse_result.status === QETWP1Rows.MatrixPredicateSatisfied
        big_boundary = QETWP1Compat.IsHermPreserving(
            Complex{BigFloat}.(boundary_choi), BigFloat(defect)
        )
        @test big_boundary.status === QETWP1Rows.MatrixPredicateUnknown

        cp_column = _wp1_cp_column([Rational{Int}[1 0; 0 1], Rational{Int}[0 1; 1 0]])
        @test QETWP1Compat.IsHermPreserving(cp_column).status ===
            QETWP1Rows.MatrixPredicateSatisfied
        hp_pair = _wp1_factor_cells([Rational{Int}[1 0; 0 -1]], [Rational{Int}[1 0; 0 -1]])
        @test QETWP1Compat.IsHermPreserving(hp_pair).status ===
            QETWP1Rows.MatrixPredicateSatisfied
        @test QETWP1Compat.IsHermPreserving(boundary_choi, defect).status ===
            QETWP1Rows.MatrixPredicateUnknown

        nonsquare = QETWP1Compat.IsHermPreserving(ones(2, 3))
        @test nonsquare.status === QETWP1Rows.MatrixPredicateViolated
        @test nonsquare.witness.kind === :nonsquare_choi

        rectangular = QETWP1Rows.OperatorSumRepresentation([ones(3, 2)], [ones(4, 3)])
        @test_throws ArgumentError QETWP1Rows.is_hermiticity_preserving(rectangular)
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(rectangular)
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(fill(NaN, 2, 2))
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(Matrix{Number}([1 0; 0 1]))
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(boundary_choi, -1)
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(boundary_choi, Inf)
        @test_throws ArgumentError QETWP1Compat.IsHermPreserving(boundary_choi, true)
        @test_throws DimensionMismatch QETWP1Compat.IsHermPreserving(
            _wp1_cp_row([ones(2, 2) for _ in 1:3])
        )
    end
end
