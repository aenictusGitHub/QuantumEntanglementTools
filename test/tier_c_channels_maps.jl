using LinearAlgebra
using Random
using SparseArrays

const QET = QuantumEntanglementTools
const CompatC = QuantumEntanglementTools.MATLABCompat

struct _TierCZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_TierCZeroBasedVector) = size(vector.storage)
Base.axes(vector::_TierCZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_TierCZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_TierCZeroBasedVector, index::Int) = vector.storage[index + 1]

struct _TierCZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_TierCZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_TierCZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_TierCZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierCZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

mutable struct _TierCMutableMatrix{T} <: AbstractMatrix{T}
    storage::Matrix{T}
    rows::Int
    columns::Int
    one_based::Bool
end

Base.size(matrix::_TierCMutableMatrix) = (matrix.rows, matrix.columns)
function Base.axes(matrix::_TierCMutableMatrix)
    if matrix.one_based
        return (Base.OneTo(matrix.rows), Base.OneTo(matrix.columns))
    end
    return (0:(matrix.rows - 1), 0:(matrix.columns - 1))
end
Base.IndexStyle(::Type{<:_TierCMutableMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierCMutableMatrix, row::Int, column::Int)
    shift = matrix.one_based ? 0 : 1
    return matrix.storage[row + shift, column + shift]
end
function Base.copy(matrix::_TierCMutableMatrix)
    return _TierCMutableMatrix(
        copy(matrix.storage), matrix.rows, matrix.columns, matrix.one_based
    )
end

@testset "Tier C channel/map representations" begin
    @testset "validated representation constructors" begin
        identity_kraus = QET.KrausRepresentation([[1.0 0.0; 0.0 1.0]])
        @test QET.input_dimension(identity_kraus) == 2
        @test QET.output_dimension(identity_kraus) == 2
        @test eltype(identity_kraus) == Float64

        source = [1.0 0.0; 0.0 1.0]
        copied = QET.KrausRepresentation([source])
        source[1, 1] = 7
        @test copied.operators[1][1, 1] == 1
        @test_throws Base.CanonicalIndexError setindex!(copied.operators, ones(1, 1), 1)
        @test size(copied.operators[1]) == (2, 2)
        @test QET.output_dimension(QET.complementary_channel(copied)) == 1

        sparse_kraus = QET.KrausRepresentation([sparse([1.0 0.0; 0.0 1.0])])
        @test issparse(sparse_kraus.operators[1])
        @test issparse(QET.choi_matrix(sparse_kraus))
        @test issparse(QET.superoperator_matrix(sparse_kraus))

        @test_throws ArgumentError QET.KrausRepresentation([])
        @test_throws DimensionMismatch QET.KrausRepresentation([ones(2, 2), ones(3, 2)])
        @test_throws ArgumentError QET.KrausRepresentation([fill(NaN, 2, 2)])
        @test_throws ArgumentError QET.KrausRepresentation([fill(Inf, 2, 2)])
        @test_throws ArgumentError QET.KrausRepresentation([fill("x", 2, 2)])
        @test_throws ArgumentError QET.KrausRepresentation([Matrix{Number}([1 0; 0 1])])

        @test_throws DimensionMismatch QET.ChoiRepresentation(ones(5, 5), 2, 3)
        @test_throws ArgumentError QET.ChoiRepresentation(ones(6, 6))
        @test_throws ArgumentError QET.ChoiRepresentation(fill(NaN, 4, 4))
        @test_throws ArgumentError QET.ChoiRepresentation(ones(4, 4); input_dim=2)
        @test_throws DimensionMismatch QET.SuperoperatorRepresentation(ones(4, 9), 2, 2)
        @test_throws ArgumentError QET.SuperoperatorRepresentation(ones(3, 4))

        forged_token = QET._ValidatedConstructorToken()
        invalid_kraus = [ones(1, 1)]
        @test_throws MethodError QET.KrausRepresentation{Float64,typeof(invalid_kraus)}(
            Val(:validated), invalid_kraus, 2, 2
        )
        @test_throws ArgumentError QET.KrausRepresentation{Float64,typeof(invalid_kraus)}(
            forged_token, invalid_kraus, 2, 2
        )
        invalid_matrix = ones(1, 1)
        @test_throws MethodError QET.ChoiRepresentation{Float64,typeof(invalid_matrix)}(
            Val(:validated), invalid_matrix, 2, 2
        )
        @test_throws ArgumentError QET.ChoiRepresentation{Float64,typeof(invalid_matrix)}(
            forged_token, invalid_matrix, 2, 2
        )
        @test_throws MethodError QET.SuperoperatorRepresentation{
            Float64,typeof(invalid_matrix)
        }(
            Val(:validated), invalid_matrix, 2, 2
        )
        @test_throws ArgumentError QET.SuperoperatorRepresentation{
            Float64,typeof(invalid_matrix)
        }(
            forged_token, invalid_matrix, 2, 2
        )

        zero_based_matrix = _TierCZeroBasedMatrix(Matrix{Float64}(I, 4, 4))
        @test_throws ArgumentError QET.ChoiRepresentation(zero_based_matrix)
        @test_throws ArgumentError QET.SuperoperatorRepresentation(zero_based_matrix)
        @test_throws ArgumentError QET.KrausRepresentation([
            _TierCZeroBasedMatrix(Matrix{Float64}(I, 2, 2))
        ])

        rectangular_super = QET.SuperoperatorRepresentation(zeros(9, 4))
        @test QET.input_dimension(rectangular_super) == 2
        @test QET.output_dimension(rectangular_super) == 3

        mutable_choi = QET.ChoiRepresentation(
            _TierCMutableMatrix(Matrix{Float64}(I, 4, 4), 4, 4, true)
        )
        mutable_choi.matrix.rows = 3
        @test_throws DimensionMismatch QET.superoperator_representation(mutable_choi)

        mutable_superoperator = QET.SuperoperatorRepresentation(
            _TierCMutableMatrix(Matrix{Float64}(I, 4, 4), 4, 4, true)
        )
        mutable_superoperator.matrix.one_based = false
        @test_throws ArgumentError QET.choi_representation(mutable_superoperator)

        mutable_kraus = QET.KrausRepresentation([
            _TierCMutableMatrix(Matrix{Float64}(I, 2, 2), 2, 2, true)
        ])
        mutable_kraus.operators[1].one_based = false
        @test_throws ArgumentError QET.apply_channel(
            Matrix{Float64}(I, 2, 2), mutable_kraus
        )
        @test_throws ArgumentError QET.is_completely_positive(mutable_kraus)

        nonfinite_choi = QET.ChoiRepresentation(Matrix{Float64}(I, 4, 4))
        nonfinite_choi.matrix[1, 1] = NaN
        @test_throws ArgumentError QET.superoperator_representation(nonfinite_choi)

        nonfinite_kraus = QET.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
        nonfinite_kraus.operators[1][1, 1] = NaN
        @test_throws ArgumentError QET.is_completely_positive(nonfinite_kraus)
    end

    @testset "Kraus, Choi, and superoperator round trips" begin
        amplitude = 0.37
        operators = [[1.0 0.0; 0.0 sqrt(1 - amplitude)], [0.0 sqrt(amplitude); 0.0 0.0]]
        kraus = QET.KrausRepresentation(operators)
        choi = QET.choi_representation(kraus)
        superoperator = QET.superoperator_representation(kraus)

        @test QET.choi_matrix(superoperator) ≈ QET.choi_matrix(choi)
        @test QET.superoperator_matrix(choi) ≈ QET.superoperator_matrix(superoperator)
        @test QET.choi_matrix(
            QET.choi_representation(QET.superoperator_representation(choi))
        ) ≈ QET.choi_matrix(choi)

        choi_data = QET.choi_matrix(choi)
        for input_column in 1:2, input_row in 1:2
            matrix_unit = zeros(2, 2)
            matrix_unit[input_row, input_column] = 1
            row_range = ((input_row - 1) * 2 + 1):(input_row * 2)
            column_range = ((input_column - 1) * 2 + 1):(input_column * 2)
            @test choi_data[row_range, column_range] ≈ QET.apply_channel(matrix_unit, kraus)
        end

        @test_throws DomainError QET.kraus_representation(choi)
        @test_throws DomainError QET.kraus_representation(superoperator)
        recovered = QET.canonical_map_decomposition(choi)
        @test QET.choi_matrix(recovered.representation) ≈ QET.choi_matrix(choi)
        @test recovered.complete_positivity.status === MatrixPredicateUnknown

        full_rank_choi = QET.ChoiRepresentation(
            Matrix(Diagonal([4.0, 3.0, 2.0, 1.0])), 2, 2
        )
        full_rank_decomposition = QET.canonical_map_decomposition(full_rank_choi)
        @test full_rank_decomposition.complete_positivity.status ===
            MatrixPredicateSatisfied
        @test all(
            >(full_rank_decomposition.threshold), full_rank_decomposition.spectral_values
        )
        @test QET.choi_matrix(QET.kraus_representation(full_rank_choi)) ≈
            QET.choi_matrix(full_rank_choi)

        tiny_positive_choi = QET.ChoiRepresentation(
            Matrix(Diagonal([1.0, 1.0, 1.0, 1.0e-10])), 2, 2
        )
        tiny_positive_decomposition = QET.canonical_map_decomposition(
            tiny_positive_choi; atol=0, rtol=0
        )
        @test tiny_positive_decomposition.complete_positivity.status ===
            MatrixPredicateSatisfied
        @test all(
            >(tiny_positive_decomposition.threshold),
            tiny_positive_decomposition.spectral_values,
        )
        tiny_positive_kraus = QET.kraus_representation(tiny_positive_choi; atol=0, rtol=0)
        @test QET.choi_matrix(tiny_positive_kraus) ≈ QET.choi_matrix(tiny_positive_choi) atol =
            1.0e-14 rtol = 1.0e-14

        # The PSD diagnostic scales from the largest matrix entry, while the
        # factorization cutoff scales from the largest eigenvalue. This map is
        # robustly CP for the former but would lose three positive modes under
        # the latter. The plain conversion must not silently change the map.
        roundoff = sqrt(eps(Float64))
        scale_mismatch_choi = QET.ChoiRepresentation(
            ones(4, 4) + 3 * roundoff * Matrix{Float64}(I, 4, 4), 2, 2
        )
        scale_mismatch = QET.canonical_map_decomposition(scale_mismatch_choi)
        @test scale_mismatch.complete_positivity.status === MatrixPredicateSatisfied
        @test scale_mismatch.retained_rank < length(scale_mismatch.spectral_values)
        @test scale_mismatch.discarded_frobenius_norm > 0
        @test scale_mismatch.reconstruction_frobenius_norm > 0
        @test_throws DomainError QET.kraus_representation(scale_mismatch_choi)

        tiny_negative_choi = QET.ChoiRepresentation(
            Matrix(Diagonal([1.0, 1.0, 1.0, -1.0e-10])), 2, 2
        )
        approximate = QET.canonical_map_decomposition(
            tiny_negative_choi; atol=1.0e-8, rtol=0
        )
        @test approximate isa QET.CanonicalMapDecompositionResult
        @test approximate.discarded_frobenius_norm == 1.0e-10
        @test approximate.reconstruction_frobenius_norm ≈ 1.0e-10
        @test_throws DomainError QET.kraus_representation(
            tiny_negative_choi; atol=1.0e-8, rtol=0
        )

        nearly_hermitian_data = Matrix{Float64}(I, 4, 4)
        nearly_hermitian_data[1, 2] = 1.0e-10
        nearly_hermitian_choi = QET.ChoiRepresentation(nearly_hermitian_data, 2, 2)
        @test_throws DomainError QET.kraus_representation(
            nearly_hermitian_choi; atol=1.0e-8, rtol=0
        )

        input = ComplexF64[0.6 0.2im; -0.2im 0.4]
        expected = sum(operator * input * operator' for operator in operators)
        @test QET.apply_channel(input, kraus) ≈ expected
        @test QET.apply_channel(input, choi) ≈ expected
        @test QET.apply_channel(input, superoperator) ≈ expected
        @test QET.apply_channel(input, operators) ≈ expected

        rng = MersenneTwister(0x43484f49)
        dense_operators = [randn(rng, ComplexF64, 4, 3) for _ in 1:5]
        dense_kraus = QET.KrausRepresentation(dense_operators)
        dense_choi_reference = sum(
            vec(operator) * adjoint(vec(operator)) for operator in dense_operators
        )
        @test QET.choi_matrix(dense_kraus) ≈ dense_choi_reference
        @test QET.choi_matrix(dense_kraus) isa Matrix{ComplexF64}

        raw_choi = QET.choi_matrix(choi)
        raw_choi[1, 1] = 99
        @test QET.choi_matrix(choi)[1, 1] != 99

        @test_throws DimensionMismatch QET.apply_channel(ones(3, 3), kraus)
        @test_throws ArgumentError QET.apply_channel(fill(NaN, 2, 2), kraus)
        @test_throws ArgumentError QET.apply_channel(
            _TierCZeroBasedMatrix(Matrix{Float64}(I, 2, 2)), kraus
        )

        zero_map = QET.ChoiRepresentation(zeros(4, 4))
        @test_throws DomainError QET.kraus_representation(zero_map)
        zero_decomposition = QET.canonical_map_decomposition(zero_map)
        @test zero_decomposition.retained_rank == 0
        @test iszero(QET.choi_matrix(zero_decomposition.representation))
    end

    @testset "rectangular maps and complex precision" begin
        isometry = ComplexF64[
            1 0
            0 1
            0 0
        ]
        map = QET.KrausRepresentation([isometry])
        @test QET.input_dimension(map) == 2
        @test QET.output_dimension(map) == 3
        @test size(QET.choi_matrix(map)) == (6, 6)
        @test size(QET.superoperator_matrix(map)) == (9, 4)
        @test QET.is_completely_positive(map).status === MatrixPredicateSatisfied
        @test QET.is_trace_preserving(map)
        @test !QET.is_unital(map)

        input = ComplexF64[0.7 0.1-0.2im; 0.1+0.2im 0.3]
        expected = isometry * input * isometry'
        for representation in
            (map, QET.choi_representation(map), QET.superoperator_representation(map))
            @test QET.apply_channel(input, representation) ≈ expected
        end

        big_operator = Complex{BigFloat}[
            big"1.0" 0
            0 cis(big"0.125")
        ]
        big_map = QET.KrausRepresentation([big_operator])
        big_input = Complex{BigFloat}[
            big"0.4" big"0.1"+big"0.2"*im
            big"0.1"-big"0.2"*im big"0.6"
        ]
        big_expected = big_operator * big_input * big_operator'
        @test eltype(QET.choi_matrix(big_map)) == Complex{BigFloat}
        @test QET.apply_channel(big_input, QET.superoperator_representation(big_map)) ≈
            big_expected
        @test QET.apply_channel(big_input, QET.choi_representation(big_map)) ≈ big_expected
    end

    @testset "physicality diagnostics" begin
        identity_map = QET.KrausRepresentation([Matrix{Float64}(I, 2, 2)])
        @test QET.is_completely_positive(identity_map).status === MatrixPredicateSatisfied
        @test QET.is_trace_preserving(identity_map)
        @test QET.is_unital(identity_map)

        transpose_super = zeros(Float64, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QET.SuperoperatorRepresentation(transpose_super, 2, 2)
        @test QET.is_completely_positive(transpose_map; atol=1e-14, rtol=1e-14).status ===
            MatrixPredicateViolated
        @test QET.is_trace_preserving(transpose_map)
        @test QET.is_unital(transpose_map)
        @test_throws DomainError QET.kraus_representation(transpose_map)

        nonhermitian = QET.ChoiRepresentation(
            ComplexF64[
                1 1 0 0
                0 0 0 0
                0 0 0 0
                0 0 0 1
            ],
            2,
            2,
        )
        @test QET.is_completely_positive(nonhermitian).status === MatrixPredicateViolated
        @test_throws DomainError QET.kraus_representation(nonhermitian)

        scaled_identity = QET.KrausRepresentation([2.0 * Matrix{Float64}(I, 2, 2)])
        @test QET.is_completely_positive(scaled_identity).status ===
            MatrixPredicateSatisfied
        @test !QET.is_trace_preserving(scaled_identity)
        @test !QET.is_unital(scaled_identity)

        @test_throws ArgumentError QET.is_completely_positive(identity_map; atol=-1)
        @test_throws ArgumentError QET.is_trace_preserving(identity_map; rtol=Inf)
    end

    @testset "analytic channel constructors" begin
        density = ComplexF64[0.7 0.2+0.1im; 0.2-0.1im 0.3]
        identity2 = Matrix{ComplexF64}(I, 2, 2)

        completely_depolarizing = QET.depolarizing_channel(2)
        @test QET.apply_channel(density, completely_depolarizing) ≈
            tr(density) * identity2 / 2
        @test QET.is_completely_positive(completely_depolarizing; allow_densify=true).status ===
            MatrixPredicateSatisfied
        @test QET.is_trace_preserving(completely_depolarizing)
        @test QET.is_unital(completely_depolarizing)
        @test issparse(QET.choi_matrix(completely_depolarizing))

        identity_depolarizing = QET.depolarizing_channel(2, 1)
        @test QET.apply_channel(density, identity_depolarizing) ≈ density
        @test QET.is_completely_positive(
            QET.depolarizing_channel(2, -1 // 3); allow_densify=true
        ).status === MatrixPredicateSatisfied
        @test_throws DomainError QET.depolarizing_channel(2, -0.34)
        @test_throws DomainError QET.depolarizing_channel(2, 1.01)

        complete_dephasing = QET.dephasing_channel(2)
        @test QET.apply_channel(density, complete_dephasing) ≈ Diagonal(diag(density))
        @test QET.apply_channel(density, QET.dephasing_channel(2, 1)) ≈ density
        @test QET.is_completely_positive(
            QET.dephasing_channel(3, -1 // 2); allow_densify=true
        ).status === MatrixPredicateSatisfied
        @test_throws DomainError QET.dephasing_channel(3, -0.51)

        bit_flip = QET.pauli_channel([0.0, 1.0, 0.0, 0.0])
        x = ComplexF64[0 1; 1 0]
        @test QET.apply_channel(density, bit_flip) ≈ x * density * x
        @test QET.is_completely_positive(bit_flip; allow_densify=true).status ===
            MatrixPredicateUnknown
        @test QET.is_trace_preserving(bit_flip)
        @test QET.is_unital(bit_flip)
        @test issparse(QET.choi_matrix(bit_flip))

        two_qubit_identity = QET.pauli_channel(reshape([1.0; zeros(15)], 4, 4))
        @test QET.apply_channel(Matrix{Float64}(I, 4, 4), two_qubit_identity) ≈
            Matrix{Float64}(I, 4, 4)
        @test_throws ArgumentError QET.pauli_channel([0.5, 0.5])
        @test_throws ArgumentError QET.pauli_channel(fill(1 / 8, 8))
        @test_throws DomainError QET.pauli_channel([0.5, 0.5, 0.1, -0.1])
        @test_throws DomainError QET.pauli_channel(fill(0.2, 4))
        @test_throws ArgumentError QET.pauli_channel([0.5, 0.5, 0.0, NaN])
        @test QET.is_trace_preserving(QET.pauli_channel(Real[1.0, 0.0, 0.0, 0.0]))
        @test_throws ArgumentError QET.pauli_channel(ComplexF64[1, 0, 0, 0])
        @test_throws ArgumentError QET.pauli_channel(
            _TierCZeroBasedVector([1.0, 0.0, 0.0, 0.0])
        )
    end

    @testset "dual, complementary, and partial maps" begin
        amplitude = 0.2
        map = QET.KrausRepresentation([
            [1.0 0.0; 0.0 sqrt(1 - amplitude)], [0.0 sqrt(amplitude); 0.0 0.0]
        ])
        dual = QET.dual_channel(map)
        x = ComplexF64[0.7 0.1im; -0.1im 0.3]
        y = ComplexF64[0.2 0.05; 0.05 0.8]
        @test tr(y' * QET.apply_channel(x, map)) ≈ tr(QET.apply_channel(y, dual)' * x)
        @test QET.choi_matrix(QET.dual_channel(QET.choi_representation(map))) ≈
            QET.choi_matrix(QET.choi_representation(dual))
        @test QET.superoperator_matrix(
            QET.dual_channel(QET.superoperator_representation(map))
        ) ≈ QET.superoperator_matrix(QET.superoperator_representation(dual))

        complement = QET.complementary_channel(map)
        @test QET.input_dimension(complement) == 2
        @test QET.output_dimension(complement) == 2
        @test QET.is_completely_positive(complement).status === MatrixPredicateSatisfied
        @test QET.is_trace_preserving(complement)
        choi_complement = QET.complementary_channel(QET.choi_representation(map))
        @test QET.input_dimension(choi_complement) == 2
        @test QET.output_dimension(choi_complement) == 2
        @test QET.is_completely_positive(choi_complement; allow_densify=true).status ===
            MatrixPredicateUnknown
        @test QET.is_trace_preserving(choi_complement)

        rho_a = ComplexF64[0.6 0.1; 0.1 0.4]
        rho_b = ComplexF64[0.3 0.05im; -0.05im 0.7]
        product = kron(rho_a, rho_b)
        dephase = QET.dephasing_channel(2)
        @test QET.partial_map(product, dephase, 1, (2, 2)) ≈
            kron(Diagonal(diag(rho_a)), rho_b)
        @test QET.partial_map(product, dephase, 2, (2, 2)) ≈
            kron(rho_a, Diagonal(diag(rho_b)))

        rectangular = QET.KrausRepresentation([ComplexF64[
            1 0
            0 1
            0 0
        ]])
        mapped = QET.partial_map(product, rectangular, 2, (2, 2))
        @test size(mapped) == (6, 6)
        @test mapped ≈
            kron(rho_a, rectangular.operators[1] * rho_b * rectangular.operators[1]')

        sparse_product = sparse(product)
        sparse_dephase = QET.dephasing_channel(2)
        @test issparse(QET.partial_map(sparse_product, sparse_dephase, 2, (2, 2)))
        @test_throws DimensionMismatch QET.partial_map(product, dephase, 1, (4, 1))
        @test_throws ArgumentError QET.partial_map(product, dephase, 0, (2, 2))
        @test_throws ArgumentError QET.partial_map(product, dephase, (1, 2), (2, 2))
    end

    @testset "positive-map fixtures" begin
        reduction = QET.reduction_map(2)
        matrix = ComplexF64[1 2im; -2im 3]
        @test QET.apply_channel(matrix, reduction) ≈
            tr(matrix) * Matrix{ComplexF64}(I, 2, 2) - matrix
        @test QET.is_completely_positive(reduction; atol=0, rtol=0, allow_densify=true).status ===
            MatrixPredicateViolated

        generalized_reduction = QET.reduction_map(3, 2)
        matrix3 = reshape(ComplexF64.(1:9), 3, 3)
        @test QET.apply_channel(matrix3, generalized_reduction) ≈
            2 * tr(matrix3) * Matrix{ComplexF64}(I, 3, 3) - matrix3

        choi = QET.choi_map()
        @test size(QET.choi_matrix(choi)) == (9, 9)
        @test QET.is_completely_positive(choi; atol=0, rtol=0, allow_densify=true).status ===
            MatrixPredicateViolated
        @test QET.choi_matrix(choi) == QET.choi_matrix(QET.choi_map(1, 1, 0))
        choi_input = reshape(ComplexF64.(1:9), 3, 3)
        @test QET.apply_channel(choi_input, choi) ≈ ComplexF64[
            choi_input[1, 1]+choi_input[2, 2] -choi_input[1, 2] -choi_input[1, 3]
            -choi_input[2, 1] choi_input[2, 2]+choi_input[3, 3] -choi_input[2, 3]
            -choi_input[3, 1] -choi_input[3, 2] choi_input[1, 1]+choi_input[3, 3]
        ]
        @test_throws ArgumentError QET.choi_map(Inf, 1, 0)
        @test_throws ArgumentError QET.choi_map(1 + im, 1, 0)
        @test_throws ArgumentError QET.reduction_map(0)
    end

    @testset "MATLAB/QETLAB compatibility wrappers" begin
        identity2 = Matrix{Float64}(I, 2, 2)
        bit_flip = Float64[0 1; 1 0]
        density = ComplexF64[0.7 0.1+0.2im; 0.1-0.2im 0.3]
        kraus = [sqrt(0.75) * identity2, sqrt(0.25) * bit_flip]

        raw_choi = CompatC.ChoiMatrix(kraus)
        @test raw_choi ≈ QET.choi_matrix(QET.KrausRepresentation(kraus))
        @test CompatC.ApplyMap(density, kraus) ≈
            sum(operator * density * operator' for operator in kraus)
        @test CompatC.ApplyMap(density, raw_choi) ≈ CompatC.ApplyMap(density, kraus)

        copied_choi = CompatC.ChoiMatrix(raw_choi, 99)
        @test copied_choi == raw_choi
        @test copied_choi !== raw_choi

        rectangular = [ComplexF64[
            1 0
            0 1
            0 0
        ]]
        rectangular_map = QET.KrausRepresentation(rectangular)
        swapped_choi = CompatC.ChoiMatrix(rectangular, 1)
        @test swapped_choi ≈ QET.permute_subsystems(
            QET.choi_matrix(rectangular_map), (2, 3); permutation=(2, 1)
        )
        @test size(CompatC.ApplyMap(density, CompatC.ChoiMatrix(rectangular))) == (3, 3)

        canonical = CompatC.KrausOperators(raw_choi)
        @test canonical isa Matrix
        canonical_map = QET.OperatorSumRepresentation(
            [canonical[index, 1] for index in axes(canonical, 1)],
            [canonical[index, 2] for index in axes(canonical, 1)],
        )
        @test QET.choi_matrix(canonical_map) ≈ raw_choi
        rectangular_factors = CompatC.KrausOperators(
            CompatC.ChoiMatrix(rectangular), (2, 3)
        )
        @test size(rectangular_factors) == (1, 2)
        @test QET.choi_matrix(
            QET.OperatorSumRepresentation(
                [rectangular_factors[1, 1]], [rectangular_factors[1, 2]]
            ),
        ) ≈ CompatC.ChoiMatrix(rectangular)

        dual_kraus = CompatC.DualMap(kraus)
        @test dual_kraus == [operator' for operator in kraus]
        @test CompatC.DualMap(raw_choi) ≈
            QET.choi_matrix(QET.dual_channel(QET.ChoiRepresentation(raw_choi)))

        complementary_kraus = CompatC.ComplementaryMap(kraus)
        @test complementary_kraus isa Vector
        @test QET.is_trace_preserving(QET.KrausRepresentation(complementary_kraus))
        complementary_choi = CompatC.ComplementaryMap(raw_choi)
        @test complementary_choi isa AbstractMatrix
        @test QET.is_trace_preserving(QET.ChoiRepresentation(complementary_choi))

        rho_a = ComplexF64[0.6 0.1im; -0.1im 0.4]
        rho_b = ComplexF64[0.3 0.05; 0.05 0.7]
        product = kron(rho_a, rho_b)
        raw_dephasing = CompatC.DephasingChannel(2)
        @test CompatC.PartialMap(product, raw_dephasing) ≈
            kron(rho_a, Diagonal(diag(rho_b)))
        @test CompatC.PartialMap(product, raw_dephasing, 1, (2, 2)) ≈
            kron(Diagonal(diag(rho_a)), rho_b)

        @test CompatC.DepolarizingChannel(3, 0.2) ≈
            QET.choi_matrix(QET.depolarizing_channel(3, 0.2))
        @test CompatC.DephasingChannel(3, 0.25) ≈
            QET.choi_matrix(QET.dephasing_channel(3, 0.25))
        @test size(CompatC.DepolarizingChannel(2, -0.5)) == (4, 4)
        @test size(CompatC.DephasingChannel(2, -2)) == (4, 4)
        @test CompatC.PauliChannel([0.0, 1.0, 0.0, 0.0]) ≈
            QET.choi_matrix(QET.pauli_channel([0.0, 1.0, 0.0, 0.0]))
        @test CompatC.ChoiMap() == QET.choi_matrix(QET.choi_map())
        @test CompatC.ReductionMap(3, 2) == QET.choi_matrix(QET.reduction_map(3, 2))

        Random.seed!(0x51c0)
        expected_global = rand(UInt)
        Random.seed!(0x51c0)
        local_rng = MersenneTwister(0xacc0)
        random_pauli = CompatC.PauliChannel(local_rng, 1)
        @test size(random_pauli) == (4, 4)
        @test rand(UInt) == expected_global
        @test !hasmethod(CompatC.PauliChannel, Tuple{Int})

        two_sided = reshape([identity2, bit_flip], 1, 2)
        @test CompatC.ApplyMap(density, two_sided) ≈ identity2 * density * adjoint(bit_flip)
        reduction_raw = Matrix{Float64}(CompatC.ReductionMap(2))
        reduction_factors = CompatC.KrausOperators(reduction_raw)
        @test reduction_factors isa Matrix
        @test QET.choi_matrix(
            QET.OperatorSumRepresentation(
                [reduction_factors[index, 1] for index in axes(reduction_factors, 1)],
                [reduction_factors[index, 2] for index in axes(reduction_factors, 1)],
            ),
        ) ≈ reduction_raw
        @test_throws ArgumentError CompatC.ChoiMatrix(kraus, 0)
        @test_throws DimensionMismatch CompatC.KrausOperators(
            CompatC.ChoiMatrix(rectangular), (2, 2)
        )
        @test_throws DimensionMismatch CompatC.PartialMap(
            product, raw_dephasing, 2, [2 2; 1 3]
        )
        @test_throws ArgumentError CompatC.KrausOperators(
            raw_choi, _TierCZeroBasedMatrix([2 2; 2 2])
        )
        @test_throws ArgumentError CompatC.PartialMap(
            product, raw_dephasing, 2, _TierCZeroBasedMatrix([2 2; 2 2])
        )
        @test_throws ArgumentError CompatC.PauliChannel(
            _TierCZeroBasedVector([1.0, 0.0, 0.0, 0.0])
        )
    end
end
