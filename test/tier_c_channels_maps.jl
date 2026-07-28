using LinearAlgebra
using Random
using SparseArrays

const QET = QuantumEntanglementTools
const CompatC = QuantumEntanglementTools.MATLABCompat

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

        rectangular_super = QET.SuperoperatorRepresentation(zeros(9, 4))
        @test QET.input_dimension(rectangular_super) == 2
        @test QET.output_dimension(rectangular_super) == 3
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

        recovered = QET.kraus_representation(choi)
        @test QET.choi_matrix(recovered) ≈ QET.choi_matrix(choi)
        @test QET.superoperator_matrix(QET.kraus_representation(superoperator)) ≈
            QET.superoperator_matrix(superoperator)

        input = ComplexF64[0.6 0.2im; -0.2im 0.4]
        expected = sum(operator * input * operator' for operator in operators)
        @test QET.apply_channel(input, kraus) ≈ expected
        @test QET.apply_channel(input, choi) ≈ expected
        @test QET.apply_channel(input, superoperator) ≈ expected
        @test QET.apply_channel(input, operators) ≈ expected

        raw_choi = QET.choi_matrix(choi)
        raw_choi[1, 1] = 99
        @test QET.choi_matrix(choi)[1, 1] != 99

        @test_throws DimensionMismatch QET.apply_channel(ones(3, 3), kraus)
        @test_throws ArgumentError QET.apply_channel(fill(NaN, 2, 2), kraus)

        zero_map = QET.ChoiRepresentation(zeros(4, 4))
        zero_kraus = QET.kraus_representation(zero_map)
        @test length(zero_kraus.operators) == 1
        @test iszero(only(zero_kraus.operators))
        @test iszero(QET.choi_matrix(zero_kraus))
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
        @test QET.is_completely_positive(map)
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
        @test QET.is_completely_positive(identity_map)
        @test QET.is_trace_preserving(identity_map)
        @test QET.is_unital(identity_map)

        transpose_super = zeros(Float64, 4, 4)
        for column in 1:2, row in 1:2
            transpose_super[column + (row - 1) * 2, row + (column - 1) * 2] = 1
        end
        transpose_map = QET.SuperoperatorRepresentation(transpose_super, 2, 2)
        @test !QET.is_completely_positive(transpose_map; atol=1e-14, rtol=1e-14)
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
        @test !QET.is_completely_positive(nonhermitian)
        @test_throws DomainError QET.kraus_representation(nonhermitian)

        scaled_identity = QET.KrausRepresentation([2.0 * Matrix{Float64}(I, 2, 2)])
        @test QET.is_completely_positive(scaled_identity)
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
        @test QET.is_completely_positive(completely_depolarizing)
        @test QET.is_trace_preserving(completely_depolarizing)
        @test QET.is_unital(completely_depolarizing)
        @test issparse(QET.choi_matrix(completely_depolarizing))

        identity_depolarizing = QET.depolarizing_channel(2, 1)
        @test QET.apply_channel(density, identity_depolarizing) ≈ density
        @test QET.is_completely_positive(QET.depolarizing_channel(2, -1 // 3))
        @test_throws DomainError QET.depolarizing_channel(2, -0.34)
        @test_throws DomainError QET.depolarizing_channel(2, 1.01)

        complete_dephasing = QET.dephasing_channel(2)
        @test QET.apply_channel(density, complete_dephasing) ≈ Diagonal(diag(density))
        @test QET.apply_channel(density, QET.dephasing_channel(2, 1)) ≈ density
        @test QET.is_completely_positive(QET.dephasing_channel(3, -1 // 2))
        @test_throws DomainError QET.dephasing_channel(3, -0.51)

        bit_flip = QET.pauli_channel([0.0, 1.0, 0.0, 0.0])
        x = ComplexF64[0 1; 1 0]
        @test QET.apply_channel(density, bit_flip) ≈ x * density * x
        @test QET.is_completely_positive(bit_flip)
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
        @test QET.is_completely_positive(complement)
        @test QET.is_trace_preserving(complement)
        choi_complement = QET.complementary_channel(QET.choi_representation(map))
        @test QET.input_dimension(choi_complement) == 2
        @test QET.output_dimension(choi_complement) == 2
        @test QET.is_completely_positive(choi_complement)
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
        @test !QET.is_completely_positive(reduction; atol=1e-14, rtol=1e-14)

        generalized_reduction = QET.reduction_map(3, 2)
        matrix3 = reshape(ComplexF64.(1:9), 3, 3)
        @test QET.apply_channel(matrix3, generalized_reduction) ≈
            2 * tr(matrix3) * Matrix{ComplexF64}(I, 3, 3) - matrix3

        choi = QET.choi_map()
        @test size(QET.choi_matrix(choi)) == (9, 9)
        @test !QET.is_completely_positive(choi; atol=1e-14, rtol=1e-14)
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
        @test QET.choi_matrix(QET.KrausRepresentation(canonical)) ≈ raw_choi
        @test length(CompatC.KrausOperators(CompatC.ChoiMatrix(rectangular), (2, 3))) == 1

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
        @test_throws ArgumentError CompatC.ApplyMap(density, two_sided)
        @test_throws DomainError CompatC.KrausOperators(CompatC.ReductionMap(2))
        @test_throws ArgumentError CompatC.ChoiMatrix(kraus, 0)
        @test_throws DimensionMismatch CompatC.KrausOperators(
            CompatC.ChoiMatrix(rectangular), (2, 2)
        )
        @test_throws ArgumentError CompatC.PartialMap(product, raw_dephasing, 2, [2 2; 1 4])
    end
end
