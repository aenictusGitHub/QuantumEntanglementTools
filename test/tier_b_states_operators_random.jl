using LinearAlgebra
using Random
using SparseArrays

const CompatB = QuantumEntanglementTools.MATLABCompat

function _minimum_eigenvalue(matrix)
    return eigmin(Hermitian(Matrix(matrix)))
end

@testset "Tier B operators" begin
    @testset "Pauli operators" begin
        identity2 = [1 0; 0 1]
        x = [0 1; 1 0]
        y = [0 -im; im 0]
        z = [1 0; 0 -1]

        @test pauli(0) == identity2
        @test pauli('X') == x
        @test pauli("y") == y
        @test pauli(:Z) == z
        @test pauli([:X, :Z]) == kron(x, z)
        @test pauli("XYZ") == kron(kron(x, y), z)
        @test issparse(pauli((1, 2); sparse_output=true))
        @test pauli(:X)^2 == identity2
        @test pauli(:Y)^2 == identity2
        @test pauli(:Z)^2 == identity2
        @test pauli(:X) * pauli(:Y) == im * pauli(:Z)

        @test_throws ArgumentError pauli(4)
        @test_throws ArgumentError pauli(:A)
        @test_throws ArgumentError pauli(())
        @test_throws ArgumentError pauli(true)
    end

    @testset "generalized Pauli operators" begin
        for dimension in 1:5
            for shift in 0:(dimension - 1), clock in 0:(dimension - 1)
                operator = generalized_pauli(shift, clock, dimension)
                @test operator' * operator ≈ I atol = 2e-14
                @test count(!iszero, operator) == dimension
                sparse_operator = generalized_pauli(
                    shift, clock, dimension; sparse_output=true
                )
                @test issparse(sparse_operator)
                @test sparse_operator == operator
                @test nnz(sparse_operator) == dimension
            end
        end

        shift = generalized_pauli(1, 0, 3)
        clock = generalized_pauli(0, 1, 3)
        omega = cis(2π / 3)
        @test clock * shift ≈ omega * shift * clock atol = 2e-14
        @test generalized_pauli(0, 0, 4) ≈ Matrix{ComplexF64}(I, 4, 4)
        @test eltype(generalized_pauli(1, 2, 3; T=BigFloat)) == Complex{BigFloat}
        @test_throws ArgumentError generalized_pauli(-1, 0, 3)
        @test_throws ArgumentError generalized_pauli(0, 3, 3)
        @test_throws ArgumentError generalized_pauli(0, 0, 0)
    end

    @testset "Gell-Mann bases" begin
        @test gell_mann(0) == Matrix{Float64}(I, 3, 3)
        @test gell_mann(1) == [0.0 1 0; 1 0 0; 0 0 0]
        @test gell_mann(2) == [0 -im 0; im 0 0; 0 0 0]
        @test gell_mann(3) == [1.0 0 0; 0 -1 0; 0 0 0]
        @test gell_mann(8) ≈ Diagonal([1 / sqrt(3), 1 / sqrt(3), -2 / sqrt(3)])

        matrices = [gell_mann(index) for index in 1:8]
        for first in 1:8, second in 1:8
            inner_product = tr(matrices[first]' * matrices[second])
            @test inner_product ≈ (first == second ? 2 : 0) atol = 2e-14
        end
        for matrix in matrices
            @test ishermitian(matrix)
            @test tr(matrix) ≈ 0 atol = 2e-14
        end

        mapping = ((0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (2, 0), (1, 2), (2, 1), (2, 2))
        for index in 0:8
            @test gell_mann(index) == generalized_gell_mann(mapping[index + 1]..., 3)
        end
        for dimension in 2:5, first in 0:(dimension - 1), second in 0:(dimension - 1)
            matrix = generalized_gell_mann(first, second, dimension)
            @test ishermitian(matrix)
            if (first, second) != (0, 0)
                @test tr(matrix) ≈ 0 atol = 2e-14
                @test tr(matrix' * matrix) ≈ 2 atol = 2e-14
            end
        end
        @test issparse(gell_mann(7; sparse_output=true))
        @test eltype(gell_mann(8; T=BigFloat)) == BigFloat
        @test eltype(generalized_gell_mann(2, 1, 4; T=BigFloat)) == Complex{BigFloat}
        @test_throws ArgumentError gell_mann(9)
        @test_throws ArgumentError generalized_gell_mann(0, 3, 3)
    end

    @testset "Fourier matrix" begin
        for dimension in 1:8
            fourier = fourier_matrix(dimension)
            @test fourier' * fourier ≈ I atol = 3e-14
            @test fourier[1, :] ≈ fill(inv(sqrt(dimension)), dimension)
            @test fourier[:, 1] ≈ fill(inv(sqrt(dimension)), dimension)
        end
        high_precision = fourier_matrix(3; T=BigFloat)
        @test eltype(high_precision) == Complex{BigFloat}
        @test norm(high_precision' * high_precision - I) < big"1e-70"
        @test_throws ArgumentError fourier_matrix(0)
        @test_throws ArgumentError fourier_matrix(true)
    end
end

@testset "Tier B named states" begin
    @testset "pure-state constructors" begin
        phi = maximally_entangled(3)
        @test phi == [inv(sqrt(3)), 0, 0, 0, inv(sqrt(3)), 0, 0, 0, inv(sqrt(3))]
        @test norm(phi) ≈ 1
        @test maximally_entangled(3; normalized=false, sparse_output=true) ==
            sparsevec([1, 5, 9], ones(3), 9)
        @test eltype(maximally_entangled(2; T=BigFloat)) == BigFloat

        expected_bells = (
            [1, 0, 0, 1] / sqrt(2),
            [1, 0, 0, -1] / sqrt(2),
            [0, 1, 1, 0] / sqrt(2),
            [0, 1, -1, 0] / sqrt(2),
        )
        bells = [bell_state(index) for index in 0:3]
        @test Tuple(bells) == expected_bells
        @test hcat(bells...)' * hcat(bells...) ≈ I atol = 2e-14
        @test bell_state(3; normalized=false, sparse_output=true) ==
            sparsevec([2, 3], [1.0, -1.0], 4)

        ghz = ghz_state(2, 3)
        @test issparse(ghz)
        @test ghz == sparsevec([1, 8], fill(inv(sqrt(2)), 2), 8)
        @test ghz_state(3, 2; coefficients=[1, 2im, -3], sparse_output=false) ==
            [1, 0, 0, 0, 2im, 0, 0, 0, -3]
        @test norm(ghz_state(5, 4)) ≈ 1

        w = w_state(4)
        @test w == sparsevec([9, 5, 3, 2], fill(0.5, 4), 16)
        @test w_state(3; coefficients=[1, 2, 3], sparse_output=false) ==
            [0, 3, 2, 0, 1, 0, 0, 0]
        @test dicke_state(4, 1) == w
        @test nnz(dicke_state(6, 3)) == binomial(6, 3)
        @test norm(dicke_state(6, 3)) ≈ 1
        @test dicke_state(4, 0) == sparsevec([1], [1.0], 16)
        @test dicke_state(4, 4) == sparsevec([16], [1.0], 16)
        @test all(==(1), nonzeros(dicke_state(5, 2; normalized=false)))

        @test_throws ArgumentError maximally_entangled(0)
        @test_throws ArgumentError bell_state(4)
        @test_throws DimensionMismatch ghz_state(3, 2; coefficients=[1, 2])
        @test_throws ArgumentError w_state(1)
        @test_throws ArgumentError dicke_state(3, 4)
        @test_throws ArgumentError dicke_state(3, -1)
    end

    @testset "mixed-state families" begin
        isotropic_exact = isotropic_state(3, 1 // 4)
        @test eltype(isotropic_exact) == Rational{Int}
        @test tr(isotropic_exact) == 1
        @test ishermitian(isotropic_exact)
        @test _minimum_eigenvalue(isotropic_exact) >= -1e-14
        @test isotropic_state(2, 1; sparse_output=false) ≈ bell_state() * bell_state()'
        @test isotropic_state(1, -12) == sparse([1], [1], [1.0], 1, 1)
        @test eltype(werner_state(3, 1 // 4)) == Rational{Int}
        @test tr(werner_state(3, 1 // 4)) == 1

        for dimension in 2:4, alpha in (-1.0, -0.2, 0.0, 0.7, 1.0)
            werner = werner_state(dimension, alpha)
            @test tr(werner) ≈ 1 atol = 2e-14
            @test ishermitian(werner)
            @test _minimum_eigenvalue(werner) >= -2e-14
            @test partial_trace(werner, (dimension, dimension); trace_out=2) ≈
                Matrix(I, dimension, dimension) / dimension
        end

        for dims in ((3, 3), (2, 4)), parameter in (0.0, 0.3, 1.0)
            horodecki = horodecki_state(parameter; dims=dims)
            @test size(horodecki) == (prod(dims), prod(dims))
            @test tr(horodecki) ≈ 1 atol = 2e-14
            @test ishermitian(horodecki)
            @test _minimum_eigenvalue(horodecki) >= -2e-14
            transposed = partial_transpose(horodecki, dims; systems=1)
            @test _minimum_eigenvalue(transposed) >= -2e-14
        end
        @test eltype(horodecki_state(big"0.5")) == BigFloat

        for weight in (0.0, 0.4, 1.0), angle in (-0.2, 0.7)
            gisin = gisin_state(weight, angle)
            @test tr(gisin) ≈ 1 atol = 2e-14
            @test ishermitian(gisin)
            @test _minimum_eigenvalue(gisin) >= -2e-14
        end

        for dimension in (2, 4, 6), weight in (0.0, 0.35, 1.0)
            breuer = breuer_state(dimension, weight)
            @test issparse(breuer)
            @test tr(breuer) ≈ 1 atol = 3e-14
            @test ishermitian(breuer)
            @test _minimum_eigenvalue(breuer) >= -3e-14
        end

        brauer = brauer_states(2, 2)
        @test issparse(brauer)
        @test size(brauer) == (16, 3)
        @test vec(sum(abs2, brauer; dims=1)) == fill(4.0, 3)
        @test sort(vec(Matrix(brauer' * brauer))) ==
            sort([4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
        @test size(brauer_states(3, 1)) == (9, 1)

        chessboard = chessboard_state(1, 2, 3, 4, 5, 6)
        @test size(chessboard) == (9, 9)
        @test tr(chessboard) ≈ 1 atol = 2e-14
        @test ishermitian(chessboard)
        @test _minimum_eigenvalue(chessboard) >= -2e-14
        complex_chessboard = chessboard_state(
            1 + im, 2 - im, 3 + 2im, -1 + im, 2, 1 - im; s=0.4im, t=-0.7
        )
        @test tr(complex_chessboard) ≈ 1 atol = 2e-14
        @test ishermitian(complex_chessboard)
        @test _minimum_eigenvalue(complex_chessboard) >= -2e-14

        @test_throws ArgumentError isotropic_state(2, -1)
        @test_throws ArgumentError isotropic_state(2, 1.1)
        @test_throws ArgumentError werner_state(2, 1.1)
        @test_throws ArgumentError werner_state(1, 0)
        @test_throws ArgumentError horodecki_state(-0.1)
        @test_throws ArgumentError horodecki_state(0.2; dims=(2, 3))
        @test_throws ArgumentError gisin_state(1.1, 0)
        @test_throws ArgumentError gisin_state(0.5, Inf)
        @test_throws ArgumentError breuer_state(3, 0.5)
        @test_throws ArgumentError breuer_state(4, -0.1)
        @test_throws ArgumentError brauer_states(2, 0)
        @test_throws ArgumentError chessboard_state(0, 0, 0, 0, 0, 1)
        @test_throws ArgumentError chessboard_state(0, 0, 0, 0, 1, 0)
        @test_throws ArgumentError chessboard_state(1, 2, 3, 4, 5, Inf)
        @test_throws ArgumentError chessboard_state(1, 2, 3, 4, 5, 6; s=NaN)
    end
end

@testset "Tier B explicit-RNG constructors" begin
    @testset "probabilities and state vectors" begin
        probabilities = random_probabilities(Xoshiro(11), 20)
        @test length(probabilities) == 20
        @test all(>(0), probabilities)
        @test sum(probabilities) ≈ 1
        @test random_probabilities(Xoshiro(11), 20) == probabilities
        @test random_probabilities(Xoshiro(8), 1) == [1.0]

        state = random_state_vector(Xoshiro(12), 9)
        @test length(state) == 9
        @test norm(state) ≈ 1
        @test random_state_vector(Xoshiro(12), 9) == state
        real_state = random_state_vector(Xoshiro(13), (2, 3); real=true)
        @test eltype(real_state) == Float64
        @test length(real_state) == 6
        @test norm(real_state) ≈ 1

        ranked = random_state_vector(Xoshiro(14), (4, 5); schmidt_rank=2)
        coefficient_matrix = Matrix(transpose(reshape(ranked, 5, 4)))
        @test rank(coefficient_matrix; atol=1e-11) == 2
        @test norm(ranked) ≈ 1
        @test eltype(random_state_vector(Xoshiro(15), 4; real=true, schmidt_rank=3)) ==
            Float64

        @test_throws ArgumentError random_probabilities(Xoshiro(1), 0)
        @test_throws DimensionMismatch random_state_vector(Xoshiro(1), (2, 3, 4))
        @test_throws ArgumentError random_state_vector(Xoshiro(1), (2, 3); schmidt_rank=3)
        @test_throws ArgumentError random_state_vector(Xoshiro(1), 3; schmidt_rank=0)
    end

    @testset "unitaries and density matrices" begin
        unitary = random_unitary(Xoshiro(20), 8)
        @test unitary' * unitary ≈ I atol = 3e-14
        @test random_unitary(Xoshiro(20), 8) == unitary
        orthogonal = random_unitary(Xoshiro(21), 7; real=true)
        @test eltype(orthogonal) == Float64
        @test orthogonal' * orthogonal ≈ I atol = 3e-14

        for real_output in (false, true), rank_bound in (1, 3, 5)
            rho = random_density_matrix(
                Xoshiro(22 + rank_bound + real_output), 5; real=real_output, rank=rank_bound
            )
            @test tr(rho) ≈ 1 atol = 3e-14
            @test ishermitian(rho)
            @test _minimum_eigenvalue(rho) >= -3e-14
            @test rank(rho; atol=2e-11) == rank_bound
            @test real_output ? eltype(rho) == Float64 : eltype(rho) == ComplexF64
        end
        bures = random_density_matrix(Xoshiro(30), 4; distribution=:bures)
        @test tr(bures) ≈ 1 atol = 3e-14
        @test _minimum_eigenvalue(bures) >= -3e-14
        @test random_density_matrix(Xoshiro(31), 3; distribution="HAAR") ==
            random_density_matrix(Xoshiro(31), 3; distribution=:hilbert_schmidt)

        @test_throws ArgumentError random_unitary(Xoshiro(1), 0)
        @test_throws ArgumentError random_density_matrix(Xoshiro(1), 3; rank=4)
        @test_throws ArgumentError random_density_matrix(
            Xoshiro(1), 3; distribution=:unknown
        )
    end

    @testset "graphs and POVMs" begin
        graph = random_graph(Xoshiro(40), 30; edge_probability=0.3)
        @test graph isa BitMatrix
        @test graph == graph'
        @test all(iszero, diag(graph))
        @test random_graph(Xoshiro(40), 30; edge_probability=0.3) == graph
        @test !any(random_graph(Xoshiro(2), 5; edge_probability=0))
        @test random_graph(Xoshiro(2), 5; edge_probability=1) == .!Matrix{Bool}(I, 5, 5)

        for real_output in (false, true), outcomes in (1, 2, 5)
            effects = random_povm(
                Xoshiro(50 + outcomes + real_output), 4, outcomes; real=real_output
            )
            @test length(effects) == outcomes
            @test sum(effects) ≈ I atol = 4e-14
            for effect in effects
                @test ishermitian(effect)
                @test _minimum_eigenvalue(effect) >= -4e-14
                @test size(effect) == (4, 4)
            end
        end
        @test random_povm(Xoshiro(60), 3, 4) == random_povm(Xoshiro(60), 3, 4)

        @test_throws ArgumentError random_graph(Xoshiro(1), 3; edge_probability=-0.1)
        @test_throws ArgumentError random_graph(Xoshiro(1), 3; edge_probability=NaN)
        @test_throws ArgumentError random_povm(Xoshiro(1), 2, 0)
    end

    @testset "global stream isolation" begin
        Random.seed!(0x5127)
        reference_first = rand()
        reference_second = rand()

        Random.seed!(0x5127)
        observed_first = rand()
        rng = Xoshiro(0x88)
        random_probabilities(rng, 5)
        random_state_vector(rng, (3, 4); schmidt_rank=2)
        random_density_matrix(rng, 4; distribution=:bures)
        random_unitary(rng, 4)
        random_graph(rng, 8)
        random_povm(rng, 3, 4)
        CompatB.RandomProbabilities(rng, 3)
        CompatB.RandomStateVector(rng, 4)
        CompatB.RandomDensityMatrix(rng, 3)
        CompatB.RandomUnitary(rng, 3)
        CompatB.RandomGraph(rng, 4)
        CompatB.RandomPOVM(rng, 2, 3)
        observed_second = rand()

        @test observed_first == reference_first
        @test observed_second == reference_second
        @test !hasmethod(random_unitary, Tuple{Int})
        @test !hasmethod(random_probabilities, Tuple{Int})
    end
end

@testset "Tier B MATLAB compatibility wrappers" begin
    @test CompatB.Pauli("XZ", 0) == pauli("XZ")
    @test issparse(CompatB.Pauli("XZ"))
    @test CompatB.GenPauli(1, 2, 3) == generalized_pauli(1, 2, 3)
    @test CompatB.GellMann(8) == gell_mann(8)
    @test CompatB.GenGellMann(2, 1, 3) == generalized_gell_mann(2, 1, 3)
    @test CompatB.FourierMatrix(4) == fourier_matrix(4)
    @test CompatB.MaxEntangled(3, 1, 0) ==
        maximally_entangled(3; sparse_output=true, normalized=false)
    @test CompatB.Bell(5) == bell_state(1)
    @test CompatB.GHZState(2, 3) == ghz_state(2, 3)
    @test CompatB.WState(4) == w_state(4)
    @test CompatB.DickeState(5, 2, 0) == dicke_state(5, 2; normalized=false)
    @test CompatB.IsotropicState(3, 0.2) == isotropic_state(3, 0.2)
    @test CompatB.WernerState(3, 0.2) == werner_state(3, 0.2)
    @test CompatB.HorodeckiState(0.2, [2, 4]) == horodecki_state(0.2; dims=(2, 4))
    @test CompatB.GisinState(0.2, 0.3) == gisin_state(0.2, 0.3)
    @test CompatB.BreuerState(4, 0.2) == breuer_state(4, 0.2)
    @test CompatB.BrauerStates(2, 2) == brauer_states(2, 2)
    @test CompatB.ChessboardState(1, 2, 3, 4, 5, 6) == chessboard_state(1, 2, 3, 4, 5, 6)

    rng_first = Xoshiro(77)
    rng_second = Xoshiro(77)
    @test CompatB.RandomProbabilities(rng_first, 4) == random_probabilities(rng_second, 4)
    @test CompatB.RandomStateVector(rng_first, 5, 1) ==
        random_state_vector(rng_second, 5; real=true)
    @test CompatB.RandomDensityMatrix(rng_first, 4, 0, 2, "hs") ==
        random_density_matrix(rng_second, 4; rank=2)
    @test CompatB.RandomUnitary(rng_first, 4, 1) == random_unitary(rng_second, 4; real=true)
    @test CompatB.RandomGraph(rng_first, 4, 0.2) ==
        Float64.(random_graph(rng_second, 4; edge_probability=0.2))
    @test CompatB.RandomPOVM(rng_first, 3, 2, 1) == random_povm(rng_second, 3, 2; real=true)

    @test_throws ArgumentError CompatB.WernerState(2, [0.1, 0.2])
    @test_throws ArgumentError CompatB.RandomStateVector(Xoshiro(1), 3, 0, -1)
    @test_throws ArgumentError CompatB.RandomGraph(Xoshiro(1), 3, 2)
end
