using LinearAlgebra
using SparseArrays

const QETCoherence = QuantumEntanglementTools

@testset "Tier E coherence measures" begin
    @testset "l1 coherence" begin
        basis_state = ComplexF64[1, 0, 0, 0]
        plus_state = fill(ComplexF64(1 / 2), 4)
        @test QETCoherence.l1_coherence(basis_state) == 0.0
        @test QETCoherence.l1_coherence(plus_state) ≈ 3.0
        @test QETCoherence.l1_coherence(im * plus_state) ≈ 3.0

        incoherent = Diagonal([0.1, 0.2, 0.3, 0.4])
        coherent = plus_state * plus_state'
        @test QETCoherence.l1_coherence(incoherent) == 0.0
        @test QETCoherence.l1_coherence(coherent) ≈ 3.0
        @test QETCoherence.l1_coherence(Float32[inv(sqrt(2.0f0)), inv(sqrt(2.0f0))]) isa
            Float32
        big_plus = BigFloat[1, 1] / sqrt(big(2))
        @test QETCoherence.l1_coherence(big_plus) isa BigFloat
        @test QETCoherence.l1_coherence(big_plus) ≈ one(BigFloat)

        sparse_basis = sparsevec([1], ComplexF64[1], 4)
        @test QETCoherence.l1_coherence(sparse_basis) == 0.0
        sparse_mixed = sparse(Matrix{Float64}(I, 4, 4) / 4)
        @test_throws ArgumentError QETCoherence.l1_coherence(sparse_mixed)
        @test QETCoherence.l1_coherence(sparse_mixed; allow_densify=true) == 0.0

        @test_throws ArgumentError QETCoherence.l1_coherence(0.9basis_state)
        @test_throws DimensionMismatch QETCoherence.l1_coherence(ones(2, 3))
        @test_throws ArgumentError QETCoherence.l1_coherence(ComplexF64[1, NaN])
    end

    @testset "relative entropy of coherence" begin
        plus_qubit = ComplexF64[1, 1] / sqrt(2)
        basis_qubit = ComplexF64[1, 0]
        @test QETCoherence.relative_entropy_coherence(plus_qubit; base=2) ≈ 1.0
        @test QETCoherence.relative_entropy_coherence(basis_qubit; base=2) == 0.0

        plus_qutrit = fill(ComplexF64(inv(sqrt(3))), 3)
        @test QETCoherence.relative_entropy_coherence(plus_qutrit; base=3) ≈ 1.0

        diagonal = Diagonal([0.2, 0.3, 0.5])
        @test QETCoherence.relative_entropy_coherence(diagonal; base=2) ≈ 0.0 atol = 5e-15

        pure_density = plus_qubit * plus_qubit'
        @test QETCoherence.relative_entropy_coherence(pure_density; base=2) ≈ 1.0
        @test QETCoherence.relative_entropy_coherence(pure_density; base=exp(1)) ≈ log(2)
        big_plus = BigFloat[1, 1] / sqrt(big(2))
        big_entropy = QETCoherence.relative_entropy_coherence(big_plus; base=big(2))
        @test big_entropy isa BigFloat
        @test big_entropy ≈ one(BigFloat)

        @test_throws UndefKeywordError QETCoherence.relative_entropy_coherence(plus_qubit)
        @test_throws ArgumentError QETCoherence.relative_entropy_coherence(
            plus_qubit; base=1
        )
        @test_throws ArgumentError QETCoherence.relative_entropy_coherence(
            0.9plus_qubit; base=2
        )
    end

    @testset "coherence rank" begin
        basis_state = ComplexF64[1, 0]
        plus_state = ComplexF64[1, 1] / sqrt(2)
        @test QETCoherence.coherence_rank(basis_state) == 1
        @test QETCoherence.coherence_rank(plus_state) == 2
        @test QETCoherence.coherence_rank(BigFloat[1, 1] / sqrt(big(2))) == 2
        @test QETCoherence.coherence_rank(
            ComplexF64[sqrt(1 - 1e-14), 1e-7]; atol=1e-6, rtol=0
        ) == 1
        @test QETCoherence.coherence_rank(
            ComplexF64[sqrt(1 - 1e-14), 1e-7]; atol=1e-8, rtol=0
        ) == 2

        hadamard = ComplexF64[1 1; 1 -1] / sqrt(2)
        @test QETCoherence.coherence_rank(plus_state; basis=hadamard) == 1
        @test QETCoherence.coherence_rank(basis_state; basis=hadamard) == 2

        sparse_basis = sparsevec([1], ComplexF64[1], 2)
        @test QETCoherence.coherence_rank(sparse_basis) == 1
        @test_throws ArgumentError QETCoherence.coherence_rank(
            sparse_basis; basis=sparse(hadamard)
        )
        @test QETCoherence.coherence_rank(
            sparse_basis; basis=sparse(hadamard), allow_densify=true
        ) == 2

        @test_throws DimensionMismatch QETCoherence.coherence_rank(
            basis_state; basis=Matrix{Float64}(I, 3, 3)
        )
        @test_throws ArgumentError QETCoherence.coherence_rank(
            basis_state; basis=ones(2, 2)
        )
        @test_throws ArgumentError QETCoherence.coherence_rank(0.5basis_state)
        @test_throws ArgumentError QETCoherence.coherence_rank(basis_state; atol=-1)
    end

    @testset "MATLAB/QETLAB coherence compatibility" begin
        plus = ComplexF64[1, 1] / sqrt(2)
        plus_row = reshape(plus, 1, :)
        basis_qutrit = ComplexF64[1, 0, 0]
        @test QETCoherence.MATLABCompat.L1NormCoherence(plus) ≈ 1
        @test QETCoherence.MATLABCompat.L1NormCoherence(plus_row) ≈ 1
        @test QETCoherence.MATLABCompat.RelEntCoherence(plus) ≈ 1
        @test QETCoherence.MATLABCompat.RelEntCoherence(plus_row; base=exp(1)) ≈ log(2)
        @test QETCoherence.MATLABCompat.CoherenceRank(plus) == 2
        @test QETCoherence.MATLABCompat.CoherenceRank(plus_row) == 2
        @test QETCoherence.MATLABCompat.CoherenceRank(basis_qutrit) == 1

        almost_basis = ComplexF64[sqrt(1 - 1e-14), 1e-7]
        @test QETCoherence.MATLABCompat.CoherenceRank(almost_basis, 1e-6) == 1
        @test QETCoherence.MATLABCompat.CoherenceRank(almost_basis, 1e-8) == 2

        hadamard = ComplexF64[1 1; 1 -1] / sqrt(2)
        @test QETCoherence.MATLABCompat.CoherenceRank(plus, 1e-10, hadamard) == 1
        @test_throws DimensionMismatch QETCoherence.MATLABCompat.CoherenceRank(
            Matrix{Float64}(I, 2, 2)
        )
        @test_throws ArgumentError QETCoherence.MATLABCompat.CoherenceRank(plus, -1)
        @test_throws ArgumentError QETCoherence.MATLABCompat.RelEntCoherence(0.9plus)
    end
end
