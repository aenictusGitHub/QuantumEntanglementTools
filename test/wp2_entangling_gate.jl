using LinearAlgebra
using SparseArrays
using Test

const EGCompat = QuantumEntanglementTools.MATLABCompat

function _entangling_gate_reconstruction(gate, result)
    local_product = tensor_product(result.local_factors...)
    permutation = permutation_operator(
        result.factor_output_dims,
        result.subsystem_permutation;
        T=eltype(gate),
        sparse_output=false,
    )
    return permutation * local_product
end

@testset "WP2 entangling-gate certificates" begin
    @testset "nonentangling local and permutation gates" begin
        identity_gate = Matrix{ComplexF64}(I, 4, 4)
        identity_result = is_entangling_gate(identity_gate)
        @test identity_result isa EntanglingGateResult
        @test identity_result.status === :not_entangling
        @test identity_result.product_witness === nothing
        @test identity_result.certificate_analysis.status === :within_tolerance
        @test identity_result.input_dims == (2, 2)
        @test identity_result.output_dims == (2, 2)
        @test isapprox(
            _entangling_gate_reconstruction(identity_gate, identity_result),
            identity_gate;
            atol=1e-13,
            rtol=1e-13,
        )

        hadamard = [1.0 1.0; 1.0 -1.0] / sqrt(2)
        phase = ComplexF64[1 0; 0 im]
        local_gate = tensor_product(hadamard, phase)
        local_copy = copy(local_gate)
        local_result = is_entangling_gate(local_gate, (2, 2))
        @test local_result.status === :not_entangling
        @test local_gate == local_copy
        @test isapprox(
            _entangling_gate_reconstruction(local_gate, local_result),
            local_gate;
            atol=1e-13,
            rtol=1e-13,
        )

        swap_gate = Matrix(swap_operator((2, 2)))
        swap_result = is_entangling_gate(swap_gate, (2, 2))
        @test swap_result.status === :not_entangling
        @test swap_result.permutations_checked == 2
        @test swap_result.subsystem_permutation == (2, 1)
        @test isapprox(
            _entangling_gate_reconstruction(swap_gate, swap_result),
            swap_gate;
            atol=1e-13,
            rtol=1e-13,
        )

        rectangular_plan = SubsystemPermutationPlan((2, 3), (2, 1))
        rectangular_swap = Matrix(permutation_operator(rectangular_plan; T=ComplexF64))
        rectangular_dims = [3 2; 2 3]
        rectangular_result = is_entangling_gate(rectangular_swap, rectangular_dims)
        @test rectangular_result.status === :not_entangling
        @test rectangular_result.input_dims == (2, 3)
        @test rectangular_result.output_dims == (3, 2)
        @test isapprox(
            _entangling_gate_reconstruction(rectangular_swap, rectangular_result),
            rectangular_swap;
            atol=2e-13,
            rtol=2e-13,
        )
    end

    @testset "explicit entangling witnesses" begin
        cnot = ComplexF64[
            1 0 0 0
            0 1 0 0
            0 0 0 1
            0 0 1 0
        ]
        cnot_result = is_entangling_gate(cnot, (2, 2))
        @test cnot_result.status === :entangling
        @test cnot_result.subsystem_permutation === nothing
        @test cnot_result.local_factors === nothing
        @test cnot_result.product_witness isa SparseVector
        @test nnz(cnot_result.product_witness) <= 2
        @test norm(cnot_result.product_witness) ≈ 1
        @test is_product_vector(Vector(cnot_result.product_witness), (2, 2)).status ===
            :within_tolerance
        @test cnot_result.certificate_analysis.status === :outside_tolerance
        @test is_product_vector(cnot * cnot_result.product_witness, (2, 2)).status ===
            :outside_tolerance
        @test cnot_result.permutations_checked == 2
        @test 1 <= cnot_result.witnesses_checked <= 12

        controlled_z = Diagonal(ComplexF64[1, 1, 1, -1]) |> Matrix
        cz_result = is_entangling_gate(controlled_z, 2)
        @test cz_result.status === :entangling
        @test nnz(cz_result.product_witness) <= 4
        @test norm(cz_result.product_witness) ≈ 1
        @test cz_result.certificate_analysis.status === :outside_tolerance

        cnot32 = Float32.(real(cnot))
        result32 = is_entangling_gate(cnot32, (2, 2))
        @test result32.status === :entangling
        @test eltype(result32.product_witness) == Float32
    end

    @testset "compatibility, sparse policy, and guards" begin
        cnot = [
            1.0 0 0 0
            0 1 0 0
            0 0 0 1
            0 0 1 0
        ]
        compatibility = EGCompat.IsEntanglingGate(cnot)
        @test compatibility isa EntanglingGateResult
        @test compatibility.status === :entangling

        sparse_cnot = sparse(cnot)
        @test_throws ArgumentError is_entangling_gate(sparse_cnot)
        @test is_entangling_gate(sparse_cnot, (2, 2); allow_densify=true).status ===
            :entangling

        @test_throws ArgumentError is_entangling_gate(cnot, (2, 2); max_permutations=1)
        @test_throws ArgumentError is_entangling_gate(cnot, (2, 2); max_witnesses=11)
        @test_throws ArgumentError is_entangling_gate(cnot, (2, 2); max_work=223)
        @test is_entangling_gate(
            cnot, (2, 2); max_permutations=nothing, max_witnesses=nothing, max_work=nothing
        ).status === :entangling
    end

    @testset "invalid inputs and boundary retention" begin
        @test_throws DimensionMismatch is_entangling_gate(ones(2, 3), (2, 2))
        @test_throws DomainError is_entangling_gate(2.0 .* Matrix(I, 4, 4))
        @test_throws ArgumentError is_entangling_gate(Matrix(I, 6, 6))
        @test_throws DimensionMismatch is_entangling_gate(Matrix(I, 6, 6), 4)
        @test_throws DimensionMismatch is_entangling_gate(Matrix(I, 4, 4), (2, 3))
        @test_throws ArgumentError is_entangling_gate(Matrix(I, 4, 4), (4,))
        @test_throws DimensionMismatch is_entangling_gate(Matrix(I, 4, 4), [2, 2, 2, 2])
        @test_throws ArgumentError is_entangling_gate(
            Matrix(I, 4, 4), (2, 2); max_permutations=-1
        )
        @test_throws ArgumentError is_entangling_gate(
            Matrix(I, 4, 4), (2, 2); max_witnesses=true
        )
        @test_throws ArgumentError is_entangling_gate(Matrix{BigFloat}(I, 4, 4), (2, 2))

        scaled = (1 + 1e-8) .* Matrix{Float64}(I, 4, 4)
        residual = norm(scaled' * scaled - Matrix{Float64}(I, 4, 4))
        boundary = is_entangling_gate(scaled, (2, 2); atol=residual, rtol=0)
        @test boundary.status === :unknown
        @test boundary.permutations_checked == 0
        @test boundary.witnesses_checked == 0
        @test boundary.unitarity_residual == boundary.unitarity_tolerance
    end
end
