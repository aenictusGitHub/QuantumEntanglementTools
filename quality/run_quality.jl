using Aqua
using JET
using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const QET = QuantumEntanglementTools

@testset "Aqua package quality" begin
    Aqua.test_all(QET)
end

@testset "representative JET inference checks" begin
    density = Matrix{ComplexF64}(I, 4, 4) / 4
    permutation_plan = QET.SubsystemPermutationPlan((2, 3, 2), (3, 1, 2))
    trace_plan = QET.PartialTracePlan((2, 2), (2,))
    transpose_plan = QET.PartialTransposePlan((2, 2), (1,))
    realignment_plan = QET.RealignmentPlan((2, 2))
    identity_channel = QET.KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    bell_density = bell * bell'
    matrix_analysis_input = Float64[
        1 2 3
        4 5 6
        7 8 10
    ]
    JET.@test_opt target_modules = (QET,) QET.tensor_product(
        ComplexF64[1, 0], ComplexF64[0, 1]
    )
    JET.@test_opt target_modules = (QET,) QET.permute_subsystems(
        collect(1:12), permutation_plan
    )
    JET.@test_opt target_modules = (QET,) QET.partial_trace(density, trace_plan)
    JET.@test_opt target_modules = (QET,) QET.partial_transpose(density, transpose_plan)
    JET.@test_opt target_modules = (QET,) QET.realign(density, realignment_plan)
    JET.@test_opt target_modules = (QET,) QET.fourier_matrix(4)
    JET.@test_opt target_modules = (QET,) QET.isotropic_state(3, 0.2; sparse_output=false)
    JET.@test_opt target_modules = (QET,) QET.random_state_vector(Xoshiro(0x514554), 8)
    JET.@test_opt target_modules = (QET,) QET.apply_channel(
        density[1:2, 1:2], identity_channel
    )
    JET.@test_opt target_modules = (QET,) QET.choi_representation(identity_channel)
    JET.@test_opt target_modules = (QET,) QET.trace_norm(density)
    JET.@test_opt target_modules = (QET,) QET.concurrence(bell)
    JET.@test_opt target_modules = (QET,) QET.ppt_criterion(
        bell_density, (2, 2); atol=1e-12, rtol=0.0
    )
    JET.@test_opt target_modules = (QET,) QET.l1_coherence(bell)
    JET.@test_opt target_modules = (QET,) QET.relative_entropy_coherence(bell; base=2)
    JET.@test_opt target_modules = (QET,) QET.coherence_rank(bell)
    JET.@test_opt target_modules = (QET,) QET.operator_schmidt_coefficients(
        bell_density, (2, 2)
    )
    JET.@test_opt target_modules = (QET,) QET.is_product_vector(bell, (2, 2))
    JET.@test_opt target_modules = (QET,) QET.entanglement_of_formation(bell, (2, 2))
    JET.@test_opt target_modules = (QET,) QET.in_separable_ball(density, (2, 2))
    JET.@test_opt target_modules = (QET,) QET.majorizes([4.0, 1.0, 1.0], [3.0, 2.0, 1.0])
    JET.@test_opt target_modules = (QET,) QET.elementary_symmetric_polynomial(
        [1.0, 2.0, 3.0, 4.0], 2
    )
    JET.@test_opt target_modules = (QET,) QET.compound_matrix(matrix_analysis_input, 2)
    JET.@test_opt target_modules = (QET,) QET.additive_compound_matrix(
        matrix_analysis_input, 2
    )
end
