using Test
using LinearAlgebra
using SparseArrays

using QuantumEntanglementTools: QuantumEntanglementTools
const QETPSDConstraint = QuantumEntanglementTools

if !isdefined(QETPSDConstraint, :positive_semidefinite_constraint)
    Base.include(
        QETPSDConstraint,
        joinpath(@__DIR__, "..", "src", "optimization", "psd_constraints.jl"),
    )
end

@testset "WP3 solver-neutral PSD constraint builder" begin
    variable = QETPSDConstraint.hermitian_variable(
        :rho, 2; variable_count=4, coefficient_type=Float64
    )
    constraint = QETPSDConstraint.positive_semidefinite_constraint(variable)

    @test constraint isa QETPSDConstraint.HermitianAffineMatrix{Float64}
    @test constraint.name === :rho
    @test constraint.dimension == 2
    @test constraint.variable_count == 4
    @test constraint.constant == variable.constant
    @test constraint.constant !== variable.constant
    @test length(constraint.terms) == length(variable.terms) == 4
    @test all(
        left.variable == right.variable && left.coefficient == right.coefficient for
        (left, right) in zip(constraint.terms, variable.terms)
    )
    @test all(
        left.coefficient !== right.coefficient for
        (left, right) in zip(constraint.terms, variable.terms)
    )

    coordinates = [1.0, 2.0, 0.5, -0.25]
    @test QETPSDConstraint.evaluate_affine(constraint, coordinates) ==
        QETPSDConstraint.evaluate_affine(variable, coordinates)

    variable.constant[1, 1] = 7
    variable.terms[1].coefficient[1, 1] = 9
    @test iszero(constraint.constant[1, 1])
    @test constraint.terms[1].coefficient[1, 1] == 1

    objective = QETPSDConstraint.AffineScalar(0.0, [1.0, 0.0, 0.0, 0.0])
    program = QETPSDConstraint.SemidefiniteProgram(
        :psd_constraint_smoke,
        :minimize,
        4,
        objective;
        psd_constraints=[constraint],
        primal_views=[constraint],
        initial_point=[1.0, 1.0, 0.0, 0.0],
        known_feasible_point=[1.0, 1.0, 0.0, 0.0],
    )
    @test length(program.psd_constraints) == 1
    @test program.psd_constraints[1].name === :rho
    @test QETPSDConstraint.primal_residual(
        program, program.known_feasible_point; allow_densify=true
    ) == 0

    unavailable = QETPSDConstraint.solve_optimization(program)
    @test unavailable.status === QETPSDConstraint.OptimizationBackendUnavailable
    @test unavailable.primal === nothing
    @test unavailable.dual === nothing
    @test !unavailable.certified

    complex_variable = QETPSDConstraint.HermitianAffineMatrix(
        :complex_affine, ComplexF32[1 im; -im 2], [1], [ComplexF32[1 0; 0 0]], 1
    )
    complex_constraint = QETPSDConstraint.positive_semidefinite_constraint(complex_variable)
    @test complex_constraint isa QETPSDConstraint.HermitianAffineMatrix{Float32}
    @test complex_constraint.constant == complex_variable.constant
    @test complex_constraint.constant !== complex_variable.constant

    @test_throws MethodError QETPSDConstraint.positive_semidefinite_constraint(
        Matrix{Float64}(I, 2, 2)
    )
end
