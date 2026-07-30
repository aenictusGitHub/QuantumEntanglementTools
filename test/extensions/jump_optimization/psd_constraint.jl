using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using Test

const PSDConstraintQET = QuantumEntanglementTools

if !isdefined(PSDConstraintQET, :positive_semidefinite_constraint)
    Base.include(
        PSDConstraintQET,
        joinpath(@__DIR__, "..", "..", "..", "src", "optimization", "psd_constraints.jl"),
    )
end

function psd_constraint_backend(; kwargs...)
    return PSDConstraintQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=1.0e-7,
        rtol=1.0e-7,
        kwargs...,
    )
end

function psd_constraint_program()
    variable = PSDConstraintQET.hermitian_variable(
        :rho, 2; variable_count=4, coefficient_type=Float64
    )
    constraint = PSDConstraintQET.positive_semidefinite_constraint(variable)
    trace_constraint = PSDConstraintQET.AffineEquality(
        PSDConstraintQET.AffineScalar(-1.0, [1.0, 1.0, 0.0, 0.0]), :unit_trace
    )
    objective = PSDConstraintQET.AffineScalar(0.0, [1.0, 0.0, 0.0, 0.0])
    return PSDConstraintQET.SemidefiniteProgram(
        :psd_constraint_extension_test,
        :minimize,
        4,
        objective;
        equalities=[trace_constraint],
        psd_constraints=[constraint],
        primal_views=[constraint],
        initial_point=[0.5, 0.5, 0.0, 0.0],
        known_feasible_point=[0.5, 0.5, 0.0, 0.0],
    )
end

@testset "WP3 PSD constraint optional extension" begin
    program = psd_constraint_program()
    explicit_model = JuMP.Model(Hypatia.Optimizer)
    explicit_coordinates = JuMP.@variable(explicit_model, [1:(program.variable_count)])
    explicit_constraint = PSDConstraintQET.positive_semidefinite_constraint(
        explicit_model,
        first(program.psd_constraints),
        explicit_coordinates;
        allow_densify=true,
    )
    @test explicit_constraint isa JuMP.ConstraintRef
    @test_throws ArgumentError PSDConstraintQET.positive_semidefinite_constraint(
        explicit_model, first(program.psd_constraints), explicit_coordinates
    )
    @test_throws DimensionMismatch PSDConstraintQET.positive_semidefinite_constraint(
        explicit_model,
        first(program.psd_constraints),
        explicit_coordinates[1:(end - 1)];
        allow_densify=true,
    )

    result = PSDConstraintQET.solve_optimization(program, psd_constraint_backend())
    @test result.status === PSDConstraintQET.OptimizationOptimal
    @test result.termination_status === :optimal
    @test result.primal_status === :feasible_point
    @test result.dual_status === :feasible_point
    @test result.objective_value ≈ 0 atol = 2e-8
    @test result.objective_bound ≈ 0 atol = 2e-8
    @test result.absolute_gap <= 1e-7
    @test result.primal_residual <= 1e-7
    @test result.dual_residual <= 1e-7
    @test result.primal !== nothing
    @test result.dual !== nothing
    @test result.primal.views.rho[1, 1] ≈ 0 atol = 2e-8
    @test result.primal.views.rho[2, 2] ≈ 1 atol = 2e-8
    @test minimum(eigvals(Hermitian(Matrix(result.primal.views.rho)))) >= -1e-7
    @test result.optimizer.configured_optimizer_name == "Hypatia"
    @test !result.certified

    limited = PSDConstraintQET.solve_optimization(
        program, psd_constraint_backend(optimizer_options=(iter_limit=0,))
    )
    @test limited.status === PSDConstraintQET.OptimizationLimit
    @test limited.termination_status === :iteration_limit
    @test !limited.certified

    malformed_backend = PSDConstraintQET.JuMPBackend(
        () -> error("intentional PSD-constraint factory failure");
        optimizer_name="intentional malformed factory",
        allow_densify=true,
    )
    malformed = PSDConstraintQET.solve_optimization(program, malformed_backend)
    @test malformed.status === PSDConstraintQET.OptimizationMalformedBackend
    @test malformed.primal === nothing
    @test malformed.dual === nothing
    @test occursin("failed", malformed.message)
end
