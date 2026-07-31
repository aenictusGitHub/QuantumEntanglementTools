#!/usr/bin/env julia

using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using SparseArrays
using Test

const QET = QuantumEntanglementTools
const MOI = JuMP.MOI
const LOAD_ORDER_PROBE = joinpath(@__DIR__, "load_order_probe.jl")

function extension_module()
    extension = Base.get_extension(QET, :QuantumEntanglementToolsJuMPExt)
    isnothing(extension) && error("JuMP extension is not loaded")
    return extension
end

function run_load_order_probe(order)
    project_file = Base.active_project()
    isnothing(project_file) && error("no active extension-test project")
    command = `$(Base.julia_cmd()) --compiled-modules=no --startup-file=no --history-file=no --project=$(dirname(project_file)) $LOAD_ORDER_PROBE $order`
    return success(pipeline(command; stdout=devnull, stderr=stderr))
end

@testset "JuMP extension discovery" begin
    @test extension_module() isa Module
    @test run_load_order_probe("core_first")
    @test run_load_order_probe("backend_first")
    @test isempty(Test.detect_ambiguities(QET, JuMP, extension_module(); recursive=false))

    status = QET.backend_status()
    @test status.optimization.installed === true
    @test status.optimization.installation_status === :confirmed_loaded
    @test status.optimization.extension_loaded
    @test status.optimization.ready_for_configuration
    @test !status.optimization.ready
    configured = QET.backend_status(
        QET.JuMPBackend(
            Hypatia.Optimizer;
            optimizer_name="Hypatia",
            optimizer_version=Base.pkgversion(Hypatia),
            allow_densify=true,
        ),
    )
    @test configured.configured
    @test configured.ready
    @test configured.optimizer_name == "Hypatia"
    @test configured.optimizer_version == Base.pkgversion(Hypatia)
    @test configured.allow_densify
    @test occursin("ready", configured.message)
end

function hypatia_backend(; options=NamedTuple(), kwargs...)
    return QET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_options=options,
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        kwargs...,
    )
end

function scs_backend(; options=(eps_abs=1.0e-7, eps_rel=1.0e-7), kwargs...)
    return QET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=options,
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=5.0e-6,
        rtol=5.0e-6,
        kwargs...,
    )
end

function complex_psd_problem()
    constant = ComplexF64[0 1 + im; 1 - im 2]
    coefficient = ComplexF64[1 0; 0 0]
    block = QET.HermitianAffineMatrix(:complex_block, constant, [1], [coefficient], 1)
    objective = QET.AffineScalar(0.0, [1.0])
    return QET.SemidefiniteProgram(
        :complex_psd_epigraph,
        :minimize,
        1,
        objective;
        psd_constraints=[block],
        primal_views=[block],
        initial_point=[2.0],
        known_feasible_point=[2.0],
        metadata=(formulation=:complex_schur_complement,),
    )
end

function scalar_problem(
    name,
    sense,
    objective_coefficient;
    lower=nothing,
    upper=nothing,
    initial_point=nothing,
    known_feasible_point=nothing,
)
    objective = QET.AffineScalar(0.0, [objective_coefficient])
    intervals = QET.AffineInterval{Float64}[]
    if lower !== nothing || upper !== nothing
        push!(
            intervals,
            QET.AffineInterval(QET.AffineScalar(0.0, [1.0]), lower, upper, :scalar_bounds),
        )
    end
    return QET.SemidefiniteProgram(
        name, sense, 1, objective; intervals, initial_point, known_feasible_point
    )
end

function infeasible_problem()
    coordinate = QET.AffineScalar(0.0, [1.0])
    return QET.SemidefiniteProgram(
        :infeasible_interval,
        :feasibility,
        1,
        QET.AffineScalar(0.0, [0.0]);
        intervals=[
            QET.AffineInterval(coordinate, 1.0, nothing, :lower),
            QET.AffineInterval(coordinate, nothing, 0.0, :upper),
        ],
    )
end

function mock_factory(status; primal=nothing, add_con_allowed=true)
    return function ()
        mock = MOI.Utilities.MockOptimizer(Float64; add_con_allowed)
        optimize_function = function (optimizer)
            if primal === nothing
                MOI.Utilities.mock_optimize!(optimizer, status)
            else
                MOI.Utilities.mock_optimize!(
                    optimizer, status, (MOI.FEASIBLE_POINT, Float64[primal...])
                )
            end
        end
        MOI.Utilities.set_mock_optimize!(mock, optimize_function)
        return mock
    end
end

function mock_backend(factory; allow_densify=false)
    return QET.JuMPBackend(
        factory;
        optimizer_name="MOI MockOptimizer",
        optimizer_version=Base.pkgversion(MOI),
        silent=false,
        allow_densify,
    )
end

@testset "solver-neutral optimization model" begin
    @testset "backend and limit validation" begin
        @test QET.NoOptimizationBackend() isa QET.AbstractOptimizationBackend
        backend = QET.JuMPBackend(
            Hypatia.Optimizer;
            optimizer_options=(iter_limit=10,),
            optimizer_name="Hypatia",
            optimizer_version=Base.pkgversion(Hypatia),
            time_limit_seconds=2,
            allow_densify=true,
            atol=1.0f-6,
            rtol=2.0f-6,
        )
        @test backend.optimizer_options.iter_limit == 10
        @test backend.time_limit_seconds === 2.0
        @test backend.atol === 1.0f-6
        @test_throws ArgumentError QET.JuMPBackend(Hypatia.Optimizer; atol=-1)
        @test_throws ArgumentError QET.JuMPBackend(
            Hypatia.Optimizer; time_limit_seconds=Inf
        )
        @test_throws ArgumentError QET.JuMPBackend(
            Hypatia.Optimizer; optimizer_options=(:x => 1, :x => 2)
        )
        @test_throws ArgumentError QET.OptimizationLimits(max_variables=0)

        problem = scalar_problem(:budget, :feasibility, 0.0)
        @test_throws ArgumentError QET.SemidefiniteProgram(
            :over_budget,
            :feasibility,
            1,
            problem.objective;
            limits=QET.OptimizationLimits(max_variables=1, max_model_entries=1),
            primal_views=[
                QET.HermitianAffineMatrix(
                    :view, zeros(2, 2), [1], [Matrix{Float64}(I, 2, 2)], 1
                ),
            ],
        )
    end

    @testset "Hermitian coordinates and maps" begin
        X = QET.hermitian_variable(:X, 2; coefficient_type=Float32)
        @test X.variable_count == 4
        @test eltype(X.constant) === ComplexF32
        coordinates = Float32[2, 3, 4, 5]
        evaluated = Matrix(QET.evaluate_affine(X, coordinates))
        @test evaluated == ComplexF32[2 4 + 5im; 4 - 5im 3]
        @test ishermitian(evaluated)
        @test QET.evaluate_affine(QET.trace_affine(X), coordinates) == 5

        reduced = QET.partial_trace_affine(
            QET.tensor_affine(X; right=sparse(Matrix{Float32}(I, 2, 2)), name=:X_tensor_I),
            (2, 2);
            trace_out=(2,),
        )
        @test Matrix(QET.evaluate_affine(reduced, coordinates)) ≈ 2evaluated

        transposed = QET.partial_transpose_affine(
            QET.tensor_affine(X; right=Float32[1 0; 0 0]), (2, 2); systems=(1,)
        )
        transposed_value = Matrix(QET.evaluate_affine(transposed, coordinates))
        @test transposed_value ≈ kron(transpose(evaluated), Float32[1 0; 0 0])

        choi = QET.hermitian_variable(:J, 4; coefficient_type=Float64)
        trace_preserving = QET.choi_trace_preserving_affine(choi, 2, 2)
        identity_channel_choi = zeros(ComplexF64, 4, 4)
        maximally_entangled = ComplexF64[1, 0, 0, 1]
        identity_channel_choi .= maximally_entangled * maximally_entangled'
        coordinates_choi = zeros(16)
        coordinates_choi[1:4] .= real.(diag(identity_channel_choi))
        position = 5
        for column in 2:4, row in 1:(column - 1)
            coordinates_choi[position] = real(identity_channel_choi[row, column])
            position += 1
        end
        for column in 2:4, row in 1:(column - 1)
            coordinates_choi[position] = imag(identity_channel_choi[row, column])
            position += 1
        end
        @test Matrix(QET.evaluate_affine(trace_preserving, coordinates_choi)) ≈
            zeros(ComplexF64, 2, 2)
        equalities = QET.hermitian_equalities(trace_preserving)
        @test length(equalities) == 4
        @test maximum(
            abs(QET.evaluate_affine(constraint.function_data, coordinates_choi)) for
            constraint in equalities
        ) ≤ 10eps()

        complex_matrix = ComplexF64[2 1 + 2im; 1 - 2im 3]
        block = QET.real_block_embedding(complex_matrix)
        @test issymmetric(block)
        @test sort(eigvals(Symmetric(block))) ≈
            sort(repeat(eigvals(Hermitian(complex_matrix)), inner=2))
        dual = Matrix{Float64}(I, 4, 4)
        mapped_dual = QET.hermitian_dual_from_real_block(dual)
        @test mapped_dual == 2Matrix{ComplexF64}(I, 2, 2)
        @test real(dot(dual, block)) ≈ real(dot(mapped_dual, complex_matrix))
    end

    @testset "model validation and residuals" begin
        @test_throws ArgumentError QET.AffineScalar(NaN, [1.0])
        @test_throws DimensionMismatch QET.AffineScalar(0.0, [1.0]; variable_count=2)
        @test_throws ArgumentError QET.AffineInterval(
            QET.AffineScalar(0.0, [1.0]), 1.0, 1.0
        )
        @test_throws ArgumentError QET.HermitianAffineMatrix(
            :bad, ComplexF64[1 1; 0 1], Int[], Matrix{ComplexF64}[], 1
        )
        @test_throws DimensionMismatch QET.SemidefiniteProgram(
            :bad, :minimize, 2, QET.AffineScalar(0.0, [1.0])
        )

        problem = complex_psd_problem()
        @test QET.primal_residual(problem, [1.0]; allow_densify=true) ≤ 20eps()
        @test QET.primal_residual(problem, [0.5]; allow_densify=true) > 0.1
        @test_throws ArgumentError QET.primal_residual(problem, [1.0])
        zero_dual = [zeros(4, 4)]
        @test QET.dual_stationarity_residual(problem, Float64[], Float64[], zero_dual) == 1
        @test QET._optimization_consistent_status(
            QET.OptimizationOptimal, 1.0, 0.0, 1.0e-6
        ) === QET.OptimizationInconsistent
    end

    @testset "missing backend and explicit densification" begin
        problem = complex_psd_problem()
        missing = QET.solve_optimization(problem)
        @test missing.status === QET.OptimizationBackendUnavailable
        @test missing.termination_status === :backend_unavailable
        @test missing.objective_value === nothing
        @test !missing.certified

        denied = QET.solve_optimization(
            problem,
            QET.JuMPBackend(
                Hypatia.Optimizer;
                optimizer_name="Hypatia",
                optimizer_version=Base.pkgversion(Hypatia),
            ),
        )
        @test denied.status === QET.OptimizationUnsupported
        @test denied.termination_status === :densification_not_allowed
    end
end

@testset "JuMP optimization extension" begin
    @testset "optimal complex SDP and second-solver cross-check" begin
        problem = complex_psd_problem()
        result = QET.solve_optimization(problem, hypatia_backend())
        @test result.status === QET.OptimizationOptimal
        @test result.termination_status === :optimal
        @test result.primal_status === :feasible_point
        @test result.objective_value ≈ 1 atol = 5.0e-7
        @test result.objective_bound ≈ 1 atol = 5.0e-7
        @test result.dual_objective_value ≈ 1 atol = 5.0e-7
        @test result.absolute_gap ≤ 5.0e-7
        @test result.primal_residual ≤ 5.0e-7
        @test result.dual_residual ≤ 5.0e-7
        @test result.primal.views.complex_block ≈ ComplexF64[1 1 + im; 1 - im 2] atol =
            5.0e-7
        @test length(result.dual.psd_real_blocks) == 1
        @test ishermitian(first(result.dual.psd_hermitian_blocks))
        @test result.optimizer.modeling_layer === :JuMP
        @test result.optimizer.reported_optimizer_name == "Hypatia"
        @test result.optimizer.optimizer_version == Base.pkgversion(Hypatia)
        @test result.optimizer.complex_psd_embedding === :real_block
        @test result.iterations isa Int
        @test result.solve_time_seconds ≥ 0
        @test !result.certified

        second = QET.solve_optimization(problem, scs_backend())
        @test second.status in (QET.OptimizationOptimal, QET.OptimizationFeasible)
        @test second.objective_value ≈ result.objective_value atol = 2.0e-5
        @test second.primal_residual ≤ 2.0e-5
        @test second.optimizer.reported_optimizer_name == "SCS"
    end

    @testset "infeasible and unbounded certificates" begin
        infeasible = QET.solve_optimization(infeasible_problem(), hypatia_backend())
        @test infeasible.status === QET.OptimizationInfeasible
        @test infeasible.termination_status === :infeasible
        @test infeasible.dual_status === :infeasibility_certificate
        @test !infeasible.certified

        unbounded_problem = scalar_problem(
            :unbounded,
            :maximize,
            1.0;
            lower=0.0,
            initial_point=[0.0],
            known_feasible_point=[0.0],
        )
        unbounded = QET.solve_optimization(unbounded_problem, hypatia_backend())
        @test unbounded.termination_status === :dual_infeasible
        @test unbounded.status === QET.OptimizationUnbounded
        @test unbounded.primal_status === :infeasibility_certificate
    end

    @testset "limits, failures, unsupported and malformed backends" begin
        limited = QET.solve_optimization(
            complex_psd_problem(), hypatia_backend(options=(iter_limit=0,))
        )
        @test limited.status === QET.OptimizationLimit
        @test limited.termination_status === :iteration_limit
        @test !limited.certified

        numerical = QET.solve_optimization(
            scalar_problem(:numerical, :feasibility, 0.0),
            mock_backend(mock_factory(MOI.NUMERICAL_ERROR)),
        )
        @test numerical.status === QET.OptimizationNumericalFailure
        @test numerical.termination_status === :numerical_error
        @test numerical.objective_value === nothing

        unsupported = QET.solve_optimization(
            scalar_problem(:unsupported, :minimize, 1.0; lower=0.0),
            mock_backend(mock_factory(MOI.OPTIMAL; primal=[0.0], add_con_allowed=false)),
        )
        @test unsupported.status === QET.OptimizationUnsupported
        @test unsupported.termination_status === :unsupported_model

        malformed = QET.solve_optimization(
            scalar_problem(:malformed, :feasibility, 0.0), mock_backend(() -> 42)
        )
        @test malformed.status === QET.OptimizationMalformedBackend
        @test malformed.termination_status === :backend_exception
        @test occursin("failed", malformed.message)
    end

    @testset "feasible/inaccurate and primal inconsistency" begin
        feasible_problem = scalar_problem(:feasible, :feasibility, 0.0)
        feasible = QET.solve_optimization(
            feasible_problem, mock_backend(mock_factory(MOI.ALMOST_OPTIMAL; primal=[0.0]))
        )
        @test feasible.status === QET.OptimizationFeasible
        @test feasible.primal_residual == 0
        @test !feasible.certified

        inconsistent_problem = scalar_problem(:inconsistent, :feasibility, 0.0; lower=1.0)
        inconsistent = QET.solve_optimization(
            inconsistent_problem, mock_backend(mock_factory(MOI.OPTIMAL; primal=[0.0]))
        )
        @test inconsistent.status === QET.OptimizationInconsistent
        @test inconsistent.termination_status === :optimal
        @test inconsistent.primal_residual == 1
        @test !isempty(inconsistent.warnings)
        @test !inconsistent.certified
    end
end

include("absolute_ppt.jl")
include("symmetric_extensions.jl")
include("psd_constraint.jl")
include("matsumoto_fidelity_model.jl")
include("top_k_p_norm_epigraph.jl")
include("sk_norms.jl")
include("state_discrimination.jl")
include("separability.jl")
include("nonlocal_games.jl")
include("channel_optimization.jl")
include("polynomial_sos.jl")
include("copositivity_clique.jl")
include("coherence_optimization.jl")
