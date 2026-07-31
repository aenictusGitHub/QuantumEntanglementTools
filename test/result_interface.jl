using LinearAlgebra
using Test
using QuantumEntanglementTools

const QETResultInterface = QuantumEntanglementTools

function _result_interface_optimization(
    status::QETResultInterface.OptimizationStatus;
    certified::Bool=false,
    certificate_kind::Union{Nothing,Symbol}=nothing,
    warnings=(),
)
    return QETResultInterface.OptimizationResult{Float64,Nothing,Nothing,Nothing}(
        status,
        :test_termination,
        :test_primal,
        :test_dual,
        status === QETResultInterface.OptimizationOptimal ? 1.0 : nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        0,
        0.0,
        certified,
        certificate_kind,
        nothing,
        nothing,
        nothing,
        Tuple(String(warning) for warning in warnings),
        "test optimization explanation",
    )
end

@testset "Common human-facing result interface" begin
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    bell_density = bell * bell'
    maximally_mixed_qubits = Matrix{Float64}(I, 4, 4) / 4

    @testset "entanglement reports" begin
        entangled = QETResultInterface.analyze_entanglement(bell, (2, 2))
        unknown = QETResultInterface.analyze_entanglement(
            Matrix{Float64}(I, 9, 9) / 9, (3, 3); atol=1e-12, rtol=0
        )

        @test QETResultInterface.conclusion(entangled) === :entangled
        @test QETResultInterface.is_conclusive(entangled)
        @test QETResultInterface.is_certified(entangled)
        @test QETResultInterface.explain(entangled) == entangled.message
        @test QETResultInterface.conclusion(unknown) === :unknown
        @test !QETResultInterface.is_conclusive(unknown)
        @test !QETResultInterface.is_certified(unknown)

        compact = sprint(show, entangled)
        rich = sprint(show, MIME"text/plain"(), entangled)
        @test !occursin('\n', compact)
        @test occursin("EntanglementReport(status=", compact)
        @test occursin("Conclusion: entangled", rich)
        @test occursin("Certified: yes", rich)
        @test occursin("Certificate: pure state schmidt rank", rich)
        @test occursin("Method: pure schmidt", rich)
        @test occursin("Reason:", rich)
    end

    @testset "necessary criteria" begin
        detected = QETResultInterface.ppt_criterion(bell_density, (2, 2))
        satisfied = QETResultInterface.ppt_criterion(maximally_mixed_qubits, (2, 2))
        boundary = QETResultInterface.ppt_criterion(Diagonal([1.0, 0.0, 0.0, 0.0]), (2, 2))

        @test QETResultInterface.conclusion(detected) === :entangled
        @test QETResultInterface.is_conclusive(detected)
        @test QETResultInterface.is_certified(detected)
        @test QETResultInterface.conclusion(satisfied) === :unknown
        @test !QETResultInterface.is_conclusive(satisfied)
        @test !QETResultInterface.is_certified(satisfied)
        @test QETResultInterface.conclusion(boundary) === :unknown
        @test !QETResultInterface.is_certified(boundary)
        @test QETResultInterface.explain(satisfied) == satisfied.message

        compact = sprint(show, detected)
        rich = sprint(show, MIME"text/plain"(), detected)
        @test !occursin('\n', compact)
        @test occursin("CriterionResult(", compact)
        @test occursin("Criterion: ppt", rich)
        @test occursin("Conclusion: entangled", rich)
        @test occursin("Criterion status: CriterionEntanglementDetected", rich)
        @test occursin("Explanation:", rich)
    end

    @testset "sufficient separable ball" begin
        certified = QETResultInterface.in_separable_ball(maximally_mixed_qubits, (2, 2))
        outside = QETResultInterface.in_separable_ball([1.0, 0.0, 0.0, 0.0], (2, 2))
        boundary = QETResultInterface.in_separable_ball([1 / 3, 1 / 3, 1 / 3, 0.0], (2, 2))

        @test QETResultInterface.conclusion(certified) === :separable
        @test QETResultInterface.is_conclusive(certified)
        @test QETResultInterface.is_certified(certified)
        @test QETResultInterface.conclusion(outside) === :unknown
        @test !QETResultInterface.is_conclusive(outside)
        @test !QETResultInterface.is_certified(outside)
        @test QETResultInterface.conclusion(boundary) === :unknown
        @test QETResultInterface.explain(outside) == outside.message

        compact = sprint(show, certified)
        rich = sprint(show, MIME"text/plain"(), certified)
        @test !occursin('\n', compact)
        @test occursin("SeparableBallResult(status=", compact)
        @test occursin("Conclusion: separable", rich)
        @test occursin("Certificate: sufficient separable ball", rich)
        @test occursin("Ball status: separable certified", rich)
        @test occursin("Explanation:", rich)
    end

    @testset "optimization results" begin
        optimal = _result_interface_optimization(
            QETResultInterface.OptimizationOptimal;
            certified=true,
            certificate_kind=:analytic_test,
            warnings=("check the retained residual",),
        )
        feasible = _result_interface_optimization(QETResultInterface.OptimizationFeasible)
        infeasible = _result_interface_optimization(
            QETResultInterface.OptimizationInfeasible
        )
        unavailable = _result_interface_optimization(
            QETResultInterface.OptimizationBackendUnavailable
        )

        @test QETResultInterface.conclusion(optimal) === :optimal
        @test QETResultInterface.is_conclusive(optimal)
        @test QETResultInterface.is_certified(optimal)
        @test QETResultInterface.conclusion(feasible) === :feasible
        @test QETResultInterface.is_conclusive(feasible)
        @test !QETResultInterface.is_certified(feasible)
        @test QETResultInterface.conclusion(infeasible) === :infeasible
        @test QETResultInterface.is_conclusive(infeasible)
        @test QETResultInterface.conclusion(unavailable) === :unknown
        @test !QETResultInterface.is_conclusive(unavailable)
        @test !QETResultInterface.is_certified(unavailable)
        @test occursin(
            "Warnings: check the retained residual", QETResultInterface.explain(optimal)
        )

        compact = sprint(show, optimal)
        rich = sprint(show, MIME"text/plain"(), optimal)
        @test !occursin('\n', compact)
        @test occursin("OptimizationResult(status=", compact)
        @test occursin("Conclusion: optimal", rich)
        @test occursin("Certificate: analytic test", rich)
        @test occursin("Objective value: 1.0", rich)
        @test occursin("Warnings:", rich)
        @test occursin("- check the retained residual", rich)
    end
end
