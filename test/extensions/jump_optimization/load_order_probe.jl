#!/usr/bin/env julia

length(ARGS) == 1 || error("expected one load-order name")
order = only(ARGS)

if order == "core_first"
    @eval using QuantumEntanglementTools
    QET = QuantumEntanglementTools
    Base.get_extension(QET, :QuantumEntanglementToolsJuMPExt) === nothing ||
        error("JuMP extension loaded before its weak dependency")

    problem = QET.SemidefiniteProgram(
        :missing_backend_probe,
        :minimize,
        1,
        QET.AffineScalar(0.0, [1.0]);
        intervals=[
            QET.AffineInterval(QET.AffineScalar(0.0, [1.0]), 0.0, nothing, :nonnegative)
        ],
        known_feasible_point=[0.0],
    )
    missing = QET.solve_optimization(
        problem, QET.JuMPBackend(() -> nothing; optimizer_name="not loaded")
    )
    missing.status === QET.OptimizationBackendUnavailable ||
        error("missing JuMP did not return OptimizationBackendUnavailable")
    missing.termination_status === :extension_unavailable ||
        error("missing JuMP did not retain its extension-unavailable status")
    @eval using JuMP
elseif order == "backend_first"
    @eval using JuMP
    @eval using QuantumEntanglementTools
else
    error("unknown load order: $order")
end

@eval using Hypatia
QET = QuantumEntanglementTools
extension = Base.get_extension(QET, :QuantumEntanglementToolsJuMPExt)
extension === nothing && error("JuMP extension did not load for $order")

problem = QET.SemidefiniteProgram(
    :load_order_probe,
    :minimize,
    1,
    QET.AffineScalar(0.0, [1.0]);
    intervals=[
        QET.AffineInterval(QET.AffineScalar(0.0, [1.0]), 0.0, nothing, :nonnegative)
    ],
    known_feasible_point=[0.0],
)
backend = QET.JuMPBackend(
    Hypatia.Optimizer; optimizer_name="Hypatia", optimizer_version=Base.pkgversion(Hypatia)
)
result = QET.solve_optimization(problem, backend)
result.status === QET.OptimizationOptimal ||
    error("load-order probe solve was not optimal: $(result.status)")
abs(result.objective_value) <= 1.0e-7 ||
    error("load-order probe returned the wrong objective")
