#!/usr/bin/env julia

using LinearAlgebra

length(ARGS) == 1 || error("expected one load-order name")
order = only(ARGS)

if order == "core_first"
    @eval using QuantumEntanglementTools
    Base.get_extension(
        QuantumEntanglementTools, :QuantumEntanglementToolsEntanglementDetectionExt
    ) === nothing || error("extension loaded before its weak dependency")
    length(QuantumEntanglementTools.available_entanglement_backends()) == 1 ||
        error("optional backend appeared before its weak dependency loaded")
    absent_error = try
        QuantumEntanglementTools.detect_entanglement(
            Matrix{Float64}(I, 4, 4) / 4,
            (2, 2),
            QuantumEntanglementTools.EntanglementDetectionSearch(; max_iteration=1),
        )
        nothing
    catch error
        error
    end
    absent_error isa ArgumentError ||
        error("missing optional dependency did not raise ArgumentError")
    occursin("is not loaded", sprint(showerror, absent_error)) ||
        error("missing optional dependency error was not actionable")
    @eval using EntanglementDetection
elseif order == "backend_first"
    @eval using EntanglementDetection
    @eval using QuantumEntanglementTools
else
    error("unknown load order: $order")
end

# Requiring both packages again must be idempotent.
@eval using EntanglementDetection
@eval using QuantumEntanglementTools

extension = Base.get_extension(
    QuantumEntanglementTools, :QuantumEntanglementToolsEntanglementDetectionExt
)
extension === nothing && error("extension did not load for $order")
backends = QuantumEntanglementTools.available_entanglement_backends()
length(backends) == 2 || error("optional backend discovery failed for $order")
last(backends) isa QuantumEntanglementTools.EntanglementDetectionBackend ||
    error("unexpected optional backend descriptor for $order")
capabilities = QuantumEntanglementTools.backend_capabilities(last(backends))
capabilities.version == v"0.2.2" ||
    error("unexpected backend version for $order: $(capabilities.version)")
capabilities.isolation === :child_process ||
    error("optional backend does not advertise child-process isolation")
