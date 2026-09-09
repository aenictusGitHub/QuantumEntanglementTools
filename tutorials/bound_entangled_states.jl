module TutorialBoundEntangledStates

using LinearAlgebra
using QuantumEntanglementTools

function _attempt_summary(report)
    return Tuple(attempt.method => attempt.status for attempt in report.attempts)
end

"""
    run(; io=stdout)

Build a Horodecki PPT-entangled state from the original bound-entangled family
and verify the full workflow: valid density operator, PPT check, and entanglement
detection by realignment.
"""
function run(; io::IO=stdout)
    local_dimensions = (3, 3)
    family_parameter = 0.3
    bound_state = horodecki_state(family_parameter; dims=local_dimensions)
    minimum_eigenvalue = eigmin(Hermitian(bound_state))
    partial_transpose_state = partial_transpose(bound_state, local_dimensions; systems=(2,))
    minimum_partial_transpose_eigenvalue = eigmin(Hermitian(partial_transpose_state))
    ppt_report = ppt_criterion(
        bound_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    entanglement_report = analyze_entanglement(
        bound_state, local_dimensions; atol=1e-12, rtol=1e-10
    )
    attempts = _attempt_summary(entanglement_report)
    realignment_evidence = entanglement_report.evidence

    @assert entanglement_report.status === :entangled
    @assert entanglement_report.certified
    @assert entanglement_report.method === :realignment
    @assert entanglement_report.certificate_kind === :realignment_cross_norm_violation
    @assert entanglement_report.evidence isa CriterionResult
    @assert ishermitian(bound_state)
    @assert isapprox(tr(bound_state), 1; atol=1e-14, rtol=0)
    @assert minimum_eigenvalue >= -1e-14
    @assert minimum_partial_transpose_eigenvalue >= -1e-14
    @assert ppt_report.status === CriterionUnknown || ppt_report.status === CriterionSatisfied
    @assert attempts[1][1] == :ppt
    @assert attempts[1][2] == :satisfied || attempts[1][2] == :unknown
    @assert attempts[2] == (:realignment => :entangled)
    @assert realignment_evidence.status === CriterionEntanglementDetected
    @assert realignment_evidence.value > realignment_evidence.threshold
    @assert realignment_evidence.value > realignment_evidence.tolerance

    println(io, "Horodecki family parameter: ", family_parameter)
    println(io, "tr(rho)= ", tr(bound_state))
    println(io, "min eig(rho)= ", minimum_eigenvalue)
    println(io, "min eig(T_B rho)= ", minimum_partial_transpose_eigenvalue)
    println(io, "PPT criterion status= ", ppt_report.status)
    println(io, "Entanglement attempts= ", attempts)

    return (
        family_parameter=family_parameter,
        local_dimensions=local_dimensions,
        minimum_eigenvalue=minimum_eigenvalue,
        minimum_partial_transpose_eigenvalue=minimum_partial_transpose_eigenvalue,
        ppt_status=ppt_report.status,
        entanglement_status=entanglement_report.status,
        entanglement_certified=entanglement_report.certified,
        entanglement_method=entanglement_report.method,
        entanglement_certificate=entanglement_report.certificate_kind,
        attempts=attempts,
        realignment_value=realignment_evidence.value,
        realignment_threshold=realignment_evidence.threshold,
        realignment_tolerance=realignment_evidence.tolerance,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
