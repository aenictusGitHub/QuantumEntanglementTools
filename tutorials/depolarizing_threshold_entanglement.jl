module TutorialDepolarizingThresholdEntanglement

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Use a one-sided depolarizing channel on a Bell state to show the entanglement
threshold in 2×2 systems (`p = 1/3`):

* `p = 0` is maximally mixed (separable),
* `p = 2/3` remains entangled (via negative partial transpose).
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2)
    bell = bell_state()
    bell_density = bell * bell'

    separable_channel = depolarizing_channel(2, 0)
    entangled_channel = depolarizing_channel(2, 2 / 3)

    separable_complete = is_completely_positive(
        separable_channel; allow_densify=true
    )
    entangled_complete = is_completely_positive(
        entangled_channel; allow_densify=true
    )

    separable_state = partial_map(bell_density, separable_channel, 2, local_dimensions)
    entangled_state = partial_map(bell_density, entangled_channel, 2, local_dimensions)

    separable_overlap = real(tr(separable_state * bell_density))
    entangled_overlap = real(tr(entangled_state * bell_density))

    separable_ppt = ppt_criterion(
        separable_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    entangled_ppt = ppt_criterion(
        entangled_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    separable_report = analyze_entanglement(
        separable_state, local_dimensions; atol=1e-12, rtol=0
    )
    entangled_report = analyze_entanglement(
        entangled_state, local_dimensions; atol=1e-12, rtol=0
    )

    separable_concurrence = concurrence(separable_state; atol=1e-12, rtol=0)
    entangled_concurrence = concurrence(entangled_state; atol=1e-12, rtol=0)

    @assert separable_complete.status === MatrixPredicateSatisfied
    @assert entangled_complete.status === MatrixPredicateSatisfied
    @assert is_trace_preserving(separable_channel)
    @assert is_trace_preserving(entangled_channel)
    @assert is_unital(separable_channel)
    @assert is_unital(entangled_channel)
    @assert isapprox(tr(separable_state), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(entangled_state), 1; atol=1e-12, rtol=0)
    @assert separable_ppt.status === CriterionSatisfied
    @assert entangled_ppt.status === CriterionEntanglementDetected
    @assert separable_report.status === :separable
    @assert entangled_report.status === :entangled
    @assert separable_concurrence <= 1e-12
    @assert entangled_concurrence > 0.3
    @assert separable_overlap <= 0.2500000000001
    @assert entangled_overlap >= 0.6666

    println(io, "Channel CP/TP/unital (separable): ", separable_complete.status)
    println(io, "Channel CP/TP/unital (entangled): ", entangled_complete.status)
    println(io, "Separable case status: ", separable_report.status)
    println(io, "Entangled case status: ", entangled_report.status)
    println(io, "Overlap with |Φ+>: ", separable_overlap, " / ", entangled_overlap)
    println(io, "Concurrence: ", separable_concurrence, " / ", entangled_concurrence)

    return (
        local_dimensions=local_dimensions,
        is_cptp_0=is_trace_preserving(separable_channel) &&
            separable_complete.status === MatrixPredicateSatisfied &&
            is_unital(separable_channel),
        is_cptp_2_3=is_trace_preserving(entangled_channel) &&
            entangled_complete.status === MatrixPredicateSatisfied && is_unital(entangled_channel),
        separable_ppt=separable_ppt.status,
        entangled_ppt=entangled_ppt.status,
        separable_status=separable_report.status,
        entangled_status=entangled_report.status,
        separable_concurrence=separable_concurrence,
        entangled_concurrence=entangled_concurrence,
        separable_overlap=separable_overlap,
        entangled_overlap=entangled_overlap,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
