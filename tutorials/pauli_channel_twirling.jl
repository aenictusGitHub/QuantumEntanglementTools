module TutorialPauliChannelTwirl

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Build three Pauli channels:

* identity probabilities to keep the Bell state entangled,
* fully symmetric Pauli noise giving the maximally mixed output,
* intermediate symmetric-noise point still above the entanglement threshold.

The example highlights Pauli-twirled noise as a practical route to an
isotropic two-qubit family.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2)
    bell = bell_state()
    bell_density = bell * bell'

    identity_channel = pauli_channel([1, 0, 0, 0])
    fully_twirled_channel = pauli_channel(fill(1 / 4, 4))
    twirled_channel = pauli_channel([0.7, 0.1, 0.1, 0.1])

    identity_complete = is_completely_positive(identity_channel; allow_densify=true)
    twirled_complete = is_completely_positive(twirled_channel; allow_densify=true)
    full_complete = is_completely_positive(fully_twirled_channel; allow_densify=true)

    identity_state = partial_map(bell_density, identity_channel, 2, local_dimensions)
    twirled_state = partial_map(bell_density, twirled_channel, 2, local_dimensions)
    fully_twirled_state = partial_map(
        bell_density, fully_twirled_channel, 2, local_dimensions
    )

    fully_mixed = Matrix{Float64}(I, 4, 4) / 4
    fully_mixed_overlap = real(tr(fully_twirled_state * bell_density))
    twirled_overlap = real(tr(twirled_state * bell_density))

    identity_report = analyze_entanglement(
        identity_state, local_dimensions; atol=1e-12, rtol=0
    )
    twirled_report = analyze_entanglement(twirled_state, local_dimensions; atol=1e-12, rtol=0)
    fully_twirled_report = analyze_entanglement(
        fully_twirled_state, local_dimensions; atol=1e-12, rtol=0
    )

    identity_ppt = ppt_criterion(
        identity_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    twirled_ppt = ppt_criterion(
        twirled_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    fully_twirled_ppt = ppt_criterion(
        fully_twirled_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )

    identity_concurrence = concurrence(identity_state; atol=1e-12, rtol=0)
    twirled_concurrence = concurrence(twirled_state; atol=1e-12, rtol=0)
    fully_twirled_concurrence = concurrence(fully_twirled_state; atol=1e-12, rtol=0)

    @assert identity_complete.status === MatrixPredicateSatisfied
    @assert twirled_complete.status === MatrixPredicateSatisfied
    @assert full_complete.status === MatrixPredicateSatisfied
    @assert is_trace_preserving(identity_channel)
    @assert is_unital(identity_channel)
    @assert is_trace_preserving(twirled_channel)
    @assert is_unital(twirled_channel)
    @assert is_trace_preserving(fully_twirled_channel)
    @assert is_unital(fully_twirled_channel)
    @assert isapprox(fully_twirled_state, fully_mixed; atol=1e-12, rtol=0)
    @assert isapprox(fully_mixed_overlap, 1 / 4; atol=1e-12, rtol=0)
    @assert identity_report.status === :entangled
    @assert twirled_report.status === :entangled
    @assert fully_twirled_report.status === :separable
    @assert twirled_ppt.status === CriterionEntanglementDetected
    @assert fully_twirled_ppt.status === CriterionSatisfied
    @assert identity_ppt.status === CriterionEntanglementDetected
    @assert isapprox(identity_concurrence, 1; atol=1e-12, rtol=0)
    @assert twirled_concurrence >= 0.4 - 10 * eps()
    @assert fully_twirled_concurrence <= 1e-12

    println(io, "Pauli identity concurrence: ", identity_concurrence)
    println(io, "Pauli twirled concurrence: ", twirled_concurrence)
    println(io, "Pauli full-twirled concurrence: ", fully_twirled_concurrence)
    println(io, "Twirled Bell overlap: ", twirled_overlap)
    println(io, "Fully twirled overlap: ", fully_mixed_overlap)

    return (
        local_dimensions=local_dimensions,
        identity_status=identity_report.status,
        twirled_status=twirled_report.status,
        fully_twirled_status=fully_twirled_report.status,
        twirled_ppt=twirled_ppt.status,
        fully_twirled_ppt=fully_twirled_ppt.status,
        identity_is_cptp=identity_complete.status === MatrixPredicateSatisfied &&
            is_trace_preserving(identity_channel) && is_unital(identity_channel),
        twirled_is_cptp=twirled_complete.status === MatrixPredicateSatisfied &&
            is_trace_preserving(twirled_channel) && is_unital(twirled_channel),
        fully_twirled_is_cptp=full_complete.status === MatrixPredicateSatisfied &&
            is_trace_preserving(fully_twirled_channel) &&
            is_unital(fully_twirled_channel),
        identity_concurrence=identity_concurrence,
        twirled_concurrence=twirled_concurrence,
        fully_twirled_concurrence=fully_twirled_concurrence,
        twirled_overlap=twirled_overlap,
        fully_mixed_overlap=fully_mixed_overlap,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
