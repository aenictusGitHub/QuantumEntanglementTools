module TutorialDephasingNoiseThreshold

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Show how phase damping preserves entanglement for partial dephasing and destroys it
at full dephasing (`p = 0`), while the channel itself remains CPTP and unital.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2)
    bell = bell_state()
    bell_density = bell * bell'

    fully_dephasing_channel = dephasing_channel(2, 0)
    partial_dephasing_channel = dephasing_channel(2, 1 / 2)

    complete_positivity = is_completely_positive(
        partial_dephasing_channel; allow_densify=true
    )
    trace_preserving = is_trace_preserving(partial_dephasing_channel)
    unital = is_unital(partial_dephasing_channel)

    fully_dephased_state = partial_map(
        bell_density, fully_dephasing_channel, 2, local_dimensions
    )
    partially_dephased_state = partial_map(
        bell_density, partial_dephasing_channel, 2, local_dimensions
    )

    fully_dephased_target = zeros(4, 4)
    fully_dephased_target[1, 1] = 1 / 2
    fully_dephased_target[4, 4] = 1 / 2

    @assert isapprox(tr(fully_dephased_state), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(partially_dephased_state), 1; atol=1e-12, rtol=0)
    @assert isapprox(
        fully_dephased_state, fully_dephased_target; atol=1e-12, rtol=0
    )
    @assert isapprox(
        abs(partially_dephased_state[1, 4]), 1 / 4; atol=1e-12, rtol=0
    )

    fully_dephased_ppt = ppt_criterion(
        fully_dephased_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    partially_dephased_ppt = ppt_criterion(
        partially_dephased_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )
    fully_dephased_report = analyze_entanglement(
        fully_dephased_state, local_dimensions; atol=1e-12, rtol=0
    )
    partially_dephased_report = analyze_entanglement(
        partially_dephased_state, local_dimensions; atol=1e-12, rtol=0
    )
    partially_dephased_negativity = negativity(
        partially_dephased_state, local_dimensions; systems=(2,), atol=1e-12, rtol=0
    )

    @assert complete_positivity.status !== MatrixPredicateViolated
    @assert trace_preserving
    @assert unital
    @assert fully_dephased_ppt.status === CriterionSatisfied ||
        fully_dephased_ppt.status === CriterionUnknown
    @assert partially_dephased_ppt.status === CriterionEntanglementDetected
    @assert fully_dephased_report.status === :separable || fully_dephased_report.status === :unknown
    @assert partially_dephased_report.status === :entangled
    @assert partially_dephased_report.certified
    @assert partially_dephased_negativity > 0
    @assert partially_dephased_report.certificate_kind ===
        :negative_partial_transpose_witness

    println(io, "Channel CP status: ", complete_positivity.status)
    println(io, "Fully dephased Bell pair separability: ", fully_dephased_report.status)
    println(
        io,
        "Partially dephased Bell pair entanglement: ",
        partially_dephased_report.status,
    )
    println(io, "Remaining off-diagonal magnitude: ", abs(partially_dephased_state[1, 4]))
    println(io, "Negativity: ", partially_dephased_negativity)

    return (
        local_dimensions=local_dimensions,
        is_cptp=trace_preserving && complete_positivity.status !== MatrixPredicateViolated,
        is_unital=unital,
        complete_positive=complete_positivity.status !== MatrixPredicateViolated,
        fully_dephased_status=fully_dephased_report.status,
        fully_dephased_ppt=fully_dephased_ppt.status,
        fully_dephased_certified=fully_dephased_report.certified,
        partial_dephased_status=partially_dephased_report.status,
        partial_dephased_ppt=partially_dephased_ppt.status,
        partial_dephased_certified=partially_dephased_report.certified,
        partial_offdiagonal=abs(partially_dephased_state[1, 4]),
        partial_negativity=partially_dephased_negativity,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
