module TutorialEntanglementMonotoneNoiseSweep

using LinearAlgebra
using QuantumEntanglementTools

function _sweep_reports(dephasing_channel_values::AbstractVector{<:Real})
    local_dimensions = (2, 2)
    bell = bell_state()
    bell_density = bell * bell'

    concurrence_values = Vector{Float64}(undef, length(dephasing_channel_values))
    negativity_values = Vector{Float64}(undef, length(dephasing_channel_values))
    log_neg_values = Vector{Float64}(undef, length(dephasing_channel_values))
    statuses = Vector{Symbol}(undef, length(dephasing_channel_values))

    for index in eachindex(dephasing_channel_values)
        p = dephasing_channel_values[index]
        channel = dephasing_channel(2, p)
        state = partial_map(bell_density, channel, 2, local_dimensions)
        concurrence_values[index] = concurrence(state; atol=1e-12, rtol=0)
        negativity_values[index] = negativity(state, local_dimensions; systems=(2,), atol=1e-12, rtol=0)
        log_neg_values[index] = logarithmic_negativity(
            state, local_dimensions; systems=(2,), atol=1e-12, rtol=0, base=2
        )
        statuses[index] = analyze_entanglement(state, local_dimensions; atol=1e-12, rtol=0).status
    end

    return (
        concurrence_values=concurrence_values,
        negativity_values=negativity_values,
        log_neg_values=log_neg_values,
        statuses=statuses,
    )
end

"""
    run(; io=stdout)

Track two monotones (concurrence and negativity) along a dephasing sweep on a
Bell state and confirm monotonic entanglement decay.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2)
    channel_points = [1, 1 / 2, 1 / 4, 0]
    sweep = _sweep_reports(channel_points)

    @assert isapprox(sweep.concurrence_values[1], 1; atol=1e-12, rtol=0)
    @assert sweep.concurrence_values[1] > sweep.concurrence_values[2]
    @assert sweep.concurrence_values[2] > sweep.concurrence_values[3]
    @assert sweep.concurrence_values[3] > sweep.concurrence_values[4]
    @assert sweep.concurrence_values[4] <= 1e-12

    @assert sweep.negativity_values[1] > sweep.negativity_values[2]
    @assert sweep.negativity_values[2] > sweep.negativity_values[3]
    @assert sweep.negativity_values[3] > sweep.negativity_values[4]
    @assert sweep.negativity_values[4] <= 1e-12

    @assert sweep.log_neg_values[1] > sweep.log_neg_values[2]
    @assert sweep.log_neg_values[2] > sweep.log_neg_values[3]
    @assert sweep.log_neg_values[3] > sweep.log_neg_values[4]
    @assert sweep.log_neg_values[4] <= 1e-12

    @assert sweep.statuses[1:3] == [:entangled, :entangled, :entangled]
    @assert sweep.statuses[4] === :separable || sweep.statuses[4] === :unknown

    println(io, "p values: ", channel_points)
    println(io, "Concurrence: ", sweep.concurrence_values)
    println(io, "Negativity: ", sweep.negativity_values)
    println(io, "Log-negativity: ", sweep.log_neg_values)

    return (
        local_dimensions=local_dimensions,
        channel_points=channel_points,
        concurrence_values=sweep.concurrence_values,
        negativity_values=sweep.negativity_values,
        log_neg_values=sweep.log_neg_values,
        statuses=sweep.statuses,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
