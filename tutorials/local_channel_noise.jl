module TutorialLocalChannelNoise

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Apply a completely depolarizing channel to one half of a Bell state and
inspect both the channel representation and certificate-aware conclusions.
"""
function run(; io::IO=stdout)
    bell = bell_state()
    bell_density = bell * bell'
    channel = depolarizing_channel(2)

    @assert is_completely_positive(channel)
    @assert is_trace_preserving(channel)
    @assert is_unital(channel)

    choi = choi_representation(channel)
    superoperator = superoperator_representation(channel)
    representation_error = maximum(abs, choi_matrix(choi) - choi_matrix(superoperator))
    @assert representation_error <= 1e-12

    noisy_state = partial_map(bell_density, channel, 2, (2, 2))
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    @assert isapprox(noisy_state, maximally_mixed; atol=1e-12, rtol=0)

    initial_report = analyze_entanglement(bell_density, (2, 2); atol=1e-12, rtol=0)
    noisy_report = analyze_entanglement(noisy_state, (2, 2); atol=1e-12, rtol=0)
    @assert initial_report.status === :entangled
    @assert initial_report.certificate_kind === :negative_partial_transpose_witness
    @assert noisy_report.status === :separable
    @assert noisy_report.certificate_kind === :ppt_low_dimension_theorem

    println(
        io,
        "Channel diagnostics: CP=",
        is_completely_positive(channel),
        ", TP=",
        is_trace_preserving(channel),
        ", unital=",
        is_unital(channel),
    )
    println(io, "Choi/superoperator conversion error: ", representation_error)
    println(
        io,
        "Bell state conclusion: ",
        initial_report.status,
        " (",
        initial_report.certificate_kind,
        ")",
    )
    println(
        io,
        "After local depolarization: ",
        noisy_report.status,
        " (",
        noisy_report.certificate_kind,
        ")",
    )

    return (
        completely_positive=is_completely_positive(channel),
        trace_preserving=is_trace_preserving(channel),
        unital=is_unital(channel),
        representation_error=representation_error,
        initial_status=initial_report.status,
        initial_certificate=initial_report.certificate_kind,
        noisy_status=noisy_report.status,
        noisy_certificate=noisy_report.certificate_kind,
        noisy_purity=purity(noisy_state),
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
