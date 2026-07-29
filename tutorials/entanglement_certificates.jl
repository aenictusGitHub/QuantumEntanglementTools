module TutorialEntanglementCertificates

using LinearAlgebra
using QuantumEntanglementTools

function _attempt_summary(report)
    return Tuple(attempt.method => attempt.status for attempt in report.attempts)
end

"""
    run(; io=stdout)

Compare certified entangled, certified separable, and inconclusive outcomes.
The higher-dimensional examples retain their ordered attempt histories.
"""
function run(; io::IO=stdout)
    bell_report = analyze_entanglement(bell_state(), (2, 2); atol=0, rtol=1e-12)
    product_report = analyze_entanglement([1.0, 0.0, 0.0, 0.0], (2, 2))

    horodecki = horodecki_state(0.3; dims=(3, 3))
    horodecki_report = analyze_entanglement(horodecki, (3, 3); atol=1e-12, rtol=1e-10)

    maximally_mixed = Matrix{Float64}(I, 9, 9) / 9
    mixed_report = analyze_entanglement(maximally_mixed, (3, 3); atol=1e-12, rtol=0)

    @assert bell_report.status === :entangled
    @assert bell_report.certificate_kind === :pure_state_schmidt_rank
    @assert product_report.status === :separable
    @assert product_report.certificate_kind === :pure_product_decomposition
    @assert horodecki_report.status === :entangled
    @assert _attempt_summary(horodecki_report) ==
        (:ppt => :unknown, :realignment => :entangled)
    @assert mixed_report.status === :unknown
    @assert !mixed_report.certified
    @assert _attempt_summary(mixed_report) ==
        (:ppt => :unknown, :realignment => :unknown, :reduction => :unknown)

    println(
        io, "Pure Bell state: ", bell_report.status, " (", bell_report.certificate_kind, ")"
    )
    println(
        io,
        "Pure product state: ",
        product_report.status,
        " (",
        product_report.certificate_kind,
        ")",
    )
    println(io, "Horodecki attempt history: ", _attempt_summary(horodecki_report))
    println(io, "Maximally mixed 3×3 attempt history: ", _attempt_summary(mixed_report))

    return (
        bell_status=bell_report.status,
        bell_certificate=bell_report.certificate_kind,
        product_status=product_report.status,
        product_certificate=product_report.certificate_kind,
        horodecki_status=horodecki_report.status,
        horodecki_certificate=horodecki_report.certificate_kind,
        horodecki_attempts=_attempt_summary(horodecki_report),
        mixed_status=mixed_report.status,
        mixed_certified=mixed_report.certified,
        mixed_attempts=_attempt_summary(mixed_report),
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
