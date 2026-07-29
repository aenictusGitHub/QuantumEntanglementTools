module TutorialSeparabilityExamples

using LinearAlgebra
using QuantumEntanglementTools

function _attempt_summary(report)
    return Tuple(attempt.method => attempt.status for attempt in report.attempts)
end

"""
    run(; io=stdout)

Work through separability certificates for a pure product state, an explicit
two-qubit product-state mixture, and a higher-dimensional product-state
mixture. Also demonstrate that failure of a sufficient separable-ball test is
not evidence of entanglement.
"""
function run(; io::IO=stdout)
    ket0 = [1.0, 0.0]
    ket1 = [0.0, 1.0]

    # |0> ⊗ |1> is separable by its explicit product construction.
    product_state = tensor_product(ket0, ket1)
    pure_report = analyze_entanglement(product_state, (2, 2); atol=0, rtol=0)

    # A full-rank convex mixture of four computational-basis product states.
    basis_states_2x2 = (
        tensor_product(ket0, ket0),
        tensor_product(ket0, ket1),
        tensor_product(ket1, ket0),
        tensor_product(ket1, ket1),
    )
    weights_2x2 = [3 / 8, 1 / 8, 1 / 8, 3 / 8]
    mixed_state_2x2 = sum(
        weight * (state * state') for (weight, state) in zip(weights_2x2, basis_states_2x2)
    )
    mixed_report = analyze_entanglement(mixed_state_2x2, (2, 2); atol=0, rtol=0)
    mixed_ball = in_separable_ball(mixed_state_2x2, (2, 2); atol=0, rtol=0)

    # This 3×3 state is also an explicit convex mixture of product-basis
    # projectors. Its binary-exact weights keep normalization auditable.
    basis3 = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
    basis_states_3x3 = Tuple(
        tensor_product(left, right) for left in basis3 for right in basis3
    )
    weights_3x3 = vcat(fill(1 / 8, 7), [1 / 16, 1 / 16])
    mixed_state_3x3 = sum(
        weight * (state * state') for (weight, state) in zip(weights_3x3, basis_states_3x3)
    )
    higher_report = analyze_entanglement(mixed_state_3x3, (3, 3); atol=0, rtol=0)
    higher_ball = in_separable_ball(mixed_state_3x3, (3, 3); atol=0, rtol=0)

    # A pure product density matrix lies outside the separable ball even though
    # its vector representation has already supplied a separability certificate.
    product_ball = in_separable_ball(product_state * product_state', (2, 2); atol=0, rtol=0)

    @assert pure_report.status === :separable
    @assert pure_report.certified
    @assert pure_report.certificate_kind === :pure_product_decomposition
    @assert sum(weights_2x2) == 1
    @assert isapprox(tr(mixed_state_2x2), 1; atol=0, rtol=0)
    @assert mixed_report.status === :separable
    @assert mixed_report.certified
    @assert mixed_report.certificate_kind === :ppt_low_dimension_theorem
    @assert mixed_ball.status === :separable_certified
    @assert sum(weights_3x3) == 1
    @assert isapprox(tr(mixed_state_3x3), 1; atol=0, rtol=0)
    @assert higher_report.status === :unknown
    @assert !higher_report.certified
    @assert _attempt_summary(higher_report) ==
        (:ppt => :unknown, :realignment => :unknown, :reduction => :unknown)
    @assert higher_ball.status === :separable_certified
    @assert product_ball.status === :outside_ball

    println(
        io, "1. Pure |0>⊗|1>: ", pure_report.status, " (", pure_report.certificate_kind, ")"
    )
    println(
        io,
        "2. Explicit 2×2 mixture: ",
        mixed_report.status,
        " (",
        mixed_report.certificate_kind,
        "); ball=",
        mixed_ball.status,
        "; purity=",
        mixed_ball.purity,
        " < boundary=",
        mixed_ball.boundary,
    )
    println(
        io,
        "3. Explicit 3×3 mixture: pipeline=",
        higher_report.status,
        "; attempts=",
        _attempt_summary(higher_report),
        "; ball=",
        higher_ball.status,
        "; purity=",
        higher_ball.purity,
        " < boundary=",
        higher_ball.boundary,
    )
    println(
        io,
        "4. Pure product density against the ball: ",
        product_ball.status,
        " (not an entanglement verdict)",
    )

    return (
        pure_status=pure_report.status,
        pure_certified=pure_report.certified,
        pure_certificate=pure_report.certificate_kind,
        mixed_status=mixed_report.status,
        mixed_certificate=mixed_report.certificate_kind,
        mixed_ball_status=mixed_ball.status,
        higher_pipeline_status=higher_report.status,
        higher_pipeline_certified=higher_report.certified,
        higher_attempts=_attempt_summary(higher_report),
        higher_ball_status=higher_ball.status,
        product_ball_status=product_ball.status,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
