# Independently written Julia rendition of the Tiles-UPB workflow on QETLAB's
# "Getting started" homepage: https://qetlab.com/
#
# The package APIs exercised here are independent Julia implementations informed
# by QETLAB UPB.m, IsUPB.m, IsPPT.m, IsSeparable.m, Tensor.m, and
# PartialTranspose.m at revision d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB is copyright Nathaniel Johnston and distributed under BSD-2-Clause; see
# licenses/QETLAB-LICENSE.txt and PROVENANCE.toml.

module TutorialQETLABIntroTiles

using LinearAlgebra
using QuantumEntanglementTools

function _attempt_summary(report)
    return Tuple(attempt.method => attempt.status for attempt in report.attempts)
end

"""
    run(; io=stdout)

Construct the normalized state on the orthogonal complement of the two-qutrit
Tiles UPB. Exact rational arithmetic proves that it is PPT, the exact UPB
certificate supplies the range-criterion entanglement proof, and the numerical
criterion pipeline independently detects entanglement by realignment for a
fixed full-rank depolarized neighbor.
"""
function run(; io::IO=stdout)
    # Each column pair is one unnormalized local product vector. Integer input
    # lets `is_upb` and the projector identities use exact arithmetic.
    left_factors = BigInt[
        1 1 0 0 1
        0 -1 0 1 1
        0 0 1 -1 1
    ]
    right_factors = BigInt[
        1 0 0 1 1
        -1 0 1 0 1
        0 1 -1 0 1
    ]

    upb_analysis = is_upb(left_factors, right_factors; normalization=:allow)

    projector = zeros(Rational{BigInt}, 9, 9)
    for k in axes(left_factors, 2)
        product_vector = Rational{BigInt}.(
            tensor_product(left_factors[:, k], right_factors[:, k])
        )
        projector .+=
            (product_vector * adjoint(product_vector)) / dot(product_vector, product_vector)
    end

    complement_projector = Matrix{Rational{BigInt}}(I, 9, 9) - projector
    rho_exact = complement_projector / tr(complement_projector)
    partial_transpose_exact = partial_transpose(rho_exact, (3, 3); systems=(2,))

    exact_projector =
        ishermitian(projector) &&
        projector * projector == projector &&
        tr(projector) == 5 &&
        complement_projector * complement_projector == complement_projector &&
        tr(complement_projector) == 4
    exact_density_operator = tr(rho_exact) == 1 && rho_exact * rho_exact == rho_exact / 4
    exact_ppt = partial_transpose_exact == rho_exact
    exact_range_entanglement =
        upb_analysis.status === :upb &&
        upb_analysis.reason === :unextendible &&
        upb_analysis.certificate_kind === :exact
    bound_entangled =
        exact_projector && exact_density_operator && exact_ppt && exact_range_entanglement

    # Connect the exact construction to the package's normalized Tiles catalog.
    catalog = upb(:tiles)
    catalog_projector = catalog.global_vectors * adjoint(catalog.global_vectors)
    catalog_projector_error = maximum(abs, catalog_projector - Float64.(projector))

    # The rank-four PPT state is on the zero-eigenvalue boundary, so the floating
    # PPT criterion deliberately returns `CriterionUnknown`. The exact equality
    # above is the PPT proof. A fixed depolarized neighbor moves the independent
    # numerical pipeline away from that eigensolver boundary while preserving
    # PPT exactly and retaining a robust realignment violation.
    rho_boundary = Float64.(rho_exact)
    numeric_ppt = ppt_criterion(rho_boundary, (3, 3); systems=(2,), atol=1e-12, rtol=0)
    depolarizing_weight = BigInt(1) // BigInt(1024)
    maximally_mixed_exact = Matrix{Rational{BigInt}}(I, 9, 9) / 9
    rho_neighbor_exact =
        (1 - depolarizing_weight) * rho_exact + depolarizing_weight * maximally_mixed_exact
    rho_neighbor = Float64.(rho_neighbor_exact)
    neighbor_minimum_eigenvalue = eigmin(Hermitian(rho_neighbor))
    neighbor_ppt = ppt_criterion(rho_neighbor, (3, 3); systems=(2,), atol=1e-12, rtol=0)
    entanglement_report = is_separable(
        rho_neighbor, (3, 3); strategies=(:ppt, :realignment), atol=1e-12, rtol=0
    )
    attempts = _attempt_summary(entanglement_report)
    @assert entanglement_report.status === :entangled
    @assert entanglement_report.certified
    @assert entanglement_report.method === :realignment
    @assert entanglement_report.evidence isa CriterionResult
    realignment_evidence = entanglement_report.evidence
    realignment_margin = realignment_evidence.value - realignment_evidence.threshold

    @assert upb_analysis.status === :upb
    @assert upb_analysis.reason === :unextendible
    @assert upb_analysis.certificate_kind === :exact
    @assert upb_analysis.partitions_examined == 20
    @assert exact_projector
    @assert exact_density_operator
    @assert exact_ppt
    @assert exact_range_entanglement
    @assert bound_entangled
    @assert tr(rho_neighbor_exact) == 1
    @assert partial_transpose(rho_neighbor_exact, (3, 3); systems=(2,)) ==
        rho_neighbor_exact
    @assert tr(rho_neighbor) == 1
    @assert neighbor_minimum_eigenvalue > 1e-5
    @assert neighbor_ppt.status === CriterionSatisfied
    @assert catalog.family === :tiles
    @assert catalog.dimensions == (3, 3)
    @assert catalog.cardinality == 5
    @assert catalog_projector_error <= 1e-14
    @assert numeric_ppt.status === CriterionUnknown
    @assert entanglement_report.status === :entangled
    @assert entanglement_report.certified
    @assert entanglement_report.method === :realignment
    @assert entanglement_report.certificate_kind === :realignment_cross_norm_violation
    @assert attempts == (:ppt => :unknown, :realignment => :entangled)
    @assert realignment_evidence.criterion === :realignment
    @assert realignment_evidence.status === CriterionEntanglementDetected
    @assert realignment_margin > realignment_evidence.tolerance

    println(io, "Tiles product vectors: ", catalog.cardinality)
    println(
        io,
        "Exact structure: rank(complement)=4, PPT=",
        exact_ppt,
        ", no product vector in its range=",
        exact_range_entanglement,
    )
    println(
        io,
        "Floating PPT test: ",
        numeric_ppt.status,
        " (expected at the zero-eigenvalue boundary)",
    )
    println(
        io,
        "Full-rank depolarized neighbor: PPT=",
        neighbor_ppt.status,
        ", separability analysis=",
        entanglement_report.status,
        " via ",
        entanglement_report.certificate_kind,
        "; realignment margin=",
        realignment_margin,
    )

    return (
        catalog_family=catalog.family,
        local_dimensions=catalog.dimensions,
        product_vector_count=catalog.cardinality,
        catalog_projector_error=catalog_projector_error,
        upb_status=upb_analysis.status,
        upb_reason=upb_analysis.reason,
        upb_certificate=upb_analysis.certificate_kind,
        partitions_examined=upb_analysis.partitions_examined,
        complement_rank=Int(tr(complement_projector)),
        exact_trace=tr(rho_exact),
        exact_projector=exact_projector,
        exact_density_operator=exact_density_operator,
        exact_ppt=exact_ppt,
        exact_range_entanglement=exact_range_entanglement,
        bound_entangled=bound_entangled,
        numeric_ppt_status=numeric_ppt.status,
        entanglement_input=:full_rank_depolarized_neighbor,
        depolarizing_weight=depolarizing_weight,
        neighbor_minimum_eigenvalue=neighbor_minimum_eigenvalue,
        neighbor_ppt_status=neighbor_ppt.status,
        entanglement_status=entanglement_report.status,
        entanglement_certified=entanglement_report.certified,
        entanglement_method=entanglement_report.method,
        entanglement_certificate=entanglement_report.certificate_kind,
        attempts=attempts,
        realignment_value=realignment_evidence.value,
        realignment_threshold=realignment_evidence.threshold,
        realignment_tolerance=realignment_evidence.tolerance,
        realignment_margin=realignment_margin,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
