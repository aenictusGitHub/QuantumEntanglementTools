module TutorialEntangledStateFamilies

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Compare canonical bipartite families (Werner and isotropic) in `2×2` across
separable and maximally entangled regime points.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2)

    isotropic_separable = isotropic_state(2, 0; sparse_output=false)
    isotropic_entangled = isotropic_state(2, 1; sparse_output=false)
    werner_separable = werner_state(2, 0; sparse_output=false)
    werner_entangled = werner_state(2, 1; sparse_output=false)

    iso_sep_marginal = partial_trace(isotropic_separable, local_dimensions; trace_out=1)
    iso_ent_marginal = partial_trace(isotropic_entangled, local_dimensions; trace_out=2)
    werner_sep_marginal = partial_trace(werner_separable, local_dimensions; trace_out=1)
    werner_ent_marginal = partial_trace(werner_entangled, local_dimensions; trace_out=2)

    iso_sep_ppt = ppt_criterion(isotropic_separable, local_dimensions; systems=(2,), atol=1e-12, rtol=0)
    iso_ent_ppt = ppt_criterion(isotropic_entangled, local_dimensions; systems=(2,), atol=1e-12, rtol=0)
    werner_sep_ppt = ppt_criterion(werner_separable, local_dimensions; systems=(2,), atol=1e-12, rtol=0)
    werner_ent_ppt = ppt_criterion(werner_entangled, local_dimensions; systems=(2,), atol=1e-12, rtol=0)

    iso_sep_report = analyze_entanglement(
        isotropic_separable, local_dimensions; atol=1e-12, rtol=0
    )
    iso_ent_report = analyze_entanglement(
        isotropic_entangled, local_dimensions; atol=1e-12, rtol=0
    )
    werner_sep_report = analyze_entanglement(
        werner_separable, local_dimensions; atol=1e-12, rtol=0
    )
    werner_ent_report = analyze_entanglement(
        werner_entangled, local_dimensions; atol=1e-12, rtol=0
    )

    iso_ent_negativity = negativity(isotropic_entangled, local_dimensions; systems=(2,), atol=0, rtol=0)
    werner_ent_negativity = negativity(werner_entangled, local_dimensions; systems=(2,), atol=0, rtol=0)

    @assert isapprox(tr(isotropic_separable), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(werner_entangled), 1; atol=1e-12, rtol=0)
    @assert isapprox(iso_sep_marginal, Matrix{Float64}(I, 2, 2) / 2; atol=1e-12, rtol=0)
    @assert isapprox(iso_ent_marginal, Matrix{Float64}(I, 2, 2) / 2; atol=1e-12, rtol=0)
    @assert isapprox(werner_sep_marginal, Matrix{Float64}(I, 2, 2) / 2; atol=1e-12, rtol=0)
    @assert isapprox(werner_ent_marginal, Matrix{Float64}(I, 2, 2) / 2; atol=1e-12, rtol=0)
    @assert iso_sep_ppt.status === CriterionSatisfied
    @assert iso_ent_ppt.status === CriterionEntanglementDetected
    @assert werner_sep_ppt.status === CriterionSatisfied
    @assert werner_ent_ppt.status === CriterionEntanglementDetected
    @assert iso_sep_report.status === :separable
    @assert iso_ent_report.status === :entangled
    @assert werner_sep_report.status === :separable
    @assert werner_ent_report.status === :entangled
    @assert iso_sep_report.certified && iso_ent_report.certified
    @assert werner_sep_report.certified && werner_ent_report.certified
    @assert iso_ent_negativity > 0.49
    @assert werner_ent_negativity > 0.49

    println(io, "Isotropic α=0,1 separability: ", iso_sep_report.status, " / ", iso_ent_report.status)
    println(
        io,
        "Werner α=0,1 separability: ",
        werner_sep_report.status,
        " / ",
        werner_ent_report.status,
    )
    println(io, "Isotropic entangled negativity: ", iso_ent_negativity)
    println(io, "Werner entangled negativity: ", werner_ent_negativity)

    return (
        local_dimensions=local_dimensions,
        iso_sep_ppt=iso_sep_ppt.status,
        iso_ent_ppt=iso_ent_ppt.status,
        werner_sep_ppt=werner_sep_ppt.status,
        werner_ent_ppt=werner_ent_ppt.status,
        iso_sep_status=iso_sep_report.status,
        iso_ent_status=iso_ent_report.status,
        werner_sep_status=werner_sep_report.status,
        werner_ent_status=werner_ent_report.status,
        iso_ent_negativity=iso_ent_negativity,
        werner_ent_negativity=werner_ent_negativity,
        iso_sep_marginal=iso_sep_marginal,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
