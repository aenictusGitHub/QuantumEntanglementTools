module TutorialMultipartiteEntanglementMotifs

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Contrast GHZ and W three-qubit states by reducing them to one- and
two-qubit margins:

* GHZ has maximally mixed single-qubit marginals and separable two-qubit cuts.
* W has mixed one-qubit marginals and entangled pairwise marginals.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2, 2)
    ghz = ghz_state(2, 3; sparse_output=false)
    w = w_state(3; sparse_output=false)
    ghz_rho = ghz * ghz'
    w_rho = w * w'

    ghz_single = partial_trace(ghz_rho, local_dimensions; trace_out=(2, 3))
    w_single = partial_trace(w_rho, local_dimensions; trace_out=(2, 3))
    ghz_single_purity = purity(ghz_single; atol=1e-12, rtol=0)
    w_single_purity = purity(w_single; atol=1e-12, rtol=0)
    ghz_single_entropy = von_neumann_entropy(ghz_single; base=2, atol=1e-12, rtol=1e-12)
    w_single_entropy = von_neumann_entropy(w_single; base=2, atol=1e-12, rtol=1e-12)

    ghz_pair = partial_trace(ghz_rho, local_dimensions; trace_out=3)
    w_pair = partial_trace(w_rho, local_dimensions; trace_out=3)
    ghz_pair_report = analyze_entanglement(ghz_pair, (2, 2); atol=1e-12, rtol=0)
    w_pair_report = analyze_entanglement(w_pair, (2, 2); atol=1e-12, rtol=0)
    ghz_pair_negativity = negativity(ghz_pair, (2, 2); systems=(1,), atol=1e-12, rtol=0)
    w_pair_negativity = negativity(w_pair, (2, 2); systems=(1,), atol=1e-12, rtol=0)

    expected_w_entropy = -(2 / 3) * log(2, 2 / 3) - (1 / 3) * log(2, 1 / 3)

    @assert isapprox(tr(ghz_single), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(w_single), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(ghz_pair), 1; atol=1e-12, rtol=0)
    @assert isapprox(tr(w_pair), 1; atol=1e-12, rtol=0)
    @assert isapprox(
        ghz_single, Matrix{Float64}(I, 2, 2) / 2; atol=1e-14, rtol=0
    )
    @assert isapprox(w_single, Diagonal([2 / 3, 1 / 3]); atol=1e-14, rtol=0)
    @assert isapprox(ghz_single_purity, 1 / 2; atol=1e-12, rtol=0)
    @assert isapprox(w_single_purity, 5 / 9; atol=1e-12, rtol=0)
    @assert isapprox(ghz_single_entropy, 1; atol=1e-12, rtol=0)
    @assert isapprox(w_single_entropy, expected_w_entropy; atol=1e-12, rtol=0)
    @assert ghz_pair_report.status === :separable || ghz_pair_report.status === :unknown
    @assert w_pair_report.status === :entangled
    @assert ghz_pair_report.status === :separable || ghz_pair_report.status === :unknown
    @assert ghz_pair_negativity <= 1e-12
    @assert w_pair_negativity > 0
    @assert w_pair_report.certified

    println(io, "GHZ pair marginals: entangled? ", ghz_pair_report.status)
    println(io, "W pair marginals status: ", w_pair_report.status)
    println(io, "W pair negativity: ", w_pair_negativity)
    println(io, "GHZ one-qubit entropy: ", ghz_single_entropy)
    println(io, "W one-qubit entropy: ", w_single_entropy)

    return (
        local_dimensions=local_dimensions,
        ghz_single_purity=ghz_single_purity,
        w_single_purity=w_single_purity,
        ghz_single_entropy=ghz_single_entropy,
        w_single_entropy=w_single_entropy,
        ghz_pair_status=ghz_pair_report.status,
        ghz_pair_certified=ghz_pair_report.certified,
        w_pair_status=w_pair_report.status,
        w_pair_certified=w_pair_report.certified,
        w_pair_negativity=w_pair_negativity,
        ghz_pair_negativity=ghz_pair_negativity,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
