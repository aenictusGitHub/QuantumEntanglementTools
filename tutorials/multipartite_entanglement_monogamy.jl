module TutorialMultipartiteEntanglementMonogamy

using LinearAlgebra
using QuantumEntanglementTools

function _one_vs_rest_concurrence(single_system_density::AbstractMatrix)
    determinant = det(single_system_density)
    determinant_real = real(determinant)
    return 2 * sqrt(max(determinant_real, 0.0))
end

"""
    run(; io=stdout)

For three-qubit pure states, compare monogamy-type residuals
tau_A = C_{A|BC}^2 - C_{AB}^2 - C_{AC}^2.

* GHZ: all pairwise two-qubit concurrence vanishes while one-vs-rest is maximal.
* W: pairwise concurrences are nonzero and saturate the monogamy identity.
"""
function run(; io::IO=stdout)
    local_dimensions = (2, 2, 2)
    ghz = ghz_state(2, 3; sparse_output=false)
    w = w_state(3; sparse_output=false)
    ghz_rho = ghz * ghz'
    w_rho = w * w'

    ghz_single = partial_trace(ghz_rho, local_dimensions; trace_out=(2, 3))
    w_single = partial_trace(w_rho, local_dimensions; trace_out=(2, 3))
    ghz_pair_ab = partial_trace(ghz_rho, local_dimensions; trace_out=3)
    w_pair_ab = partial_trace(w_rho, local_dimensions; trace_out=3)
    ghz_pair_ac = partial_trace(ghz_rho, local_dimensions; trace_out=2)
    w_pair_ac = partial_trace(w_rho, local_dimensions; trace_out=2)

    ghz_c_ab = concurrence(ghz_pair_ab; atol=1e-12, rtol=0)
    ghz_c_ac = concurrence(ghz_pair_ac; atol=1e-12, rtol=0)
    w_c_ab = concurrence(w_pair_ab; atol=1e-12, rtol=0)
    w_c_ac = concurrence(w_pair_ac; atol=1e-12, rtol=0)

    ghz_c_a_bc = _one_vs_rest_concurrence(ghz_single)
    w_c_a_bc = _one_vs_rest_concurrence(w_single)

    ghz_tau = ghz_c_a_bc^2 - ghz_c_ab^2 - ghz_c_ac^2
    w_tau = w_c_a_bc^2 - w_c_ab^2 - w_c_ac^2

    ghz_entropy = von_neumann_entropy(ghz_single; base=2, atol=1e-12, rtol=0)
    w_entropy = von_neumann_entropy(w_single; base=2, atol=1e-12, rtol=0)

    @assert isapprox(ghz_c_ab, 0; atol=1e-12, rtol=0)
    @assert isapprox(ghz_c_ac, 0; atol=1e-12, rtol=0)
    @assert ghz_tau > 0.999
    @assert ghz_tau <= 1 + 1e-12
    @assert w_tau >= -1e-12
    @assert w_tau <= 1e-12
    @assert w_c_ab > 0.6
    @assert w_c_ac > 0.6
    @assert w_c_ab <= 1
    @assert w_c_ac <= 1
    @assert isapprox(ghz_entropy, 1; atol=1e-12, rtol=0)
    @assert w_entropy > 0.9
    @assert w_entropy < 1

    println(io, "GHZ C_A(BC)=", ghz_c_a_bc, ", C_AB=", ghz_c_ab, ", C_AC=", ghz_c_ac)
    println(io, "GHZ tau=", ghz_tau)
    println(io, "W C_A(BC)=", w_c_a_bc, ", C_AB=", w_c_ab, ", C_AC=", w_c_ac)
    println(io, "W tau=", w_tau)

    return (
        local_dimensions=local_dimensions,
        ghz_pair_concurrence_ab=ghz_c_ab,
        ghz_pair_concurrence_ac=ghz_c_ac,
        w_pair_concurrence_ab=w_c_ab,
        w_pair_concurrence_ac=w_c_ac,
        ghz_one_vs_rest_concurrence=ghz_c_a_bc,
        w_one_vs_rest_concurrence=w_c_a_bc,
        ghz_tau=ghz_tau,
        w_tau=w_tau,
        ghz_entropy=ghz_entropy,
        w_entropy=w_entropy,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
