# Source-informed independent Julia compatibility wrappers based on the
# executable contracts of QETLAB NPAHierarchy.m, NonlocalGameLB.m,
# XORGameValue.m, BellInequalityMax.m, BellInequalityMaxQubits.m, BCSGameLB.m,
# and BCSGameValue.m at d8589610f00cff106537268dee2e2a1153f3a601.
# Upstream source authors named in those files include Nathaniel Johnston,
# Vincent Russo, and Mateus Araújo. QETLAB: Copyright 2014 Nathaniel Johnston,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# This file is included inside
# `MATLABCompat` after the native bindings are imported by the parent module.
# Structured results are the default; positional scalar output is available
# only through an explicit `structured=false` request and only when the
# corresponding exact value or independently certified bound exists.

function _compat_nonlocal_regime(value, allowed, name::AbstractString)
    value isa Union{Symbol,AbstractString} ||
        throw(ArgumentError("$name must be a Symbol or string"))
    normalized = Symbol(lowercase(replace(String(value), '-' => '_')))
    normalized === :nosignal && (normalized = :no_signalling)
    choices = join(string.(allowed), ", ")
    normalized in allowed || throw(ArgumentError("$name must be one of $choices"))
    return normalized
end

function _compat_nonlocal_scalar(result, name::AbstractString)
    result.exact && result.value !== nothing && return result.value
    result.certified_lower && result.lower_bound !== nothing && return result.lower_bound
    result.certified_upper && result.upper_bound !== nothing && return result.upper_bound
    return throw(
        ErrorException(
            "$name has no exact value or certified scalar bound; inspect the " *
            "structured result and optimization statuses",
        ),
    )
end

"""
    NPAHierarchy(cg, desc, k=1; backend=NoOptimizationBackend(),
                 structured=true, kwargs...)

QETLAB-compatible argument order with status-safe output. Numerical NPA
feasibility is never rounded to `0` or `1`; `structured=false` is accepted only
for the exact level-zero basic-condition branch.
"""
function NPAHierarchy(
    cg,
    desc,
    k=1;
    backend=NoOptimizationBackend(),
    structured::Bool=true,
    atol=nothing,
    rtol=nothing,
    kwargs...,
)
    scenario = BellScenario(desc)
    behavior = CollinsGisinBehavior(cg, scenario; atol=atol, rtol=rtol)
    result = npa_membership(behavior; level=k, backend=backend, kwargs...)
    structured && return result
    result.verdict === nothing && throw(
        ErrorException(
            "NPAHierarchy numerical output is inconclusive as an exact Boolean; " *
            "use structured=true",
        ),
    )
    return Int(result.verdict)
end

"""
    NonlocalGameLB(rng, d, p, V, verbose=1; structured=true, kwargs...)

Mandatory explicit-RNG replacement for QETLAB's global-RNG heuristic.
`verbose` is validated for compatibility but library code never prints.
"""
function NonlocalGameLB(rng, d, p, V, verbose=1; structured::Bool=true, kwargs...)
    verbose in (0, 1, false, true) || throw(ArgumentError("VERBOSE must be zero or one"))
    result = nonlocal_game_lower_bound(rng, d, p, V; kwargs...)
    structured && return result
    (!result.certified_lower || result.lower_bound === nothing) && throw(
        ErrorException(
            "NonlocalGameLB produced only a numerical see-saw candidate, not a " *
            "certified lower bound; " *
            "inspect the structured status",
        ),
    )
    return result.lower_bound
end

"""
    XORGameValue(p, f, vtype="classical"; structured=true, kwargs...)

Compatibility wrapper for exact classical and numerical quantum-bound routes.
"""
function XORGameValue(p, f, vtype="classical"; structured::Bool=true, kwargs...)
    regime = _compat_nonlocal_regime(vtype, (:classical, :quantum), "VTYPE")
    result = xor_game_value(p, f; regime=regime, kwargs...)
    return structured ? result : _compat_nonlocal_scalar(result, "XORGameValue")
end

"""
    BellInequalityMax(coefficients, desc, notation, mtype="classical", k=1;
                      structured=true, kwargs...)
"""
function BellInequalityMax(
    coefficients, desc, notation, mtype="classical", k=1; structured::Bool=true, kwargs...
)
    regime = _compat_nonlocal_regime(mtype, (:classical, :quantum, :no_signalling), "MTYPE")
    notation_symbol = _compat_nonlocal_regime(
        notation,
        (:fp, :fc, :cg, :full_probability, :full_correlator, :collins_gisin),
        "NOTATION",
    )
    result = bell_inequality_bound(
        coefficients,
        BellScenario(desc);
        notation=notation_symbol,
        regime=regime,
        level=k,
        kwargs...,
    )
    return structured ? result : _compat_nonlocal_scalar(result, "BellInequalityMax")
end

"""
    BellInequalityMaxQubits(joint, alice, bob, alice_values, bob_values;
                            structured=true, kwargs...)

The native corrected rectangular-setting formulation is always used. With
`structured=false`, return `(upper_bound, relaxation_state)` only after a
certified upper bound exists.
"""
function BellInequalityMaxQubits(
    joint, alice, bob, alice_values, bob_values; structured::Bool=true, kwargs...
)
    result = bell_inequality_qubit_bound(
        joint, alice, bob, alice_values, bob_values; kwargs...
    )
    structured && return result
    (!result.certified_upper || result.upper_bound === nothing) && throw(
        ErrorException(
            "BellInequalityMaxQubits produced no certified upper bound; " *
            "inspect the structured status",
        ),
    )
    return (result.upper_bound, result.relaxation_state)
end

"""
    BCSGameLB(rng, d, constraints, verbose=1; structured=true, kwargs...)

Mandatory explicit-RNG BCS lower-bound wrapper.
"""
function BCSGameLB(rng, d, constraints, verbose=1; structured::Bool=true, kwargs...)
    verbose in (0, 1, false, true) || throw(ArgumentError("VERBOSE must be zero or one"))
    result = bcs_game_lower_bound(rng, d, constraints; kwargs...)
    structured && return result
    (!result.certified_lower || result.lower_bound === nothing) && throw(
        ErrorException(
            "BCSGameLB produced only a numerical see-saw candidate, not a " *
            "certified lower bound; inspect " *
            "the structured status",
        ),
    )
    return result.lower_bound
end

"""
    BCSGameValue(constraints, mtype="classical", k=1;
                 structured=true, kwargs...)

Compatibility mapping backed by the independently specified
`nonlocal_game_value` dispatcher that replaces the dependency absent from the
pinned QETLAB source tree.
"""
function BCSGameValue(constraints, mtype="classical", k=1; structured::Bool=true, kwargs...)
    regime = _compat_nonlocal_regime(mtype, (:classical, :quantum, :no_signalling), "MTYPE")
    result = bcs_game_value(constraints; regime=regime, level=k, kwargs...)
    return structured ? result : _compat_nonlocal_scalar(result, "BCSGameValue")
end
