# Compatibility wrappers for QETLAB IsCopositive.m and CliqueNumber.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

function _compat_copositivity_hierarchy(mode, function_name::AbstractString)
    normalized = if mode isa Symbol
        Symbol(lowercase(String(mode)))
    elseif mode isa AbstractString
        Symbol(lowercase(strip(mode)))
    else
        throw(ArgumentError("$function_name MODE must be \"sos\", \"nosdp\", :sos, or :nosdp"))
    end
    normalized in (:sos, :nosdp) || throw(
        ArgumentError(
            "$function_name MODE must be \"sos\" or \"nosdp\"; got $(repr(mode))"
        ),
    )
    return normalized
end

"""
    IsCopositive(
        rng::AbstractRNG, C, K=0, MODE="sos";
        structured=true, backend=NoOptimizationBackend(), ...
    )

Compatibility spelling and positional order for QETLAB `IsCopositive`.
The mandatory explicit RNG replaces upstream global, wall-clock-dependent
sampling. By default this returns the certificate-aware
[`QuantumEntanglementTools.CopositivityResult`](@ref).

With `structured=false`, a certified positive or negative conclusion maps to
QETLAB's `1` or `0`. Numerical boundaries, missing or failed backends,
resource limits, and hierarchy-only evidence throw a `DomainError` carrying
the structured result; they are never mapped to `false`.
"""
function IsCopositive(
    rng::AbstractRNG,
    matrix::AbstractMatrix,
    level=0,
    mode="sos";
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    result = copositivity_criterion(
        rng,
        matrix;
        level=level,
        hierarchy=_compat_copositivity_hierarchy(mode, "IsCopositive"),
        backend=backend,
        kwargs...,
    )
    structured && return result
    result.verdict === nothing && throw(
        DomainError(
            result,
            "IsCopositive is inconclusive; request structured=true to inspect its bounds and status",
        ),
    )
    return result.verdict ? 1 : 0
end

"""
    CliqueNumber(
        rng::AbstractRNG, A, K=0, MODE="sos";
        structured=true, backend=NoOptimizationBackend(), ...
    )

Compatibility spelling and positional order for QETLAB `CliqueNumber`.
The default preserves the native certified graph interval and all hierarchy
diagnostics in a
[`QuantumEntanglementTools.CliqueNumberResult`](@ref).

With `structured=false`, return `(ub, lb)` in QETLAB output order, but only
from the graph-certified integer fields. A floating hierarchy candidate is
never substituted for `ub`.
"""
function CliqueNumber(
    rng::AbstractRNG,
    adjacency::AbstractMatrix,
    level=0,
    mode="sos";
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    result = clique_number_bounds(
        rng,
        adjacency;
        level=level,
        hierarchy=_compat_copositivity_hierarchy(mode, "CliqueNumber"),
        backend=backend,
        kwargs...,
    )
    structured && return result
    result.bounds_certified || throw(
        DomainError(
            result,
            "CliqueNumber has no certified integer interval; request structured=true",
        ),
    )
    return (ub=result.upper_bound, lb=result.lower_bound)
end
