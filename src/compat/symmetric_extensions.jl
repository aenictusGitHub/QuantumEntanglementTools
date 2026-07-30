# Source-informed Julia compatibility wrappers based on QETLAB
# SymmetricExtension.m, SymmetricInnerExtension.m, and RandomPPTState.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

function _compat_symmetric_dimensions(operator, dims)
    dims === nothing || return dims
    total_dimension = max(size(operator)...)
    total_dimension > 0 || throw(ArgumentError("X must be nonempty"))
    return round(Int, sqrt(total_dimension))
end

function _compat_symmetric_tolerance(value)
    tolerance = value === nothing ? eps(Float64)^(1 / 4) : value
    tolerance isa Real && !(tolerance isa Bool) ||
        throw(ArgumentError("TOL must be a finite nonnegative real number"))
    isfinite(tolerance) && tolerance >= zero(tolerance) ||
        throw(ArgumentError("TOL must be a finite nonnegative real number"))
    return tolerance
end

function _compat_symmetric_output(result, include_witness::Bool)
    result.verdict === nothing && throw(
        DomainError(
            result,
            "the symmetric-extension computation ended with $(result.status); " *
            "pass structured=true to retain its inconclusive diagnostics or " *
            "supply an explicit optimization backend",
        ),
    )
    value = result.verdict ? 1 : 0
    include_witness || return value
    certificate = if result.verdict
        result.extension
    elseif result.witness === nothing
        nothing
    else
        result.witness.operator
    end
    certificate === nothing && throw(
        DomainError(
            result,
            "the computation was conclusive but did not reconstruct the positional WIT output",
        ),
    )
    return (ex=value, wit=certificate)
end

"""
    SymmetricExtension(
        X, K=2, DIM=nothing, PPT=0, BOS=0, TOL=nothing;
        return_witness=false, structured=false,
        backend=NoOptimizationBackend(), ...
    )

Compatibility spelling and positional argument order for QETLAB's outer
symmetric-extension hierarchy. The one-output form returns integer `1` or `0`
only after a validated native certificate. Set `return_witness=true` to
receive the named tuple `(ex, wit)`, which can be destructured positionally.
That form disables solver-free shortcuts that do not reconstruct `WIT`, just
as requesting QETLAB's second output selects its CVX path.

An unavailable backend, resource limit, numerical boundary, or invalid
certificate raises `DomainError` in the positional forms instead of becoming
`0`. Set `structured=true` to receive the complete
[`QuantumEntanglementTools.SymmetricExtensionResult`](@ref), including
`verdict === nothing`. `PPT` and `BOS` accept only `Bool`, `0`, or `1`.
`TOL` maps to native absolute tolerance with zero relative tolerance and
defaults to QETLAB's `eps(Float64)^(1/4)`.
"""
function SymmetricExtension(
    operator::AbstractMatrix{<:Number},
    order=2,
    dims=nothing,
    ppt=0,
    bosonic=0,
    tolerance=nothing;
    return_witness=false,
    structured=false,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    include_witness = _flag(return_witness, "return_witness")
    return_structured = _flag(structured, "structured")
    include_witness &&
        return_structured &&
        throw(ArgumentError("return_witness and structured cannot both be true"))
    result = symmetric_extension(
        operator;
        order=_positive_dimension(order, "K"),
        dims=_compat_symmetric_dimensions(operator, dims),
        ppt=_flag(ppt, "PPT"),
        bosonic=_flag(bosonic, "BOS"),
        backend,
        prefer_analytic=(!include_witness),
        atol=_compat_symmetric_tolerance(tolerance),
        rtol=0,
        allow_densify,
        max_dense_entries,
        max_order,
        limits,
    )
    return_structured && return result
    return _compat_symmetric_output(result, include_witness)
end

"""
    SymmetricInnerExtension(
        X, K=2, DIM=nothing, PPT=0, TOL=nothing;
        return_witness=false, structured=false,
        backend=NoOptimizationBackend(), ...
    )

Compatibility spelling and positional argument order for QETLAB's inner
hierarchy. Conclusive positional output follows [`SymmetricExtension`](@ref).
General calls require an explicit optimization backend.

When `return_witness=true` and `ex == 0`, the returned `wit` is only a
separator from the selected inner approximation. It is **not automatically an
entanglement witness**. Use `structured=true` to retain the explicit
`entanglement_witness=false` marker and warning.
"""
function SymmetricInnerExtension(
    operator::AbstractMatrix{<:Number},
    order=2,
    dims=nothing,
    ppt=0,
    tolerance=nothing;
    return_witness=false,
    structured=false,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    include_witness = _flag(return_witness, "return_witness")
    return_structured = _flag(structured, "structured")
    include_witness &&
        return_structured &&
        throw(ArgumentError("return_witness and structured cannot both be true"))
    result = symmetric_inner_extension(
        operator;
        order=_positive_dimension(order, "K"),
        dims=_compat_symmetric_dimensions(operator, dims),
        ppt=_flag(ppt, "PPT"),
        backend,
        atol=_compat_symmetric_tolerance(tolerance),
        rtol=0,
        allow_densify,
        max_dense_entries,
        max_order,
        limits,
    )
    return_structured && return result
    return _compat_symmetric_output(result, include_witness)
end

function _compat_random_ppt_iterations(value)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("MAX_ITS must be a nonnegative integer or Inf"))
    value == Inf && return 1
    isfinite(value) && isinteger(value) && value >= zero(value) ||
        throw(ArgumentError("MAX_ITS must be a nonnegative integer or Inf"))
    value <= typemax(Int) || throw(ArgumentError("MAX_ITS is too large for Int"))
    return Int(value)
end

"""
    RandomPPTState(
        rng, DIM, RNK=nothing, TOL=1e-12, MAX_ITS=1;
        structured=false, ...
    )

Compatibility spelling for QETLAB's random PPT constructor, with a mandatory
leading `rng::AbstractRNG`. A verified construction returns the density
matrix. Set `structured=true` to retain construction, rank, tolerance, work,
and verification diagnostics. Any failed or limited plain call raises
`DomainError`; an unverified candidate is never returned as a state.

The Julia-native low-rank route is a single bounded random separable mixture,
so `MAX_ITS=Inf` is accepted only as a compatibility sentinel and maps to that
one bounded step. The output distribution and random stream do not match
MATLAB or Octave.
"""
function RandomPPTState(
    rng::AbstractRNG,
    dims,
    ranks=nothing,
    tolerance=1.0e-12,
    max_iterations=1;
    structured=false,
    construction::Symbol=:auto,
    real::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    max_dense_entries=1_000_000,
    max_work=100_000_000,
)
    return_structured = _flag(structured, "structured")
    result = random_ppt_state(
        rng,
        dims;
        ranks,
        construction,
        real,
        T,
        atol=_compat_symmetric_tolerance(tolerance),
        rtol=0,
        max_iterations=_compat_random_ppt_iterations(max_iterations),
        max_dense_entries,
        max_work,
    )
    return_structured && return result
    result.verified && result.state !== nothing || throw(
        DomainError(
            result,
            "RandomPPTState ended with $(result.status); pass structured=true " *
            "to inspect the bounded construction diagnostics",
        ),
    )
    return copy(result.state)
end
