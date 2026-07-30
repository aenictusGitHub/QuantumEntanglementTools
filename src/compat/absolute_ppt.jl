# Source-informed Julia compatibility wrappers based on QETLAB
# AbsPPTConstraints.m and IsAbsPPT.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

function _compat_abs_ppt_limit(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("LIM must be a nonnegative integer"))
    0 <= value <= typemax(Int) ||
        throw(ArgumentError("LIM must be a nonnegative integer representable as Int"))
    return Int(value)
end

function _compat_abs_ppt_constraint_dims(dims)
    dims isa Integer || return dims
    local_dimension = _positive_dimension(dims, "DIM")
    return (local_dimension, local_dimension)
end

function _compat_is_abs_ppt_dims(input, dims)
    dims === nothing || return dims
    total = input isa AbstractVector ? length(input) : max(size(input)...)
    total > 0 || throw(ArgumentError("RHO must be nonempty"))
    return round(Int, sqrt(total))
end

function _compat_abs_ppt_matrices(
    family, escape_on_negative::Bool, requested_limit::Int, return_structured::Bool
)
    return_structured && return family
    expected_partial =
        (escape_on_negative && family.status === AbsPPTEnumerationEarlyViolation) ||
        (requested_limit > 0 && family.status === AbsPPTEnumerationConstraintLimit)
    if !family.exhaustive && !expected_partial
        throw(
            DomainError(
                family,
                "AbsPPTConstraints reached $(family.status) before completing the " *
                "requested family; pass structured=true to inspect its limits or " *
                "raise the explicit resource budget",
            ),
        )
    end
    return [constraint.matrix for constraint in family.constraints]
end

"""
    AbsPPTConstraints(
        LAM, DIM=nothing, ESC_IF_NPOS=0, LIM=0;
        structured=false, max_constraints=2612,
        max_work=500_000_000, max_entries=1_000_000,
        sparse_output=false, allow_densify=false,
        atol=nothing, rtol=nothing
    )

Compatibility spelling and positional argument order for QETLAB's finite
absolute-PPT LMI builder. Numeric input returns a vector of package-owned
matrices. `ESC_IF_NPOS=1` stops at the first robust negative matrix, and
positive `LIM` returns at most that many matrices. Set `structured=true` to
receive the complete
[`QuantumEntanglementTools.AbsPPTConstraintFamily`](@ref).

For safety, `LIM=0` retains the native finite `max_constraints=2612` default;
if that guard, `max_work`, or `max_entries` prevents completion, the plain
matrix form raises a `DomainError` instead of silently presenting a partial
family as complete. Pass `max_constraints=nothing` only after explicitly
accepting unbounded enumeration.

This wrapper preserves pinned `AbsPPTConstraints` scalar-DIM behavior:
`DIM=d` means `(d,d)`. The native [`abs_ppt_constraints`](@ref) instead uses
the consistent rule `(d,length(LAM)÷d)`.
"""
function AbsPPTConstraints(
    input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dims=nothing,
    escape_if_nonpositive=0,
    limit=0;
    structured=false,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
    sparse_output::Bool=false,
    allow_densify::Bool=false,
    atol=nothing,
    rtol=nothing,
)
    escape_on_negative = _flag(escape_if_nonpositive, "ESC_IF_NPOS")
    requested_limit = _compat_abs_ppt_limit(limit)
    return_structured = _flag(structured, "structured")
    dimensions = _compat_abs_ppt_constraint_dims(dims)
    constraint_limit = requested_limit > 0 ? requested_limit : max_constraints
    family = abs_ppt_constraints(
        input;
        dims=dimensions,
        stop_on_violation=escape_on_negative,
        max_constraints=constraint_limit,
        max_work,
        max_entries,
        sparse_output,
        allow_densify,
        atol,
        rtol,
    )
    return _compat_abs_ppt_matrices(
        family, escape_on_negative, requested_limit, return_structured
    )
end

"""
    AbsPPTConstraints(
        LAM::AbstractVector{<:AffineScalar},
        DIM=nothing, ESC_IF_NPOS=0, LIM=0;
        assume_ordered=true, structured=false, limits...
    )

Build package-owned affine Hermitian matrices for a symbolic eigenvalue
vector. Unlike a CVX expression, each returned value remains solver
independent. `ESC_IF_NPOS` must be zero because symbolic PSD cannot be decided
without a model. The surrounding optimization model must impose
`LAM[1] >= ... >= LAM[end] >= 0`.
"""
function AbsPPTConstraints(
    eigenvalues::AbstractVector{<:AffineScalar},
    dims=nothing,
    escape_if_nonpositive=0,
    limit=0;
    assume_ordered::Bool=true,
    structured=false,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
)
    escape_on_negative = _flag(escape_if_nonpositive, "ESC_IF_NPOS")
    escape_on_negative &&
        throw(ArgumentError("ESC_IF_NPOS=1 is unavailable for affine eigenvalues"))
    requested_limit = _compat_abs_ppt_limit(limit)
    return_structured = _flag(structured, "structured")
    dimensions = _compat_abs_ppt_constraint_dims(dims)
    constraint_limit = requested_limit > 0 ? requested_limit : max_constraints
    family = abs_ppt_constraints(
        eigenvalues;
        dims=dimensions,
        assume_ordered,
        max_constraints=constraint_limit,
        max_work,
        max_entries,
    )
    return _compat_abs_ppt_matrices(family, false, requested_limit, return_structured)
end

"""
    IsAbsPPT(RHO, DIM=nothing; structured=false, kwargs...)

Compatibility spelling for [`is_abs_ppt`](@ref). The plain form preserves
QETLAB's tri-state values: `1` is certified absolutely PPT, `0` is certified
not absolutely PPT, and `-1` is inconclusive. Capped work, floating numerical
boundaries, and backend failures map to `-1`, never to a mathematical
negative. Set `structured=true` to retain the native result and all evidence.

When `DIM` is omitted, this wrapper preserves QETLAB `IsAbsPPT` dimension
guessing: it rounds the square root of the total dimension and infers the
second factor. This notably treats a length-six spectrum as `(2,3)`. A scalar
`DIM=d` likewise means `(d,length(RHO)÷d)`.
"""
function IsAbsPPT(
    input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dims=nothing;
    structured=false,
    kwargs...,
)
    return_structured = _flag(structured, "structured")
    dimensions = _compat_is_abs_ppt_dims(input, dims)
    result = is_abs_ppt(input; dims=dimensions, kwargs...)
    return_structured && return result
    result.verdict === true && return 1
    result.verdict === false && return 0
    return -1
end
