# Source-informed Julia compatibility wrappers based on QETLAB
# LocalDistinguishability.m, IsSeparable.m, and UPBSepDistinguishable.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and Alessandro Cosentino,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

using ..QuantumEntanglementTools:
    _is_package_validated_entanglement_report, _localdisc_legacy_output_valid

const _SEPARABILITY_OPTIMIZATION_COMPAT_LOADED = true

const _COMPAT_IS_SEPARABLE_STRENGTH_ZERO = (
    :ppt,
    :low_rank_ppt,
    :realignment,
    :centered_realignment,
    :qubit_qudit,
    :rank4_chow,
    :separable_ball,
    :rank_one_identity,
    :operator_schmidt_rank,
    :positive_maps,
)

function _compat_separability_tolerance(value, default, name::AbstractString)
    tolerance = value === nothing ? default : value
    tolerance isa Real && !(tolerance isa Bool) && isfinite(tolerance) ||
        throw(ArgumentError("$name must be a finite nonnegative real number"))
    tolerance >= zero(tolerance) ||
        throw(ArgumentError("$name must be a finite nonnegative real number"))
    return tolerance
end

function _compat_local_discrimination_dimensions(states, dims)
    dimension = if states isa AbstractMatrix{<:Number}
        size(states, 1)
    elseif states isa Tuple || states isa AbstractVector
        isempty(states) && throw(ArgumentError("X must contain at least one state"))
        first_state = first(states)
        first_state isa AbstractMatrix{<:Number} ||
            throw(ArgumentError("a state collection must contain numeric matrices"))
        size(first_state, 1)
    else
        throw(
            ArgumentError(
                "X must be a numeric pure-state matrix or a tuple/vector of density matrices",
            ),
        )
    end
    dimension > 0 || throw(ArgumentError("X must have positive state dimension"))
    dimensions = dims === nothing ? round(Int, sqrt(dimension)) : dims
    return if dimensions isa Integer
        _expand_scalar_dimension(dimensions, dimension, "DIM")
    else
        dimensions
    end
end

function _compat_separability_dimensions(operator, dims)
    size(operator, 1) > 0 || throw(ArgumentError("X must be nonempty"))
    dimensions = dims === nothing ? round(Int, sqrt(size(operator, 1))) : dims
    return if dimensions isa Integer
        _expand_scalar_dimension(dimensions, size(operator, 1), "DIM")
    else
        dimensions
    end
end

function _compat_local_distinguishability_output(result)
    _localdisc_legacy_output_valid(result) || throw(
        DomainError(
            result,
            "LocalDistinguishability result did not pass package-owned structural validation",
        ),
    )
    result.relaxation_value === nothing && throw(
        DomainError(
            result,
            "LocalDistinguishability did not establish a single relaxation value; " *
            "request structured=true to inspect bounds and solver status",
        ),
    )
    result.measurement === nothing && throw(
        DomainError(
            result,
            "LocalDistinguishability has no residual-checked hierarchy POVM; " *
            "request structured=true to inspect the result",
        ),
    )
    hasproperty(result.residuals, :valid) && result.residuals.valid || throw(
        DomainError(
            result,
            "LocalDistinguishability has no validated hierarchy POVM; " *
            "request structured=true to inspect residuals",
        ),
    )
    dual = result.dual_solution
    dual !== nothing &&
    hasproperty(dual, :completeness_operator) &&
    dual.completeness_operator !== nothing || throw(
        DomainError(
            result,
            "LocalDistinguishability has no reconstructed completeness dual; " *
            "request structured=true to inspect the optimization result",
        ),
    )
    return (
        dist=result.relaxation_value,
        meas=result.measurement,
        dual_sol=dual.completeness_operator,
    )
end

"""
    LocalDistinguishability(
        X, P=nothing, DIM=nothing, COPIES=2, PPT=1, BOS=1,
        TOL=eps(Float64)^(1/4);
        structured=true, backend=NoOptimizationBackend(), kwargs...
    )

Compatibility spelling and positional argument order for QETLAB's local
state-discrimination hierarchy. The default returns the status-rich native
result. `P`, `DIM`, `COPIES`, `PPT`, `BOS`, and `TOL` map to the native prior,
dimension, hierarchy-order, PPT, bosonic, and absolute-tolerance arguments.
States and priors are validated but never normalized; their represented
normalization must be exact.

With `structured=false`, the three QETLAB output positions are returned as
`(dist, meas, dual_sol)` only when matching checked bounds establish one
relaxation value, a residual-checked hierarchy POVM is available, and the
completeness dual was reconstructed. A missing backend, resource limit,
unmatched bounds, invalid primal, or missing dual raises `DomainError` carrying
the structured result.
"""
function LocalDistinguishability(
    states,
    priors=nothing,
    dims=nothing,
    copies=2,
    ppt=1,
    bosonic=1,
    tolerance=eps(Float64)^(1 / 4);
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    dimensions = _compat_local_discrimination_dimensions(states, dims)
    checked_tolerance = _compat_separability_tolerance(
        tolerance, eps(Float64)^(1 / 4), "TOL"
    )
    result = local_distinguishability(
        states,
        dimensions;
        priors,
        order=_positive_dimension(copies, "COPIES"),
        ppt=_flag(ppt, "PPT"),
        bosonic=_flag(bosonic, "BOS"),
        backend,
        atol=checked_tolerance,
        rtol=zero(checked_tolerance),
        kwargs...,
    )
    structured && return result
    return _compat_local_distinguishability_output(result)
end

function _compat_is_separable_strength(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("STR must be a finite nonnegative integer"))
    value == -1 && throw(
        ArgumentError(
            "STR=-1 requests an unbounded search and is unsupported; choose a finite STR and explicit resource limits",
        ),
    )
    value >= 0 || throw(ArgumentError("STR must be a finite nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("STR is too large for Int"))
    return Int(value)
end

function _compat_is_separable_strategies(strength::Int)
    strategies = _COMPAT_IS_SEPARABLE_STRENGTH_ZERO
    strength >= 1 &&
        (strategies = (strategies..., :filter_covariance, :randomized_subtraction))
    strength >= 2 &&
        (strategies = (strategies..., :symmetric_extension, :symmetric_inner_extension))
    return strategies
end

function _compat_is_separable_output(result)
    _is_package_validated_entanglement_report(result) &&
    result.certified &&
    result.status in (:separable, :entangled) || throw(
        DomainError(
            result,
            "IsSeparable is inconclusive; request structured=true to inspect " *
            "necessary tests, heuristics, hierarchy results, and resource status",
        ),
    )
    return result.status === :separable ? 1 : 0
end

"""
    IsSeparable(
        rng::AbstractRNG, X, DIM=nothing, STR=2, VERBOSE=1,
        TOL=eps(Float64)^(3/8);
        structured=true, backend=NoOptimizationBackend(),
        max_extension_order=6, kwargs...
    )

Compatibility spelling for QETLAB's composite separability routine. The
mandatory leading RNG replaces the pinned routine's implicit global random
stream. `STR=0` runs the deterministic strength-zero family; finite
`STR >= 1` adds bounded product subtraction, and finite `STR >= 2` also
requests outer and inner hierarchy orders `2:STR`. The unbounded `STR=-1`
sentinel is rejected, and `STR` may not exceed the explicit
`max_extension_order` allocation guard.

The default returns an `EntanglementReport`. With `structured=false`, only a
validated separability certificate maps to integer `1`, and only a validated
entanglement certificate maps to integer `0`. An unknown, boundary, exhausted
heuristic, missing backend, or solver failure raises `DomainError`; it is
never fabricated into QETLAB's scalar `-1` or a Boolean answer.
"""
function IsSeparable(
    rng::AbstractRNG,
    operator::AbstractMatrix{<:Number},
    dims=nothing,
    strength=2,
    verbose=1,
    tolerance=eps(Float64)^(3 / 8);
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_extension_order=6,
    kwargs...,
)
    checked_strength = _compat_is_separable_strength(strength)
    checked_max_extension = _positive_dimension(max_extension_order, "max_extension_order")
    checked_strength <= checked_max_extension || throw(
        ArgumentError(
            "STR=$checked_strength exceeds max_extension_order=$checked_max_extension; raise the explicit cap only with corresponding optimization limits",
        ),
    )
    emit_message = _flag(verbose, "VERBOSE")
    checked_tolerance = _compat_separability_tolerance(
        tolerance, eps(Float64)^(3 / 8), "TOL"
    )
    orders = checked_strength >= 2 ? Tuple(2:checked_strength) : (2,)
    result = is_separable(
        rng,
        operator,
        _compat_separability_dimensions(operator, dims);
        strategies=_compat_is_separable_strategies(checked_strength),
        backend,
        extension_orders=orders,
        extension_ppt=true,
        extension_bosonic=true,
        max_extension_order=checked_max_extension,
        atol=checked_tolerance,
        rtol=zero(checked_tolerance),
        kwargs...,
    )
    emit_message && result.status !== :unknown && println(result.message)
    structured && return result
    return _compat_is_separable_output(result)
end

function _compat_upb_separable_discrimination_output(result)
    return throw(
        DomainError(
            result,
            "UPBSepDistinguishable cannot collapse the current floating-only " *
            "implementation to a Boolean; request structured=true to inspect " *
            "numerical feasibility, backend status, residuals, and inconclusive " *
            "floating-cone Farkas evidence",
        ),
    )
end

"""
    UPBSepDistinguishable(
        U, V, W...;
        structured=true, backend=NoOptimizationBackend(), kwargs...
    )

Compatibility spelling for QETLAB's perfect separable-discrimination test.
The default preserves the native status-rich result. The current public input
domain and replacement generator are floating-point only, so there is no exact
linked positive-certificate path; the floating-cone Farkas evidence likewise
does not provide a rigorous negative conclusion. Therefore
`structured=false` always raises `DomainError`.
"""
function UPBSepDistinguishable(
    first_factor::AbstractMatrix,
    second_factor::AbstractMatrix,
    remaining_factors::AbstractMatrix...;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    result = upb_sep_distinguishable(
        first_factor, second_factor, remaining_factors...; backend, kwargs...
    )
    structured && return result
    return _compat_upb_separable_discrimination_output(result)
end
