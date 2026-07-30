# Source-informed independent Julia implementation based on the specification
# and QETLAB OperatorSinkhorn.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export OperatorSinkhornResult, operator_sinkhorn

const _OPERATOR_SINKHORN_DEFAULT_MAX_ITERATIONS = 1_000
const _OPERATOR_SINKHORN_DEFAULT_MAX_ENTRIES = 10_000_000
const _OPERATOR_SINKHORN_DEFAULT_MAX_WORK = 1_000_000_000

"""
    OperatorSinkhornResult

Structured result returned by [`operator_sinkhorn`](@ref).

`scaled_operator` and `local_filters` obey

```math
\\mathtt{scaled\\_operator} \\approx
\\left(\\bigotimes_j F_j\\right)\\rho
\\left(\\bigotimes_j F_j\\right)^\\dagger.
```

`left_filter` and `right_filter` are aliases for the first two entries of
`local_filters`; `local_filters` is authoritative for multipartite input.
`status` is one of `:converged`, `:max_iterations`, `:singular_marginal`,
`:ill_conditioned`, `:work_limit`, or `:numerical_failure`. Only
`:converged` sets `converged=true`.

The aggregate `residual_history` contains the sum of the Frobenius distances
between every one-party marginal of the trace-one work matrix and its target
identity divided by local dimension. `final_marginal_residuals` exposes the
individual terms. Conditioning, minimum marginal eigenvalues, the largest
checked marginal Hermiticity residual, the failed subsystem, and conservative
work accounting make unsuccessful bounded iterations auditable.
"""
struct OperatorSinkhornResult{
    O<:AbstractMatrix,
    F<:Tuple,
    LF<:AbstractMatrix,
    RF<:AbstractMatrix,
    RH<:AbstractVector,
    MR<:AbstractVector,
    C<:AbstractVector,
    MC,
    E<:AbstractVector,
    H,
    T,
    CL,
    TR,
}
    scaled_operator::O
    local_filters::F
    left_filter::LF
    right_filter::RF
    iterations::Int
    status::Symbol
    converged::Bool
    residual_history::RH
    final_marginal_residuals::MR
    filter_condition_numbers::C
    maximum_filter_condition::MC
    minimum_marginal_eigenvalues::E
    maximum_marginal_hermiticity_residual::H
    convergence_threshold::T
    condition_limit::CL
    input_trace::TR
    failed_subsystem::Union{Nothing,Int}
    work_used::BigInt
    message::String
end

function Base.show(io::IO, result::OperatorSinkhornResult)
    return print(
        io,
        "OperatorSinkhornResult(status=",
        result.status,
        ", iterations=",
        result.iterations,
        ", residual=",
        last(result.residual_history),
        ", threshold=",
        result.convergence_threshold,
        ", max_condition=",
        result.maximum_filter_condition,
        ")",
    )
end

function _operator_sinkhorn_nonnegative_integer(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a nonnegative integer, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    return try
        Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as Int"))
    end
end

function _operator_sinkhorn_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _operator_sinkhorn_positive_real(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a finite real number, not Bool"))
    value isa Real || throw(ArgumentError("$name must be a finite real number"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value > zero(value) ||
        throw(ArgumentError("$name must be strictly positive; got $value"))
    return value
end

function _operator_sinkhorn_budget(
    dimension::Int, layout::SubsystemLayout, max_entries, max_work
)
    entry_limit = _operator_sinkhorn_limit(max_entries, "max_entries")
    work_limit = _operator_sinkhorn_limit(max_work, "max_work")
    dense_entries = BigInt(dimension)^2
    entry_limit !== nothing &&
        dense_entries > entry_limit &&
        throw(
            ArgumentError(
                "operator_sinkhorn requires a dense $dimension × $dimension work " *
                "matrix ($dense_entries entries), exceeding max_entries=$entry_limit",
            ),
        )
    snapshot_work =
        BigInt(length(layout)) * dense_entries +
        sum(BigInt(local_dimension)^3 for local_dimension in layout)
    setup_work = BigInt(dimension)^3 + snapshot_work
    work_limit !== nothing &&
        setup_work > work_limit &&
        throw(
            ArgumentError(
                "operator_sinkhorn validation and initial marginal work estimate " *
                "$setup_work exceeds max_work=$work_limit",
            ),
        )
    sweep_work =
        snapshot_work + sum(
            2 * BigInt(dimension)^3 + 2 * BigInt(local_dimension)^3 for
            local_dimension in layout
        )
    return (; work_limit, setup_work, sweep_work)
end

function _operator_sinkhorn_roundoff_tolerance(
    marginal::AbstractMatrix, ::Type{R}
) where {R<:AbstractFloat}
    scale = maximum(abs, marginal; init=zero(R))
    return 16 * size(marginal, 1) * eps(R) * max(one(R), scale)
end

function _operator_sinkhorn_marginal_data(
    marginal::AbstractMatrix, ::Type{R}
) where {R<:AbstractFloat}
    hermiticity_residual = maximum(abs, marginal - adjoint(marginal); init=zero(R))
    roundoff_tolerance = _operator_sinkhorn_roundoff_tolerance(marginal, R)
    hermiticity_residual <= roundoff_tolerance || return (
        ok=false,
        hermiticity_residual=hermiticity_residual,
        roundoff_tolerance=roundoff_tolerance,
        decomposition=nothing,
        minimum_eigenvalue=R(NaN),
        maximum_eigenvalue=R(NaN),
    )

    # The averaging is confined to a documented eigensolver work matrix after
    # checking a derived roundoff boundary. It never alters the input or the
    # returned scaled operator.
    work_matrix = (marginal + adjoint(marginal)) / 2
    decomposition = try
        eigen(Hermitian(work_matrix))
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        return (
            ok=false,
            hermiticity_residual=hermiticity_residual,
            roundoff_tolerance=roundoff_tolerance,
            decomposition=nothing,
            minimum_eigenvalue=R(NaN),
            maximum_eigenvalue=R(NaN),
        )
    end
    minimum_eigenvalue = minimum(decomposition.values)
    maximum_eigenvalue = maximum(decomposition.values)
    return (
        ok=true,
        hermiticity_residual=hermiticity_residual,
        roundoff_tolerance=roundoff_tolerance,
        decomposition=decomposition,
        minimum_eigenvalue=minimum_eigenvalue,
        maximum_eigenvalue=maximum_eigenvalue,
    )
end

function _operator_sinkhorn_snapshot(working, plans, targets, ::Type{R}) where {R}
    count = length(plans)
    residuals = Vector{R}(undef, count)
    minimum_eigenvalues = Vector{R}(undef, count)
    maximum_hermiticity_residual = zero(R)
    for subsystem in 1:count
        marginal = partial_trace(working, plans[subsystem])
        data = _operator_sinkhorn_marginal_data(marginal, R)
        maximum_hermiticity_residual = max(
            maximum_hermiticity_residual, data.hermiticity_residual
        )
        data.ok || return (
            ok=false,
            failed_subsystem=subsystem,
            residuals=fill(R(NaN), count),
            aggregate_residual=R(NaN),
            minimum_eigenvalues=fill(R(NaN), count),
            maximum_hermiticity_residual=maximum_hermiticity_residual,
        )
        residuals[subsystem] = norm(marginal - targets[subsystem])
        minimum_eigenvalues[subsystem] = data.minimum_eigenvalue
    end
    return (
        ok=true,
        failed_subsystem=nothing,
        residuals=residuals,
        aggregate_residual=sum(residuals; init=zero(R)),
        minimum_eigenvalues=minimum_eigenvalues,
        maximum_hermiticity_residual=maximum_hermiticity_residual,
    )
end

function _operator_sinkhorn_message(status::Symbol, failed_subsystem)
    status === :converged &&
        return "all one-party marginals satisfy the requested balancing tolerance"
    status === :max_iterations &&
        return "the bounded iteration reached max_iterations before convergence"
    status === :singular_marginal &&
        return "subsystem $failed_subsystem has a non-positive marginal eigenvalue"
    status === :ill_conditioned &&
        return "subsystem $failed_subsystem would exceed the filter condition limit"
    status === :work_limit &&
        return "the next complete Sinkhorn sweep would exceed max_work"
    return "a marginal or filter operation crossed a checked numerical roundoff boundary"
end

function _operator_sinkhorn_result(
    working,
    filters,
    input_trace,
    iterations,
    status,
    residual_history,
    snapshot,
    filter_conditions,
    maximum_filter_condition,
    encountered_minimum_eigenvalues,
    maximum_marginal_hermiticity_residual,
    convergence_threshold,
    condition_limit,
    failed_subsystem,
    work_used,
)
    final_minimum_eigenvalues = min.(
        encountered_minimum_eigenvalues, snapshot.minimum_eigenvalues
    )
    local_filters = Tuple(filters)
    return OperatorSinkhornResult(
        input_trace .* working,
        local_filters,
        local_filters[1],
        local_filters[2],
        iterations,
        status,
        status === :converged,
        copy(residual_history),
        copy(snapshot.residuals),
        copy(filter_conditions),
        maximum_filter_condition,
        final_minimum_eigenvalues,
        max(maximum_marginal_hermiticity_residual, snapshot.maximum_hermiticity_residual),
        convergence_threshold,
        condition_limit,
        input_trace,
        failed_subsystem,
        work_used,
        _operator_sinkhorn_message(status, failed_subsystem),
    )
end

"""
    operator_sinkhorn(
        rho, dims;
        atol=0,
        rtol=nothing,
        max_iterations=1000,
        max_condition_number=nothing,
        allow_densify=false,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    ) -> OperatorSinkhornResult

Run a bounded multipartite operator Sinkhorn iteration. `rho` must be a
finite, exactly Hermitian positive-semidefinite matrix with positive trace,
and `prod(dims) == size(rho, 1)`. At least two subsystems are required. The
input is never mutated, symmetrized, spectrally clipped, or returned with a
silently changed trace.

The algorithm explicitly forms the trace-one work matrix `rho / tr(rho)`,
records the original trace, and restores that trace algebraically in
`scaled_operator`. Consequently, the returned filters act directly on the
original `rho`. This internal gauge choice is documented rather than hidden
normalization.

Each complete sweep scales the current one-party marginal `M_j` with

```math
T_j = \\frac{M_j^{-1/2}}{\\sqrt{d_j}},
```

using a Hermitian eigendecomposition and scalar reciprocal square roots; no
matrix inverse or unrestricted matrix square root is formed. A non-positive
marginal eigenvalue returns `:singular_marginal`. A cumulative local-filter
condition number at or above `max_condition_number` returns
`:ill_conditioned`. The default condition limit is
`1 / sqrt(eps(R))`, matching the scale of the pinned routine's default
heuristic without consulting or changing global warning state.

The convergence residual is the sum of all one-party marginal Frobenius
distances from `I/d_j`. The threshold is
`atol + rtol * sum(1/sqrt(d_j))`; `rtol` defaults to `sqrt(eps(R))`.
`max_iterations` counts complete multipartite sweeps and may be zero.
`max_entries` guards every dense full-operator work matrix, while `max_work`
limits a conservative eigendecomposition, partial-trace, and dense
congruence estimate. Pass `nothing` to disable either budget.

Only `Float32`, `Float64`, `ComplexF32`, and `ComplexF64` matrices are
supported by the dependency-free spectral path. Sparse input requires
`allow_densify=true`; output and filters are dense. `BigFloat` is rejected
rather than down-converted.
"""
function operator_sinkhorn(
    rho::AbstractMatrix{<:Number},
    dims;
    atol=0,
    rtol=nothing,
    max_iterations=_OPERATOR_SINKHORN_DEFAULT_MAX_ITERATIONS,
    max_condition_number=nothing,
    allow_densify::Bool=false,
    max_entries=_OPERATOR_SINKHORN_DEFAULT_MAX_ENTRIES,
    max_work=_OPERATOR_SINKHORN_DEFAULT_MAX_WORK,
)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("operator_sinkhorn requires a square matrix"))
    size(rho, 1) > 0 || throw(ArgumentError("rho must have positive dimension"))
    layout = _as_layout(dims)
    length(layout) >= 2 ||
        throw(ArgumentError("operator_sinkhorn requires at least two subsystems"))
    _validate_matrix_dimension(rho, layout)

    checked_iterations = _operator_sinkhorn_nonnegative_integer(
        max_iterations, "max_iterations"
    )
    budget = _operator_sinkhorn_budget(size(rho, 1), layout, max_entries, max_work)
    dense, real_type = _tierd_dense_matrix(
        rho; allow_densify=allow_densify, operation="operator_sinkhorn"
    )
    ishermitian(dense) || throw(
        ArgumentError(
            "operator_sinkhorn requires an exactly Hermitian rho; the input " *
            "is never symmetrized",
        ),
    )
    input_decomposition = eigen(Hermitian(dense))
    input_minimum_eigenvalue = minimum(input_decomposition.values)
    input_minimum_eigenvalue >= zero(real_type) || throw(
        DomainError(
            input_minimum_eigenvalue,
            "rho must be positive semidefinite; negative eigenvalues are never clipped",
        ),
    )
    input_trace = real(tr(dense))
    isfinite(input_trace) ||
        throw(ArgumentError("rho has a non-finite trace after accumulation"))
    input_trace > zero(input_trace) ||
        throw(DomainError(input_trace, "rho must have positive trace"))

    checked_atol = _tierd_validate_tolerance(atol, "atol")
    checked_rtol =
        rtol === nothing ? sqrt(eps(real_type)) : _tierd_validate_tolerance(rtol, "rtol")
    absolute_tolerance = checked_atol === nothing ? zero(real_type) : checked_atol
    target_scale = sum(
        one(real_type) / sqrt(real_type(local_dimension)) for local_dimension in layout;
        init=zero(real_type),
    )
    convergence_threshold = absolute_tolerance + checked_rtol * target_scale
    isfinite(convergence_threshold) ||
        throw(ArgumentError("atol and rtol produce a non-finite convergence threshold"))
    condition_limit = if max_condition_number === nothing
        one(real_type) / sqrt(eps(real_type))
    else
        _operator_sinkhorn_positive_real(max_condition_number, "max_condition_number")
    end

    plans = [
        PartialTracePlan(
            layout,
            Tuple(
                subsystem for subsystem in 1:length(layout) if subsystem != kept_subsystem
            ),
        ) for kept_subsystem in 1:length(layout)
    ]
    targets = [
        Matrix{eltype(dense)}(I, local_dimension, local_dimension) /
        real_type(local_dimension) for local_dimension in layout
    ]
    identities = [
        Matrix{eltype(dense)}(I, local_dimension, local_dimension) for
        local_dimension in layout
    ]
    filters = copy(identities)
    filter_conditions = fill(one(real_type), length(layout))
    maximum_filter_condition = one(real_type)
    encountered_minimum_eigenvalues = fill(real_type(Inf), length(layout))
    maximum_marginal_hermiticity_residual = zero(real_type)
    working = dense / input_trace
    work_used = budget.setup_work

    snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
    maximum_marginal_hermiticity_residual = max(
        maximum_marginal_hermiticity_residual, snapshot.maximum_hermiticity_residual
    )
    snapshot.ok || return _operator_sinkhorn_result(
        working,
        filters,
        input_trace,
        0,
        :numerical_failure,
        [snapshot.aggregate_residual],
        snapshot,
        filter_conditions,
        maximum_filter_condition,
        encountered_minimum_eigenvalues,
        maximum_marginal_hermiticity_residual,
        convergence_threshold,
        condition_limit,
        snapshot.failed_subsystem,
        work_used,
    )
    encountered_minimum_eigenvalues .= min.(
        encountered_minimum_eigenvalues, snapshot.minimum_eigenvalues
    )
    residual_history = [snapshot.aggregate_residual]
    snapshot.aggregate_residual <= convergence_threshold &&
        return _operator_sinkhorn_result(
            working,
            filters,
            input_trace,
            0,
            :converged,
            residual_history,
            snapshot,
            filter_conditions,
            maximum_filter_condition,
            encountered_minimum_eigenvalues,
            maximum_marginal_hermiticity_residual,
            convergence_threshold,
            condition_limit,
            nothing,
            work_used,
        )

    completed_sweeps = 0
    while completed_sweeps < checked_iterations
        if budget.work_limit !== nothing &&
            work_used + budget.sweep_work > budget.work_limit
            return _operator_sinkhorn_result(
                working,
                filters,
                input_trace,
                completed_sweeps,
                :work_limit,
                residual_history,
                snapshot,
                filter_conditions,
                maximum_filter_condition,
                encountered_minimum_eigenvalues,
                maximum_marginal_hermiticity_residual,
                convergence_threshold,
                condition_limit,
                nothing,
                work_used,
            )
        end
        work_used += budget.sweep_work

        for subsystem in 1:length(layout)
            marginal = partial_trace(working, plans[subsystem])
            marginal_data = _operator_sinkhorn_marginal_data(marginal, real_type)
            maximum_marginal_hermiticity_residual = max(
                maximum_marginal_hermiticity_residual, marginal_data.hermiticity_residual
            )
            if !marginal_data.ok
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :numerical_failure,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end
            encountered_minimum_eigenvalues[subsystem] = min(
                encountered_minimum_eigenvalues[subsystem], marginal_data.minimum_eigenvalue
            )
            if marginal_data.minimum_eigenvalue <= zero(real_type)
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :singular_marginal,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end

            inverse_square_roots = map(
                value -> one(value) / sqrt(value), marginal_data.decomposition.values
            )
            if !all(isfinite, inverse_square_roots)
                maximum_filter_condition = real_type(Inf)
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :ill_conditioned,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end
            local_filter =
                (
                    marginal_data.decomposition.vectors *
                    Diagonal(inverse_square_roots) *
                    adjoint(marginal_data.decomposition.vectors)
                ) / sqrt(real_type(layout[subsystem]))
            candidate_filter = local_filter * filters[subsystem]
            all(isfinite, candidate_filter) || begin
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :numerical_failure,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end
            candidate_condition = try
                cond(candidate_filter)
            catch err
                err isa LinearAlgebra.LAPACKException || rethrow()
                real_type(Inf)
            end
            maximum_filter_condition = max(maximum_filter_condition, candidate_condition)
            if !isfinite(candidate_condition) || candidate_condition >= condition_limit
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :ill_conditioned,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end
            full_filter = tensor_product(
                ntuple(
                    index -> index == subsystem ? local_filter : identities[index],
                    length(layout),
                )...,
            )
            candidate_working = full_filter * working * adjoint(full_filter)
            all(isfinite, candidate_working) || begin
                snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
                return _operator_sinkhorn_result(
                    working,
                    filters,
                    input_trace,
                    completed_sweeps,
                    :numerical_failure,
                    residual_history,
                    snapshot,
                    filter_conditions,
                    maximum_filter_condition,
                    encountered_minimum_eigenvalues,
                    maximum_marginal_hermiticity_residual,
                    convergence_threshold,
                    condition_limit,
                    subsystem,
                    work_used,
                )
            end
            filters[subsystem] = candidate_filter
            filter_conditions[subsystem] = candidate_condition
            working = candidate_working
        end

        completed_sweeps += 1
        snapshot = _operator_sinkhorn_snapshot(working, plans, targets, real_type)
        maximum_marginal_hermiticity_residual = max(
            maximum_marginal_hermiticity_residual, snapshot.maximum_hermiticity_residual
        )
        snapshot.ok || return _operator_sinkhorn_result(
            working,
            filters,
            input_trace,
            completed_sweeps,
            :numerical_failure,
            residual_history,
            snapshot,
            filter_conditions,
            maximum_filter_condition,
            encountered_minimum_eigenvalues,
            maximum_marginal_hermiticity_residual,
            convergence_threshold,
            condition_limit,
            snapshot.failed_subsystem,
            work_used,
        )
        encountered_minimum_eigenvalues .= min.(
            encountered_minimum_eigenvalues, snapshot.minimum_eigenvalues
        )
        push!(residual_history, snapshot.aggregate_residual)
        snapshot.aggregate_residual <= convergence_threshold &&
            return _operator_sinkhorn_result(
                working,
                filters,
                input_trace,
                completed_sweeps,
                :converged,
                residual_history,
                snapshot,
                filter_conditions,
                maximum_filter_condition,
                encountered_minimum_eigenvalues,
                maximum_marginal_hermiticity_residual,
                convergence_threshold,
                condition_limit,
                nothing,
                work_used,
            )
    end

    return _operator_sinkhorn_result(
        working,
        filters,
        input_trace,
        completed_sweeps,
        :max_iterations,
        residual_history,
        snapshot,
        filter_conditions,
        maximum_filter_condition,
        encountered_minimum_eigenvalues,
        maximum_marginal_hermiticity_residual,
        convergence_threshold,
        condition_limit,
        nothing,
        work_used,
    )
end
