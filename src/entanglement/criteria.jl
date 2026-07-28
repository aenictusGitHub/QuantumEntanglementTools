# Source-informed independent Julia implementations based on the specifications
# in QETLAB IsPPT.m, Realignment.m, ReductionMap.m, and the corresponding
# necessary-condition checks in IsSeparable.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    CriterionStatus

Three-valued outcome vocabulary for necessary entanglement criteria:

- [`CriterionEntanglementDetected`](@ref): a robust violation certifies
  entanglement.
- [`CriterionSatisfied`](@ref): the necessary condition holds with a numerical
  margin, but this alone is not a separability certificate.
- [`CriterionUnknown`](@ref): the criterion or input lies within its numerical
  tolerance boundary.
"""
@enum CriterionStatus::UInt8 begin
    CriterionEntanglementDetected
    CriterionSatisfied
    CriterionUnknown
end

@doc "A robust necessary-criterion violation that certifies entanglement." CriterionEntanglementDetected
@doc "A necessary condition holds with margin; this is not a separability certificate." CriterionSatisfied
@doc "The criterion or input lies inside its numerical tolerance boundary." CriterionUnknown

"""
    CriterionResult

Structured result returned by [`ppt_criterion`](@ref),
[`realignment_criterion`](@ref), and [`reduction_criterion`](@ref).

`value` is the measured extremal eigenvalue or norm, `threshold` is its exact
mathematical boundary, and `tolerance` is the absolute numerical band used for
classification.  `witness` is criterion-specific and is `nothing` when no
explicit vector witness is supplied.  Unless supplied explicitly, `atol`
defaults to zero and `rtol` to `sqrt(eps(R))` for the input's real floating
type.
"""
struct CriterionResult{V,B,T,W}
    criterion::Symbol
    status::CriterionStatus
    value::V
    threshold::B
    tolerance::T
    witness::W
    message::String
end

function Base.show(io::IO, result::CriterionResult)
    return print(
        io,
        "CriterionResult(",
        result.criterion,
        ", ",
        result.status,
        ", value=",
        result.value,
        ", threshold=",
        result.threshold,
        ", tolerance=",
        result.tolerance,
        ")",
    )
end

function _tierd_criterion_status_lower(value, boundary, tolerance)
    if value < boundary - tolerance
        return CriterionEntanglementDetected
    elseif value > boundary + tolerance
        return CriterionSatisfied
    end
    return CriterionUnknown
end

function _tierd_criterion_status_upper(value, boundary, tolerance)
    if value > boundary + tolerance
        return CriterionEntanglementDetected
    elseif value < boundary - tolerance
        return CriterionSatisfied
    end
    return CriterionUnknown
end

function _tierd_criterion_message(
    criterion::Symbol, status::CriterionStatus; input_boundary_uncertain::Bool=false
)
    input_boundary_uncertain &&
        return "the input has a nonzero Hermiticity, trace, or positivity residual inside the state-validation tolerance, so no certificate is reported"
    status === CriterionEntanglementDetected &&
        return "$criterion violation certifies entanglement"
    status === CriterionSatisfied &&
        return "$criterion necessary condition holds with a tolerance margin; this test alone gives no entanglement certificate"
    return "$criterion value lies within the numerical tolerance boundary"
end

function _tierd_criterion_layout(dims, dimension::Int; bipartite::Bool=false)
    layout = _as_layout(dims)
    layout.total_dimension == dimension || throw(
        DimensionMismatch(
            "prod(dims)=$(layout.total_dimension) does not match state dimension $dimension",
        ),
    )
    if bipartite
        length(layout) == 2 || throw(
            ArgumentError("dims must describe exactly two subsystems; got $(layout.dims)"),
        )
    else
        length(layout) >= 2 ||
            throw(ArgumentError("the criterion requires at least two subsystems"))
    end
    return layout
end

function _tierd_criterion_analysis(
    rho::AbstractMatrix{<:Number};
    atol,
    rtol,
    allow_densify::Bool,
    operation::AbstractString,
)
    return _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation=operation,
    )
end

"""
    ppt_criterion(rho, dims; systems=(2,), atol=nothing, rtol=nothing,
                  allow_densify=false) -> CriterionResult

Apply the positive-partial-transpose necessary condition to a validated
multipartite density matrix.  `systems` must be a nonempty proper subset of
the subsystem indices.  A minimum partial-transpose eigenvalue below
`-tolerance` returns `CriterionEntanglementDetected`; one above `tolerance`
returns `CriterionSatisfied`; the boundary band returns `CriterionUnknown`.

A detected result includes the minimum-eigenvalue vector as `witness`, so
`witness' * partial_transpose(rho, dims; systems=systems) * witness < 0`
up to numerical error.  Passing this necessary condition is never reported as
a separability conclusion.  Sparse inputs require explicit
`allow_densify=true`.  Dense state validation and the partial-transpose
eigendecomposition cost `O(n^3)` time and `O(n^2)` workspace.
"""
function ppt_criterion(
    rho::AbstractMatrix{<:Number},
    dims;
    systems=(2,),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    analysis = _tierd_criterion_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="ppt_criterion"
    )
    layout = _tierd_criterion_layout(dims, size(analysis.matrix, 1))
    selected = _normalize_systems(
        systems,
        length(layout);
        name="systems",
        allow_empty=false,
        proper=true,
        sort_result=true,
    )
    transpose_plan = PartialTransposePlan(layout, selected)
    transposed = partial_transpose(analysis.matrix, transpose_plan)
    work_matrix = (transposed + adjoint(transposed)) / 2
    decomposition = eigen(Hermitian(work_matrix))
    index = argmin(decomposition.values)
    value = decomposition.values[index]
    scale = max(one(value), maximum(abs, decomposition.values; init=zero(value)))
    tolerance = _tierd_threshold(scale, analysis.atol, analysis.rtol)
    status = _tierd_criterion_status_lower(value, zero(value), tolerance)
    analysis.structural_boundary_uncertain && (status = CriterionUnknown)
    witness = if status === CriterionEntanglementDetected || status === CriterionUnknown
        decomposition.vectors[:, index]
    else
        nothing
    end
    return CriterionResult(
        :ppt,
        status,
        value,
        zero(value),
        tolerance,
        witness,
        _tierd_criterion_message(
            :ppt, status; input_boundary_uncertain=analysis.structural_boundary_uncertain
        ),
    )
end

"""
    realignment_criterion(rho, dims; systems=(1,), atol=nothing,
                          rtol=nothing, allow_densify=false)

Apply the computable cross-norm/realignment necessary condition.  For a
positive operator the exact boundary is `trace(rho)`; a realignment trace norm
larger than that boundary by more than the numerical tolerance certifies
entanglement.  A value below the boundary with margin only returns
`CriterionSatisfied`, never a separability conclusion.  Values in the
tolerance band return `CriterionUnknown`.  The full realignment SVD costs
`O(n^3)` time for balanced square bipartitions and uses dense `O(n^2)`
workspace.
"""
function realignment_criterion(
    rho::AbstractMatrix{<:Number},
    dims;
    systems=(1,),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    analysis = _tierd_criterion_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="realignment_criterion",
    )
    layout = _tierd_criterion_layout(dims, size(analysis.matrix, 1))
    selected = _normalize_systems(
        systems, length(layout); name="systems", allow_empty=false, proper=true
    )
    realignment_plan = RealignmentPlan(layout; systems=selected)
    realigned = realign(analysis.matrix, realignment_plan)
    value = sum(svdvals(realigned))
    boundary = analysis.trace_value
    scale = max(one(value), abs(value), abs(boundary))
    tolerance = _tierd_threshold(scale, analysis.atol, analysis.rtol)
    status = _tierd_criterion_status_upper(value, boundary, tolerance)
    analysis.structural_boundary_uncertain && (status = CriterionUnknown)
    return CriterionResult(
        :realignment,
        status,
        value,
        boundary,
        tolerance,
        nothing,
        _tierd_criterion_message(
            :realignment,
            status;
            input_boundary_uncertain=analysis.structural_boundary_uncertain,
        ),
    )
end

function _tierd_reduction_operator(
    rho::AbstractMatrix, layout::SubsystemLayout, side::Symbol
)
    dimension_a, dimension_b = layout.dims
    if side === :a
        rho_b = partial_trace(rho, layout; trace_out=(1,))
        identity_a = Matrix{eltype(rho)}(I, dimension_a, dimension_a)
        return tensor_product(identity_a, rho_b) - rho
    elseif side === :b
        rho_a = partial_trace(rho, layout; trace_out=(2,))
        identity_b = Matrix{eltype(rho)}(I, dimension_b, dimension_b)
        return tensor_product(rho_a, identity_b) - rho
    end
    return throw(ArgumentError("side must be :a or :b; got $(repr(side))"))
end

"""
    reduction_criterion(rho, dims; side=:both, atol=nothing, rtol=nothing,
                        allow_densify=false)

Apply the reduction-map necessary condition to a bipartite density matrix.
`side=:a` tests `I_A ⊗ rho_B - rho`, `side=:b` tests
`rho_A ⊗ I_B - rho`, and `side=:both` tests both.  A negative minimum
eigenvalue beyond tolerance certifies entanglement.  If every tested operator
is positive with a tolerance margin the result is `CriterionSatisfied`; a
boundary value returns `CriterionUnknown`.

The `witness` of a detected or boundary result is a named tuple with fields
`side` and `vector`.  Passing the reduction condition is not a separability
certificate.  Each selected side requires a dense `n x n` Hermitian
eigendecomposition, hence `O(n^3)` time and `O(n^2)` workspace.
"""
function reduction_criterion(
    rho::AbstractMatrix{<:Number},
    dims;
    side::Symbol=:both,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    side in (:a, :b, :both) ||
        throw(ArgumentError("side must be :a, :b, or :both; got $(repr(side))"))
    analysis = _tierd_criterion_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="reduction_criterion",
    )
    layout = _tierd_criterion_layout(dims, size(analysis.matrix, 1); bipartite=true)
    sides = side === :both ? (:a, :b) : (side,)

    best_value = nothing
    best_side = first(sides)
    best_vector = nothing
    spectral_scale = one(eltype(analysis.eigenvalues))
    for selected_side in sides
        reduction = _tierd_reduction_operator(analysis.matrix, layout, selected_side)
        work_matrix = (reduction + adjoint(reduction)) / 2
        decomposition = eigen(Hermitian(work_matrix))
        index = argmin(decomposition.values)
        value = decomposition.values[index]
        spectral_scale = max(
            spectral_scale, maximum(abs, decomposition.values; init=zero(value))
        )
        if best_value === nothing || value < best_value
            best_value = value
            best_side = selected_side
            best_vector = decomposition.vectors[:, index]
        end
    end

    tolerance = _tierd_threshold(spectral_scale, analysis.atol, analysis.rtol)
    status = _tierd_criterion_status_lower(best_value, zero(best_value), tolerance)
    analysis.structural_boundary_uncertain && (status = CriterionUnknown)
    witness = if status === CriterionEntanglementDetected || status === CriterionUnknown
        (side=best_side, vector=best_vector)
    else
        nothing
    end
    return CriterionResult(
        :reduction,
        status,
        best_value,
        zero(best_value),
        tolerance,
        witness,
        _tierd_criterion_message(
            :reduction,
            status;
            input_boundary_uncertain=analysis.structural_boundary_uncertain,
        ),
    )
end
