# Source-informed independent Julia implementation based on the specification
# and QETLAB PartialTrace.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    PartialTracePlan(dims, trace_out)

Precompute index projections for repeated partial traces with the same
subsystem dimensions and traced systems.

The order of `trace_out` is immaterial and is canonicalized.  Repeated or
out-of-range subsystem indices are rejected.
"""
struct PartialTracePlan{L<:SubsystemLayout,TO<:Tuple,K<:Tuple,O<:SubsystemLayout}
    layout::L
    trace_out::TO
    keep::K
    output_layout::O
    keep_index::_ReadOnlyPlanVector
    trace_index::_ReadOnlyPlanVector
    trace_dimension::Int
    source_indices::_ReadOnlyPlanMatrix

    function PartialTracePlan(
        token::_ValidatedConstructorToken,
        layout::L,
        trace_out::TO,
        keep::K,
        output_layout::O,
        keep_index::Vector{Int},
        trace_index::Vector{Int},
        trace_dimension::Int,
        source_indices::Matrix{Int},
    ) where {L<:SubsystemLayout,TO<:Tuple,K<:Tuple,O<:SubsystemLayout}
        _require_validated_constructor_token(token)
        length(keep_index) == layout.total_dimension ||
            throw(ArgumentError("validated keep-index lookup length is inconsistent"))
        length(trace_index) == layout.total_dimension ||
            throw(ArgumentError("validated trace-index lookup length is inconsistent"))
        trace_dimension > 0 ||
            throw(ArgumentError("validated trace dimension must be positive"))
        size(source_indices) == (output_layout.total_dimension, trace_dimension) ||
            throw(ArgumentError("validated source-index lookup shape is inconsistent"))
        for old_index in 1:layout.total_dimension
            keep_position = keep_index[old_index]
            trace_position = trace_index[old_index]
            1 <= keep_position <= output_layout.total_dimension ||
                throw(ArgumentError("validated keep-index lookup is out of bounds"))
            1 <= trace_position <= trace_dimension ||
                throw(ArgumentError("validated trace-index lookup is out of bounds"))
            source_indices[keep_position, trace_position] == old_index ||
                throw(ArgumentError("validated source-index lookup is inconsistent"))
        end
        return new{L,TO,K,O}(
            layout,
            trace_out,
            keep,
            output_layout,
            _read_only_plan_array(keep_index),
            _read_only_plan_array(trace_index),
            trace_dimension,
            _read_only_plan_array(source_indices),
        )
    end
end

function PartialTracePlan(dims, trace_out)
    layout = _as_layout(dims)
    checked_trace = _normalize_systems(
        trace_out, length(layout); name="trace_out", sort_result=true
    )
    keep = _complement_systems(checked_trace, layout)
    keep_index, output_layout = _projected_indices(layout, keep)
    trace_index, trace_layout = _projected_indices(layout, checked_trace)
    source_indices = Matrix{Int}(
        undef, output_layout.total_dimension, trace_layout.total_dimension
    )
    @inbounds for old_index in 1:layout.total_dimension
        source_indices[keep_index[old_index], trace_index[old_index]] = old_index
    end
    return PartialTracePlan(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        layout,
        checked_trace,
        keep,
        output_layout,
        keep_index,
        trace_index,
        trace_layout.total_dimension,
        source_indices,
    )
end

function _partial_trace_vector(vector::AbstractVector, plan::PartialTracePlan)
    _validate_vector_dimension(vector, plan.layout)
    kept_dimension = plan.output_layout.total_dimension
    coefficients = Matrix{eltype(vector)}(undef, kept_dimension, plan.trace_dimension)
    @inbounds for old_index in eachindex(vector)
        coefficients[plan.keep_index[old_index], plan.trace_index[old_index]] = vector[old_index]
    end
    return coefficients * adjoint(coefficients)
end

function _partial_trace_vector(
    vector::SparseArrays.AbstractSparseVector, plan::PartialTracePlan
)
    _validate_vector_dimension(vector, plan.layout)
    indices, values = findnz(vector)
    coefficients = sparse(
        plan.keep_index[indices],
        plan.trace_index[indices],
        copy(values),
        plan.output_layout.total_dimension,
        plan.trace_dimension,
    )
    return coefficients * adjoint(coefficients)
end

function _partial_trace_matrix(matrix::AbstractMatrix, plan::PartialTracePlan)
    _validate_matrix_dimension(matrix, plan.layout)
    kept_dimension = plan.output_layout.total_dimension
    result = Matrix{eltype(matrix)}(undef, kept_dimension, kept_dimension)
    fill!(result, zero(first(matrix)))
    source = plan.source_indices
    @inbounds for traced_index in 1:plan.trace_dimension
        for kept_column in 1:kept_dimension
            old_column = source[kept_column, traced_index]
            for kept_row in 1:kept_dimension
                old_row = source[kept_row, traced_index]
                result[kept_row, kept_column] += matrix[old_row, old_column]
            end
        end
    end
    return result
end

function _partial_trace_matrix(matrix::SparseMatrixCSC, plan::PartialTracePlan)
    _validate_matrix_dimension(matrix, plan.layout)
    rows, columns, values = findnz(matrix)
    matching = plan.trace_index[rows] .== plan.trace_index[columns]
    return sparse(
        plan.keep_index[rows[matching]],
        plan.keep_index[columns[matching]],
        values[matching],
        plan.output_layout.total_dimension,
        plan.output_layout.total_dimension,
    )
end

"""
    partial_trace(x, plan)
    partial_trace(x, dims; trace_out)

Trace the selected subsystems from a pure-state vector or square operator.

For a vector `ψ`, this computes the reduction of `ψ * ψ'` without forming
that full outer product.  A nontrivial reduction always returns a matrix.
Tracing every subsystem returns an explicit `1 × 1` matrix (not a scalar),
and tracing no subsystems returns `ψ * ψ'` for a vector or a copy-equivalent
matrix for an operator.  Sparse inputs produce sparse outputs.
"""
function partial_trace(vector::AbstractVector, plan::PartialTracePlan)
    return _partial_trace_vector(vector, plan)
end

function partial_trace(matrix::AbstractMatrix, plan::PartialTracePlan)
    return _partial_trace_matrix(matrix, plan)
end

function partial_trace(x::Union{AbstractVector,AbstractMatrix}, dims; trace_out=nothing)
    layout = _as_layout(dims)
    selected = if trace_out === nothing
        isempty(layout.dims) ? () : (length(layout),)
    else
        trace_out
    end
    return partial_trace(x, PartialTracePlan(layout, selected))
end
