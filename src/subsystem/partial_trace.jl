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
    if issparse(vector)
        return _partial_trace_vector(sparse(vector), plan)
    end
    if isconcretetype(eltype(vector)) && eltype(vector) <: Base.BitInteger
        return _partial_trace_bitinteger_vector(vector, plan)
    end
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
    if isconcretetype(eltype(vector)) && eltype(vector) <: Base.BitInteger
        return _partial_trace_sparse_bitinteger_vector(vector, plan)
    end
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

function _partial_trace_bitinteger_vector(
    vector::AbstractVector{T}, plan::PartialTracePlan
) where {T<:Base.BitInteger}
    kept_dimension = plan.output_layout.total_dimension
    result = zeros(T, kept_dimension, kept_dimension)
    source = plan.source_indices
    @inbounds for kept_column in 1:kept_dimension
        for kept_row in 1:kept_dimension
            value = BigInt(0)
            for traced_index in 1:plan.trace_dimension
                old_column = source[kept_column, traced_index]
                old_row = source[kept_row, traced_index]
                value += BigInt(vector[old_row]) * BigInt(vector[old_column])
            end
            result[kept_row, kept_column] = _narrow_bitinteger(T, value, "partial_trace")
        end
    end
    return result
end

function _partial_trace_sparse_bitinteger_vector(
    vector::SparseArrays.AbstractSparseVector{T}, plan::PartialTracePlan
) where {T<:Base.BitInteger}
    indices, values = findnz(vector)
    entries_by_trace = Dict{Int,Vector{Tuple{Int,T}}}()
    for position in eachindex(values)
        value = values[position]
        iszero(value) && continue
        traced_index = plan.trace_index[indices[position]]
        entries = get!(entries_by_trace, traced_index, Tuple{Int,T}[])
        push!(entries, (plan.keep_index[indices[position]], value))
    end

    accumulated = Dict{Tuple{Int,Int},BigInt}()
    for traced_index in sort!(collect(keys(entries_by_trace)))
        entries = entries_by_trace[traced_index]
        for (kept_column, column_value) in entries
            for (kept_row, row_value) in entries
                key = (kept_row, kept_column)
                accumulated[key] =
                    get(accumulated, key, BigInt(0)) +
                    BigInt(row_value) * BigInt(column_value)
            end
        end
    end

    keys_in_column_order = sort!(collect(keys(accumulated)); by=key -> (key[2], key[1]))
    rows = Int[]
    columns = Int[]
    output_values = T[]
    for key in keys_in_column_order
        exact_value = accumulated[key]
        iszero(exact_value) && continue
        push!(rows, key[1])
        push!(columns, key[2])
        push!(output_values, _narrow_bitinteger(T, exact_value, "partial_trace"))
    end
    output_dimension = plan.output_layout.total_dimension
    return sparse(rows, columns, output_values, output_dimension, output_dimension)
end

# Boolean addition is not closed in `Bool`: Julia's arithmetic result type for
# `true + true` is `Int`. A partial trace is an additive reduction, so use that
# codomain for Boolean operators in both the dense and sparse kernels. All
# other element types retain the package's existing type-preserving policy.
_partial_trace_matrix_output_type(::Type{Bool}) = Int
_partial_trace_matrix_output_type(::Type{T}) where {T} = T

function _partial_trace_matrix(matrix::AbstractMatrix, plan::PartialTracePlan)
    _validate_matrix_dimension(matrix, plan.layout)
    if issparse(matrix)
        return _partial_trace_matrix(sparse(matrix), plan)
    end
    if isconcretetype(eltype(matrix)) && eltype(matrix) <: Base.BitInteger
        return _partial_trace_bitinteger_matrix(matrix, plan)
    end
    kept_dimension = plan.output_layout.total_dimension
    output_type = _partial_trace_matrix_output_type(eltype(matrix))
    result = Matrix{output_type}(undef, kept_dimension, kept_dimension)
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

function _partial_trace_bitinteger_matrix(
    matrix::AbstractMatrix{T}, plan::PartialTracePlan
) where {T<:Base.BitInteger}
    kept_dimension = plan.output_layout.total_dimension
    result = zeros(T, kept_dimension, kept_dimension)
    source = plan.source_indices
    @inbounds for kept_column in 1:kept_dimension
        for kept_row in 1:kept_dimension
            value = BigInt(0)
            for traced_index in 1:plan.trace_dimension
                old_column = source[kept_column, traced_index]
                old_row = source[kept_row, traced_index]
                value += BigInt(matrix[old_row, old_column])
            end
            result[kept_row, kept_column] = _narrow_bitinteger(T, value, "partial_trace")
        end
    end
    return result
end

function _partial_trace_matrix(matrix::SparseMatrixCSC, plan::PartialTracePlan)
    _validate_matrix_dimension(matrix, plan.layout)
    if isconcretetype(eltype(matrix)) && eltype(matrix) <: Base.BitInteger
        return _partial_trace_sparse_bitinteger_matrix(matrix, plan)
    end
    rows, columns, values = findnz(matrix)
    matching = plan.trace_index[rows] .== plan.trace_index[columns]
    output_type = _partial_trace_matrix_output_type(eltype(matrix))
    matching_values = values[matching]
    if output_type !== eltype(matrix)
        matching_values = output_type.(matching_values)
    end
    return sparse(
        plan.keep_index[rows[matching]],
        plan.keep_index[columns[matching]],
        matching_values,
        plan.output_layout.total_dimension,
        plan.output_layout.total_dimension,
    )
end

function _partial_trace_sparse_bitinteger_matrix(
    matrix::SparseMatrixCSC{T}, plan::PartialTracePlan
) where {T<:Base.BitInteger}
    rows, columns, values = findnz(matrix)
    entries_by_output = Dict{Tuple{Int,Int},Vector{Tuple{Int,T}}}()
    for position in eachindex(values)
        row = rows[position]
        column = columns[position]
        traced_index = plan.trace_index[row]
        traced_index == plan.trace_index[column] || continue
        key = (plan.keep_index[row], plan.keep_index[column])
        entries = get!(entries_by_output, key, Tuple{Int,T}[])
        push!(entries, (traced_index, values[position]))
    end

    keys_in_column_order = sort!(
        collect(keys(entries_by_output)); by=key -> (key[2], key[1])
    )
    output_rows = Int[]
    output_columns = Int[]
    output_values = T[]
    for key in keys_in_column_order
        entries = sort!(entries_by_output[key]; by=first)
        exact_value = BigInt(0)
        for (_, entry) in entries
            exact_value += BigInt(entry)
        end
        iszero(exact_value) && continue
        push!(output_rows, key[1])
        push!(output_columns, key[2])
        push!(output_values, _narrow_bitinteger(T, exact_value, "partial_trace"))
    end
    output_dimension = plan.output_layout.total_dimension
    return sparse(
        output_rows, output_columns, output_values, output_dimension, output_dimension
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
matrix for an operator.  Sparse inputs produce sparse outputs. Real
fixed-width integer reductions accumulate exactly, retain their input element
type when every final entry is representable, and raise `OverflowError`
instead of wrapping otherwise.
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
