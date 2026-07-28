# Source-informed independent Julia implementation based on the specification
# and QETLAB Realignment.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    RealignmentPlan(dims; systems=(1,))
    RealignmentPlan(row_dims, column_dims; systems=(1,))

Precompute a generalized operator-realignment map.  `systems` defines the
first party and its complement defines the second.  The order written in
`systems` is preserved and is observable for a party containing multiple
subsystems; complementary systems retain their original relative order.
For bipartite basis operators the map is

```math
|i j\\rangle\\langle k l| \\mapsto |i k\\rangle\\langle j l|.
```

Separate row and column layouts permit rectangular operators.  At least two
subsystems and a nonempty proper partition are required.
"""
struct RealignmentPlan{
    R<:SubsystemLayout,
    C<:SubsystemLayout,
    A<:Tuple,
    B<:Tuple,
    RAL<:SubsystemLayout,
    RBL<:SubsystemLayout,
    CAL<:SubsystemLayout,
    CBL<:SubsystemLayout,
}
    row_layout::R
    column_layout::C
    systems::A
    complement::B
    row_a_layout::RAL
    row_b_layout::RBL
    column_a_layout::CAL
    column_b_layout::CBL
    row_a_index::Vector{Int}
    row_b_index::Vector{Int}
    column_a_index::Vector{Int}
    column_b_index::Vector{Int}
    row_from_groups::Matrix{Int}
    column_from_groups::Matrix{Int}
    output_size::Tuple{Int,Int}
end

function RealignmentPlan(dims; systems=(1,))
    return RealignmentPlan(dims, dims; systems=systems)
end

function RealignmentPlan(row_dims, column_dims; systems=(1,))
    row_layout = _as_layout(row_dims)
    column_layout = _as_layout(column_dims)
    length(row_layout) == length(column_layout) || throw(
        DimensionMismatch(
            "row dims has $(length(row_layout)) subsystems but column dims has $(length(column_layout))",
        ),
    )
    length(row_layout) >= 2 ||
        throw(ArgumentError("realignment requires at least two subsystems"))
    selected = _normalize_systems(
        systems, length(row_layout); name="systems", allow_empty=false, proper=true
    )
    complement = _complement_systems(selected, row_layout)

    row_a_index, row_a_layout = _projected_indices(row_layout, selected)
    row_b_index, row_b_layout = _projected_indices(row_layout, complement)
    column_a_index, column_a_layout = _projected_indices(column_layout, selected)
    column_b_index, column_b_layout = _projected_indices(column_layout, complement)

    row_from_groups = Matrix{Int}(
        undef, row_a_layout.total_dimension, row_b_layout.total_dimension
    )
    for row in 1:row_layout.total_dimension
        row_from_groups[row_a_index[row], row_b_index[row]] = row
    end
    column_from_groups = Matrix{Int}(
        undef, column_a_layout.total_dimension, column_b_layout.total_dimension
    )
    for column in 1:column_layout.total_dimension
        column_from_groups[column_a_index[column], column_b_index[column]] = column
    end

    output_size = (
        Base.checked_mul(row_a_layout.total_dimension, column_a_layout.total_dimension),
        Base.checked_mul(row_b_layout.total_dimension, column_b_layout.total_dimension),
    )
    return RealignmentPlan(
        row_layout,
        column_layout,
        selected,
        complement,
        row_a_layout,
        row_b_layout,
        column_a_layout,
        column_b_layout,
        row_a_index,
        row_b_index,
        column_a_index,
        column_b_index,
        row_from_groups,
        column_from_groups,
        output_size,
    )
end

@inline function _realigned_row(plan::RealignmentPlan, old_row::Int, old_column::Int)
    return (plan.row_a_index[old_row] - 1) * plan.column_a_layout.total_dimension +
           plan.column_a_index[old_column]
end

@inline function _realigned_column(plan::RealignmentPlan, old_row::Int, old_column::Int)
    return (plan.row_b_index[old_row] - 1) * plan.column_b_layout.total_dimension +
           plan.column_b_index[old_column]
end

function _realign(matrix::AbstractMatrix, plan::RealignmentPlan)
    _validate_matrix_dimensions(matrix, plan.row_layout, plan.column_layout)
    result = Matrix{eltype(matrix)}(undef, plan.output_size)
    @inbounds for old_column in axes(matrix, 2), old_row in axes(matrix, 1)
        result[_realigned_row(plan, old_row, old_column), _realigned_column(plan, old_row, old_column)] = matrix[
            old_row, old_column
        ]
    end
    return result
end

function _realign(matrix::SparseMatrixCSC, plan::RealignmentPlan)
    _validate_matrix_dimensions(matrix, plan.row_layout, plan.column_layout)
    rows, columns, values = findnz(matrix)
    output_rows = similar(rows)
    output_columns = similar(columns)
    @inbounds for index in eachindex(values)
        output_rows[index] = _realigned_row(plan, rows[index], columns[index])
        output_columns[index] = _realigned_column(plan, rows[index], columns[index])
    end
    return sparse(output_rows, output_columns, copy(values), plan.output_size...)
end

function _inverse_realign(matrix::AbstractMatrix, plan::RealignmentPlan)
    size(matrix) == plan.output_size || throw(
        DimensionMismatch(
            "realigned matrix size $(size(matrix)) must be $(plan.output_size) for this plan",
        ),
    )
    result = Matrix{eltype(matrix)}(
        undef, plan.row_layout.total_dimension, plan.column_layout.total_dimension
    )
    @inbounds for old_column in 1:plan.column_layout.total_dimension,
        old_row in 1:plan.row_layout.total_dimension

        result[old_row, old_column] = matrix[
            _realigned_row(plan, old_row, old_column),
            _realigned_column(plan, old_row, old_column),
        ]
    end
    return result
end

function _inverse_realign(matrix::SparseMatrixCSC, plan::RealignmentPlan)
    size(matrix) == plan.output_size || throw(
        DimensionMismatch(
            "realigned matrix size $(size(matrix)) must be $(plan.output_size) for this plan",
        ),
    )
    rows, columns, values = findnz(matrix)
    inverse_rows = Vector{Int}(undef, length(values))
    inverse_columns = Vector{Int}(undef, length(values))

    column_a_dimension = plan.column_a_layout.total_dimension
    column_b_dimension = plan.column_b_layout.total_dimension
    @inbounds for index in eachindex(values)
        row_zero_based = rows[index] - 1
        column_zero_based = columns[index] - 1
        row_a = div(row_zero_based, column_a_dimension) + 1
        column_a = mod(row_zero_based, column_a_dimension) + 1
        row_b = div(column_zero_based, column_b_dimension) + 1
        column_b = mod(column_zero_based, column_b_dimension) + 1
        inverse_rows[index] = plan.row_from_groups[row_a, row_b]
        inverse_columns[index] = plan.column_from_groups[column_a, column_b]
    end
    return sparse(
        inverse_rows,
        inverse_columns,
        copy(values),
        plan.row_layout.total_dimension,
        plan.column_layout.total_dimension,
    )
end

"""
    realign(matrix, plan)
    realign(matrix, dims; systems=(1,))

Realign an operator according to [`RealignmentPlan`](@ref).  Sparse inputs
produce sparse outputs.
"""
realign(matrix::AbstractMatrix, plan::RealignmentPlan) = _realign(matrix, plan)

function realign(matrix::AbstractMatrix, dims; systems=(1,))
    return realign(matrix, RealignmentPlan(dims; systems=systems))
end

function realign(matrix::AbstractMatrix, row_dims, column_dims; systems=(1,))
    return realign(matrix, RealignmentPlan(row_dims, column_dims; systems=systems))
end

"""
    inverse_realign(matrix, plan)
    inverse_realign(matrix, dims; systems=(1,))

Invert [`realign`](@ref) exactly for the same dimensions and subsystem
partition.
"""
function inverse_realign(matrix::AbstractMatrix, plan::RealignmentPlan)
    return _inverse_realign(matrix, plan)
end

function inverse_realign(matrix::AbstractMatrix, dims; systems=(1,))
    return inverse_realign(matrix, RealignmentPlan(dims; systems=systems))
end

function inverse_realign(matrix::AbstractMatrix, row_dims, column_dims; systems=(1,))
    return inverse_realign(matrix, RealignmentPlan(row_dims, column_dims; systems=systems))
end

"""Alias for [`realign`](@ref)."""
realignment(args...; kwargs...) = realign(args...; kwargs...)

"""Alias for [`realign`](@ref), commonly called operator reshuffling."""
reshuffle(args...; kwargs...) = realign(args...; kwargs...)

"""Alias for [`inverse_realign`](@ref)."""
inverse_realignment(args...; kwargs...) = inverse_realign(args...; kwargs...)

"""Inverse of [`reshuffle`](@ref)."""
inverse_reshuffle(args...; kwargs...) = inverse_realign(args...; kwargs...)
