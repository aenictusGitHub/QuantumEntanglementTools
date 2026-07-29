# Source-informed independent Julia implementation based on the specification
# and QETLAB PartialTranspose.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    PartialTransposePlan(dims, systems)
    PartialTransposePlan(row_dims, column_dims, systems)

Precompute the index components for a partial transpose.  Separate row and
column dimensions support rectangular operators and unequal local
dimensions; selected local row and column dimensions exchange places in the
output.
"""
struct PartialTransposePlan{
    R<:SubsystemLayout,C<:SubsystemLayout,S<:Tuple,OR<:SubsystemLayout,OC<:SubsystemLayout
}
    row_layout::R
    column_layout::C
    systems::S
    output_row_layout::OR
    output_column_layout::OC
    row_to_row::_ReadOnlyPlanVector
    row_to_column::_ReadOnlyPlanVector
    column_to_row::_ReadOnlyPlanVector
    column_to_column::_ReadOnlyPlanVector

    function PartialTransposePlan(
        token::_ValidatedConstructorToken,
        row_layout::R,
        column_layout::C,
        systems::S,
        output_row_layout::OR,
        output_column_layout::OC,
        row_to_row::Vector{Int},
        row_to_column::Vector{Int},
        column_to_row::Vector{Int},
        column_to_column::Vector{Int},
    ) where {
        R<:SubsystemLayout,
        C<:SubsystemLayout,
        S<:Tuple,
        OR<:SubsystemLayout,
        OC<:SubsystemLayout,
    }
        _require_validated_constructor_token(token)
        length(row_to_row) == row_layout.total_dimension ||
            throw(ArgumentError("validated row lookup length is inconsistent"))
        length(row_to_column) == row_layout.total_dimension ||
            throw(ArgumentError("validated row-column lookup length is inconsistent"))
        length(column_to_row) == column_layout.total_dimension ||
            throw(ArgumentError("validated column-row lookup length is inconsistent"))
        length(column_to_column) == column_layout.total_dimension ||
            throw(ArgumentError("validated column lookup length is inconsistent"))
        minimum(row_to_row) >= 0 &&
        minimum(column_to_row) >= 0 &&
        maximum(row_to_row) + maximum(column_to_row) < output_row_layout.total_dimension ||
            throw(ArgumentError("validated output-row lookups are out of bounds"))
        minimum(row_to_column) >= 0 &&
        minimum(column_to_column) >= 0 &&
        maximum(row_to_column) + maximum(column_to_column) <
        output_column_layout.total_dimension ||
            throw(ArgumentError("validated output-column lookups are out of bounds"))
        return new{R,C,S,OR,OC}(
            row_layout,
            column_layout,
            systems,
            output_row_layout,
            output_column_layout,
            _read_only_plan_array(row_to_row),
            _read_only_plan_array(row_to_column),
            _read_only_plan_array(column_to_row),
            _read_only_plan_array(column_to_column),
        )
    end
end

PartialTransposePlan(dims, systems) = PartialTransposePlan(dims, dims, systems)

function PartialTransposePlan(row_dims, column_dims, systems)
    row_layout = _as_layout(row_dims)
    column_layout = _as_layout(column_dims)
    length(row_layout) == length(column_layout) || throw(
        DimensionMismatch(
            "row dims has $(length(row_layout)) subsystems but column dims has $(length(column_layout))",
        ),
    )
    selected = _normalize_systems(
        systems, length(row_layout); name="systems", sort_result=true
    )
    selected_set = Set(selected)
    output_row_layout = SubsystemLayout(
        ntuple(system -> if system in selected_set
            column_layout.dims[system]
        else
            row_layout.dims[system]
        end, length(row_layout))
    )
    output_column_layout = SubsystemLayout(
        ntuple(system -> if system in selected_set
            row_layout.dims[system]
        else
            column_layout.dims[system]
        end, length(row_layout))
    )

    row_to_row = zeros(Int, row_layout.total_dimension)
    row_to_column = zeros(Int, row_layout.total_dimension)
    @inbounds for row in 1:row_layout.total_dimension
        zero_based = row - 1
        for system in 1:length(row_layout)
            coordinate = _coordinate(zero_based, row_layout, system) - 1
            if system in selected_set
                row_to_column[row] += coordinate * output_column_layout.strides[system]
            else
                row_to_row[row] += coordinate * output_row_layout.strides[system]
            end
        end
    end

    column_to_row = zeros(Int, column_layout.total_dimension)
    column_to_column = zeros(Int, column_layout.total_dimension)
    @inbounds for column in 1:column_layout.total_dimension
        zero_based = column - 1
        for system in 1:length(column_layout)
            coordinate = _coordinate(zero_based, column_layout, system) - 1
            if system in selected_set
                column_to_row[column] += coordinate * output_row_layout.strides[system]
            else
                column_to_column[column] +=
                    coordinate * output_column_layout.strides[system]
            end
        end
    end

    return PartialTransposePlan(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        row_layout,
        column_layout,
        selected,
        output_row_layout,
        output_column_layout,
        row_to_row,
        row_to_column,
        column_to_row,
        column_to_column,
    )
end

function _partial_transpose(matrix::AbstractMatrix, plan::PartialTransposePlan)
    _validate_matrix_dimensions(matrix, plan.row_layout, plan.column_layout)
    result = Matrix{eltype(matrix)}(
        undef,
        plan.output_row_layout.total_dimension,
        plan.output_column_layout.total_dimension,
    )
    @inbounds for old_column in axes(matrix, 2), old_row in axes(matrix, 1)
        new_row = plan.row_to_row[old_row] + plan.column_to_row[old_column] + 1
        new_column = plan.row_to_column[old_row] + plan.column_to_column[old_column] + 1
        result[new_row, new_column] = matrix[old_row, old_column]
    end
    return result
end

function _partial_transpose(matrix::SparseMatrixCSC, plan::PartialTransposePlan)
    _validate_matrix_dimensions(matrix, plan.row_layout, plan.column_layout)
    rows, columns, values = findnz(matrix)
    output_rows = plan.row_to_row[rows] .+ plan.column_to_row[columns] .+ 1
    output_columns = plan.row_to_column[rows] .+ plan.column_to_column[columns] .+ 1
    return sparse(
        output_rows,
        output_columns,
        copy(values),
        plan.output_row_layout.total_dimension,
        plan.output_column_layout.total_dimension,
    )
end

"""
    partial_transpose(matrix, plan)
    partial_transpose(matrix, dims; systems)
    partial_transpose(matrix, row_dims, column_dims; systems)

Transpose the local row/column indices of the selected subsystems.  No
complex conjugation is performed.  Applying the operation twice to the same
systems, with the first output dimensions used for the second plan, returns
the original matrix exactly.
"""
function partial_transpose(matrix::AbstractMatrix, plan::PartialTransposePlan)
    return _partial_transpose(matrix, plan)
end

function partial_transpose(matrix::AbstractMatrix, dims; systems)
    return partial_transpose(matrix, PartialTransposePlan(dims, systems))
end

function partial_transpose(matrix::AbstractMatrix, row_dims, column_dims; systems)
    return partial_transpose(matrix, PartialTransposePlan(row_dims, column_dims, systems))
end
