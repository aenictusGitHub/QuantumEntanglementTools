# Source-informed independent Julia implementation based on the specification
# and QETLAB PermuteSystems.m, Swap.m, PermutationOperator.m, and
# SwapOperator.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    SubsystemPermutationPlan(dims, permutation; inverse=false)

Precompute the basis-index map for repeatedly permuting subsystems.

`permutation[k]` names the old subsystem that occupies output position `k`,
matching `permutedims` and QETLAB's `PermuteSystems` convention.  The plan
stores only `prod(dims)` integer indices and is immutable apart from its
private lookup vector.
"""
struct SubsystemPermutationPlan{L<:SubsystemLayout,P<:Tuple,O<:SubsystemLayout}
    layout::L
    permutation::P
    output_layout::O
    forward::_ReadOnlyPlanVector

    function SubsystemPermutationPlan(
        token::_ValidatedConstructorToken,
        layout::L,
        permutation::P,
        output_layout::O,
        forward::Vector{Int},
    ) where {L<:SubsystemLayout,P<:Tuple,O<:SubsystemLayout}
        _require_validated_constructor_token(token)
        length(permutation) == length(layout) ||
            throw(ArgumentError("validated permutation length is inconsistent"))
        output_layout.dims == Tuple(layout.dims[system] for system in permutation) ||
            throw(ArgumentError("validated permutation output layout is inconsistent"))
        length(forward) == layout.total_dimension ||
            throw(ArgumentError("validated permutation lookup length is inconsistent"))
        sort(forward) == collect(1:output_layout.total_dimension) ||
            throw(ArgumentError("validated permutation lookup is not a bijection"))
        return new{L,P,O}(
            layout, permutation, output_layout, _read_only_plan_array(forward)
        )
    end
end

function SubsystemPermutationPlan(dims, permutation; inverse::Bool=false)
    layout = _as_layout(dims)
    checked_permutation = _normalize_permutation(permutation, length(layout))
    if inverse
        checked_permutation = Tuple(invperm(collect(checked_permutation)))
    end
    output_layout = SubsystemLayout(
        Tuple(layout.dims[system] for system in checked_permutation)
    )
    forward = Vector{Int}(undef, layout.total_dimension)
    for old_index in 1:layout.total_dimension
        old_basis = linear_to_basis(old_index, layout)
        new_basis = ntuple(
            position -> old_basis[checked_permutation[position]], length(layout)
        )
        forward[old_index] = basis_to_linear(new_basis, output_layout)
    end
    return SubsystemPermutationPlan(
        _VALIDATED_CONSTRUCTOR_TOKEN, layout, checked_permutation, output_layout, forward
    )
end

function _permute_vector(vector::AbstractVector, plan::SubsystemPermutationPlan)
    _validate_vector_dimension(vector, plan.layout)
    if issparse(vector)
        return _permute_vector(sparse(vector), plan)
    end
    result = Vector{eltype(vector)}(undef, length(vector))
    @inbounds for old_index in eachindex(plan.forward)
        result[plan.forward[old_index]] = vector[old_index]
    end
    return result
end

function _permute_vector(
    vector::SparseArrays.AbstractSparseVector, plan::SubsystemPermutationPlan
)
    _validate_vector_dimension(vector, plan.layout)
    indices, values = findnz(vector)
    return sparsevec(plan.forward[indices], copy(values), length(vector))
end

function _permute_matrix_rows(matrix::AbstractMatrix, row_plan::SubsystemPermutationPlan)
    Base.require_one_based_indexing(matrix)
    size(matrix, 1) == row_plan.layout.total_dimension || throw(
        DimensionMismatch(
            "matrix has $(size(matrix, 1)) rows; expected $(row_plan.layout.total_dimension) for row dims=$(row_plan.layout.dims)",
        ),
    )
    if issparse(matrix)
        return _permute_matrix_rows(sparse(matrix), row_plan)
    end
    result = Matrix{eltype(matrix)}(undef, size(matrix))
    @inbounds for column in axes(matrix, 2), old_row in axes(matrix, 1)
        result[row_plan.forward[old_row], column] = matrix[old_row, column]
    end
    return result
end

function _permute_matrix_rows(matrix::SparseMatrixCSC, row_plan::SubsystemPermutationPlan)
    Base.require_one_based_indexing(matrix)
    size(matrix, 1) == row_plan.layout.total_dimension || throw(
        DimensionMismatch(
            "matrix has $(size(matrix, 1)) rows; expected $(row_plan.layout.total_dimension) for row dims=$(row_plan.layout.dims)",
        ),
    )
    rows, columns, values = findnz(matrix)
    return sparse(
        row_plan.forward[rows], columns, copy(values), size(matrix, 1), size(matrix, 2)
    )
end

function _permute_matrix(
    matrix::AbstractMatrix,
    row_plan::SubsystemPermutationPlan,
    column_plan::SubsystemPermutationPlan,
)
    _validate_matrix_dimensions(matrix, row_plan.layout, column_plan.layout)
    if issparse(matrix)
        return _permute_matrix(sparse(matrix), row_plan, column_plan)
    end
    result = Matrix{eltype(matrix)}(undef, size(matrix))
    @inbounds for old_column in axes(matrix, 2), old_row in axes(matrix, 1)
        result[row_plan.forward[old_row], column_plan.forward[old_column]] = matrix[
            old_row, old_column
        ]
    end
    return result
end

function _permute_matrix(
    matrix::SparseMatrixCSC,
    row_plan::SubsystemPermutationPlan,
    column_plan::SubsystemPermutationPlan,
)
    _validate_matrix_dimensions(matrix, row_plan.layout, column_plan.layout)
    rows, columns, values = findnz(matrix)
    return sparse(
        row_plan.forward[rows],
        column_plan.forward[columns],
        copy(values),
        size(matrix, 1),
        size(matrix, 2),
    )
end

"""
    permute_subsystems(x, plan; rows_only=false)
    permute_subsystems(x, dims; permutation, rows_only=false, inverse=false)

Permute tensor factors of a state vector or operator.  For matrices, rows
and columns are permuted together unless `rows_only=true`.  Sparse vectors
and matrices produce sparse outputs.
"""
function permute_subsystems(
    vector::AbstractVector, plan::SubsystemPermutationPlan; rows_only::Bool=false
)
    rows_only && throw(ArgumentError("rows_only is not meaningful for an AbstractVector"))
    return _permute_vector(vector, plan)
end

function permute_subsystems(
    matrix::AbstractMatrix, plan::SubsystemPermutationPlan; rows_only::Bool=false
)
    return if rows_only
        _permute_matrix_rows(matrix, plan)
    else
        _permute_matrix(matrix, plan, plan)
    end
end

"""
    permute_subsystems(matrix, row_plan, column_plan)

Permute row and column tensor factors independently.  This method supports
rectangular operators and unequal local row/column dimensions.
"""
function permute_subsystems(
    matrix::AbstractMatrix,
    row_plan::SubsystemPermutationPlan,
    column_plan::SubsystemPermutationPlan,
)
    return _permute_matrix(matrix, row_plan, column_plan)
end

function permute_subsystems(
    x::Union{AbstractVector,AbstractMatrix},
    dims;
    permutation,
    rows_only::Bool=false,
    inverse::Bool=false,
)
    plan = SubsystemPermutationPlan(dims, permutation; inverse=inverse)
    return permute_subsystems(x, plan; rows_only=rows_only)
end

function permute_subsystems(
    x::Union{AbstractVector,AbstractMatrix},
    permutation,
    dims;
    rows_only::Bool=false,
    inverse::Bool=false,
)
    return permute_subsystems(
        x, dims; permutation=permutation, rows_only=rows_only, inverse=inverse
    )
end

"""
    swap_subsystems(x, dims, first, second; rows_only=false)
    swap_subsystems(x, dims; systems=(1, 2), rows_only=false)

Exchange two subsystem positions in a vector or operator.
"""
function swap_subsystems(
    x::Union{AbstractVector,AbstractMatrix}, dims, first, second; rows_only::Bool=false
)
    layout = _as_layout(dims)
    systems = _normalize_systems(
        (first, second), length(layout); name="systems", allow_empty=false
    )
    permutation = collect(1:length(layout))
    permutation[systems[1]], permutation[systems[2]] = permutation[systems[2]],
    permutation[systems[1]]
    return permute_subsystems(
        x, layout; permutation=Tuple(permutation), rows_only=rows_only
    )
end

"""
    permutation_operator(plan; T=Float64, sparse_output=true)
    permutation_operator(dims, permutation;
                         inverse=false, T=Float64, sparse_output=true)

Construct the linear operator `P` satisfying
`P * x == permute_subsystems(x, plan)`.  The sparse-first representation has
exactly `prod(dims)` stored entries.
"""
function permutation_operator(
    plan::SubsystemPermutationPlan; T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    dimension = plan.layout.total_dimension
    operator = sparse(
        plan.forward, collect(1:dimension), fill(one(T), dimension), dimension, dimension
    )
    return sparse_output ? operator : Matrix(operator)
end

function permutation_operator(
    dims,
    permutation;
    inverse::Bool=false,
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
)
    plan = SubsystemPermutationPlan(dims, permutation; inverse=inverse)
    return permutation_operator(plan; T=T, sparse_output=sparse_output)
end

"""
    swap_operator(dims; systems=(1, 2), T=Float64, sparse_output=true)

Construct the operator that swaps two subsystem positions.  Unlike QETLAB's
two-party `SwapOperator`, the native method also supports a selected pair
inside a larger multipartite layout.
"""
function swap_operator(
    dims; systems=(1, 2), T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    layout = _as_layout(dims)
    checked = _normalize_systems(systems, length(layout); name="systems", allow_empty=false)
    length(checked) == 2 ||
        throw(ArgumentError("systems must contain exactly two distinct subsystem indices"))
    permutation = collect(1:length(layout))
    permutation[checked[1]], permutation[checked[2]] = permutation[checked[2]],
    permutation[checked[1]]
    return permutation_operator(
        layout, Tuple(permutation); T=T, sparse_output=sparse_output
    )
end

function swap_subsystems(
    x::Union{AbstractVector,AbstractMatrix}, dims; systems=(1, 2), rows_only::Bool=false
)
    checked = if systems isa Tuple || systems isa AbstractVector
        Tuple(systems)
    else
        throw(ArgumentError("systems must contain exactly two subsystem indices"))
    end
    length(checked) == 2 || throw(
        ArgumentError(
            "systems must contain exactly two subsystem indices; got $(length(checked))"
        ),
    )
    return swap_subsystems(x, dims, checked[1], checked[2]; rows_only=rows_only)
end
