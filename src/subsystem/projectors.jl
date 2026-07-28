# Source-informed independent Julia implementation based on the specification
# and QETLAB SymmetricProjection.m and AntisymmetricProjection.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

function _foreach_nondecreasing_sequence(visitor, local_dimension::Int, copies::Int)
    sequence = Vector{Int}(undef, copies)
    function visit(position::Int, minimum::Int)
        if position > copies
            visitor(copy(sequence))
            return nothing
        end
        for value in minimum:local_dimension
            sequence[position] = value
            visit(position + 1, value)
        end
    end
    visit(1, 1)
    return nothing
end

function _foreach_strictly_increasing_sequence(visitor, local_dimension::Int, copies::Int)
    sequence = Vector{Int}(undef, copies)
    function visit(position::Int, minimum::Int)
        if position > copies
            visitor(copy(sequence))
            return nothing
        end
        maximum = local_dimension - (copies - position)
        for value in minimum:maximum
            sequence[position] = value
            visit(position + 1, value + 1)
        end
    end
    visit(1, 1)
    return nothing
end

function _unique_permutations(values::Vector{Int})
    count = length(values)
    permutations = Vector{Vector{Int}}()
    current = Vector{Int}(undef, count)
    used = falses(count)
    function visit(position::Int)
        if position > count
            push!(permutations, copy(current))
            return nothing
        end
        previous = 0
        has_previous = false
        for index in eachindex(values)
            if !used[index] && (!has_previous || values[index] != previous)
                used[index] = true
                current[position] = values[index]
                visit(position + 1)
                used[index] = false
                previous = values[index]
                has_previous = true
            end
        end
    end
    visit(1)
    return permutations
end

function _permutation_sign(values::Vector{Int})
    inversions = 0
    for right in 2:length(values), left in 1:(right - 1)
        inversions += values[left] > values[right]
    end
    return iseven(inversions) ? 1 : -1
end

function _projection_arguments(local_dimension, copies)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _positive_int(copies, "copies")
    total_dimension = _checked_power(dimension, count, "local_dimension^copies")
    layout = SubsystemLayout(ntuple(_ -> dimension, count))
    return dimension, count, total_dimension, layout
end

function _finish_projection(matrix::SparseMatrixCSC, sparse_output::Bool)
    return sparse_output ? matrix : Matrix(matrix)
end

"""
    symmetric_subspace_basis(local_dimension, copies=2;
                             T=Float64, sparse_output=true)

Return a matrix whose columns are an orthonormal occupation-number basis for
the symmetric subspace of `copies` equal local systems.  The result has
`local_dimension^copies` rows and `binomial(local_dimension + copies - 1,
copies)` columns.
"""
function symmetric_subspace_basis(
    local_dimension, copies=2; T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    dimension, count, total_dimension, layout = _projection_arguments(
        local_dimension, copies
    )
    value_type = typeof(inv(sqrt(convert(T, 1))))
    rows = Int[]
    columns = Int[]
    values = value_type[]
    column = Ref(0)

    _foreach_nondecreasing_sequence(dimension, count) do sequence
        column[] += 1
        orbit = _unique_permutations(sequence)
        coefficient = convert(value_type, inv(sqrt(convert(T, length(orbit)))))
        for basis in orbit
            push!(rows, basis_to_linear(Tuple(basis), layout))
            push!(columns, column[])
            push!(values, coefficient)
        end
    end
    basis = sparse(rows, columns, values, total_dimension, column[])
    return _finish_projection(basis, sparse_output)
end

"""
    antisymmetric_subspace_basis(local_dimension, copies=2;
                                 T=Float64, sparse_output=true)

Return an orthonormal signed-permutation basis for the antisymmetric
subspace.  When `copies > local_dimension`, the result has zero columns.
"""
function antisymmetric_subspace_basis(
    local_dimension, copies=2; T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    dimension, count, total_dimension, layout = _projection_arguments(
        local_dimension, copies
    )
    value_type = typeof(inv(sqrt(convert(T, 1))))
    rows = Int[]
    columns = Int[]
    values = value_type[]
    column = Ref(0)

    if count <= dimension
        _foreach_strictly_increasing_sequence(dimension, count) do sequence
            column[] += 1
            orbit = _unique_permutations(sequence)
            coefficient = convert(value_type, inv(sqrt(convert(T, length(orbit)))))
            for basis in orbit
                push!(rows, basis_to_linear(Tuple(basis), layout))
                push!(columns, column[])
                push!(values, _permutation_sign(basis) * coefficient)
            end
        end
    end
    basis = sparse(rows, columns, values, total_dimension, column[])
    return _finish_projection(basis, sparse_output)
end

"""
    symmetric_projector(local_dimension, copies=2;
                        T=Float64, sparse_output=true)

Construct the orthogonal projector onto the symmetric subspace of
`copies` equal local systems.  The occupation-orbit construction avoids an
explicit sum over permutation matrices.  Choosing `T=Rational{Int}` yields
exact rational projector entries.
"""
function symmetric_projector(
    local_dimension, copies=2; T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    dimension, count, total_dimension, layout = _projection_arguments(
        local_dimension, copies
    )
    value_type = typeof(one(T) / 1)
    rows = Int[]
    columns = Int[]
    values = value_type[]

    _foreach_nondecreasing_sequence(dimension, count) do sequence
        orbit = _unique_permutations(sequence)
        basis_indices = [basis_to_linear(Tuple(basis), layout) for basis in orbit]
        coefficient = convert(value_type, one(T) / length(orbit))
        for column in basis_indices, row in basis_indices
            push!(rows, row)
            push!(columns, column)
            push!(values, coefficient)
        end
    end
    projector = sparse(rows, columns, values, total_dimension, total_dimension)
    return _finish_projection(projector, sparse_output)
end

"""
    antisymmetric_projector(local_dimension, copies=2;
                            T=Float64, sparse_output=true)

Construct the orthogonal projector onto the antisymmetric subspace of
`copies` equal local systems.  If `copies > local_dimension`, the
antisymmetric subspace is empty and a zero projector is returned.
"""
function antisymmetric_projector(
    local_dimension, copies=2; T::Type{<:Number}=Float64, sparse_output::Bool=true
)
    dimension, count, total_dimension, layout = _projection_arguments(
        local_dimension, copies
    )
    value_type = typeof(one(T) / 1)
    rows = Int[]
    columns = Int[]
    values = value_type[]

    if count <= dimension
        _foreach_strictly_increasing_sequence(dimension, count) do sequence
            orbit = _unique_permutations(sequence)
            basis_indices = [basis_to_linear(Tuple(basis), layout) for basis in orbit]
            signs = [_permutation_sign(basis) for basis in orbit]
            coefficient = convert(value_type, one(T) / length(orbit))
            for column_position in eachindex(basis_indices),
                row_position in eachindex(basis_indices)

                push!(rows, basis_indices[row_position])
                push!(columns, basis_indices[column_position])
                push!(values, signs[row_position] * signs[column_position] * coefficient)
            end
        end
    end
    projector = sparse(rows, columns, values, total_dimension, total_dimension)
    return _finish_projection(projector, sparse_output)
end

"""Alias for [`symmetric_projector`](@ref)."""
symmetric_projection(args...; kwargs...) = symmetric_projector(args...; kwargs...)

"""Alias for [`antisymmetric_projector`](@ref)."""
antisymmetric_projection(args...; kwargs...) = antisymmetric_projector(args...; kwargs...)
