# Source-informed independent Julia implementation based on the specification
# and QETLAB SymmetricProjection.m and AntisymmetricProjection.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

const _PROJECTION_DEFAULT_MAX_COLUMNS = 100_000
const _PROJECTION_DEFAULT_MAX_NONZEROS = 5_000_000
const _PROJECTION_DEFAULT_MAX_DENSE_ENTRIES = 10_000_000
const _PROJECTION_DEFAULT_MAX_WORK = 100_000_000

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

function _projection_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _projection_check_resource(
    planned::BigInt,
    limit,
    name::AbstractString,
    resource::AbstractString,
    function_name::AbstractString,
)
    limit === nothing && return nothing
    planned <= limit || throw(
        ArgumentError(
            "$function_name requires $planned $resource, exceeding $name=$limit; " *
            "raise the explicit guard only after reviewing the combinatorial cost",
        ),
    )
    return nothing
end

function _projection_symmetric_projector_nonzeros(dimension::Int, copies::Int)
    dimension == 1 && return BigInt(1)
    permutation_count = factorial(BigInt(copies))
    nonzeros = Ref(BigInt(0))
    _foreach_nondecreasing_sequence(dimension, copies) do sequence
        denominator = BigInt(1)
        run_length = 1
        for position in 2:copies
            if sequence[position] == sequence[position - 1]
                run_length += 1
            else
                denominator *= factorial(BigInt(run_length))
                run_length = 1
            end
        end
        denominator *= factorial(BigInt(run_length))
        orbit_size = div(permutation_count, denominator)
        return nonzeros[] += orbit_size * orbit_size
    end
    return nonzeros[]
end

function _projection_plan(
    local_dimension,
    copies;
    antisymmetric::Bool,
    projector::Bool,
    sparse_output::Bool,
    max_columns,
    max_nonzeros,
    max_dense_entries,
    max_work,
)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _positive_int(copies, "copies")
    column_limit = _projection_limit(max_columns, "max_columns")
    nonzero_limit = _projection_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _projection_limit(max_dense_entries, "max_dense_entries")
    work_limit = _projection_limit(max_work, "max_work")
    function_name = if antisymmetric
        projector ? "antisymmetric_projector" : "antisymmetric_subspace_basis"
    else
        projector ? "symmetric_projector" : "symmetric_subspace_basis"
    end

    # `_checked_power` performs `count` checked multiplications. Guard that
    # deterministic work before it starts, including the dimension-one case
    # where an arbitrarily large exponent would otherwise remain Int-sized.
    _projection_check_resource(
        BigInt(count),
        work_limit,
        "max_work",
        "checked-power multiplications",
        function_name,
    )
    # The one-dimensional tensor power is exactly one. Avoid a linear checked-
    # multiplication loop before the later quadratic-work preflight.
    total_dimension =
        dimension == 1 ? 1 : _checked_power(dimension, count, "local_dimension^copies")

    column_count, basis_nonzeros, permutation_factor = if antisymmetric
        if count > dimension
            (BigInt(0), BigInt(0), BigInt(0))
        else
            columns = binomial(BigInt(dimension), count)
            permutations = factorial(BigInt(count))
            (columns, columns * permutations, permutations)
        end
    else
        columns = binomial(BigInt(dimension) + BigInt(count) - 1, min(count, dimension - 1))
        (columns, BigInt(total_dimension), BigInt(0))
    end

    _projection_check_resource(
        column_count,
        column_limit,
        "max_columns",
        antisymmetric ? "antisymmetric basis columns" : "symmetric occupation columns",
        function_name,
    )
    if projector
        _projection_check_resource(
            BigInt(total_dimension),
            column_limit,
            "max_columns",
            "projector output columns",
            function_name,
        )
    end

    dense_entries = if projector
        BigInt(total_dimension)^2
    else
        BigInt(total_dimension) * column_count
    end
    if !sparse_output
        dense_entries <= typemax(Int) || throw(
            ArgumentError(
                "$function_name dense output requires $dense_entries entries, which " *
                "cannot be represented as an array length",
            ),
        )
        _projection_check_resource(
            dense_entries,
            dense_limit,
            "max_dense_entries",
            "dense output entries",
            function_name,
        )
    end

    # `_unique_permutations` scans up to `count` positions at every recursion
    # depth, and the antisymmetric path additionally computes every sign.
    permutation_work = basis_nonzeros * BigInt(count)^2
    antisymmetric && (permutation_work *= 2)
    preliminary_work = BigInt(count) + permutation_work + column_count * BigInt(count)
    _projection_check_resource(
        preliminary_work,
        work_limit,
        "max_work",
        "estimated scalar operations",
        function_name,
    )

    output_nonzeros = if !projector
        basis_nonzeros
    elseif antisymmetric
        basis_nonzeros * permutation_factor
    else
        _projection_symmetric_projector_nonzeros(dimension, count)
    end
    output_nonzeros <= typemax(Int) || throw(
        ArgumentError(
            "$function_name would generate $output_nonzeros stored entries, which " *
            "cannot be represented as an array length",
        ),
    )
    _projection_check_resource(
        output_nonzeros, nonzero_limit, "max_nonzeros", "stored entries", function_name
    )

    total_work = preliminary_work + (projector ? output_nonzeros : BigInt(0))
    _projection_check_resource(
        total_work, work_limit, "max_work", "estimated scalar operations", function_name
    )

    return (
        dimension=dimension,
        copies=count,
        total_dimension=total_dimension,
        columns=Int(column_count),
        nonzeros=Int(output_nonzeros),
        work=total_work,
    )
end

function _finish_projection(matrix::SparseMatrixCSC, sparse_output::Bool)
    return sparse_output ? matrix : Matrix(matrix)
end

"""
    symmetric_subspace_basis(
        local_dimension, copies=2;
        T=Float64, sparse_output=true,
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

Return a matrix whose columns are an orthonormal occupation-number basis for
the symmetric subspace of `copies` equal local systems.  The result has
`local_dimension^copies` rows and `binomial(local_dimension + copies - 1,
copies)` columns.

`max_columns`, `max_nonzeros`, and `max_work` bound the occupation sectors,
stored sparse entries, and a conservative permutation-work estimate before
enumeration. `max_dense_entries` additionally bounds a requested dense result.
Set an individual guard to `nothing` only after independently reviewing the
requested allocation; representability checks remain mandatory.
"""
function symmetric_subspace_basis(
    local_dimension,
    copies=2;
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
    max_columns=_PROJECTION_DEFAULT_MAX_COLUMNS,
    max_nonzeros=_PROJECTION_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_PROJECTION_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_PROJECTION_DEFAULT_MAX_WORK,
)
    plan = _projection_plan(
        local_dimension,
        copies;
        antisymmetric=false,
        projector=false,
        sparse_output=sparse_output,
        max_columns=max_columns,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
    dimension, count, total_dimension = plan.dimension, plan.copies, plan.total_dimension
    value_type = typeof(inv(sqrt(convert(T, 1))))
    if dimension == 1
        result = sparse([1], [1], [one(value_type)], 1, 1)
        return _finish_projection(result, sparse_output)
    end
    layout = SubsystemLayout(ntuple(_ -> dimension, count))
    rows = sizehint!(Int[], plan.nonzeros)
    columns = sizehint!(Int[], plan.nonzeros)
    values = sizehint!(value_type[], plan.nonzeros)
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
    antisymmetric_subspace_basis(
        local_dimension, copies=2;
        T=Float64, sparse_output=true,
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

Return an orthonormal signed-permutation basis for the antisymmetric
subspace.  When `copies > local_dimension`, the result has zero columns.
The four resource keywords have the same preflight and explicit-`nothing`
semantics as [`symmetric_subspace_basis`](@ref).
"""
function antisymmetric_subspace_basis(
    local_dimension,
    copies=2;
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
    max_columns=_PROJECTION_DEFAULT_MAX_COLUMNS,
    max_nonzeros=_PROJECTION_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_PROJECTION_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_PROJECTION_DEFAULT_MAX_WORK,
)
    plan = _projection_plan(
        local_dimension,
        copies;
        antisymmetric=true,
        projector=false,
        sparse_output=sparse_output,
        max_columns=max_columns,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
    dimension, count, total_dimension = plan.dimension, plan.copies, plan.total_dimension
    value_type = typeof(inv(sqrt(convert(T, 1))))
    rows = sizehint!(Int[], plan.nonzeros)
    columns = sizehint!(Int[], plan.nonzeros)
    values = sizehint!(value_type[], plan.nonzeros)
    column = Ref(0)

    if count <= dimension
        layout = SubsystemLayout(ntuple(_ -> dimension, count))
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
    symmetric_projector(
        local_dimension, copies=2;
        T=Float64, sparse_output=true,
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

Construct the orthogonal projector onto the symmetric subspace of
`copies` equal local systems.  The occupation-orbit construction avoids an
explicit sum over permutation matrices.  Choosing `T=Rational{Int}` yields
exact rational projector entries.

The resource guards bound occupation sectors, sparse output columns and
nonzeros, dense entries, and estimated enumeration work before result
allocation. Set an individual guard to `nothing` only after reviewing the
combinatorial cost.
"""
function symmetric_projector(
    local_dimension,
    copies=2;
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
    max_columns=_PROJECTION_DEFAULT_MAX_COLUMNS,
    max_nonzeros=_PROJECTION_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_PROJECTION_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_PROJECTION_DEFAULT_MAX_WORK,
)
    plan = _projection_plan(
        local_dimension,
        copies;
        antisymmetric=false,
        projector=true,
        sparse_output=sparse_output,
        max_columns=max_columns,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
    dimension, count, total_dimension = plan.dimension, plan.copies, plan.total_dimension
    value_type = typeof(one(T) / 1)
    if dimension == 1
        result = sparse([1], [1], [one(value_type)], 1, 1)
        return _finish_projection(result, sparse_output)
    end
    layout = SubsystemLayout(ntuple(_ -> dimension, count))
    rows = sizehint!(Int[], plan.nonzeros)
    columns = sizehint!(Int[], plan.nonzeros)
    values = sizehint!(value_type[], plan.nonzeros)

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
    antisymmetric_projector(
        local_dimension, copies=2;
        T=Float64, sparse_output=true,
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

Construct the orthogonal projector onto the antisymmetric subspace of
`copies` equal local systems.  If `copies > local_dimension`, the
antisymmetric subspace is empty and a zero projector is returned.
The four resource keywords have the same preflight and explicit-`nothing`
semantics as [`symmetric_projector`](@ref).
"""
function antisymmetric_projector(
    local_dimension,
    copies=2;
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
    max_columns=_PROJECTION_DEFAULT_MAX_COLUMNS,
    max_nonzeros=_PROJECTION_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_PROJECTION_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_PROJECTION_DEFAULT_MAX_WORK,
)
    plan = _projection_plan(
        local_dimension,
        copies;
        antisymmetric=true,
        projector=true,
        sparse_output=sparse_output,
        max_columns=max_columns,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
    dimension, count, total_dimension = plan.dimension, plan.copies, plan.total_dimension
    value_type = typeof(one(T) / 1)
    rows = sizehint!(Int[], plan.nonzeros)
    columns = sizehint!(Int[], plan.nonzeros)
    values = sizehint!(value_type[], plan.nonzeros)

    if count <= dimension
        layout = SubsystemLayout(ntuple(_ -> dimension, count))
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
