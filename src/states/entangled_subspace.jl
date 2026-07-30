# Source-informed independent Julia implementation based on the specification
# and QETLAB EntangledSubspace.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2022 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

const _ENTANGLED_SUBSPACE_DEFAULT_MAX_NONZEROS = 1_000_000
const _ENTANGLED_SUBSPACE_DEFAULT_MAX_WORK = 5_000_000

function _entangled_subspace_local_dims(local_dims)
    if local_dims isa Integer
        dimension = _positive_int(local_dims, "local_dims")
        return (dimension, dimension)
    elseif local_dims isa Tuple
        length(local_dims) == 2 ||
            throw(ArgumentError("local_dims must be a scalar or a two-entry collection"))
        return (
            _positive_int(local_dims[1], "local_dims[1]"),
            _positive_int(local_dims[2], "local_dims[2]"),
        )
    elseif local_dims isa AbstractVector
        Base.require_one_based_indexing(local_dims)
        length(local_dims) == 2 ||
            throw(ArgumentError("local_dims must be a scalar or a two-entry collection"))
        return (
            _positive_int(local_dims[1], "local_dims[1]"),
            _positive_int(local_dims[2], "local_dims[2]"),
        )
    end
    return throw(
        ArgumentError(
            "local_dims must be a positive integer or a two-entry tuple/vector; " *
            "got $(typeof(local_dims))",
        ),
    )
end

function _entangled_subspace_limit(value, name::AbstractString)
    value === nothing && return nothing
    return _positive_int(value, name)
end

function _entangled_subspace_value(::Type{T}, base::Int, exponent::Int) where {T}
    exact_value = big(base)^exponent
    value = try
        convert(T, exact_value)
    catch err
        err isa InexactError || err isa OverflowError || rethrow()
        throw(
            ArgumentError(
                "Vandermonde coefficient $base^$exponent cannot be represented " *
                "as $T; choose a wider coefficient_type",
            ),
        )
    end
    isfinite(value) || throw(
        ArgumentError(
            "Vandermonde coefficient is not finite in $T; choose a wider " *
            "coefficient_type",
        ),
    )
    return value
end

function _entangled_subspace_specs(dimensions::NTuple{2,Int}, r::Int, column_count::Int)
    first_dimension, second_dimension = dimensions
    smaller_dimension = min(dimensions...)
    specs = Vector{NTuple{3,Int}}()
    sizehint!(specs, column_count)
    for vandermonde_column in 1:(smaller_dimension - r)
        for diagonal in (r + 1 - second_dimension):(first_dimension - r - 1)
            diagonal_length = if diagonal >= 0
                min(second_dimension, first_dimension - diagonal)
            else
                min(second_dimension + diagonal, first_dimension)
            end
            if vandermonde_column <= diagonal_length - r
                push!(specs, (vandermonde_column, diagonal, diagonal_length))
                length(specs) == column_count && return specs
            end
        end
    end
    return error("internal entangled-subspace enumeration produced too few columns")
end

"""
    entangled_subspace(
        subspace_dimension,
        local_dims;
        r=1,
        coefficient_type=Float64,
        max_nonzeros=1_000_000,
        max_work=5_000_000,
    )

Construct a sparse basis for a bipartite `r`-entangled subspace.

`local_dims == (d₁, d₂)` uses subsystem `1` as the slowest-varying tensor
factor. Every nonzero vector in the returned column span has Schmidt rank at
least `r + 1`. Such a subspace can have at most
`(d₁-r) * (d₂-r)` dimensions, and the requested positive
`subspace_dimension` may be any value up to that sharp bound. A scalar
`local_dims` selects equal local dimensions. `r=0` is allowed and constructs a
subspace with no entanglement restriction.

The columns reproduce the pinned diagonal Vandermonde construction and are
linearly independent, but they are not normalized or orthogonal. The default
coefficient type is `Float64`, matching QETLAB's numeric output. Pass an
explicit concrete numeric type such as `BigInt`, `Rational{BigInt}`,
`Float32`, or `BigFloat` when exact or alternate-precision coefficients are
needed. No normalization or numerical rank repair is performed.

The output is a `SparseMatrixCSC` of size
`(d₁*d₂, subspace_dimension)`. `max_nonzeros` guards its exact stored-entry
count and `max_work` guards the same construction plus coefficient-power
work; either guard may be disabled explicitly with `nothing`.

# Examples

```jldoctest
julia> using LinearAlgebra, SparseArrays

julia> basis = entangled_subspace(2, (2, 3); r=1, coefficient_type=BigInt);

julia> size(basis), issparse(basis)
((6, 2), true)

julia> rank(Matrix{Float64}(reshape(basis[:, 1], 3, 2)))
2
```

The construction takes `O(max_nonzeros_required)` storage and at most
`O(max_work_required)` scalar multiplications.
"""
function entangled_subspace(
    subspace_dimension,
    local_dims;
    r=1,
    coefficient_type::Type{T}=Float64,
    max_nonzeros=_ENTANGLED_SUBSPACE_DEFAULT_MAX_NONZEROS,
    max_work=_ENTANGLED_SUBSPACE_DEFAULT_MAX_WORK,
) where {T}
    column_count = _positive_int(subspace_dimension, "subspace_dimension")
    dimensions = _entangled_subspace_local_dims(local_dims)
    entanglement_order = _nonnegative_int(r, "r")
    entanglement_order < min(dimensions...) || throw(
        ArgumentError(
            "r must be smaller than both local dimensions; got r=$entanglement_order " *
            "and local_dims=$dimensions",
        ),
    )
    T <: Number && isconcretetype(T) && T !== Bool || throw(
        ArgumentError(
            "coefficient_type must be a concrete numeric type other than Bool; got $T"
        ),
    )

    maximum_dimension = try
        Base.checked_mul(dimensions[1] - entanglement_order, dimensions[2] - entanglement_order)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the maximal subspace dimension exceeds typemax(Int)"))
    end
    column_count <= maximum_dimension || throw(
        ArgumentError(
            "no r-entangled subspace of dimension $column_count exists for " *
            "local_dims=$dimensions and r=$entanglement_order; the sharp maximum " *
            "is $maximum_dimension",
        ),
    )
    total_dimension = _checked_product(dimensions, "local_dims")
    nonzero_limit = _entangled_subspace_limit(max_nonzeros, "max_nonzeros")
    work_limit = _entangled_subspace_limit(max_work, "max_work")
    minimum_nonzeros = try
        Base.checked_mul(column_count, entanglement_order + 1)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the entangled-subspace minimum work exceeds Int"))
    end
    if nonzero_limit !== nothing && minimum_nonzeros > nonzero_limit
        throw(
            ArgumentError(
                "entangled_subspace requires at least $minimum_nonzeros stored " *
                "entries, exceeding max_nonzeros=$nonzero_limit",
            ),
        )
    end
    if work_limit !== nothing && minimum_nonzeros > work_limit
        throw(
            ArgumentError(
                "entangled_subspace requires at least $minimum_nonzeros scalar " *
                "operations, exceeding max_work=$work_limit",
            ),
        )
    end
    specs = _entangled_subspace_specs(dimensions, entanglement_order, column_count)
    nonzero_count = 0
    power_work = 0
    for (vandermonde_column, _, diagonal_length) in specs
        nonzero_count = try
            Base.checked_add(nonzero_count, diagonal_length)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("the entangled-subspace nonzero count exceeds Int"))
        end
        power_work = try
            Base.checked_add(
                power_work, Base.checked_mul(diagonal_length, vandermonde_column - 1)
            )
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("the entangled-subspace work estimate exceeds Int"))
        end
    end
    total_work = try
        Base.checked_add(nonzero_count, power_work)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the entangled-subspace work estimate exceeds Int"))
    end

    if nonzero_limit !== nothing && nonzero_count > nonzero_limit
        throw(
            ArgumentError(
                "entangled_subspace requires $nonzero_count stored entries, " *
                "exceeding max_nonzeros=$nonzero_limit",
            ),
        )
    end
    if work_limit !== nothing && total_work > work_limit
        throw(
            ArgumentError(
                "entangled_subspace requires an estimated $total_work scalar " *
                "operations, exceeding max_work=$work_limit",
            ),
        )
    end

    row_indices = Vector{Int}(undef, nonzero_count)
    column_indices = Vector{Int}(undef, nonzero_count)
    values = Vector{T}(undef, nonzero_count)
    cursor = 1
    second_dimension = dimensions[2]
    for (basis_column, spec) in pairs(specs)
        vandermonde_column, diagonal, diagonal_length = spec
        exponent = vandermonde_column - 1
        for position in 1:diagonal_length
            row, column = if diagonal >= 0
                (position, position + diagonal)
            else
                (position - diagonal, position)
            end
            row_indices[cursor] = row + (column - 1) * second_dimension
            column_indices[cursor] = basis_column
            values[cursor] = _entangled_subspace_value(T, position, exponent)
            cursor += 1
        end
    end
    return sparse(row_indices, column_indices, values, total_dimension, column_count)
end
