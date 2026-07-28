# Source-informed independent Julia implementation based on the specification
# and QETLAB PermuteSystems.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    SubsystemLayout(dims)

Validated dimensions and column-major basis strides for a multipartite space.

`dims` is a tuple or vector of positive integers.  Subsystems are written in
the same order as tensor-product factors: subsystem `1` is the
most-significant (slowest-varying) basis index, while the final subsystem is
the fastest-varying one.  Thus, for dimensions `(d₁, d₂)`,
`basis_to_linear((i, j), (d₁, d₂)) == (i - 1) * d₂ + j`.

The empty layout `SubsystemLayout(())` represents the one-dimensional scalar
tensor product and is used internally when every subsystem has been traced
out.
"""
struct SubsystemLayout{N}
    dims::NTuple{N,Int}
    strides::NTuple{N,Int}
    total_dimension::Int
end

function _positive_int(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer; got Bool $value"))
    value isa Integer ||
        throw(ArgumentError("$name must be a positive integer; got $(repr(value))"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    try
        return Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as Int"))
    end
end

function _nonnegative_int(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a nonnegative integer; got Bool $value"))
    value isa Integer ||
        throw(ArgumentError("$name must be a nonnegative integer; got $(repr(value))"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    try
        return Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as Int"))
    end
end

function _checked_product(values, name::AbstractString)
    result = 1
    for value in values
        try
            result = Base.checked_mul(result, value)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("$name has a product that exceeds typemax(Int)"))
        end
    end
    return result
end

function _checked_power(base::Int, exponent::Int, name::AbstractString)
    result = 1
    for _ in 1:exponent
        try
            result = Base.checked_mul(result, base)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("$name=$base^$exponent exceeds typemax(Int)"))
        end
    end
    return result
end

function SubsystemLayout(dims::Tuple)
    all(dim -> dim isa Integer && !(dim isa Bool), dims) ||
        throw(ArgumentError("dims must contain only positive integers; got $(repr(dims))"))
    checked_dims = ntuple(i -> _positive_int(dims[i], "dims[$i]"), length(dims))
    total = _checked_product(checked_dims, "dims")

    stride_values = Vector{Int}(undef, length(checked_dims))
    stride = 1
    for i in length(checked_dims):-1:1
        stride_values[i] = stride
        try
            stride = Base.checked_mul(stride, checked_dims[i])
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("dims has a product that exceeds typemax(Int)"))
        end
    end
    strides = ntuple(i -> stride_values[i], length(stride_values))
    return SubsystemLayout(checked_dims, strides, total)
end

SubsystemLayout(dims::AbstractVector{<:Integer}) = SubsystemLayout(Tuple(dims))
SubsystemLayout(layout::SubsystemLayout) = layout

Base.length(layout::SubsystemLayout) = length(layout.dims)
Base.getindex(layout::SubsystemLayout, index::Integer) = layout.dims[index]
Base.iterate(layout::SubsystemLayout, state...) = iterate(layout.dims, state...)
Base.Tuple(layout::SubsystemLayout) = layout.dims

function Base.show(io::IO, layout::SubsystemLayout)
    return print(io, "SubsystemLayout(", layout.dims, ")")
end

_as_layout(dims) = SubsystemLayout(dims)

function _tuple_from_int_vector(values::Vector{Int}, ::NTuple{N,Any}) where {N}
    return ntuple(index -> values[index], Val(N))
end

function _normalize_systems(
    systems,
    count::Int;
    name::AbstractString="systems",
    allow_empty::Bool=true,
    proper::Bool=false,
    sort_result::Bool=false,
)
    raw = if systems isa Integer
        (systems,)
    elseif systems isa Tuple || systems isa AbstractVector
        Tuple(systems)
    else
        throw(
            ArgumentError(
                "$name must be an integer or a tuple/vector of integers; got $(typeof(systems))",
            ),
        )
    end

    !allow_empty &&
        isempty(raw) &&
        throw(ArgumentError("$name must contain at least one subsystem"))
    proper &&
        length(raw) == count &&
        throw(ArgumentError("$name must select a proper subset of the $count subsystems"))

    normalized = Vector{Int}(undef, length(raw))
    seen = Set{Int}()
    for (position, system) in pairs(raw)
        system isa Bool &&
            throw(ArgumentError("$name[$position] must be an integer, not Bool"))
        system isa Integer ||
            throw(ArgumentError("$name[$position] must be an integer; got $(repr(system))"))
        index = try
            Int(system)
        catch err
            err isa InexactError || rethrow()
            throw(ArgumentError("$name[$position]=$system cannot be represented as Int"))
        end
        1 <= index <= count || throw(
            ArgumentError("$name[$position]=$index is outside the valid range 1:$count")
        )
        index in seen && throw(ArgumentError("$name contains repeated subsystem $index"))
        push!(seen, index)
        normalized[position] = index
    end
    sort_result && sort!(normalized)
    return _tuple_from_int_vector(normalized, raw)
end

function _normalize_permutation(permutation, count::Int)
    raw = if permutation isa Tuple || permutation isa AbstractVector
        Tuple(permutation)
    else
        throw(
            ArgumentError(
                "permutation must be a tuple or vector of integers; got $(typeof(permutation))",
            ),
        )
    end
    length(raw) == count || throw(
        ArgumentError(
            "permutation has length $(length(raw)); expected $count for the supplied dimensions",
        ),
    )
    normalized = _normalize_systems(raw, count; name="permutation", allow_empty=count == 0)
    length(normalized) == count ||
        throw(ArgumentError("permutation must contain every subsystem exactly once"))
    return normalized
end

function _complement_systems(systems::Tuple, count::Int)
    selected = Set(systems)
    return Tuple(system for system in 1:count if !(system in selected))
end

function _complement_systems(systems::NTuple{K,Int}, ::SubsystemLayout{N}) where {K,N}
    selected = Set(systems)
    values = Int[system for system in 1:N if !(system in selected)]
    return ntuple(index -> values[index], Val(N - K))
end

function _validate_vector_dimension(vector::AbstractVector, layout::SubsystemLayout)
    length(vector) == layout.total_dimension || throw(
        DimensionMismatch(
            "vector length $(length(vector)) does not match prod(dims)=$(layout.total_dimension) for dims=$(layout.dims)",
        ),
    )
    return nothing
end

function _validate_matrix_dimension(matrix::AbstractMatrix, layout::SubsystemLayout)
    expected = layout.total_dimension
    size(matrix) == (expected, expected) || throw(
        DimensionMismatch(
            "matrix size $(size(matrix)) must be ($expected, $expected) for dims=$(layout.dims)",
        ),
    )
    return nothing
end

function _validate_matrix_dimensions(
    matrix::AbstractMatrix, row_layout::SubsystemLayout, column_layout::SubsystemLayout
)
    expected = (row_layout.total_dimension, column_layout.total_dimension)
    size(matrix) == expected || throw(
        DimensionMismatch(
            "matrix size $(size(matrix)) must be $expected for row dims=$(row_layout.dims) and column dims=$(column_layout.dims)",
        ),
    )
    return nothing
end

function _coordinate(linear_zero_based::Int, layout::SubsystemLayout, system::Int)
    return mod(div(linear_zero_based, layout.strides[system]), layout.dims[system]) + 1
end

function _projected_indices(layout::SubsystemLayout, systems::NTuple{K,Int}) where {K}
    projected_dims = ntuple(position -> layout.dims[systems[position]], Val(K))
    projected_layout = SubsystemLayout(projected_dims)
    indices = Vector{Int}(undef, layout.total_dimension)
    for linear in 1:layout.total_dimension
        zero_based = linear - 1
        projected_zero_based = 0
        for (position, system) in pairs(systems)
            coordinate = _coordinate(zero_based, layout, system)
            projected_zero_based += (coordinate - 1) * projected_layout.strides[position]
        end
        indices[linear] = projected_zero_based + 1
    end
    return indices, projected_layout
end

"""
    basis_to_linear(indices, dims) -> Int

Convert one-based subsystem basis indices to a one-based linear index.

The first subsystem is slowest-varying.  Both `indices` and `dims` may be
tuples or vectors; a prevalidated [`SubsystemLayout`](@ref) may be supplied
in place of `dims`.
"""
function basis_to_linear(indices, dims)
    layout = _as_layout(dims)
    raw = if indices isa Tuple || indices isa AbstractVector
        Tuple(indices)
    else
        throw(
            ArgumentError(
                "indices must be a tuple or vector of integers; got $(typeof(indices))"
            ),
        )
    end
    length(raw) == length(layout) || throw(
        DimensionMismatch(
            "indices has length $(length(raw)); expected $(length(layout)) for dims=$(layout.dims)",
        ),
    )

    linear_zero_based = 0
    for position in eachindex(raw)
        index = raw[position]
        index isa Bool &&
            throw(ArgumentError("indices[$position] must be an integer, not Bool"))
        index isa Integer || throw(
            ArgumentError("indices[$position] must be an integer; got $(repr(index))")
        )
        coordinate = try
            Int(index)
        catch err
            err isa InexactError || rethrow()
            throw(ArgumentError("indices[$position]=$index cannot be represented as Int"))
        end
        1 <= coordinate <= layout.dims[position] || throw(
            ArgumentError(
                "indices[$position]=$coordinate is outside 1:$(layout.dims[position])"
            ),
        )
        linear_zero_based += (coordinate - 1) * layout.strides[position]
    end
    return linear_zero_based + 1
end

"""
    linear_to_basis(index, dims) -> Tuple

Convert a one-based linear basis index to one-based subsystem indices using
the tensor-factor order documented by [`SubsystemLayout`](@ref).
"""
function linear_to_basis(index, dims)
    layout = _as_layout(dims)
    index isa Bool && throw(ArgumentError("index must be an integer, not Bool"))
    index isa Integer ||
        throw(ArgumentError("index must be an integer; got $(repr(index))"))
    linear = try
        Int(index)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("index=$index cannot be represented as Int"))
    end
    1 <= linear <= layout.total_dimension || throw(
        ArgumentError(
            "index=$linear is outside 1:$(layout.total_dimension) for dims=$(layout.dims)",
        ),
    )

    zero_based = linear - 1
    return ntuple(system -> _coordinate(zero_based, layout, system), length(layout))
end
