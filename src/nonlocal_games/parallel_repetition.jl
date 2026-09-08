# Source-informed independent Julia implementation based on the specification
# and QETLAB ParallelRepetition.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2022 Nathaniel Johnston and Mateus Araújo, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

const _PARALLEL_REPETITION_DEFAULT_MAX_ENTRIES = 10_000_000
const _PARALLEL_REPETITION_DEFAULT_MAX_WORK = 100_000_000

function _parallel_repetition_power(base::Int, exponent::Int, name::AbstractString)
    base == 1 && return 1
    return _checked_power(base, exponent, name)
end

function _parallel_repetition_limit(value, name::AbstractString)
    value === nothing && return nothing
    return _positive_int(value, name)
end

function _parallel_repetition_shape(
    game::AbstractArray{<:Number,4}, repetitions::Int; max_entries, max_work
)
    input_dims = size(game)
    all(>(0), input_dims) ||
        throw(ArgumentError("every output and setting dimension must be positive"))

    output_dims = ntuple(
        axis -> _parallel_repetition_power(
            input_dims[axis], repetitions, "size(game, $axis)^repetitions"
        ),
        4,
    )
    entries = _checked_product(output_dims, "parallel-repetition output dimensions")

    entry_limit = _parallel_repetition_limit(max_entries, "max_entries")
    if entry_limit !== nothing && entries > entry_limit
        throw(
            ArgumentError(
                "parallel repetition would allocate $entries entries, exceeding " *
                "max_entries=$entry_limit",
            ),
        )
    end

    work = try
        Base.checked_mul(entries, repetitions)
    catch err
        err isa OverflowError || rethrow()
        throw(
            ArgumentError(
                "parallel-repetition work estimate exceeds typemax(Int); " *
                "reduce repetitions or the game dimensions",
            ),
        )
    end
    work_limit = _parallel_repetition_limit(max_work, "max_work")
    if work_limit !== nothing && work > work_limit
        throw(
            ArgumentError(
                "parallel repetition requires an estimated $work scalar products, " *
                "exceeding max_work=$work_limit",
            ),
        )
    end

    return input_dims, output_dims
end

function _parallel_repetition_strides(input_dims::NTuple{4,Int}, repetitions::Int)
    strides = Matrix{Int}(undef, 4, repetitions)
    for axis in 1:4
        stride = 1
        for copy_index in repetitions:-1:1
            strides[axis, copy_index] = stride
            stride = try
                Base.checked_mul(stride, input_dims[axis])
            catch err
                err isa OverflowError || rethrow()
                throw(ArgumentError("size(game, $axis)^repetitions exceeds typemax(Int)"))
            end
        end
    end
    return strides
end

"""
    parallel_repetition(game, repetitions;
                        max_entries=10_000_000,
                        max_work=100_000_000)

Return the coefficient tensor for `repetitions` parallel copies of a bipartite
nonlocal game in full-probability notation.

`game[a, b, x, y]` is the coefficient for Alice output `a`, Bob output `b`,
Alice setting `x`, and Bob setting `y`. If the four one-copy dimensions are
`(oₐ, oᵦ, mₐ, mᵦ)`, the result has dimensions
`(oₐ^r, oᵦ^r, mₐ^r, mᵦ^r)`, where `r == repetitions`. Within each packed
index, copy `1` is the most-significant digit and copy `r` the
least-significant digit. Thus

```math
V^{\\otimes r}_{\\mathbf a,\\mathbf b,\\mathbf x,\\mathbf y}
= \\prod_{t=1}^r V_{a_t,b_t,x_t,y_t}.
```

The input must be a finite, one-based, four-dimensional numeric array with
nonempty axes, and `repetitions` must be positive. The output is dense because
Julia's standard sparse arrays are two-dimensional. `max_entries` guards the
allocation and `max_work` guards the estimated number of scalar
multiplications; either guard may be disabled explicitly with `nothing`.
Inputs are never normalized or otherwise repaired.
Fixed-width integer coefficients are multiplied with checked arithmetic and
raise `OverflowError` instead of wrapping.

# Examples

```jldoctest
julia> game = reshape(1:16, 2, 2, 2, 2);

julia> repeated = parallel_repetition(game, 2; max_entries=1_000);

julia> size(repeated)
(4, 4, 4, 4)

julia> repeated[2, 3, 1, 4] == game[1, 2, 1, 2] * game[2, 1, 1, 2]
true
```

The algorithm uses `O(r prod(size(game).^r))` scalar work and exactly one
dense output allocation of `prod(size(game).^r)` elements, in addition to
`O(r)` integer indexing data.

This is an independent Julia implementation of the public contract of QETLAB
`ParallelRepetition` at the pinned upstream revision. Unlike the pinned
routine's accidental behavior for `repetitions <= 0`, this API rejects
nonpositive copy counts.
"""
function parallel_repetition(
    game::AbstractArray{T,4},
    repetitions::Integer;
    max_entries=_PARALLEL_REPETITION_DEFAULT_MAX_ENTRIES,
    max_work=_PARALLEL_REPETITION_DEFAULT_MAX_WORK,
) where {T<:Number}
    Base.require_one_based_indexing(game)
    isconcretetype(T) ||
        throw(ArgumentError("game must have a concrete numeric element type; got $T"))
    all(isfinite, game) ||
        throw(ArgumentError("game must contain only finite coefficients"))

    copy_count = _positive_int(repetitions, "repetitions")
    input_dims, output_dims = _parallel_repetition_shape(
        game, copy_count; max_entries=max_entries, max_work=max_work
    )
    copy_count == 1 && return copy(game)

    packed_strides = _parallel_repetition_strides(input_dims, copy_count)
    repeated = Array{T}(undef, output_dims)
    for packed_index in CartesianIndices(repeated)
        coefficient = one(T)
        @inbounds for copy_index in 1:copy_count
            a =
                mod(
                    div(packed_index[1] - 1, packed_strides[1, copy_index]), input_dims[1]
                ) + 1
            b =
                mod(
                    div(packed_index[2] - 1, packed_strides[2, copy_index]), input_dims[2]
                ) + 1
            x =
                mod(
                    div(packed_index[3] - 1, packed_strides[3, copy_index]), input_dims[3]
                ) + 1
            y =
                mod(
                    div(packed_index[4] - 1, packed_strides[4, copy_index]), input_dims[4]
                ) + 1
            value = game[a, b, x, y]
            coefficient = if T <: Base.BitInteger
                _checked_bitinteger_mul(coefficient, value, "parallel_repetition")
            else
                coefficient * value
            end
        end
        @inbounds repeated[packed_index] = coefficient
    end
    return repeated
end

function parallel_repetition(game::AbstractArray{<:Number}, repetitions::Integer; kwargs...)
    return throw(
        DimensionMismatch(
            "game must be a four-dimensional full-probability tensor; got ndims=$(ndims(game))",
        ),
    )
end
