# Project-native occupation-coordinate tools for permutation-symmetric states.
# SPDX-License-Identifier: BSD-3-Clause

export symmetric_subspace_dimension,
    symmetric_occupations,
    symmetric_basis_index,
    symmetric_basis_occupation,
    generalized_dicke_state,
    symmetric_product_coordinates,
    symmetric_collective_operator,
    symmetric_split_isometry,
    symmetric_reduced_state,
    symmetric_maximally_mixed_state

const _SYMMETRIC_DEFAULT_MAX_OCCUPATIONS = 100_000
const _SYMMETRIC_DEFAULT_MAX_NONZEROS = 5_000_000
const _SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES = 10_000_000
const _SYMMETRIC_DEFAULT_MAX_WORK = 100_000_000
const _SYMMETRIC_MAX_TUPLE_ARITY = 1_000

function _symmetric_require_tuple_arity(local_dimension::Int, operation::AbstractString)
    local_dimension <= _SYMMETRIC_MAX_TUPLE_ARITY || throw(
        ArgumentError(
            "$operation requires occupation tuples with $local_dimension fields; " *
            "the supported safety limit is $_SYMMETRIC_MAX_TUPLE_ARITY",
        ),
    )
    return nothing
end

function _symmetric_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or `nothing`, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be a positive integer or `nothing`"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _symmetric_check_resource(
    planned::BigInt,
    limit,
    name::AbstractString,
    resource::AbstractString,
    operation::AbstractString,
)
    limit === nothing && return nothing
    planned <= limit || throw(
        ArgumentError(
            "$operation requires $planned $resource, exceeding $name=$limit; " *
            "raise the explicit guard only after reviewing the combinatorial cost",
        ),
    )
    return nothing
end

function _symmetric_array_length(value::BigInt, resource::AbstractString)
    value <= typemax(Int) || throw(
        ArgumentError("$resource=$value cannot be represented as a Julia array dimension"),
    )
    return Int(value)
end

function _symmetric_checked_occupation(occupation)
    occupation isa Tuple ||
        occupation isa AbstractVector ||
        throw(ArgumentError("occupation must be a nonempty tuple or vector of integers"))
    occupation isa AbstractArray && Base.require_one_based_indexing(occupation)
    isempty(occupation) &&
        throw(ArgumentError("occupation must contain at least one level"))
    _symmetric_require_tuple_arity(length(occupation), "occupation validation")

    checked = Vector{Int}(undef, length(occupation))
    total = BigInt(0)
    for index in eachindex(occupation)
        value = occupation[index]
        value isa Bool && throw(
            ArgumentError("occupation[$index] must be a nonnegative integer, not Bool")
        )
        value isa Integer || throw(
            ArgumentError(
                "occupation[$index] must be a nonnegative integer; got $(repr(value))"
            ),
        )
        value >= 0 ||
            throw(ArgumentError("occupation[$index] must be nonnegative; got $value"))
        try
            checked[index] = Int(value)
        catch err
            err isa InexactError || rethrow()
            throw(ArgumentError("occupation[$index]=$value cannot be represented as Int"))
        end
        total += value
    end
    total <= typemax(Int) ||
        throw(ArgumentError("the total occupation $total cannot be represented as Int"))
    return Tuple(checked), Int(total)
end

@inline function _symmetric_dimension(local_dimension::Int, parties::Int)
    return binomial(BigInt(local_dimension) + BigInt(parties) - 1, parties)
end

function _symmetric_foreach_occupation(visitor, local_dimension::Int, parties::Int)
    occupation = zeros(Int, local_dimension)
    occupation[1] = parties
    while true
        visitor(Tuple(occupation))

        pivot = 0
        for level in (local_dimension - 1):-1:1
            if occupation[level] > 0
                pivot = level
                break
            end
        end
        pivot == 0 && break

        redistributed = 1
        for level in (pivot + 1):local_dimension
            redistributed += occupation[level]
            occupation[level] = 0
        end
        occupation[pivot] -= 1
        occupation[pivot + 1] = redistributed
    end
    return nothing
end

function _symmetric_occupations_unchecked(local_dimension::Int, parties::Int)
    count = _symmetric_array_length(
        _symmetric_dimension(local_dimension, parties), "symmetric occupation count"
    )
    tuple_type = NTuple{local_dimension,Int}
    occupations = sizehint!(Vector{tuple_type}(), count)
    _symmetric_foreach_occupation(local_dimension, parties) do occupation
        return push!(occupations, occupation)
    end
    return occupations
end

"""
    symmetric_subspace_dimension(local_dimension, parties) -> BigInt

Return the exact dimension `binomial(local_dimension + parties - 1, parties)`
of `parties` identical `local_dimension`-level systems in the permutation-
symmetric subspace. `parties=0` uses the one-dimensional vacuum convention.

This project-native helper allocates no state or basis and cannot overflow an
`Int`; convert its result to `Int` only when an actual array is required.

# Example
```jldoctest
julia> symmetric_subspace_dimension(3, 4)
15
```

Complexity is that of one exact binomial coefficient, with `BigInt` output.
For astronomically large balanced arguments, the exact integer itself can be
expensive to compute and store; this scalar helper has no allocation override.
"""
function symmetric_subspace_dimension(local_dimension, parties)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    return _symmetric_dimension(dimension, count)
end

"""
    symmetric_occupations(
        local_dimension, parties;
        max_occupations=100_000, max_entries=10_000_000,
        max_work=100_000_000,
    )

Enumerate occupation tuples `(n₁, ..., n_d)` with nonnegative entries summing
to `parties`. Entry `n_j` counts local level `|j-1⟩`. The order exactly matches
the columns of [`symmetric_subspace_basis`](@ref): earlier occupations descend
first. For example, `(3, 2)` begins `(2,0,0)`, `(1,1,0)`, `(1,0,1)`.

The project-native implementation preflights the exact sector count, the total
number of stored tuple fields, and a conservative enumeration-work estimate.
Set a guard to `nothing` only after reviewing the requested allocation.

# Example
```jldoctest
julia> symmetric_occupations(3, 2)
6-element Vector{Tuple{Int64, Int64, Int64}}:
 (2, 0, 0)
 (1, 1, 0)
 (1, 0, 1)
 (0, 2, 0)
 (0, 1, 1)
 (0, 0, 2)
```

The result stores `binomial(d+N-1,N)` tuples and uses `O(d)` work per tuple.
Occupation-tuple arity is capped at 1,000 local levels for supported-Julia
compiler stability.
"""
function symmetric_occupations(
    local_dimension,
    parties;
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    _symmetric_require_tuple_arity(dimension, "symmetric_occupations")
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    entry_limit = _symmetric_limit(max_entries, "max_entries")
    work_limit = _symmetric_limit(max_work, "max_work")
    _symmetric_check_resource(
        BigInt(dimension) + BigInt(count) + 1,
        work_limit,
        "max_work",
        "preflight scalar operations",
        "symmetric_occupations",
    )
    sector_count = _symmetric_dimension(dimension, count)
    _symmetric_check_resource(
        sector_count,
        occupation_limit,
        "max_occupations",
        "occupation sectors",
        "symmetric_occupations",
    )
    _symmetric_check_resource(
        sector_count * BigInt(dimension),
        entry_limit,
        "max_entries",
        "stored occupation entries",
        "symmetric_occupations",
    )
    work = sector_count * (BigInt(dimension) + 1) + BigInt(count)
    _symmetric_check_resource(
        work, work_limit, "max_work", "estimated scalar operations", "symmetric_occupations"
    )
    return _symmetric_occupations_unchecked(dimension, count)
end

function _symmetric_basis_index_checked(occupation::Tuple, parties::Int)
    local_dimension = length(occupation)
    rank = BigInt(1)
    remaining = parties
    for level in 1:(local_dimension - 1)
        count = occupation[level]
        remaining_after = remaining - count
        remaining_levels = local_dimension - level
        top = BigInt(remaining_after) + BigInt(remaining_levels) - 1
        if top >= remaining_levels
            rank += binomial(top, remaining_levels)
        end
        remaining = remaining_after
    end
    return rank
end

"""
    symmetric_basis_index(occupation) -> BigInt

Return the exact one-based symmetric-basis column for an occupation tuple.
The ordering is the one returned by [`symmetric_occupations`](@ref) and used by
[`symmetric_subspace_basis`](@ref). Trailing zeros remain significant because
they determine the local dimension.

No enumeration is performed, and invalid, negative, nonintegral, or Boolean
occupations are rejected rather than repaired.

# Example
```jldoctest
julia> symmetric_basis_index((1, 0, 1))
3
```

Complexity is `O(d)` exact-binomial operations and `O(d)` validation storage.
"""
function symmetric_basis_index(occupation)
    checked, parties = _symmetric_checked_occupation(occupation)
    return _symmetric_basis_index_checked(checked, parties)
end

"""
    symmetric_basis_occupation(
        index, local_dimension, parties;
        max_entries=100_000, max_work=100_000_000,
    ) -> Tuple

Invert [`symmetric_basis_index`](@ref) without enumerating earlier sectors.
`index` is one-based and must not exceed
[`symmetric_subspace_dimension`](@ref)`(local_dimension, parties)`.

# Example
```jldoctest
julia> symmetric_basis_occupation(5, 3, 2)
(0, 1, 1)
```

The project-native unranking algorithm uses exact `BigInt` block counts. Its
worst-case work is `O(dN)` and it allocates only the returned tuple. The
`max_entries` and `max_work` guards are checked before that tuple is allocated.
"""
function symmetric_basis_occupation(
    index,
    local_dimension,
    parties;
    max_entries=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    index isa Bool && throw(ArgumentError("index must be a positive integer, not Bool"))
    index isa Integer ||
        throw(ArgumentError("index must be a positive integer; got $(repr(index))"))
    index > 0 || throw(ArgumentError("index must be positive; got $index"))
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    _symmetric_require_tuple_arity(dimension, "symmetric_basis_occupation")
    _symmetric_check_resource(
        BigInt(dimension),
        _symmetric_limit(max_entries, "max_entries"),
        "max_entries",
        "returned occupation entries",
        "symmetric_basis_occupation",
    )
    work = BigInt(dimension) * (BigInt(count) + 1) + 1
    _symmetric_check_resource(
        work,
        _symmetric_limit(max_work, "max_work"),
        "max_work",
        "estimated unranking operations",
        "symmetric_basis_occupation",
    )
    sector_count = _symmetric_dimension(dimension, count)
    exact_index = BigInt(index)
    exact_index <= sector_count || throw(BoundsError(BigInt(1):sector_count, exact_index))

    residual = exact_index - 1
    occupation = zeros(Int, dimension)
    remaining = count
    for level in 1:(dimension - 1)
        remaining_levels = dimension - level
        selected = -1
        for candidate in remaining:-1:0
            tail = remaining - candidate
            top = BigInt(tail) + BigInt(remaining_levels) - 1
            block = binomial(top, remaining_levels - 1)
            if residual < block
                selected = candidate
                break
            end
            residual -= block
        end
        selected >= 0 || error("internal symmetric-basis unranking failure")
        occupation[level] = selected
        remaining -= selected
    end
    occupation[end] = remaining
    return Tuple(occupation)
end

function _symmetric_multinomial(occupation::Tuple, parties::Int)
    result = BigInt(1)
    remaining = parties
    for level in 1:(length(occupation) - 1)
        count = occupation[level]
        result *= binomial(BigInt(remaining), count)
        remaining -= count
    end
    return result
end

function _symmetric_root_type(::Type{T}) where {T<:Number}
    return promote_type(T, typeof(sqrt(one(T))))
end

function _symmetric_ratio_sqrt(
    numerator::BigInt, denominator::BigInt, ::Type{T}
) where {T<:Number}
    denominator > 0 || throw(ArgumentError("the square-root denominator must be positive"))
    numerator >= 0 || throw(ArgumentError("the square-root numerator must be nonnegative"))
    output_type = _symmetric_root_type(T)
    value = sqrt(BigFloat(numerator) / BigFloat(denominator))
    converted = convert(output_type, value)
    if !isfinite(converted) || (numerator > 0 && iszero(converted))
        throw(
            OverflowError(
                "a square-root coefficient is not representable as $output_type; " *
                "use a wider numeric type",
            ),
        )
    end
    return converted
end

function _symmetric_inverse_sqrt(exact::BigInt, ::Type{T}) where {T<:Number}
    try
        root = _symmetric_ratio_sqrt(exact, BigInt(1), T)
        value = inv(root)
        if isfinite(value) && !iszero(value)
            return value
        end
    catch err
        err isa OverflowError || rethrow()
    end
    return _symmetric_ratio_sqrt(BigInt(1), exact, T)
end

function _symmetric_promoted_eltype(array)
    isempty(array) && throw(ArgumentError("input must be nonempty"))
    value_type = eltype(array)
    if value_type <: Number && isconcretetype(value_type)
        return value_type
    end
    values = _symmetric_stored_values(array)
    isempty(values) && return Float64
    first_value = first(values)
    first_value isa Number || throw(ArgumentError("input must contain only numbers"))
    promoted = typeof(first_value)
    for value in Iterators.drop(values, 1)
        value isa Number || throw(ArgumentError("input must contain only numbers"))
        promoted = promote_type(promoted, typeof(value))
    end
    return promoted
end

function _symmetric_is_diagonal_storage(array)
    return array isa Diagonal ||
           (array isa Union{Hermitian,Symmetric} && parent(array) isa Diagonal)
end

function _symmetric_stored_values(array)
    if _symmetric_is_diagonal_storage(array)
        return diag(array)
    elseif !issparse(array)
        return array
    elseif array isa Union{SparseVector,SparseMatrixCSC}
        return nonzeros(array)
    end
    return nonzeros(sparse(array))
end

function _symmetric_stored_count(array)
    return length(_symmetric_stored_values(array))
end

function _symmetric_input_entries(array)
    stored_layout = issparse(array) || _symmetric_is_diagonal_storage(array)
    return stored_layout ? BigInt(_symmetric_stored_count(array)) : BigInt(length(array))
end

function _symmetric_require_finite(array, operation::AbstractString)
    values = _symmetric_stored_values(array)
    all(value -> value isa Number && isfinite(value), values) ||
        throw(ArgumentError("$operation requires finite numeric input"))
    return nothing
end

function _symmetric_stored_bigfloat_precision(value)
    return 0
end

function _symmetric_stored_bigfloat_precision(value::BigFloat)
    return precision(value)
end

function _symmetric_stored_bigfloat_precision(value::Complex{BigFloat})
    return max(precision(real(value)), precision(imag(value)))
end

function _symmetric_precision_guard(value)
    return 0
end

function _symmetric_precision_guard(value::Integer)
    bits = iszero(value) ? 1 : ndigits(abs(BigInt(value)); base=2)
    return bits + 128
end

function _symmetric_precision_guard(value::Rational)
    numerator_bits = _symmetric_precision_guard(numerator(value))
    denominator_bits = _symmetric_precision_guard(denominator(value))
    return numerator_bits + denominator_bits + 128
end

function _symmetric_precision_guard(value::T) where {T<:AbstractFloat}
    iszero(value) && return precision(T) + 128
    return precision(T) + abs(exponent(value)) + 128
end

function _symmetric_precision_guard(value::BigFloat)
    return iszero(value) ? 128 : abs(exponent(value)) + 128
end

function _symmetric_precision_guard(value::Complex{T}) where {T<:Number}
    return max(
        _symmetric_precision_guard(real(value)), _symmetric_precision_guard(imag(value))
    )
end

function _symmetric_target_precision(array)
    result = precision(BigFloat)
    for value in _symmetric_stored_values(array)
        result = max(result, _symmetric_stored_bigfloat_precision(value))
    end
    return result
end

function _symmetric_binary_exponent(value::Real)
    iszero(value) && return nothing
    if value isa AbstractFloat
        return exponent(value)
    elseif value isa Integer
        return ndigits(abs(BigInt(value)); base=2) - 1
    elseif value isa Rational
        numerator_exponent = ndigits(abs(BigInt(numerator(value))); base=2) - 1
        denominator_exponent = ndigits(abs(BigInt(denominator(value))); base=2) - 1
        return numerator_exponent - denominator_exponent
    end
    return nothing
end

function _symmetric_binary_exponents(value::Real)
    exponent_value = _symmetric_binary_exponent(value)
    return exponent_value === nothing ? () : (exponent_value,)
end

function _symmetric_binary_exponents(value::Complex)
    return (
        _symmetric_binary_exponents(real(value))...,
        _symmetric_binary_exponents(imag(value))...,
    )
end

function _symmetric_work_precision(array)
    target = _symmetric_target_precision(array)
    guard = 0
    minimum_exponent = nothing
    maximum_exponent = nothing
    for value in _symmetric_stored_values(array)
        guard = max(guard, _symmetric_precision_guard(value))
        for exponent_value in _symmetric_binary_exponents(value)
            minimum_exponent = if minimum_exponent === nothing
                exponent_value
            else
                min(minimum_exponent, exponent_value)
            end
            maximum_exponent = if maximum_exponent === nothing
                exponent_value
            else
                max(maximum_exponent, exponent_value)
            end
        end
    end
    if minimum_exponent !== nothing
        exponent_span = maximum_exponent - minimum_exponent
        guard = max(guard, 2 * exponent_span + 128)
    end
    return target + guard
end

function _symmetric_wide_value(value::Real)
    return BigFloat(value)
end

function _symmetric_wide_value(value::BigFloat)
    return value
end

function _symmetric_wide_value(value::Complex)
    return Complex{BigFloat}(BigFloat(real(value)), BigFloat(imag(value)))
end

function _symmetric_wide_value(value::Complex{BigFloat})
    return value
end

function _symmetric_next_permutation!(word::Vector{Int})
    pivot = length(word) - 1
    while pivot >= 1 && word[pivot] >= word[pivot + 1]
        pivot -= 1
    end
    pivot == 0 && return false

    successor = length(word)
    while word[successor] <= word[pivot]
        successor -= 1
    end
    word[pivot], word[successor] = word[successor], word[pivot]
    left = pivot + 1
    right = length(word)
    while left < right
        word[left], word[right] = word[right], word[left]
        left += 1
        right -= 1
    end
    return true
end

function _symmetric_typed_sparse_vector(vector::AbstractVector, ::Type{T}) where {T}
    indices = Int[]
    values = T[]
    if issparse(vector)
        canonical = vector isa SparseVector ? vector : sparsevec(vector)
        stored_indices, stored_values = findnz(canonical)
        for index in eachindex(stored_values)
            iszero(stored_values[index]) && continue
            push!(indices, stored_indices[index])
            push!(values, convert(T, stored_values[index]))
        end
    else
        for index in eachindex(vector)
            value = vector[index]
            iszero(value) && continue
            push!(indices, index)
            push!(values, convert(T, value))
        end
    end
    return sparsevec(indices, values, length(vector))
end

function _symmetric_typed_matrix(matrix::AbstractMatrix, ::Type{T}) where {T}
    if issparse(matrix) || _symmetric_is_diagonal_storage(matrix)
        canonical = matrix isa SparseMatrixCSC ? matrix : sparse(matrix)
        rows, columns, values = findnz(canonical)
        converted = T[convert(T, value) for value in values]
        result = sparse(rows, columns, converted, size(matrix, 1), size(matrix, 2))
        dropzeros!(result)
        return result
    end
    return Matrix{T}(matrix)
end

"""
    generalized_dicke_state(
        occupation; normalized=true, sparse_output=true, T=Float64,
        max_nonzeros=1_000_000, max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

Construct the multiqudit Dicke state for occupation `(n₁,...,n_d)`:
`sum |x₁...x_N⟩ / sqrt(N! / prod(n_j!))`, where every computational word has
the requested type. Subsystem 1 is the slowest-varying tensor factor. With
`normalized=false`, every nonzero coefficient is exactly `one(T)`.

The occupation fixes both local dimension and particle number; trailing zeros
are therefore significant. The zero-particle occupation returns the explicit
one-entry vacuum vector. Inputs are never normalized or altered. All ambient,
nonzero, dense-entry, and work costs are checked before enumeration.

# Example
```jldoctest
julia> using SparseArrays

julia> psi = generalized_dicke_state((1, 1, 1));

julia> (length(psi), nnz(psi), isapprox(sum(abs2, psi), 1))
(27, 6, true)
```

The output has ambient length `d^N` and exactly `N!/prod(n_j!)` stored entries;
construction work is `O(N)` per stored entry. This is a project-native
multiqudit generalization, not a QETLAB parity claim.
"""
function generalized_dicke_state(
    occupation;
    normalized::Bool=true,
    sparse_output::Bool=true,
    T::Type{<:Number}=Float64,
    max_nonzeros=_DICKE_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _symmetric_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symmetric_limit(max_work, "max_work")
    occupation isa Tuple ||
        occupation isa AbstractVector ||
        throw(ArgumentError("occupation must be a nonempty tuple or vector of integers"))
    occupation isa AbstractArray && Base.require_one_based_indexing(occupation)
    isempty(occupation) &&
        throw(ArgumentError("occupation must contain at least one level"))
    _symmetric_check_resource(
        BigInt(length(occupation)) + 1,
        work_limit,
        "max_work",
        "occupation-validation operations",
        "generalized_dicke_state",
    )
    checked, parties = _symmetric_checked_occupation(occupation)
    local_dimension = length(checked)

    _symmetric_check_resource(
        BigInt(parties) + BigInt(local_dimension),
        work_limit,
        "max_work",
        "preflight scalar operations",
        "generalized_dicke_state",
    )
    ambient = if local_dimension == 1
        1
    else
        _checked_power(local_dimension, parties, "local_dimension^parties")
    end
    orbit_size = _symmetric_multinomial(checked, parties)
    stored = _symmetric_array_length(orbit_size, "Dicke-state nonzero count")
    _symmetric_check_resource(
        orbit_size,
        nonzero_limit,
        "max_nonzeros",
        "stored entries",
        "generalized_dicke_state",
    )
    if !sparse_output
        _symmetric_check_resource(
            BigInt(ambient),
            dense_limit,
            "max_dense_entries",
            "dense output entries",
            "generalized_dicke_state",
        )
    end
    work =
        BigInt(parties) +
        BigInt(local_dimension) +
        orbit_size * (2 * BigInt(parties) + BigInt(local_dimension) + 2)
    _symmetric_check_resource(
        work,
        work_limit,
        "max_work",
        "estimated scalar operations",
        "generalized_dicke_state",
    )

    indices = sizehint!(Int[], stored)
    if orbit_size == 1
        one_based = 1
        if parties > 0 && local_dimension > 1
            occupied_level = findfirst(==(parties), checked)
            occupied_level === nothing && error("internal single-orbit occupation failure")
            zero_based = 0
            for _ in 1:parties
                zero_based = Base.checked_add(
                    Base.checked_mul(zero_based, local_dimension), occupied_level - 1
                )
            end
            one_based = zero_based + 1
        end
        push!(indices, one_based)
    else
        word = Vector{Int}(undef, parties)
        position = 1
        for level in 1:local_dimension
            for _ in 1:checked[level]
                word[position] = level - 1
                position += 1
            end
        end
        while true
            zero_based = 0
            for label in word
                zero_based = Base.checked_add(
                    Base.checked_mul(zero_based, local_dimension), label
                )
            end
            push!(indices, zero_based + 1)
            _symmetric_next_permutation!(word) || break
        end
    end

    value = if normalized
        orbit_size == 1 ? one(T) : _symmetric_inverse_sqrt(orbit_size, T)
    else
        one(T)
    end
    state = sparsevec(indices, fill(value, stored), ambient)
    return sparse_output ? state : Vector(state)
end

"""
    symmetric_product_coordinates(
        local_state, parties;
        sparse_output=issparse(local_state),
        max_occupations=100_000, max_work=100_000_000,
    )

Return occupation-basis coordinates of `local_state^⊗parties`. For occupation
`n`, the coefficient is `sqrt(N!/prod(n_j!)) * prod(local_state[j]^n_j)`.
There is no conjugation and the supplied local vector is never normalized.
Consequently `norm(result) == norm(local_state)^parties` up to arithmetic
roundoff. `parties=0` returns the one-entry vacuum coordinate.

The output order matches [`symmetric_occupations`](@ref). Floating conversion
uses widened `BigFloat` intermediates so large multinomial factors are not
formed as `Inf/Inf`; the final element type still follows the input precision.

# Example
```jldoctest
julia> c = symmetric_product_coordinates([1 / sqrt(2), 1 / sqrt(2)], 2);

julia> isapprox(c, [0.5, inv(sqrt(2)), 0.5])
true
```

Dense output stores one value per symmetric sector and uses `O(d+N)` work per
sector. Sparse output enumerates only sectors supported by the nonzero local
levels. The function rejects nonfinite input and excessive sector/work
requests before either enumeration.
"""
function symmetric_product_coordinates(
    local_state::AbstractVector{<:Number},
    parties;
    sparse_output::Bool=issparse(local_state),
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    Base.require_one_based_indexing(local_state)
    local_dimension = length(local_state)
    local_dimension > 0 || throw(ArgumentError("local_state must be nonempty"))
    count = _nonnegative_int(parties, "parties")
    count > 1 &&
        _symmetric_require_tuple_arity(local_dimension, "symmetric_product_coordinates")
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    work_limit = _symmetric_limit(max_work, "max_work")
    _symmetric_check_resource(
        BigInt(local_dimension) + BigInt(count) + 1,
        work_limit,
        "max_work",
        "preflight scalar operations",
        "symmetric_product_coordinates",
    )
    sector_exact = _symmetric_dimension(local_dimension, count)
    _symmetric_check_resource(
        sector_exact,
        occupation_limit,
        "max_occupations",
        "occupation sectors",
        "symmetric_product_coordinates",
    )
    sector_count = _symmetric_array_length(sector_exact, "symmetric coordinate count")
    _symmetric_require_finite(local_state, "symmetric_product_coordinates")
    input_type = _symmetric_promoted_eltype(local_state)
    input_type <: Union{Real,Complex} || throw(
        ArgumentError(
            "symmetric_product_coordinates supports real or complex numeric vectors"
        ),
    )
    if count == 0
        return setprecision(BigFloat, _symmetric_work_precision(local_state)) do
            value = one(input_type)
            return sparse_output ? sparsevec([1], [value], 1) : [value]
        end
    elseif count == 1
        _symmetric_check_resource(
            2 * BigInt(local_dimension) + 1,
            work_limit,
            "max_work",
            "estimated endpoint operations",
            "symmetric_product_coordinates",
        )
        return setprecision(BigFloat, _symmetric_work_precision(local_state)) do
            if sparse_output
                return _symmetric_typed_sparse_vector(local_state, input_type)
            end
            return Vector{input_type}(local_state)
        end
    end

    active_dimension = if issparse(local_state)
        Base.count(!iszero, _symmetric_stored_values(local_state))
    else
        Base.count(!iszero, local_state)
    end
    active_exact = if active_dimension == 0
        BigInt(0)
    else
        _symmetric_dimension(active_dimension, count)
    end
    enumerated = sparse_output ? active_exact : sector_exact
    work =
        BigInt(local_dimension) +
        enumerated *
        (BigInt(local_dimension) + BigInt(active_dimension) + BigInt(count) + 3)
    _symmetric_check_resource(
        work,
        work_limit,
        "max_work",
        "estimated scalar operations",
        "symmetric_product_coordinates",
    )

    active_levels = if issparse(local_state)
        canonical = local_state isa SparseVector ? local_state : sparsevec(local_state)
        stored_indices, stored_values = findnz(canonical)
        [
            stored_indices[index] for
            index in eachindex(stored_values) if !iszero(stored_values[index])
        ]
    else
        findall(!iszero, local_state)
    end
    length(active_levels) == active_dimension ||
        error("internal product-coordinate support-count mismatch")

    output_type = _symmetric_root_type(input_type)
    function coordinate(occupation)
        multinomial = _symmetric_multinomial(occupation, count)
        wide_product = one(_symmetric_wide_value(first(local_state)))
        for level in 1:local_dimension
            wide_product *= _symmetric_wide_value(local_state[level])^occupation[level]
        end
        wide_value = sqrt(BigFloat(multinomial)) * wide_product
        value = convert(output_type, wide_value)
        (isfinite(value) && (iszero(wide_value) || !iszero(value))) || throw(
            OverflowError(
                "a nonzero product-state coordinate is not representable as a finite " *
                "$output_type; use BigFloat or Complex{BigFloat} input",
            ),
        )
        return value
    end

    return setprecision(BigFloat, _symmetric_work_precision(local_state)) do
        if !sparse_output
            values = Vector{output_type}(undef, sector_count)
            index = 0
            _symmetric_foreach_occupation(local_dimension, count) do occupation
                index += 1
                return values[index] = coordinate(occupation)
            end
            return values
        elseif active_dimension == 0
            return spzeros(output_type, sector_count)
        end

        active_count = _symmetric_array_length(
            active_exact, "supported symmetric coordinate count"
        )
        indices = Vector{Int}(undef, active_count)
        values = Vector{output_type}(undef, active_count)
        full_occupation = zeros(Int, local_dimension)
        index = 0
        _symmetric_foreach_occupation(active_dimension, count) do active_occupation
            index += 1
            fill!(full_occupation, 0)
            for active_index in 1:active_dimension
                full_occupation[active_levels[active_index]] = active_occupation[active_index]
            end
            occupation = Tuple(full_occupation)
            indices[index] = _symmetric_array_length(
                _symmetric_basis_index_checked(occupation, count), "symmetric basis index"
            )
            return values[index] = coordinate(occupation)
        end
        return sparsevec(indices, values, sector_count)
    end
end

function _symmetric_matrix_entries(matrix::AbstractMatrix, ::Type{T}) where {T}
    rows = Int[]
    columns = Int[]
    values = T[]
    if issparse(matrix) || !(matrix isa StridedMatrix)
        canonical = matrix isa SparseMatrixCSC ? matrix : sparse(matrix)
        sparse_rows, sparse_columns, sparse_values = findnz(canonical)
        for index in eachindex(sparse_values)
            iszero(sparse_values[index]) && continue
            push!(rows, sparse_rows[index])
            push!(columns, sparse_columns[index])
            push!(values, convert(T, sparse_values[index]))
        end
    else
        for column in axes(matrix, 2), row in axes(matrix, 1)
            value = matrix[row, column]
            iszero(value) && continue
            push!(rows, row)
            push!(columns, column)
            push!(values, convert(T, value))
        end
    end
    return rows, columns, values
end

function _symmetric_matrix_profile(matrix::AbstractMatrix)
    stored = 0
    has_offdiagonal = false
    if _symmetric_is_diagonal_storage(matrix)
        for value in diag(matrix)
            iszero(value) || (stored += 1)
        end
    elseif matrix isa SparseMatrixCSC
        values = nonzeros(matrix)
        row_indices = rowvals(matrix)
        for column in axes(matrix, 2), pointer in nzrange(matrix, column)
            iszero(values[pointer]) && continue
            stored += 1
            has_offdiagonal |= row_indices[pointer] != column
        end
    elseif issparse(matrix)
        return _symmetric_matrix_profile(sparse(matrix))
    else
        for column in axes(matrix, 2), row in axes(matrix, 1)
            iszero(matrix[row, column]) && continue
            stored += 1
            has_offdiagonal |= row != column
        end
    end
    return BigInt(stored), has_offdiagonal
end

function _symmetric_checked_scale(value, factor::Int, operation::AbstractString)
    result = value * factor
    isfinite(result) || throw(
        OverflowError(
            "$operation produced a nonfinite entry; use a wider numeric element type"
        ),
    )
    return result
end

function _symmetric_checked_scale(
    value::T, factor::Int, operation::AbstractString
) where {T<:Base.BitInteger}
    converted_factor = try
        convert(T, factor)
    catch err
        err isa InexactError || rethrow()
        throw(OverflowError("$operation exceeds the representable range of $T"))
    end
    return try
        Base.checked_mul(value, converted_factor)
    catch err
        err isa OverflowError || rethrow()
        throw(OverflowError("$operation exceeds the representable range of $T"))
    end
end

function _symmetric_checked_add(left, right, operation::AbstractString)
    result = left + right
    isfinite(result) || throw(
        OverflowError(
            "$operation produced a nonfinite entry; use a wider numeric element type"
        ),
    )
    return result
end

function _symmetric_checked_add(
    left::T, right::T, operation::AbstractString
) where {T<:Base.BitInteger}
    return try
        Base.checked_add(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(OverflowError("$operation exceeds the representable range of $T"))
    end
end

function _symmetric_accumulation_type(::Type{T}) where {T}
    return T
end

function _symmetric_accumulation_type(::Type{T}) where {T<:Base.BitInteger}
    return BigInt
end

function _symmetric_accumulation_type(::Type{Rational{T}}) where {T<:Base.BitInteger}
    return Rational{BigInt}
end

function _symmetric_accumulation_type(::Type{T}) where {T<:AbstractFloat}
    return T <: BigFloat ? BigFloat : BigFloat
end

function _symmetric_accumulation_type(::Type{Complex{T}}) where {T}
    return Complex{_symmetric_accumulation_type(T)}
end

function _symmetric_narrow_collective_value(value, ::Type{T}) where {T}
    converted = try
        convert(T, value)
    catch err
        err isa InexactError || err isa OverflowError || rethrow()
        throw(
            OverflowError(
                "a collective-operator entry is not representable as $T; " *
                "use a wider numeric type",
            ),
        )
    end
    if !isfinite(converted) || (!iszero(value) && iszero(converted))
        throw(
            OverflowError(
                "a nonzero collective-operator entry is not representable as a finite " *
                "$T; use a wider numeric type",
            ),
        )
    end
    return converted
end

"""
    symmetric_collective_operator(
        local_operator, parties;
        sparse_output=true, max_occupations=100_000,
        max_nonzeros=5_000_000, max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

Represent `sum_{j=1}^N local_operator^(j)` directly in symmetric occupation
coordinates. A local matrix element `A[a,b]` moves one particle from level
`b` to `a` with bosonic factor `sqrt(n_b*(n_a+1))`; diagonal elements act as
`A[a,a]*n_a`. The input is never Hermitized, normalized, or densified.

The result has `D × D` shape for `D=binomial(N+d-1,N)`. Sparse output is built
without constructing the ambient `d^N` operator or basis. `parties=0` returns
the explicit `1 × 1` zero operator.

# Example
```jldoctest
julia> Z = [1 0; 0 -1];

julia> Matrix(symmetric_collective_operator(Z, 2))
3×3 Matrix{Int64}:
 2  0   0
 0  0   0
 0  0  -2
```

For `q` stored local entries, candidate work is `O(D*q*d)` and candidate
storage is at most `D*q`, both checked with exact `BigInt` estimates.
"""
function symmetric_collective_operator(
    local_operator::AbstractMatrix{<:Number},
    parties;
    sparse_output::Bool=true,
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_nonzeros=_SYMMETRIC_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    Base.require_one_based_indexing(local_operator)
    size(local_operator, 1) == size(local_operator, 2) || throw(
        DimensionMismatch(
            "local_operator must be square; got size $(size(local_operator))"
        ),
    )
    local_dimension = size(local_operator, 1)
    local_dimension > 0 || throw(ArgumentError("local_operator must be nonempty"))
    count = _nonnegative_int(parties, "parties")
    count > 1 &&
        _symmetric_require_tuple_arity(local_dimension, "symmetric_collective_operator")
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _symmetric_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symmetric_limit(max_work, "max_work")

    _symmetric_check_resource(
        BigInt(local_dimension) + BigInt(count) + 1,
        work_limit,
        "max_work",
        "preflight scalar operations",
        "symmetric_collective_operator",
    )
    output_exact = _symmetric_dimension(local_dimension, count)
    _symmetric_check_resource(
        output_exact,
        occupation_limit,
        "max_occupations",
        "occupation sectors",
        "symmetric_collective_operator",
    )
    output_dimension = _symmetric_array_length(
        output_exact, "collective-operator dimension"
    )
    input_entries = _symmetric_input_entries(local_operator)
    _symmetric_check_resource(
        input_entries + BigInt(local_dimension) + 1,
        work_limit,
        "max_work",
        "input-validation operations",
        "symmetric_collective_operator",
    )
    input_type = _symmetric_promoted_eltype(local_operator)
    _symmetric_require_finite(local_operator, "symmetric_collective_operator")
    if count == 0
        if !sparse_output
            _symmetric_check_resource(
                BigInt(1),
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_collective_operator",
            )
        end
        return setprecision(BigFloat, _symmetric_work_precision(local_operator)) do
            zero_operator = spzeros(input_type, 1, 1)
            return sparse_output ? zero_operator : Matrix(zero_operator)
        end
    end

    stored_entries, has_offdiagonal = _symmetric_matrix_profile(local_operator)
    if count == 1
        _symmetric_check_resource(
            stored_entries,
            nonzero_limit,
            "max_nonzeros",
            "stored entries",
            "symmetric_collective_operator",
        )
        if !sparse_output
            _symmetric_check_resource(
                output_exact^2,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_collective_operator",
            )
        end
        _symmetric_check_resource(
            input_entries + stored_entries + BigInt(local_dimension),
            work_limit,
            "max_work",
            "estimated endpoint operations",
            "symmetric_collective_operator",
        )
        return setprecision(BigFloat, _symmetric_work_precision(local_operator)) do
            local_rows, local_columns, local_values = _symmetric_matrix_entries(
                local_operator, input_type
            )
            exact_operator = sparse(
                local_rows, local_columns, local_values, local_dimension, local_dimension
            )
            return sparse_output ? exact_operator : Matrix(exact_operator)
        end
    end

    if iszero(stored_entries)
        if !sparse_output
            _symmetric_check_resource(
                output_exact^2,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_collective_operator",
            )
        end
        _symmetric_check_resource(
            input_entries + output_exact,
            work_limit,
            "max_work",
            "zero-output construction operations",
            "symmetric_collective_operator",
        )
        return setprecision(BigFloat, _symmetric_work_precision(local_operator)) do
            zero_operator = spzeros(input_type, output_dimension, output_dimension)
            return sparse_output ? zero_operator : Matrix(zero_operator)
        end
    end

    candidate_work = output_exact * stored_entries
    candidate_entries = if has_offdiagonal
        candidate_work
    else
        output_exact
    end
    _symmetric_check_resource(
        candidate_entries,
        nonzero_limit,
        "max_nonzeros",
        "candidate stored entries",
        "symmetric_collective_operator",
    )
    if !sparse_output
        _symmetric_check_resource(
            output_exact^2,
            dense_limit,
            "max_dense_entries",
            "dense output entries",
            "symmetric_collective_operator",
        )
    end
    work =
        input_entries +
        output_exact * (BigInt(local_dimension) + 1) +
        candidate_work * (BigInt(local_dimension) + 3)
    _symmetric_check_resource(
        work,
        work_limit,
        "max_work",
        "estimated scalar operations",
        "symmetric_collective_operator",
    )

    output_type = if has_offdiagonal
        _symmetric_root_type(input_type)
    elseif input_type <: Bool
        Int
    elseif input_type <: Complex{<:Base.BitInteger}
        Complex{BigInt}
    else
        input_type
    end
    return setprecision(BigFloat, _symmetric_work_precision(local_operator)) do
        local_rows, local_columns, local_values = _symmetric_matrix_entries(
            local_operator, input_type
        )
        BigInt(length(local_values)) == stored_entries ||
            error("internal collective-operator nonzero-count mismatch")
        converted_local_values = convert.(output_type, local_values)

        accumulation_type = _symmetric_accumulation_type(output_type)
        diagonal = Vector{output_type}(undef, output_dimension)
        column = 0
        _symmetric_foreach_occupation(local_dimension, count) do occupation
            column += 1
            value = zero(accumulation_type)
            for entry in eachindex(converted_local_values)
                local_rows[entry] == local_columns[entry] || continue
                factor = occupation[local_columns[entry]]
                factor == 0 && continue
                term = _symmetric_checked_scale(
                    convert(accumulation_type, local_values[entry]),
                    factor,
                    "symmetric_collective_operator",
                )
                value = _symmetric_checked_add(value, term, "symmetric_collective_operator")
            end
            return diagonal[column] = _symmetric_narrow_collective_value(value, output_type)
        end
        diagonal_operator = spdiagm(0 => diagonal)
        dropzeros!(diagonal_operator)
        if !has_offdiagonal
            return sparse_output ? diagonal_operator : Matrix(diagonal_operator)
        end

        rows = Int[]
        columns = Int[]
        values = output_type[]
        hint = min(
            _symmetric_array_length(candidate_entries, "candidate entries"), 1_000_000
        )
        sizehint!(rows, hint)
        sizehint!(columns, hint)
        sizehint!(values, hint)
        column = 0
        _symmetric_foreach_occupation(local_dimension, count) do occupation
            column += 1
            for entry in eachindex(converted_local_values)
                target = local_rows[entry]
                source = local_columns[entry]
                target == source && continue
                source_count = occupation[source]
                source_count == 0 && continue
                target_occupation = collect(occupation)
                target_occupation[source] -= 1
                target_occupation[target] += 1
                row = _symmetric_array_length(
                    symmetric_basis_index(target_occupation), "symmetric basis index"
                )
                factor = _symmetric_ratio_sqrt(
                    BigInt(source_count) * (BigInt(occupation[target]) + 1),
                    BigInt(1),
                    output_type,
                )
                value = converted_local_values[entry] * factor
                (
                    isfinite(value) &&
                    (iszero(converted_local_values[entry]) || !iszero(value))
                ) || throw(
                    OverflowError(
                        "a nonzero collective-operator entry is not representable " *
                        "as a finite $output_type; use a wider numeric type",
                    ),
                )
                push!(rows, row)
                push!(columns, column)
                push!(values, value)
            end
        end
        offdiagonal_operator = sparse(
            rows, columns, values, output_dimension, output_dimension
        )
        operator = diagonal_operator + offdiagonal_operator
        all(isfinite, nonzeros(operator)) || throw(
            OverflowError(
                "assembled collective-operator entries are not representable as " *
                "$output_type; use a wider numeric type",
            ),
        )
        return sparse_output ? operator : Matrix(operator)
    end
end

function _symmetric_split_counts(local_dimension::Int, parties::Int, keep::Int)
    kept_dimension = _symmetric_dimension(local_dimension, keep)
    traced_dimension = _symmetric_dimension(local_dimension, parties - keep)
    global_dimension = _symmetric_dimension(local_dimension, parties)
    pair_count = kept_dimension * traced_dimension
    return kept_dimension, traced_dimension, global_dimension, pair_count
end

function _symmetric_split_preflight(
    local_dimension::Int,
    parties::Int,
    keep::Int,
    operation::AbstractString;
    max_occupations,
    max_nonzeros,
    max_work,
)
    _symmetric_require_tuple_arity(local_dimension, operation)
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    work_limit = _symmetric_limit(max_work, "max_work")
    _symmetric_check_resource(
        BigInt(local_dimension) + BigInt(parties) + BigInt(keep) + 1,
        work_limit,
        "max_work",
        "preflight scalar operations",
        operation,
    )
    kept_exact, traced_exact, global_exact, pair_count = _symmetric_split_counts(
        local_dimension, parties, keep
    )
    for (planned, label) in (
        (kept_exact, "kept occupation sectors"),
        (traced_exact, "traced occupation sectors"),
        (global_exact, "global occupation sectors"),
    )
        _symmetric_check_resource(
            planned, occupation_limit, "max_occupations", label, operation
        )
    end
    _symmetric_check_resource(
        pair_count, nonzero_limit, "max_nonzeros", "split coefficients", operation
    )
    work =
        pair_count * (3 * BigInt(local_dimension) + 4) +
        global_exact * (BigInt(local_dimension) + 1)
    _symmetric_check_resource(
        work, work_limit, "max_work", "estimated split-construction operations", operation
    )

    kept_dimension = _symmetric_array_length(kept_exact, "kept symmetric dimension")
    traced_dimension = _symmetric_array_length(traced_exact, "traced symmetric dimension")
    global_dimension = _symmetric_array_length(global_exact, "global symmetric dimension")
    _symmetric_array_length(pair_count, "symmetric split row count")
    return (
        kept_exact=kept_exact,
        traced_exact=traced_exact,
        global_exact=global_exact,
        pair_count=pair_count,
        kept_dimension=kept_dimension,
        traced_dimension=traced_dimension,
        global_dimension=global_dimension,
        work=work,
    )
end

function _symmetric_split_data(
    local_dimension::Int, parties::Int, keep::Int, ::Type{T}, plan
) where {T<:Number}
    kept_occupations = _symmetric_occupations_unchecked(local_dimension, keep)
    traced_occupations = _symmetric_occupations_unchecked(local_dimension, parties - keep)
    global_occupations = _symmetric_occupations_unchecked(local_dimension, parties)
    global_multinomials = [
        _symmetric_multinomial(occupation, parties) for occupation in global_occupations
    ]
    kept_multinomials = [
        _symmetric_multinomial(occupation, keep) for occupation in kept_occupations
    ]
    traced_multinomials = [
        _symmetric_multinomial(occupation, parties - keep) for
        occupation in traced_occupations
    ]

    root_type = _symmetric_root_type(T)
    indices = Matrix{Int}(undef, plan.kept_dimension, plan.traced_dimension)
    weights = Matrix{root_type}(undef, plan.kept_dimension, plan.traced_dimension)
    for traced_index in 1:plan.traced_dimension, kept_index in 1:plan.kept_dimension
        global_occupation = ntuple(
            level ->
                kept_occupations[kept_index][level] +
                traced_occupations[traced_index][level],
            local_dimension,
        )
        global_index = _symmetric_array_length(
            _symmetric_basis_index_checked(global_occupation, parties),
            "global symmetric basis index",
        )
        indices[kept_index, traced_index] = global_index
        numerator = kept_multinomials[kept_index] * traced_multinomials[traced_index]
        weights[kept_index, traced_index] = _symmetric_ratio_sqrt(
            numerator, global_multinomials[global_index], T
        )
    end
    return (
        kept_dimension=plan.kept_dimension,
        traced_dimension=plan.traced_dimension,
        global_dimension=plan.global_dimension,
        indices=indices,
        weights=weights,
    )
end

"""
    symmetric_split_isometry(
        local_dimension, parties, keep;
        T=Float64, sparse_output=true, max_occupations=100_000,
        max_nonzeros=5_000_000, max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

Return the isometry from global `N`-particle occupation coordinates into
`Sym^keep(C^d) ⊗ Sym^(N-keep)(C^d)`. The kept factor is the slow (left) tensor
factor. For occupations `a+c=n`, its nonzero coefficient is
`sqrt(M_a*M_c/M_n)`, where `M_n=N!/prod(n_j!)`.

The matrix has `D_keep*D_trace` rows, `D_N` columns, exactly one stored entry
per row, and satisfies `V'V ≈ I` up to the precision of `T`. `keep=0` and
`keep=N` are exact identity embeddings and preserve exact `T` values. Costs are
preflighted with exact combinatorics and nontrivial ratio square roots use
widened intermediates.

# Example
```jldoctest
julia> using LinearAlgebra

julia> V = symmetric_split_isometry(2, 3, 1);

julia> isapprox(Matrix(V' * V), Matrix{Float64}(I, 4, 4))
true
```

Sparse construction uses `O(D_keep*D_trace*d)` work and exactly
`D_keep*D_trace` nonzeros; dense output additionally requires explicit room
under `max_dense_entries`.
"""
function symmetric_split_isometry(
    local_dimension,
    parties,
    keep;
    T::Type{<:Number}=Float64,
    sparse_output::Bool=true,
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_nonzeros=_SYMMETRIC_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    kept = _nonnegative_int(keep, "keep")
    kept <= count || throw(ArgumentError("keep=$kept must not exceed parties=$count"))
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _symmetric_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symmetric_limit(max_work, "max_work")
    if kept == 0 || kept == count
        preliminary_work = BigInt(dimension) + BigInt(count) + BigInt(kept) + 1
        _symmetric_check_resource(
            preliminary_work,
            work_limit,
            "max_work",
            "preflight scalar operations",
            "symmetric_split_isometry",
        )
        global_exact = _symmetric_dimension(dimension, count)
        _symmetric_check_resource(
            global_exact,
            occupation_limit,
            "max_occupations",
            "global occupation sectors",
            "symmetric_split_isometry",
        )
        _symmetric_check_resource(
            global_exact,
            nonzero_limit,
            "max_nonzeros",
            "identity entries",
            "symmetric_split_isometry",
        )
        _symmetric_check_resource(
            preliminary_work + global_exact,
            work_limit,
            "max_work",
            "estimated identity-construction operations",
            "symmetric_split_isometry",
        )
        if !sparse_output
            _symmetric_check_resource(
                global_exact^2,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_split_isometry",
            )
        end
        global_dimension = _symmetric_array_length(
            global_exact, "global symmetric dimension"
        )
        identity = spdiagm(0 => fill(one(T), global_dimension))
        return sparse_output ? identity : Matrix(identity)
    end

    plan = _symmetric_split_preflight(
        dimension,
        count,
        kept,
        "symmetric_split_isometry";
        max_occupations=max_occupations,
        max_nonzeros=max_nonzeros,
        max_work=max_work,
    )
    dense_entries = plan.pair_count * plan.global_exact
    if !sparse_output
        _symmetric_check_resource(
            dense_entries,
            dense_limit,
            "max_dense_entries",
            "dense output entries",
            "symmetric_split_isometry",
        )
    end

    row_count = _symmetric_array_length(plan.pair_count, "symmetric split row count")
    split = _symmetric_split_data(dimension, count, kept, T, plan)
    rows = collect(1:row_count)
    columns = Vector{Int}(undef, row_count)
    values = Vector{eltype(split.weights)}(undef, row_count)
    for kept_index in 1:split.kept_dimension, traced_index in 1:split.traced_dimension
        row = (kept_index - 1) * split.traced_dimension + traced_index
        columns[row] = split.indices[kept_index, traced_index]
        values[row] = split.weights[kept_index, traced_index]
    end
    isometry = sparse(rows, columns, values, row_count, plan.global_dimension)
    return sparse_output ? isometry : Matrix(isometry)
end

function _symmetric_sparse_reduction(state::AbstractVector, split, ::Type{T}) where {T}
    accumulated = Dict{Tuple{Int,Int},T}()
    for traced_index in 1:split.traced_dimension
        active = Tuple{Int,T}[]
        for kept_index in 1:split.kept_dimension
            amplitude =
                convert(T, split.weights[kept_index, traced_index]) *
                convert(T, state[split.indices[kept_index, traced_index]])
            iszero(amplitude) || push!(active, (kept_index, amplitude))
        end
        for (row, row_value) in active, (column, column_value) in active
            key = (row, column)
            accumulated[key] =
                get(accumulated, key, zero(T)) + row_value * conj(column_value)
        end
    end
    rows = Int[]
    columns = Int[]
    values = T[]
    for ((row, column), value) in accumulated
        iszero(value) && continue
        push!(rows, row)
        push!(columns, column)
        push!(values, value)
    end
    return sparse(rows, columns, values, split.kept_dimension, split.kept_dimension)
end

function _symmetric_sparse_reduction(state::AbstractMatrix, split, ::Type{T}) where {T}
    accumulated = Dict{Tuple{Int,Int},T}()
    for traced_index in 1:split.traced_dimension,
        column in 1:split.kept_dimension,
        row in 1:split.kept_dimension

        value = state[split.indices[row, traced_index], split.indices[column, traced_index]]
        iszero(value) && continue
        contribution =
            convert(T, split.weights[row, traced_index]) *
            convert(T, value) *
            convert(T, split.weights[column, traced_index])
        key = (row, column)
        accumulated[key] = get(accumulated, key, zero(T)) + contribution
    end
    rows = Int[]
    columns = Int[]
    values = T[]
    for ((row, column), value) in accumulated
        iszero(value) && continue
        push!(rows, row)
        push!(columns, column)
        push!(values, value)
    end
    return sparse(rows, columns, values, split.kept_dimension, split.kept_dimension)
end

function _symmetric_wide_numeric_type(::Type{T}) where {T<:Real}
    return BigFloat
end

function _symmetric_wide_numeric_type(::Type{Complex{T}}) where {T}
    return Complex{BigFloat}
end

function _symmetric_reduction_weight_plan(
    state, input_type::Type, local_dimension::Int, parties::Int, traced_dimension::Int
)
    output_type = _symmetric_root_type(input_type)
    real_type = _symmetric_real_type(output_type)
    input_real_type = _symmetric_real_type(input_type)
    input_real_type <: AbstractFloat || return BigFloat, true
    real_type <: AbstractFloat || return BigFloat, true
    real_type <: BigFloat && return BigFloat, false

    minimum_state = nothing
    maximum_state = BigFloat(0)
    for value in _symmetric_stored_values(state)
        iszero(value) && continue
        magnitude = abs(_symmetric_wide_value(value))
        minimum_state =
            minimum_state === nothing ? magnitude : min(minimum_state, magnitude)
        maximum_state = max(maximum_state, magnitude)
    end
    minimum_state === nothing && return real_type, false

    term_count = BigFloat(max(traced_dimension, 1))
    smallest_positive = BigFloat(nextfloat(zero(real_type)))
    largest_finite = BigFloat(floatmax(real_type))
    roundoff = BigFloat(eps(real_type))
    minimum_weight_squared = exp2(-BigFloat(parties) * log2(BigFloat(local_dimension)))
    minimum_contribution = if state isa AbstractVector
        minimum_state^2 * minimum_weight_squared
    else
        minimum_state * minimum_weight_squared
    end
    maximum_contribution = state isa AbstractVector ? maximum_state^2 : maximum_state
    underflow_risk = minimum_contribution < smallest_positive
    overflow_risk = maximum_contribution * term_count > largest_finite
    cancellation_risk = minimum_state < 64 * roundoff * maximum_state * term_count
    use_wide = underflow_risk || overflow_risk || cancellation_risk
    return use_wide ? BigFloat : real_type, use_wide
end

function _symmetric_needs_wide_reduction(state, split, output_type::Type)
    real_type = _symmetric_real_type(output_type)
    real_type <: AbstractFloat || return false
    real_type <: BigFloat && return false

    minimum_contribution, maximum_contribution = if state isa AbstractVector
        minimum_amplitude = nothing
        maximum_amplitude = BigFloat(0)
        for traced_index in 1:split.traced_dimension, kept_index in 1:split.kept_dimension

            value = state[split.indices[kept_index, traced_index]]
            iszero(value) && continue
            amplitude =
                abs(_symmetric_wide_value(value)) *
                BigFloat(abs(split.weights[kept_index, traced_index]))
            minimum_amplitude = if minimum_amplitude === nothing
                amplitude
            else
                min(minimum_amplitude, amplitude)
            end
            maximum_amplitude = max(maximum_amplitude, amplitude)
        end
        minimum_amplitude === nothing && return false
        (minimum_amplitude^2, maximum_amplitude^2)
    else
        minimum_state = nothing
        maximum_state = BigFloat(0)
        for value in _symmetric_stored_values(state)
            iszero(value) && continue
            magnitude = abs(_symmetric_wide_value(value))
            minimum_state =
                minimum_state === nothing ? magnitude : min(minimum_state, magnitude)
            maximum_state = max(maximum_state, magnitude)
        end
        minimum_state === nothing && return false
        minimum_weight = minimum(
            BigFloat(abs(value)) for value in split.weights if !iszero(value)
        )
        maximum_weight = maximum(
            BigFloat(abs(value)) for value in split.weights if !iszero(value)
        )
        (minimum_state * minimum_weight^2, maximum_state * maximum_weight^2)
    end
    term_count = BigFloat(max(split.traced_dimension, 1))
    smallest_positive = BigFloat(nextfloat(zero(real_type)))
    largest_finite = BigFloat(floatmax(real_type))
    roundoff = BigFloat(eps(real_type))
    underflow_risk = minimum_contribution < smallest_positive
    overflow_risk = maximum_contribution * term_count > largest_finite
    cancellation_risk =
        minimum_contribution < 64 * roundoff * maximum_contribution * term_count
    return underflow_risk || overflow_risk || cancellation_risk
end

function _symmetric_narrow_value(value, ::Type{T}) where {T}
    converted = try
        convert(T, value)
    catch err
        err isa InexactError || err isa OverflowError || rethrow()
        throw(
            OverflowError(
                "a reduced-state entry is not representable as $T; " *
                "use a wider numeric element type",
            ),
        )
    end
    if !isfinite(converted) || (!iszero(value) && iszero(converted))
        throw(
            OverflowError(
                "a nonzero reduced-state entry is not representable as a finite $T; " *
                "use a wider numeric element type",
            ),
        )
    end
    return converted
end

function _symmetric_narrow_matrix(matrix::AbstractMatrix, ::Type{T}) where {T}
    if issparse(matrix)
        rows, columns, values = findnz(matrix isa SparseMatrixCSC ? matrix : sparse(matrix))
        narrowed = T[_symmetric_narrow_value(value, T) for value in values]
        result = sparse(rows, columns, narrowed, size(matrix, 1), size(matrix, 2))
        dropzeros!(result)
        return result
    end
    return map(value -> _symmetric_narrow_value(value, T), matrix)
end

function _symmetric_exact_quotient_type(::Type{T}) where {T}
    return T <: Rational
end

function _symmetric_exact_quotient_type(::Type{Complex{T}}) where {T}
    return T <: Rational
end

function _symmetric_reciprocal(exact::BigInt, ::Type{T}) where {T<:Number}
    exact > 0 || throw(ArgumentError("the reciprocal denominator must be positive"))
    output_type = typeof(one(T) / one(T))
    value = try
        if _symmetric_exact_quotient_type(output_type)
            convert(output_type, BigInt(1) // exact)
        else
            convert(output_type, inv(BigFloat(exact)))
        end
    catch err
        err isa InexactError || err isa OverflowError || rethrow()
        throw(
            OverflowError(
                "1/$exact is not representable as $output_type; use a wider numeric type",
            ),
        )
    end
    if !isfinite(value) || iszero(value)
        throw(
            OverflowError(
                "1/$exact is not representable as a finite nonzero $output_type; " *
                "use a wider numeric type",
            ),
        )
    end
    return value
end

function _symmetric_endpoint_type(::Type{T}) where {T}
    return T
end

function _symmetric_endpoint_type(::Type{T}) where {T<:Integer}
    return BigInt
end

function _symmetric_endpoint_type(::Type{Complex{T}}) where {T<:Integer}
    return Complex{BigInt}
end

function _symmetric_widen_endpoint_type(::Type{T}) where {T}
    return T
end

function _symmetric_widen_endpoint_type(::Type{Rational{T}}) where {T<:Base.BitInteger}
    return Rational{BigInt}
end

function _symmetric_widen_endpoint_type(::Type{Complex{T}}) where {T}
    widened = _symmetric_widen_endpoint_type(T)
    return widened === T ? Complex{T} : Complex{widened}
end

function _symmetric_validate_finite_output(result, operation::AbstractString)
    values = issparse(result) ? nonzeros(result) : result
    all(isfinite, values) || throw(
        OverflowError(
            "$operation produced a nonfinite output; use a wider numeric element type"
        ),
    )
    return result
end

function _symmetric_real_type(::Type{T}) where {T}
    return T
end

function _symmetric_real_type(::Type{Complex{T}}) where {T}
    return T
end

function _symmetric_validate_reduction_trace(
    state, result, input_type::Type, target_precision::Int
)
    input_real_type = _symmetric_real_type(input_type)
    result_real_type = _symmetric_real_type(eltype(result))
    real_type = if input_real_type <: AbstractFloat
        input_real_type
    elseif result_real_type <: AbstractFloat
        result_real_type
    else
        return result
    end

    expected = if state isa AbstractVector
        total = BigFloat(0)
        for value in _symmetric_stored_values(state)
            total += abs2(_symmetric_wide_value(value))
        end
        total
    else
        total = zero(_symmetric_wide_value(state[firstindex(state, 1), firstindex(state, 2)]))
        for index in axes(state, 1)
            total += _symmetric_wide_value(state[index, index])
        end
        total
    end
    actual = zero(_symmetric_wide_value(result[1, 1]))
    for index in axes(result, 1)
        actual += _symmetric_wide_value(result[index, index])
    end
    scale = if state isa AbstractVector
        abs(expected)
    else
        sum(abs(_symmetric_wide_value(state[index, index])) for index in axes(state, 1))
    end
    unit_roundoff = if real_type <: BigFloat
        exp2(BigFloat(1 - target_precision))
    else
        BigFloat(eps(real_type))
    end
    tolerance_scale = if input_real_type <: AbstractFloat
        max(scale, abs(expected), abs(actual))
    else
        max(BigFloat(1), abs(expected), abs(actual))
    end
    tolerance = 64 * unit_roundoff * tolerance_scale
    abs(actual - expected) <= tolerance || throw(
        OverflowError(
            "symmetric_reduced_state cannot preserve the trace at the input " *
            "precision; use a wider numeric element type",
        ),
    )
    return result
end

function _symmetric_finish_reduction(state, result, input_type::Type, target_precision::Int)
    _symmetric_validate_finite_output(result, "symmetric_reduced_state")
    return _symmetric_validate_reduction_trace(state, result, input_type, target_precision)
end

function _symmetric_reduction_endpoint(
    state,
    input_type::Type,
    kept::Int,
    count::Int,
    sparse_output::Bool,
    target_precision::Int,
)
    try
        return setprecision(BigFloat, _symmetric_work_precision(state)) do
            safe_type = _symmetric_endpoint_type(input_type)
            if kept == 0
                accumulation_type = _symmetric_accumulation_type(safe_type)
                wide_value = if state isa AbstractVector
                    total = zero(typeof(abs2(one(accumulation_type))))
                    for entry in _symmetric_stored_values(state)
                        total += abs2(convert(accumulation_type, entry))
                    end
                    total
                else
                    total = zero(accumulation_type)
                    for index in axes(state, 1)
                        total += convert(accumulation_type, state[index, index])
                    end
                    total
                end
                target_type =
                    state isa AbstractVector ? typeof(abs2(one(safe_type))) : safe_type
                value = _symmetric_narrow_value(wide_value, target_type)
                isfinite(value) || throw(
                    OverflowError(
                        "symmetric_reduced_state produced a nonfinite scalar; " *
                        "use a wider numeric element type",
                    ),
                )
                result = if iszero(value)
                    spzeros(typeof(value), 1, 1)
                else
                    sparse([1], [1], [value], 1, 1)
                end
                finished = sparse_output ? result : Matrix(result)
                return _symmetric_finish_reduction(
                    state, finished, input_type, target_precision
                )
            end

            kept == count || error("internal symmetric-reduction endpoint failure")
            result = if state isa AbstractVector
                accumulation_type = _symmetric_accumulation_type(safe_type)
                safe_state = if issparse(state)
                    _symmetric_typed_sparse_vector(state, accumulation_type)
                else
                    Vector{accumulation_type}(state)
                end
                wide_result = safe_state * adjoint(safe_state)
                if accumulation_type === safe_type
                    wide_result
                else
                    _symmetric_narrow_matrix(wide_result, safe_type)
                end
            else
                _symmetric_typed_matrix(state, input_type)
            end
            finished = sparse_output ? sparse(result) : Matrix(result)
            return _symmetric_finish_reduction(
                state, finished, input_type, target_precision
            )
        end
    catch err
        err isa OverflowError || rethrow()
        widened = _symmetric_widen_endpoint_type(input_type)
        widened === input_type && rethrow()
        return _symmetric_reduction_endpoint(
            state, widened, kept, count, sparse_output, target_precision
        )
    end
end

"""
    symmetric_reduced_state(
        state, local_dimension, parties; keep,
        sparse_output=issparse(state), max_occupations=100_000,
        max_nonzeros=5_000_000, max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

Trace `parties-keep` systems directly in occupation coordinates. A vector input
is interpreted as a pure state and returns the reduction of `state*state'`; a
matrix input is treated linearly as an operator. Inputs are never normalized,
Hermitized, projected, or repaired.

The input dimension is `D_N`; the output is a `D_keep × D_keep` operator in
the same ordered occupation convention. `keep=0` returns an explicit `1 × 1`
norm-squared or trace, and `keep=parties` returns the pure projector or an
equivalent copy of the operator. These endpoint branches preserve exact input
arithmetic. Sparse inputs remain sparse by default.

# Example
```jldoctest
julia> psi = generalized_dicke_state((2, 1));

julia> coordinates = symmetric_subspace_basis(2, 3)' * psi;

julia> Matrix(symmetric_reduced_state(coordinates, 2, 3; keep=1))
2×2 Matrix{Float64}:
 0.666667  0.0
 0.0       0.333333
```

The direct algorithm uses `O(D_trace*D_keep^2)` arithmetic and avoids the
ambient `d^N × d^N` density matrix. Exact combinatorial sizes and conservative
work/output limits are checked before allocation.
"""
function symmetric_reduced_state(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    local_dimension,
    parties;
    keep,
    sparse_output::Bool=issparse(state),
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_nonzeros=_SYMMETRIC_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    Base.require_one_based_indexing(state)
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    kept = _nonnegative_int(keep, "keep")
    kept <= count || throw(ArgumentError("keep=$kept must not exceed parties=$count"))
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _symmetric_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symmetric_limit(max_work, "max_work")
    preliminary_work = BigInt(dimension) + BigInt(count) + BigInt(kept) + 1
    _symmetric_check_resource(
        preliminary_work,
        work_limit,
        "max_work",
        "preflight scalar operations",
        "symmetric_reduced_state",
    )
    global_exact = _symmetric_dimension(dimension, count)
    _symmetric_check_resource(
        global_exact,
        occupation_limit,
        "max_occupations",
        "global occupation sectors",
        "symmetric_reduced_state",
    )
    global_dimension = _symmetric_array_length(global_exact, "global symmetric dimension")
    if state isa AbstractVector
        length(state) == global_dimension || throw(
            DimensionMismatch(
                "state has length $(length(state)); expected $global_dimension symmetric coordinates",
            ),
        )
    else
        size(state) == (global_dimension, global_dimension) || throw(
            DimensionMismatch(
                "operator has size $(size(state)); expected " *
                "($global_dimension, $global_dimension)",
            ),
        )
    end
    input_entries = _symmetric_input_entries(state)

    if kept == 0 || kept == count
        output_exact = kept == 0 ? BigInt(1) : global_exact
        output_entries = output_exact^2
        planned_stored = if kept == 0
            BigInt(1)
        elseif state isa AbstractVector && issparse(state)
            input_entries^2
        elseif state isa AbstractMatrix &&
            (issparse(state) || _symmetric_is_diagonal_storage(state))
            input_entries
        else
            output_entries
        end
        endpoint_work = preliminary_work + planned_stored + global_exact + input_entries
        _symmetric_check_resource(
            endpoint_work,
            work_limit,
            "max_work",
            "estimated endpoint operations",
            "symmetric_reduced_state",
        )
        if sparse_output
            _symmetric_check_resource(
                planned_stored,
                nonzero_limit,
                "max_nonzeros",
                "possible output entries",
                "symmetric_reduced_state",
            )
        else
            _symmetric_check_resource(
                output_entries,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_reduced_state",
            )
        end
        input_type = _symmetric_promoted_eltype(state)
        _symmetric_require_finite(state, "symmetric_reduced_state")
        target_precision = _symmetric_target_precision(state)
        return _symmetric_reduction_endpoint(
            state, input_type, kept, count, sparse_output, target_precision
        )
    end

    _symmetric_check_resource(
        preliminary_work + input_entries,
        work_limit,
        "max_work",
        "input-validation operations",
        "symmetric_reduced_state",
    )
    input_type = _symmetric_promoted_eltype(state)
    _symmetric_require_finite(state, "symmetric_reduced_state")
    target_precision = _symmetric_target_precision(state)
    output_type = _symmetric_root_type(input_type)
    if all(iszero, _symmetric_stored_values(state))
        kept_exact = _symmetric_dimension(dimension, kept)
        _symmetric_check_resource(
            kept_exact,
            occupation_limit,
            "max_occupations",
            "kept occupation sectors",
            "symmetric_reduced_state",
        )
        kept_dimension = _symmetric_array_length(kept_exact, "kept symmetric dimension")
        zero_work = preliminary_work + input_entries + kept_exact + 1
        if !sparse_output
            output_entries = kept_exact^2
            _symmetric_check_resource(
                output_entries,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_reduced_state",
            )
            zero_work += output_entries
        end
        _symmetric_check_resource(
            zero_work,
            work_limit,
            "max_work",
            "zero-output construction operations",
            "symmetric_reduced_state",
        )
        return setprecision(BigFloat, _symmetric_work_precision(state)) do
            zero_result = spzeros(output_type, kept_dimension, kept_dimension)
            finished = sparse_output ? zero_result : Matrix(zero_result)
            return _symmetric_finish_reduction(
                state, finished, input_type, target_precision
            )
        end
    end

    plan = _symmetric_split_preflight(
        dimension,
        count,
        kept,
        "symmetric_reduced_state";
        max_occupations=max_occupations,
        max_nonzeros=max_nonzeros,
        max_work=max_work,
    )
    output_entries = plan.kept_exact^2
    reduction_work = plan.traced_exact * output_entries + plan.pair_count
    work = plan.work + reduction_work + input_entries
    _symmetric_check_resource(
        work,
        work_limit,
        "max_work",
        "estimated reduction operations",
        "symmetric_reduced_state",
    )
    if sparse_output
        _symmetric_check_resource(
            output_entries,
            nonzero_limit,
            "max_nonzeros",
            "possible output entries",
            "symmetric_reduced_state",
        )
    else
        _symmetric_check_resource(
            output_entries,
            dense_limit,
            "max_dense_entries",
            "dense output entries",
            "symmetric_reduced_state",
        )
        if state isa AbstractVector && kept != 0 && kept != count
            _symmetric_check_resource(
                plan.pair_count,
                dense_limit,
                "max_dense_entries",
                "dense amplitude-workspace entries",
                "symmetric_reduced_state",
            )
        end
    end

    weight_type, force_wide = _symmetric_reduction_weight_plan(
        state, input_type, dimension, count, plan.traced_dimension
    )
    return setprecision(BigFloat, _symmetric_work_precision(state)) do
        split = _symmetric_split_data(dimension, count, kept, weight_type, plan)
        use_wide = force_wide || _symmetric_needs_wide_reduction(state, split, output_type)
        working_type = use_wide ? _symmetric_wide_numeric_type(output_type) : output_type
        if sparse_output
            result = _symmetric_sparse_reduction(state, split, working_type)
            finished = use_wide ? _symmetric_narrow_matrix(result, output_type) : result
            return _symmetric_finish_reduction(
                state, finished, input_type, target_precision
            )
        end
        if state isa AbstractVector
            amplitudes = Matrix{working_type}(
                undef, split.kept_dimension, split.traced_dimension
            )
            for traced_index in 1:split.traced_dimension,
                kept_index in 1:split.kept_dimension

                amplitudes[kept_index, traced_index] =
                    convert(working_type, split.weights[kept_index, traced_index]) *
                    convert(working_type, state[split.indices[kept_index, traced_index]])
            end
            result = amplitudes * adjoint(amplitudes)
            finished = use_wide ? _symmetric_narrow_matrix(result, output_type) : result
            return _symmetric_finish_reduction(
                state, finished, input_type, target_precision
            )
        end
        result = zeros(working_type, split.kept_dimension, split.kept_dimension)
        for traced_index in 1:split.traced_dimension,
            column in 1:split.kept_dimension,
            row in 1:split.kept_dimension

            result[row, column] +=
                convert(working_type, split.weights[row, traced_index]) *
                convert(
                    working_type,
                    state[
                        split.indices[row, traced_index],
                        split.indices[column, traced_index],
                    ],
                ) *
                convert(working_type, split.weights[column, traced_index])
        end
        finished = use_wide ? _symmetric_narrow_matrix(result, output_type) : result
        return _symmetric_finish_reduction(state, finished, input_type, target_precision)
    end
end

"""
    symmetric_maximally_mixed_state(
        local_dimension, parties; representation=:coordinates,
        T=Rational{BigInt}, sparse_output=true,
        max_occupations=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

Construct the maximally mixed state on the symmetric subspace. With
`representation=:coordinates` (default), return `I_D/D` in occupation
coordinates. With `representation=:ambient`, return `P_sym/D` in the full
`local_dimension^parties` tensor space. Rational output has exact unit trace;
floating output uses the representable rounding of `1/D`. Neither
representation is a claim about separability.

`T=Rational{BigInt}` gives exact coordinate entries and exact ambient projector
entries. Accumulate a large low-precision floating diagonal in a wider type
when checking its approximate trace. Ambient output is explicitly selected
because it can be exponentially larger and is guarded before construction.
`parties=0` returns the exact `1 × 1` vacuum state.

# Example
```jldoctest
julia> using LinearAlgebra

julia> rho = symmetric_maximally_mixed_state(2, 3);

julia> (size(rho), tr(rho))
((4, 4), 1//1)
```

Coordinate construction costs `O(D)` sparse storage. Ambient construction has
the same guarded complexity as [`symmetric_projector`](@ref).
"""
function symmetric_maximally_mixed_state(
    local_dimension,
    parties;
    representation::Symbol=:coordinates,
    T::Type{<:Number}=Rational{BigInt},
    sparse_output::Bool=true,
    max_occupations=_SYMMETRIC_DEFAULT_MAX_OCCUPATIONS,
    max_nonzeros=_SYMMETRIC_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_SYMMETRIC_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_SYMMETRIC_DEFAULT_MAX_WORK,
)
    representation in (:coordinates, :ambient) ||
        throw(ArgumentError("representation must be :coordinates or :ambient"))
    dimension = _positive_int(local_dimension, "local_dimension")
    count = _nonnegative_int(parties, "parties")
    occupation_limit = _symmetric_limit(max_occupations, "max_occupations")
    nonzero_limit = _symmetric_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _symmetric_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symmetric_limit(max_work, "max_work")
    _symmetric_check_resource(
        BigInt(dimension) + BigInt(count) + 1,
        work_limit,
        "max_work",
        "preflight scalar operations",
        "symmetric_maximally_mixed_state",
    )
    sector_exact = _symmetric_dimension(dimension, count)
    _symmetric_check_resource(
        sector_exact,
        occupation_limit,
        "max_occupations",
        "symmetric occupation sectors",
        "symmetric_maximally_mixed_state",
    )
    sector_count = _symmetric_array_length(sector_exact, "symmetric subspace dimension")
    work = sector_exact + BigInt(count) + 1
    _symmetric_check_resource(
        work,
        work_limit,
        "max_work",
        "estimated scalar operations",
        "symmetric_maximally_mixed_state",
    )
    value = _symmetric_reciprocal(sector_exact, T)

    if representation == :coordinates
        _symmetric_check_resource(
            sector_exact,
            nonzero_limit,
            "max_nonzeros",
            "stored diagonal entries",
            "symmetric_maximally_mixed_state",
        )
        if !sparse_output
            _symmetric_check_resource(
                sector_exact^2,
                dense_limit,
                "max_dense_entries",
                "dense output entries",
                "symmetric_maximally_mixed_state",
            )
        end
        result = spdiagm(0 => fill(value, sector_count))
        return sparse_output ? result : Matrix(result)
    end
    if count == 0
        result = sparse([1], [1], [value], 1, 1)
        return sparse_output ? result : Matrix(result)
    end

    ambient_dimension =
        dimension == 1 ? 1 : _checked_power(dimension, count, "local_dimension^parties")
    if !sparse_output
        _symmetric_check_resource(
            BigInt(ambient_dimension)^2,
            dense_limit,
            "max_dense_entries",
            "dense ambient output entries",
            "symmetric_maximally_mixed_state",
        )
    end
    projector = symmetric_projector(
        dimension,
        count;
        T=T,
        sparse_output=true,
        max_columns=nothing,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
    result = projector * value
    return sparse_output ? result : Matrix(result)
end
