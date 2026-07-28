# Source-informed independent Julia implementation based on the specification
# and QETLAB Tensor.m, TensorSum.m, and KroneckerSum.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    tensor_product(A, B, ...)
    tensor_product(A; copies=1)

Compute a Kronecker tensor product without changing element types or
densifying sparse factors.

Factors are ordered conventionally: in `tensor_product(A, B)`, `A` is
subsystem 1 and `B` is subsystem 2.  Consequently the basis index of `B`
varies fastest.  Use [`tensor_power`](@ref) when the copy count is computed
dynamically.
"""
function tensor_product(
    first_factor::Union{AbstractVector,AbstractMatrix},
    remaining_factors::Union{AbstractVector,AbstractMatrix}...;
    copies=nothing,
)
    if copies !== nothing
        isempty(remaining_factors) || throw(
            ArgumentError(
                "copies may only be used when one tensor-product factor is supplied"
            ),
        )
        copies isa Integer || throw(ArgumentError("copies must be a nonnegative integer"))
        return tensor_power(first_factor, copies)
    end

    result = copy(first_factor)
    for factor in remaining_factors
        result = kron(result, factor)
    end
    return result
end

"""
    tensor_power(A, copies)

Return `A` tensored with itself `copies` times.  Zero copies returns the
scalar multiplicative identity of `eltype(A)`, matching the empty tensor
product convention.
"""
function tensor_power(factor::Union{AbstractVector,AbstractMatrix}, copies::Integer)
    count = _nonnegative_int(copies, "copies")
    count == 0 && return one(eltype(factor))
    result = copy(factor)
    for _ in 2:count
        result = kron(result, factor)
    end
    return result
end

function _kronecker_identity(::Type{T}, dimension::Int, sparse_output::Bool) where {T}
    if sparse_output
        return spdiagm(0 => fill(one(T), dimension))
    end
    return Matrix{T}(I, dimension, dimension)
end

"""
    kronecker_sum(A, B, ...)
    kronecker_sum(A; copies=1)

Compute the Kronecker sum of square matrices:

```math
A_1 \\oplus \\cdots \\oplus A_n =
\\sum_k I \\otimes \\cdots \\otimes A_k \\otimes \\cdots \\otimes I.
```

If every factor is sparse, sparse identities and sparse outputs are used.
`copies=0` is rejected because a zero-factor Kronecker sum has no useful
matrix shape.
"""
function kronecker_sum(
    first_factor::AbstractMatrix, remaining_factors::AbstractMatrix...; copies=nothing
)
    factors = if copies === nothing
        (first_factor, remaining_factors...)
    else
        isempty(remaining_factors) || throw(
            ArgumentError(
                "copies may only be used when one Kronecker-sum factor is supplied"
            ),
        )
        copies isa Integer || throw(ArgumentError("copies must be a positive integer"))
        count = _positive_int(copies, "copies")
        ntuple(_ -> first_factor, count)
    end

    for (position, factor) in pairs(factors)
        size(factor, 1) == size(factor, 2) || throw(
            DimensionMismatch(
                "Kronecker-sum factor $position has size $(size(factor)); every factor must be square",
            ),
        )
    end
    length(factors) == 1 && return copy(factors[1])

    promoted_type = promote_type(map(eltype, factors)...)
    sparse_output = all(issparse, factors)
    identities = ntuple(
        position ->
            _kronecker_identity(promoted_type, size(factors[position], 1), sparse_output),
        length(factors),
    )

    result = nothing
    for active in eachindex(factors)
        local_factors = ntuple(
            position -> position == active ? factors[position] : identities[position],
            length(factors),
        )
        term = tensor_product(local_factors...)
        result = result === nothing ? term : result + term
    end
    return result
end

function _tensor_sum_term_count(factor, position::Int)
    if factor isa AbstractMatrix
        return size(factor, 2)
    elseif factor isa AbstractVector{<:Number}
        return 1
    elseif factor isa Tuple || (factor isa AbstractVector && !(eltype(factor) <: Number))
        isempty(factor) &&
            throw(ArgumentError("factor $position contains no decomposition terms"))
        all(term -> term isa AbstractVector || term isa AbstractMatrix, factor) ||
            throw(ArgumentError("factor $position must contain only vectors or matrices"))
        return length(factor)
    end
    return throw(
        ArgumentError(
            "factor $position must be a matrix of column vectors, one vector, or a collection of arrays; got $(typeof(factor))",
        ),
    )
end

function _tensor_sum_term(factor::AbstractMatrix, term::Int)
    return factor[:, term]
end

function _tensor_sum_term(factor::AbstractVector{<:Number}, term::Int)
    term == 1 || throw(BoundsError(factor, term))
    return factor
end

_tensor_sum_term(factor::Union{Tuple,AbstractVector}, term::Int) = factor[term]

function _tensor_sum_weights(weights, term_count::Int)
    if weights === nothing
        return nothing
    elseif weights isa Number
        term_count == 1 || throw(
            DimensionMismatch(
                "a scalar weight is only valid for one term; decomposition has $term_count terms",
            ),
        )
        return (weights,)
    elseif weights isa Tuple || weights isa AbstractVector
        length(weights) == term_count || throw(
            DimensionMismatch(
                "weights has length $(length(weights)); expected $term_count"
            ),
        )
        all(weight -> weight isa Number, weights) ||
            throw(ArgumentError("weights must contain only numbers"))
        return weights
    end
    return throw(
        ArgumentError(
            "weights must be `nothing`, a number, or a tuple/vector of numbers; got $(typeof(weights))",
        ),
    )
end

"""
    tensor_sum(factors...; weights=nothing)

Reconstruct a vector or operator from a tensor decomposition.

Every factor supplies the same number of terms.  A matrix supplies its
columns as vector terms; a numeric vector supplies one term; and a tuple or
vector of arrays supplies those arrays as terms.  The result is

```math
\\sum_k w_k A_{1,k} \\otimes A_{2,k} \\otimes \\cdots.
```

This is the Julia-native counterpart of QETLAB's `TensorSum`; it is not the
Kronecker sum of square matrices.  Sparse terms remain sparse when all
arithmetic operands support sparse addition.
"""
function tensor_sum(first_factor, remaining_factors...; weights=nothing)
    factors = (first_factor, remaining_factors...)
    counts = ntuple(
        position -> _tensor_sum_term_count(factors[position], position), length(factors)
    )
    term_count = counts[1]
    term_count > 0 ||
        throw(ArgumentError("a tensor decomposition must have at least one term"))
    all(==(term_count), counts) || throw(
        DimensionMismatch(
            "all factors must supply the same number of terms; got counts $counts"
        ),
    )
    checked_weights = _tensor_sum_weights(weights, term_count)

    expected_sizes = Vector{Tuple}(undef, length(factors))
    for position in eachindex(factors)
        term = _tensor_sum_term(factors[position], 1)
        expected_sizes[position] = size(term)
    end

    result = nothing
    for term_index in 1:term_count
        local_terms = ntuple(
            position -> begin
                term = _tensor_sum_term(factors[position], term_index)
                size(term) == expected_sizes[position] || throw(
                    DimensionMismatch(
                        "term $term_index of factor $position has size $(size(term)); expected $(expected_sizes[position])",
                    ),
                )
                term
            end,
            length(factors),
        )
        product_term = tensor_product(local_terms...)
        weighted_term = if checked_weights === nothing
            product_term
        else
            checked_weights[term_index] * product_term
        end
        result = result === nothing ? weighted_term : result + weighted_term
    end
    return result
end
