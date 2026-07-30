# Source-informed independent Julia implementation based on the specification
# and QETLAB MinUPBSize.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2013 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    MinimumUPBSizeResult

Exact-or-unknown result returned by [`minimum_upb_size`](@ref).

`status == :known` means `size` is an exact theorem-backed minimum for
`dimensions`. `status == :unknown` means only that the reviewed pinned theorem
table has no exact entry; `size` is then `nothing`. `lower_bound` is always the
Alon--Lovász counting bound
`1 + sum(dimension - 1 for dimension in dimensions)`.

For known cases, `reference_key`, `reference`, and `reference_url` identify
the primary source used by the pinned table. The dimensions are stored in
canonical sorted order because subsystem permutations do not change the
minimum.
"""
struct MinimumUPBSizeResult{S,D,R,U}
    status::Symbol
    size::S
    lower_bound::Int
    dimensions::D
    reference_key::Symbol
    reference::R
    reference_url::U
    message::String
end

function Base.show(io::IO, result::MinimumUPBSizeResult)
    return print(
        io,
        "MinimumUPBSizeResult(status=",
        result.status,
        ", size=",
        result.size,
        ", lower_bound=",
        result.lower_bound,
        ", dimensions=",
        result.dimensions,
        ")",
    )
end

function _upb_reference(reference_key::Symbol)
    reference_key === :divincenzo_mor_shor_smolin_terhal_2003 && return (
        "D. P. DiVincenzo, T. Mor, P. W. Shor, J. A. Smolin, and " *
        "B. M. Terhal, Commun. Math. Phys. 238, 379–410 (2003).",
        "https://arxiv.org/abs/quant-ph/9908070",
    )
    reference_key === :alon_lovasz_2001 && return (
        "N. Alon and L. Lovász, J. Combin. Theory Ser. A 95, 169–179 (2001).",
        "https://doi.org/10.1006/jcta.2000.3122",
    )
    reference_key === :feng_2006 && return (
        "K. Feng, Discrete Appl. Math. 154, 942–949 (2006).",
        "https://doi.org/10.1016/j.dam.2005.10.011",
    )
    reference_key === :chen_johnston_2015 && return (
        "J. Chen and N. Johnston, Commun. Math. Phys. 333, 351–365 (2015).",
        "https://arxiv.org/abs/1301.1406",
    )
    reference_key === :johnston_2013 && return (
        "N. Johnston, Proc. TQC 2013, LIPIcs 22, 93–105 (2013).",
        "https://arxiv.org/abs/1302.1604",
    )
    return error("internal minimum-UPB reference key is unknown: $reference_key")
end

function _minimum_upb_dimensions(dims)
    dims isa Tuple ||
        dims isa AbstractVector ||
        throw(
            ArgumentError(
                "dims must be a tuple or vector of at least two local dimensions; " *
                "got $(typeof(dims))",
            ),
        )
    dims isa AbstractVector && Base.require_one_based_indexing(dims)
    length(dims) >= 2 ||
        throw(ArgumentError("minimum_upb_size requires at least two subsystems"))
    checked = [
        _positive_int(dimension, "dims[$index]") for (index, dimension) in pairs(dims)
    ]
    all(dimension -> dimension >= 2, checked) || throw(
        ArgumentError(
            "every local dimension must be at least 2; one-dimensional factors " *
            "must be removed before requesting a UPB minimum",
        ),
    )
    sort!(checked)
    return Tuple(checked)
end

function _minimum_upb_lower_bound(dimensions::Tuple)
    total = 1
    for dimension in dimensions
        total = try
            Base.checked_add(total, dimension - 1)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("the UPB counting lower bound exceeds typemax(Int)"))
        end
    end
    return total
end

function _known_minimum_upb_result(
    size::Int, lower_bound::Int, dimensions::Tuple, reference_key::Symbol
)
    reference, reference_url = _upb_reference(reference_key)
    return MinimumUPBSizeResult(
        :known,
        size,
        lower_bound,
        dimensions,
        reference_key,
        reference,
        reference_url,
        "the reviewed theorem table gives an exact minimum",
    )
end

"""
    minimum_upb_size(dims) -> MinimumUPBSizeResult

Return the exact minimum cardinality of an orthogonal unextendible product
basis when that value is covered by the theorem and exception table reviewed
for the pinned QETLAB revision.

`dims` contains at least two local Hilbert-space dimensions, each at least
two. Subsystem order is irrelevant. Known results include the complete
bipartite table, the Alon--Lovász lower-bound cases, the reviewed qubit
families, and the pinned `(2,2,d)` exceptions. An uncovered multipartite
family returns `status == :unknown`, `size === nothing`, and the universal
counting lower bound. It does not throw or guess an exact minimum.

This is a constant-work arithmetic/table lookup after sorting `length(dims)`
integers. Checked integer arithmetic rejects an unrepresentable bound. The
result records a primary-source citation and URL for every exact branch.

# Examples

```jldoctest
julia> minimum_upb_size((3, 3)).size
5

julia> unresolved = minimum_upb_size((2, 3, 4));

julia> (unresolved.status, unresolved.size, unresolved.lower_bound)
(:unknown, nothing, 7)
```
"""
function minimum_upb_size(dims)
    dimensions = _minimum_upb_dimensions(dims)
    subsystem_count = length(dimensions)
    lower_bound = _minimum_upb_lower_bound(dimensions)
    largest_dimension = dimensions[end]

    if subsystem_count == 2 && dimensions[1] <= 2
        exact_size = _checked_product(dimensions, "dims")
        return _known_minimum_upb_result(
            exact_size, lower_bound, dimensions, :divincenzo_mor_shor_smolin_terhal_2003
        )
    elseif iseven(lower_bound) || all(isodd, dimensions)
        return _known_minimum_upb_result(
            lower_bound, lower_bound, dimensions, :alon_lovasz_2001
        )
    elseif subsystem_count == 2 &&
        largest_dimension - 1 >= lower_bound - largest_dimension &&
        lower_bound - largest_dimension >= 3
        return _known_minimum_upb_result(
            Base.checked_add(lower_bound, 1), lower_bound, dimensions, :chen_johnston_2015
        )
    elseif (subsystem_count == 4 && dimensions == (2, 2, 2, 2)) ||
        (rem(subsystem_count, 4) == 2 && all(==(2), dimensions))
        return _known_minimum_upb_result(
            Base.checked_add(lower_bound, 1), lower_bound, dimensions, :feng_2006
        )
    elseif subsystem_count == 8 && dimensions == ntuple(_ -> 2, 8)
        return _known_minimum_upb_result(11, lower_bound, dimensions, :johnston_2013)
    elseif all(==(2), dimensions)
        return _known_minimum_upb_result(
            Base.checked_add(lower_bound, 3), lower_bound, dimensions, :johnston_2013
        )
    elseif subsystem_count == 3 && (dimensions == (2, 2, 3) || dimensions == (2, 2, 5))
        return _known_minimum_upb_result(
            Base.checked_add(lower_bound, 1), lower_bound, dimensions, :feng_2006
        )
    elseif subsystem_count == 3 &&
        dimensions[1] == 2 &&
        dimensions[2] == 2 &&
        rem(dimensions[3], 4) == 1
        return _known_minimum_upb_result(
            Base.checked_add(lower_bound, 1), lower_bound, dimensions, :chen_johnston_2015
        )
    end

    return MinimumUPBSizeResult(
        :unknown,
        nothing,
        lower_bound,
        dimensions,
        :unresolved_by_reviewed_table,
        nothing,
        nothing,
        "the reviewed pinned theorem table does not determine this exact minimum",
    )
end
