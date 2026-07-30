# Source-informed independent Julia implementation based on the specification
# and QETLAB IsEntanglingGate.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2012 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

const _ENTANGLING_GATE_DEFAULT_MAX_PERMUTATIONS = 40_320
const _ENTANGLING_GATE_DEFAULT_MAX_WITNESSES = 1_000_000
const _ENTANGLING_GATE_DEFAULT_MAX_WORK = 100_000_000

"""
    EntanglingGateResult

Certificate-aware result returned by [`is_entangling_gate`](@ref).

`status` is one of:

- `:not_entangling`: `local_factors` and `subsystem_permutation` reconstruct
  the gate as a subsystem permutation composed with a tensor product;
- `:entangling`: `product_witness` is a normalized product input whose output
  has the reported outside-tolerance `certificate_analysis`;
- `:unknown`: a unitarity or product-structure decision landed on a numerical
  boundary, or the exhaustive sparse-witness search was numerically
  inconclusive.

For a `:not_entangling` result, the checked reconstruction is

```julia
permutation_operator(
    result.factor_output_dims,
    result.subsystem_permutation;
    T=eltype(gate),
    sparse_output=false,
) * tensor_product(result.local_factors...)
```

`input_dims` describe the gate input and `output_dims` its output. The
permutation acts on `factor_output_dims` and restores `output_dims`.
"""
struct EntanglingGateResult{P,F,D,W,A,R,T,ID,OD}
    status::Symbol
    subsystem_permutation::P
    local_factors::F
    factor_output_dims::D
    product_witness::W
    certificate_analysis::A
    unitarity_residual::R
    unitarity_tolerance::T
    permutations_checked::Int
    witnesses_checked::Int
    input_dims::ID
    output_dims::OD
    message::String
end

function Base.show(io::IO, result::EntanglingGateResult)
    return print(
        io,
        "EntanglingGateResult(status=",
        result.status,
        ", permutations_checked=",
        result.permutations_checked,
        ", witnesses_checked=",
        result.witnesses_checked,
        ")",
    )
end

function _entangling_gate_layouts(gate::AbstractMatrix, dims)
    Base.require_one_based_indexing(gate)
    size(gate, 1) == size(gate, 2) || throw(
        DimensionMismatch(
            "is_entangling_gate requires a square gate; got size $(size(gate))"
        ),
    )
    dimension = size(gate, 1)
    row_layout, column_layout = if dims === nothing
        local_dimension = isqrt(dimension)
        local_dimension^2 == dimension || throw(
            ArgumentError(
                "the default layout requires a perfect-square gate dimension; " *
                "provide dims explicitly for a $dimension-dimensional gate",
            ),
        )
        layout = SubsystemLayout((local_dimension, local_dimension))
        (layout, layout)
    elseif dims isa Integer
        first_dimension = _positive_int(dims, "dims")
        rem(dimension, first_dimension) == 0 || throw(
            DimensionMismatch(
                "scalar dims=$first_dimension does not evenly divide gate " *
                "dimension $dimension",
            ),
        )
        layout = SubsystemLayout((first_dimension, div(dimension, first_dimension)))
        (layout, layout)
    elseif dims isa AbstractMatrix
        Base.require_one_based_indexing(dims)
        size(dims, 1) == 2 || throw(
            ArgumentError(
                "matrix dims must have two rows (output dimensions above input " *
                "dimensions); got size $(size(dims))",
            ),
        )
        size(dims, 2) >= 2 ||
            throw(ArgumentError("is_entangling_gate requires at least two subsystems"))
        output_dimensions = ntuple(
            index -> _positive_int(dims[1, index], "dims[1,$index]"), size(dims, 2)
        )
        input_dimensions = ntuple(
            index -> _positive_int(dims[2, index], "dims[2,$index]"), size(dims, 2)
        )
        (SubsystemLayout(output_dimensions), SubsystemLayout(input_dimensions))
    else
        layout = _as_layout(dims)
        (layout, layout)
    end
    length(row_layout) == length(column_layout) || throw(
        DimensionMismatch(
            "output dims has $(length(row_layout)) subsystems but input dims has " *
            "$(length(column_layout))",
        ),
    )
    length(row_layout) >= 2 ||
        throw(ArgumentError("is_entangling_gate requires at least two subsystems"))
    row_layout.total_dimension == dimension || throw(
        DimensionMismatch(
            "output dims $(row_layout.dims) have product " *
            "$(row_layout.total_dimension), expected $dimension",
        ),
    )
    column_layout.total_dimension == dimension || throw(
        DimensionMismatch(
            "input dims $(column_layout.dims) have product " *
            "$(column_layout.total_dimension), expected $dimension",
        ),
    )
    return row_layout, column_layout
end

function _entangling_gate_limit(value, name::AbstractString)
    value === nothing && return nothing
    return _nonnegative_int(value, name)
end

function _entangling_gate_checked_product(left::Int, right::Int, label::AbstractString)
    return try
        Base.checked_mul(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("$label exceeds typemax(Int)"))
    end
end

function _entangling_gate_checked_sum(left::Int, right::Int, label::AbstractString)
    return try
        Base.checked_add(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("$label exceeds typemax(Int)"))
    end
end

function _entangling_gate_permutation_count(subsystem_count::Int)
    count = 1
    for factor in 2:subsystem_count
        count = _entangling_gate_checked_product(
            count, factor, "the subsystem-permutation count"
        )
    end
    return count
end

function _entangling_gate_witness_count(layout::SubsystemLayout)
    count = 0
    total_dimension = layout.total_dimension
    for local_dimension in layout.dims
        local_choices =
            _entangling_gate_checked_product(
                local_dimension,
                local_dimension + 1,
                "the local sparse-witness choice count",
            ) ÷ 2
        assignments = div(total_dimension, local_dimension)
        count = _entangling_gate_checked_sum(
            count,
            _entangling_gate_checked_product(
                assignments, local_choices, "the sparse-witness count"
            ),
            "the sparse-witness count",
        )
    end
    return count
end

function _entangling_gate_next_permutation!(permutation::Vector{Int})
    pivot = length(permutation) - 1
    while pivot >= 1 && permutation[pivot] >= permutation[pivot + 1]
        pivot -= 1
    end
    pivot == 0 && return false
    successor = length(permutation)
    while permutation[successor] <= permutation[pivot]
        successor -= 1
    end
    permutation[pivot], permutation[successor] = permutation[successor], permutation[pivot]
    reverse!(permutation, pivot + 1, length(permutation))
    return true
end

function _entangling_gate_candidate(
    ::Type{T},
    total_dimension::Int,
    prefix_index::Int,
    local_dimension::Int,
    first_local_index::Int,
    second_local_index::Int,
    suffix_dimension::Int,
    suffix_index::Int,
) where {T<:LinearAlgebra.BlasFloat}
    first_index =
        ((prefix_index - 1) * local_dimension + first_local_index - 1) * suffix_dimension +
        suffix_index
    if first_local_index == second_local_index
        return sparsevec([first_index], T[one(T)], total_dimension)
    end
    second_index =
        ((prefix_index - 1) * local_dimension + second_local_index - 1) * suffix_dimension +
        suffix_index
    real_type = typeof(real(zero(T)))
    coefficient = convert(T, inv(sqrt(real_type(2))))
    return sparsevec(
        [first_index, second_index], T[coefficient, coefficient], total_dimension
    )
end

function _entangling_gate_local_candidates(
    ::Type{T}, local_dimension::Int
) where {T<:LinearAlgebra.BlasFloat}
    witness_type = T <: Real ? Complex{T} : T
    candidates = SparseVector{witness_type,Int}[]
    for index in 1:local_dimension
        push!(
            candidates, sparsevec([index], witness_type[one(witness_type)], local_dimension)
        )
    end
    real_type = typeof(real(zero(witness_type)))
    coefficient = convert(witness_type, inv(sqrt(real_type(2))))
    phases = (
        one(witness_type),
        -one(witness_type),
        convert(witness_type, im),
        convert(witness_type, -im),
    )
    for phase in phases
        for first_index in 1:(local_dimension - 1)
            for second_index in (first_index + 1):local_dimension
                push!(
                    candidates,
                    sparsevec(
                        [first_index, second_index],
                        witness_type[coefficient, coefficient * phase],
                        local_dimension,
                    ),
                )
            end
        end
    end
    return candidates
end

function _entangling_gate_advance_odometer!(indices::Vector{Int}, candidate_sets::Vector)
    position = length(indices)
    while position >= 1
        if indices[position] < length(candidate_sets[position])
            indices[position] += 1
            return true
        end
        indices[position] = 1
        position -= 1
    end
    return false
end

function _entangling_gate_result(
    status::Symbol,
    permutation,
    factors,
    factor_output_dims,
    witness,
    certificate,
    unitarity_residual,
    unitarity_tolerance,
    permutations_checked::Int,
    witnesses_checked::Int,
    column_layout::SubsystemLayout,
    row_layout::SubsystemLayout,
    message::String,
)
    return EntanglingGateResult(
        status,
        permutation,
        factors,
        factor_output_dims,
        witness,
        certificate,
        unitarity_residual,
        unitarity_tolerance,
        permutations_checked,
        witnesses_checked,
        column_layout.dims,
        row_layout.dims,
        message,
    )
end

"""
    is_entangling_gate(
        gate,
        dims=nothing;
        atol=nothing,
        rtol=nothing,
        allow_densify=false,
        max_permutations=40_320,
        max_witnesses=1_000_000,
        max_work=100_000_000,
    ) -> EntanglingGateResult

Determine whether a multipartite unitary can entangle a product input.

`dims` is either one shared input/output layout, a scalar first bipartite
dimension, or a two-row matrix whose first row gives output dimensions and
second row gives input dimensions. `nothing` selects two equal subsystems and
therefore requires a perfect-square gate dimension. Subsystem `1` is the
slowest-varying tensor factor.

The method first checks every output-subsystem permutation for a local tensor
factorization. A successful factorization certifies `:not_entangling`. If none
is found, it first tests the pinned family of normalized product inputs having
one or two nonzero coefficients. Because that upstream search misses gates
such as controlled phase, the corrected path then exhausts a finite local
phase grid made from basis vectors and normalized
`eᵢ + z*eⱼ`, `z in (1, -1, im, -im)`, product inputs. An outside-tolerance
output product analysis supplies an explicit `:entangling` witness. Numerical
boundary or inconclusive cases return `:unknown`; they are never collapsed to
a Boolean.

The gate must be finite, square, and unitary within
`atol + rtol*sqrt(size(gate,1))`. The same tolerances govern product analyses;
defaults are zero absolute tolerance and `sqrt(eps(T))` relative tolerance.
Sparse input requires `allow_densify=true`. Only BLAS floating element types
are accepted, and the input is never normalized, projected, or mutated.

All factorial permutation work and sparse-witness work are counted before
their respective phases. The three guards may be raised or disabled explicitly
with `nothing`. The dense algorithm has worst-case spectral cost polynomial in
the gate dimension per subsystem permutation, plus one dense matrix-vector
product per witness.
"""
function is_entangling_gate(
    gate::AbstractMatrix{<:Number},
    dims=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_permutations=_ENTANGLING_GATE_DEFAULT_MAX_PERMUTATIONS,
    max_witnesses=_ENTANGLING_GATE_DEFAULT_MAX_WITNESSES,
    max_work=_ENTANGLING_GATE_DEFAULT_MAX_WORK,
)
    row_layout, column_layout = _entangling_gate_layouts(gate, dims)
    dense, real_type = _tierd_dense_matrix(
        gate; allow_densify=allow_densify, operation="is_entangling_gate"
    )
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    subsystem_count = length(row_layout)
    permutation_count = _entangling_gate_permutation_count(subsystem_count)
    permutation_limit = _entangling_gate_limit(max_permutations, "max_permutations")
    if permutation_limit !== nothing && permutation_count > permutation_limit
        throw(
            ArgumentError(
                "is_entangling_gate requires $permutation_count subsystem " *
                "permutations, exceeding max_permutations=$permutation_limit",
            ),
        )
    end
    matrix_entries = _entangling_gate_checked_product(
        size(dense, 1), size(dense, 2), "the gate entry count"
    )
    spectral_work = _entangling_gate_checked_product(
        size(dense, 1), matrix_entries, "one dense spectral-work estimate"
    )
    permutation_spectral_work = _entangling_gate_checked_product(
        permutation_count, spectral_work, "the permutation-phase spectral-work estimate"
    )
    permutation_work = _entangling_gate_checked_sum(
        spectral_work,
        permutation_spectral_work,
        "the unitarity-plus-permutation work estimate",
    )
    work_limit = _entangling_gate_limit(max_work, "max_work")
    if work_limit !== nothing && permutation_work > work_limit
        throw(
            ArgumentError(
                "the unitarity and permutation phases require estimated work " *
                "$permutation_work, exceeding max_work=$work_limit",
            ),
        )
    end

    identity_matrix = Matrix{eltype(dense)}(I, size(dense, 1), size(dense, 1))
    unitarity_residual = norm(adjoint(dense) * dense - identity_matrix)
    unitarity_tolerance = absolute + relative * sqrt(real_type(size(dense, 1)))
    unitarity_residual > unitarity_tolerance && throw(
        DomainError(
            unitarity_residual,
            "gate is not unitary within tolerance $unitarity_tolerance; " *
            "is_entangling_gate never projects an input onto a unitary",
        ),
    )
    if !iszero(unitarity_residual) && unitarity_residual == unitarity_tolerance
        return _entangling_gate_result(
            :unknown,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            unitarity_residual,
            unitarity_tolerance,
            0,
            0,
            column_layout,
            row_layout,
            "gate unitarity lies exactly on the requested numerical boundary",
        )
    end

    permutation = collect(1:subsystem_count)
    permutations_checked = 0
    saw_product_boundary = false
    while true
        permutations_checked += 1
        row_plan = SubsystemPermutationPlan(row_layout, Tuple(permutation))
        permuted_gate = permute_subsystems(dense, row_plan; rows_only=true)
        analysis = is_product_operator(
            permuted_gate,
            row_plan.output_layout,
            column_layout;
            atol=absolute,
            rtol=relative,
            allow_densify=false,
        )
        if analysis.status === :within_tolerance
            inverse_permutation = Tuple(invperm(permutation))
            return _entangling_gate_result(
                :not_entangling,
                inverse_permutation,
                analysis.factors,
                row_plan.output_layout.dims,
                nothing,
                analysis,
                unitarity_residual,
                unitarity_tolerance,
                permutations_checked,
                0,
                column_layout,
                row_layout,
                "gate is a subsystem permutation composed with local factors",
            )
        end
        saw_product_boundary |= analysis.status === :boundary
        _entangling_gate_next_permutation!(permutation) || break
    end

    witness_count = _entangling_gate_witness_count(column_layout)
    witness_limit = _entangling_gate_limit(max_witnesses, "max_witnesses")
    if witness_limit !== nothing && witness_count > witness_limit
        throw(
            ArgumentError(
                "is_entangling_gate requires $witness_count sparse product-state " *
                "candidates, exceeding max_witnesses=$witness_limit",
            ),
        )
    end
    witness_work = _entangling_gate_checked_product(
        witness_count, matrix_entries, "the witness-phase work estimate"
    )
    total_work = _entangling_gate_checked_sum(
        permutation_work, witness_work, "the total work estimate"
    )
    if work_limit !== nothing && total_work > work_limit
        throw(
            ArgumentError(
                "is_entangling_gate requires estimated total work $total_work, " *
                "exceeding max_work=$work_limit",
            ),
        )
    end

    witnesses_checked = 0
    saw_witness_boundary = false
    total_dimension = column_layout.total_dimension
    for subsystem in 1:length(column_layout)
        local_dimension = column_layout[subsystem]
        prefix_dimension = _checked_product(
            column_layout.dims[1:(subsystem - 1)], "the product-witness prefix dimensions"
        )
        suffix_dimension = _checked_product(
            column_layout.dims[(subsystem + 1):end], "the product-witness suffix dimensions"
        )
        for first_local_index in 1:local_dimension
            for prefix_index in 1:prefix_dimension, suffix_index in 1:suffix_dimension
                witness = _entangling_gate_candidate(
                    eltype(dense),
                    total_dimension,
                    prefix_index,
                    local_dimension,
                    first_local_index,
                    first_local_index,
                    suffix_dimension,
                    suffix_index,
                )
                witnesses_checked += 1
                output_analysis = is_product_vector(
                    dense * witness,
                    row_layout;
                    atol=absolute,
                    rtol=relative,
                    allow_densify=false,
                )
                if output_analysis.status === :outside_tolerance
                    return _entangling_gate_result(
                        :entangling,
                        nothing,
                        nothing,
                        nothing,
                        witness,
                        output_analysis,
                        unitarity_residual,
                        unitarity_tolerance,
                        permutations_checked,
                        witnesses_checked,
                        column_layout,
                        row_layout,
                        "a normalized one-support product input has an entangled output",
                    )
                end
                saw_witness_boundary |= output_analysis.status === :boundary
            end
        end
        for first_local_index in 1:(local_dimension - 1)
            for second_local_index in (first_local_index + 1):local_dimension
                for prefix_index in 1:prefix_dimension, suffix_index in 1:suffix_dimension
                    witness = _entangling_gate_candidate(
                        eltype(dense),
                        total_dimension,
                        prefix_index,
                        local_dimension,
                        first_local_index,
                        second_local_index,
                        suffix_dimension,
                        suffix_index,
                    )
                    witnesses_checked += 1
                    output_analysis = is_product_vector(
                        dense * witness,
                        row_layout;
                        atol=absolute,
                        rtol=relative,
                        allow_densify=false,
                    )
                    if output_analysis.status === :outside_tolerance
                        return _entangling_gate_result(
                            :entangling,
                            nothing,
                            nothing,
                            nothing,
                            witness,
                            output_analysis,
                            unitarity_residual,
                            unitarity_tolerance,
                            permutations_checked,
                            witnesses_checked,
                            column_layout,
                            row_layout,
                            "a normalized two-support product input has an entangled output",
                        )
                    end
                    saw_witness_boundary |= output_analysis.status === :boundary
                end
            end
        end
    end

    grid_count = 1
    for local_dimension in column_layout.dims
        twice_dimension = _entangling_gate_checked_product(
            2, local_dimension, "the extended local witness count"
        )
        local_candidate_count = _entangling_gate_checked_product(
            local_dimension, twice_dimension - 1, "the extended local witness count"
        )
        grid_count = _entangling_gate_checked_product(
            grid_count, local_candidate_count, "the extended witness-grid count"
        )
    end
    total_candidate_count = _entangling_gate_checked_sum(
        witness_count, grid_count, "the total witness candidate count"
    )
    if witness_limit !== nothing && total_candidate_count > witness_limit
        throw(
            ArgumentError(
                "the corrected extended search requires up to " *
                "$total_candidate_count product-state candidates, exceeding " *
                "max_witnesses=$witness_limit",
            ),
        )
    end
    extended_witness_work = _entangling_gate_checked_product(
        grid_count, matrix_entries, "the extended witness-grid work estimate"
    )
    extended_total_work = _entangling_gate_checked_sum(
        total_work, extended_witness_work, "the extended total work estimate"
    )
    if work_limit !== nothing && extended_total_work > work_limit
        throw(
            ArgumentError(
                "the corrected extended search requires estimated total work " *
                "$extended_total_work, exceeding max_work=$work_limit",
            ),
        )
    end

    local_candidate_sets = [
        _entangling_gate_local_candidates(eltype(dense), local_dimension) for
        local_dimension in column_layout.dims
    ]
    candidate_indices = ones(Int, length(local_candidate_sets))
    while true
        local_factors = ntuple(
            subsystem -> local_candidate_sets[subsystem][candidate_indices[subsystem]],
            length(local_candidate_sets),
        )
        witness = tensor_product(local_factors...)
        witnesses_checked += 1
        output_analysis = is_product_vector(
            dense * witness, row_layout; atol=absolute, rtol=relative, allow_densify=false
        )
        if output_analysis.status === :outside_tolerance
            return _entangling_gate_result(
                :entangling,
                nothing,
                nothing,
                nothing,
                witness,
                output_analysis,
                unitarity_residual,
                unitarity_tolerance,
                permutations_checked,
                witnesses_checked,
                column_layout,
                row_layout,
                "a normalized finite-phase-grid product input has an entangled output",
            )
        end
        saw_witness_boundary |= output_analysis.status === :boundary
        _entangling_gate_advance_odometer!(candidate_indices, local_candidate_sets) || break
    end
    boundary_message = if saw_product_boundary || saw_witness_boundary
        "a product-factorization or witness output lies on a numerical boundary"
    else
        "the exhaustive sparse-witness family was numerically inconclusive"
    end
    return _entangling_gate_result(
        :unknown,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        unitarity_residual,
        unitarity_tolerance,
        permutations_checked,
        witnesses_checked,
        column_layout,
        row_layout,
        boundary_message,
    )
end
