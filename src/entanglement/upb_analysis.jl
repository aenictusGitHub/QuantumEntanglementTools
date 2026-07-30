# Source-informed independent Julia implementation based on the specification
# and QETLAB IsUPB.m and helpers/vec_partitions.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The mathematical definition is also checked against:
# C. H. Bennett et al., Phys. Rev. Lett. 82, 5385 (1999),
# https://arxiv.org/abs/quant-ph/9808030.

const _UPB_DEFAULT_MAX_PARTITIONS = 100_000
const _UPB_DEFAULT_MAX_WORK = 25_000_000
const _UPB_DEFAULT_MAX_DENSE_ENTRIES = 1_000_000
const _UPB_DEFAULT_MAX_STATES = 256

"""
    UPBAnalysisResult

Certificate-aware result returned by [`is_upb`](@ref).

`status` is `:upb`, `:not_upb`, or `:unknown`, and `is_upb` is respectively
`true`, `false`, or `nothing`. Important `reason` values are:

- `:unextendible`: exhaustive partition analysis found no product extension;
- `:extension_witness`: `witness_factors` form a product vector orthogonal to
  every supplied state;
- `:not_product`, `:not_orthogonal`, or `:complete_basis`: the input does not
  satisfy the definition of an incomplete orthogonal product basis;
- `:numerical_boundary`: a product, orthogonality, rank, normalization, or
  witness check fell inside the declared tolerance boundary band;
- `:state_limit`, `:partition_limit`, or `:work_limit`: an explicit
  deterministic size/search budget stopped the analysis.

`certificate_kind == :exact` means that exact arithmetic established the
result. `:tolerance_robust` means that a floating-point result stayed outside
the boundary band. An `:unknown` result always has `certificate_kind == :none`.

The witness is represented primarily by local factors, so it remains compact
even when the tensor-product dimension is large. `witness_vector` is populated
only when `materialize_witness=true` was requested. `witness_partition`
records which supplied state indices were assigned to each local orthogonal
complement. `product_residual` records the largest tensor-factorization
residual for global-vector input, while `witness_residual` records the checked
orthogonality of a returned extension.
"""
struct UPBAnalysisResult{B,F,V,P,D,OR,PR,WR,T,O}
    status::Symbol
    reason::Symbol
    is_upb::B
    witness_factors::F
    witness_vector::V
    witness_partition::P
    dimensions::D
    state_count::Int
    input_form::Symbol
    certificate_kind::Symbol
    orthogonality_residual::OR
    product_residual::PR
    witness_residual::WR
    atol::T
    rtol::T
    boundary_factor::T
    partitions_examined::Int
    work_used::Int
    max_partitions::Union{Nothing,Int}
    max_work::Union{Nothing,Int}
    max_states::Union{Nothing,Int}
    normalization::Symbol
    input_normalized::Union{Bool,Nothing}
    analysis_rescaled::Bool
    densified::Bool
    offending_state::O
    message::String
end

function Base.show(io::IO, result::UPBAnalysisResult)
    return print(
        io,
        "UPBAnalysisResult(status=",
        result.status,
        ", reason=",
        result.reason,
        ", dimensions=",
        result.dimensions,
        ", state_count=",
        result.state_count,
        ", partitions_examined=",
        result.partitions_examined,
        ")",
    )
end

struct _UPBOptions{T}
    atol::T
    rtol::T
    boundary_factor::T
    max_partitions::Union{Nothing,Int}
    max_work::Union{Nothing,Int}
    max_dense_entries::Union{Nothing,Int}
    max_states::Union{Nothing,Int}
    normalization::Symbol
    allow_densify::Bool
    materialize_witness::Bool
    exact::Bool
end

mutable struct _UPBSearchState
    partitions_examined::Int
    work_used::Int
    boundary_seen::Bool
    stopped_reason::Union{Nothing,Symbol}
end

function _upb_checked_nonnegative(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a nonnegative integer or nothing"))
    value isa Integer ||
        throw(ArgumentError("$name must be a nonnegative integer or nothing"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative"))
    value <= typemax(Int) || throw(ArgumentError("$name exceeds typemax(Int)"))
    return Int(value)
end

function _upb_limit(value, name::AbstractString)
    value === nothing && return nothing
    return _upb_checked_nonnegative(value, name)
end

function _upb_checked_add(left::Int, right::Int, label::AbstractString)
    return try
        Base.checked_add(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("$label exceeds typemax(Int)"))
    end
end

function _upb_checked_mul(left::Int, right::Int, label::AbstractString)
    return try
        Base.checked_mul(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("$label exceeds typemax(Int)"))
    end
end

function _upb_checked_product(values, label::AbstractString)
    product_value = 1
    for value in values
        product_value = _upb_checked_mul(product_value, value, label)
    end
    return product_value
end

_upb_exact_real_type(::Type{T}) where {T<:Integer} = Rational{BigInt}
_upb_exact_real_type(::Type{<:Rational}) = Rational{BigInt}
_upb_exact_work_type(::Type{T}) where {T<:Integer} = Rational{BigInt}
_upb_exact_work_type(::Type{<:Rational}) = Rational{BigInt}
function _upb_exact_work_type(::Type{Complex{T}}) where {T}
    return Complex{_upb_exact_real_type(T)}
end

_upb_is_exact_type(::Type{T}) where {T<:Union{Integer,Rational}} = true
_upb_is_exact_type(::Type{Complex{T}}) where {T} = _upb_is_exact_type(T)
_upb_is_exact_type(::Type) = false

_upb_real_type(::Type{Complex{T}}) where {T} = T
_upb_real_type(::Type{T}) where {T} = T

function _upb_supported_scalar_type(::Type{T}) where {T}
    T === Bool && return false
    T <: Union{Integer,Rational,AbstractFloat} && return isconcretetype(T)
    T <: Complex || return false
    R = _upb_real_type(T)
    R === Bool && return false
    return R <: Union{Integer,Rational,AbstractFloat} && isconcretetype(T)
end

function _upb_options(
    ::Type{W};
    atol=0,
    rtol=nothing,
    boundary_factor=8,
    max_partitions=_UPB_DEFAULT_MAX_PARTITIONS,
    max_work=_UPB_DEFAULT_MAX_WORK,
    max_dense_entries=_UPB_DEFAULT_MAX_DENSE_ENTRIES,
    max_states=_UPB_DEFAULT_MAX_STATES,
    normalization=:require,
    allow_densify=false,
    materialize_witness=false,
) where {W}
    _upb_supported_scalar_type(W) || throw(
        ArgumentError(
            "UPB inputs must have a concrete integer, rational, floating, or " *
            "complex counterpart element type; got $W",
        ),
    )
    normalization in (:require, :allow) ||
        throw(ArgumentError("normalization must be :require or :allow; got $normalization"))
    allow_densify isa Bool || throw(ArgumentError("allow_densify must be Bool"))
    materialize_witness isa Bool || throw(ArgumentError("materialize_witness must be Bool"))

    exact = _upb_is_exact_type(W)
    R = exact ? _upb_real_type(_upb_exact_work_type(W)) : _upb_real_type(W)
    checked_atol = try
        convert(R, atol)
    catch
        throw(ArgumentError("atol must be representable as $R"))
    end
    checked_rtol = if rtol === nothing
        exact ? zero(R) : convert(R, 8) * eps(R)
    else
        try
            convert(R, rtol)
        catch
            throw(ArgumentError("rtol must be representable as $R"))
        end
    end
    checked_boundary = try
        convert(R, boundary_factor)
    catch
        throw(ArgumentError("boundary_factor must be representable as $R"))
    end
    isfinite(checked_atol) && checked_atol >= zero(R) ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    isfinite(checked_rtol) && checked_rtol >= zero(R) ||
        throw(ArgumentError("rtol must be finite and nonnegative"))
    isfinite(checked_boundary) && checked_boundary > one(R) ||
        throw(ArgumentError("boundary_factor must be finite and greater than one"))
    exact &&
        (!iszero(checked_atol) || !iszero(checked_rtol)) &&
        throw(
            ArgumentError(
                "exact integer/rational inputs use exact decisions and require zero " *
                "atol and rtol",
            ),
        )

    return _UPBOptions(
        checked_atol,
        checked_rtol,
        checked_boundary,
        _upb_limit(max_partitions, "max_partitions"),
        _upb_limit(max_work, "max_work"),
        _upb_limit(max_dense_entries, "max_dense_entries"),
        _upb_limit(max_states, "max_states"),
        normalization,
        allow_densify,
        materialize_witness,
        exact,
    )
end

function _upb_threshold(options::_UPBOptions, scale::Int)
    return options.atol + options.rtol * max(scale, 1)
end

function _upb_boundary_threshold(options::_UPBOptions, scale::Int)
    return options.boundary_factor * _upb_threshold(options, scale)
end

function _upb_classify_zero(value, options::_UPBOptions, scale::Int)
    iszero(value) && return :zero
    lower = _upb_threshold(options, scale)
    value < lower && return :zero
    value <= _upb_boundary_threshold(options, scale) && return :boundary
    return :nonzero
end

function _upb_make_result(
    metadata;
    status,
    reason,
    is_upb_value,
    certificate_kind=:none,
    witness_factors=nothing,
    witness_vector=nothing,
    witness_partition=nothing,
    orthogonality_residual=nothing,
    product_residual=nothing,
    witness_residual=nothing,
    partitions_examined=0,
    work_used=0,
    offending_state=nothing,
    message,
)
    options = metadata.options
    recorded_product_residual =
        if product_residual === nothing && hasproperty(metadata, :product_residual)
            metadata.product_residual
        else
            product_residual
        end
    return UPBAnalysisResult(
        status,
        reason,
        is_upb_value,
        witness_factors,
        witness_vector,
        witness_partition,
        metadata.dimensions,
        metadata.state_count,
        metadata.input_form,
        certificate_kind,
        orthogonality_residual,
        recorded_product_residual,
        witness_residual,
        options.atol,
        options.rtol,
        options.boundary_factor,
        partitions_examined,
        work_used,
        options.max_partitions,
        options.max_work,
        options.max_states,
        options.normalization,
        metadata.input_normalized,
        metadata.analysis_rescaled,
        metadata.densified,
        offending_state,
        message,
    )
end

function _upb_dimensions(dims)
    dims isa Tuple ||
        dims isa AbstractVector ||
        throw(ArgumentError("dims must be a tuple or vector of local dimensions"))
    dims isa AbstractVector && Base.require_one_based_indexing(dims)
    length(dims) >= 2 || throw(ArgumentError("is_upb requires at least two tensor factors"))
    dimensions = ntuple(length(dims)) do index
        value = dims[index]
        value isa Bool && throw(ArgumentError("dims[$index] must be an integer at least 2"))
        value isa Integer ||
            throw(ArgumentError("dims[$index] must be an integer at least 2"))
        2 <= value <= typemax(Int) ||
            throw(ArgumentError("dims[$index] must be an integer at least 2"))
        return Int(value)
    end
    _upb_checked_product(dimensions, "the global Hilbert-space dimension")
    return dimensions
end

function _upb_matrix_entry_count(matrices)
    count = 0
    for matrix in matrices
        count = _upb_checked_add(
            count,
            _upb_checked_mul(
                size(matrix, 1), size(matrix, 2), "the UPB matrix entry count"
            ),
            "the UPB matrix entry count",
        )
    end
    return count
end

function _upb_check_dense_budget(required::Int, options::_UPBOptions, action)
    limit = options.max_dense_entries
    limit !== nothing &&
        required > limit &&
        throw(
            ArgumentError(
                "$action requires $required dense entries, exceeding " *
                "max_dense_entries=$limit",
            ),
        )
    return nothing
end

function _upb_isfinite_matrix(matrix)
    values = issparse(matrix) ? nonzeros(matrix) : matrix
    all(isfinite, values) ||
        throw(ArgumentError("UPB inputs must contain only finite values"))
    return nothing
end

function _upb_convert_local_factors(local_factors, options::_UPBOptions)
    length(local_factors) >= 2 ||
        throw(ArgumentError("is_upb requires at least two local-factor matrices"))
    all(factor -> factor isa AbstractMatrix, local_factors) || throw(
        ArgumentError(
            "local-factor input must contain only matrices, one matrix per subsystem"
        ),
    )
    for factor in local_factors
        Base.require_one_based_indexing(factor)
    end
    state_count = size(local_factors[1], 2)
    dimensions = ntuple(length(local_factors)) do party
        factor = local_factors[party]
        size(factor, 2) == state_count || throw(
            DimensionMismatch(
                "all local-factor matrices must have the same number of columns; " *
                "factor 1 has $state_count and factor $party has $(size(factor, 2))",
            ),
        )
        size(factor, 1) >= 2 || throw(
            ArgumentError(
                "local factor $party has dimension $(size(factor, 1)); remove " *
                "one-dimensional tensor factors before calling is_upb",
            ),
        )
        return size(factor, 1)
    end
    _upb_checked_product(dimensions, "the global Hilbert-space dimension")
    scalar_types = map(eltype, local_factors)
    all(_upb_supported_scalar_type, scalar_types) || throw(
        ArgumentError(
            "every local-factor matrix must have a supported concrete numeric " *
            "element type; got $(Tuple(scalar_types))",
        ),
    )
    promoted = promote_type(scalar_types...)
    _upb_supported_scalar_type(promoted) || throw(
        ArgumentError("the promoted local-factor element type $promoted is unsupported")
    )
    sparse_input = any(issparse, local_factors)
    sparse_input &&
        !options.allow_densify &&
        throw(
            ArgumentError(
                "sparse local factors require allow_densify=true because rank and " *
                "orthogonal-complement analysis uses dense work vectors",
            ),
        )
    dense_entries = _upb_matrix_entry_count(local_factors)
    _upb_check_dense_budget(dense_entries, options, "local-factor analysis")
    foreach(_upb_isfinite_matrix, local_factors)

    W = options.exact ? _upb_exact_work_type(promoted) : promoted
    converted = ntuple(party -> Matrix{W}(local_factors[party]), length(local_factors))
    return converted, dimensions, state_count, sparse_input
end

function _upb_column_norm_squared(column)
    value = real(dot(column, column))
    value > zero(value) ||
        throw(ArgumentError("every supplied local and global state must be nonzero"))
    isfinite(value) || throw(ArgumentError("state norms must be finite"))
    return value
end

function _upb_normalize_local_factors(factors, options::_UPBOptions)
    state_count = size(factors[1], 2)
    party_count = length(factors)
    input_normalized = true
    normalization_boundary = false
    analysis_rescaled = false

    if options.exact
        for state in 1:state_count
            global_norm_squared = prod(
                _upb_column_norm_squared(@view factors[party][:, state]) for
                party in 1:party_count
            )
            normalized = global_norm_squared == one(global_norm_squared)
            input_normalized &= normalized
            if options.normalization === :require && !normalized
                throw(
                    ArgumentError(
                        "product state $state has squared norm " *
                        "$global_norm_squared; normalization=:require never " *
                        "rescales caller input. Pass normalization=:allow to " *
                        "accept arbitrary nonzero scaling.",
                    ),
                )
            end
        end
        return factors, input_normalized, false, false
    end

    R = _upb_real_type(eltype(factors[1]))
    norms = [Vector{R}(undef, state_count) for _ in 1:party_count]
    for party in 1:party_count, state in 1:state_count
        norms[party][state] = sqrt(_upb_column_norm_squared(@view factors[party][:, state]))
    end
    for state in 1:state_count
        global_norm = prod(norms[party][state] for party in 1:party_count)
        isfinite(global_norm) || throw(
            ArgumentError(
                "the product norm of state $state is not finite; rescale the " *
                "local factors explicitly",
            ),
        )
        delta = abs(global_norm - one(global_norm))
        normalized_class = _upb_classify_zero(
            delta, options, sum(size(factors[party], 1) for party in 1:party_count)
        )
        input_normalized &= normalized_class === :zero
        normalization_boundary |= normalized_class === :boundary
        if options.normalization === :require && normalized_class === :nonzero
            throw(
                ArgumentError(
                    "product state $state has norm $global_norm; " *
                    "normalization=:require never rescales caller input. Pass " *
                    "normalization=:allow to accept arbitrary nonzero scaling.",
                ),
            )
        end
    end

    normalized = ntuple(party_count) do party
        matrix = copy(factors[party])
        for state in 1:state_count
            matrix[:, state] ./= norms[party][state]
            analysis_rescaled |= norms[party][state] != one(norms[party][state])
        end
        return matrix
    end
    return normalized, input_normalized, normalization_boundary, analysis_rescaled
end

function _upb_factor_exact_product(vector, dimensions)
    factors = Vector{Vector{eltype(vector)}}()
    current = copy(vector)
    for party in 1:(length(dimensions) - 1)
        local_dimension = dimensions[party]
        remaining_dimension = div(length(current), local_dimension)
        unfolding = reshape(current, remaining_dimension, local_dimension)
        pivot = findfirst(!iszero, unfolding)
        pivot === nothing && return :not_product, nothing, nothing
        pivot_row, pivot_column = Tuple(pivot)
        remaining_factor = copy(@view unfolding[:, pivot_column])
        local_factor = [
            unfolding[pivot_row, column] / unfolding[pivot_row, pivot_column] for
            column in 1:local_dimension
        ]
        remaining_factor * transpose(local_factor) == unfolding ||
            return :not_product, nothing, nothing
        push!(factors, local_factor)
        current = remaining_factor
    end
    push!(factors, copy(current))
    return :product, Tuple(factors), zero(_upb_real_type(eltype(vector)))
end

function _upb_factor_floating_product(vector, dimensions, options::_UPBOptions)
    factors = Vector{Vector{eltype(vector)}}()
    current = copy(vector)
    maximum_residual = zero(_upb_real_type(eltype(vector)))
    boundary_seen = false
    for party in 1:(length(dimensions) - 1)
        local_dimension = dimensions[party]
        remaining_dimension = div(length(current), local_dimension)
        unfolding = reshape(current, remaining_dimension, local_dimension)
        pivot = first(CartesianIndices(unfolding))
        pivot_magnitude = abs(unfolding[pivot])
        for index in CartesianIndices(unfolding)
            magnitude = abs(unfolding[index])
            if magnitude > pivot_magnitude
                pivot = index
                pivot_magnitude = magnitude
            end
        end
        iszero(pivot_magnitude) && return :not_product, nothing, maximum_residual
        pivot_row, pivot_column = Tuple(pivot)
        remaining_factor = copy(@view unfolding[:, pivot_column])
        local_factor = [
            unfolding[pivot_row, column] / unfolding[pivot_row, pivot_column] for
            column in 1:local_dimension
        ]
        reconstructed = remaining_factor * transpose(local_factor)
        residual = norm(unfolding - reconstructed) / norm(unfolding)
        maximum_residual = max(maximum_residual, residual)
        classification = _upb_classify_zero(residual, options, length(unfolding))
        classification === :nonzero && return :not_product, nothing, maximum_residual
        boundary_seen |= classification === :boundary
        push!(factors, local_factor)
        current = remaining_factor
    end
    push!(factors, copy(current))
    return boundary_seen ? :boundary : :product, Tuple(factors), maximum_residual
end

function _upb_convert_global_vectors(global_vectors, dims, option_keywords)
    dimensions = _upb_dimensions(dims)
    total_dimension = _upb_checked_product(dimensions, "the global Hilbert-space dimension")
    matrix = if global_vectors isa AbstractVector
        Base.require_one_based_indexing(global_vectors)
        reshape(global_vectors, :, 1)
    elseif global_vectors isa AbstractMatrix
        Base.require_one_based_indexing(global_vectors)
        global_vectors
    else
        throw(ArgumentError("global product states must be a vector or matrix"))
    end
    size(matrix, 1) == total_dimension || throw(
        DimensionMismatch(
            "global vectors have length $(size(matrix, 1)), but dims=$dimensions " *
            "have product $total_dimension",
        ),
    )
    scalar_type = eltype(matrix)
    _upb_supported_scalar_type(scalar_type) || throw(
        ArgumentError(
            "global vectors must have a supported concrete numeric element type; " *
            "got $scalar_type",
        ),
    )
    options = _upb_options(scalar_type; option_keywords...)
    sparse_input = issparse(global_vectors)
    sparse_input &&
        !options.allow_densify &&
        throw(
            ArgumentError(
                "sparse global vectors require allow_densify=true because product " *
                "factorization uses dense tensor unfoldings",
            ),
        )
    factor_entries = _upb_checked_mul(
        sum(dimensions), size(matrix, 2), "the global-factor entry count"
    )
    dense_entries = _upb_checked_add(
        _upb_checked_mul(size(matrix, 1), size(matrix, 2), "the global-vector entry count"),
        factor_entries,
        "the global-vector analysis entry count",
    )
    _upb_check_dense_budget(dense_entries, options, "global-vector analysis")
    _upb_isfinite_matrix(matrix)
    state_count = size(matrix, 2)
    if options.max_states !== nothing && state_count > options.max_states
        return (;
            early_status=:state_limit,
            offending_state=nothing,
            product_residual=nothing,
            factors=nothing,
            dimensions,
            state_count,
            options,
            densified=sparse_input,
        )
    end
    W = options.exact ? _upb_exact_work_type(scalar_type) : scalar_type
    dense = Matrix{W}(matrix)
    local_factors = ntuple(
        party -> Matrix{W}(undef, dimensions[party], state_count), length(dimensions)
    )
    maximum_product_residual = zero(_upb_real_type(W))
    for state in 1:state_count
        vector = @view dense[:, state]
        _upb_column_norm_squared(vector)
        factor_status, factors, residual = if options.exact
            _upb_factor_exact_product(vector, dimensions)
        else
            _upb_factor_floating_product(vector, dimensions, options)
        end
        residual !== nothing &&
            (maximum_product_residual = max(maximum_product_residual, residual))
        if factor_status !== :product
            return (;
                early_status=factor_status,
                offending_state=state,
                product_residual=maximum_product_residual,
                factors=nothing,
                dimensions,
                state_count,
                options,
                densified=sparse_input,
            )
        end
        for party in eachindex(factors)
            local_factors[party][:, state] .= factors[party]
        end
    end
    return (;
        early_status=:ready,
        offending_state=nothing,
        product_residual=maximum_product_residual,
        factors=local_factors,
        dimensions,
        state_count,
        options,
        densified=sparse_input,
    )
end

function _upb_orthogonality(factors, options::_UPBOptions)
    state_count = size(factors[1], 2)
    maximum_residual = zero(_upb_real_type(eltype(factors[1])))
    boundary_pair = nothing
    for first in 1:(state_count - 1), second in (first + 1):state_count
        overlap = prod(
            dot(@view(factors[party][:, first]), @view(factors[party][:, second])) for
            party in eachindex(factors)
        )
        residual = abs(overlap)
        maximum_residual = max(maximum_residual, residual)
        if options.exact
            iszero(overlap) || return :not_orthogonal, maximum_residual, (first, second)
        else
            classification = _upb_classify_zero(
                residual, options, sum(size(factor, 1) for factor in factors)
            )
            classification === :nonzero &&
                return :not_orthogonal, maximum_residual, (first, second)
            classification === :boundary &&
                boundary_pair === nothing &&
                (boundary_pair = (first, second))
        end
    end
    boundary_pair === nothing || return :boundary, maximum_residual, boundary_pair
    return :orthogonal, maximum_residual, nothing
end

function _upb_rref!(matrix)
    rows, columns = size(matrix)
    pivot_columns = Int[]
    pivot_row = 1
    for column in 1:columns
        selected = findfirst(row -> !iszero(matrix[row, column]), pivot_row:rows)
        selected === nothing && continue
        selected_row = pivot_row + selected - 1
        if selected_row != pivot_row
            temporary = copy(@view matrix[pivot_row, :])
            matrix[pivot_row, :] .= @view matrix[selected_row, :]
            matrix[selected_row, :] .= temporary
        end
        pivot_value = matrix[pivot_row, column]
        matrix[pivot_row, :] ./= pivot_value
        for row in 1:rows
            row == pivot_row && continue
            factor = matrix[row, column]
            iszero(factor) && continue
            matrix[row, :] .-= factor .* @view(matrix[pivot_row, :])
        end
        push!(pivot_columns, column)
        pivot_row += 1
        pivot_row > rows && break
    end
    return pivot_columns
end

function _upb_exact_span(matrix, indices)
    dimension = size(matrix, 1)
    isempty(indices) && begin
        witness = zeros(eltype(matrix), dimension)
        witness[1] = one(eltype(matrix))
        return :deficient, witness
    end
    equations = Matrix(adjoint(matrix[:, indices]))
    pivot_columns = _upb_rref!(equations)
    length(pivot_columns) == dimension && return :full, nothing
    free_column = findfirst(column -> column ∉ pivot_columns, 1:dimension)
    free_column === nothing && error("internal UPB nullspace construction failed")
    witness = zeros(eltype(matrix), dimension)
    witness[free_column] = one(eltype(matrix))
    for (row, pivot_column) in pairs(pivot_columns)
        witness[pivot_column] = -equations[row, free_column]
    end
    all(iszero, adjoint(matrix[:, indices]) * witness) ||
        error("internal exact UPB witness failed its orthogonality check")
    return :deficient, witness
end

function _upb_project_residual(vector, basis)
    residual = copy(vector)
    for _ in 1:2
        for direction in basis
            residual .-= direction .* dot(direction, residual)
        end
    end
    return residual
end

function _upb_float_span(matrix, indices, options::_UPBOptions)
    dimension = size(matrix, 1)
    isempty(indices) && begin
        witness = zeros(eltype(matrix), dimension)
        witness[1] = one(eltype(matrix))
        return :deficient, witness
    end
    remaining = [copy(@view matrix[:, index]) for index in indices]
    basis = Vector{Vector{eltype(matrix)}}()
    boundary_seen = false
    while !isempty(remaining) && length(basis) < dimension
        best_position = 1
        best_residual = _upb_project_residual(remaining[1], basis)
        best_norm = norm(best_residual)
        for position in 2:length(remaining)
            residual = _upb_project_residual(remaining[position], basis)
            residual_norm = norm(residual)
            if residual_norm > best_norm
                best_position = position
                best_residual = residual
                best_norm = residual_norm
            end
        end
        classification = _upb_classify_zero(
            best_norm, options, max(dimension, length(indices))
        )
        classification === :zero && break
        boundary_seen |= classification === :boundary
        best_residual ./= best_norm
        push!(basis, best_residual)
        deleteat!(remaining, best_position)
    end
    length(basis) == dimension &&
        return boundary_seen ? (:boundary, nothing) : (:full, nothing)
    boundary_seen && return :boundary, nothing

    best_witness = nothing
    best_norm = zero(_upb_real_type(eltype(matrix)))
    for coordinate in 1:dimension
        candidate = zeros(eltype(matrix), dimension)
        candidate[coordinate] = one(eltype(matrix))
        residual = _upb_project_residual(candidate, basis)
        residual_norm = norm(residual)
        if residual_norm > best_norm
            best_witness = residual
            best_norm = residual_norm
        end
    end
    best_witness === nothing &&
        error("internal floating UPB complement construction failed")
    complement_class = _upb_classify_zero(best_norm, options, dimension)
    complement_class === :zero &&
        error("internal floating UPB complement unexpectedly vanished")
    complement_class === :boundary && return :boundary, nothing
    best_witness ./= best_norm
    witness_residual = maximum(
        index -> abs(dot(best_witness, @view matrix[:, index])),
        indices;
        init=zero(best_norm),
    )
    witness_class = _upb_classify_zero(
        witness_residual, options, max(dimension, length(indices))
    )
    witness_class === :zero && return :deficient, best_witness
    return :boundary, nothing
end

function _upb_rank_work(dimension::Int, column_count::Int)
    common = min(dimension, column_count)
    return try
        Base.checked_mul(
            4,
            Base.checked_mul(
                dimension, Base.checked_mul(max(column_count, 1), max(common, 1))
            ),
        )
    catch err
        err isa OverflowError || rethrow()
        typemax(Int)
    end
end

function _upb_reserve_work!(state::_UPBSearchState, amount::Int, options::_UPBOptions)
    new_work = try
        Base.checked_add(state.work_used, amount)
    catch err
        err isa OverflowError || rethrow()
        state.stopped_reason = :work_limit
        return false
    end
    if options.max_work !== nothing && new_work > options.max_work
        state.stopped_reason = :work_limit
        return false
    end
    state.work_used = new_work
    return true
end

function _upb_witness_residual(factors, witnesses, options::_UPBOptions)
    state_count = size(factors[1], 2)
    maximum_residual = zero(_upb_real_type(eltype(factors[1])))
    for state in 1:state_count
        overlap = prod(
            dot(witnesses[party], @view(factors[party][:, state])) for
            party in eachindex(factors)
        )
        residual = abs(overlap)
        maximum_residual = max(maximum_residual, residual)
        if options.exact
            iszero(overlap) || return :invalid, maximum_residual
        end
    end
    options.exact && return :valid, maximum_residual
    classification = _upb_classify_zero(
        maximum_residual, options, sum(size(factor, 1) for factor in factors)
    )
    classification === :zero && return :valid, maximum_residual
    return :boundary, maximum_residual
end

function _upb_materialize_witness(witnesses, dimensions, options::_UPBOptions)
    options.materialize_witness || return nothing
    total_dimension = _upb_checked_product(dimensions, "the product-witness dimension")
    _upb_check_dense_budget(
        total_dimension, options, "materializing the product-vector witness"
    )
    return foldl(kron, witnesses)
end

function _upb_evaluate_partition!(state, factors, groups, options)
    if options.max_partitions !== nothing &&
        state.partitions_examined >= options.max_partitions
        state.stopped_reason = :partition_limit
        return nothing
    end
    state.partitions_examined += 1
    witnesses = Vector{Vector{eltype(factors[1])}}(undef, length(factors))
    for party in eachindex(factors)
        rank_work = _upb_rank_work(size(factors[party], 1), length(groups[party]))
        _upb_reserve_work!(state, rank_work, options) || return nothing
        span_status, witness = if options.exact
            _upb_exact_span(factors[party], groups[party])
        else
            _upb_float_span(factors[party], groups[party], options)
        end
        span_status === :full && return nothing
        if span_status === :boundary
            state.boundary_seen = true
            return nothing
        end
        witnesses[party] = witness
    end
    witness_tuple = Tuple(witnesses)
    validity, residual = _upb_witness_residual(factors, witness_tuple, options)
    if validity === :valid
        partition = Tuple(Tuple(group) for group in groups)
        return witness_tuple, partition, residual
    end
    state.boundary_seen = true
    return nothing
end

function _upb_search_partitions!(
    state, factors, minimum_sizes, groups, counts, state_index, state_count, options
)
    state.stopped_reason === nothing || return nothing
    if state_index > state_count
        all(counts .>= minimum_sizes) || return nothing
        return _upb_evaluate_partition!(state, factors, groups, options)
    end
    remaining_after = state_count - state_index
    for party in eachindex(groups)
        counts[party] += 1
        push!(groups[party], state_index)
        feasible = true
        for check_party in eachindex(groups)
            if counts[check_party] + remaining_after < minimum_sizes[check_party]
                feasible = false
                break
            end
        end
        result = if feasible
            _upb_search_partitions!(
                state,
                factors,
                minimum_sizes,
                groups,
                counts,
                state_index + 1,
                state_count,
                options,
            )
        else
            nothing
        end
        pop!(groups[party])
        counts[party] -= 1
        result === nothing || return result
        state.stopped_reason === nothing || return nothing
    end
    return nothing
end

function _upb_small_family_witness(factors, dimensions, state_count, options)
    groups = [Int[] for _ in dimensions]
    next_state = 1
    for party in eachindex(dimensions)
        capacity = dimensions[party] - 1
        while length(groups[party]) < capacity && next_state <= state_count
            push!(groups[party], next_state)
            next_state += 1
        end
    end
    next_state == state_count + 1 ||
        error("internal small-family UPB assignment did not cover every state")
    state = _UPBSearchState(0, 0, false, nothing)
    result = _upb_evaluate_partition!(state, factors, groups, options)
    return result, state
end

function _upb_analyze_local(
    factors, dimensions, state_count, input_form, densified, options, product_residual
)
    analyzed_factors, input_normalized, normalization_boundary, analysis_rescaled = _upb_normalize_local_factors(
        factors, options
    )
    metadata = (;
        dimensions,
        state_count,
        input_form,
        options,
        input_normalized,
        analysis_rescaled,
        densified,
        product_residual,
    )
    if options.normalization === :require && normalization_boundary
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:numerical_boundary,
            is_upb_value=nothing,
            offending_state=nothing,
            message="normalization lies inside the declared numerical boundary band",
        )
    end

    orthogonality_status, orthogonality_residual, offending_pair = _upb_orthogonality(
        analyzed_factors, options
    )
    if orthogonality_status === :not_orthogonal
        return _upb_make_result(
            metadata;
            status=:not_upb,
            reason=:not_orthogonal,
            is_upb_value=false,
            certificate_kind=options.exact ? :exact : :tolerance_robust,
            orthogonality_residual,
            offending_state=offending_pair,
            message="the supplied product states are not mutually orthogonal",
        )
    elseif orthogonality_status === :boundary
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:numerical_boundary,
            is_upb_value=nothing,
            orthogonality_residual,
            offending_state=offending_pair,
            message="a pairwise orthogonality check lies inside the boundary band",
        )
    end

    total_dimension = _upb_checked_product(dimensions, "the global Hilbert-space dimension")
    if state_count >= total_dimension
        reason = state_count == total_dimension ? :complete_basis : :too_many_states
        return _upb_make_result(
            metadata;
            status=:not_upb,
            reason,
            is_upb_value=false,
            certificate_kind=options.exact ? :exact : :tolerance_robust,
            orthogonality_residual,
            message=if state_count == total_dimension
                "a UPB must be incomplete; these states form a complete basis"
            else
                "an orthogonal family cannot contain more states than the ambient dimension"
            end,
        )
    end

    minimum_sizes = [dimension - 1 for dimension in dimensions]
    minimum_total = sum(minimum_sizes)
    search_result, search_state = if state_count < minimum_total
        _upb_small_family_witness(analyzed_factors, dimensions, state_count, options)
    else
        groups = [Int[] for _ in dimensions]
        counts = zeros(Int, length(dimensions))
        state = _UPBSearchState(0, 0, false, nothing)
        result = _upb_search_partitions!(
            state,
            analyzed_factors,
            minimum_sizes,
            groups,
            counts,
            1,
            state_count,
            options,
        )
        result, state
    end

    if search_result !== nothing
        witness_factors, witness_partition, witness_residual = search_result
        witness_vector = _upb_materialize_witness(witness_factors, dimensions, options)
        return _upb_make_result(
            metadata;
            status=:not_upb,
            reason=:extension_witness,
            is_upb_value=false,
            certificate_kind=options.exact ? :exact : :tolerance_robust,
            witness_factors,
            witness_vector,
            witness_partition,
            orthogonality_residual,
            witness_residual,
            partitions_examined=search_state.partitions_examined,
            work_used=search_state.work_used,
            message="an explicit product vector extends the supplied orthogonal family",
        )
    elseif search_state.stopped_reason !== nothing
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=search_state.stopped_reason,
            is_upb_value=nothing,
            orthogonality_residual,
            partitions_examined=search_state.partitions_examined,
            work_used=search_state.work_used,
            message=if search_state.stopped_reason === :partition_limit
                "the exhaustive search reached max_partitions"
            else
                "the exhaustive search reached max_work"
            end,
        )
    elseif search_state.boundary_seen
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:numerical_boundary,
            is_upb_value=nothing,
            orthogonality_residual,
            partitions_examined=search_state.partitions_examined,
            work_used=search_state.work_used,
            message="at least one local-span decision lies inside the boundary band",
        )
    end

    return _upb_make_result(
        metadata;
        status=:upb,
        reason=:unextendible,
        is_upb_value=true,
        certificate_kind=options.exact ? :exact : :tolerance_robust,
        orthogonality_residual,
        partitions_examined=search_state.partitions_examined,
        work_used=search_state.work_used,
        message="all admissible local-span partitions were exhausted without an extension",
    )
end

"""
    is_upb(
        local_factors;
        normalization=:require,
        atol=0,
        rtol=nothing,
        boundary_factor=8,
        max_partitions=100_000,
        max_work=25_000_000,
        max_dense_entries=1_000_000,
        max_states=256,
        allow_densify=false,
        materialize_witness=false,
    ) -> UPBAnalysisResult

Determine whether local-factor columns form an unextendible product basis.

`local_factors` is a tuple or vector of matrices. Matrix `j` has size
`dⱼ × s`, and its `k`th column is the local factor of product state `k`.
At least two local dimensions are required. The routine verifies nonzero,
finite, one-based input, equal column counts, normalization according to the
explicit policy, mutual orthogonality, incompleteness, and unextendibility.

The default `normalization=:require` rejects clearly non-unit global product
columns and never modifies caller input. `normalization=:allow` explicitly
accepts arbitrary nonzero scaling. Floating-point analysis uses normalized
private work copies for scale-independent rank decisions, reported through
`analysis_rescaled`; exact integer/rational analysis performs no normalization
and uses exact Gaussian elimination.

Unextendibility is decided by a lazy, ordered version of the finite partition
criterion: every supplied state is assigned to a subsystem, and an extension
exists exactly when every assigned local family has rank below its local
dimension. `max_states` guards input-dependent pair checks and recursion depth;
`max_partitions` and `max_work` stop the combinatorial search. Each produces
`status == :unknown` and may be disabled explicitly with `nothing`.

Sparse inputs are rejected unless `allow_densify=true`; the dense workspace is
then guarded by `max_dense_entries` and reported in the result. No input is
silently normalized, repaired, conjugated, or projected.
"""
function is_upb(
    local_factors::Tuple;
    atol=0,
    rtol=nothing,
    boundary_factor=8,
    max_partitions=_UPB_DEFAULT_MAX_PARTITIONS,
    max_work=_UPB_DEFAULT_MAX_WORK,
    max_dense_entries=_UPB_DEFAULT_MAX_DENSE_ENTRIES,
    max_states=_UPB_DEFAULT_MAX_STATES,
    normalization=:require,
    allow_densify=false,
    materialize_witness=false,
)
    isempty(local_factors) &&
        throw(ArgumentError("local_factors must contain at least two matrices"))
    scalar_types = map(
        factor -> factor isa AbstractMatrix ? eltype(factor) : Bool, local_factors
    )
    all(factor -> factor isa AbstractMatrix, local_factors) ||
        throw(ArgumentError("local_factors must contain only matrices"))
    promoted = promote_type(scalar_types...)
    options = _upb_options(
        promoted;
        atol,
        rtol,
        boundary_factor,
        max_partitions,
        max_work,
        max_dense_entries,
        max_states,
        normalization,
        allow_densify,
        materialize_witness,
    )
    factors, dimensions, state_count, densified = _upb_convert_local_factors(
        local_factors, options
    )
    if options.max_states !== nothing && state_count > options.max_states
        metadata = (;
            dimensions,
            state_count,
            input_form=:local_factors,
            options,
            input_normalized=nothing,
            analysis_rescaled=false,
            densified,
            product_residual=nothing,
        )
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:state_limit,
            is_upb_value=nothing,
            message="the input has $state_count states, exceeding " *
                    "max_states=$(options.max_states)",
        )
    end
    return _upb_analyze_local(
        factors, dimensions, state_count, :local_factors, densified, options, nothing
    )
end

function is_upb(local_factors::AbstractVector; kwargs...)
    Base.require_one_based_indexing(local_factors)
    all(factor -> factor isa AbstractMatrix, local_factors) || throw(
        ArgumentError(
            "a one-argument vector input must contain local-factor matrices; " *
            "pass global vectors together with dims",
        ),
    )
    return is_upb(Tuple(local_factors); kwargs...)
end

function is_upb(
    first_factor::AbstractMatrix,
    second_factor::AbstractMatrix,
    remaining_factors::AbstractMatrix...;
    kwargs...,
)
    return is_upb((first_factor, second_factor, remaining_factors...); kwargs...)
end

"""
    is_upb(global_vectors, dims; kwargs...) -> UPBAnalysisResult

Analyze global product-vector columns in the tensor layout `dims`.

Each column is factorized deterministically in the package tensor convention
(subsystem 1 slowest). A column whose rank-one tensor residual is outside the
declared tolerance returns `status == :not_upb` and `reason == :not_product`;
a residual inside the boundary band returns `:unknown`. The resulting local
factors are then analyzed by the same certificate-aware algorithm as the
local-factor form.

All keywords and resource semantics are the same as for
`is_upb(local_factors)`. A single global vector is accepted as a one-column
family.
"""
function is_upb(
    global_vectors::Union{AbstractVector,AbstractMatrix},
    dims::Union{Tuple,AbstractVector};
    kwargs...,
)
    converted = _upb_convert_global_vectors(global_vectors, dims, kwargs)
    metadata = (;
        dimensions=converted.dimensions,
        state_count=converted.state_count,
        input_form=:global_vectors,
        options=converted.options,
        input_normalized=nothing,
        analysis_rescaled=false,
        densified=converted.densified,
        product_residual=converted.product_residual,
    )
    if converted.early_status === :state_limit
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:state_limit,
            is_upb_value=nothing,
            message="the input has $(converted.state_count) states, exceeding " *
                    "max_states=$(converted.options.max_states)",
        )
    elseif converted.early_status === :not_product
        return _upb_make_result(
            metadata;
            status=:not_upb,
            reason=:not_product,
            is_upb_value=false,
            certificate_kind=converted.options.exact ? :exact : :tolerance_robust,
            offending_state=converted.offending_state,
            message="a supplied global column is not a product vector",
        )
    elseif converted.early_status === :boundary
        return _upb_make_result(
            metadata;
            status=:unknown,
            reason=:numerical_boundary,
            is_upb_value=nothing,
            offending_state=converted.offending_state,
            message="a global product-factorization residual lies inside the boundary band",
        )
    end
    return _upb_analyze_local(
        converted.factors,
        converted.dimensions,
        converted.state_count,
        :global_vectors,
        converted.densified,
        converted.options,
        converted.product_residual,
    )
end
