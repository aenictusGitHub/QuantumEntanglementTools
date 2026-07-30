# Source-informed independent Julia implementation based on the specifications
# in QETLAB AbsPPTConstraints.m and IsAbsPPT.m at
# d8589610f00cff106537268dee2e2a1153f3a601 and on R. Hildebrand,
# "Positive partial transpose from spectra", Phys. Rev. A 76, 052325
# (2007), arXiv:quant-ph/0502170.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    AbsPPTEnumerationStatus

Completion state of a Hildebrand/QETLAB absolute-PPT LMI enumeration.
`AbsPPTEnumerationExhaustive` means that every ordering admitted by the
documented monotonicity and criss-cross rules was visited. A limit status is a
partial family and must not be treated as the full absolute-PPT criterion.
"""
@enum AbsPPTEnumerationStatus::UInt8 begin
    AbsPPTEnumerationExhaustive
    AbsPPTEnumerationConstraintLimit
    AbsPPTEnumerationWorkLimit
    AbsPPTEnumerationEntryLimit
    AbsPPTEnumerationEarlyViolation
end

"""
    AbsolutePPTStatus

Certificate-aware status returned by [`is_abs_ppt`](@ref).
"""
@enum AbsolutePPTStatus::UInt8 begin
    AbsolutePPTAnalyticCertified
    AbsolutePPTSufficientTestPassed
    AbsolutePPTExhaustiveCertified
    AbsolutePPTCertifiedNot
    AbsolutePPTCappedUnknown
    AbsolutePPTNumericalBoundary
    AbsolutePPTBackendUnavailable
    AbsolutePPTBackendFailure
end

const _AbsPPTPair = NTuple{2,Int}
const _AbsPPTReadOnlyPairVector = _ReadOnlyPlanArray{_AbsPPTPair,1,Vector{_AbsPPTPair}}

"""
    AbsPPTOrdering

One explicit ordering of the products `xᵢxⱼ` in Hildebrand's absolute-PPT
criterion. `positive_pairs[r] == (i,j)` means that `xᵢxⱼ` has rank `r`
(largest first) among pairs with `i <= j`. `negative_pairs` is the induced
relative ordering of the strictly off-diagonal pairs.

The pair arrays are owned, read-only copies. They expose QETLAB's permutation
and criss-cross semantics without storing a solver expression.
"""
struct AbsPPTOrdering
    positive_pairs::_AbsPPTReadOnlyPairVector
    negative_pairs::_AbsPPTReadOnlyPairVector
end

function AbsPPTOrdering(
    positive_pairs::AbstractVector{<:Tuple}, negative_pairs::AbstractVector{<:Tuple}
)
    Base.require_one_based_indexing(positive_pairs, negative_pairs)
    converted_positive = _AbsPPTPair[]
    converted_negative = _AbsPPTPair[]
    for (name, source, target) in (
        ("positive_pairs", positive_pairs, converted_positive),
        ("negative_pairs", negative_pairs, converted_negative),
    )
        for pair in source
            length(pair) == 2 || throw(ArgumentError("$name must contain two-index pairs"))
            left, right = pair
            left isa Integer && !(left isa Bool) ||
                throw(ArgumentError("$name indices must be integers"))
            right isa Integer && !(right isa Bool) ||
                throw(ArgumentError("$name indices must be integers"))
            left <= right || throw(ArgumentError("$name pairs must be stored with i <= j"))
            left >= 1 ||
                throw(ArgumentError("$name indices must be one-based and positive"))
            push!(target, (Int(left), Int(right)))
        end
    end
    isempty(converted_positive) &&
        !isempty(converted_negative) &&
        throw(ArgumentError("negative_pairs cannot be nonempty without positive_pairs"))
    if !isempty(converted_positive)
        p = maximum(last, converted_positive)
        expected_positive = p * (p + 1) ÷ 2
        expected_negative = p * (p - 1) ÷ 2
        length(converted_positive) == expected_positive || throw(
            ArgumentError(
                "positive_pairs has length $(length(converted_positive)); " *
                "expected $expected_positive for local dimension $p",
            ),
        )
        length(converted_negative) == expected_negative || throw(
            ArgumentError(
                "negative_pairs has length $(length(converted_negative)); " *
                "expected $expected_negative for local dimension $p",
            ),
        )
        Set(converted_positive) == Set((row, column) for row in 1:p for column in row:p) ||
            throw(ArgumentError("positive_pairs is not a permutation of i <= j pairs"))
        converted_negative ==
        filter(pair -> first(pair) != last(pair), converted_positive) || throw(
            ArgumentError(
                "negative_pairs must be the off-diagonal positive pairs in " *
                "the same relative order",
            ),
        )
    end
    return AbsPPTOrdering(
        _read_only_plan_array(converted_positive), _read_only_plan_array(converted_negative)
    )
end

"""
    AbsPPTConstraint

One numbered absolute-PPT LMI together with the product ordering that generated
it. `matrix` is either a concrete numeric Hermitian matrix or a
[`HermitianAffineMatrix`](@ref); it never contains JuMP or solver-owned values.
"""
struct AbsPPTConstraint{M}
    index::Int
    ordering::AbsPPTOrdering
    matrix::M
end

"""
    AbsPPTConstraintFamily

Typed, solver-independent output of [`abs_ppt_constraints`](@ref).

Besides the generated constraints, the family records whether enumeration was
exhaustive, combinatorial estimates, actual work and storage, explicit limits,
and the first robust PSD violation or numerical boundary when PSD checking was
requested.
"""
struct AbsPPTConstraintFamily{S,C,V,B,L}
    spectrum::S
    dimensions::NTuple{2,Int}
    local_dimension::Int
    constraints::C
    status::AbsPPTEnumerationStatus
    exhaustive::Bool
    monotone_ordering_count::BigInt
    known_criss_cross_count::Union{Nothing,BigInt}
    candidate_orderings_checked::BigInt
    work_used::BigInt
    entries_stored::BigInt
    limits::L
    first_violation::V
    first_boundary::B
    input_kind::Symbol
    symbolic::Bool
    message::String
end

function Base.show(io::IO, family::AbsPPTConstraintFamily)
    return print(
        io,
        "AbsPPTConstraintFamily(status=",
        family.status,
        ", constraints=",
        length(family.constraints),
        ", exhaustive=",
        family.exhaustive,
        ", dims=",
        family.dimensions,
        ")",
    )
end

"""
    AbsPPTOrderingCertificate

Exact rational realization of an ordering. `log_gaps[i]` represents
`log(xᵢ) - log(xᵢ₊₁) > 0`. Every adjacent product-order inequality is
rechecked with rational arithmetic before this value is constructed.
"""
struct AbsPPTOrderingCertificate{R<:Rational}
    log_gaps::Vector{R}
    log_coordinates::Vector{R}
    minimum_gap::R
    minimum_product_margin::R
    source::Symbol
end

"""
    IsAbsPPTResult

Status-rich absolute-PPT analysis. `verdict` is `true` or `false` only for a
mathematically certified branch; otherwise it is `nothing`. A negative result
retains the violating LMI/PSD diagnostic and an exact ordering-realization
certificate. Optional backend evidence is retained verbatim.
"""
struct IsAbsPPTResult{S,F,V,O,B,T,M,U}
    status::AbsolutePPTStatus
    verdict::Union{Nothing,Bool}
    certificate_kind::Union{Nothing,Symbol}
    spectrum::S
    dimensions::NTuple{2,Int}
    trace_value::T
    family::F
    violating_constraint::V
    ordering_certificate::O
    backend_result::B
    checked_constraints::BigInt
    planned_orderings::BigInt
    margin::M
    tolerance::U
    message::String
end

function Base.show(io::IO, result::IsAbsPPTResult)
    return print(
        io,
        "IsAbsPPTResult(status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", checked=",
        result.checked_constraints,
        "/",
        result.planned_orderings,
        ")",
    )
end

function _abs_ppt_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _abs_ppt_nonnegative_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _abs_ppt_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be a positive integer or nothing"))
    return BigInt(value)
end

function _abs_ppt_dimensions(dims, dimension::Int)
    if dims === nothing
        root = isqrt(dimension)
        root * root == dimension || throw(
            ArgumentError(
                "dims must be supplied when the spectrum length $dimension is not " *
                "a perfect square",
            ),
        )
        return (root, root)
    elseif dims isa Integer
        first_dimension = _abs_ppt_positive_integer(dims, "dims")
        dimension % first_dimension == 0 || throw(
            DimensionMismatch(
                "scalar dims=$first_dimension does not divide spectrum length $dimension",
            ),
        )
        return (first_dimension, dimension ÷ first_dimension)
    elseif dims isa Union{Tuple,AbstractVector}
        length(dims) == 2 ||
            throw(ArgumentError("dims must contain exactly two local dimensions"))
        first_dimension = _abs_ppt_positive_integer(dims[1], "dims[1]")
        second_dimension = _abs_ppt_positive_integer(dims[2], "dims[2]")
        BigInt(first_dimension) * BigInt(second_dimension) == dimension || throw(
            DimensionMismatch(
                "dims=($first_dimension, $second_dimension) has product " *
                "$(BigInt(first_dimension) * second_dimension), not spectrum " *
                "length $dimension",
            ),
        )
        return (first_dimension, second_dimension)
    end
    return throw(
        ArgumentError(
            "dims must be nothing, a positive integer, or a two-element tuple/vector"
        ),
    )
end

function _abs_ppt_real_component_type(::Type{T}) where {T<:Real}
    return T
end

function _abs_ppt_real_component_type(::Type{Complex{T}}) where {T<:Real}
    return T
end

function _abs_ppt_tolerance(real_type::Type; atol, rtol, scale)
    exact = real_type <: Integer || real_type <: Rational
    absolute = atol === nothing ? zero(real_type) : atol
    relative = if rtol === nothing
        exact ? zero(real_type) : sqrt(eps(real_type))
    else
        rtol
    end
    for (name, value) in (("atol", absolute), ("rtol", relative))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    exact &&
        (!iszero(absolute) || !iszero(relative)) &&
        throw(
            ArgumentError(
                "integer and rational spectra use exact comparisons; atol and rtol " *
                "must be zero",
            ),
        )
    promoted_absolute, promoted_relative, promoted_scale = promote(
        absolute, relative, scale
    )
    return (
        promoted_absolute + promoted_relative * max(one(promoted_scale), promoted_scale)
    )
end

function _abs_ppt_numeric_spectrum(spectrum::AbstractVector{<:Number}; allow_densify::Bool)
    Base.require_one_based_indexing(spectrum)
    isempty(spectrum) && throw(ArgumentError("the spectrum must be nonempty"))
    if SparseArrays.issparse(spectrum) && !allow_densify
        throw(
            ArgumentError(
                "sorting a sparse spectrum requires dense storage; pass " *
                "allow_densify=true after reviewing the length-$(length(spectrum)) " *
                "allocation",
            ),
        )
    end
    values = Real[]
    for (index, value) in pairs(spectrum)
        value isa Real && !(value isa Bool) || throw(
            ArgumentError(
                "spectrum entry $index must be real and non-Boolean; got $(repr(value))"
            ),
        )
        isfinite(value) || throw(ArgumentError("spectrum entry $index must be finite"))
        push!(values, value)
    end
    T = promote_type(map(typeof, values)...)
    T <: Real || throw(ArgumentError("spectrum entries must promote to a real type"))
    converted = T[value for value in values]
    sort!(converted; rev=true)
    return converted
end

function _abs_ppt_matrix_spectrum(
    matrix::AbstractMatrix{<:Number};
    allow_densify::Bool,
    max_entries,
    atol,
    rtol,
    boundary_allowed::Bool,
)
    Base.require_one_based_indexing(matrix)
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("the input matrix must be square; got $(size(matrix))"))
    dimension = size(matrix, 1)
    dimension > 0 || throw(ArgumentError("the input matrix must be nonempty"))
    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Number || throw(
        ArgumentError(
            "the input matrix must have a concrete numeric element type; got " *
            "$value_type",
        ),
    )
    for (index, value) in pairs(matrix)
        value isa Number && !(value isa Bool) ||
            throw(ArgumentError("matrix entry $index must be numeric and non-Boolean"))
        isfinite(value) || throw(ArgumentError("matrix entry $index must be finite"))
    end

    real_type = _abs_ppt_real_component_type(value_type)
    (real_type <: Integer || real_type <: Rational || real_type <: AbstractFloat) || throw(
        ArgumentError(
            "matrix real components must be integer, rational, or floating point; " *
            "got $real_type",
        ),
    )
    exact = real_type <: Integer || real_type <: Rational
    scale = maximum(abs, matrix; init=one(real_type))
    tolerance = _abs_ppt_tolerance(real_type; atol, rtol, scale)
    defect = maximum(abs, matrix - adjoint(matrix); init=zero(real_type))
    if exact
        iszero(defect) || throw(ArgumentError("the input matrix must be exactly Hermitian"))
    elseif defect > tolerance
        throw(
            ArgumentError(
                "the input matrix is non-Hermitian outside tolerance; maximum " *
                "defect $defect exceeds $tolerance",
            ),
        )
    elseif !iszero(defect)
        boundary_allowed || throw(
            ArgumentError(
                "the input matrix has a nonzero Hermiticity defect inside tolerance; " *
                "it is not symmetrized implicitly",
            ),
        )
        return (
            spectrum=nothing,
            boundary=true,
            tolerance=tolerance,
            message="the input Hermiticity defect lies inside the numerical boundary",
        )
    end

    if isdiag(matrix)
        diagonal = real.(diag(matrix))
        sort!(diagonal; rev=true)
        return (
            spectrum=diagonal,
            boundary=false,
            tolerance=tolerance,
            message="diagonal spectrum extracted without densification",
        )
    end

    dense_entries = BigInt(dimension)^2
    checked_entries = _abs_ppt_limit(max_entries, "max_entries")
    checked_entries !== nothing &&
        dense_entries > checked_entries &&
        throw(
            ArgumentError(
                "the eigendecomposition needs $dense_entries dense entries, exceeding " *
                "max_entries=$checked_entries",
            ),
        )
    if SparseArrays.issparse(matrix) && !allow_densify
        throw(
            ArgumentError(
                "the matrix eigendecomposition would densify sparse input; pass " *
                "allow_densify=true after reviewing the $(dimension)x$(dimension) " *
                "allocation",
            ),
        )
    end
    dense = SparseArrays.issparse(matrix) ? Matrix(matrix) : matrix
    eigenvalues = try
        eigvals(Hermitian(dense))
    catch error
        if real_type <: BigFloat
            throw(
                ArgumentError(
                    "the standard-library eigensolver does not support this " *
                    "BigFloat matrix; pass its real eigenvalue vector instead",
                ),
            )
        end
        rethrow(error)
    end
    sort!(eigenvalues; rev=true)
    return (
        spectrum=eigenvalues,
        boundary=false,
        tolerance=tolerance,
        message=if exact
            "a non-diagonal exact matrix was spectrally promoted by the " *
            "standard-library eigensolver"
        else
            "matrix eigenvalues computed without modifying the input"
        end,
    )
end

function _abs_ppt_monotone_ordering_count(p::Int)
    p <= 1 && return BigInt(1)
    q = p * (p + 1) ÷ 2
    numerator = factorial(BigInt(q))
    for index in 1:(p - 1)
        numerator *= factorial(BigInt(index))
    end
    denominator = BigInt(1)
    for index in 1:p
        denominator *= factorial(BigInt(2index - 1))
    end
    return numerator ÷ denominator
end

function _abs_ppt_known_criss_cross_count(p::Int)
    counts = (
        BigInt(1),
        BigInt(1),
        BigInt(2),
        BigInt(10),
        BigInt(114),
        BigInt(2612),
        BigInt(108664),
    )
    return p <= length(counts) ? counts[p] : nothing
end

@inline function _abs_ppt_rank(ordering_matrix::Matrix{Int}, left::Int, right::Int)
    return ordering_matrix[min(left, right), max(left, right)]
end

function _abs_ppt_charge!(state, amount::Integer)
    requested = BigInt(amount)
    limit = state.max_work
    if limit !== nothing && state.work_used + requested > limit
        state.halted = true
        state.status = AbsPPTEnumerationWorkLimit
        return false
    end
    state.work_used += requested
    return true
end

function _abs_ppt_criss_cross_valid!(ordering_matrix::Matrix{Int}, p::Int, state)
    # Collapse the innermost `n` loop into exact bit-set intersections. This
    # preserves QETLAB's six-index predicate but avoids allocating temporary
    # arrays and makes exhaustive p=6 enumeration practical. Dimensions above
    # 64 fall back to the literal predicate; their combinatorial preflight is
    # already prohibitive under the default limits.
    if p <= 64
        q = p * (p + 1) ÷ 2
        greater_masks = zeros(UInt64, p, q)
        less_masks = zeros(UInt64, p, q)
        for row in 1:p, threshold in 1:q, column in 1:p
            _abs_ppt_charge!(state, 2) || return nothing
            rank = _abs_ppt_rank(ordering_matrix, row, column)
            bit = UInt64(1) << (column - 1)
            rank > threshold && (greater_masks[row, threshold] |= bit)
            rank < threshold && (less_masks[row, threshold] |= bit)
        end
        for i in 1:p, j in 1:p, k in 1:p, l in 1:p
            _abs_ppt_charge!(state, 1) || return nothing
            _abs_ppt_rank(ordering_matrix, i, j) > _abs_ppt_rank(ordering_matrix, k, l) ||
                continue
            for m in 1:p
                _abs_ppt_charge!(state, 1) || return nothing
                possible_n =
                    greater_masks[l, _abs_ppt_rank(ordering_matrix, j, m)] &
                    less_masks[i, _abs_ppt_rank(ordering_matrix, k, m)]
                iszero(possible_n) || return false
            end
        end
        return true
    end

    for i in 1:p, j in 1:p, k in 1:p, l in 1:p, m in 1:p, n in 1:p
        _abs_ppt_charge!(state, 1) || return nothing
        if _abs_ppt_rank(ordering_matrix, i, j) > _abs_ppt_rank(ordering_matrix, k, l) &&
            _abs_ppt_rank(ordering_matrix, l, n) > _abs_ppt_rank(ordering_matrix, j, m) &&
            _abs_ppt_rank(ordering_matrix, i, n) < _abs_ppt_rank(ordering_matrix, k, m)
            return false
        end
    end
    return true
end

function _abs_ppt_ordering(ordering_matrix::Matrix{Int}, p::Int)
    q = p * (p + 1) ÷ 2
    positive_pairs = fill((0, 0), q)
    for row in 1:p, column in row:p
        positive_pairs[ordering_matrix[row, column]] = (row, column)
    end
    negative_pairs = filter(pair -> first(pair) != last(pair), positive_pairs)
    return AbsPPTOrdering(positive_pairs, negative_pairs)
end

function _abs_ppt_hard_coded_ordering(p::Int)
    if p == 1
        return AbsPPTOrdering(_AbsPPTPair[], _AbsPPTPair[])
    elseif p == 2
        return AbsPPTOrdering(_AbsPPTPair[(1, 1), (1, 2), (2, 2)], _AbsPPTPair[(1, 2)])
    end
    return throw(ArgumentError("hard-coded orderings exist only for p <= 2"))
end

function _abs_ppt_lmi_dense(
    ordering::AbsPPTOrdering, spectrum::AbstractVector{T}
) where {T<:Real}
    p = isempty(ordering.positive_pairs) ? 1 : maximum(last, ordering.positive_pairs)
    matrix = zeros(T, p, p)
    dimension = length(spectrum)
    for (rank, pair) in enumerate(ordering.positive_pairs)
        row, column = pair
        value = spectrum[dimension + 1 - rank]
        if row == column
            matrix[row, row] += 2value
        else
            matrix[row, column] += value
            matrix[column, row] += value
        end
    end
    for (rank, pair) in enumerate(ordering.negative_pairs)
        row, column = pair
        value = spectrum[rank]
        matrix[row, column] -= value
        matrix[column, row] -= value
    end
    return matrix
end

"""
    abs_ppt_lmi_matrix(ordering, spectrum; sparse_output=false)

Evaluate one explicit Hildebrand ordering at a finite real spectrum. The
spectrum must already be sorted in nonincreasing order. No normalization,
sorting, or clipping is performed by this low-level evaluator.
"""
function abs_ppt_lmi_matrix(
    ordering::AbsPPTOrdering, spectrum::AbstractVector{<:Number}; sparse_output::Bool=false
)
    values = _abs_ppt_numeric_spectrum(spectrum; allow_densify=true)
    original = collect(spectrum)
    values == original ||
        throw(ArgumentError("abs_ppt_lmi_matrix requires an already sorted spectrum"))
    required = length(ordering.positive_pairs)
    length(values) >= required || throw(
        DimensionMismatch(
            "the ordering needs at least $required eigenvalues; got $(length(values))"
        ),
    )
    matrix = _abs_ppt_lmi_dense(ordering, values)
    return sparse_output ? sparse(matrix) : matrix
end

function _abs_ppt_affine_lmi(
    ordering::AbsPPTOrdering, spectrum::AbstractVector{<:AffineScalar}, index::Int
)
    variable_count = length(first(spectrum).coefficients)
    constants = [entry.constant for entry in spectrum]
    constant_matrix = _abs_ppt_lmi_dense(ordering, constants)
    active_variables = Int[]
    for entry in spectrum
        length(entry.coefficients) == variable_count || throw(
            DimensionMismatch("all affine eigenvalues must use the same variable count")
        )
        append!(active_variables, entry.coefficients.nzind)
    end
    sort!(unique!(active_variables))
    coefficient_matrices = Vector{typeof(constant_matrix)}()
    for variable in active_variables
        coefficient_spectrum = [entry.coefficients[variable] for entry in spectrum]
        push!(coefficient_matrices, _abs_ppt_lmi_dense(ordering, coefficient_spectrum))
    end
    return HermitianAffineMatrix(
        Symbol(:abs_ppt_lmi_, index),
        constant_matrix,
        active_variables,
        coefficient_matrices,
        variable_count,
    )
end

mutable struct _AbsPPTEnumerationState
    max_constraints::Union{Nothing,BigInt}
    max_work::Union{Nothing,BigInt}
    max_entries::Union{Nothing,BigInt}
    entries_per_constraint::BigInt
    work_used::BigInt
    entries_stored::BigInt
    candidate_orderings_checked::BigInt
    halted::Bool
    status::AbsPPTEnumerationStatus
    first_violation
    first_boundary
end

function _abs_ppt_process_ordering!(
    ordering::AbsPPTOrdering,
    constraints,
    state::_AbsPPTEnumerationState,
    builder,
    psd_checker,
    stop_on_violation::Bool,
    p::Int,
)
    if state.max_constraints !== nothing &&
        BigInt(length(constraints)) >= state.max_constraints
        state.halted = true
        state.status = AbsPPTEnumerationConstraintLimit
        return false
    end
    if state.max_entries !== nothing &&
        state.entries_stored + state.entries_per_constraint > state.max_entries
        state.halted = true
        state.status = AbsPPTEnumerationEntryLimit
        return false
    end
    _abs_ppt_charge!(state, p^2) || return false
    index = length(constraints) + 1
    constraint = AbsPPTConstraint(index, ordering, builder(ordering, index))
    push!(constraints, constraint)
    state.entries_stored += state.entries_per_constraint
    psd_checker === nothing && return true

    _abs_ppt_charge!(state, p^3) || return false
    diagnostic = psd_checker(constraint)
    if diagnostic.status === MatrixPredicateViolated
        state.first_violation = (constraint=constraint, diagnostic=diagnostic)
        if stop_on_violation
            state.halted = true
            state.status = AbsPPTEnumerationEarlyViolation
            return false
        end
    elseif diagnostic.status === MatrixPredicateUnknown && state.first_boundary === nothing
        state.first_boundary = (constraint=constraint, diagnostic=diagnostic)
    end
    return true
end

function _abs_ppt_enumerate!(
    constraints,
    state::_AbsPPTEnumerationState,
    p::Int,
    builder,
    psd_checker,
    stop_on_violation::Bool,
)
    if p == 1
        return nothing
    elseif p == 2
        state.candidate_orderings_checked += 1
        _abs_ppt_process_ordering!(
            _abs_ppt_hard_coded_ordering(2),
            constraints,
            state,
            builder,
            psd_checker,
            stop_on_violation,
            p,
        )
        return nothing
    end

    q = p * (p + 1) ÷ 2
    ordering_matrix = fill(q + 1, p, p)
    used = falses(q)
    ordering_matrix[1, 1] = 1
    ordering_matrix[1, 2] = 2
    ordering_matrix[p, p] = q
    ordering_matrix[p - 1, p] = q - 1
    used[1] = true
    used[2] = true
    used[q - 1] = true
    used[q] = true

    function fill_ordering!(row::Int, column::Int, lower::Int)
        state.halted && return nothing
        upper = min(column * (column + 1) ÷ 2 + row * (p - column), q - 2)
        for rank in lower:upper
            state.halted && return nothing
            used[rank] && continue
            _abs_ppt_charge!(state, 1) || return nothing
            ordering_matrix[row, column] = rank
            used[rank] = true
            column_ordered = row == 1 || ordering_matrix[row - 1, column] < rank
            if column_ordered
                if row == p - 1 && column == p - 1
                    state.candidate_orderings_checked += 1
                    cross_valid = _abs_ppt_criss_cross_valid!(ordering_matrix, p, state)
                    cross_valid === nothing && return nothing
                    if cross_valid
                        ordering = _abs_ppt_ordering(ordering_matrix, p)
                        _abs_ppt_process_ordering!(
                            ordering,
                            constraints,
                            state,
                            builder,
                            psd_checker,
                            stop_on_violation,
                            p,
                        )
                    end
                elseif column == p
                    fill_ordering!(row + 1, row + 1, 3)
                else
                    fill_ordering!(row, column + 1, rank + 1)
                end
            end
            used[rank] = false
        end
        ordering_matrix[row, column] = q + 1
        return nothing
    end

    fill_ordering!(1, 3, 3)
    return nothing
end

function _abs_ppt_family_message(status::AbsPPTEnumerationStatus, count::Integer)
    status === AbsPPTEnumerationExhaustive &&
        return "enumerated the full monotone, criss-cross-filtered LMI family"
    status === AbsPPTEnumerationConstraintLimit &&
        return "stopped after storing $count constraints at max_constraints"
    status === AbsPPTEnumerationWorkLimit &&
        return "stopped before exceeding the deterministic max_work budget"
    status === AbsPPTEnumerationEntryLimit &&
        return "stopped before exceeding the max_entries storage budget"
    return "stopped at the first robust non-PSD constraint"
end

function _abs_ppt_build_family(
    spectrum,
    dimensions::NTuple{2,Int},
    builder;
    max_constraints,
    max_work,
    max_entries,
    entries_per_constraint,
    psd_checker,
    stop_on_violation::Bool,
    input_kind::Symbol,
    symbolic::Bool,
)
    p = min(dimensions...)
    checked_constraints = _abs_ppt_limit(max_constraints, "max_constraints")
    checked_work = _abs_ppt_limit(max_work, "max_work")
    checked_entries = _abs_ppt_limit(max_entries, "max_entries")
    per_constraint = BigInt(entries_per_constraint)
    per_constraint > 0 || throw(ArgumentError("entries_per_constraint must be positive"))
    state = _AbsPPTEnumerationState(
        checked_constraints,
        checked_work,
        checked_entries,
        per_constraint,
        BigInt(0),
        BigInt(0),
        BigInt(0),
        false,
        AbsPPTEnumerationExhaustive,
        nothing,
        nothing,
    )
    constraints = AbsPPTConstraint[]
    _abs_ppt_enumerate!(constraints, state, p, builder, psd_checker, stop_on_violation)
    status = state.halted ? state.status : AbsPPTEnumerationExhaustive
    exhaustive = status === AbsPPTEnumerationExhaustive
    limits = (
        max_constraints=checked_constraints,
        max_work=checked_work,
        max_entries=checked_entries,
        entries_per_constraint=per_constraint,
    )
    return AbsPPTConstraintFamily(
        _read_only_plan_array(collect(spectrum)),
        dimensions,
        p,
        _read_only_plan_array(constraints),
        status,
        exhaustive,
        _abs_ppt_monotone_ordering_count(p),
        _abs_ppt_known_criss_cross_count(p),
        state.candidate_orderings_checked,
        state.work_used,
        state.entries_stored,
        limits,
        state.first_violation,
        state.first_boundary,
        input_kind,
        symbolic,
        _abs_ppt_family_message(status, length(constraints)),
    )
end

"""
    abs_ppt_constraints(
        spectrum; dims=nothing, stop_on_violation=false,
        max_constraints=2612, max_work=500_000_000,
        max_entries=1_000_000, sparse_output=false,
        atol=nothing, rtol=nothing, allow_densify=false
    ) -> AbsPPTConstraintFamily

Construct the finite Hildebrand/QETLAB eigenvalue LMIs for a real numeric
spectrum. Numeric spectra are sorted in nonincreasing order. A matrix input is
also accepted; it must be finite and exactly Hermitian. General sparse matrix
inputs and sparse spectra require explicit `allow_densify=true`.

For an ordering `(σ₊,σ₋)`, the returned matrix is
`Λ(σ₊,σ₋) + transpose(Λ(σ₊,σ₋))`. Its upper triangle receives the smallest
eigenvalues in reverse ordering rank, and the off-diagonal lower contribution
subtracts the largest eigenvalues in the induced relative order. The
criss-cross test is QETLAB's documented exclusion rule.

`max_constraints`, `max_work`, and `max_entries` are checked before exceeding
their budgets; pass `nothing` only after explicitly accepting unbounded work.
Work units count recursive placements, criss-cross comparisons, matrix entries,
and (when requested) cubic PSD-check work. A limit returns a partial family
with a non-exhaustive status rather than pretending to be the full criterion.

The default `dims=nothing` is accepted only for a square bipartite dimension.
A scalar `dims=d` consistently means `(d, length(spectrum) ÷ d)`, correcting
the inconsistent scalar-dimension behavior of the two pinned MATLAB routines.
No eigenvalue is normalized, clipped, or otherwise repaired.
"""
function abs_ppt_constraints(
    spectrum::AbstractVector{<:Number};
    dims=nothing,
    stop_on_violation::Bool=false,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
    sparse_output::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    values = _abs_ppt_numeric_spectrum(spectrum; allow_densify)
    dimensions = _abs_ppt_dimensions(dims, length(values))
    p = min(dimensions...)
    real_type = eltype(values)
    scale = maximum(abs, values; init=one(real_type))
    _abs_ppt_tolerance(real_type; atol, rtol, scale)
    builder = function (ordering, _)
        matrix = _abs_ppt_lmi_dense(ordering, values)
        return sparse_output ? sparse(matrix) : matrix
    end
    psd_checker = if stop_on_violation
        constraint -> is_positive_semidefinite(
            constraint.matrix; atol, rtol, allow_densify=sparse_output
        )
    else
        nothing
    end
    return _abs_ppt_build_family(
        values,
        dimensions,
        builder;
        max_constraints,
        max_work,
        max_entries,
        entries_per_constraint=p^2,
        psd_checker,
        stop_on_violation,
        input_kind=:spectrum,
        symbolic=false,
    )
end

function abs_ppt_constraints(
    matrix::AbstractMatrix{<:Number};
    dims=nothing,
    stop_on_violation::Bool=false,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
    sparse_output::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    analysis = _abs_ppt_matrix_spectrum(
        matrix; allow_densify, max_entries, atol, rtol, boundary_allowed=false
    )
    return abs_ppt_constraints(
        analysis.spectrum;
        dims,
        stop_on_violation,
        max_constraints,
        max_work,
        max_entries,
        sparse_output,
        atol,
        rtol,
        allow_densify=true,
    )
end

"""
    abs_ppt_constraints(
        eigenvalues::AbstractVector{<:AffineScalar};
        dims=nothing, assume_ordered=false, limits...
    )

Construct package-owned affine Hermitian LMIs for scalar affine eigenvalue
expressions. Symbolic spectra cannot be sorted, so callers must pass
`assume_ordered=true` and separately impose
`λ₁ >= λ₂ >= ... >= λₙ >= 0`. The family metadata records this requirement.
No JuMP value is stored in the result.
"""
function abs_ppt_constraints(
    eigenvalues::AbstractVector{<:AffineScalar};
    dims=nothing,
    assume_ordered::Bool=false,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
)
    Base.require_one_based_indexing(eigenvalues)
    isempty(eigenvalues) &&
        throw(ArgumentError("the affine eigenvalue vector must be nonempty"))
    assume_ordered || throw(
        ArgumentError(
            "affine eigenvalues cannot be sorted; pass assume_ordered=true and " *
            "impose ordering/nonnegativity constraints in the surrounding model",
        ),
    )
    variable_count = length(first(eigenvalues).coefficients)
    all(entry -> length(entry.coefficients) == variable_count, eigenvalues) ||
        throw(DimensionMismatch("all affine eigenvalues need the same variable count"))
    dimensions = _abs_ppt_dimensions(dims, length(eigenvalues))
    p = min(dimensions...)
    entries_per_constraint = BigInt(p)^2 * BigInt(variable_count + 1)
    builder = (ordering, index) -> _abs_ppt_affine_lmi(ordering, eigenvalues, index)
    return _abs_ppt_build_family(
        eigenvalues,
        dimensions,
        builder;
        max_constraints,
        max_work,
        max_entries,
        entries_per_constraint,
        psd_checker=nothing,
        stop_on_violation=false,
        input_kind=:affine_spectrum,
        symbolic=true,
    )
end

function _abs_ppt_ordering_inequalities(ordering::AbsPPTOrdering)
    p = isempty(ordering.positive_pairs) ? 1 : maximum(last, ordering.positive_pairs)
    p <= 1 && return Vector{Vector{Int}}()
    rows = Vector{Vector{Int}}()
    for index in 1:(p - 1)
        coefficients = zeros(Int, p - 1)
        coefficients[index] = 1
        push!(rows, coefficients)
    end
    pair_coefficients = function (pair)
        coefficients = zeros(Int, p - 1)
        for endpoint in pair
            for gap in endpoint:(p - 1)
                coefficients[gap] += 1
            end
        end
        return coefficients
    end
    for rank in 1:(length(ordering.positive_pairs) - 1)
        higher = pair_coefficients(ordering.positive_pairs[rank])
        lower = pair_coefficients(ordering.positive_pairs[rank + 1])
        push!(rows, higher - lower)
    end
    return rows
end

function _abs_ppt_exact_rational(value::Integer)
    return Rational{BigInt}(BigInt(value), BigInt(1))
end

function _abs_ppt_exact_rational(value::Rational)
    return Rational{BigInt}(BigInt(numerator(value)), BigInt(denominator(value)))
end

function _abs_ppt_exact_rational(value::AbstractFloat)
    return rationalize(BigInt, value; tol=zero(value))
end

function _abs_ppt_verify_ordering_gaps(
    ordering::AbsPPTOrdering, gaps::AbstractVector{<:Real}, source::Symbol
)
    rows = _abs_ppt_ordering_inequalities(ordering)
    p = isempty(ordering.positive_pairs) ? 1 : maximum(last, ordering.positive_pairs)
    length(gaps) == max(p - 1, 0) ||
        throw(DimensionMismatch("expected $(max(p - 1, 0)) log gaps; got $(length(gaps))"))
    exact_gaps = Rational{BigInt}[_abs_ppt_exact_rational(value) for value in gaps]
    any(value -> value <= 0, exact_gaps) && return nothing
    product_margins = Rational{BigInt}[]
    for coefficients in rows[(p):end]
        margin = sum(
            Rational{BigInt}(coefficient) * exact_gaps[index] for
            (index, coefficient) in enumerate(coefficients);
            init=Rational{BigInt}(0),
        )
        margin > 0 || return nothing
        push!(product_margins, margin)
    end
    coordinates = fill(Rational{BigInt}(0), p)
    for index in (p - 1):-1:1
        coordinates[index] = coordinates[index + 1] + exact_gaps[index]
    end
    minimum_product = if isempty(product_margins)
        minimum(exact_gaps; init=Rational{BigInt}(1))
    else
        minimum(product_margins)
    end
    return AbsPPTOrderingCertificate(
        exact_gaps,
        coordinates,
        minimum(exact_gaps; init=Rational{BigInt}(1)),
        minimum_product,
        source,
    )
end

function _abs_ppt_deterministic_ordering_certificate(
    ordering::AbsPPTOrdering; max_iterations::Int
)
    p = isempty(ordering.positive_pairs) ? 1 : maximum(last, ordering.positive_pairs)
    p <= 1 && return AbsPPTOrderingCertificate(
        Rational{BigInt}[], Rational{BigInt}[0 // 1], 1 // 1, 1 // 1, :analytic
    )
    max_iterations == 0 && return nothing
    rows = _abs_ppt_ordering_inequalities(ordering)
    gaps = ones(Float64, p - 1)
    for _ in 1:max_iterations
        changed = false
        for coefficients in rows
            value = dot(coefficients, gaps)
            if value < 1.0
                norm_squared = sum(abs2, coefficients)
                norm_squared > 0 || return nothing
                gaps .+= ((1.0 - value) / norm_squared) .* coefficients
                changed = true
            end
        end
        certificate = _abs_ppt_verify_ordering_gaps(
            ordering, gaps, :deterministic_projection
        )
        certificate === nothing || return certificate
        changed || break
    end
    return nothing
end

"""
    abs_ppt_ordering_program(ordering; limits=OptimizationLimits())

Return a solver-neutral linear feasibility model that realizes an ordering in
log-gap coordinates. Feasibility of

```math
\\log x_i-\\log x_{i+1} > 0,\\qquad
\\log(x_ix_j)-\\log(x_kx_l) > 0
```

is made closed by homogeneous rescaling and represented with right-hand side
one. A solver point is never trusted directly: [`is_abs_ppt`](@ref) converts it
to exact rationals and rechecks every strict inequality before using it as a
negative absolute-PPT certificate.
"""
function abs_ppt_ordering_program(
    ordering::AbsPPTOrdering; limits::OptimizationLimits=OptimizationLimits()
)
    p = isempty(ordering.positive_pairs) ? 1 : maximum(last, ordering.positive_pairs)
    p >= 2 ||
        throw(ArgumentError("the one-dimensional ordering needs no feasibility model"))
    rows = _abs_ppt_ordering_inequalities(ordering)
    variable_count = p - 1
    intervals = AffineInterval[]
    for (index, coefficients) in enumerate(rows)
        push!(
            intervals,
            AffineInterval(
                AffineScalar(0.0, Float64.(coefficients)),
                1.0,
                nothing,
                Symbol(:ordering_margin_, index),
            ),
        )
    end
    return SemidefiniteProgram(
        :abs_ppt_ordering_realization,
        :feasibility,
        variable_count,
        AffineScalar(0.0, zeros(variable_count));
        intervals,
        initial_point=ones(variable_count),
        limits,
        metadata=(
            formulation=:hildebrand_log_product_ordering,
            local_dimension=p,
            homogeneous_margin=1,
        ),
    )
end

function _abs_ppt_backend_ordering_certificate(
    ordering::AbsPPTOrdering,
    backend::AbstractOptimizationBackend,
    optimization_limits::OptimizationLimits,
)
    problem = try
        abs_ppt_ordering_program(ordering; limits=optimization_limits)
    catch error
        if error isa ArgumentError
            return nothing, nothing, sprint(showerror, error)
        end
        rethrow(error)
    end
    result = solve_optimization(problem, backend)
    if result.status in (OptimizationOptimal, OptimizationFeasible) &&
        result.primal !== nothing
        certificate = _abs_ppt_verify_ordering_gaps(
            ordering, result.primal.coordinates, :verified_backend_primal
        )
        return certificate, result, nothing
    end
    return nothing, result, nothing
end

function _abs_ppt_result(
    status::AbsolutePPTStatus,
    verdict,
    certificate_kind,
    spectrum,
    dimensions,
    trace_value,
    family,
    violating_constraint,
    ordering_certificate,
    backend_result,
    margin,
    tolerance,
    message,
)
    checked = family === nothing ? BigInt(0) : BigInt(length(family.constraints))
    p = min(dimensions...)
    known_count = _abs_ppt_known_criss_cross_count(p)
    planned = if family === nothing
        known_count === nothing ? _abs_ppt_monotone_ordering_count(p) : known_count
    else
        if family.known_criss_cross_count === nothing
            family.monotone_ordering_count
        else
            family.known_criss_cross_count
        end
    end
    return IsAbsPPTResult(
        status,
        verdict,
        certificate_kind,
        spectrum,
        dimensions,
        trace_value,
        family,
        violating_constraint,
        ordering_certificate,
        backend_result,
        checked,
        planned,
        margin,
        tolerance,
        String(message),
    )
end

function _abs_ppt_sufficient_margins(values, p::Int)
    dimension = length(values)
    trace_value = sum(values)
    ball_margin = trace_value^2 - (dimension - 1) * sum(abs2, values)
    gershgorin_margin = if p <= 1
        zero(trace_value)
    else
        2values[end] + sum(@view(values[(dimension - p + 1):(dimension - 1)])) -
        sum(@view(values[1:(p - 1)]))
    end
    return trace_value, ball_margin, gershgorin_margin
end

function _abs_ppt_analyze_spectrum(
    values::AbstractVector{<:Real},
    dimensions::NTuple{2,Int};
    use_sufficient_tests::Bool,
    max_constraints,
    max_work,
    max_entries,
    sparse_output::Bool,
    atol,
    rtol,
    backend::AbstractOptimizationBackend,
    optimization_limits::OptimizationLimits,
    max_realization_iterations,
)
    p = min(dimensions...)
    checked_max_constraints = _abs_ppt_limit(max_constraints, "max_constraints")
    checked_max_work = _abs_ppt_limit(max_work, "max_work")
    checked_max_entries = _abs_ppt_limit(max_entries, "max_entries")
    checked_realization_iterations = _abs_ppt_nonnegative_integer(
        max_realization_iterations, "max_realization_iterations"
    )
    real_type = eltype(values)
    scale = maximum(abs, values; init=one(real_type))
    tolerance = _abs_ppt_tolerance(real_type; atol, rtol, scale)
    exact = real_type <: Integer || real_type <: Rational
    minimum_value = minimum(values)
    if exact
        minimum_value < 0 && throw(
            DomainError(minimum_value, "an absolute-PPT spectrum must be nonnegative")
        )
    elseif minimum_value < -tolerance
        throw(
            DomainError(
                minimum_value,
                "the spectrum has a negative eigenvalue outside tolerance $tolerance",
            ),
        )
    elseif minimum_value < 0
        return _abs_ppt_result(
            AbsolutePPTNumericalBoundary,
            nothing,
            nothing,
            _read_only_plan_array(collect(values)),
            dimensions,
            sum(values),
            nothing,
            nothing,
            nothing,
            nothing,
            minimum_value,
            tolerance,
            "a negative input eigenvalue lies inside the numerical PSD boundary; " *
            "the spectrum was not clipped",
        )
    end

    trace_value, ball_margin, gershgorin_margin = _abs_ppt_sufficient_margins(values, p)
    if trace_value < zero(trace_value)
        throw(DomainError(trace_value, "the spectrum trace must be nonnegative"))
    elseif iszero(trace_value)
        return _abs_ppt_result(
            AbsolutePPTAnalyticCertified,
            true,
            :zero_operator,
            _read_only_plan_array(collect(values)),
            dimensions,
            trace_value,
            nothing,
            nothing,
            nothing,
            nothing,
            zero(trace_value),
            tolerance,
            "the zero operator remains positive under every unitary and partial transpose",
        )
    elseif p == 1
        return _abs_ppt_result(
            AbsolutePPTAnalyticCertified,
            true,
            :one_dimensional_factor,
            _read_only_plan_array(collect(values)),
            dimensions,
            trace_value,
            nothing,
            nothing,
            nothing,
            nothing,
            minimum_value,
            tolerance,
            "a one-dimensional local factor makes partial transposition an ordinary " *
            "global transpose up to subsystem convention",
        )
    end

    if use_sufficient_tests
        if exact ? ball_margin >= 0 : ball_margin > tolerance
            return _abs_ppt_result(
                AbsolutePPTSufficientTestPassed,
                true,
                :gurvits_barnum_separable_ball,
                _read_only_plan_array(collect(values)),
                dimensions,
                trace_value,
                nothing,
                nothing,
                nothing,
                nothing,
                ball_margin,
                tolerance,
                "the scale-invariant Gurvits--Barnum separable-ball inequality passed",
            )
        elseif exact ? gershgorin_margin >= 0 : gershgorin_margin > tolerance
            return _abs_ppt_result(
                AbsolutePPTSufficientTestPassed,
                true,
                :gershgorin_hildebrand_lmis,
                _read_only_plan_array(collect(values)),
                dimensions,
                trace_value,
                nothing,
                nothing,
                nothing,
                nothing,
                gershgorin_margin,
                tolerance,
                "the Gershgorin sufficient condition for every Hildebrand LMI passed",
            )
        end
    end

    family = abs_ppt_constraints(
        values;
        dims=dimensions,
        stop_on_violation=true,
        max_constraints=checked_max_constraints,
        max_work=checked_max_work,
        max_entries=checked_max_entries,
        sparse_output,
        atol,
        rtol,
        allow_densify=true,
    )
    violation = family.first_violation
    if violation !== nothing
        ordering_certificate = _abs_ppt_deterministic_ordering_certificate(
            violation.constraint.ordering; max_iterations=checked_realization_iterations
        )
        backend_result = nothing
        program_limit_message = nothing
        if ordering_certificate === nothing
            ordering_certificate, backend_result, program_limit_message = _abs_ppt_backend_ordering_certificate(
                violation.constraint.ordering, backend, optimization_limits
            )
        end
        if ordering_certificate !== nothing
            return _abs_ppt_result(
                AbsolutePPTCertifiedNot,
                false,
                :realized_hildebrand_lmi_violation,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                violation,
                ordering_certificate,
                backend_result,
                violation.diagnostic.value,
                violation.diagnostic.tolerance,
                "a robust negative LMI eigenvalue and an exactly verified product-" *
                "ordering realization certify that the operator is not absolutely PPT",
            )
        elseif program_limit_message !== nothing
            return _abs_ppt_result(
                AbsolutePPTCappedUnknown,
                nothing,
                nothing,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                violation,
                nothing,
                nothing,
                violation.diagnostic.value,
                violation.diagnostic.tolerance,
                "a candidate LMI is negative, but the ordering-realization model " *
                "exceeded its deterministic pre-allocation limits: " *
                program_limit_message,
            )
        elseif backend_result !== nothing &&
            backend_result.status === OptimizationBackendUnavailable
            return _abs_ppt_result(
                AbsolutePPTBackendUnavailable,
                nothing,
                nothing,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                violation,
                nothing,
                backend_result,
                violation.diagnostic.value,
                violation.diagnostic.tolerance,
                "a candidate LMI is negative, but no exact ordering realization was " *
                "found before the bounded search and no backend is available",
            )
        elseif backend_result !== nothing &&
            backend_result.status in (OptimizationLimit, OptimizationInfeasible)
            reason = if backend_result.status === OptimizationInfeasible
                "the optional backend reported this candidate ordering infeasible; " *
                "enumeration had stopped at that candidate and therefore did not " *
                "establish a mathematical verdict"
            else
                "the optional ordering-realization backend reached its configured limit"
            end
            return _abs_ppt_result(
                AbsolutePPTCappedUnknown,
                nothing,
                nothing,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                violation,
                nothing,
                backend_result,
                violation.diagnostic.value,
                violation.diagnostic.tolerance,
                reason,
            )
        elseif backend_result !== nothing
            return _abs_ppt_result(
                AbsolutePPTBackendFailure,
                nothing,
                nothing,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                violation,
                nothing,
                backend_result,
                violation.diagnostic.value,
                violation.diagnostic.tolerance,
                "the optional ordering-realization backend did not produce a primal " *
                "point that passed exact rational verification",
            )
        end
        return _abs_ppt_result(
            AbsolutePPTCappedUnknown,
            nothing,
            nothing,
            family.spectrum,
            dimensions,
            trace_value,
            family,
            violation,
            nothing,
            nothing,
            violation.diagnostic.value,
            violation.diagnostic.tolerance,
            "a candidate LMI is negative, but its product ordering was not certified " *
            "realizable within the deterministic realization budget",
        )
    end

    if family.exhaustive
        if family.first_boundary !== nothing
            boundary = family.first_boundary
            return _abs_ppt_result(
                AbsolutePPTNumericalBoundary,
                nothing,
                nothing,
                family.spectrum,
                dimensions,
                trace_value,
                family,
                nothing,
                nothing,
                nothing,
                boundary.diagnostic.value,
                boundary.diagnostic.tolerance,
                "every robust constraint passed, but at least one exhaustive LMI " *
                "lies on the floating-point PSD boundary",
            )
        end
        minimum_margin = minimum(
            constraint -> begin
                diagnostic = is_positive_semidefinite(
                    constraint.matrix; atol, rtol, allow_densify=sparse_output
                )
                diagnostic.value
            end,
            family.constraints;
            init=typemax(real_type),
        )
        return _abs_ppt_result(
            AbsolutePPTExhaustiveCertified,
            true,
            :exhaustive_hildebrand_lmis,
            family.spectrum,
            dimensions,
            trace_value,
            family,
            nothing,
            nothing,
            nothing,
            minimum_margin,
            tolerance,
            "every LMI in the exhaustive criss-cross-filtered family is positive " *
            "semidefinite outside the numerical boundary",
        )
    end

    return _abs_ppt_result(
        AbsolutePPTCappedUnknown,
        nothing,
        nothing,
        family.spectrum,
        dimensions,
        trace_value,
        family,
        nothing,
        nothing,
        nothing,
        nothing,
        tolerance,
        "the generated LMI subfamily contains no robust violation, but enumeration " *
        "was capped and therefore cannot certify absolute PPT",
    )
end

"""
    is_abs_ppt(
        rho; dims=nothing, use_sufficient_tests=true,
        max_constraints=2612, max_work=500_000_000,
        max_entries=1_000_000, max_realization_iterations=10_000,
        backend=NoOptimizationBackend(), optimization_limits=OptimizationLimits(),
        sparse_output=false, atol=nothing, rtol=nothing,
        allow_densify=false
    ) -> IsAbsPPTResult

Analyze whether a positive operator is PPT after every global unitary
conjugation. `rho` may be a finite real eigenvalue vector or a finite Hermitian
matrix. Inputs are never normalized, clipped, or symmetrized; the property and
both sufficient tests are homogeneous in the spectrum.

The result distinguishes:

- analytic and sufficient-test certificates;
- a complete all-LMI certificate;
- a robust negative LMI with an exactly verified realizable ordering;
- capped enumeration;
- floating PSD/Hermiticity boundaries;
- missing optional backend; and
- backend failure or an unverifiable backend point.

The optional backend solves only the package-owned log-ordering feasibility
model after the bounded deterministic realization search fails. A solver point
becomes certificate evidence only after exact rational re-verification.
Passing a capped LMI family is never a positive certificate.
"""
function is_abs_ppt(
    spectrum::AbstractVector{<:Number};
    dims=nothing,
    use_sufficient_tests::Bool=true,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
    max_realization_iterations=10_000,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    optimization_limits::OptimizationLimits=OptimizationLimits(),
    sparse_output::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    values = _abs_ppt_numeric_spectrum(spectrum; allow_densify)
    dimensions = _abs_ppt_dimensions(dims, length(values))
    return _abs_ppt_analyze_spectrum(
        values,
        dimensions;
        use_sufficient_tests,
        max_constraints,
        max_work,
        max_entries,
        sparse_output,
        atol,
        rtol,
        backend,
        optimization_limits,
        max_realization_iterations,
    )
end

function is_abs_ppt(
    matrix::AbstractMatrix{<:Number};
    dims=nothing,
    use_sufficient_tests::Bool=true,
    max_constraints=2612,
    max_work=500_000_000,
    max_entries=1_000_000,
    max_realization_iterations=10_000,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    optimization_limits::OptimizationLimits=OptimizationLimits(),
    sparse_output::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    dimension = size(matrix, 1)
    dimensions = _abs_ppt_dimensions(dims, dimension)
    analysis = _abs_ppt_matrix_spectrum(
        matrix; allow_densify, max_entries, atol, rtol, boundary_allowed=true
    )
    if analysis.boundary
        real_type = _abs_ppt_real_component_type(eltype(matrix))
        zero_value = zero(real_type)
        return _abs_ppt_result(
            AbsolutePPTNumericalBoundary,
            nothing,
            nothing,
            nothing,
            dimensions,
            real(tr(matrix)),
            nothing,
            nothing,
            nothing,
            nothing,
            zero_value,
            analysis.tolerance,
            analysis.message * "; the matrix was not symmetrized",
        )
    end
    return _abs_ppt_analyze_spectrum(
        analysis.spectrum,
        dimensions;
        use_sufficient_tests,
        max_constraints,
        max_work,
        max_entries,
        sparse_output,
        atol,
        rtol,
        backend,
        optimization_limits,
        max_realization_iterations,
    )
end
