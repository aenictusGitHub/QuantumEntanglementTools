# Source-informed independent Julia implementations based on the specifications
# in QETLAB IsPSD.m, IsLocallyPSD.m, IsTotallyPositive.m, and
# IsTotallyNonsingular.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2012-2024 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    MatrixPredicateStatus

Three-valued outcome vocabulary for numerical matrix predicates:

- [`MatrixPredicateSatisfied`](@ref): every requested condition holds with a
  margin outside the numerical tolerance band, or holds exactly.
- [`MatrixPredicateViolated`](@ref): a violation has been established.
- [`MatrixPredicateUnknown`](@ref): the input lies on a floating-point
  tolerance boundary, so a Boolean answer would hide numerical uncertainty.
"""
@enum MatrixPredicateStatus::UInt8 begin
    MatrixPredicateSatisfied
    MatrixPredicateViolated
    MatrixPredicateUnknown
end

@doc "The matrix predicate is satisfied with a numerical margin or exact proof." MatrixPredicateSatisfied
@doc "A matrix-predicate violation has been established." MatrixPredicateViolated
@doc "The matrix predicate is numerically inconclusive at the tolerance boundary." MatrixPredicateUnknown

"""
    MatrixPredicateResult

Structured result returned by the native matrix predicates. `value` is the
decisive or boundary quantity when one is available, and `tolerance` is the
absolute band used to classify it. `witness` identifies a failing or
inconclusive entry, principal submatrix, or minor. `checked` and `planned`
record how many subproblems were inspected and how many were requested.
"""
struct MatrixPredicateResult{V,T,W}
    predicate::Symbol
    status::MatrixPredicateStatus
    value::V
    tolerance::T
    witness::W
    checked::BigInt
    planned::BigInt
    message::String
end

function Base.show(io::IO, result::MatrixPredicateResult)
    return print(
        io,
        "MatrixPredicateResult(",
        result.predicate,
        ", ",
        result.status,
        ", checked=",
        result.checked,
        "/",
        result.planned,
        ")",
    )
end

function _matrix_predicate_order(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be an integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be an integer; got $(repr(value))"))
    try
        return Int(value)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as Int"))
    end
end

function _matrix_predicate_tolerance(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a finite nonnegative real number, not Bool"))
    value isa Real || throw(
        ArgumentError("$name must be a finite nonnegative real number; got $(repr(value))"),
    )
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= zero(value) || throw(ArgumentError("$name must be nonnegative; got $value"))
    return value
end

function _matrix_predicate_type_information(matrix)
    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Number || throw(
        ArgumentError(
            "matrix predicates require a concrete numeric element type; " *
            "got $value_type",
        ),
    )
    real_type = typeof(real(zero(value_type)))
    supported_real =
        real_type <: Integer || real_type <: Rational || real_type <: AbstractFloat
    supported_real || throw(
        ArgumentError(
            "matrix predicates support integer, rational, and floating-point " *
            "real components; got $real_type",
        ),
    )
    exact = real_type <: Integer || real_type <: Rational
    return value_type, real_type, exact
end

function _matrix_predicate_tolerances(real_type::Type, exact::Bool; atol, rtol)
    absolute =
        atol === nothing ? zero(real_type) : _matrix_predicate_tolerance(atol, "atol")
    relative = if rtol === nothing
        (exact ? zero(real_type) : sqrt(eps(real_type)))
    else
        _matrix_predicate_tolerance(rtol, "rtol")
    end
    if exact && (!iszero(absolute) || !iszero(relative))
        throw(
            ArgumentError(
                "exact integer and rational inputs use exact comparisons; " *
                "atol and rtol must both be zero",
            ),
        )
    end
    return promote(absolute, relative)
end

_matrix_predicate_magnitude(value::Real) = abs(value)
_matrix_predicate_magnitude(value::Complex) = max(abs(real(value)), abs(imag(value)))

function _matrix_predicate_require_finite(matrix, operation::AbstractString)
    Base.require_one_based_indexing(matrix)
    for (index, value) in pairs(matrix)
        isfinite(value) || throw(
            ArgumentError(
                "$operation requires finite entries; entry $index is " * "$(repr(value))",
            ),
        )
    end
    return nothing
end

function _matrix_predicate_prepare(
    matrix::AbstractMatrix{<:Number}, operation::AbstractString; allow_densify::Bool
)
    _matrix_predicate_require_finite(matrix, operation)
    if SparseArrays.issparse(matrix)
        allow_densify || throw(
            ArgumentError(
                "$operation does not implicitly densify sparse input; pass " *
                "allow_densify=true after reviewing the $(size(matrix)) " *
                "dense allocation",
            ),
        )
        return Matrix(matrix)
    end
    return matrix
end

function _matrix_predicate_scale(matrix, real_type::Type, exact::Bool)
    exact && return one(real_type)
    scale = one(real_type)
    for value in matrix
        scale = max(scale, _matrix_predicate_magnitude(value))
    end
    return scale
end

function _matrix_predicate_threshold(scale, absolute, relative)
    relative_term = if iszero(relative)
        zero(promote_type(typeof(scale), typeof(relative)))
    else
        relative * scale
    end
    return absolute + relative_term
end

function _matrix_predicate_result(
    predicate::Symbol,
    status::MatrixPredicateStatus;
    value=nothing,
    tolerance=nothing,
    witness=nothing,
    checked=1,
    planned=1,
    message::AbstractString,
)
    return MatrixPredicateResult(
        predicate,
        status,
        value,
        tolerance,
        witness,
        BigInt(checked),
        BigInt(planned),
        String(message),
    )
end

function _matrix_predicate_hermiticity_defect(matrix)
    dimension = size(matrix, 1)
    first_value = matrix[1, 1] - conj(matrix[1, 1])
    maximum_defect = _matrix_predicate_magnitude(first_value)
    maximum_pair = (1, 1)
    maximum_difference = first_value
    for column in 1:dimension
        for row in 1:column
            difference = matrix[row, column] - conj(matrix[column, row])
            defect = _matrix_predicate_magnitude(difference)
            if defect > maximum_defect
                maximum_defect = defect
                maximum_pair = (row, column)
                maximum_difference = difference
            end
        end
    end
    return maximum_defect, maximum_pair, maximum_difference
end

function _matrix_predicate_psd_message(status::MatrixPredicateStatus)
    status === MatrixPredicateSatisfied &&
        return "the matrix is positive semidefinite with an exact proof or a margin outside tolerance"
    status === MatrixPredicateViolated && return "the matrix is not positive semidefinite"
    return "positive semidefiniteness is numerically inconclusive at the tolerance boundary"
end

function _matrix_predicate_psd_diagonal(matrix::Diagonal, exact::Bool, tolerance)
    diagonal = matrix.diag
    value = real(diagonal[1])
    index = 1
    for position in 2:length(diagonal)
        candidate = real(diagonal[position])
        if candidate < value
            value = candidate
            index = position
        end
    end
    status = if exact
        value < zero(value) ? MatrixPredicateViolated : MatrixPredicateSatisfied
    elseif value < -tolerance
        MatrixPredicateViolated
    elseif value > tolerance
        MatrixPredicateSatisfied
    else
        MatrixPredicateUnknown
    end
    witness =
        status === MatrixPredicateSatisfied ? nothing : (kind=:diagonal_entry, index=index)
    return _matrix_predicate_result(
        :positive_semidefinite,
        status;
        value=value,
        tolerance=tolerance,
        witness=witness,
        message=_matrix_predicate_psd_message(status),
    )
end

function _matrix_predicate_exact_scalar(value::Integer)
    return Rational{BigInt}(BigInt(value), BigInt(1))
end

function _matrix_predicate_exact_scalar(value::Rational)
    return Rational{BigInt}(BigInt(numerator(value)), BigInt(denominator(value)))
end

function _matrix_predicate_exact_scalar(value::AbstractFloat)
    return rationalize(BigInt, value; tol=zero(value))
end

function _matrix_predicate_exact_scalar(value::Complex)
    return complex(
        _matrix_predicate_exact_scalar(real(value)),
        _matrix_predicate_exact_scalar(imag(value)),
    )
end

function _matrix_predicate_exact_matrix(matrix)
    return map(_matrix_predicate_exact_scalar, matrix)
end

function _matrix_predicate_psd_ldl(matrix, exact::Bool, tolerance)
    work = exact ? _matrix_predicate_exact_matrix(matrix) : Matrix(matrix)
    dimension = size(work, 1)
    minimum_pivot = nothing

    for index in 1:dimension
        pivot = real(work[index, index])
        if minimum_pivot === nothing || pivot < minimum_pivot
            minimum_pivot = pivot
        end

        if exact
            if pivot < zero(pivot)
                return _matrix_predicate_result(
                    :positive_semidefinite,
                    MatrixPredicateViolated;
                    value=pivot,
                    tolerance=tolerance,
                    witness=(kind=:negative_ldl_pivot, index=index),
                    message=_matrix_predicate_psd_message(MatrixPredicateViolated),
                )
            elseif iszero(pivot)
                nonzero_index = findfirst(
                    row -> !iszero(work[row, index]), (index + 1):dimension
                )
                if nonzero_index !== nothing
                    row = index + nonzero_index
                    off_diagonal = work[row, index]
                    return _matrix_predicate_result(
                        :positive_semidefinite,
                        MatrixPredicateViolated;
                        value=(-abs2(off_diagonal)),
                        tolerance=tolerance,
                        witness=(
                            kind=:zero_pivot_nonzero_schur_entry,
                            indices=(index, row),
                            entry=off_diagonal,
                        ),
                        message=_matrix_predicate_psd_message(MatrixPredicateViolated),
                    )
                end
                continue
            end
        elseif pivot < -tolerance
            return _matrix_predicate_result(
                :positive_semidefinite,
                MatrixPredicateViolated;
                value=pivot,
                tolerance=tolerance,
                witness=(kind=:negative_ldl_pivot, index=index),
                message=_matrix_predicate_psd_message(MatrixPredicateViolated),
            )
        elseif pivot <= tolerance
            return _matrix_predicate_result(
                :positive_semidefinite,
                MatrixPredicateUnknown;
                value=pivot,
                tolerance=tolerance,
                witness=(kind=:ldl_pivot_boundary, index=index),
                message=_matrix_predicate_psd_message(MatrixPredicateUnknown),
            )
        end

        for column in (index + 1):dimension
            for row in column:dimension
                updated =
                    work[row, column] - work[row, index] * conj(work[column, index]) / pivot
                work[row, column] = updated
                work[column, row] = conj(updated)
            end
        end
    end

    return _matrix_predicate_result(
        :positive_semidefinite,
        MatrixPredicateSatisfied;
        value=minimum_pivot,
        tolerance=tolerance,
        message=_matrix_predicate_psd_message(MatrixPredicateSatisfied),
    )
end

function _matrix_predicate_positive_semidefinite(
    matrix, real_type::Type, exact::Bool, absolute, relative
)
    scale = _matrix_predicate_scale(matrix, real_type, exact)
    tolerance = _matrix_predicate_threshold(scale, absolute, relative)
    defect, pair, difference = _matrix_predicate_hermiticity_defect(matrix)

    if exact
        if !iszero(defect)
            return _matrix_predicate_result(
                :positive_semidefinite,
                MatrixPredicateViolated;
                value=defect,
                tolerance=zero(defect),
                witness=(
                    kind=:nonhermitian, row=pair[1], column=pair[2], difference=difference
                ),
                message="the matrix is not Hermitian; no Hermitian part was substituted",
            )
        end
    elseif defect > tolerance
        return _matrix_predicate_result(
            :positive_semidefinite,
            MatrixPredicateViolated;
            value=defect,
            tolerance=tolerance,
            witness=(
                kind=:nonhermitian, row=pair[1], column=pair[2], difference=difference
            ),
            message="the matrix is non-Hermitian outside tolerance; no Hermitian part was substituted",
        )
    elseif !iszero(defect)
        return _matrix_predicate_result(
            :positive_semidefinite,
            MatrixPredicateUnknown;
            value=defect,
            tolerance=tolerance,
            witness=(
                kind=:hermiticity_boundary,
                row=pair[1],
                column=pair[2],
                difference=difference,
            ),
            message="the Hermiticity defect lies inside tolerance; the input was not symmetrized",
        )
    end

    matrix isa Diagonal && return _matrix_predicate_psd_diagonal(matrix, exact, tolerance)

    if eltype(matrix) <: LinearAlgebra.BlasFloat
        decomposition = eigen(Hermitian(Matrix(matrix)))
        index = argmin(decomposition.values)
        value = decomposition.values[index]
        status = if value < -tolerance
            MatrixPredicateViolated
        elseif value > tolerance
            MatrixPredicateSatisfied
        else
            MatrixPredicateUnknown
        end
        witness =
            status === MatrixPredicateSatisfied ? nothing : decomposition.vectors[:, index]
        return _matrix_predicate_result(
            :positive_semidefinite,
            status;
            value=value,
            tolerance=tolerance,
            witness=witness,
            message=_matrix_predicate_psd_message(status),
        )
    end

    return _matrix_predicate_psd_ldl(matrix, exact, tolerance)
end

"""
    is_positive_semidefinite(
        matrix; atol=nothing, rtol=nothing, allow_densify=false
    ) -> MatrixPredicateResult

Determine whether a nonempty square matrix is Hermitian positive
semidefinite. Integer and rational matrices are decided exactly. For
floating-point matrices, a smallest eigenvalue or LDL pivot below
`-tolerance` violates the predicate, one above `tolerance` satisfies it, and
the intervening boundary returns `MatrixPredicateUnknown`. The default is
`atol=0` and `rtol=sqrt(eps(R))`.

The input is never symmetrized, clipped, or mutated. A nonzero Hermiticity
defect inside tolerance is `unknown`; one outside tolerance is a violation.
Sparse inputs require `allow_densify=true`. BLAS floating types use a dense
Hermitian eigendecomposition; exact and other supported types use a generic
LDL congruence calculation. Both cost `O(n^3)` time and dense `O(n^2)`
workspace, except that `Diagonal` inputs are handled in `O(n)` storage and
time.
"""
function is_positive_semidefinite(
    matrix::AbstractMatrix{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    size(matrix, 1) == size(matrix, 2) || throw(
        DimensionMismatch(
            "is_positive_semidefinite requires a square matrix; got " *
            "size $(size(matrix))",
        ),
    )
    !isempty(matrix) ||
        throw(ArgumentError("is_positive_semidefinite requires positive dimension"))
    prepared = _matrix_predicate_prepare(
        matrix, "is_positive_semidefinite"; allow_densify=allow_densify
    )
    _, real_type, exact = _matrix_predicate_type_information(prepared)
    absolute, relative = _matrix_predicate_tolerances(
        real_type, exact; atol=atol, rtol=rtol
    )
    return _matrix_predicate_positive_semidefinite(
        prepared, real_type, exact, absolute, relative
    )
end

function _matrix_predicate_advance_combination!(combination::Vector{Int}, dimension::Int)
    order = length(combination)
    position = order
    while position >= 1 && combination[position] == dimension - order + position
        position -= 1
    end
    position == 0 && return false
    combination[position] += 1
    for later in (position + 1):order
        combination[later] = combination[later - 1] + 1
    end
    return true
end

function _matrix_predicate_guard(
    planned::BigInt, limit, keyword::AbstractString, operation::AbstractString
)
    limit === nothing && return nothing
    checked_limit = _matrix_predicate_order(limit, keyword)
    checked_limit >= 1 || throw(ArgumentError("$keyword must be positive or `nothing`"))
    planned <= checked_limit || throw(
        ArgumentError(
            "$operation would inspect $planned submatrices/minors, exceeding " *
            "$keyword=$checked_limit; raise the guard only after reviewing " *
            "the combinatorial cost",
        ),
    )
    return nothing
end

"""
    is_locally_positive_semidefinite(
        matrix, order;
        atol=nothing, rtol=nothing, allow_densify=false,
        max_submatrices=100_000
    ) -> MatrixPredicateResult

Determine whether every `order × order` principal submatrix is positive
semidefinite. A robust violation takes precedence over earlier inconclusive
submatrices; otherwise any boundary submatrix makes the aggregate result
`unknown`. The witness records its one-based principal indices and nested PSD
result.

The matrix must be nonempty and square, and `1 <= order <= size(matrix, 1)`.
The routine inspects `binomial(n, order)` submatrices. `max_submatrices`
guards that combinatorial cost; setting it to `nothing` explicitly disables
the guard. Sparse inputs require `allow_densify=true`, and no submatrix is
silently symmetrized or repaired.
"""
function is_locally_positive_semidefinite(
    matrix::AbstractMatrix{<:Number},
    order;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_submatrices=100_000,
)
    size(matrix, 1) == size(matrix, 2) || throw(
        DimensionMismatch(
            "is_locally_positive_semidefinite requires a square matrix; " *
            "got size $(size(matrix))",
        ),
    )
    dimension = size(matrix, 1)
    dimension >= 1 ||
        throw(ArgumentError("is_locally_positive_semidefinite requires positive dimension"))
    checked_order = _matrix_predicate_order(order, "order")
    1 <= checked_order <= dimension ||
        throw(ArgumentError("order must lie in 1:$dimension; got $checked_order"))
    planned = binomial(BigInt(dimension), checked_order)
    _matrix_predicate_guard(
        planned, max_submatrices, "max_submatrices", "is_locally_positive_semidefinite"
    )

    prepared = _matrix_predicate_prepare(
        matrix, "is_locally_positive_semidefinite"; allow_densify=allow_densify
    )
    _, real_type, exact = _matrix_predicate_type_information(prepared)
    absolute, relative = _matrix_predicate_tolerances(
        real_type, exact; atol=atol, rtol=rtol
    )

    indices = collect(1:checked_order)
    checked = BigInt(0)
    first_unknown = nothing
    while true
        principal = prepared[indices, indices]
        subresult = _matrix_predicate_positive_semidefinite(
            principal, real_type, exact, absolute, relative
        )
        checked += 1
        witness = (indices=Tuple(indices), subresult=subresult)
        if subresult.status === MatrixPredicateViolated
            return _matrix_predicate_result(
                :locally_positive_semidefinite,
                MatrixPredicateViolated;
                value=subresult.value,
                tolerance=subresult.tolerance,
                witness=witness,
                checked=checked,
                planned=planned,
                message="a requested principal submatrix is not positive semidefinite",
            )
        elseif subresult.status === MatrixPredicateUnknown && first_unknown === nothing
            first_unknown = witness
        end
        _matrix_predicate_advance_combination!(indices, dimension) || break
    end

    if first_unknown !== nothing
        subresult = first_unknown.subresult
        return _matrix_predicate_result(
            :locally_positive_semidefinite,
            MatrixPredicateUnknown;
            value=subresult.value,
            tolerance=subresult.tolerance,
            witness=first_unknown,
            checked=checked,
            planned=planned,
            message="at least one requested principal submatrix lies on a PSD tolerance boundary",
        )
    end
    return _matrix_predicate_result(
        :locally_positive_semidefinite,
        MatrixPredicateSatisfied;
        checked=checked,
        planned=planned,
        message="every requested principal submatrix is positive semidefinite with an exact proof or numerical margin",
    )
end

function _matrix_predicate_orders(orders, maximum_order::Int)
    raw_orders = if orders === nothing
        collect(1:maximum_order)
    elseif orders isa Integer
        [orders]
    else
        try
            collect(orders)
        catch error
            error isa MethodError || rethrow()
            throw(ArgumentError("orders must be an integer or iterable of integers"))
        end
    end
    isempty(raw_orders) && throw(ArgumentError("orders must contain at least one order"))
    normalized = Int[]
    for raw_order in raw_orders
        checked = _matrix_predicate_order(raw_order, "minor order")
        1 <= checked <= maximum_order || throw(
            ArgumentError("every minor order must lie in 1:$maximum_order; got $checked"),
        )
        push!(normalized, checked)
    end
    sort!(unique!(normalized))
    return normalized
end

function _matrix_predicate_minor_count(rows::Int, columns::Int, orders::Vector{Int})
    count = BigInt(0)
    for order in orders
        count += binomial(BigInt(rows), order) * binomial(BigInt(columns), order)
    end
    return count
end

function _matrix_predicate_minor(matrix, rows::Vector{Int}, columns::Vector{Int})
    order = length(rows)
    minor = Matrix{eltype(matrix)}(undef, order, order)
    for column in 1:order
        for row in 1:order
            minor[row, column] = matrix[rows[row], columns[column]]
        end
    end
    return minor
end

function _matrix_predicate_determinant(work, exact::Bool)
    dimension = size(work, 1)
    dimension == 1 && return work[1, 1]
    determinant = one(eltype(work))

    for column in 1:dimension
        pivot_row = 0
        if exact
            for row in column:dimension
                if !iszero(work[row, column])
                    pivot_row = row
                    break
                end
            end
        else
            pivot_magnitude = zero(
                typeof(_matrix_predicate_magnitude(work[column, column]))
            )
            for row in column:dimension
                magnitude = _matrix_predicate_magnitude(work[row, column])
                if pivot_row == 0 || magnitude > pivot_magnitude
                    pivot_row = row
                    pivot_magnitude = magnitude
                end
            end
            iszero(pivot_magnitude) && (pivot_row = 0)
        end
        pivot_row == 0 && return zero(eltype(work))

        if pivot_row != column
            for entry_column in 1:dimension
                work[column, entry_column], work[pivot_row, entry_column] = work[
                    pivot_row, entry_column
                ],
                work[column, entry_column]
            end
            determinant = -determinant
        end
        pivot = work[column, column]
        determinant *= pivot
        column == dimension && continue

        for row in (column + 1):dimension
            factor = work[row, column] / pivot
            for entry_column in (column + 1):dimension
                work[row, entry_column] -= factor * work[column, entry_column]
            end
        end
    end
    return determinant
end

function _matrix_predicate_minor_determinant(minor, exact::Bool)
    work = exact ? _matrix_predicate_exact_matrix(minor) : copy(minor)
    return _matrix_predicate_determinant(work, exact)
end

function _matrix_predicate_minor_tolerance(
    minor, order::Int, real_type::Type, exact::Bool, absolute, relative
)
    exact && return zero(real_type)
    scale = _matrix_predicate_scale(minor, real_type, false)
    determinant_scale = one(scale)
    for _ in 1:order
        determinant_scale *= scale
    end
    return _matrix_predicate_threshold(determinant_scale, absolute, relative)
end

function _matrix_predicate_exact_zero_determinant(minor)
    exact_minor = _matrix_predicate_exact_matrix(minor)
    return iszero(_matrix_predicate_determinant(exact_minor, true))
end

function _matrix_predicate_positive_minor_status(minor, determinant, tolerance, exact::Bool)
    if exact
        return if determinant > zero(determinant)
            MatrixPredicateSatisfied
        else
            MatrixPredicateViolated
        end
    elseif !isfinite(determinant) || !isfinite(tolerance)
        return MatrixPredicateUnknown
    elseif determinant > tolerance
        return MatrixPredicateSatisfied
    elseif determinant < -tolerance
        return MatrixPredicateViolated
    elseif _matrix_predicate_exact_zero_determinant(minor)
        return MatrixPredicateViolated
    end
    return MatrixPredicateUnknown
end

function _matrix_predicate_nonsingular_minor_status(
    minor, determinant, tolerance, exact::Bool
)
    if exact
        return iszero(determinant) ? MatrixPredicateViolated : MatrixPredicateSatisfied
    elseif !isfinite(determinant) || !isfinite(tolerance)
        return MatrixPredicateUnknown
    elseif abs(determinant) > tolerance
        return MatrixPredicateSatisfied
    elseif _matrix_predicate_exact_zero_determinant(minor)
        return MatrixPredicateViolated
    end
    return MatrixPredicateUnknown
end

function _matrix_predicate_all_minors(
    matrix,
    predicate::Symbol,
    orders::Vector{Int},
    planned::BigInt,
    real_type::Type,
    exact::Bool,
    absolute,
    relative,
)
    checked = BigInt(0)
    first_unknown = nothing
    row_count, column_count = size(matrix)

    for order in orders
        rows = collect(1:order)
        while true
            columns = collect(1:order)
            while true
                minor = _matrix_predicate_minor(matrix, rows, columns)
                determinant = _matrix_predicate_minor_determinant(minor, exact)
                tolerance = _matrix_predicate_minor_tolerance(
                    minor, order, real_type, exact, absolute, relative
                )
                status = if predicate === :totally_positive
                    _matrix_predicate_positive_minor_status(
                        minor, determinant, tolerance, exact
                    )
                else
                    _matrix_predicate_nonsingular_minor_status(
                        minor, determinant, tolerance, exact
                    )
                end
                checked += 1
                witness = (
                    rows=Tuple(rows),
                    columns=Tuple(columns),
                    order=order,
                    determinant=determinant,
                )
                if status === MatrixPredicateViolated
                    noun = predicate === :totally_positive ? "positive" : "nonsingular"
                    return _matrix_predicate_result(
                        predicate,
                        MatrixPredicateViolated;
                        value=determinant,
                        tolerance=tolerance,
                        witness=witness,
                        checked=checked,
                        planned=planned,
                        message="a requested minor is not $noun",
                    )
                elseif status === MatrixPredicateUnknown && first_unknown === nothing
                    first_unknown = (
                        value=determinant, tolerance=tolerance, witness=witness
                    )
                end
                _matrix_predicate_advance_combination!(columns, column_count) || break
            end
            _matrix_predicate_advance_combination!(rows, row_count) || break
        end
    end

    if first_unknown !== nothing
        return _matrix_predicate_result(
            predicate,
            MatrixPredicateUnknown;
            value=first_unknown.value,
            tolerance=first_unknown.tolerance,
            witness=first_unknown.witness,
            checked=checked,
            planned=planned,
            message="at least one requested minor lies on the numerical tolerance boundary",
        )
    end
    description = predicate === :totally_positive ? "positive" : "nonsingular"
    return _matrix_predicate_result(
        predicate,
        MatrixPredicateSatisfied;
        checked=checked,
        planned=planned,
        message="every requested minor is $description with an exact proof or numerical margin",
    )
end

function _matrix_predicate_prepare_all_minors(
    matrix, operation::AbstractString, orders, max_minors, allow_densify, atol, rtol
)
    rows, columns = size(matrix)
    rows >= 1 && columns >= 1 || throw(
        ArgumentError(
            "$operation requires a matrix with positive row and column dimensions"
        ),
    )
    normalized_orders = _matrix_predicate_orders(orders, min(rows, columns))
    planned = _matrix_predicate_minor_count(rows, columns, normalized_orders)
    _matrix_predicate_guard(planned, max_minors, "max_minors", operation)
    prepared = _matrix_predicate_prepare(matrix, operation; allow_densify=allow_densify)
    _, real_type, exact = _matrix_predicate_type_information(prepared)
    absolute, relative = _matrix_predicate_tolerances(
        real_type, exact; atol=atol, rtol=rtol
    )
    return (prepared, normalized_orders, planned, real_type, exact, absolute, relative)
end

"""
    is_totally_positive(
        matrix; orders=nothing, atol=nothing, rtol=nothing,
        allow_densify=false, max_minors=100_000
    ) -> MatrixPredicateResult

Determine whether every requested square minor of a nonempty real matrix has
strictly positive determinant. `orders=nothing` selects
`1:min(size(matrix)...)`. Exact integer and rational inputs use exact
determinants. For floating-point inputs, determinants above their tolerance
band pass, those below its negative edge fail, and nonzero determinants inside
the band are `unknown`. An exactly singular represented minor is a violation
of strict positivity.

For an order-`k` minor, the band is
`atol + rtol * max(1, maximum(abs, minor))^k`; defaults are zero and
`sqrt(eps(R))`. The work is combinatorial:
`sum(binomial(m,k) * binomial(n,k) for k in orders)` determinants, each taking
`O(k^3)` arithmetic. `max_minors` guards this cost; `nothing` explicitly
disables the guard. Sparse input requires `allow_densify=true`. The input is
never normalized, sign-adjusted, or mutated.
"""
function is_totally_positive(
    matrix::AbstractMatrix{<:Number};
    orders=nothing,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_minors=100_000,
)
    eltype(matrix) <: Real || throw(
        ArgumentError(
            "is_totally_positive is defined for real matrices; got " *
            "eltype $(eltype(matrix))",
        ),
    )
    prepared, normalized_orders, planned, real_type, exact, absolute, relative = _matrix_predicate_prepare_all_minors(
        matrix, "is_totally_positive", orders, max_minors, allow_densify, atol, rtol
    )
    return _matrix_predicate_all_minors(
        prepared,
        :totally_positive,
        normalized_orders,
        planned,
        real_type,
        exact,
        absolute,
        relative,
    )
end

"""
    is_totally_nonsingular(
        matrix; orders=nothing, atol=nothing, rtol=nothing,
        allow_densify=false, max_minors=100_000
    ) -> MatrixPredicateResult

Determine whether every requested square minor of a nonempty real or complex
matrix is nonsingular. Exact integer and rational inputs are decided exactly.
A floating determinant whose magnitude exceeds its tolerance passes. A
boundary determinant returns `unknown`, unless exact rational reconstruction
of the represented floating entries proves that the minor is singular.

`orders=nothing` selects every order through `min(size(matrix)...)`. The
order-`k` determinant band and combinatorial `max_minors` guard are the same as
for [`is_totally_positive`](@ref). Each checked minor costs `O(k^3)`
arithmetic. Sparse input requires `allow_densify=true`, and the input is never
mutated or implicitly repaired.
"""
function is_totally_nonsingular(
    matrix::AbstractMatrix{<:Number};
    orders=nothing,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_minors=100_000,
)
    prepared, normalized_orders, planned, real_type, exact, absolute, relative = _matrix_predicate_prepare_all_minors(
        matrix, "is_totally_nonsingular", orders, max_minors, allow_densify, atol, rtol
    )
    return _matrix_predicate_all_minors(
        prepared,
        :totally_nonsingular,
        normalized_orders,
        planned,
        real_type,
        exact,
        absolute,
        relative,
    )
end
