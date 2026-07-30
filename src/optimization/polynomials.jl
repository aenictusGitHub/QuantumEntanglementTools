# Source-informed independent Julia implementations based on the specifications
# and QETLAB CopositivePolynomial.m, PolynomialAsMatrix.m,
# PolynomialOptimize.m, and their polynomial-index helpers at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2012 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

using LinearAlgebra
using Random: AbstractRNG, randn
using SparseArrays

export HomogeneousPolynomial,
    PolynomialOptimizationResult,
    monomial_exponents,
    evaluate_polynomial,
    copositive_polynomial,
    polynomial_as_matrix,
    polynomial_bounds

const _POLYNOMIAL_DEFAULT_MAX_TERMS = 100_000
const _POLYNOMIAL_DEFAULT_MAX_DEGREE = 256
const _POLYNOMIAL_DEFAULT_MAX_DIMENSION = 10_000
const _POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES = 2_000_000
const _POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES = 10_000_000
const _POLYNOMIAL_DEFAULT_MAX_NONZEROS = 2_000_000
const _POLYNOMIAL_DEFAULT_MAX_WORK = 100_000_000
const _POLYNOMIAL_DEFAULT_MAX_SAMPLES = 100_000

function _polynomial_nonnegative_integer(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a nonnegative integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be a nonnegative integer; got $(repr(value))"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    return try
        Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as Int"))
    end
end

function _polynomial_positive_integer(value, name::AbstractString)
    checked = _polynomial_nonnegative_integer(value, name)
    checked > 0 || throw(ArgumentError("$name must be positive; got $checked"))
    return checked
end

function _polynomial_limit(value, name::AbstractString; allow_zero::Bool=false)
    value === nothing && return nothing
    checked = if allow_zero
        _polynomial_nonnegative_integer(value, name)
    else
        _polynomial_positive_integer(value, name)
    end
    return BigInt(checked)
end

function _polynomial_checked_sum(left::Int, right::Int, label::AbstractString)
    return try
        Base.checked_add(left, right)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("$label exceeds typemax(Int)"))
    end
end

function _polynomial_monomial_count(
    variables::Int, degree::Int; max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS
)
    count = binomial(BigInt(variables) + degree - 1, degree)
    count <= typemax(Int) || throw(
        ArgumentError(
            "the $variables-variable degree-$degree monomial count $count " *
            "cannot be represented as an array dimension",
        ),
    )
    limit = _polynomial_limit(max_terms, "max_terms")
    if limit !== nothing && count > limit
        throw(
            ArgumentError(
                "the $variables-variable degree-$degree polynomial has $count " *
                "monomials, exceeding max_terms=$limit",
            ),
        )
    end
    return Int(count)
end

function _polynomial_require_finite(value, label::AbstractString)
    value isa Number || throw(ArgumentError("$label must be numeric; got $(repr(value))"))
    finite = try
        isfinite(value)
    catch err
        err isa MethodError || rethrow()
        throw(
            ArgumentError(
                "$label has numeric type $(typeof(value)) whose finiteness " *
                "cannot be validated",
            ),
        )
    end
    finite || throw(ArgumentError("$label must be finite; got $(repr(value))"))
    return nothing
end

"""
    HomogeneousPolynomial(coefficients, variables, degree;
                          coefficient_order=:qetlab_lexicographic,
                          max_terms=100_000)

An owned coefficient representation of a homogeneous polynomial. The only
accepted coefficient order is `:qetlab_lexicographic`: row `j` corresponds to
row `j` of `monomial_exponents(variables, degree)`. Equivalently, exponents of
`x₁` decrease first, followed recursively by those of `x₂`, and so on.

Dense input coefficients remain dense and sparse vectors remain sparse.
Numeric element types are not widened by this constructor. The input is
copied, and every consuming routine revalidates its contents because Julia
array fields remain mutable.
"""
struct HomogeneousPolynomial{T,V<:AbstractVector{T}}
    variables::Int
    degree::Int
    coefficients::V

    function HomogeneousPolynomial{T,V}(
        variables::Int, degree::Int, coefficients::V
    ) where {T,V<:AbstractVector{T}}
        variables > 0 || throw(ArgumentError("variables must be positive; got $variables"))
        degree >= 0 || throw(ArgumentError("degree must be nonnegative; got $degree"))
        expected = _polynomial_monomial_count(variables, degree; max_terms=nothing)
        length(coefficients) == expected || throw(
            DimensionMismatch(
                "a $variables-variable degree-$degree polynomial requires " *
                "$expected coefficients; got $(length(coefficients))",
            ),
        )
        Base.require_one_based_indexing(coefficients)
        for index in eachindex(coefficients)
            _polynomial_require_finite(coefficients[index], "coefficient $index")
        end
        return new{T,V}(variables, degree, coefficients)
    end
end

function HomogeneousPolynomial(
    coefficients::AbstractVector,
    variables,
    degree;
    coefficient_order=:qetlab_lexicographic,
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
)
    Base.require_one_based_indexing(coefficients)
    coefficient_order === :qetlab_lexicographic || throw(
        ArgumentError(
            "coefficient_order must be :qetlab_lexicographic; got " *
            "$(repr(coefficient_order))",
        ),
    )
    variable_count = _polynomial_positive_integer(variables, "variables")
    checked_degree = _polynomial_nonnegative_integer(degree, "degree")
    expected = _polynomial_monomial_count(
        variable_count, checked_degree; max_terms=max_terms
    )
    length(coefficients) == expected || throw(
        DimensionMismatch(
            "a $variable_count-variable degree-$checked_degree polynomial " *
            "requires $expected coefficients; got $(length(coefficients))",
        ),
    )
    eltype(coefficients) <: Number || throw(
        ArgumentError(
            "coefficient storage must have a numeric element type; got " *
            "$(eltype(coefficients))",
        ),
    )
    eltype(coefficients) === Bool &&
        throw(ArgumentError("Bool is not a supported polynomial coefficient type"))
    owned = copy(coefficients)
    return HomogeneousPolynomial{eltype(owned),typeof(owned)}(
        variable_count, checked_degree, owned
    )
end

Base.eltype(::Type{<:HomogeneousPolynomial{T}}) where {T} = T
Base.eltype(::HomogeneousPolynomial{T}) where {T} = T
Base.length(polynomial::HomogeneousPolynomial) = length(polynomial.coefficients)

function Base.copy(polynomial::HomogeneousPolynomial)
    return HomogeneousPolynomial(
        copy(polynomial.coefficients),
        polynomial.variables,
        polynomial.degree;
        max_terms=nothing,
    )
end

function Base.show(io::IO, polynomial::HomogeneousPolynomial)
    storage = polynomial.coefficients isa SparseVector ? "sparse" : "dense"
    return print(
        io,
        "HomogeneousPolynomial(",
        polynomial.variables,
        " variables, degree ",
        polynomial.degree,
        ", ",
        storage,
        " ",
        eltype(polynomial),
        " coefficients)",
    )
end

function _validate_polynomial(
    polynomial::HomogeneousPolynomial; max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS
)
    expected = _polynomial_monomial_count(
        polynomial.variables, polynomial.degree; max_terms=max_terms
    )
    length(polynomial.coefficients) == expected || throw(
        DimensionMismatch(
            "the polynomial coefficient storage has length " *
            "$(length(polynomial.coefficients)); expected $expected",
        ),
    )
    Base.require_one_based_indexing(polynomial.coefficients)
    for index in eachindex(polynomial.coefficients)
        _polynomial_require_finite(polynomial.coefficients[index], "coefficient $index")
    end
    return expected
end

function _qetlab_monomial_exponents(
    variables::Int,
    degree::Int;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
)
    count = _polynomial_monomial_count(variables, degree; max_terms=max_terms)
    entry_count = BigInt(count) * variables
    entry_count <= typemax(Int) || throw(
        ArgumentError(
            "the $count × $variables monomial-exponent table has " *
            "$entry_count entries, which cannot be represented as one Julia array",
        ),
    )
    entry_limit = _polynomial_limit(
        max_exponent_entries, "max_exponent_entries"; allow_zero=true
    )
    if entry_limit !== nothing && entry_count > entry_limit
        throw(
            ArgumentError(
                "the $count × $variables monomial-exponent table has " *
                "$entry_count entries, exceeding " *
                "max_exponent_entries=$entry_limit",
            ),
        )
    end
    exponents = Matrix{Int}(undef, count, variables)
    current = zeros(Int, variables)
    current[1] = degree
    for row in 1:count
        @inbounds exponents[row, :] .= current
        row == count && break
        amount = Base.checked_add(current[end], 1)
        current[end] = 0
        for column in (variables - 1):-1:1
            if current[column] > 0
                current[column] -= 1
                current[column + 1] = Base.checked_add(current[column + 1], amount)
                break
            end
        end
    end
    return exponents
end

"""
    monomial_exponents(polynomial;
                       max_terms=100_000,
                       max_exponent_entries=2_000_000)
    monomial_exponents(variables, degree; ...)

Return the exponent rows defining the package's QETLAB-compatible coefficient
order. For two variables and degree two the rows are `[2 0; 1 1; 0 2]`, so
the corresponding coefficient vector represents
`c₁*x₁^2 + c₂*x₁*x₂ + c₃*x₂^2`.
"""
function monomial_exponents(
    variables,
    degree;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
)
    variable_count = _polynomial_positive_integer(variables, "variables")
    checked_degree = _polynomial_nonnegative_integer(degree, "degree")
    return _qetlab_monomial_exponents(
        variable_count,
        checked_degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
end

function monomial_exponents(
    polynomial::HomogeneousPolynomial;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
)
    _validate_polynomial(polynomial; max_terms=max_terms)
    return _qetlab_monomial_exponents(
        polynomial.variables,
        polynomial.degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
end

function _qetlab_monomial_index_unchecked(
    exponents::AbstractVector{<:Integer}, variables::Int, degree::Int
)
    rank = BigInt(1)
    remaining = degree
    @inbounds for variable in 1:(variables - 1)
        exponent = Int(exponents[variable])
        skipped = remaining - exponent
        if skipped > 0
            trailing_variables = variables - variable
            rank += binomial(BigInt(skipped + trailing_variables - 1), trailing_variables)
        end
        remaining -= exponent
    end
    rank <= typemax(Int) ||
        throw(ArgumentError("the monomial index cannot be represented as Int"))
    return Int(rank)
end

function _evaluate_polynomial_with_exponents(
    polynomial::HomogeneousPolynomial, point::AbstractVector, exponents::AbstractMatrix{Int}
)
    result_type = promote_type(eltype(polynomial), eltype(point))
    result_type <: Number ||
        throw(ArgumentError("polynomial and point values do not promote to a numeric type"))
    result = zero(result_type)
    @inbounds for term_index in axes(exponents, 1)
        coefficient = convert(result_type, polynomial.coefficients[term_index])
        iszero(coefficient) && continue
        term = coefficient
        for variable in axes(exponents, 2)
            exponent = exponents[term_index, variable]
            exponent == 0 && continue
            term *= convert(result_type, point[variable])^exponent
        end
        result += term
    end
    _polynomial_require_finite(result, "polynomial value")
    return result
end

"""
    evaluate_polynomial(polynomial, point;
                        max_terms=100_000,
                        max_exponent_entries=2_000_000,
                        max_work=100_000_000)

Evaluate a [`HomogeneousPolynomial`](@ref) without changing its coefficient
type or normalizing `point`. The point must contain exactly
`polynomial.variables` finite numeric entries.
"""
function evaluate_polynomial(
    polynomial::HomogeneousPolynomial,
    point::AbstractVector;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    Base.require_one_based_indexing(point)
    term_count = _validate_polynomial(polynomial; max_terms=max_terms)
    length(point) == polynomial.variables || throw(
        DimensionMismatch(
            "point has length $(length(point)); expected " * "$(polynomial.variables)"
        ),
    )
    eltype(point) <: Number ||
        throw(ArgumentError("point storage must have a numeric element type"))
    for index in eachindex(point)
        _polynomial_require_finite(point[index], "point entry $index")
    end
    work = BigInt(term_count) * (BigInt(polynomial.variables) + max(polynomial.degree, 1))
    work_limit = _polynomial_limit(max_work, "max_work")
    if work_limit !== nothing && work > work_limit
        throw(
            ArgumentError(
                "polynomial evaluation work estimate $work exceeds " *
                "max_work=$work_limit",
            ),
        )
    end
    exponents = _qetlab_monomial_exponents(
        polynomial.variables,
        polynomial.degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    return _evaluate_polynomial_with_exponents(polynomial, point, exponents)
end

function _polynomial_checked_double(value)
    doubled = if value isa Base.BitInteger
        try
            Base.checked_add(value, value)
        catch err
            err isa OverflowError || rethrow()
            throw(
                OverflowError(
                    "doubling an off-diagonal matrix coefficient overflowed " *
                    "$(typeof(value))",
                ),
            )
        end
    else
        value + value
    end
    _polynomial_require_finite(doubled, "doubled matrix coefficient")
    return doubled
end

"""
    copositive_polynomial(matrix; max_terms=100_000)

Construct the sparse degree-four polynomial

`p(x) = y' * matrix * y`, where `yᵢ = xᵢ^2`.

The matrix must be nonempty, square, finite, real, and exactly symmetric.
Unlike the pinned MATLAB routine, this function never silently replaces a
nonsymmetric input by `(matrix + matrix') / 2`. Exact coefficient types are
preserved; fixed-width integer overflow is reported.
"""
function copositive_polynomial(
    matrix::AbstractMatrix; max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS
)
    Base.require_one_based_indexing(matrix)
    rows, columns = size(matrix)
    rows == columns || throw(
        DimensionMismatch(
            "copositive_polynomial requires a square matrix; got " * "$(size(matrix))"
        ),
    )
    rows > 0 || throw(ArgumentError("copositive_polynomial requires a nonempty matrix"))
    eltype(matrix) <: Real || throw(
        ArgumentError(
            "copositive_polynomial requires a real matrix; got element type " *
            "$(eltype(matrix))",
        ),
    )
    eltype(matrix) === Bool &&
        throw(ArgumentError("Bool is not a supported matrix element type"))
    for index in eachindex(matrix)
        _polynomial_require_finite(matrix[index], "matrix entry $index")
    end
    issymmetric(matrix) || throw(
        ArgumentError(
            "matrix must be exactly symmetric; symmetrize explicitly before " *
            "calling if that is the intended mathematical input",
        ),
    )

    term_count = _polynomial_monomial_count(rows, 4; max_terms=max_terms)
    indices = Int[]
    values = Vector{eltype(matrix)}()
    exponent = zeros(Int, rows)
    for row in 1:rows
        fill!(exponent, 0)
        exponent[row] = 4
        diagonal = matrix[row, row]
        if !iszero(diagonal)
            push!(indices, _qetlab_monomial_index_unchecked(exponent, rows, 4))
            push!(values, diagonal)
        end
        for column in (row + 1):rows
            fill!(exponent, 0)
            exponent[row] = 2
            exponent[column] = 2
            off_diagonal = _polynomial_checked_double(matrix[row, column])
            if !iszero(off_diagonal)
                push!(indices, _qetlab_monomial_index_unchecked(exponent, rows, 4))
                push!(values, off_diagonal)
            end
        end
    end
    coefficients = sparsevec(indices, values, term_count)
    return HomogeneousPolynomial(coefficients, rows, 4; max_terms=max_terms)
end

function _polynomial_real_scaling_type(::Type{T}) where {T<:Number}
    floating_unit = try
        float(one(T))
    catch err
        (err isa MethodError || err isa ArgumentError) || rethrow()
        throw(
            ArgumentError(
                "cannot choose a floating scaling type for polynomial " *
                "coefficient type $T",
            ),
        )
    end
    real_unit = real(floating_unit)
    scaling_type = typeof(sqrt(real_unit))
    scaling_type <: AbstractFloat || throw(
        ArgumentError(
            "polynomial scaling requires a real floating type; got " *
            "$scaling_type from coefficient type $T",
        ),
    )
    return scaling_type
end

function _polynomial_factorials(maximum::Int)
    factorials = Vector{BigInt}(undef, maximum + 1)
    factorials[1] = BigInt(1)
    for value in 1:maximum
        factorials[value + 1] = factorials[value] * value
    end
    return factorials
end

function _polynomial_matrix_scale(
    residual_exponents,
    left_degree_exponents,
    right_degree_exponents,
    row_exponents,
    column_exponents,
    half_degree::Int,
    lifted_degree::Int,
    level::Int,
    factorials,
    ::Type{R},
) where {R<:AbstractFloat}
    numerator = factorials[level + 1] * factorials[half_degree + 1]^2
    denominator = factorials[2 * half_degree + 1] * factorials[lifted_degree + 1]
    square_root_numerator = BigInt(1)
    @inbounds for variable in eachindex(residual_exponents)
        combined = left_degree_exponents[variable] + right_degree_exponents[variable]
        numerator *= factorials[combined + 1]
        denominator *=
            factorials[residual_exponents[variable] + 1] *
            factorials[left_degree_exponents[variable] + 1] *
            factorials[right_degree_exponents[variable] + 1]
        square_root_numerator *=
            factorials[row_exponents[variable] + 1] *
            factorials[column_exponents[variable] + 1]
    end
    precision_bits = max(256, precision(BigFloat))
    high_precision_value = setprecision(BigFloat, precision_bits) do
        return BigFloat(numerator) / BigFloat(denominator) *
               sqrt(BigFloat(square_root_numerator))
    end
    value = convert(R, high_precision_value)
    isfinite(value) || throw(
        OverflowError(
            "a polynomial matrix scaling coefficient is not representable " *
            "as $R; reduce the degree or use BigFloat coefficients",
        ),
    )
    return value
end

function _polynomial_pair_rows(
    polynomial::HomogeneousPolynomial, degree_exponents::Matrix{Int}
)
    half_basis_count = size(degree_exponents, 1)
    variables = polynomial.variables
    degree = polynomial.degree
    coefficient_rows = Vector{Vector{Tuple{Int,eltype(polynomial)}}}(
        undef, half_basis_count
    )
    combined = zeros(Int, variables)
    pair_nonzeros = 0
    for left in 1:half_basis_count
        row = Tuple{Int,eltype(polynomial)}[]
        for right in 1:half_basis_count
            @inbounds for variable in 1:variables
                combined[variable] =
                    degree_exponents[left, variable] + degree_exponents[right, variable]
            end
            coefficient_index = _qetlab_monomial_index_unchecked(
                combined, variables, degree
            )
            coefficient = polynomial.coefficients[coefficient_index]
            if !iszero(coefficient)
                push!(row, (right, coefficient))
                pair_nonzeros += 1
            end
        end
        coefficient_rows[left] = row
    end
    return coefficient_rows, pair_nonzeros
end

function _polynomial_add_upper_entry!(
    entries::Dict{Int,T},
    row::Int,
    column::Int,
    value::T,
    dimension::Int,
    output_nonzeros::Int,
    max_nonzeros,
) where {T}
    iszero(value) && return output_nonzeros
    key = row + (column - 1) * dimension
    present = haskey(entries, key)
    previous = present ? entries[key] : zero(T)
    updated = previous + value
    _polynomial_require_finite(updated, "polynomial matrix entry")
    if iszero(updated)
        if present
            delete!(entries, key)
            return output_nonzeros - (row == column ? 1 : 2)
        end
        return output_nonzeros
    end
    if !present
        added = row == column ? 1 : 2
        projected = output_nonzeros + added
        if max_nonzeros !== nothing && projected > max_nonzeros
            throw(
                ArgumentError(
                    "the polynomial matrix exceeded max_nonzeros=" *
                    "$max_nonzeros during guarded sparse construction",
                ),
            )
        end
        entries[key] = updated
        return projected
    end
    entries[key] = updated
    return output_nonzeros
end

"""
    polynomial_as_matrix(polynomial; level=0, sparse_output=true,
                         max_degree=256, max_dimension=10_000,
                         max_exponent_entries=2_000_000,
                         max_dense_entries=10_000_000,
                         max_nonzeros=2_000_000, max_work=100_000_000)
    polynomial_as_matrix(coefficients, variables, half_degree; kwargs...)

Return the compact fully symmetric QETLAB matrix for an even-degree
homogeneous polynomial. If `polynomial` has degree `2d`, `level=K` constructs
the representation of
`(x₁^2 + ... + xₙ^2)^K * polynomial(x)` in the normalized degree-`d+K`
symmetric-tensor monomial basis.

The numeric construction is solver-independent. Symbolic/CVX coefficient
objects are deliberately not accepted. Sparse output is the default;
`sparse_output=false` is an explicit, guarded densification request. Integer
and rational polynomial coefficients remain exact in
[`HomogeneousPolynomial`](@ref), while the square-root normalization of this
matrix necessarily promotes them to a real floating scaling type.
"""
function polynomial_as_matrix(
    polynomial::HomogeneousPolynomial;
    level=0,
    sparse_output::Bool=true,
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    _validate_polynomial(polynomial; max_terms=max_terms)
    iseven(polynomial.degree) || throw(
        ArgumentError(
            "polynomial_as_matrix requires an even polynomial degree; got " *
            "$(polynomial.degree)",
        ),
    )
    checked_level = _polynomial_nonnegative_integer(level, "level")
    half_degree = div(polynomial.degree, 2)
    lifted_degree = _polynomial_checked_sum(
        half_degree, checked_level, "half_degree + level"
    )
    degree_limit = _polynomial_limit(max_degree, "max_degree"; allow_zero=true)
    largest_degree = max(polynomial.degree, lifted_degree)
    if degree_limit !== nothing && largest_degree > degree_limit
        throw(
            ArgumentError(
                "required degree $largest_degree exceeds max_degree=$degree_limit"
            ),
        )
    end

    half_basis_count = _polynomial_monomial_count(
        polynomial.variables, half_degree; max_terms=max_terms
    )
    matrix_dimension = _polynomial_monomial_count(
        polynomial.variables, lifted_degree; max_terms=max_terms
    )
    dimension_limit = _polynomial_limit(max_dimension, "max_dimension")
    if dimension_limit !== nothing && matrix_dimension > dimension_limit
        throw(
            ArgumentError(
                "polynomial matrix dimension $matrix_dimension exceeds " *
                "max_dimension=$dimension_limit",
            ),
        )
    end
    dense_entry_limit = _polynomial_limit(
        max_dense_entries, "max_dense_entries"; allow_zero=true
    )
    dense_entries = BigInt(matrix_dimension)^2
    dense_entries <= typemax(Int) || throw(
        ArgumentError(
            "the $matrix_dimension × $matrix_dimension polynomial matrix " *
            "has $dense_entries addressable entries, which exceeds typemax(Int)",
        ),
    )
    if !sparse_output && dense_entry_limit !== nothing && dense_entries > dense_entry_limit
        throw(
            ArgumentError(
                "dense polynomial matrix output would contain $dense_entries " *
                "entries, exceeding max_dense_entries=$dense_entry_limit",
            ),
        )
    end
    nonzero_limit = _polynomial_limit(max_nonzeros, "max_nonzeros")
    work_limit = _polynomial_limit(max_work, "max_work")
    pair_work = BigInt(half_basis_count)^2
    exponent_entry_work =
        (BigInt(half_basis_count) + matrix_dimension) * polynomial.variables
    initial_work = pair_work + exponent_entry_work
    if work_limit !== nothing && initial_work > work_limit
        throw(
            ArgumentError(
                "polynomial exponent-table and coefficient-map work estimate " *
                "$initial_work exceeds " *
                "max_work=$work_limit before matrix allocation",
            ),
        )
    end

    degree_exponents = _qetlab_monomial_exponents(
        polynomial.variables,
        half_degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    coefficient_rows, pair_nonzeros = _polynomial_pair_rows(polynomial, degree_exponents)
    construction_work = initial_work + BigInt(matrix_dimension) * pair_nonzeros
    if work_limit !== nothing && construction_work > work_limit
        throw(
            ArgumentError(
                "polynomial matrix work estimate $construction_work exceeds " *
                "max_work=$work_limit before output allocation",
            ),
        )
    end

    lifted_exponents = _qetlab_monomial_exponents(
        polynomial.variables,
        lifted_degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    maximum_factorial = max(polynomial.degree, lifted_degree)
    factorials = _polynomial_factorials(maximum_factorial)
    scaling_type = _polynomial_real_scaling_type(eltype(polynomial))
    output_type = promote_type(eltype(polynomial), scaling_type)
    entries = Dict{Int,output_type}()
    output_nonzeros = 0
    residual = zeros(Int, polynomial.variables)
    column_exponents = zeros(Int, polynomial.variables)

    for row_index in 1:matrix_dimension
        row_exponents = @view lifted_exponents[row_index, :]
        for left_index in 1:half_basis_count
            left_exponents = @view degree_exponents[left_index, :]
            admissible = true
            @inbounds for variable in 1:polynomial.variables
                residual[variable] = row_exponents[variable] - left_exponents[variable]
                if residual[variable] < 0 || residual[variable] > checked_level
                    admissible = false
                    break
                end
            end
            admissible || continue
            for (right_index, coefficient) in coefficient_rows[left_index]
                right_exponents = @view degree_exponents[right_index, :]
                @inbounds for variable in 1:polynomial.variables
                    column_exponents[variable] =
                        residual[variable] + right_exponents[variable]
                end
                column_index = _qetlab_monomial_index_unchecked(
                    column_exponents, polynomial.variables, lifted_degree
                )
                column_index < row_index && continue
                scale = _polynomial_matrix_scale(
                    residual,
                    left_exponents,
                    right_exponents,
                    row_exponents,
                    column_exponents,
                    half_degree,
                    lifted_degree,
                    checked_level,
                    factorials,
                    scaling_type,
                )
                contribution =
                    convert(output_type, coefficient) * convert(output_type, scale)
                _polynomial_require_finite(contribution, "polynomial matrix contribution")
                output_nonzeros = _polynomial_add_upper_entry!(
                    entries,
                    row_index,
                    column_index,
                    contribution,
                    matrix_dimension,
                    output_nonzeros,
                    nonzero_limit,
                )
            end
        end
    end

    row_indices = Vector{Int}(undef, output_nonzeros)
    column_indices = Vector{Int}(undef, output_nonzeros)
    values = Vector{output_type}(undef, output_nonzeros)
    cursor = 1
    for (key, value) in entries
        column = div(key - 1, matrix_dimension) + 1
        row = key - (column - 1) * matrix_dimension
        row_indices[cursor] = row
        column_indices[cursor] = column
        values[cursor] = value
        cursor += 1
        if row != column
            row_indices[cursor] = column
            column_indices[cursor] = row
            values[cursor] = value
            cursor += 1
        end
    end
    matrix = sparse(row_indices, column_indices, values, matrix_dimension, matrix_dimension)
    return sparse_output ? matrix : Matrix(matrix)
end

function polynomial_as_matrix(
    coefficients::AbstractVector,
    variables,
    half_degree;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    kwargs...,
)
    checked_half_degree = _polynomial_nonnegative_integer(half_degree, "half_degree")
    degree = try
        Base.checked_mul(2, checked_half_degree)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("2 * half_degree exceeds typemax(Int)"))
    end
    polynomial = HomogeneousPolynomial(coefficients, variables, degree; max_terms=max_terms)
    return polynomial_as_matrix(polynomial; max_terms=max_terms, kwargs...)
end

"""
    PolynomialOptimizationResult

Structured output from [`polynomial_bounds`](@ref). `outer_kind` is
`:hierarchy_outer_bound`: it records the theorem-level generalized-eigenvalue
outer bound together with a numerical reconstruction uncertainty.
`inner_kind` is `:sampled_feasible`, `:not_requested`, or
`:skipped_by_outer_target`. A sampled inner bound is an attained feasible
value, but its quality is heuristic.
"""
struct PolynomialOptimizationResult{T,V,U}
    sense::Symbol
    outer_bound::Union{Nothing,T}
    inner_bound::Union{Nothing,T}
    outer_kind::Symbol
    outer_status::Symbol
    inner_kind::Symbol
    best_point::Union{Nothing,V}
    hierarchy_level::Int
    matrix_dimension::Int
    samples_requested::Int
    samples_evaluated::Int
    outer_uncertainty::T
    sphere_residual::T
    target::U
    target_status::Symbol
    message::String
end

function Base.show(io::IO, result::PolynomialOptimizationResult)
    return print(
        io,
        "PolynomialOptimizationResult(sense=",
        result.sense,
        ", outer_bound=",
        result.outer_bound,
        ", inner_bound=",
        result.inner_bound,
        ", samples=",
        result.samples_evaluated,
        ")",
    )
end

function _polynomial_sense(sense)
    sense isa Symbol ||
        throw(ArgumentError("sense must be :min or :max; got $(repr(sense))"))
    sense in (:min, :max) ||
        throw(ArgumentError("sense must be :min or :max; got $(repr(sense))"))
    return sense
end

function _sphere_power_polynomial(
    variables::Int, half_degree::Int, ::Type{T}; max_terms, max_exponent_entries
) where {T<:AbstractFloat}
    total_degree = 2 * half_degree
    term_count = _polynomial_monomial_count(variables, total_degree; max_terms=max_terms)
    half_exponents = _qetlab_monomial_exponents(
        variables,
        half_degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    factorials = _polynomial_factorials(half_degree)
    indices = Vector{Int}(undef, size(half_exponents, 1))
    values = Vector{T}(undef, size(half_exponents, 1))
    doubled = zeros(Int, variables)
    for row in axes(half_exponents, 1)
        denominator = BigInt(1)
        @inbounds for variable in 1:variables
            exponent = half_exponents[row, variable]
            doubled[variable] = 2 * exponent
            denominator *= factorials[exponent + 1]
        end
        coefficient = factorials[half_degree + 1] ÷ denominator
        converted = convert(T, coefficient)
        isfinite(converted) || throw(
            OverflowError("a sphere-polynomial coefficient is not representable as $T")
        )
        indices[row] = _qetlab_monomial_index_unchecked(doubled, variables, total_degree)
        values[row] = converted
    end
    return HomogeneousPolynomial(
        sparsevec(indices, values, term_count), variables, total_degree; max_terms=max_terms
    )
end

function _polynomial_solver_type(::Type{T}) where {T<:Real}
    solver_type = typeof(float(zero(T)))
    solver_type <: Union{Float32,Float64} || throw(
        ArgumentError(
            "polynomial_bounds uses the standard-library Float32/Float64 " *
            "generalized symmetric eigensolver; coefficient type $T promotes " *
            "to unsupported solver type $solver_type",
        ),
    )
    return solver_type
end

function _polynomial_target(target, ::Type{T}) where {T<:AbstractFloat}
    target === nothing && return nothing
    target isa Real ||
        throw(ArgumentError("target must be a finite real number or nothing"))
    _polynomial_require_finite(target, "target")
    converted = convert(T, target)
    isfinite(converted) || throw(ArgumentError("target is not representable as $T"))
    return converted
end

function _polynomial_dense_outer_bound(
    matrix::SparseMatrixCSC{T,Int}, sphere_matrix::SparseMatrixCSC{T,Int}, sense::Symbol
) where {T<:Union{Float32,Float64}}
    dense_matrix = Matrix(matrix)
    dense_sphere = Matrix(sphere_matrix)
    factor = try
        cholesky(Symmetric(dense_sphere); check=true)
    catch err
        err isa PosDefException || rethrow()
        throw(
            ErrorException(
                "the normalized sphere polynomial matrix was not positive " *
                "definite; the hierarchy representation cannot be validated",
            ),
        )
    end
    lower = factor.L
    whitened = Matrix(lower \ dense_matrix / adjoint(lower))
    all(isfinite, whitened) || throw(
        ErrorException("the generalized-eigenvalue reduction produced non-finite entries"),
    )
    scale = max(opnorm(whitened, Inf), one(T))
    asymmetry = opnorm(whitened - transpose(whitened), Inf)
    symmetry_tolerance = T(64) * eps(T) * scale * size(whitened, 1)
    asymmetry <= symmetry_tolerance || throw(
        ErrorException(
            "the generalized-eigenvalue reduction lost symmetry: residual " *
            "$asymmetry exceeds tolerance $symmetry_tolerance",
        ),
    )
    symmetric_whitened = (whitened + transpose(whitened)) / T(2)
    decomposition = eigen(Symmetric(symmetric_whitened))
    reconstructed =
        decomposition.vectors *
        Diagonal(decomposition.values) *
        transpose(decomposition.vectors)
    reconstruction_residual = opnorm(symmetric_whitened - reconstructed, Inf)
    rounding_allowance = T(64) * eps(T) * scale * max(size(whitened, 1), 1)
    uncertainty = reconstruction_residual + asymmetry / T(2) + rounding_allowance
    raw = sense === :max ? maximum(decomposition.values) : minimum(decomposition.values)
    outer = if sense === :max
        nextfloat(raw + uncertainty)
    else
        prevfloat(raw - uncertainty)
    end
    isfinite(outer) ||
        throw(OverflowError("the polynomial hierarchy outer bound is non-finite"))
    return outer, uncertainty
end

function _polynomial_target_status(sense::Symbol, outer_bound, target)
    target === nothing && return :not_requested
    if sense === :max && outer_bound <= target
        return :outer_proves_at_most
    elseif sense === :min && outer_bound >= target
        return :outer_proves_at_least
    end
    return :not_proven
end

"""
    polynomial_bounds(rng, polynomial; level=0, sense=:max, target=nothing,
                      inner_samples=0, allow_densify=false, ...)
    polynomial_bounds(rng, coefficients, variables, half_degree, level; ...)

Compute the solver-free generalized-eigenvalue hierarchy bound for a real,
even-degree homogeneous polynomial on the real unit sphere. For `sense=:max`
the outer value is an upper bound and each sampled feasible value is a lower
bound. For `sense=:min` these directions are reversed.

The outer result is tagged `:hierarchy_outer_bound` and includes the numerical
reconstruction allowance used to pad the extremal generalized eigenvalue. An
inner result is tagged `:sampled_feasible`: it is an attained point and hence
a valid inner value, but its quality is heuristic. The routine never uses a
solver or a symbolic/CVX branch.

`rng::AbstractRNG` is mandatory. `inner_samples` is an exact deterministic
count, never derived from elapsed time. If the outer bound already establishes
the requested `target`, sampling is skipped without consuming `rng`. The
routine densifies two guarded matrices for the standard-library eigensolver,
so callers must set `allow_densify=true`.

Maximization is handled directly. This intentionally corrects the pinned
routine's recursion defect, which negates the polynomial but not `target`.
"""
function polynomial_bounds(
    rng::AbstractRNG,
    polynomial::HomogeneousPolynomial;
    level=0,
    sense=:max,
    target=nothing,
    inner_samples=0,
    allow_densify::Bool=false,
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_samples=_POLYNOMIAL_DEFAULT_MAX_SAMPLES,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    term_count = _validate_polynomial(polynomial; max_terms=max_terms)
    eltype(polynomial) <: Real || throw(
        ArgumentError(
            "polynomial_bounds requires real coefficients; got " * "$(eltype(polynomial))",
        ),
    )
    iseven(polynomial.degree) || throw(
        ArgumentError(
            "polynomial_bounds requires an even degree; got " * "$(polynomial.degree)"
        ),
    )
    checked_sense = _polynomial_sense(sense)
    checked_level = _polynomial_nonnegative_integer(level, "level")
    checked_samples = _polynomial_nonnegative_integer(inner_samples, "inner_samples")
    sample_limit = _polynomial_limit(max_samples, "max_samples"; allow_zero=true)
    if sample_limit !== nothing && checked_samples > sample_limit
        throw(
            ArgumentError(
                "inner_samples=$checked_samples exceeds max_samples=$sample_limit"
            ),
        )
    end
    allow_densify || throw(
        ArgumentError(
            "polynomial_bounds requires guarded dense generalized " *
            "eigenvalue work; pass allow_densify=true explicitly",
        ),
    )

    half_degree = div(polynomial.degree, 2)
    lifted_degree = _polynomial_checked_sum(
        half_degree, checked_level, "half_degree + level"
    )
    degree_limit = _polynomial_limit(max_degree, "max_degree"; allow_zero=true)
    largest_degree = max(polynomial.degree, lifted_degree)
    if degree_limit !== nothing && largest_degree > degree_limit
        throw(
            ArgumentError(
                "required degree $largest_degree exceeds max_degree=$degree_limit"
            ),
        )
    end
    matrix_dimension = _polynomial_monomial_count(
        polynomial.variables, lifted_degree; max_terms=max_terms
    )
    half_basis_count = _polynomial_monomial_count(
        polynomial.variables, half_degree; max_terms=max_terms
    )
    dimension_limit = _polynomial_limit(max_dimension, "max_dimension")
    if dimension_limit !== nothing && matrix_dimension > dimension_limit
        throw(
            ArgumentError(
                "polynomial matrix dimension $matrix_dimension exceeds " *
                "max_dimension=$dimension_limit",
            ),
        )
    end
    dense_entry_limit = _polynomial_limit(
        max_dense_entries, "max_dense_entries"; allow_zero=true
    )
    dense_entries = BigInt(matrix_dimension)^2
    dense_entries <= typemax(Int) || throw(
        ArgumentError(
            "the $matrix_dimension × $matrix_dimension hierarchy matrix " *
            "has $dense_entries entries, which exceeds typemax(Int)",
        ),
    )
    if dense_entry_limit !== nothing && dense_entries > dense_entry_limit
        throw(
            ArgumentError(
                "each dense hierarchy matrix would contain $dense_entries " *
                "entries, exceeding max_dense_entries=$dense_entry_limit",
            ),
        )
    end
    matrix_exponent_entries =
        (BigInt(half_basis_count) + matrix_dimension) * polynomial.variables
    full_exponent_entries = BigInt(term_count) * polynomial.variables
    exponent_entry_limit = _polynomial_limit(
        max_exponent_entries, "max_exponent_entries"; allow_zero=true
    )
    largest_exponent_table = max(
        BigInt(half_basis_count) * polynomial.variables,
        BigInt(matrix_dimension) * polynomial.variables,
        full_exponent_entries,
    )
    if exponent_entry_limit !== nothing && largest_exponent_table > exponent_entry_limit
        throw(
            ArgumentError(
                "the largest required monomial-exponent table has " *
                "$largest_exponent_table entries, exceeding " *
                "max_exponent_entries=$exponent_entry_limit",
            ),
        )
    end
    construction_work =
        BigInt(2) * BigInt(half_basis_count)^2 * (BigInt(matrix_dimension) + 1)
    eigen_work = BigInt(matrix_dimension)^3
    sampling_work =
        BigInt(checked_samples) * (
            polynomial.variables +
            BigInt(term_count) * (polynomial.variables + max(polynomial.degree, 1))
        )
    total_work =
        construction_work +
        BigInt(2) * matrix_exponent_entries +
        full_exponent_entries +
        eigen_work +
        sampling_work
    work_limit = _polynomial_limit(max_work, "max_work")
    if work_limit !== nothing && total_work > work_limit
        throw(
            ArgumentError(
                "polynomial_bounds work estimate $total_work exceeds " *
                "max_work=$work_limit before matrix construction or sampling",
            ),
        )
    end

    solver_type = _polynomial_solver_type(eltype(polynomial))
    converted_coefficients = solver_type.(polynomial.coefficients)
    converted_polynomial = HomogeneousPolynomial(
        converted_coefficients, polynomial.variables, polynomial.degree; max_terms=max_terms
    )
    sphere_polynomial = _sphere_power_polynomial(
        polynomial.variables,
        half_degree,
        solver_type;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    polynomial_matrix = polynomial_as_matrix(
        converted_polynomial;
        level=checked_level,
        sparse_output=true,
        max_terms=max_terms,
        max_degree=max_degree,
        max_dimension=max_dimension,
        max_exponent_entries=max_exponent_entries,
        max_dense_entries=max_dense_entries,
        max_nonzeros=max_nonzeros,
        max_work=nothing,
    )
    sphere_matrix = polynomial_as_matrix(
        sphere_polynomial;
        level=checked_level,
        sparse_output=true,
        max_terms=max_terms,
        max_degree=max_degree,
        max_dimension=max_dimension,
        max_exponent_entries=max_exponent_entries,
        max_dense_entries=max_dense_entries,
        max_nonzeros=max_nonzeros,
        max_work=nothing,
    )
    outer_bound, outer_uncertainty = _polynomial_dense_outer_bound(
        polynomial_matrix, sphere_matrix, checked_sense
    )

    checked_target = _polynomial_target(target, solver_type)
    target_status = _polynomial_target_status(checked_sense, outer_bound, checked_target)
    skip_sampling = target_status in (:outer_proves_at_most, :outer_proves_at_least)
    if checked_samples == 0 || skip_sampling
        inner_kind = skip_sampling ? :skipped_by_outer_target : :not_requested
        message = if skip_sampling
            "the hierarchy outer bound establishes the requested target; " *
            "no random samples were drawn"
        else
            "no inner samples were requested"
        end
        return PolynomialOptimizationResult{
            solver_type,Vector{solver_type},typeof(checked_target)
        }(
            checked_sense,
            outer_bound,
            nothing,
            :hierarchy_outer_bound,
            :validated_numerical,
            inner_kind,
            nothing,
            checked_level,
            matrix_dimension,
            checked_samples,
            0,
            outer_uncertainty,
            zero(solver_type),
            checked_target,
            target_status,
            message,
        )
    end

    exponents = _qetlab_monomial_exponents(
        polynomial.variables,
        polynomial.degree;
        max_terms=max_terms,
        max_exponent_entries=max_exponent_entries,
    )
    best_value = nothing
    best_point = nothing
    maximum_sphere_residual = zero(solver_type)
    for _ in 1:checked_samples
        point = randn(rng, solver_type, polynomial.variables)
        point_norm = norm(point)
        isfinite(point_norm) && !iszero(point_norm) || throw(
            ArgumentError("the supplied RNG produced a zero or non-finite Gaussian point"),
        )
        point ./= point_norm
        sphere_residual = abs(norm(point) - one(solver_type))
        maximum_sphere_residual = max(maximum_sphere_residual, sphere_residual)
        value = _evaluate_polynomial_with_exponents(converted_polynomial, point, exponents)
        if best_value === nothing ||
            (checked_sense === :max ? value > best_value : value < best_value)
            best_value = value
            best_point = copy(point)
        end
    end
    message =
        "the outer value is a hierarchy bound; the inner value is the " *
        "best of $checked_samples explicit feasible samples"
    return PolynomialOptimizationResult{
        solver_type,Vector{solver_type},typeof(checked_target)
    }(
        checked_sense,
        outer_bound,
        best_value,
        :hierarchy_outer_bound,
        :validated_numerical,
        :sampled_feasible,
        best_point,
        checked_level,
        matrix_dimension,
        checked_samples,
        checked_samples,
        outer_uncertainty,
        maximum_sphere_residual,
        checked_target,
        target_status,
        message,
    )
end

function polynomial_bounds(
    rng::AbstractRNG,
    coefficients::AbstractVector,
    variables,
    half_degree,
    level;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    kwargs...,
)
    checked_half_degree = _polynomial_nonnegative_integer(half_degree, "half_degree")
    degree = try
        Base.checked_mul(2, checked_half_degree)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("2 * half_degree exceeds typemax(Int)"))
    end
    polynomial = HomogeneousPolynomial(coefficients, variables, degree; max_terms=max_terms)
    return polynomial_bounds(rng, polynomial; level=level, max_terms=max_terms, kwargs...)
end
