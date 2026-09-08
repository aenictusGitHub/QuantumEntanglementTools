# Source-informed independent Julia implementations based on the specifications
# and QETLAB Majorizes.m, ElemSymPoly.m, CompoundMatrix.m, and
# AdditiveCompoundMatrix.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# ElemSymPoly.m and the compound-matrix routines also credit Benjamin Talbot.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export majorizes, elementary_symmetric_polynomial, compound_matrix, additive_compound_matrix

const _COMPOUND_DEFAULT_MAX_ENTRIES = 10_000_000
const _COMPOUND_DEFAULT_MAX_WORK = 100_000_000

function _matrix_analysis_order(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a nonnegative integer, not Bool"))
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

function _matrix_analysis_tolerance(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a nonnegative finite real number, not Bool"))
    value isa Real || throw(
        ArgumentError("$name must be a nonnegative finite real number; got $(repr(value))"),
    )
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    return value
end

function _reject_sparse_matrix_analysis_input(array, function_name::AbstractString)
    SparseArrays.issparse(array) || return nothing
    return throw(
        ArgumentError(
            "$function_name does not implicitly densify sparse inputs; " *
            "convert to a dense array explicitly before calling it",
        ),
    )
end

function _majorization_values(vector::AbstractVector)
    Base.require_one_based_indexing(vector)
    _reject_sparse_matrix_analysis_input(vector, "majorizes")
    values = collect(vector)
    for (index, value) in pairs(values)
        value isa Real || throw(
            ArgumentError(
                "vector entry $index must be real for majorization; got $(repr(value))"
            ),
        )
        isfinite(value) ||
            throw(ArgumentError("vector entry $index must be finite; got $(repr(value))"))
    end
    return values
end

function _majorization_values(matrix::AbstractMatrix)
    Base.require_one_based_indexing(matrix)
    _reject_sparse_matrix_analysis_input(matrix, "majorizes")
    for (index, value) in pairs(matrix)
        value isa Number ||
            throw(ArgumentError("matrix entry $index must be numeric; got $(repr(value))"))
        isfinite(value) ||
            throw(ArgumentError("matrix entry $index must be finite; got $(repr(value))"))
    end
    isempty(matrix) && return Float64[]
    try
        return collect(LinearAlgebra.svdvals(matrix))
    catch err
        err isa MethodError || rethrow()
        throw(
            ArgumentError(
                "majorizes requires a matrix element type supported by the " *
                "standard-library singular value decomposition; got $(eltype(matrix))",
            ),
        )
    end
end

function _majorization_zero(values)
    if !isempty(values)
        return zero(first(values))
    end
    try
        return zero(eltype(values))
    catch
        return 0
    end
end

function _majorization_pad_and_sort(values, length_out::Int)
    padded = copy(values)
    zero_value = _majorization_zero(values)
    while length(padded) < length_out
        push!(padded, zero_value)
    end
    sort!(padded; rev=true)
    return padded
end

function _majorization_default_rtol(first_values, second_values)
    default_rtol = 0
    for value in Iterators.flatten((first_values, second_values))
        (value isa Integer || value isa Rational) && continue
        float_type = typeof(float(value))
        default_rtol = max(default_rtol, sqrt(eps(float_type)))
    end
    return default_rtol
end

function _majorization_dynamic_accumulator_type(first_values, second_values)
    values = collect(Iterators.flatten((first_values, second_values)))
    isempty(values) && return BigInt
    exact = all(value -> value isa Integer || value isa Rational, values)
    types = if exact
        typeof.(_matrix_analysis_widen_exact.(values))
    else
        typeof.(values)
    end
    return foldl(promote_type, types)
end

_majorization_widened_type(::Type{T}) where {T<:Integer} = BigInt
_majorization_widened_type(::Type{T}) where {T<:Rational} = Rational{BigInt}

function _majorization_accumulator_type(
    first_values::AbstractVector{T1}, second_values::AbstractVector{T2}
) where {T1,T2}
    if isconcretetype(T1) && isconcretetype(T2)
        first_exact = T1 <: Integer || T1 <: Rational
        second_exact = T2 <: Integer || T2 <: Rational
        if first_exact && second_exact
            return promote_type(
                _majorization_widened_type(T1), _majorization_widened_type(T2)
            )
        elseif T1 <: Real && T2 <: Real
            return promote_type(T1, T2)
        end
    end
    return _majorization_dynamic_accumulator_type(first_values, second_values)
end

function _majorization_accumulator_value(value, ::Type{T}) where {T}
    widened = if value isa Integer || value isa Rational
        _matrix_analysis_widen_exact(value)
    else
        value
    end
    return convert(T, widened)
end

function _majorization_scale(first_values, second_values, ::Type{T}) where {T}
    first_scale = zero(T)
    second_scale = zero(T)
    for value in first_values
        first_scale += abs(_majorization_accumulator_value(value, T))
    end
    for value in second_values
        second_scale += abs(_majorization_accumulator_value(value, T))
    end
    return max(first_scale, second_scale, one(T))
end

_matrix_analysis_widen_exact(value::Integer) = BigInt(value)
_matrix_analysis_widen_exact(value::Rational) = Rational{BigInt}(value)
_matrix_analysis_widen_exact(value) = value

"""
    majorizes(a, b; atol=0, rtol=nothing) -> Bool

Return whether `a` majorizes `b` in the standard strong sense.

Vectors must have finite real entries. They are padded with zeros to a common
length, sorted in descending order, and compared by prefix sums; their total
sums must agree within the same tolerance. Matrix inputs may be real or complex
and are compared through their singular-value vectors. A vector and a matrix
may be compared, but a one-row `AbstractMatrix` is treated as a matrix rather
than as a vector. Matrix element types must be supported by the standard-library
singular value decomposition.

`atol` defaults to zero. If `rtol` is omitted, it is zero for exact vector
arithmetic and `sqrt(eps(T))` when floating-point values or singular values are
present. The comparison allowance is
`atol + rtol * max(sum(abs, a), sum(abs, b), 1)`.

Sparse inputs are rejected rather than silently densified. Sorting costs
`O(n log n)` for vectors; matrix inputs additionally require a dense singular
value decomposition.
"""
function majorizes(
    first::Union{AbstractVector,AbstractMatrix},
    second::Union{AbstractVector,AbstractMatrix};
    atol=0,
    rtol=nothing,
)
    checked_atol = _matrix_analysis_tolerance(atol, "atol")
    first_values = _majorization_values(first)
    second_values = _majorization_values(second)
    checked_rtol = if rtol === nothing
        _majorization_default_rtol(first_values, second_values)
    else
        _matrix_analysis_tolerance(rtol, "rtol")
    end

    common_length = max(length(first_values), length(second_values))
    sorted_first = _majorization_pad_and_sort(first_values, common_length)
    sorted_second = _majorization_pad_and_sort(second_values, common_length)
    accumulator_type = _majorization_accumulator_type(sorted_first, sorted_second)
    tolerance =
        checked_atol +
        checked_rtol * _majorization_scale(sorted_first, sorted_second, accumulator_type)

    first_prefix = zero(accumulator_type)
    second_prefix = zero(accumulator_type)
    for position in 1:common_length
        first_prefix += _majorization_accumulator_value(
            sorted_first[position], accumulator_type
        )
        second_prefix += _majorization_accumulator_value(
            sorted_second[position], accumulator_type
        )
        if position < common_length && first_prefix + tolerance < second_prefix
            return false
        end
    end
    return abs(first_prefix - second_prefix) <= tolerance
end

function _matrix_analysis_numeric_type(array)
    declared_type = eltype(array)
    if isconcretetype(declared_type) && declared_type <: Number
        return declared_type
    end

    value_type = nothing
    for (index, value) in pairs(array)
        value isa Number ||
            throw(ArgumentError("entry $index must be numeric; got $(repr(value))"))
        value_type =
            value_type === nothing ? typeof(value) : promote_type(value_type, typeof(value))
    end
    return value_type === nothing ? Int : value_type
end

function _matrix_analysis_require_finite(array, operation::AbstractString)
    Base.require_one_based_indexing(array)
    for (index, value) in pairs(array)
        value isa Number || throw(
            ArgumentError("$operation entry $index must be numeric; got $(repr(value))")
        )
        finite = try
            isfinite(value)
        catch err
            err isa MethodError || rethrow()
            throw(
                ArgumentError(
                    "$operation does not support finite-value validation for " *
                    "entry type $(typeof(value))",
                ),
            )
        end
        finite || throw(
            ArgumentError("$operation entry $index must be finite; got $(repr(value))")
        )
    end
    return nothing
end

function _matrix_analysis_result_type(value_type::Type)
    if value_type === Bool
        return Int
    elseif value_type <: Unsigned
        return BigInt
    elseif value_type <: Complex
        component_type = typeof(real(zero(value_type)))
        if component_type === Bool || component_type <: Unsigned
            return Complex{BigInt}
        end
    end
    return value_type
end

function _matrix_analysis_work_type(value_type::Type)
    if value_type <: Integer
        return BigInt
    elseif value_type <: Rational
        return Rational{BigInt}
    elseif value_type <: Complex
        component_type = typeof(real(zero(value_type)))
        if component_type <: Integer || component_type <: Rational
            return Complex{Rational{BigInt}}
        end
    end
    return value_type
end

@inline function _matrix_analysis_narrow(target_type::Type, value, context::AbstractString)
    try
        return convert(target_type, value)
    catch err
        if err isa InexactError || err isa OverflowError
            throw(OverflowError("$context is not representable as $target_type"))
        end
        rethrow()
    end
end

"""
    elementary_symmetric_polynomial(values, order)

Evaluate the elementary symmetric polynomial

```math
e_k(x_1,\\ldots,x_n) =
\\sum_{1 \\le i_1 < \\cdots < i_k \\le n}
x_{i_1}\\cdots x_{i_k}.
```

`order` must be an integer in `0:length(values)`. The empty-product convention
gives `e₀ = 1`. A descending dynamic program uses `O(length(values) * order)`
arithmetic operations and `O(order)` storage, without enumerating combinations
or densifying a sparse vector. Integer, rational, complex, and other compatible
numeric element types retain their arithmetic.
"""
function elementary_symmetric_polynomial(values::AbstractVector, order)
    checked_order = _matrix_analysis_order(order, "order")
    checked_order <= length(values) || throw(
        ArgumentError("order=$checked_order exceeds length(values)=$(length(values))")
    )
    _matrix_analysis_require_finite(values, "elementary_symmetric_polynomial")
    input_type = _matrix_analysis_numeric_type(values)
    result_type = input_type === Bool ? Int : input_type
    work_type = _matrix_analysis_work_type(result_type)
    coefficients = fill(zero(work_type), checked_order + 1)
    coefficients[1] = one(work_type)

    processed = 0
    for raw_value in values
        raw_value isa Number ||
            throw(ArgumentError("every value must be numeric; got $(repr(raw_value))"))
        value = convert(work_type, raw_value)
        maximum_degree = min(checked_order, processed + 1)
        for degree in maximum_degree:-1:1
            coefficients[degree + 1] += value * coefficients[degree]
        end
        processed += 1
    end
    return _matrix_analysis_narrow(
        result_type, coefficients[checked_order + 1], "elementary symmetric polynomial"
    )
end

function _matrix_analysis_binomial(dimension::Int, order::Int)
    order > dimension && return 0
    count = binomial(BigInt(dimension), order)
    count <= typemax(Int) || throw(
        ArgumentError(
            "binomial($dimension, $order)=$count cannot be represented " *
            "as an array dimension",
        ),
    )
    return Int(count)
end

function _matrix_analysis_resource_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _matrix_analysis_check_resource(
    function_name::AbstractString,
    planned::BigInt,
    limit,
    name::AbstractString,
    resource::AbstractString,
)
    limit === nothing && return nothing
    planned <= limit || throw(
        ArgumentError(
            "$function_name requires $planned $resource, exceeding $name=$limit; " *
            "raise the explicit guard only after reviewing the resource cost",
        ),
    )
    return nothing
end

function _compound_matrix_resource_plan(
    row_count::Int,
    column_count::Int,
    order::Int,
    sparse_output::Bool;
    max_entries,
    max_work,
)
    entry_limit = _matrix_analysis_resource_limit(max_entries, "max_entries")
    work_limit = _matrix_analysis_resource_limit(max_work, "max_work")
    rows = BigInt(row_count)
    columns = BigInt(column_count)
    checked_order = BigInt(order)
    minors = rows * columns
    minors <= typemax(Int) ||
        throw(ArgumentError("the compound-matrix output dimensions overflow Int"))

    row_workspace = row_count == 0 || column_count == 0 ? BigInt(0) : rows * checked_order
    column_workspace =
        row_count == 0 || column_count == 0 ? BigInt(0) : columns * checked_order
    row_workspace <= typemax(Int) || throw(
        ArgumentError("compound_matrix row-combination workspace exceeds typemax(Int)")
    )
    column_workspace <= typemax(Int) || throw(
        ArgumentError("compound_matrix column-combination workspace exceeds typemax(Int)"),
    )
    combination_workspace = row_workspace + column_workspace
    minor_workspace = minors == 0 ? BigInt(0) : checked_order^2
    workspace = combination_workspace + minor_workspace
    stored_entries = if sparse_output
        # Three temporary coordinate arrays coexist with the returned CSC
        # row/value arrays and column pointer while `sparse` is materialized.
        5 * minors + columns + 1 + workspace
    else
        minors + workspace
    end
    stored_entries <= typemax(Int) || throw(
        ArgumentError(
            "compound_matrix requires $stored_entries output/workspace entries, " *
            "which cannot be represented as an array length",
        ),
    )
    _matrix_analysis_check_resource(
        "compound_matrix",
        stored_entries,
        entry_limit,
        "max_entries",
        "output/workspace entries",
    )

    determinant_work = max(BigInt(1), checked_order^3)
    work = workspace + minors * determinant_work
    _matrix_analysis_check_resource(
        "compound_matrix", work, work_limit, "max_work", "estimated scalar operations"
    )
    return (entries=stored_entries, work=work, minors=minors)
end

function _additive_compound_resource_plan(
    dimension::Int, count::Int, order::Int, sparse_output::Bool; max_entries, max_work
)
    entry_limit = _matrix_analysis_resource_limit(max_entries, "max_entries")
    work_limit = _matrix_analysis_resource_limit(max_work, "max_work")
    n = BigInt(dimension)
    columns = BigInt(count)
    checked_order = BigInt(order)
    replacements = max(BigInt(0), n - checked_order)
    candidate_entries = if count == 0
        BigInt(0)
    else
        columns * (1 + checked_order * replacements)
    end
    candidate_entries <= typemax(Int) || throw(
        ArgumentError("additive_compound_matrix sparse-entry plan exceeds typemax(Int)")
    )

    # The combination matrix and copied dictionary keys each contain
    # `count * order` indices.
    workspace = 2 * columns * checked_order
    workspace <= typemax(Int) || throw(
        ArgumentError(
            "additive_compound_matrix combination workspace exceeds typemax(Int)"
        ),
    )
    stored_entries = if sparse_output
        # Three temporary coordinate arrays coexist with the returned CSC
        # row/value arrays and column pointer while `sparse` is materialized.
        5 * candidate_entries + columns + 1 + workspace
    else
        columns^2 + workspace
    end
    stored_entries <= typemax(Int) || throw(
        ArgumentError(
            "additive_compound_matrix requires $stored_entries output/workspace " *
            "entries, which cannot be represented as an array length",
        ),
    )
    _matrix_analysis_check_resource(
        "additive_compound_matrix",
        stored_entries,
        entry_limit,
        "max_entries",
        "output/workspace entries",
    )

    # Membership scans are linear in `order`; every accepted replacement also
    # copies, sorts, and hashes an order-sized key. This intentionally
    # overestimates the small fixed-order paths.
    membership_work = columns * checked_order * n * max(checked_order, BigInt(1))
    replacement_work =
        3 * columns * checked_order * replacements * max(checked_order, BigInt(1))
    work = columns + workspace + membership_work + replacement_work
    _matrix_analysis_check_resource(
        "additive_compound_matrix",
        work,
        work_limit,
        "max_work",
        "estimated scalar operations",
    )
    return (entries=stored_entries, work=work, candidates=candidate_entries)
end

function _matrix_analysis_combinations(dimension::Int, order::Int)
    count = _matrix_analysis_binomial(dimension, order)
    combinations = Matrix{Int}(undef, count, order)
    (count == 0 || order == 0) && return combinations

    current = collect(1:order)
    for row in 1:count
        @inbounds for column in 1:order
            combinations[row, column] = current[column]
        end
        row == count && break
        pivot = order
        while current[pivot] == dimension - order + pivot
            pivot -= 1
        end
        current[pivot] += 1
        @inbounds for column in (pivot + 1):order
            current[column] = current[column - 1] + 1
        end
    end
    return combinations
end

function _matrix_analysis_exact_work_type(::Type{T}) where {T}
    if T <: Integer || T <: Rational
        return true
    elseif T <: Complex
        component_type = typeof(real(zero(T)))
        return component_type <: Integer || component_type <: Rational
    end
    return false
end

function _matrix_analysis_pivoted_determinant_2x2(a11, a12, a21, a22)
    if abs(a21) > abs(a11)
        iszero(a21) && return zero(a21)
        return -a21 * (a12 - (a11 / a21) * a22)
    end
    iszero(a11) && return zero(a11)
    return a11 * (a22 - (a21 / a11) * a12)
end

function _matrix_analysis_pivoted_determinant_3x3(
    a11, a12, a13, a21, a22, a23, a31, a32, a33
)
    sign = one(a11)
    if abs(a21) > abs(a11) && abs(a21) >= abs(a31)
        a11, a21 = a21, a11
        a12, a22 = a22, a12
        a13, a23 = a23, a13
        sign = -sign
    elseif abs(a31) > abs(a11) && abs(a31) > abs(a21)
        a11, a31 = a31, a11
        a12, a32 = a32, a12
        a13, a33 = a33, a13
        sign = -sign
    end
    iszero(a11) && return zero(a11)

    multiplier21 = a21 / a11
    multiplier31 = a31 / a11
    a22 -= multiplier21 * a12
    a23 -= multiplier21 * a13
    a32 -= multiplier31 * a12
    a33 -= multiplier31 * a13

    if abs(a32) > abs(a22)
        a22, a32 = a32, a22
        a23, a33 = a33, a23
        sign = -sign
    end
    iszero(a22) && return zero(a22)

    a33 -= (a32 / a22) * a23
    return sign * a11 * a22 * a33
end

function _matrix_analysis_minor_determinant(
    matrix::AbstractMatrix, row_indices, column_indices, ::Type{R}
) where {R}
    order = length(row_indices)
    order == 0 && return one(R)
    work_type = _matrix_analysis_work_type(R)
    if order == 1
        value = convert(work_type, matrix[row_indices[1], column_indices[1]])
        return _matrix_analysis_narrow(R, value, "minor determinant")
    elseif order == 2 && (
        work_type <: LinearAlgebra.BlasFloat || _matrix_analysis_exact_work_type(work_type)
    )
        @inbounds begin
            row1, row2 = row_indices[1], row_indices[2]
            column1, column2 = column_indices[1], column_indices[2]
            a11 = convert(work_type, matrix[row1, column1])
            a12 = convert(work_type, matrix[row1, column2])
            a21 = convert(work_type, matrix[row2, column1])
            a22 = convert(work_type, matrix[row2, column2])
            value = if work_type <: LinearAlgebra.BlasFloat
                _matrix_analysis_pivoted_determinant_2x2(a11, a12, a21, a22)
            else
                a11 * a22 - a12 * a21
            end
        end
        return _matrix_analysis_narrow(R, value, "minor determinant")
    elseif order == 3 && (
        work_type <: LinearAlgebra.BlasFloat || _matrix_analysis_exact_work_type(work_type)
    )
        @inbounds begin
            row1, row2, row3 = row_indices[1], row_indices[2], row_indices[3]
            column1, column2, column3 = column_indices[1],
            column_indices[2],
            column_indices[3]
            a11 = convert(work_type, matrix[row1, column1])
            a12 = convert(work_type, matrix[row1, column2])
            a13 = convert(work_type, matrix[row1, column3])
            a21 = convert(work_type, matrix[row2, column1])
            a22 = convert(work_type, matrix[row2, column2])
            a23 = convert(work_type, matrix[row2, column3])
            a31 = convert(work_type, matrix[row3, column1])
            a32 = convert(work_type, matrix[row3, column2])
            a33 = convert(work_type, matrix[row3, column3])
            value = if work_type <: LinearAlgebra.BlasFloat
                _matrix_analysis_pivoted_determinant_3x3(
                    a11, a12, a13, a21, a22, a23, a31, a32, a33
                )
            else
                a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) +
                a13 * (a21 * a32 - a22 * a31)
            end
        end
        return _matrix_analysis_narrow(R, value, "minor determinant")
    end
    minor = Matrix{work_type}(undef, order, order)
    @inbounds for column in 1:order, row in 1:order
        minor[row, column] = convert(
            work_type, matrix[row_indices[row], column_indices[column]]
        )
    end
    determinant = LinearAlgebra.det(minor)
    return _matrix_analysis_narrow(R, determinant, "minor determinant")
end

@inline function _additive_compound_pair_index(first::Int, second::Int, dimension::Int)
    # Pairs are ordered as (1,2), (1,3), ..., (n-1,n). The resource plan has
    # already established that the dense result and this index fit in `Int`.
    return ((first - 1) * (2 * dimension - first)) ÷ 2 + second - first
end

function _additive_compound_order_two_dense(
    matrix::AbstractMatrix, dimension::Int, count::Int, ::Type{R}
) where {R}
    work_type = _matrix_analysis_work_type(R)
    result = fill(zero(R), count, count)

    @inbounds for first in 1:(dimension - 1)
        for second in (first + 1):dimension
            column = _additive_compound_pair_index(first, second, dimension)
            # Match the generic exterior-action accumulation order, including
            # IEEE signed-zero behavior: start from +zero, then add each term.
            diagonal_value = zero(work_type)
            diagonal_value += convert(work_type, matrix[first, first])
            diagonal_value += convert(work_type, matrix[second, second])
            result[column, column] = _matrix_analysis_narrow(
                R, diagonal_value, "additive-compound diagonal entry"
            )

            for replacement in 1:dimension
                (replacement == first || replacement == second) && continue

                first_value = convert(work_type, matrix[replacement, first])
                if !iszero(first_value)
                    row = if replacement < second
                        _additive_compound_pair_index(replacement, second, dimension)
                    else
                        first_value = -first_value
                        _additive_compound_pair_index(second, replacement, dimension)
                    end
                    result[row, column] = _matrix_analysis_narrow(
                        R, first_value, "additive-compound off-diagonal entry"
                    )
                end

                second_value = convert(work_type, matrix[replacement, second])
                if !iszero(second_value)
                    row = if replacement < first
                        second_value = -second_value
                        _additive_compound_pair_index(replacement, first, dimension)
                    else
                        _additive_compound_pair_index(first, replacement, dimension)
                    end
                    result[row, column] = _matrix_analysis_narrow(
                        R, second_value, "additive-compound off-diagonal entry"
                    )
                end
            end
        end
    end
    return result
end

"""
    compound_matrix(
        matrix, order;
        sparse_output=issparse(matrix),
        max_entries=10_000_000,
        max_work=100_000_000,
    )

Construct the `order`-th multiplicative compound matrix. Rows and columns use
lexicographically ordered increasing subsets. For an `m × n` input, the output
has size `binomial(m, order) × binomial(n, order)` and entry `(I,J)` is
`det(matrix[I,J])`. Thus order zero returns a `1 × 1` multiplicative identity,
and an order exceeding only one input dimension retains the mathematically
meaningful zero-by-nonzero shape.

Dense inputs produce dense outputs by default, while sparse inputs produce
sparse outputs; `sparse_output` may select either representation explicitly.
Order one returns the selected entry directly. Orders two and three use
allocation-light determinant kernels for standard BLAS floating-point and
exact element types; floating-point kernels retain partial pivoting. Other
numeric types at those orders, and all higher-order selected minors, use the
standard-library determinant on an independently materialized minor. The full
input is never implicitly densified. Integer and rational minors are evaluated
in widened exact arithmetic and narrowed only after a representability check.
The number of determinants is
`binomial(m, order) * binomial(n, order)`.

Before combination or result allocation, `max_entries` bounds a conservative
count of output, one selected-minor workspace, and combination-workspace
slots. Sparse output includes three worst-case coordinate arrays together with
the returned CSC row/value arrays and column pointer. `max_work` bounds the
determinant count weighted by a cubic minor-order estimate. Either guard may be
set to `nothing` only after independent review; array-length representability
checks cannot be disabled.
"""
function compound_matrix(
    matrix::AbstractMatrix,
    order;
    sparse_output::Bool=SparseArrays.issparse(matrix),
    max_entries=_COMPOUND_DEFAULT_MAX_ENTRIES,
    max_work=_COMPOUND_DEFAULT_MAX_WORK,
)
    checked_order = _matrix_analysis_order(order, "order")
    row_count = _matrix_analysis_binomial(size(matrix, 1), checked_order)
    column_count = _matrix_analysis_binomial(size(matrix, 2), checked_order)
    _matrix_analysis_require_finite(matrix, "compound_matrix")
    input_type = _matrix_analysis_numeric_type(matrix)
    result_type = _matrix_analysis_result_type(input_type)
    _compound_matrix_resource_plan(
        row_count,
        column_count,
        checked_order,
        sparse_output;
        max_entries=max_entries,
        max_work=max_work,
    )
    if row_count == 0 || column_count == 0
        return if sparse_output
            spzeros(result_type, row_count, column_count)
        else
            Matrix{result_type}(undef, row_count, column_count)
        end
    end

    row_indices = _matrix_analysis_combinations(size(matrix, 1), checked_order)
    column_indices = _matrix_analysis_combinations(size(matrix, 2), checked_order)
    if sparse_output
        result_rows = Int[]
        result_columns = Int[]
        result_values = result_type[]
        for column in 1:column_count, row in 1:row_count
            @views value = _matrix_analysis_minor_determinant(
                matrix, row_indices[row, :], column_indices[column, :], result_type
            )
            iszero(value) && continue
            push!(result_rows, row)
            push!(result_columns, column)
            push!(result_values, value)
        end
        return sparse(result_rows, result_columns, result_values, row_count, column_count)
    end

    result = Matrix{result_type}(undef, row_count, column_count)
    for column in 1:column_count, row in 1:row_count
        @views result[row, column] = _matrix_analysis_minor_determinant(
            matrix, row_indices[row, :], column_indices[column, :], result_type
        )
    end
    return result
end

"""
    additive_compound_matrix(
        matrix, order;
        sparse_output=issparse(matrix),
        max_entries=10_000_000,
        max_work=100_000_000,
    )

Construct the `order`-th additive compound of a square matrix. It is the
derivative at zero of the multiplicative compound,

```math
A^{[k]} =
\\left.\\frac{d}{dt}(I+tA)^{(k)}\\right|_{t=0},
```

in the lexicographically ordered exterior-power basis. Order zero returns the
`1 × 1` zero matrix, order one returns a dense copy of `matrix`, and an order
larger than the matrix dimension returns a `0 × 0` matrix. Its eigenvalues,
with algebraic multiplicity, are the `order`-fold sums of eigenvalues of
`matrix`.

The direct exterior-action construction takes
`O(binomial(n, order) * order * (n-order))` off-diagonal updates. Dense inputs
produce dense outputs by default, sparse inputs remain sparse, and
`sparse_output` can select either representation explicitly. Exact integer and
rational arithmetic is widened internally and checked when converted to the
result type. Dense order-two output uses direct lexicographic pair indexing;
higher orders and sparse output use the general exterior-action construction.

`max_entries` bounds conservative result, coordinate-array, combination-table,
and lookup-key storage before allocation. `max_work` also accounts for the
membership scans and order-sized key updates in the direct construction. Set a
guard to `nothing` only after reviewing the cost; mandatory array-length checks
remain active.
"""
function additive_compound_matrix(
    matrix::AbstractMatrix,
    order;
    sparse_output::Bool=SparseArrays.issparse(matrix),
    max_entries=_COMPOUND_DEFAULT_MAX_ENTRIES,
    max_work=_COMPOUND_DEFAULT_MAX_WORK,
)
    size(matrix, 1) == size(matrix, 2) || throw(
        DimensionMismatch(
            "additive_compound_matrix requires a square matrix; got size $(size(matrix))",
        ),
    )
    checked_order = _matrix_analysis_order(order, "order")
    dimension = size(matrix, 1)
    count = _matrix_analysis_binomial(dimension, checked_order)
    _matrix_analysis_require_finite(matrix, "additive_compound_matrix")
    input_type = _matrix_analysis_numeric_type(matrix)
    result_type = _matrix_analysis_result_type(input_type)
    _additive_compound_resource_plan(
        dimension,
        count,
        checked_order,
        sparse_output;
        max_entries=max_entries,
        max_work=max_work,
    )
    if count == 0
        return if sparse_output
            spzeros(result_type, count, count)
        else
            Matrix{result_type}(undef, count, count)
        end
    end
    if !sparse_output && checked_order == 2
        return _additive_compound_order_two_dense(matrix, dimension, count, result_type)
    end
    result = sparse_output ? nothing : fill(zero(result_type), count, count)
    result_rows = Int[]
    result_columns = Int[]
    result_values = result_type[]

    combinations = _matrix_analysis_combinations(dimension, checked_order)
    # Runtime-sized tuples force an abstract `Tuple` key type.  Stable copied
    # vectors retain content-based hashing while keeping the lookup concrete.
    lookup = Dict{Vector{Int},Int}()
    for row in 1:count
        lookup[collect(@view combinations[row, :])] = row
    end

    work_type = _matrix_analysis_work_type(result_type)
    for column in 1:count
        diagonal_value = zero(work_type)
        @inbounds for position in 1:checked_order
            index = combinations[column, position]
            diagonal_value += convert(work_type, matrix[index, index])
        end
        narrowed_diagonal = _matrix_analysis_narrow(
            result_type, diagonal_value, "additive-compound diagonal entry"
        )
        if sparse_output
            if !iszero(narrowed_diagonal)
                push!(result_rows, column)
                push!(result_columns, column)
                push!(result_values, narrowed_diagonal)
            end
        else
            result[column, column] = narrowed_diagonal
        end

        selected = collect(@view combinations[column, :])
        for removed_position in 1:checked_order
            removed = selected[removed_position]
            for replacement in 1:dimension
                replacement in selected && continue
                target = copy(selected)
                target[removed_position] = replacement
                sort!(target)
                inserted_position = searchsortedfirst(target, replacement)
                row = lookup[target]
                value = convert(work_type, matrix[replacement, removed])
                iszero(value) && continue
                if isodd(inserted_position + removed_position)
                    value = -value
                end
                narrowed_value = _matrix_analysis_narrow(
                    result_type, value, "additive-compound off-diagonal entry"
                )
                if sparse_output
                    push!(result_rows, row)
                    push!(result_columns, column)
                    push!(result_values, narrowed_value)
                else
                    result[row, column] = narrowed_value
                end
            end
        end
    end
    return if sparse_output
        sparse(result_rows, result_columns, result_values, count, count)
    else
        result
    end
end
