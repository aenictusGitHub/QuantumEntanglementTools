# Source-informed independent Julia implementations based on the specifications
# and QETLAB Majorizes.m, ElemSymPoly.m, CompoundMatrix.m, and
# AdditiveCompoundMatrix.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# ElemSymPoly.m and the compound-matrix routines also credit Benjamin Talbot.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export majorizes, elementary_symmetric_polynomial, compound_matrix, additive_compound_matrix

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

function _majorization_scale(first_values, second_values)
    first_scale = BigInt(0)
    second_scale = BigInt(0)
    for value in first_values
        first_scale += abs(_matrix_analysis_widen_exact(value))
    end
    for value in second_values
        second_scale += abs(_matrix_analysis_widen_exact(value))
    end
    return max(first_scale, second_scale, 1)
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
    tolerance =
        checked_atol + checked_rtol * _majorization_scale(sorted_first, sorted_second)

    first_prefix = BigInt(0)
    second_prefix = BigInt(0)
    for position in 1:common_length
        first_prefix += _matrix_analysis_widen_exact(sorted_first[position])
        second_prefix += _matrix_analysis_widen_exact(sorted_second[position])
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

function _matrix_analysis_narrow(target_type::Type, value, context::AbstractString)
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

function _matrix_analysis_minor_determinant(
    matrix::AbstractMatrix, row_indices, column_indices, result_type::Type
)
    order = length(row_indices)
    order == 0 && return one(result_type)
    work_type = _matrix_analysis_work_type(result_type)
    minor = Matrix{work_type}(undef, order, order)
    @inbounds for column in 1:order, row in 1:order
        minor[row, column] = convert(
            work_type, matrix[row_indices[row], column_indices[column]]
        )
    end
    determinant = LinearAlgebra.det(minor)
    return _matrix_analysis_narrow(result_type, determinant, "minor determinant")
end

"""
    compound_matrix(matrix, order; sparse_output=issparse(matrix))

Construct the `order`-th multiplicative compound matrix. Rows and columns use
lexicographically ordered increasing subsets. For an `m × n` input, the output
has size `binomial(m, order) × binomial(n, order)` and entry `(I,J)` is
`det(matrix[I,J])`. Thus order zero returns a `1 × 1` multiplicative identity,
and an order exceeding only one input dimension retains the mathematically
meaningful zero-by-nonzero shape.

Dense inputs produce dense outputs by default, while sparse inputs produce
sparse outputs; `sparse_output` may select either representation explicitly.
Each selected `order × order` minor is materialized independently, but the full
input is never implicitly densified. Integer and rational minors are evaluated
in widened exact arithmetic and narrowed only after a representability check.
The number of determinants is
`binomial(m, order) * binomial(n, order)`.
"""
function compound_matrix(
    matrix::AbstractMatrix, order; sparse_output::Bool=SparseArrays.issparse(matrix)
)
    checked_order = _matrix_analysis_order(order, "order")
    row_count = _matrix_analysis_binomial(size(matrix, 1), checked_order)
    column_count = _matrix_analysis_binomial(size(matrix, 2), checked_order)
    try
        Base.checked_mul(row_count, column_count)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the compound-matrix output dimensions overflow Int"))
    end

    _matrix_analysis_require_finite(matrix, "compound_matrix")
    input_type = _matrix_analysis_numeric_type(matrix)
    result_type = _matrix_analysis_result_type(input_type)
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
    additive_compound_matrix(matrix, order;
                             sparse_output=issparse(matrix))

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
result type.
"""
function additive_compound_matrix(
    matrix::AbstractMatrix, order; sparse_output::Bool=SparseArrays.issparse(matrix)
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
    if count == 0
        return if sparse_output
            spzeros(result_type, count, count)
        else
            Matrix{result_type}(undef, count, count)
        end
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
