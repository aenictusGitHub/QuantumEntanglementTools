# Source-informed independent Julia implementation based on the specifications
# of QETLAB ApplyMap.m, ChoiMatrix.m, KrausOperators.m, and
# IsHermPreserving.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export AbstractMapRepresentation,
    OperatorSpace,
    KrausRepresentation,
    OperatorSumRepresentation,
    OperatorSumDecompositionResult,
    CanonicalMapDecompositionResult,
    ChoiRepresentation,
    SuperoperatorRepresentation,
    operator_space,
    input_size,
    output_size,
    input_dimension,
    output_dimension,
    operator_sum_factors,
    operator_sum_representation,
    operator_sum_decomposition,
    canonical_map_decomposition,
    kraus_representation,
    choi_representation,
    superoperator_representation,
    kraus_operators,
    choi_matrix,
    superoperator_matrix,
    apply_channel,
    is_hermiticity_preserving,
    is_completely_positive,
    is_trace_preserving,
    is_unital,
    dual_channel,
    complementary_channel,
    partial_map,
    depolarizing_channel,
    dephasing_channel,
    pauli_channel,
    choi_map,
    reduction_map

"""
    OperatorSpace(input_rows, input_columns, output_rows, output_columns)
    OperatorSpace(input_size, output_size)

Immutable descriptor for a linear map from `input_rows × input_columns`
matrices to `output_rows × output_columns` matrices. All dimensions must be
positive integers and all products needed by the Choi and superoperator
representations must fit in `Int`.

Use [`input_size`](@ref) and [`output_size`](@ref) for every operator space.
The scalar [`input_dimension`](@ref) and [`output_dimension`](@ref) accessors
are intentionally defined only when the corresponding matrix space is square.
"""
struct OperatorSpace
    input_rows::Int
    input_columns::Int
    output_rows::Int
    output_columns::Int

    function OperatorSpace(
        token::_ValidatedConstructorToken,
        input_rows::Int,
        input_columns::Int,
        output_rows::Int,
        output_columns::Int,
    )
        _require_validated_constructor_token(token)
        all(>(0), (input_rows, input_columns, output_rows, output_columns)) ||
            throw(ArgumentError("validated operator-space dimensions must be positive"))
        _checked_product((input_rows, output_rows), "operator-space Choi row dimensions")
        _checked_product(
            (input_columns, output_columns), "operator-space Choi column dimensions"
        )
        _checked_product((output_rows, output_columns), "operator-space output dimensions")
        _checked_product((input_rows, input_columns), "operator-space input dimensions")
        return new(input_rows, input_columns, output_rows, output_columns)
    end
end

function OperatorSpace(input_rows, input_columns, output_rows, output_columns)
    checked_input_rows = _positive_int(input_rows, "input_rows")
    checked_input_columns = _positive_int(input_columns, "input_columns")
    checked_output_rows = _positive_int(output_rows, "output_rows")
    checked_output_columns = _positive_int(output_columns, "output_columns")
    return try
        OperatorSpace(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            checked_input_rows,
            checked_input_columns,
            checked_output_rows,
            checked_output_columns,
        )
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("operator-space representation dimensions exceed typemax(Int)"))
    end
end

function OperatorSpace(
    input_shape::Union{Tuple,AbstractVector}, output_shape::Union{Tuple,AbstractVector}
)
    input_shape isa AbstractVector && Base.require_one_based_indexing(input_shape)
    output_shape isa AbstractVector && Base.require_one_based_indexing(output_shape)
    length(input_shape) == 2 ||
        throw(ArgumentError("input_size must contain exactly two dimensions"))
    length(output_shape) == 2 ||
        throw(ArgumentError("output_size must contain exactly two dimensions"))
    return OperatorSpace(input_shape[1], input_shape[2], output_shape[1], output_shape[2])
end

"""Return `(rows, columns)` for the input matrix space."""
input_size(space::OperatorSpace) = (space.input_rows, space.input_columns)

"""Return `(rows, columns)` for the output matrix space."""
output_size(space::OperatorSpace) = (space.output_rows, space.output_columns)

"""Return the common input dimension, rejecting nonsquare input spaces."""
function input_dimension(space::OperatorSpace)
    space.input_rows == space.input_columns || throw(
        ArgumentError(
            "input_dimension is defined only for square input spaces; " *
            "use input_size instead (got $(input_size(space)))",
        ),
    )
    return space.input_rows
end

"""Return the common output dimension, rejecting nonsquare output spaces."""
function output_dimension(space::OperatorSpace)
    space.output_rows == space.output_columns || throw(
        ArgumentError(
            "output_dimension is defined only for square output spaces; " *
            "use output_size instead (got $(output_size(space)))",
        ),
    )
    return space.output_rows
end

"""
    AbstractMapRepresentation{T}

Supertype for finite-dimensional linear-map representations with scalar
element type `T`. Every representation is associated with an
[`OperatorSpace`](@ref) descriptor. Use `input_size(map)` and
`output_size(map)` for rectangular operator spaces.

All representations in this package use column-major vectorization and the
unnormalized Choi convention
``J = sum(E_ij ⊗ Φ(E_ij), i, j)``. Thus
`vec(Φ(X)) == superoperator_matrix(map) * vec(X)`, and a Kraus operator `K`
contributes `vec(K) * vec(K)'` to `J`.
"""
abstract type AbstractMapRepresentation{T} end

"""
    KrausRepresentation(operators)

Represent the completely positive map
``Φ(X) = sum(K * X * K', K in operators)``.

`operators` must be a nonempty collection of finite numeric matrices with one
common `(output_dimension, input_dimension)` size. The matrices are copied;
sparse matrices remain sparse, and the stored collection does not permit
entry replacement. Use [`kraus_operators`](@ref) to obtain mutable copies.
Trace preservation is deliberately not imposed by the constructor.
"""
struct KrausRepresentation{T,V<:AbstractVector} <: AbstractMapRepresentation{T}
    operators::V
    input_dim::Int
    output_dim::Int

    function KrausRepresentation{T,V}(
        token::_ValidatedConstructorToken, operators::V, input_dim::Int, output_dim::Int
    ) where {T,V<:AbstractVector}
        _require_validated_constructor_token(token)
        input_dim > 0 && output_dim > 0 ||
            throw(ArgumentError("validated map dimensions must be positive"))
        isempty(operators) &&
            throw(ArgumentError("validated Kraus operators must be nonempty"))
        all(
            operator ->
                operator isa AbstractMatrix && size(operator) == (output_dim, input_dim),
            operators,
        ) || throw(ArgumentError("validated Kraus operator shapes are inconsistent"))
        return new{T,V}(operators, input_dim, output_dim)
    end
end

"""
    OperatorSumRepresentation(left_operators, right_operators)

Represent the general linear map
``Φ(X) = sum(A * X * B', (A, B) in zip(left_operators, right_operators))``.
For a map from `m × n` matrices to `p × q` matrices, every left factor has
size `p × m` and every right factor has size `q × n`.

The two collections must be nonempty and have equal length. Factors are
validated and copied; sparse factors remain sparse, and the stored collections
do not permit entry replacement. Unlike [`KrausRepresentation`](@ref), this
type does not imply complete positivity.
"""
struct OperatorSumRepresentation{T,VL<:AbstractVector,VR<:AbstractVector} <:
       AbstractMapRepresentation{T}
    left_operators::VL
    right_operators::VR
    space::OperatorSpace

    function OperatorSumRepresentation{T,VL,VR}(
        token::_ValidatedConstructorToken,
        left_operators::VL,
        right_operators::VR,
        space::OperatorSpace,
    ) where {T,VL<:AbstractVector,VR<:AbstractVector}
        _require_validated_constructor_token(token)
        isempty(left_operators) &&
            throw(ArgumentError("validated operator-sum factors must be nonempty"))
        length(left_operators) == length(right_operators) ||
            throw(ArgumentError("validated operator-sum factor counts must agree"))
        expected_left = (space.output_rows, space.input_rows)
        expected_right = (space.output_columns, space.input_columns)
        all(
            operator -> operator isa AbstractMatrix && size(operator) == expected_left,
            left_operators,
        ) || throw(ArgumentError("validated left-factor shapes are inconsistent"))
        all(
            operator -> operator isa AbstractMatrix && size(operator) == expected_right,
            right_operators,
        ) || throw(ArgumentError("validated right-factor shapes are inconsistent"))
        return new{T,VL,VR}(left_operators, right_operators, space)
    end
end

"""
    ChoiRepresentation(matrix, space)
    ChoiRepresentation(matrix, input_dim, output_dim)
    ChoiRepresentation(matrix)

Wrap a finite numeric Choi matrix using the input-first tensor ordering
``J[(i,a),(j,b)] = <a|Φ(|i><j|)|b>``. The explicit dimensions are required
for rectangular operator spaces. For a map from `m × n` to `p × q`, the
matrix has size `(m*p, n*q)`. The legacy `(input_dim, output_dim)` form
constructs the square-algebra space `(input_dim, input_dim) →
(output_dim, output_dim)`. The one-argument form infers an equal
input/output dimension from a square matrix of size `d^2`.

The matrix is copied, and sparse storage is preserved. Construction validates
shape and finite entries but does not require positivity, trace preservation,
or Hermiticity, so this type also represents nonphysical linear maps.
"""
struct ChoiRepresentation{T,M<:AbstractMatrix{T}} <: AbstractMapRepresentation{T}
    matrix::M
    space::OperatorSpace

    function ChoiRepresentation{T,M}(
        token::_ValidatedConstructorToken, matrix::M, space::OperatorSpace
    ) where {T,M<:AbstractMatrix{T}}
        _require_validated_constructor_token(token)
        expected = (
            _checked_product(
                (space.input_rows, space.output_rows), "validated Choi row dimensions"
            ),
            _checked_product(
                (space.input_columns, space.output_columns),
                "validated Choi column dimensions",
            ),
        )
        size(matrix) == expected ||
            throw(ArgumentError("validated Choi matrix shape is inconsistent"))
        return new{T,M}(matrix, space)
    end
end

"""
    SuperoperatorRepresentation(matrix, space)
    SuperoperatorRepresentation(matrix, input_dim, output_dim)
    SuperoperatorRepresentation(matrix)

Wrap a transfer matrix `S` satisfying `vec(Φ(X)) = S * vec(X)`, with Julia's
column-major `vec`. For a map from `m × n` to `p × q`, its size is
`(p*q, m*n)`. The legacy scalar-dimension form constructs a map between
square matrix algebras. The one-argument form infers both square-algebra
dimensions when the row and column counts are perfect squares.

The matrix is copied, and sparse storage is preserved. No physicality
condition is imposed.
"""
struct SuperoperatorRepresentation{T,M<:AbstractMatrix{T}} <: AbstractMapRepresentation{T}
    matrix::M
    space::OperatorSpace

    function SuperoperatorRepresentation{T,M}(
        token::_ValidatedConstructorToken, matrix::M, space::OperatorSpace
    ) where {T,M<:AbstractMatrix{T}}
        _require_validated_constructor_token(token)
        expected = (
            _checked_product(
                (space.output_rows, space.output_columns),
                "validated superoperator output dimensions",
            ),
            _checked_product(
                (space.input_rows, space.input_columns),
                "validated superoperator input dimensions",
            ),
        )
        size(matrix) == expected ||
            throw(ArgumentError("validated superoperator matrix shape is inconsistent"))
        return new{T,M}(matrix, space)
    end
end

# Preserve the legacy validated-constructor arity for invariant tests and
# internal square-algebra callers. The unforgeable token remains mandatory.
function ChoiRepresentation{T,M}(
    token::_ValidatedConstructorToken, matrix::M, input_dim::Int, output_dim::Int
) where {T,M<:AbstractMatrix{T}}
    _require_validated_constructor_token(token)
    return ChoiRepresentation{T,M}(
        token, matrix, OperatorSpace(input_dim, input_dim, output_dim, output_dim)
    )
end

function SuperoperatorRepresentation{T,M}(
    token::_ValidatedConstructorToken, matrix::M, input_dim::Int, output_dim::Int
) where {T,M<:AbstractMatrix{T}}
    _require_validated_constructor_token(token)
    return SuperoperatorRepresentation{T,M}(
        token, matrix, OperatorSpace(input_dim, input_dim, output_dim, output_dim)
    )
end

Base.eltype(::Type{<:AbstractMapRepresentation{T}}) where {T} = T
Base.eltype(map::AbstractMapRepresentation) = eltype(typeof(map))

"""
    operator_space(map)

Return the complete immutable [`OperatorSpace`](@ref) descriptor associated
with a map representation.
"""
function operator_space(map::KrausRepresentation)
    return OperatorSpace(map.input_dim, map.input_dim, map.output_dim, map.output_dim)
end
operator_space(map::OperatorSumRepresentation) = getfield(map, :space)
operator_space(map::ChoiRepresentation) = getfield(map, :space)
operator_space(map::SuperoperatorRepresentation) = getfield(map, :space)

input_size(map::AbstractMapRepresentation) = input_size(operator_space(map))
output_size(map::AbstractMapRepresentation) = output_size(operator_space(map))
input_dimension(map::AbstractMapRepresentation) = input_dimension(operator_space(map))
output_dimension(map::AbstractMapRepresentation) = output_dimension(operator_space(map))

function _has_square_operator_algebras(space::OperatorSpace)
    return space.input_rows == space.input_columns &&
           space.output_rows == space.output_columns
end

function _require_square_operator_algebras(
    map::AbstractMapRepresentation, operation::AbstractString
)
    space = operator_space(map)
    _has_square_operator_algebras(space) || throw(
        ArgumentError(
            "$operation is defined only for maps between square matrix algebras; " *
            "got $(input_size(space)) → $(output_size(space))",
        ),
    )
    return space.input_rows, space.output_rows
end

function Base.getproperty(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation}, name::Symbol
)
    name === :input_dim && return input_dimension(map)
    name === :output_dim && return output_dimension(map)
    return getfield(map, name)
end

function Base.propertynames(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation}, private::Bool=false
)
    return (fieldnames(typeof(map))..., :input_dim, :output_dim)
end

_all_finite(matrix::AbstractMatrix) = all(isfinite, matrix)
_all_finite(matrix::SparseMatrixCSC) = all(isfinite, nonzeros(matrix))

function _validate_numeric_matrix(matrix::AbstractMatrix, name::AbstractString)
    Base.require_one_based_indexing(matrix)
    eltype(matrix) <: Number || throw(
        ArgumentError("$name must have a numeric element type; got $(eltype(matrix))")
    )
    isconcretetype(eltype(matrix)) || throw(
        ArgumentError(
            "$name must have a concrete numeric element type; got $(eltype(matrix))"
        ),
    )
    _all_finite(matrix) || throw(ArgumentError("$name must contain only finite values"))
    return nothing
end

function _map_dimension(value, name::AbstractString)
    return _positive_int(value, name)
end

function _squared_dimension(dimension::Int, name::AbstractString)
    return _checked_power(dimension, 2, name)
end

function KrausRepresentation(operators)
    (operators isa Tuple || operators isa AbstractVector) ||
        throw(ArgumentError("operators must be a nonempty tuple or vector of matrices"))
    isempty(operators) &&
        throw(ArgumentError("operators must contain at least one Kraus operator"))
    all(operator -> operator isa AbstractMatrix, operators) ||
        throw(ArgumentError("every Kraus operator must be an AbstractMatrix"))

    first_operator = first(operators)
    _validate_numeric_matrix(first_operator, "operators[1]")
    output_dim, input_dim = size(first_operator)
    output_dim > 0 && input_dim > 0 ||
        throw(ArgumentError("Kraus operators must have nonzero dimensions"))

    for (index, operator) in enumerate(operators)
        _validate_numeric_matrix(operator, "operators[$index]")
        size(operator) == (output_dim, input_dim) || throw(
            DimensionMismatch(
                "operators[$index] has size $(size(operator)); expected " *
                "($output_dim, $input_dim)",
            ),
        )
    end

    scalar_type = foldl(promote_type, (eltype(operator) for operator in operators))
    copied = [copy(operator) for operator in operators]
    read_only = _read_only_plan_array(copied)
    return KrausRepresentation{scalar_type,typeof(read_only)}(
        _VALIDATED_CONSTRUCTOR_TOKEN, read_only, input_dim, output_dim
    )
end

function OperatorSumRepresentation(left_operators, right_operators)
    (left_operators isa Tuple || left_operators isa AbstractVector) || throw(
        ArgumentError("left_operators must be a nonempty tuple or vector of matrices")
    )
    (right_operators isa Tuple || right_operators isa AbstractVector) || throw(
        ArgumentError("right_operators must be a nonempty tuple or vector of matrices")
    )
    isempty(left_operators) &&
        throw(ArgumentError("left_operators must contain at least one matrix"))
    length(left_operators) == length(right_operators) || throw(
        DimensionMismatch(
            "left_operators and right_operators must contain the same number of factors"
        ),
    )
    all(operator -> operator isa AbstractMatrix, left_operators) ||
        throw(ArgumentError("every left operator must be an AbstractMatrix"))
    all(operator -> operator isa AbstractMatrix, right_operators) ||
        throw(ArgumentError("every right operator must be an AbstractMatrix"))

    first_left = first(left_operators)
    first_right = first(right_operators)
    _validate_numeric_matrix(first_left, "left_operators[1]")
    _validate_numeric_matrix(first_right, "right_operators[1]")
    output_rows, input_rows = size(first_left)
    output_columns, input_columns = size(first_right)
    space = OperatorSpace(input_rows, input_columns, output_rows, output_columns)

    for (index, operator) in enumerate(left_operators)
        _validate_numeric_matrix(operator, "left_operators[$index]")
        size(operator) == (output_rows, input_rows) || throw(
            DimensionMismatch(
                "left_operators[$index] has size $(size(operator)); expected " *
                "($output_rows, $input_rows)",
            ),
        )
    end
    for (index, operator) in enumerate(right_operators)
        _validate_numeric_matrix(operator, "right_operators[$index]")
        size(operator) == (output_columns, input_columns) || throw(
            DimensionMismatch(
                "right_operators[$index] has size $(size(operator)); expected " *
                "($output_columns, $input_columns)",
            ),
        )
    end

    scalar_type = foldl(
        promote_type,
        Iterators.flatten((
            (eltype(operator) for operator in left_operators),
            (eltype(operator) for operator in right_operators),
        )),
    )
    copied_left = _read_only_plan_array([copy(operator) for operator in left_operators])
    copied_right = _read_only_plan_array([copy(operator) for operator in right_operators])
    return OperatorSumRepresentation{scalar_type,typeof(copied_left),typeof(copied_right)}(
        _VALIDATED_CONSTRUCTOR_TOKEN, copied_left, copied_right, space
    )
end

"""
    operator_sum_factors(map)

Return `(left=[...], right=[...])`, with mutable copies of both factor
collections. Mutating the returned matrices does not change `map`.
"""
function operator_sum_factors(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return (
        left=[copy(operator) for operator in left_operators],
        right=[copy(operator) for operator in right_operators],
    )
end

function ChoiRepresentation(matrix::AbstractMatrix, space::OperatorSpace)
    _validate_numeric_matrix(matrix, "Choi matrix")
    expected = (
        _checked_product((space.input_rows, space.output_rows), "Choi row dimensions"),
        _checked_product(
            (space.input_columns, space.output_columns), "Choi column dimensions"
        ),
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "Choi matrix has size $(size(matrix)); expected $expected for " *
            "operator space $(input_size(space)) → $(output_size(space))",
        ),
    )
    copied = copy(matrix)
    return ChoiRepresentation{eltype(copied),typeof(copied)}(
        _VALIDATED_CONSTRUCTOR_TOKEN, copied, space
    )
end

function ChoiRepresentation(matrix::AbstractMatrix, input_dim, output_dim)
    checked_input = _map_dimension(input_dim, "input_dim")
    checked_output = _map_dimension(output_dim, "output_dim")
    space = OperatorSpace(checked_input, checked_input, checked_output, checked_output)
    return ChoiRepresentation(matrix, space)
end

function ChoiRepresentation(matrix::AbstractMatrix, input_shape::Tuple, output_shape::Tuple)
    return ChoiRepresentation(matrix, OperatorSpace(input_shape, output_shape))
end

function ChoiRepresentation(matrix::AbstractMatrix; input_dim=nothing, output_dim=nothing)
    if input_dim === nothing && output_dim === nothing
        size(matrix, 1) == size(matrix, 2) ||
            throw(DimensionMismatch("a Choi matrix must be square"))
        dimension = isqrt(size(matrix, 1))
        dimension > 0 && dimension^2 == size(matrix, 1) || throw(
            ArgumentError(
                "cannot infer equal input/output dimensions from Choi size $(size(matrix))",
            ),
        )
        return ChoiRepresentation(matrix, dimension, dimension)
    elseif input_dim === nothing || output_dim === nothing
        throw(ArgumentError("input_dim and output_dim must be supplied together"))
    end
    return ChoiRepresentation(matrix, input_dim, output_dim)
end

function SuperoperatorRepresentation(matrix::AbstractMatrix, space::OperatorSpace)
    _validate_numeric_matrix(matrix, "superoperator matrix")
    expected = (
        _checked_product(
            (space.output_rows, space.output_columns), "superoperator output dimensions"
        ),
        _checked_product(
            (space.input_rows, space.input_columns), "superoperator input dimensions"
        ),
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "superoperator matrix has size $(size(matrix)); expected $expected for " *
            "operator space $(input_size(space)) → $(output_size(space))",
        ),
    )
    copied = copy(matrix)
    return SuperoperatorRepresentation{eltype(copied),typeof(copied)}(
        _VALIDATED_CONSTRUCTOR_TOKEN, copied, space
    )
end

function SuperoperatorRepresentation(matrix::AbstractMatrix, input_dim, output_dim)
    checked_input = _map_dimension(input_dim, "input_dim")
    checked_output = _map_dimension(output_dim, "output_dim")
    space = OperatorSpace(checked_input, checked_input, checked_output, checked_output)
    return SuperoperatorRepresentation(matrix, space)
end

function SuperoperatorRepresentation(
    matrix::AbstractMatrix, input_shape::Tuple, output_shape::Tuple
)
    return SuperoperatorRepresentation(matrix, OperatorSpace(input_shape, output_shape))
end

function SuperoperatorRepresentation(
    matrix::AbstractMatrix; input_dim=nothing, output_dim=nothing
)
    if input_dim === nothing && output_dim === nothing
        inferred_output = isqrt(size(matrix, 1))
        inferred_input = isqrt(size(matrix, 2))
        inferred_output > 0 &&
        inferred_input > 0 &&
        inferred_output^2 == size(matrix, 1) &&
        inferred_input^2 == size(matrix, 2) || throw(
            ArgumentError(
                "superoperator row and column counts must be nonzero perfect squares; " *
                "got size $(size(matrix))",
            ),
        )
        return SuperoperatorRepresentation(matrix, inferred_input, inferred_output)
    elseif input_dim === nothing || output_dim === nothing
        throw(ArgumentError("input_dim and output_dim must be supplied together"))
    end
    return SuperoperatorRepresentation(matrix, input_dim, output_dim)
end

function _validated_representation_matrix(map::ChoiRepresentation)
    matrix = getfield(map, :matrix)
    Base.require_one_based_indexing(matrix)
    space = operator_space(map)
    expected = (
        _checked_product(
            (space.input_rows, space.output_rows), "stored Choi row dimensions"
        ),
        _checked_product(
            (space.input_columns, space.output_columns), "stored Choi column dimensions"
        ),
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "stored Choi matrix has size $(size(matrix)); expected $expected for " *
            "operator space $(input_size(space)) → $(output_size(space))",
        ),
    )
    _all_finite(matrix) ||
        throw(ArgumentError("stored Choi matrix must contain only finite values"))
    return matrix
end

function _validated_representation_matrix(map::SuperoperatorRepresentation)
    matrix = getfield(map, :matrix)
    Base.require_one_based_indexing(matrix)
    space = operator_space(map)
    expected = (
        _checked_product(
            (space.output_rows, space.output_columns),
            "stored superoperator output dimensions",
        ),
        _checked_product(
            (space.input_rows, space.input_columns), "stored superoperator input dimensions"
        ),
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "stored superoperator matrix has size $(size(matrix)); expected $expected " *
            "for operator space $(input_size(space)) → $(output_size(space))",
        ),
    )
    _all_finite(matrix) ||
        throw(ArgumentError("stored superoperator matrix must contain only finite values"))
    return matrix
end

function _validated_kraus_operators(map::KrausRepresentation)
    operators = getfield(map, :operators)
    for (index, operator) in enumerate(operators)
        Base.require_one_based_indexing(operator)
        size(operator) == (map.output_dim, map.input_dim) || throw(
            DimensionMismatch(
                "stored Kraus operator $index has size $(size(operator)); expected " *
                "($(map.output_dim), $(map.input_dim))",
            ),
        )
        _all_finite(operator) || throw(
            ArgumentError("stored Kraus operator $index must contain only finite values"),
        )
    end
    return operators
end

function _validated_operator_sum_factors(map::OperatorSumRepresentation)
    left_operators = getfield(map, :left_operators)
    right_operators = getfield(map, :right_operators)
    length(left_operators) == length(right_operators) ||
        throw(DimensionMismatch("stored operator-sum factor counts do not agree"))
    isempty(left_operators) &&
        throw(ArgumentError("stored operator-sum factors must be nonempty"))
    space = operator_space(map)
    expected_left = (space.output_rows, space.input_rows)
    expected_right = (space.output_columns, space.input_columns)
    for (index, operator) in enumerate(left_operators)
        Base.require_one_based_indexing(operator)
        size(operator) == expected_left || throw(
            DimensionMismatch(
                "stored left operator $index has size $(size(operator)); " *
                "expected $expected_left",
            ),
        )
        _all_finite(operator) || throw(
            ArgumentError("stored left operator $index must contain only finite values")
        )
    end
    for (index, operator) in enumerate(right_operators)
        Base.require_one_based_indexing(operator)
        size(operator) == expected_right || throw(
            DimensionMismatch(
                "stored right operator $index has size $(size(operator)); " *
                "expected $expected_right",
            ),
        )
        _all_finite(operator) || throw(
            ArgumentError("stored right operator $index must contain only finite values"),
        )
    end
    return left_operators, right_operators
end

function _choi_to_superoperator(matrix::SparseMatrixCSC, space::OperatorSpace)
    input_rows, input_columns = input_size(space)
    output_rows, output_columns = output_size(space)
    rows, columns, values = findnz(matrix)
    super_rows = Vector{Int}(undef, length(rows))
    super_columns = Vector{Int}(undef, length(rows))
    @inbounds for index in eachindex(rows)
        output_row = mod(rows[index] - 1, output_rows) + 1
        input_row = div(rows[index] - 1, output_rows) + 1
        output_column = mod(columns[index] - 1, output_columns) + 1
        input_column = div(columns[index] - 1, output_columns) + 1
        super_rows[index] = output_row + (output_column - 1) * output_rows
        super_columns[index] = input_row + (input_column - 1) * input_rows
    end
    return sparse(
        super_rows,
        super_columns,
        copy(values),
        output_rows * output_columns,
        input_rows * input_columns,
    )
end

function _choi_to_superoperator(matrix::AbstractMatrix, space::OperatorSpace)
    input_rows, input_columns = input_size(space)
    output_rows, output_columns = output_size(space)
    result = Matrix{eltype(matrix)}(
        undef, output_rows * output_columns, input_rows * input_columns
    )
    @inbounds for input_column in 1:input_columns,
        input_row in 1:input_rows, output_column in 1:output_columns,
        output_row in 1:output_rows

        result[output_row + (output_column - 1) * output_rows, input_row + (input_column - 1) * input_rows] = matrix[
            output_row + (input_row - 1) * output_rows,
            output_column + (input_column - 1) * output_columns,
        ]
    end
    return result
end

function _superoperator_to_choi(matrix::SparseMatrixCSC, space::OperatorSpace)
    input_rows, input_columns = input_size(space)
    output_rows, output_columns = output_size(space)
    rows, columns, values = findnz(matrix)
    choi_rows = Vector{Int}(undef, length(rows))
    choi_columns = Vector{Int}(undef, length(rows))
    @inbounds for index in eachindex(rows)
        output_row = mod(rows[index] - 1, output_rows) + 1
        output_column = div(rows[index] - 1, output_rows) + 1
        input_row = mod(columns[index] - 1, input_rows) + 1
        input_column = div(columns[index] - 1, input_rows) + 1
        choi_rows[index] = output_row + (input_row - 1) * output_rows
        choi_columns[index] = output_column + (input_column - 1) * output_columns
    end
    return sparse(
        choi_rows,
        choi_columns,
        copy(values),
        input_rows * output_rows,
        input_columns * output_columns,
    )
end

function _superoperator_to_choi(matrix::AbstractMatrix, space::OperatorSpace)
    input_rows, input_columns = input_size(space)
    output_rows, output_columns = output_size(space)
    result = Matrix{eltype(matrix)}(
        undef, input_rows * output_rows, input_columns * output_columns
    )
    @inbounds for input_column in 1:input_columns,
        input_row in 1:input_rows, output_column in 1:output_columns,
        output_row in 1:output_rows

        result[output_row + (input_row - 1) * output_rows, output_column + (input_column - 1) * output_columns] = matrix[
            output_row + (output_column - 1) * output_rows,
            input_row + (input_column - 1) * input_rows,
        ]
    end
    return result
end

_column_vector(operator::AbstractMatrix) = vec(operator)

function _column_vector(operator::SparseMatrixCSC)
    rows, columns, values = findnz(operator)
    indices = rows .+ (columns .- 1) .* size(operator, 1)
    return sparsevec(indices, copy(values), length(operator))
end

function _kraus_to_choi(operators)
    terms = Base.map(operators) do operator
        vector = _column_vector(operator)
        return vector * adjoint(vector)
    end
    return reduce(+, terms)
end

function _kraus_to_superoperator(operators)
    terms = Base.map(operator -> kron(conj(operator), operator), operators)
    return reduce(+, terms)
end

function _operator_sum_to_choi(left_operators, right_operators)
    terms = Base.map(left_operators, right_operators) do left, right
        return _column_vector(left) * adjoint(_column_vector(right))
    end
    return reduce(+, terms)
end

function _operator_sum_to_superoperator(left_operators, right_operators)
    terms = Base.map(left_operators, right_operators) do left, right
        return kron(conj(right), left)
    end
    return reduce(+, terms)
end

"""
    choi_representation(map)
    superoperator_representation(map)
    kraus_representation(map; atol, rtol)

Convert between map representations. Choi/superoperator conversion is an
exact index reshuffle and preserves sparse storage. Kraus-to-matrix conversion
preserves sparse storage when every Kraus operator is sparse.

Choi-to-Kraus conversion is defined only when complete positivity is
established by the structured diagnostic. A numerical boundary is rejected,
not projected onto the positive semidefinite cone. The conversion performs a
dense Hermitian eigendecomposition, so sparse Choi data require
`allow_densify=true`. Use [`canonical_map_decomposition`](@ref) when a stable
paired result is required for Hermiticity-preserving, general, or
numerically-inconclusive maps.
"""
function choi_representation(map::ChoiRepresentation)
    return ChoiRepresentation(_validated_representation_matrix(map), operator_space(map))
end

function choi_representation(map::KrausRepresentation)
    operators = _validated_kraus_operators(map)
    return ChoiRepresentation(_kraus_to_choi(operators), operator_space(map))
end

function choi_representation(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return ChoiRepresentation(
        _operator_sum_to_choi(left_operators, right_operators), operator_space(map)
    )
end

function choi_representation(map::SuperoperatorRepresentation)
    return ChoiRepresentation(
        _superoperator_to_choi(_validated_representation_matrix(map), operator_space(map)),
        operator_space(map),
    )
end

function superoperator_representation(map::SuperoperatorRepresentation)
    return SuperoperatorRepresentation(
        _validated_representation_matrix(map), operator_space(map)
    )
end

function superoperator_representation(map::KrausRepresentation)
    operators = _validated_kraus_operators(map)
    return SuperoperatorRepresentation(
        _kraus_to_superoperator(operators), operator_space(map)
    )
end

function superoperator_representation(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return SuperoperatorRepresentation(
        _operator_sum_to_superoperator(left_operators, right_operators), operator_space(map)
    )
end

function superoperator_representation(map::ChoiRepresentation)
    return SuperoperatorRepresentation(
        _choi_to_superoperator(_validated_representation_matrix(map), operator_space(map)),
        operator_space(map),
    )
end

function kraus_representation(
    map::KrausRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    _checked_tolerances(map, atol, rtol)
    return KrausRepresentation(_validated_kraus_operators(map))
end

function _real_float_type(::Type{T}) where {T}
    return typeof(float(real(zero(T))))
end

function _default_rtol(map::AbstractMapRepresentation)
    return sqrt(eps(_real_float_type(eltype(map))))
end

function _checked_tolerances(map, atol, rtol)
    !(atol isa Bool) && atol isa Real && isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("atol must be a finite nonnegative real number"))
    !(rtol isa Bool) && rtol isa Real && isfinite(rtol) && rtol >= 0 ||
        throw(ArgumentError("rtol must be a finite nonnegative real number"))
    return atol, rtol
end

"""
    OperatorSumDecompositionResult

Diagnostics returned by [`operator_sum_decomposition`](@ref). `threshold` is
the singular-value cutoff `atol + rtol * max(σ₁, 1)`;
`discarded_frobenius_norm` is the Euclidean norm of the discarded singular
values. The represented map therefore differs from the input Choi matrix by
that Frobenius norm, up to the numerical SVD error.
"""
struct OperatorSumDecompositionResult{R,V,T,N}
    representation::R
    singular_values::V
    threshold::T
    retained_rank::Int
    discarded_frobenius_norm::N
end

"""
    CanonicalMapDecompositionResult

Structured output of [`canonical_map_decomposition`](@ref). `representation`
is always an [`OperatorSumRepresentation`](@ref), so paired factors have a
stable native type even when the map is not completely positive.
`classification` is one of:

- `:completely_positive`, when complete positivity is established;
- `:hermiticity_preserving`, when Choi Hermiticity is established but
  complete positivity is violated or numerically inconclusive;
- `:general`, when an SVD is required without symmetrizing the Choi matrix.

`spectral_values`, `threshold`, `retained_rank`, and
`discarded_frobenius_norm` describe the eigenvalue or singular-value cutoff.
`reconstruction_frobenius_norm` measures the actual Choi reconstruction
residual. `complete_positivity` and `hermiticity_preservation` retain the
structured diagnostics used to select the branch; they are `nothing` when
those properties are not applicable to a genuinely rectangular operator
space.
"""
struct CanonicalMapDecompositionResult{R,V,T,N,C,H}
    representation::R
    classification::Symbol
    spectral_values::V
    threshold::T
    retained_rank::Int
    discarded_frobenius_norm::N
    reconstruction_frobenius_norm::N
    complete_positivity::C
    hermiticity_preservation::H
end

"""
    operator_sum_representation(map; kwargs...)

Return a general paired-factor representation. Conversion from an existing
operator sum or a Kraus representation is exact. Conversion from a Choi or
superoperator representation uses [`operator_sum_decomposition`](@ref) and
accepts its tolerance and densification keywords.
"""
function operator_sum_representation(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return OperatorSumRepresentation(left_operators, right_operators)
end

function operator_sum_representation(map::KrausRepresentation)
    operators = _validated_kraus_operators(map)
    return OperatorSumRepresentation(operators, operators)
end

"""
    operator_sum_decomposition(map; atol=0, rtol=sqrt(eps(T)),
                               allow_densify=false)

Compute a general paired-factor decomposition of a Choi or superoperator
representation using a thin SVD. If
`J = U * Diagonal(σ) * V'`, each retained singular triplet yields
`A = reshape(sqrt(σ) * U[:,i], output_rows, input_rows)` and
`B = reshape(sqrt(σ) * V[:,i], output_columns, input_columns)`.

Only `Float32`, `Float64`, and their complex counterparts are accepted,
matching the SVD scalar types supported by Julia's standard linear algebra.
Exact and arbitrary-precision inputs are rejected rather than silently
downcast. Sparse input is rejected unless `allow_densify=true`. Singular
values at or below `atol + rtol * max(σ₁, 1)` are discarded and reported in
the returned [`OperatorSumDecompositionResult`](@ref). A zero map is
represented by one explicitly zero factor pair.
"""
function operator_sum_decomposition(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation};
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    choi = map isa ChoiRepresentation ? map : choi_representation(map)
    matrix = _validated_representation_matrix(choi)
    issparse(matrix) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "operator-sum SVD would densify sparse Choi data; " *
                "pass allow_densify=true to permit this",
            ),
        )
    eltype(matrix) <: LinearAlgebra.BlasFloat || throw(
        ArgumentError(
            "operator-sum SVD requires Float32, Float64, ComplexF32, or " *
            "ComplexF64 data; got $(eltype(matrix))",
        ),
    )

    decomposition = svd(Matrix(matrix); full=false)
    singular_values = copy(decomposition.S)
    scale = max(first(singular_values), one(eltype(singular_values)))
    threshold = atol + rtol * scale
    isfinite(threshold) || throw(
        ArgumentError("atol + rtol * max(σ₁, 1) must be finite; got threshold=$threshold"),
    )
    retained = findall(>(threshold), singular_values)
    discarded = findall(<=(threshold), singular_values)
    discarded_norm = norm(view(singular_values, discarded))

    space = operator_space(choi)
    left_operators = Base.map(retained) do index
        return copy(
            reshape(
                sqrt(singular_values[index]) * decomposition.U[:, index],
                space.output_rows,
                space.input_rows,
            ),
        )
    end
    right_operators = Base.map(retained) do index
        return copy(
            reshape(
                sqrt(singular_values[index]) * decomposition.V[:, index],
                space.output_columns,
                space.input_columns,
            ),
        )
    end
    if isempty(retained)
        scalar_type = eltype(decomposition.U)
        push!(left_operators, zeros(scalar_type, space.output_rows, space.input_rows))
        push!(
            right_operators, zeros(scalar_type, space.output_columns, space.input_columns)
        )
    end

    representation = OperatorSumRepresentation(left_operators, right_operators)
    return OperatorSumDecompositionResult(
        representation, singular_values, threshold, length(retained), discarded_norm
    )
end

function operator_sum_representation(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation}; kwargs...
)
    return operator_sum_decomposition(map; kwargs...).representation
end

function _decomposition_threshold(values, atol, rtol)
    scale = max(maximum(abs, values), one(eltype(values)))
    threshold = atol + rtol * scale
    isfinite(threshold) || throw(
        ArgumentError(
            "atol + rtol * max(maximum(abs, values), 1) must be finite; " *
            "got threshold=$threshold",
        ),
    )
    return threshold
end

function _zero_operator_sum(space::OperatorSpace, ::Type{T}) where {T}
    return OperatorSumRepresentation(
        [zeros(T, space.output_rows, space.input_rows)],
        [zeros(T, space.output_columns, space.input_columns)],
    )
end

function _canonical_hermitian_operator_sum(
    matrix, space::OperatorSpace, classification::Symbol, atol, rtol
)
    decomposition = eigen(Hermitian(Matrix(matrix)))
    order = reverse(eachindex(decomposition.values))
    spectral_values = copy(decomposition.values[order])
    vectors = decomposition.vectors[:, order]
    threshold = _decomposition_threshold(spectral_values, atol, rtol)
    retained = if classification === :completely_positive
        findall(>(threshold), spectral_values)
    else
        findall(value -> abs(value) > threshold, spectral_values)
    end
    discarded = if classification === :completely_positive
        findall(<=(threshold), spectral_values)
    else
        findall(value -> abs(value) <= threshold, spectral_values)
    end
    discarded_norm = norm(view(spectral_values, discarded))

    left_operators = Base.map(retained) do index
        vector = sqrt(abs(spectral_values[index])) * vectors[:, index]
        return copy(reshape(vector, space.output_rows, space.input_rows))
    end
    right_operators = Base.map(retained) do index
        vector = sqrt(abs(spectral_values[index])) * vectors[:, index]
        signed_vector = if classification === :completely_positive
            vector
        else
            sign(spectral_values[index]) * vector
        end
        return copy(reshape(signed_vector, space.output_columns, space.input_columns))
    end
    representation = if isempty(retained)
        _zero_operator_sum(space, eltype(vectors))
    else
        OperatorSumRepresentation(left_operators, right_operators)
    end
    return representation, spectral_values, threshold, length(retained), discarded_norm
end

"""
    canonical_map_decomposition(map; atol=0, rtol=sqrt(eps(T)),
                                allow_densify=false)
        -> CanonicalMapDecompositionResult

Compute QETLAB's three canonical factor branches without calling a general
two-sided factorization “Kraus”:

1. a map whose complete positivity is established uses a Choi
   eigendecomposition with equal left and right factors;
2. an exactly Choi-Hermitian map uses positive factor pairs first and signed
   negative pairs second;
3. every other map uses a thin SVD of the unmodified Choi matrix.

The returned representation always has stable paired-factor type. Numerical
boundary outcomes never enter the completely-positive branch. No Choi matrix
is symmetrized, clipped, or projected. A nonzero cutoff is an explicit
approximation and both its discarded Frobenius norm and actual reconstruction
residual are reported.

Canonical spectral work supports `Float32`, `Float64`, and their complex
counterparts. Exact and arbitrary-precision inputs are rejected rather than
downcast. Sparse Choi data require `allow_densify=true`.
"""
function canonical_map_decomposition(
    map::AbstractMapRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    choi = choi_representation(map)
    matrix = _validated_representation_matrix(choi)
    issparse(matrix) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "canonical map decomposition would densify sparse Choi data; " *
                "pass allow_densify=true to permit this",
            ),
        )
    eltype(matrix) <: LinearAlgebra.BlasFloat || throw(
        ArgumentError(
            "canonical map decomposition requires Float32, Float64, ComplexF32, " *
            "or ComplexF64 data; got $(eltype(matrix))",
        ),
    )

    space = operator_space(choi)
    square_algebras = _has_square_operator_algebras(space)
    hermiticity = if square_algebras
        is_hermiticity_preserving(choi; atol=atol, rtol=rtol)
    else
        nothing
    end
    complete_positivity = if map isa KrausRepresentation
        is_completely_positive(map; atol=atol, rtol=rtol)
    elseif square_algebras
        is_completely_positive(choi; atol=atol, rtol=rtol, allow_densify=true)
    else
        nothing
    end

    classification =
        if complete_positivity !== nothing &&
            complete_positivity.status === MatrixPredicateSatisfied
            :completely_positive
        elseif hermiticity !== nothing && hermiticity.status === MatrixPredicateSatisfied
            :hermiticity_preserving
        else
            :general
        end

    representation, spectral_values, threshold, retained_rank, discarded_norm =
        if classification === :general
            general = operator_sum_decomposition(
                choi; atol=atol, rtol=rtol, allow_densify=true
            )
            (
                general.representation,
                general.singular_values,
                general.threshold,
                general.retained_rank,
                general.discarded_frobenius_norm,
            )
        else
            _canonical_hermitian_operator_sum(matrix, space, classification, atol, rtol)
        end
    reconstruction_residual = norm(Matrix(matrix) - Matrix(choi_matrix(representation)))
    return CanonicalMapDecompositionResult(
        representation,
        classification,
        spectral_values,
        threshold,
        retained_rank,
        discarded_norm,
        reconstruction_residual,
        complete_positivity,
        hermiticity,
    )
end

function kraus_representation(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation};
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    result = canonical_map_decomposition(
        map; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    result.complete_positivity !== nothing &&
    result.complete_positivity.status === MatrixPredicateSatisfied || throw(
        DomainError(
            result.complete_positivity,
            "complete positivity was not established; inspect " *
            "canonical_map_decomposition for a paired reconstruction",
        ),
    )
    factors = operator_sum_factors(result.representation)
    return KrausRepresentation(factors.left)
end

function kraus_representation(
    map::OperatorSumRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    _require_square_operator_algebras(map, "kraus_representation")
    return kraus_representation(
        choi_representation(map); atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

"""
    choi_matrix(map)

Return a mutable copy of the ordinary-array Choi matrix, converting the map
representation when necessary.
"""
choi_matrix(map::AbstractMapRepresentation) = copy(choi_representation(map).matrix)

"""
    superoperator_matrix(map)

Return a mutable copy of the ordinary-array superoperator matrix, converting
the map representation when necessary.
"""
function superoperator_matrix(map::AbstractMapRepresentation)
    return copy(superoperator_representation(map).matrix)
end

"""
    kraus_operators(map; atol, rtol)

Return mutable copies of the map's Kraus operators. Conversion follows the
explicit dense eigendecomposition and tolerance behavior of
`kraus_representation`.
"""
function kraus_operators(map::AbstractMapRepresentation; kwargs...)
    return [copy(operator) for operator in kraus_representation(map; kwargs...).operators]
end

function _validate_channel_input(input::AbstractMatrix, map::AbstractMapRepresentation)
    _validate_numeric_matrix(input, "channel input")
    expected = input_size(map)
    size(input) == expected || throw(
        DimensionMismatch("channel input has size $(size(input)); expected $expected")
    )
    return nothing
end

"""
    apply_channel(input, map)
    apply_channel(input, kraus_operators)

Apply a finite-dimensional map to a matrix of size `input_size(map)`. `map`
may be a `KrausRepresentation`, `OperatorSumRepresentation`,
`ChoiRepresentation`, or `SuperoperatorRepresentation`. A vector/tuple of
matrices is interpreted as Kraus operators.

Kraus application evaluates `sum(K * input * K')`. Operator-sum application
evaluates `sum(A * input * B')` directly without constructing a
superoperator. Choi application uses the exact Choi-to-superoperator
reshuffle. Sparse output is retained when sparse linear algebra permits it;
no input is normalized, symmetrized, or repaired.
"""
function apply_channel(input::AbstractMatrix, map::KrausRepresentation)
    _validate_channel_input(input, map)
    operators = _validated_kraus_operators(map)
    terms = Base.map(operators) do operator
        return operator * input * adjoint(operator)
    end
    return reduce(+, terms)
end

function _can_use_strided_operator_sum_kernel(input, left_operators, right_operators)
    scalar_type = promote_type(
        eltype(input),
        (eltype(operator) for operator in left_operators)...,
        (eltype(operator) for operator in right_operators)...,
    )
    compatible =
        input isa StridedMatrix{scalar_type} &&
        all(operator -> operator isa StridedMatrix{scalar_type}, left_operators) &&
        all(operator -> operator isa StridedMatrix{scalar_type}, right_operators)
    return compatible, scalar_type
end

function _apply_strided_operator_sum(
    input::StridedMatrix{T}, left_operators, right_operators
) where {T}
    output_rows = size(first(left_operators), 1)
    output_columns = size(first(right_operators), 1)
    input_columns = size(input, 2)
    intermediate = Matrix{T}(undef, output_rows, input_columns)
    result = Matrix{T}(undef, output_rows, output_columns)

    mul!(intermediate, first(left_operators), input)
    mul!(result, intermediate, adjoint(first(right_operators)))
    length(left_operators) == 1 && return result

    for index in 2:length(left_operators)
        mul!(intermediate, left_operators[index], input)
        mul!(result, intermediate, adjoint(right_operators[index]), one(T), one(T))
    end
    return result
end

function _apply_generic_operator_sum(input, left_operators, right_operators)
    result = first(left_operators) * input * adjoint(first(right_operators))
    for index in 2:length(left_operators)
        result = result + left_operators[index] * input * adjoint(right_operators[index])
    end
    return result
end

function apply_channel(input::AbstractMatrix, map::OperatorSumRepresentation)
    _validate_channel_input(input, map)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    use_strided_kernel, scalar_type = _can_use_strided_operator_sum_kernel(
        input, left_operators, right_operators
    )
    if use_strided_kernel
        return _apply_strided_operator_sum(
            input::StridedMatrix{scalar_type}, left_operators, right_operators
        )
    end
    return _apply_generic_operator_sum(input, left_operators, right_operators)
end

function apply_channel(input::AbstractMatrix, map::SuperoperatorRepresentation)
    _validate_channel_input(input, map)
    output_vector = _validated_representation_matrix(map) * vec(input)
    return reshape(output_vector, output_size(map))
end

function apply_channel(input::AbstractMatrix, map::ChoiRepresentation)
    _validate_channel_input(input, map)
    return apply_channel(input, superoperator_representation(map))
end

function apply_channel(input::AbstractMatrix, operators::Union{Tuple,AbstractVector})
    return apply_channel(input, KrausRepresentation(operators))
end

function _map_hermiticity_scale(matrix, ::Type{R}) where {R}
    scale = one(R)
    values = matrix isa SparseMatrixCSC ? nonzeros(matrix) : matrix
    for value in values
        scale = max(scale, abs(value))
    end
    return scale
end

function _maximum_hermiticity_defect(matrix::SparseMatrixCSC)
    difference = matrix - adjoint(matrix)
    rows, columns, values = findnz(difference)
    isempty(values) &&
        return zero(typeof(abs(zero(eltype(matrix))))), (1, 1), zero(eltype(matrix))
    selected = firstindex(values)
    residual = abs(values[selected])
    for index in Iterators.drop(eachindex(values), 1)
        candidate = abs(values[index])
        if candidate > residual
            residual = candidate
            selected = index
        end
    end
    return residual, (rows[selected], columns[selected]), values[selected]
end

function _maximum_hermiticity_defect(matrix::AbstractMatrix)
    selected = (1, 1)
    difference = matrix[1, 1] - conj(matrix[1, 1])
    residual = abs(difference)
    for column in axes(matrix, 2), row in first(axes(matrix, 1)):column
        candidate_difference = matrix[row, column] - conj(matrix[column, row])
        candidate = abs(candidate_difference)
        if candidate > residual
            residual = candidate
            selected = (row, column)
            difference = candidate_difference
        end
    end
    return residual, selected, difference
end

function _choi_hermiticity_result(
    matrix::AbstractMatrix; atol=nothing, rtol=nothing, nonsquare_status::Symbol=:reject
)
    _validate_numeric_matrix(matrix, "Choi matrix")
    _, real_type, exact = _matrix_predicate_type_information(matrix)
    absolute, relative = _matrix_predicate_tolerances(
        real_type, exact; atol=atol, rtol=rtol
    )
    scale = _map_hermiticity_scale(matrix, real_type)
    tolerance = _matrix_predicate_threshold(scale, absolute, relative)
    isfinite(tolerance) ||
        throw(ArgumentError("the combined Hermiticity tolerance must be finite"))

    if size(matrix, 1) != size(matrix, 2)
        nonsquare_status === :reject && throw(
            ArgumentError(
                "Hermiticity preservation is defined only for maps between " *
                "square matrix algebras; the Choi matrix has size $(size(matrix))",
            ),
        )
        nonsquare_status === :violated ||
            throw(ArgumentError("nonsquare_status must be :reject or :violated"))
        return _matrix_predicate_result(
            :hermiticity_preserving,
            MatrixPredicateViolated;
            value=nothing,
            tolerance=tolerance,
            witness=(kind=:nonsquare_choi, size=size(matrix)),
            message="the Choi matrix is nonsquare, so the QETLAB-compatible Choi Hermiticity test is false",
        )
    end

    residual, indices, difference = _maximum_hermiticity_defect(matrix)
    status = if exact
        iszero(residual) ? MatrixPredicateSatisfied : MatrixPredicateViolated
    elseif iszero(residual)
        MatrixPredicateSatisfied
    elseif residual > tolerance
        MatrixPredicateViolated
    else
        MatrixPredicateUnknown
    end
    witness = if status === MatrixPredicateSatisfied
        nothing
    else
        (
            kind=(
                if status === MatrixPredicateUnknown
                    :hermiticity_boundary
                else
                    :nonhermitian
                end
            ),
            row=indices[1],
            column=indices[2],
            difference=difference,
            residual=residual,
        )
    end
    message = if status === MatrixPredicateSatisfied
        "the Choi matrix is exactly Hermitian"
    elseif status === MatrixPredicateViolated
        "the Choi matrix is non-Hermitian outside the requested tolerance"
    else
        "the nonzero Choi Hermiticity defect lies inside the numerical tolerance boundary; the map was not repaired"
    end
    return _matrix_predicate_result(
        :hermiticity_preserving,
        status;
        value=residual,
        tolerance=tolerance,
        witness=witness,
        message=message,
    )
end

"""
    is_hermiticity_preserving(map; atol=nothing, rtol=nothing)
        -> MatrixPredicateResult

Test Hermiticity preservation through the map's Choi matrix. This property is
defined only for maps between square matrix algebras; genuinely rectangular
operator spaces are rejected explicitly.

Exact integer and rational maps are decided exactly. For floating-point maps,
an exactly zero entrywise Choi Hermiticity defect returns
`MatrixPredicateSatisfied`, a defect above
`atol + rtol * max(maximum(abs, J), 1)` returns
`MatrixPredicateViolated`, and a nonzero defect inside that boundary returns
`MatrixPredicateUnknown`. The Choi matrix is never symmetrized or repaired,
and sparse representations are inspected without dense factorization.
"""
function is_hermiticity_preserving(
    map::AbstractMapRepresentation; atol=nothing, rtol=nothing
)
    _require_square_operator_algebras(map, "is_hermiticity_preserving")
    choi = choi_representation(map)
    return _choi_hermiticity_result(
        _validated_representation_matrix(choi); atol=atol, rtol=rtol
    )
end

function _complete_positivity_result_from_choi(
    matrix::AbstractMatrix;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    nonsquare_status::Symbol=:reject,
)
    hermiticity = _choi_hermiticity_result(
        matrix; atol=atol, rtol=rtol, nonsquare_status=nonsquare_status
    )
    if hermiticity.status !== MatrixPredicateSatisfied
        message = if hermiticity.status === MatrixPredicateViolated
            "complete positivity is violated because the Choi matrix is not Hermitian"
        else
            "complete positivity is unknown because Choi Hermiticity lies on a numerical boundary"
        end
        return _matrix_predicate_result(
            :completely_positive,
            hermiticity.status;
            value=hermiticity.value,
            tolerance=hermiticity.tolerance,
            witness=(kind=:choi_hermiticity, result=hermiticity),
            message=message,
        )
    end

    semidefiniteness = is_positive_semidefinite(
        matrix; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    message = if semidefiniteness.status === MatrixPredicateSatisfied
        "the Choi matrix is positive semidefinite"
    elseif semidefiniteness.status === MatrixPredicateViolated
        "the Choi matrix has a negative-semidefiniteness witness"
    else
        "Choi positive semidefiniteness lies on a numerical boundary"
    end
    return _matrix_predicate_result(
        :completely_positive,
        semidefiniteness.status;
        value=semidefiniteness.value,
        tolerance=semidefiniteness.tolerance,
        witness=(
            if semidefiniteness.status === MatrixPredicateSatisfied
                nothing
            else
                (kind=:choi_positive_semidefiniteness, result=semidefiniteness)
            end
        ),
        message=message,
    )
end

function _kraus_complete_positivity_result(map, atol, rtol)
    value_type = eltype(map)
    real_type = typeof(real(zero(value_type)))
    exact = real_type <: Integer || real_type <: Rational
    absolute, relative = _matrix_predicate_tolerances(
        real_type, exact; atol=atol, rtol=rtol
    )
    tolerance = _matrix_predicate_threshold(one(real_type), absolute, relative)
    return _matrix_predicate_result(
        :completely_positive,
        MatrixPredicateSatisfied;
        value=zero(real_type),
        tolerance=tolerance,
        witness=nothing,
        message="the supplied Kraus representation is completely positive by construction",
    )
end

"""
    is_completely_positive(map; atol=nothing, rtol=nothing,
                           allow_densify=false) -> MatrixPredicateResult
    is_trace_preserving(map; atol, rtol)
    is_unital(map; atol, rtol)

Diagnose complete positivity without modifying the representation. A Kraus
representation is satisfied by construction. Other representations first
require exact Choi Hermiticity outside any numerical boundary, then use the
structured positive-semidefiniteness diagnostic:

- `MatrixPredicateSatisfied` establishes complete positivity;
- `MatrixPredicateViolated` records a Choi-Hermiticity or negativity witness;
- `MatrixPredicateUnknown` preserves a nonzero tolerance-boundary defect.

The Choi matrix is never symmetrized or projected. Spectral work on sparse
Choi data requires `allow_densify=true`. Exact integer and rational data use
exact comparisons and require zero tolerances. Floating data default to the
matrix-predicate relative tolerance.

Trace preservation checks `Tr_output(J) == I_input`; unitality checks
`Φ(I_input) == I_output` and retain their Boolean APIs.
"""
function is_completely_positive(
    map::KrausRepresentation; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    _validated_kraus_operators(map)
    return _kraus_complete_positivity_result(map, atol, rtol)
end

function is_completely_positive(
    map::Union{OperatorSumRepresentation,ChoiRepresentation,SuperoperatorRepresentation};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    _require_square_operator_algebras(map, "is_completely_positive")
    choi = map isa ChoiRepresentation ? map : choi_representation(map)
    return _complete_positivity_result_from_choi(
        _validated_representation_matrix(choi);
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
    )
end

function _trace_output(matrix, input_dim::Int, output_dim::Int)
    result = zeros(eltype(matrix), input_dim, input_dim)
    @inbounds for input_column in 1:input_dim,
        input_row in 1:input_dim,
        output_index in 1:output_dim

        result[input_row, input_column] += matrix[
            output_index + (input_row - 1) * output_dim,
            output_index + (input_column - 1) * output_dim,
        ]
    end
    return result
end

function is_trace_preserving(
    map::AbstractMapRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    input_dim, output_dim = _require_square_operator_algebras(map, "is_trace_preserving")
    choi = choi_representation(map)
    reduced = _trace_output(choi.matrix, input_dim, output_dim)
    identity = Matrix{eltype(reduced)}(I, input_dim, input_dim)
    return isapprox(reduced, identity; atol=atol, rtol=rtol)
end

function is_unital(
    map::AbstractMapRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    input_dim, output_dim = _require_square_operator_algebras(map, "is_unital")
    input_identity = Matrix{eltype(map)}(I, input_dim, input_dim)
    output_identity = Matrix{eltype(map)}(I, output_dim, output_dim)
    output = apply_channel(input_identity, map)
    return isapprox(output, output_identity; atol=atol, rtol=rtol)
end
