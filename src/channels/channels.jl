# Source-informed independent Julia implementation based on the specifications
# of QETLAB ApplyMap.m, ChoiMatrix.m, and KrausOperators.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export AbstractMapRepresentation,
    KrausRepresentation,
    ChoiRepresentation,
    SuperoperatorRepresentation,
    input_dimension,
    output_dimension,
    kraus_representation,
    choi_representation,
    superoperator_representation,
    kraus_operators,
    choi_matrix,
    superoperator_matrix,
    apply_channel,
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
    AbstractMapRepresentation{T}

Supertype for finite-dimensional linear-map representations with scalar
element type `T`. A map sends `input_dimension(map)`-square matrices to
`output_dimension(map)`-square matrices.

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
sparse matrices remain sparse. Trace preservation is deliberately not imposed
by the constructor.
"""
struct KrausRepresentation{T,V<:AbstractVector} <: AbstractMapRepresentation{T}
    operators::V
    input_dim::Int
    output_dim::Int

    function KrausRepresentation{T,V}(
        ::Val{:validated}, operators::V, input_dim::Int, output_dim::Int
    ) where {T,V<:AbstractVector}
        return new{T,V}(operators, input_dim, output_dim)
    end
end

"""
    ChoiRepresentation(matrix, input_dim, output_dim)
    ChoiRepresentation(matrix)

Wrap a finite numeric Choi matrix using the input-first tensor ordering
``J[(i,a),(j,b)] = <a|Φ(|i><j|)|b>``. The explicit dimensions are required
for unequal input and output dimensions. The one-argument form infers an
equal input/output dimension from a square matrix of size `d^2`.

The matrix is copied, and sparse storage is preserved. Construction validates
shape and finite entries but does not require positivity, trace preservation,
or Hermiticity, so this type also represents nonphysical linear maps.
"""
struct ChoiRepresentation{T,M<:AbstractMatrix{T}} <: AbstractMapRepresentation{T}
    matrix::M
    input_dim::Int
    output_dim::Int

    function ChoiRepresentation{T,M}(
        ::Val{:validated}, matrix::M, input_dim::Int, output_dim::Int
    ) where {T,M<:AbstractMatrix{T}}
        return new{T,M}(matrix, input_dim, output_dim)
    end
end

"""
    SuperoperatorRepresentation(matrix, input_dim, output_dim)
    SuperoperatorRepresentation(matrix)

Wrap a transfer matrix `S` satisfying `vec(Φ(X)) = S * vec(X)`, with Julia's
column-major `vec`. Its size is `(output_dim^2, input_dim^2)`. The
one-argument form infers both dimensions when the row and column counts are
perfect squares.

The matrix is copied, and sparse storage is preserved. No physicality
condition is imposed.
"""
struct SuperoperatorRepresentation{T,M<:AbstractMatrix{T}} <: AbstractMapRepresentation{T}
    matrix::M
    input_dim::Int
    output_dim::Int

    function SuperoperatorRepresentation{T,M}(
        ::Val{:validated}, matrix::M, input_dim::Int, output_dim::Int
    ) where {T,M<:AbstractMatrix{T}}
        return new{T,M}(matrix, input_dim, output_dim)
    end
end

Base.eltype(::Type{<:AbstractMapRepresentation{T}}) where {T} = T
Base.eltype(map::AbstractMapRepresentation) = eltype(typeof(map))

input_dimension(map::AbstractMapRepresentation) = map.input_dim
output_dimension(map::AbstractMapRepresentation) = map.output_dim

_all_finite(matrix::AbstractMatrix) = all(isfinite, matrix)
_all_finite(matrix::SparseMatrixCSC) = all(isfinite, nonzeros(matrix))

function _validate_numeric_matrix(matrix::AbstractMatrix, name::AbstractString)
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
    return KrausRepresentation{scalar_type,typeof(copied)}(
        Val(:validated), copied, input_dim, output_dim
    )
end

function ChoiRepresentation(matrix::AbstractMatrix, input_dim, output_dim)
    _validate_numeric_matrix(matrix, "Choi matrix")
    checked_input = _map_dimension(input_dim, "input_dim")
    checked_output = _map_dimension(output_dim, "output_dim")
    total = try
        Base.checked_mul(checked_input, checked_output)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("input_dim * output_dim exceeds typemax(Int)"))
    end
    size(matrix) == (total, total) || throw(
        DimensionMismatch(
            "Choi matrix has size $(size(matrix)); expected ($total, $total) " *
            "for input_dim=$checked_input and output_dim=$checked_output",
        ),
    )
    copied = copy(matrix)
    return ChoiRepresentation{eltype(copied),typeof(copied)}(
        Val(:validated), copied, checked_input, checked_output
    )
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

function SuperoperatorRepresentation(matrix::AbstractMatrix, input_dim, output_dim)
    _validate_numeric_matrix(matrix, "superoperator matrix")
    checked_input = _map_dimension(input_dim, "input_dim")
    checked_output = _map_dimension(output_dim, "output_dim")
    expected = (
        _squared_dimension(checked_output, "output_dim"),
        _squared_dimension(checked_input, "input_dim"),
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "superoperator matrix has size $(size(matrix)); expected $expected " *
            "for input_dim=$checked_input and output_dim=$checked_output",
        ),
    )
    copied = copy(matrix)
    return SuperoperatorRepresentation{eltype(copied),typeof(copied)}(
        Val(:validated), copied, checked_input, checked_output
    )
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

function _choi_to_superoperator(matrix::SparseMatrixCSC, input_dim::Int, output_dim::Int)
    rows, columns, values = findnz(matrix)
    super_rows = Vector{Int}(undef, length(rows))
    super_columns = Vector{Int}(undef, length(rows))
    @inbounds for index in eachindex(rows)
        output_row = mod(rows[index] - 1, output_dim) + 1
        input_row = div(rows[index] - 1, output_dim) + 1
        output_column = mod(columns[index] - 1, output_dim) + 1
        input_column = div(columns[index] - 1, output_dim) + 1
        super_rows[index] = output_row + (output_column - 1) * output_dim
        super_columns[index] = input_row + (input_column - 1) * input_dim
    end
    return sparse(super_rows, super_columns, copy(values), output_dim^2, input_dim^2)
end

function _choi_to_superoperator(matrix::AbstractMatrix, input_dim::Int, output_dim::Int)
    result = Matrix{eltype(matrix)}(undef, output_dim^2, input_dim^2)
    @inbounds for input_column in 1:input_dim,
        input_row in 1:input_dim, output_column in 1:output_dim,
        output_row in 1:output_dim

        result[output_row + (output_column - 1) * output_dim, input_row + (input_column - 1) * input_dim] = matrix[
            output_row + (input_row - 1) * output_dim,
            output_column + (input_column - 1) * output_dim,
        ]
    end
    return result
end

function _superoperator_to_choi(matrix::SparseMatrixCSC, input_dim::Int, output_dim::Int)
    rows, columns, values = findnz(matrix)
    choi_rows = Vector{Int}(undef, length(rows))
    choi_columns = Vector{Int}(undef, length(rows))
    @inbounds for index in eachindex(rows)
        output_row = mod(rows[index] - 1, output_dim) + 1
        output_column = div(rows[index] - 1, output_dim) + 1
        input_row = mod(columns[index] - 1, input_dim) + 1
        input_column = div(columns[index] - 1, input_dim) + 1
        choi_rows[index] = output_row + (input_row - 1) * output_dim
        choi_columns[index] = output_column + (input_column - 1) * output_dim
    end
    total = input_dim * output_dim
    return sparse(choi_rows, choi_columns, copy(values), total, total)
end

function _superoperator_to_choi(matrix::AbstractMatrix, input_dim::Int, output_dim::Int)
    total = input_dim * output_dim
    result = Matrix{eltype(matrix)}(undef, total, total)
    @inbounds for input_column in 1:input_dim,
        input_row in 1:input_dim, output_column in 1:output_dim,
        output_row in 1:output_dim

        result[output_row + (input_row - 1) * output_dim, output_column + (input_column - 1) * output_dim] = matrix[
            output_row + (output_column - 1) * output_dim,
            input_row + (input_column - 1) * input_dim,
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

"""
    choi_representation(map)
    superoperator_representation(map)
    kraus_representation(map; atol, rtol)

Convert between map representations. Choi/superoperator conversion is an
exact index reshuffle and preserves sparse storage. Kraus-to-matrix conversion
preserves sparse storage when every Kraus operator is sparse.

Choi-to-Kraus conversion is defined only for completely positive maps. It
performs a dense Hermitian eigendecomposition and therefore explicitly
densifies sparse Choi matrices. Eigenvalues whose magnitude is at most
`atol + rtol * max(norm(J, Inf), 1)` are treated as numerical zero; a negative
eigenvalue below that threshold raises `DomainError`.
"""
function choi_representation(map::ChoiRepresentation)
    return ChoiRepresentation(map.matrix, map.input_dim, map.output_dim)
end

function choi_representation(map::KrausRepresentation)
    return ChoiRepresentation(_kraus_to_choi(map.operators), map.input_dim, map.output_dim)
end

function choi_representation(map::SuperoperatorRepresentation)
    return ChoiRepresentation(
        _superoperator_to_choi(map.matrix, map.input_dim, map.output_dim),
        map.input_dim,
        map.output_dim,
    )
end

function superoperator_representation(map::SuperoperatorRepresentation)
    return SuperoperatorRepresentation(map.matrix, map.input_dim, map.output_dim)
end

function superoperator_representation(map::KrausRepresentation)
    return SuperoperatorRepresentation(
        _kraus_to_superoperator(map.operators), map.input_dim, map.output_dim
    )
end

function superoperator_representation(map::ChoiRepresentation)
    return SuperoperatorRepresentation(
        _choi_to_superoperator(map.matrix, map.input_dim, map.output_dim),
        map.input_dim,
        map.output_dim,
    )
end

function kraus_representation(
    map::KrausRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    _checked_tolerances(map, atol, rtol)
    return KrausRepresentation(map.operators)
end

function _real_float_type(::Type{T}) where {T}
    return typeof(float(real(zero(T))))
end

function _default_rtol(map::AbstractMapRepresentation)
    return sqrt(eps(_real_float_type(eltype(map))))
end

function _checked_tolerances(map, atol, rtol)
    atol isa Real && isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("atol must be a finite nonnegative real number"))
    rtol isa Real && isfinite(rtol) && rtol >= 0 ||
        throw(ArgumentError("rtol must be a finite nonnegative real number"))
    return atol, rtol
end

function _hermitian_part_for_diagnostic(matrix::AbstractMatrix)
    dense = Matrix(matrix)
    return Hermitian((dense + adjoint(dense)) / 2)
end

function kraus_representation(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation};
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    choi = map isa ChoiRepresentation ? map : choi_representation(map)
    matrix = choi.matrix
    isapprox(matrix, adjoint(matrix); atol=atol, rtol=rtol) ||
        throw(DomainError(matrix, "the map's Choi matrix is not Hermitian"))

    scale = max(norm(matrix, Inf), one(_real_float_type(eltype(matrix))))
    threshold = atol + rtol * scale
    decomposition = eigen(_hermitian_part_for_diagnostic(matrix))
    minimum(decomposition.values) >= -threshold || throw(
        DomainError(
            minimum(decomposition.values),
            "the map is not completely positive at the requested tolerances",
        ),
    )

    selected = findall(value -> value > threshold, decomposition.values)
    operators = Base.map(reverse(selected)) do index
        value = decomposition.values[index]
        return reshape(
            sqrt(value) * decomposition.vectors[:, index], choi.output_dim, choi.input_dim
        )
    end
    if isempty(operators)
        scalar_type = eltype(decomposition.vectors)
        push!(operators, zeros(scalar_type, choi.output_dim, choi.input_dim))
    end
    return KrausRepresentation(operators)
end

"""
    choi_matrix(map)
    superoperator_matrix(map)
    kraus_operators(map; atol, rtol)

Return copied ordinary-array data for a representation. Conversion is
performed when necessary. `kraus_operators` follows the explicit dense
eigendecomposition behavior of `kraus_representation`.
"""
choi_matrix(map::AbstractMapRepresentation) = copy(choi_representation(map).matrix)

function superoperator_matrix(map::AbstractMapRepresentation)
    return copy(superoperator_representation(map).matrix)
end

function kraus_operators(map::AbstractMapRepresentation; kwargs...)
    return [copy(operator) for operator in kraus_representation(map; kwargs...).operators]
end

function _validate_channel_input(input::AbstractMatrix, map::AbstractMapRepresentation)
    _validate_numeric_matrix(input, "channel input")
    expected = input_dimension(map)
    size(input) == (expected, expected) || throw(
        DimensionMismatch(
            "channel input has size $(size(input)); expected ($expected, $expected)"
        ),
    )
    return nothing
end

"""
    apply_channel(input, map)
    apply_channel(input, kraus_operators)

Apply a finite-dimensional map to a square matrix. `map` may be a
`KrausRepresentation`, `ChoiRepresentation`, or
`SuperoperatorRepresentation`. A vector/tuple of matrices is interpreted as
Kraus operators.

Kraus application evaluates `sum(K * input * K')`. Choi application uses the
exact Choi-to-superoperator reshuffle. Sparse output is retained when sparse
linear algebra permits it; no input is normalized, symmetrized, or repaired.
"""
function apply_channel(input::AbstractMatrix, map::KrausRepresentation)
    _validate_channel_input(input, map)
    terms = Base.map(map.operators) do operator
        return operator * input * adjoint(operator)
    end
    return reduce(+, terms)
end

function apply_channel(input::AbstractMatrix, map::SuperoperatorRepresentation)
    _validate_channel_input(input, map)
    output_vector = map.matrix * vec(input)
    return reshape(output_vector, map.output_dim, map.output_dim)
end

function apply_channel(input::AbstractMatrix, map::ChoiRepresentation)
    _validate_channel_input(input, map)
    return apply_channel(input, superoperator_representation(map))
end

function apply_channel(input::AbstractMatrix, operators::Union{Tuple,AbstractVector})
    return apply_channel(input, KrausRepresentation(operators))
end

"""
    is_completely_positive(map; atol, rtol)
    is_trace_preserving(map; atol, rtol)
    is_unital(map; atol, rtol)

Numerically test standard map properties without modifying the representation.
The default relative tolerance is `sqrt(eps(T))` for the real floating type
associated with `eltype(map)`; the default absolute tolerance is zero.

Complete positivity first requires the Choi matrix to be Hermitian within the
requested tolerances, then tests the smallest eigenvalue of its Hermitian part.
The Hermitian part is used only for the diagnostic eigenvalue calculation and
is never returned as a repaired representation. This eigenvalue calculation
explicitly densifies sparse Choi data. Trace preservation checks
`Tr_output(J) == I_input`; unitality checks `Φ(I_input) == I_output`.
"""
function is_completely_positive(
    map::KrausRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    _checked_tolerances(map, atol, rtol)
    return true
end

function is_completely_positive(
    map::Union{ChoiRepresentation,SuperoperatorRepresentation};
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    choi = map isa ChoiRepresentation ? map : choi_representation(map)
    matrix = choi.matrix
    isapprox(matrix, adjoint(matrix); atol=atol, rtol=rtol) || return false
    scale = max(norm(matrix, Inf), one(_real_float_type(eltype(matrix))))
    threshold = atol + rtol * scale
    values = eigvals(_hermitian_part_for_diagnostic(matrix))
    return minimum(values) >= -threshold
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
    choi = choi_representation(map)
    reduced = _trace_output(choi.matrix, choi.input_dim, choi.output_dim)
    identity = Matrix{eltype(reduced)}(I, choi.input_dim, choi.input_dim)
    return isapprox(reduced, identity; atol=atol, rtol=rtol)
end

function is_unital(
    map::AbstractMapRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    atol, rtol = _checked_tolerances(map, atol, rtol)
    input_identity = Matrix{eltype(map)}(I, map.input_dim, map.input_dim)
    output_identity = Matrix{eltype(map)}(I, map.output_dim, map.output_dim)
    output = apply_channel(input_identity, map)
    return isapprox(output, output_identity; atol=atol, rtol=rtol)
end
