# Source-informed independent Julia implementations based on the specifications
# and QETLAB Pauli.m, GenPauli.m, GellMann.m, GenGellMann.m, and
# FourierMatrix.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export pauli, generalized_pauli, gell_mann, generalized_gell_mann, fourier_matrix

function _operator_index(value, upper::Int, name::AbstractString)
    index = _nonnegative_int(value, name)
    index <= upper ||
        throw(ArgumentError("$name=$index is outside the valid range 0:$upper"))
    return index
end

function _pauli_index(value)
    if value isa Integer
        return _operator_index(value, 3, "index")
    elseif value isa Char
        label = uppercase(value)
    elseif value isa Symbol
        label = uppercase(String(value))
    elseif value isa AbstractString
        ncodeunits(value) == 1 || throw(
            ArgumentError(
                "a Pauli string label must be one of I, X, Y, or Z; got $(repr(value))"
            ),
        )
        label = uppercase(value)
    else
        throw(ArgumentError("a Pauli label must be 0:3 or I, X, Y, Z; got $(repr(value))"))
    end

    if label == 'I' || label == "I"
        0
    elseif label == 'X' || label == "X"
        1
    elseif label == 'Y' || label == "Y"
        2
    elseif label == 'Z' || label == "Z"
        3
    else
        throw(ArgumentError("a Pauli label must be I, X, Y, or Z; got $(repr(value))"))
    end
end

function _single_pauli(index::Int)
    if index == 0
        return [1 0; 0 1]
    elseif index == 1
        return [0 1; 1 0]
    elseif index == 2
        return [0 -im; im 0]
    end
    return [1 0; 0 -1]
end

"""
    pauli(index; sparse_output=false)

Construct a one- or multi-qubit Pauli operator.

`index` may be an integer in `0:3`, one of `I`, `X`, `Y`, `Z` as a
character, string, or symbol, or a nonempty tuple/vector of such labels.  For
a collection, the first label acts on subsystem 1 (the slowest-varying tensor
factor).  Invalid labels are rejected rather than being treated as identity.
"""
function pauli(index; sparse_output::Bool=false)
    labels = if index isa Tuple || (index isa AbstractVector && !(index isa AbstractString))
        isempty(index) && throw(ArgumentError("a multi-qubit Pauli label cannot be empty"))
        Tuple(index)
    elseif index isa AbstractString && ncodeunits(index) > 1
        Tuple(index)
    else
        nothing
    end

    result = if labels === nothing
        factor = _single_pauli(_pauli_index(index))
        sparse_output ? sparse(factor) : factor
    else
        factors = map(labels) do label
            factor = _single_pauli(_pauli_index(label))
            return sparse_output ? sparse(factor) : factor
        end
        tensor_product(factors...)
    end
    return sparse_output ? result : Matrix(result)
end

"""
    generalized_pauli(shift, clock, dim; sparse_output=false, T=Float64)

Construct the Weyl operator ``X^{shift} Z^{clock}`` in dimension `dim`.

Both indices are zero-based and must lie in `0:dim-1`.  The shift satisfies
``X|j\\rangle=|j+shift \\pmod {dim}\\rangle`` and the clock phase is
``\\exp(2\\pi i\\,clock\\,j/dim)``.  The sparse representation has exactly
`dim` stored entries.
"""
function generalized_pauli(
    shift, clock, dim; sparse_output::Bool=false, T::Type{<:AbstractFloat}=Float64
)
    dimension = _positive_int(dim, "dim")
    shift_index = _operator_index(shift, dimension - 1, "shift")
    clock_index = _operator_index(clock, dimension - 1, "clock")

    columns = collect(1:dimension)
    rows = [mod(column - 1 + shift_index, dimension) + 1 for column in columns]
    angle_scale = T(2) * T(π) / T(dimension)
    phases = Complex{T}[
        cis(angle_scale * T(clock_index) * T(column - 1)) for column in columns
    ]
    result = sparse(rows, columns, phases, dimension, dimension)
    return sparse_output ? result : Matrix(result)
end

"""
    generalized_gell_mann(i, j, dim; sparse_output=false, T=Float64)

Construct a Hermitian generalized Gell-Mann basis matrix.

Indices are zero-based and belong to `0:dim-1`.  `(0,0)` gives identity,
`i < j` gives a symmetric off-diagonal matrix, `i > j` gives its imaginary
antisymmetric partner, and `(i,i)` for `i > 0` gives a traceless diagonal
matrix.  Nonidentity elements have Hilbert--Schmidt norm `sqrt(2)`.
"""
function generalized_gell_mann(
    first, second, dim; sparse_output::Bool=false, T::Type{<:AbstractFloat}=Float64
)
    dimension = _positive_int(dim, "dim")
    i = _operator_index(first, dimension - 1, "first index")
    j = _operator_index(second, dimension - 1, "second index")

    result = if i == j
        if i == 0
            spdiagm(0 => ones(T, dimension))
        else
            scale = sqrt(T(2) / (T(i) * T(i + 1)))
            diagonal = zeros(T, dimension)
            diagonal[1:i] .= scale
            diagonal[i + 1] = -i * scale
            spdiagm(0 => diagonal)
        end
    elseif i < j
        sparse([i + 1, j + 1], [j + 1, i + 1], T[1, 1], dimension, dimension)
    else
        sparse(
            [i + 1, j + 1],
            [j + 1, i + 1],
            Complex{T}[complex(zero(T), one(T)), complex(zero(T), -one(T))],
            dimension,
            dimension,
        )
    end

    return sparse_output ? result : Matrix(result)
end

const _GELL_MANN_INDICES = (
    (0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (2, 0), (1, 2), (2, 1), (2, 2)
)

"""
    gell_mann(index; sparse_output=false, T=Float64)

Construct the identity (`index == 0`) or one of the eight conventional
three-dimensional Gell-Mann matrices (`index in 1:8`).
"""
function gell_mann(index; sparse_output::Bool=false, T::Type{<:AbstractFloat}=Float64)
    checked_index = _operator_index(index, 8, "index")
    first, second = _GELL_MANN_INDICES[checked_index + 1]
    return generalized_gell_mann(first, second, 3; sparse_output=sparse_output, T=T)
end

"""
    fourier_matrix(dim; T=Float64)

Return the unitary `dim`-dimensional quantum Fourier matrix with entries
``F[j,k] = \\exp(2\\pi i(j-1)(k-1)/dim)/\\sqrt{dim}``.

`T` selects the real floating-point precision of the returned
`Matrix{Complex{T}}`.
"""
function fourier_matrix(dim; T::Type{<:AbstractFloat}=Float64)
    dimension = _positive_int(dim, "dim")
    scale = inv(sqrt(T(dimension)))
    angle_scale = T(2) * T(π) / T(dimension)
    return Complex{T}[
        scale * cis(angle_scale * T(row - 1) * T(column - 1)) for
        row in 1:dimension, column in 1:dimension
    ]
end
