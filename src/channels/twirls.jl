# Source-informed independent Julia implementation based on the specification
# and QETLAB Twirl.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

export twirl

const _TWIRL_KINDS = (:werner, :isotropic, :real, :pauli)
const _TWIRL_DEFAULT_MAX_BASIS_SIZE = 256
const _TWIRL_DEFAULT_MAX_NONZEROS = 5_000_000
const _TWIRL_DEFAULT_MAX_DENSE_ENTRIES = 1_000_000
const _TWIRL_DEFAULT_MAX_WORK = 1_000_000_000

function _twirl_kind(kind)
    kind isa Symbol || throw(
        ArgumentError(
            "kind must be one of :werner, :isotropic, :real, or :pauli; " *
            "got $(repr(kind))",
        ),
    )
    kind in _TWIRL_KINDS || throw(
        ArgumentError(
            "kind must be one of :werner, :isotropic, :real, or :pauli; " *
            "got $(repr(kind))",
        ),
    )
    return kind
end

function _twirl_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _twirl_check_resource(
    planned::BigInt, limit, name::AbstractString, resource::AbstractString
)
    limit === nothing && return nothing
    planned <= limit || throw(
        ArgumentError(
            "twirl requires $planned $resource, exceeding $name=$limit; " *
            "raise the explicit guard only after reviewing the resource cost",
        ),
    )
    return nothing
end

function _twirl_exact_local_dimension(total::Int, copies::Int)
    total == 1 && return 1
    copies <= 8 * sizeof(Int) || throw(
        DimensionMismatch(
            "matrix dimension $total is not an exact $copies-th power of an " *
            "integer local dimension",
        ),
    )

    target = BigInt(total)
    lower = 2
    upper = total
    while lower <= upper
        midpoint = lower + div(upper - lower, 2)
        power = BigInt(midpoint)^copies
        if power == target
            return midpoint
        elseif power < target
            lower = midpoint + 1
        else
            upper = midpoint - 1
        end
    end
    return throw(
        DimensionMismatch(
            "matrix dimension $total is not an exact $copies-th power of an " *
            "integer local dimension",
        ),
    )
end

function _twirl_coefficient_type(::Type{T}) where {T<:Number}
    T === Bool && throw(ArgumentError("operator entries must not have type Bool"))
    if T <: Integer || T <: Rational
        return Rational{BigInt}
    elseif T <: Complex
        real_type = typeof(real(zero(T)))
        if real_type <: Integer || real_type <: Rational
            return Complex{Rational{BigInt}}
        end
    end
    coefficient_type = try
        typeof(one(T) / one(T))
    catch err
        (err isa MethodError || err isa ArgumentError) || rethrow()
        throw(
            ArgumentError("twirl cannot construct scalar coefficients for element type $T"),
        )
    end
    isconcretetype(coefficient_type) && coefficient_type <: Number || throw(
        ArgumentError(
            "twirl requires a concrete numeric coefficient type; inferred " *
            "$coefficient_type from input element type $T",
        ),
    )
    return coefficient_type
end

function _twirl_convert_operator(operator::AbstractMatrix, coefficient_type::Type)
    converted = try
        coefficient_type.(operator)
    catch err
        (err isa InexactError || err isa MethodError || err isa ArgumentError) ||
            rethrow()
        throw(
            ArgumentError(
                "operator entries cannot be represented by twirl coefficient type " *
                "$coefficient_type",
            ),
        )
    end
    return converted
end

function _twirl_validate(
    operator::AbstractMatrix,
    kind,
    copies;
    sparse_output::Bool,
    allow_densify::Bool,
    max_basis_size,
    max_nonzeros,
    max_dense_entries,
    max_work,
)
    _validate_numeric_matrix(operator, "operator")
    rows, columns = size(operator)
    rows == columns || throw(
        DimensionMismatch("twirl requires a square operator; got size $(size(operator))"),
    )
    rows > 0 || throw(ArgumentError("twirl requires a positive matrix dimension"))

    checked_kind = _twirl_kind(kind)
    copy_count = _positive_int(copies, "copies")
    if checked_kind in (:isotropic, :pauli)
        copy_count == 2 || throw(
            ArgumentError("kind=$checked_kind requires copies=2; got copies=$copy_count"),
        )
    end
    local_dimension = _twirl_exact_local_dimension(rows, copy_count)
    if checked_kind === :pauli
        ispow2(local_dimension) || throw(
            DimensionMismatch(
                "Pauli twirling requires a power-of-two local dimension; " *
                "got $local_dimension",
            ),
        )
    end

    basis_limit = _twirl_limit(max_basis_size, "max_basis_size")
    nonzero_limit = _twirl_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _twirl_limit(max_dense_entries, "max_dense_entries")
    work_limit = _twirl_limit(max_work, "max_work")
    dense_entries = BigInt(rows)^2
    if !sparse_output
        issparse(operator) &&
            !allow_densify &&
            throw(
                ArgumentError(
                    "sparse input requires allow_densify=true when sparse_output=false"
                ),
            )
        _twirl_check_resource(
            dense_entries, dense_limit, "max_dense_entries", "dense output entries"
        )
    end

    return (;
        kind=checked_kind,
        copies=copy_count,
        local_dimension,
        total_dimension=rows,
        coefficient_type=_twirl_coefficient_type(eltype(operator)),
        basis_limit,
        nonzero_limit,
        dense_limit,
        work_limit,
    )
end

function _twirl_next_permutation!(permutation::Vector{Int})
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

function _twirl_permutations(copies::Int)
    permutation = collect(1:copies)
    permutations = NTuple{copies,Int}[Tuple(permutation)]
    while _twirl_next_permutation!(permutation)
        push!(permutations, Tuple(permutation))
    end
    return permutations
end

function _twirl_permutation_count(copies::Int, basis_limit)
    count = BigInt(1)
    for factor in 2:copies
        count *= factor
        _twirl_check_resource(
            count, basis_limit, "max_basis_size", "Werner permutation operators"
        )
    end
    return count
end

function _twirl_matching_count(copies::Int, basis_limit)
    count = BigInt(1)
    for factor in 1:2:(2 * copies - 1)
        count *= factor
        _twirl_check_resource(
            count, basis_limit, "max_basis_size", "real-twirl Brauer operators"
        )
    end
    return count
end

function _twirl_relative_cycles(first::Tuple, second::Tuple)
    count = length(first)
    inverse_first = invperm(collect(first))
    relative = [inverse_first[second[index]] for index in 1:count]
    visited = falses(count)
    cycles = 0
    for start in 1:count
        visited[start] && continue
        cycles += 1
        current = start
        while !visited[current]
            visited[current] = true
            current = relative[current]
        end
    end
    return cycles
end

function _twirl_permutation_gram(permutations, local_dimension::Int)
    count = length(permutations)
    gram = Matrix{Int}(undef, count, count)
    for column in 1:count, row in 1:column
        cycles = _twirl_relative_cycles(permutations[row], permutations[column])
        value = _checked_power(local_dimension, cycles, "twirl permutation Gram entry")
        gram[row, column] = value
        gram[column, row] = value
    end
    return gram
end

function _twirl_matching_components(first, second, vertex_count::Int)
    parent = collect(1:vertex_count)
    function root(vertex::Int)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        return vertex
    end
    function join!(left::Int, right::Int)
        left_root = root(left)
        right_root = root(right)
        left_root == right_root || (parent[right_root] = left_root)
        return nothing
    end
    for (left, right) in first
        join!(left, right)
    end
    for (left, right) in second
        join!(left, right)
    end
    return length(Set(root(vertex) for vertex in 1:vertex_count))
end

function _twirl_matching_gram(matchings, local_dimension::Int, copies::Int)
    count = length(matchings)
    gram = Matrix{Int}(undef, count, count)
    for column in 1:count, row in 1:column
        components = _twirl_matching_components(
            matchings[row], matchings[column], 2 * copies
        )
        value = _checked_power(local_dimension, components, "twirl Brauer Gram entry")
        gram[row, column] = value
        gram[column, row] = value
    end
    return gram
end

function _twirl_independent_indices(gram::Matrix{Int})
    rows, columns = size(gram)
    rows == columns || error("internal twirl Gram matrix must be square")
    echelon = Rational{BigInt}.(gram)
    pivot_indices = Int[]
    pivot_row = 1
    for column in 1:columns
        candidate = findfirst(row -> !iszero(echelon[row, column]), pivot_row:rows)
        candidate === nothing && continue
        selected_row = pivot_row + candidate - 1
        if selected_row != pivot_row
            echelon[pivot_row, :], echelon[selected_row, :] = copy(
                echelon[selected_row, :]
            ),
            copy(echelon[pivot_row, :])
        end
        pivot = echelon[pivot_row, column]
        for row in (pivot_row + 1):rows
            iszero(echelon[row, column]) && continue
            factor = echelon[row, column] / pivot
            for trailing_column in column:columns
                echelon[row, trailing_column] -=
                    factor * echelon[pivot_row, trailing_column]
            end
        end
        push!(pivot_indices, column)
        pivot_row += 1
        pivot_row > rows && break
    end
    isempty(pivot_indices) && error("internal twirl spanning family has zero rank")
    return pivot_indices
end

function _twirl_sparse_inner(
    basis_operator::SparseMatrixCSC, operator::AbstractMatrix, coefficient_type::Type
)
    rows, columns, values = findnz(basis_operator)
    result = zero(coefficient_type)
    for index in eachindex(values)
        result +=
            conj(convert(coefficient_type, values[index])) *
            operator[rows[index], columns[index]]
    end
    return result
end

function _twirl_project_onto_basis(
    operator::AbstractMatrix,
    gram::Matrix{Int},
    independent_indices,
    basis_builder,
    coefficient_type::Type,
    dimension::Int,
)
    basis = [basis_builder(index, coefficient_type) for index in independent_indices]
    right_hand_side = coefficient_type[
        _twirl_sparse_inner(matrix, operator, coefficient_type) for matrix in basis
    ]
    reduced_gram = coefficient_type.(gram[independent_indices, independent_indices])
    coefficients = reduced_gram \ right_hand_side
    result = spzeros(coefficient_type, dimension, dimension)
    for (coefficient, matrix) in zip(coefficients, basis)
        result += coefficient * matrix
    end
    dropzeros!(result)
    return result
end

function _twirl_brauer_operator(
    matching, local_dimension::Int, copies::Int, coefficient_type::Type
)
    total_dimension = _checked_power(local_dimension, copies, "twirl local dimension")
    full_dims = ntuple(_ -> local_dimension, 2 * copies)
    digits = Vector{Int}(undef, 2 * copies)
    assignment_digits = Vector{Int}(undef, copies)
    rows = Vector{Int}(undef, total_dimension)
    columns = Vector{Int}(undef, total_dimension)
    values = fill(one(coefficient_type), total_dimension)
    for assignment in 0:(total_dimension - 1)
        residual = assignment
        for pair in copies:-1:1
            assignment_digits[pair] = mod(residual, local_dimension) + 1
            residual = div(residual, local_dimension)
        end
        for (pair, (first_vertex, second_vertex)) in pairs(matching)
            digits[first_vertex] = assignment_digits[pair]
            digits[second_vertex] = assignment_digits[pair]
        end
        linear = basis_to_linear(Tuple(digits), full_dims)
        rows[assignment + 1] = mod(linear - 1, total_dimension) + 1
        columns[assignment + 1] = div(linear - 1, total_dimension) + 1
    end
    return sparse(rows, columns, values, total_dimension, total_dimension)
end

function _twirl_finish(result::SparseMatrixCSC, sparse_output::Bool)
    return sparse_output ? result : Matrix(result)
end

function _twirl_check_span_budget(
    input_entries::BigInt,
    basis_count::BigInt,
    dimension::Int,
    cubic_gram::Bool,
    nonzero_limit,
    dense_limit,
    work_limit,
)
    stored_entries = basis_count * dimension
    _twirl_check_resource(
        stored_entries, nonzero_limit, "max_nonzeros", "spanning-family stored entries"
    )
    if cubic_gram
        compact_dense_entries = 3 * basis_count^2
        _twirl_check_resource(
            compact_dense_entries,
            dense_limit,
            "max_dense_entries",
            "compact Gram/elimination entries",
        )
    end
    work =
        input_entries +
        2 * stored_entries +
        basis_count^2 +
        (cubic_gram ? 2 * basis_count^3 : BigInt(0))
    _twirl_check_resource(work, work_limit, "max_work", "estimated scalar operations")
    return nothing
end

function _twirl_werner(operator, options)
    basis_count = _twirl_permutation_count(options.copies, options.basis_limit)
    _twirl_check_span_budget(
        BigInt(issparse(operator) ? nnz(operator) : length(operator)),
        basis_count,
        options.total_dimension,
        true,
        options.nonzero_limit,
        options.dense_limit,
        options.work_limit,
    )
    permutations = _twirl_permutations(options.copies)
    length(permutations) == basis_count ||
        error("internal Werner permutation count is inconsistent")
    gram = _twirl_permutation_gram(permutations, options.local_dimension)
    independent = _twirl_independent_indices(gram)
    dims = ntuple(_ -> options.local_dimension, options.copies)
    result = _twirl_project_onto_basis(
        operator,
        gram,
        independent,
        (index, coefficient_type) -> permutation_operator(
            dims, permutations[index]; T=coefficient_type, sparse_output=true
        ),
        options.coefficient_type,
        options.total_dimension,
    )
    return result
end

function _twirl_isotropic(operator, options)
    dimension = options.local_dimension
    total = options.total_dimension
    coefficient_type = options.coefficient_type
    input_entries = BigInt(issparse(operator) ? nnz(operator) : length(operator))
    _twirl_check_resource(
        BigInt(total == 1 ? 1 : 2),
        options.basis_limit,
        "max_basis_size",
        "isotropic projectors",
    )
    _twirl_check_resource(
        BigInt(2) * total,
        options.nonzero_limit,
        "max_nonzeros",
        "isotropic-projector stored entries",
    )
    _twirl_check_resource(
        input_entries + 4 * BigInt(total),
        options.work_limit,
        "max_work",
        "estimated scalar operations",
    )
    total == 1 && return sparse(operator)

    entangled_indices = [1 + (position - 1) * (dimension + 1) for position in 1:dimension]
    rows = repeat(entangled_indices; outer=dimension)
    columns = repeat(entangled_indices; inner=dimension)
    unnormalized_projector = sparse(
        rows, columns, fill(one(coefficient_type), total), total, total
    )
    overlap =
        _twirl_sparse_inner(unnormalized_projector, operator, coefficient_type) /
        coefficient_type(dimension)
    projector = unnormalized_projector / coefficient_type(dimension)
    complement_weight = (tr(operator) - overlap) / coefficient_type(total - 1)
    identity_matrix = spdiagm(0 => fill(one(coefficient_type), total))
    result = complement_weight * identity_matrix + (overlap - complement_weight) * projector
    dropzeros!(result)
    return result
end

function _twirl_real(operator, options)
    basis_count = _twirl_matching_count(options.copies, options.basis_limit)
    _twirl_check_span_budget(
        BigInt(issparse(operator) ? nnz(operator) : length(operator)),
        basis_count,
        options.total_dimension,
        true,
        options.nonzero_limit,
        options.dense_limit,
        options.work_limit,
    )
    matchings = _perfect_matchings(2 * options.copies)
    length(matchings) == basis_count ||
        error("internal real-twirl matching count is inconsistent")
    gram = _twirl_matching_gram(matchings, options.local_dimension, options.copies)
    independent = _twirl_independent_indices(gram)
    result = _twirl_project_onto_basis(
        operator,
        gram,
        independent,
        (index, coefficient_type) -> _twirl_brauer_operator(
            matchings[index], options.local_dimension, options.copies, coefficient_type
        ),
        options.coefficient_type,
        options.total_dimension,
    )
    return result
end

function _twirl_pauli_labels(index::Int, qubits::Int)
    labels = Vector{Int}(undef, qubits)
    residual = index
    for position in qubits:-1:1
        labels[position] = mod(residual, 4)
        residual = div(residual, 4)
    end
    return labels
end

function _twirl_pauli_vector(
    index::Int, qubits::Int, coefficient_type::Type, local_dimension::Int
)
    qubits == 0 && return sparsevec([1], [one(coefficient_type)], 1)
    labels = _twirl_pauli_labels(index, qubits)
    operator = pauli(labels; sparse_output=true)
    indices, values = findnz(sparsevec(vec(operator)))
    phase = (-im)^count(==(2), labels)
    real_values = map(values) do value
        phased = phase * value
        iszero(imag(phased)) || error("internal phase-adjusted Pauli vector is not real")
        return convert(coefficient_type, real(phased))
    end
    length(indices) == local_dimension ||
        error("internal Pauli vector support has unexpected size")
    return sparsevec(indices, real_values, local_dimension^2)
end

function _twirl_pauli_projector(
    index::Int, qubits::Int, coefficient_type::Type, local_dimension::Int
)
    vector = _twirl_pauli_vector(index, qubits, coefficient_type, local_dimension)
    return vector * adjoint(vector)
end

function _twirl_pauli(operator, options)
    basis_count = BigInt(options.total_dimension)
    _twirl_check_resource(
        basis_count, options.basis_limit, "max_basis_size", "Bell/Pauli projectors"
    )
    _twirl_check_span_budget(
        BigInt(issparse(operator) ? nnz(operator) : length(operator)),
        basis_count,
        options.total_dimension,
        false,
        options.nonzero_limit,
        options.dense_limit,
        options.work_limit,
    )
    qubits = 0
    residual = options.local_dimension
    while residual > 1
        residual = div(residual, 2)
        qubits += 1
    end

    coefficient_type = options.coefficient_type
    total = options.total_dimension
    result = spzeros(coefficient_type, total, total)
    for index in 0:(total - 1)
        projector = _twirl_pauli_projector(
            index, qubits, coefficient_type, options.local_dimension
        )
        coefficient =
            _twirl_sparse_inner(projector, operator, coefficient_type) /
            coefficient_type(total)
        result += coefficient * projector
    end
    dropzeros!(result)
    return result
end

@doc raw"""
    twirl(
        operator;
        kind=:werner,
        copies=2,
        sparse_output=issparse(operator),
        allow_densify=false,
        max_basis_size=256,
        max_nonzeros=5_000_000,
        max_dense_entries=1_000_000,
        max_work=1_000_000_000,
    )

Apply one of the four deterministic group twirls supported by the pinned
QETLAB `Twirl` entry point:

- `:werner` is the Haar average of
  `U^⊗copies * operator * adjoint(U^⊗copies)` and projects onto the span of
  the `copies!` subsystem-permutation operators.
- `:isotropic` requires `copies=2` and projects onto the maximally entangled
  projector and its orthogonal complement.
- `:real` is the analogous orthogonal-group average and projects onto the
  Brauer span indexed by the `(2copies-1)!!` perfect matchings.
- `:pauli` requires `copies=2` and a power-of-two local dimension; it projects
  onto the `4^q` Bell projectors obtained by vectorizing the `q`-qubit Pauli
  basis.

The input must be a finite, one-based, nonempty square numeric matrix of size
`d^copies` for an exact integer `d`. No entry, trace, Hermiticity, positivity,
or normalization is repaired. In particular, this linear map accepts general
operators, not only density matrices. `kind` is a strict lowercase `Symbol`,
and both isotropic and Pauli twirls require exactly two copies, correcting
the pinned routine's under-validation of `copies < 2`.

Sparse inputs stay sparse by default. The implementation forms only sparse
permutation/Brauer/projector operators and a compact dense integer Gram matrix;
it never densifies the input. A dense result from sparse input requires both
`sparse_output=false` and `allow_densify=true`. `max_basis_size` bounds the
factorial, double-factorial, or exponential spanning family before
enumeration. `max_nonzeros`, `max_dense_entries`, and `max_work` bound
conservative storage, output densification, and scalar-work estimates before
the corresponding allocations. Pass `nothing` to disable a guard explicitly.

Linearly dependent permutation or Brauer families are reduced using exact
integer-Gram row elimination before the coefficient solve. Integer and
rational inputs consequently produce exact `Rational{BigInt}` results.
Floating and complex floating types retain their precision when supported by
Julia's generic dense linear solve. No pseudoinverse tolerance or post-hoc
trace correction is used.

The Werner and real branches use `O(b^2)` Gram storage, `O(b^3)` compact
elimination/solve work, and at most `O(b*d^copies)` sparse construction work,
where `b` is respectively `copies!` or `(2copies-1)!!`. The Pauli branch uses
`b=4^q=d^2` projectors and `O(b*d^2)` work. The isotropic branch is linear in
the output dimension apart from scanning a dense input.

# Example

```jldoctest
julia> using LinearAlgebra, SparseArrays

julia> rho = sparse([0.4 0 0 0.1; 0 0.2 0 0; 0 0 0.1 0; 0.1 0 0 0.3]);

julia> twirled = twirl(rho; kind=:isotropic);

julia> tr(twirled) ≈ tr(rho)
true

julia> twirl(twirled; kind=:isotropic) ≈ twirled
true
```

This is an independent Julia implementation informed by QETLAB `Twirl.m` at
revision `d8589610f00cff106537268dee2e2a1153f3a601` (BSD-2-Clause).
""" function twirl(
    operator::AbstractMatrix;
    kind=:werner,
    copies=2,
    sparse_output::Bool=issparse(operator),
    allow_densify::Bool=false,
    max_basis_size=_TWIRL_DEFAULT_MAX_BASIS_SIZE,
    max_nonzeros=_TWIRL_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_TWIRL_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_TWIRL_DEFAULT_MAX_WORK,
)
    options = _twirl_validate(
        operator,
        kind,
        copies;
        sparse_output,
        allow_densify,
        max_basis_size,
        max_nonzeros,
        max_dense_entries,
        max_work,
    )
    converted = _twirl_convert_operator(operator, options.coefficient_type)
    result = if options.total_dimension == 1
        sparse(converted)
    elseif options.kind === :werner
        _twirl_werner(converted, options)
    elseif options.kind === :isotropic
        _twirl_isotropic(converted, options)
    elseif options.kind === :real
        _twirl_real(converted, options)
    else
        _twirl_pauli(converted, options)
    end
    return _twirl_finish(result, sparse_output)
end
