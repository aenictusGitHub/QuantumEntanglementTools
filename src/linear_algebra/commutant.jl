# Source-informed independent Julia implementation based on the specification
# and QETLAB Commutant.m and helpers/spnull.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# spnull: Copyright 2010 Bruno Luong, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt and
# licenses/QETLAB-BRUNO-LUONG-HELPERS-LICENSE.txt.

export commutant

const _COMMUTANT_DEFAULT_MAX_ENTRIES = 10_000_000
const _COMMUTANT_DEFAULT_MAX_WORK = 1_000_000_000

function _commutant_generators(matrix::AbstractMatrix)
    return AbstractMatrix[matrix]
end

function _commutant_generators(matrices::Union{Tuple,AbstractVector})
    matrices isa AbstractVector && Base.require_one_based_indexing(matrices)
    isempty(matrices) &&
        throw(ArgumentError("commutant requires at least one generator matrix"))
    generators = AbstractMatrix[]
    for (index, matrix) in pairs(matrices)
        matrix isa AbstractMatrix || throw(
            ArgumentError(
                "commutant generator $index must be a matrix; got $(typeof(matrix))"
            ),
        )
        push!(generators, matrix)
    end
    return generators
end

function _commutant_generators(input)
    return throw(
        ArgumentError(
            "commutant expects a square matrix or a nonempty tuple/vector of " *
            "equally sized square matrices; got $(typeof(input))",
        ),
    )
end

function _commutant_validate_generators(generators)
    first_generator = first(generators)
    Base.require_one_based_indexing(first_generator)
    rows, columns = size(first_generator)
    rows == columns || throw(
        DimensionMismatch(
            "commutant generators must be square; generator 1 has size " *
            "$(size(first_generator))",
        ),
    )
    rows > 0 || throw(ArgumentError("commutant generators must have positive dimension"))

    value_type = Union{}
    for (generator_index, generator) in pairs(generators)
        Base.require_one_based_indexing(generator)
        size(generator) == (rows, rows) || throw(
            DimensionMismatch(
                "commutant generators must have the same square size; generator " *
                "$generator_index has size $(size(generator)), expected ($rows, $rows)",
            ),
        )
        for (entry_index, value) in pairs(generator)
            value isa Number || throw(
                ArgumentError(
                    "commutant generator $generator_index entry $entry_index must " *
                    "be numeric; got $(repr(value))",
                ),
            )
            finite = try
                isfinite(value)
            catch err
                err isa MethodError || rethrow()
                throw(
                    ArgumentError(
                        "commutant cannot validate finiteness for generator " *
                        "$generator_index entry type $(typeof(value))",
                    ),
                )
            end
            finite || throw(
                ArgumentError(
                    "commutant generator $generator_index entry $entry_index must " *
                    "be finite; got $(repr(value))",
                ),
            )
            value_type = if value_type === Union{}
                typeof(value)
            else
                promote_type(value_type, typeof(value))
            end
        end
    end
    return rows, value_type
end

function _commutant_solver_type(value_type::Type)
    float_value = try
        float(zero(value_type))
    catch err
        (err isa MethodError || err isa ArgumentError) || rethrow()
        throw(
            ArgumentError(
                "commutant cannot select a floating SVD type for input value type " *
                "$value_type",
            ),
        )
    end
    solver_type = typeof(float_value)
    if solver_type <: Union{BigFloat,Complex{BigFloat}}
        throw(
            ArgumentError(
                "commutant does not down-convert BigFloat input: the dependency-free " *
                "core has no generic BigFloat SVD; convert explicitly to Float32 or " *
                "Float64, or use an external generic linear-algebra backend",
            ),
        )
    end
    solver_type <: Union{Float32,Float64,ComplexF32,ComplexF64} || throw(
        ArgumentError(
            "commutant requires an input type representable by the standard-library " *
            "Float32/Float64 dense SVD; got promoted floating type $solver_type",
        ),
    )
    return solver_type
end

function _commutant_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _commutant_check_budget(
    generator_count::Int, dimension::Int; max_entries, max_work
)
    columns = BigInt(dimension)^2
    rows = BigInt(generator_count) * columns
    entries = rows * columns
    work = rows * columns^2
    entry_limit = _commutant_limit(max_entries, "max_entries")
    work_limit = _commutant_limit(max_work, "max_work")
    if entry_limit !== nothing && entries > entry_limit
        throw(
            ArgumentError(
                "commutant would densify a $rows × $columns commutator " *
                "($entries entries), exceeding max_entries=$entry_limit",
            ),
        )
    end
    if work_limit !== nothing && work > work_limit
        throw(
            ArgumentError(
                "commutant dense SVD work estimate $work exceeds " *
                "max_work=$work_limit; reduce the problem or raise the explicit budget",
            ),
        )
    end
    return nothing
end

function _commutant_tolerance(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a nonnegative finite real number, not Bool"))
    value isa Real || throw(ArgumentError("$name must be a nonnegative finite real number"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    return value
end

function _commutant_sparse_operator(generators, solver_type::Type, dimension::Int)
    identity_matrix = sparse(
        1:dimension, 1:dimension, fill(one(solver_type), dimension), dimension, dimension
    )
    blocks = map(generators) do generator
        converted = sparse(solver_type.(generator))
        return kron(identity_matrix, converted) -
               kron(sparse(transpose(converted)), identity_matrix)
    end
    return reduce(vcat, blocks)
end

function _commutant_dense_operator(generators, solver_type::Type, dimension::Int)
    identity_matrix = Matrix{solver_type}(I, dimension, dimension)
    block_size = dimension^2
    operator = Matrix{solver_type}(undef, length(generators) * block_size, block_size)
    for (index, generator) in pairs(generators)
        converted = Matrix{solver_type}(generator)
        row_range = ((index - 1) * block_size + 1):(index * block_size)
        @views operator[row_range, :] .=
            kron(identity_matrix, converted) - kron(transpose(converted), identity_matrix)
    end
    return operator
end

"""
    commutant(generators;
              atol=0,
              rtol=nothing,
              allow_densify=false,
              max_entries=10_000_000,
              max_work=1_000_000_000)

Return a vector of matrices forming a Hilbert--Schmidt-orthonormal basis for
the matrices that commute with every supplied generator. `generators` may be
one square matrix or a nonempty tuple/vector of equally sized square matrices.
All matrices must use one-based axes and contain finite numeric values.

For an `n × n` matrix `X`, column-major vectorization is used. Each generator
`A` contributes the null-space equation

```math
\\left(I_n \\otimes A - A^\\mathsf{T} \\otimes I_n\\right)
\\mathrm{vec}(X) = 0.
```

Consequently, each returned matrix has size `n × n`, and vectorizing the
returned matrices gives orthonormal columns. A basis of a degenerate null
space is not unique; compare dimensions, residuals, or the orthogonal
projectors onto two returned spans rather than comparing basis entries.

The numerical rank threshold is
`max(atol, rtol * largest_singular_value)`. `atol` defaults to zero, and the
default `rtol` is `max(size(K)...) * eps(T)`, where `K` is the stacked
commutator and `T` is its real floating component type. `Float32`, `Float64`,
and their complex counterparts are supported. Integer and rational inputs are
converted to `Float64`. `BigFloat` is rejected explicitly because the
dependency-free core has no generic-precision SVD and never down-converts
precision silently.

Sparse generators are assembled into a sparse Kronecker commutator, but the
dependency-free null-space calculation is a dense SVD. Any sparse generator
therefore requires `allow_densify=true`. Before allocation, `max_entries`
limits the dense commutator size and `max_work` limits the conservative
`rows * columns^2` SVD work estimate; pass `nothing` to disable either guard.
For `g` generators, dense storage is `O(g n^4)` and the SVD estimate is
`O(g n^6)`.

# Examples

```jldoctest
julia> A = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0];

julia> basis = commutant(A);

julia> length(basis)
3

julia> maximum(maximum(abs, A * X - X * A) for X in basis) < 1e-12
true
```
"""
function commutant(
    input;
    atol=0,
    rtol=nothing,
    allow_densify::Bool=false,
    max_entries=_COMMUTANT_DEFAULT_MAX_ENTRIES,
    max_work=_COMMUTANT_DEFAULT_MAX_WORK,
)
    generators = _commutant_generators(input)
    dimension, value_type = _commutant_validate_generators(generators)
    solver_type = _commutant_solver_type(value_type)
    checked_atol = _commutant_tolerance(atol, "atol")
    generator_count = length(generators)

    contains_sparse = any(SparseArrays.issparse, generators)
    if contains_sparse && !allow_densify
        throw(
            ArgumentError(
                "commutant requires a dense null-space SVD; pass " *
                "allow_densify=true to permit bounded conversion of sparse generators",
            ),
        )
    end
    _commutant_check_budget(
        generator_count, dimension; max_entries=max_entries, max_work=max_work
    )
    column_count = try
        Base.checked_mul(dimension, dimension)
    catch err
        err isa OverflowError || rethrow()
        throw(
            ArgumentError(
                "commutant vectorized dimension $dimension^2 cannot be represented as Int",
            ),
        )
    end
    row_count = try
        Base.checked_mul(generator_count, column_count)
    catch err
        err isa OverflowError || rethrow()
        throw(
            ArgumentError(
                "commutant stacked row count cannot be represented as Int for " *
                "$generator_count generators of dimension $dimension",
            ),
        )
    end
    real_type = typeof(real(zero(solver_type)))
    checked_rtol = if rtol === nothing
        max(row_count, column_count) * eps(real_type)
    else
        _commutant_tolerance(rtol, "rtol")
    end

    all_sparse = all(SparseArrays.issparse, generators)
    structured_operator = if all_sparse
        _commutant_sparse_operator(generators, solver_type, dimension)
    else
        _commutant_dense_operator(generators, solver_type, dimension)
    end
    dense_operator = all_sparse ? Matrix(structured_operator) : structured_operator
    decomposition = svd(dense_operator; full=true)
    largest_singular_value = maximum(decomposition.S; init=zero(real_type))
    threshold = max(checked_atol, checked_rtol * largest_singular_value)
    numerical_rank = count(value -> value > threshold, decomposition.S)
    null_vectors = decomposition.V[:, (numerical_rank + 1):end]

    basis = Matrix{solver_type}[]
    sizehint!(basis, size(null_vectors, 2))
    for column in axes(null_vectors, 2)
        push!(basis, copy(reshape(@view(null_vectors[:, column]), dimension, dimension)))
    end
    return basis
end
