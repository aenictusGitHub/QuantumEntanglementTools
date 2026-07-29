# Source-informed independent Julia implementations based on the specifications
# and QETLAB RandomProbabilities.m, RandomStateVector.m,
# RandomDensityMatrix.m, RandomUnitary.m, RandomGraph.m, and RandomPOVM.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

using Random

export random_probabilities,
    random_state_vector, random_density_matrix, random_unitary, random_graph, random_povm

function _ginibre(rng::AbstractRNG, rows::Int, columns::Int, real_output::Bool)
    real_part = randn(rng, rows, columns)
    return real_output ? real_part : complex.(real_part, randn(rng, rows, columns))
end

"""
    random_probabilities(rng, n)

Draw a probability vector uniformly from the `(n-1)`-simplex (the
Dirichlet distribution with every concentration equal to one).

An explicit `rng::AbstractRNG` is mandatory.  The global random stream is
never read or mutated.
"""
function random_probabilities(rng::AbstractRNG, n)
    sample_count = _positive_int(n, "n")
    samples = randexp(rng, sample_count)
    total = sum(samples)
    isfinite(total) && !iszero(total) ||
        throw(ArgumentError("the RNG produced a non-finite probability draw"))
    return samples / total
end

function _random_dimension_tuple(dimension)
    if dimension isa Integer
        checked = _positive_int(dimension, "dimension")
        return (checked,)
    elseif dimension isa Tuple || dimension isa AbstractVector
        dimension isa AbstractVector && Base.require_one_based_indexing(dimension)
        length(dimension) == 2 || throw(
            DimensionMismatch("a bipartite dimension must contain exactly two entries")
        )
        return (
            _positive_int(dimension[1], "dimension[1]"),
            _positive_int(dimension[2], "dimension[2]"),
        )
    end
    return throw(
        ArgumentError("dimension must be a positive integer or a two-entry tuple/vector")
    )
end

"""
    random_state_vector(rng, dimension; real=false, schmidt_rank=nothing)

Draw a normalized Gaussian (Haar) state vector using only `rng`.

If `schmidt_rank` is supplied, `dimension` describes one local dimension
(used for both parties) or a two-entry local-dimension tuple.  The returned
bipartite vector has Schmidt rank at most that value and, almost surely, that
rank exactly.  The first subsystem is the slowest-varying basis index.
"""
function random_state_vector(
    rng::AbstractRNG, dimension; real::Bool=false, schmidt_rank=nothing
)
    dims = _random_dimension_tuple(dimension)
    if schmidt_rank === nothing
        total = _checked_product(dims, "dimension")
        vector = vec(_ginibre(rng, total, 1, real))
        vector_norm = norm(vector)
        isfinite(vector_norm) && !iszero(vector_norm) ||
            throw(ArgumentError("the RNG produced a zero or non-finite vector"))
        return vector / vector_norm
    end

    local_dims = length(dims) == 1 ? (dims[1], dims[1]) : dims
    rank_bound = _positive_int(schmidt_rank, "schmidt_rank")
    rank_bound <= min(local_dims...) || throw(
        ArgumentError(
            "schmidt_rank=$rank_bound exceeds the smaller local dimension $(min(local_dims...))",
        ),
    )
    left = _ginibre(rng, local_dims[1], rank_bound, real)
    right = _ginibre(rng, local_dims[2], rank_bound, real)
    coefficient_matrix = left * transpose(right)
    vector = collect(vec(transpose(coefficient_matrix)))
    vector_norm = norm(vector)
    isfinite(vector_norm) && !iszero(vector_norm) ||
        throw(ArgumentError("the RNG produced a zero or non-finite vector"))
    return vector / vector_norm
end

"""
    random_unitary(rng, dim; real=false)

Draw a Haar unitary, or a Haar orthogonal matrix when `real=true`, by
phase-correcting the QR factorization of a Ginibre matrix.

An explicit RNG is mandatory.
"""
function random_unitary(rng::AbstractRNG, dim; real::Bool=false)
    dimension = _positive_int(dim, "dim")
    ginibre = _ginibre(rng, dimension, dimension, real)
    factorization = qr(ginibre)
    q = Matrix(factorization.Q)
    diagonal = diag(factorization.R)
    phases = map(diagonal) do value
        return iszero(value) ? one(value) : value / abs(value)
    end
    return q * Diagonal(phases)
end

function _density_distribution(distribution)
    normalized = if distribution isa Symbol
        distribution
    elseif distribution isa AbstractString
        Symbol(lowercase(distribution))
    else
        throw(
            ArgumentError(
                "distribution must be a symbol or string; got $(repr(distribution))"
            ),
        )
    end
    normalized in (:hilbert_schmidt, :hs, :haar) && return :hilbert_schmidt
    normalized == :bures && return :bures
    return throw(
        ArgumentError(
            "distribution must be :hilbert_schmidt, :hs, :haar, or :bures; got $(repr(distribution))",
        ),
    )
end

"""
    random_density_matrix(rng, dim;
                          real=false, rank=dim,
                          distribution=:hilbert_schmidt)

Draw a positive semidefinite, trace-one matrix from an induced Ginibre
ensemble.  `rank` is the number of Ginibre columns and bounds the matrix
rank.  `distribution=:bures` applies the standard `(I+U)G` construction.

Only the explicit `rng` is used, including for the Bures unitary.
"""
function random_density_matrix(
    rng::AbstractRNG, dim; real::Bool=false, rank=dim, distribution=:hilbert_schmidt
)
    dimension = _positive_int(dim, "dim")
    rank_bound = _positive_int(rank, "rank")
    rank_bound <= dimension ||
        throw(ArgumentError("rank=$rank_bound must not exceed dim=$dimension"))
    selected_distribution = _density_distribution(distribution)
    ginibre = _ginibre(rng, dimension, rank_bound, real)
    if selected_distribution == :bures
        unitary = random_unitary(rng, dimension; real=real)
        ginibre = (Matrix{eltype(unitary)}(I, dimension, dimension) + unitary) * ginibre
    end
    rho = ginibre * ginibre'
    normalization = Base.real(tr(rho))
    isfinite(normalization) && normalization > zero(normalization) ||
        throw(ArgumentError("the RNG produced a zero or non-finite density-matrix trace"))
    return rho / normalization
end

"""
    random_graph(rng, n; edge_probability=0.5)

Draw the loop-free undirected Erdős--Rényi graph `G(n,p)` and return its
symmetric `BitMatrix` adjacency matrix.  Every unordered edge uses one draw
from the explicit RNG.
"""
function random_graph(rng::AbstractRNG, n; edge_probability=0.5)
    vertex_count = _positive_int(n, "n")
    probability = _finite_real_parameter(edge_probability, "edge_probability")
    zero(probability) <= probability <= one(probability) ||
        throw(ArgumentError("edge_probability must lie in [0, 1]; got $probability"))
    adjacency = falses(vertex_count, vertex_count)
    for column in 2:vertex_count
        for row in 1:(column - 1)
            present = rand(rng) < probability
            adjacency[row, column] = present
            adjacency[column, row] = present
        end
    end
    return adjacency
end

"""
    random_povm(rng, dim, outcomes; real=false)

Draw a POVM from a Haar isometry.  If `Q` is the first `dim` columns of an
`outcomes*dim` dimensional Haar unitary and `Q_j` are its row blocks, the
returned effects are `Q_j'Q_j`.  They are positive semidefinite and sum to
identity up to factorization roundoff.

The return value is a vector of dense matrices and the explicit RNG is the
only source of randomness.
"""
function random_povm(rng::AbstractRNG, dim, outcomes; real::Bool=false)
    dimension = _positive_int(dim, "dim")
    outcome_count = _positive_int(outcomes, "outcomes")
    row_count = _checked_product((dimension, outcome_count), "dim*outcomes")
    ginibre = _ginibre(rng, row_count, dimension, real)
    factorization = qr(ginibre)
    selector = Matrix{eltype(ginibre)}(I, row_count, dimension)
    isometry = factorization.Q * selector

    effects = Vector{Matrix{eltype(ginibre)}}(undef, outcome_count)
    for outcome in 1:outcome_count
        rows = ((outcome - 1) * dimension + 1):(outcome * dimension)
        block = isometry[rows, :]
        effects[outcome] = block' * block
    end
    return effects
end
