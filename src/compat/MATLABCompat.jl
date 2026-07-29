# Source-informed independent Julia compatibility wrappers based on the
# corresponding QETLAB Tier-A and Tier-B entry points at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
MATLAB/QETLAB-style entry points for migration of existing scripts.

The wrappers in this namespace preserve QETLAB argument order and defaults
where the corresponding primitive is implemented.  Randomized entry points
add a mandatory leading `rng` argument to prevent global RNG mutation. New
Julia code should prefer the lowercase native API in
`QuantumEntanglementTools`.
"""
module MATLABCompat

using ..QuantumEntanglementTools:
    AbstractMapRepresentation,
    ChoiRepresentation,
    CriterionEntanglementDetected,
    CriterionResult,
    CriterionSatisfied,
    CriterionUnknown,
    KrausRepresentation,
    SubsystemPermutationPlan,
    SuperoperatorRepresentation,
    additive_compound_matrix,
    antisymmetric_projector,
    antisymmetric_subspace_basis,
    apply_channel,
    basis_to_linear,
    bell_state,
    brauer_states,
    breuer_state,
    chessboard_state,
    choi_map,
    choi_matrix,
    complementary_channel,
    compound_matrix,
    concurrence,
    dicke_state,
    dual_channel,
    elementary_symmetric_polynomial,
    entanglement_of_formation,
    fidelity,
    fourier_matrix,
    generalized_gell_mann,
    generalized_pauli,
    gell_mann,
    ghz_state,
    gisin_state,
    horodecki_state,
    in_separable_ball,
    is_locally_positive_semidefinite,
    is_positive_semidefinite,
    input_dimension,
    inverse_realign,
    isotropic_state,
    is_product_operator,
    is_product_vector,
    is_totally_nonsingular,
    is_totally_positive,
    ky_fan_norm,
    kraus_operators,
    kronecker_sum,
    l1_coherence,
    linear_to_basis,
    majorizes,
    maximally_entangled,
    negativity,
    operator_schmidt_decomposition,
    operator_schmidt_rank,
    output_dimension,
    pauli,
    partial_trace,
    partial_transpose,
    partial_map,
    pauli_channel,
    permutation_operator,
    permute_subsystems,
    random_density_matrix,
    random_graph,
    random_povm,
    random_probabilities,
    random_state_vector,
    random_unitary,
    realign,
    reduction_map,
    relative_entropy_coherence,
    schatten_norm,
    coherence_rank,
    schmidt_decomposition,
    schmidt_rank,
    swap_operator,
    symmetric_projector,
    symmetric_subspace_basis,
    tensor_power,
    tensor_product,
    tensor_sum,
    trace_norm,
    von_neumann_entropy,
    werner_state,
    w_state
using Random: AbstractRNG
using LinearAlgebra: Hermitian, diag, eigen, norm, svdvals, tr
using SparseArrays: AbstractSparseVector, findnz, issparse, sparse, spdiagm, sparsevec

export Tensor,
    TensorSum,
    KroneckerSum,
    PermuteSystems,
    PermutationOperator,
    Swap,
    SwapOperator,
    PartialTrace,
    PartialTranspose,
    Realignment,
    InverseRealignment,
    SymmetricProjection,
    AntisymmetricProjection,
    BasisToLinear,
    LinearToBasis,
    Pauli,
    GenPauli,
    GellMann,
    GenGellMann,
    FourierMatrix,
    MaxEntangled,
    Bell,
    GHZState,
    WState,
    DickeState,
    IsotropicState,
    WernerState,
    HorodeckiState,
    GisinState,
    BreuerState,
    BrauerStates,
    ChessboardState,
    RandomProbabilities,
    RandomStateVector,
    RandomDensityMatrix,
    RandomUnitary,
    RandomGraph,
    RandomPOVM,
    ApplyMap,
    ChoiMatrix,
    KrausOperators,
    ComplementaryMap,
    DualMap,
    PartialMap,
    DepolarizingChannel,
    DephasingChannel,
    PauliChannel,
    ChoiMap,
    ReductionMap,
    TraceNorm,
    SchattenNorm,
    KyFanNorm,
    Purity,
    Entropy,
    Fidelity,
    Negativity,
    SchmidtDecomposition,
    SchmidtRank,
    OperatorSchmidtDecomposition,
    OperatorSchmidtRank,
    IsProductVector,
    IsProductOperator,
    Concurrence,
    EntFormation,
    InSeparableBall,
    IsPPT,
    L1NormCoherence,
    RelEntCoherence,
    CoherenceRank,
    Majorizes,
    ElemSymPoly,
    CompoundMatrix,
    AdditiveCompoundMatrix,
    IsPSD,
    IsLocallyPSD,
    IsTotallyPositive,
    IsTotallyNonsingular

function _flag(value, name::AbstractString)
    if value isa Bool
        return value
    elseif value isa Integer && (value == 0 || value == 1)
        return value == 1
    end
    return throw(ArgumentError("$name must be Bool, 0, or 1; got $(repr(value))"))
end

function _positive_dimension(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a positive integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be a positive integer; got $(repr(value))"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return Int(value)
end

function _nonnegative_dimension(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a nonnegative integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be a nonnegative integer; got $(repr(value))"))
    value >= 0 || throw(ArgumentError("$name must be nonnegative; got $value"))
    return Int(value)
end

function _expand_scalar_dimension(value, total::Int, name::AbstractString)
    first_dimension = _positive_dimension(value, name)
    rem(total, first_dimension) == 0 ||
        throw(ArgumentError("$name=$first_dimension must evenly divide dimension $total"))
    return (first_dimension, div(total, first_dimension))
end

function _dimension_tuple(dim, name::AbstractString="DIM")
    if dim isa Tuple || dim isa AbstractVector
        return Tuple(dim)
    end
    return throw(
        ArgumentError("$name must be an integer, tuple/vector, or a 2-row dimension matrix")
    )
end

function _matrix_dimension_rows(dim::AbstractMatrix, subsystem_count=nothing)
    Base.require_one_based_indexing(dim)
    size(dim, 1) == 2 || throw(
        DimensionMismatch(
            "a row/column DIM matrix must have exactly 2 rows; got size $(size(dim))"
        ),
    )
    subsystem_count === nothing ||
        size(dim, 2) == subsystem_count ||
        throw(
            DimensionMismatch(
                "DIM has $(size(dim, 2)) subsystems; expected $subsystem_count"
            ),
        )
    return Tuple(dim[1, :]), Tuple(dim[2, :])
end

function _equal_default_dimensions(total::Int, subsystem_count::Int)
    subsystem_count > 0 ||
        throw(ArgumentError("a subsystem permutation must contain at least one system"))
    local_dimension = round(Int, total^(1 / subsystem_count))
    return ntuple(_ -> local_dimension, subsystem_count)
end

function _state_vector(matrix::AbstractMatrix)
    Base.require_one_based_indexing(matrix)
    min(size(matrix)...) == 1 ||
        throw(ArgumentError("matrix is not a row or column pure-state array"))
    if issparse(matrix)
        rows, columns, values = findnz(matrix)
        indices = size(matrix, 1) == 1 ? columns : rows
        return sparsevec(indices, copy(values), length(matrix))
    end
    return vec(matrix)
end

function _restore_state_orientation(vector::AbstractVector, original::AbstractMatrix)
    shaped = if size(original, 1) == 1
        reshape(vector, 1, length(vector))
    else
        reshape(vector, length(vector), 1)
    end
    return issparse(original) ? sparse(shaped) : shaped
end

"""QETLAB-compatible Kronecker tensor product."""
function Tensor(first, rest...)
    if (first isa Tuple || (first isa AbstractVector && !(eltype(first) <: Number))) &&
        isempty(rest)
        isempty(first) && throw(ArgumentError("Tensor requires at least one factor"))
        return tensor_product(first...)
    elseif length(rest) == 1 && rest[1] isa Number
        copies = rest[1]
        copies isa Integer ||
            throw(ArgumentError("the Tensor copy count must be an integer"))
        return tensor_power(first, copies)
    end
    return tensor_product(first, rest...)
end

function _compat_term_count(factor)
    if factor isa AbstractMatrix
        return size(factor, 2)
    elseif factor isa AbstractVector{<:Number}
        return 1
    elseif factor isa Tuple || factor isa AbstractVector
        return length(factor)
    end
    return throw(ArgumentError("invalid TensorSum factor $(typeof(factor))"))
end

function _looks_like_weights(first, second)
    term_count = _compat_term_count(second)
    term_count == 1 && return false
    if first isa AbstractVector{<:Number}
        return length(first) == term_count
    elseif first isa AbstractMatrix{<:Number} && min(size(first)...) == 1
        return length(first) == term_count
    end
    return false
end

"""
    TensorSum(A1, A2, ...)
    TensorSum(weights, A1, A2, ...)

QETLAB-compatible reconstruction from a tensor decomposition.
"""
function TensorSum(first, second, rest...)
    if _looks_like_weights(first, second)
        first isa AbstractArray && Base.require_one_based_indexing(first)
        weights = vec(first)
        return tensor_sum(second, rest...; weights=weights)
    end
    return tensor_sum(first, second, rest...)
end

"""QETLAB-compatible Kronecker sum."""
function KroneckerSum(first, rest...)
    if (first isa Tuple || (first isa AbstractVector && !(eltype(first) <: Number))) &&
        isempty(rest)
        isempty(first) && throw(ArgumentError("KroneckerSum requires at least one factor"))
        return kronecker_sum(first...)
    elseif length(rest) == 1 && rest[1] isa Number
        copies = rest[1]
        copies isa Integer ||
            throw(ArgumentError("the KroneckerSum copy count must be an integer"))
        return kronecker_sum(first; copies=copies)
    end
    return kronecker_sum(first, rest...)
end

"""
    PermuteSystems(X, PERM, DIM=nothing, ROW_ONLY=0, INV_PERM=0)

Compatibility wrapper for QETLAB `PermuteSystems`, including independent
row/column dimensions and row-only rectangular operation.
"""
function PermuteSystems(
    x::AbstractVector, permutation, dim=nothing, row_only=0, inverse_permutation=0
)
    _flag(row_only, "ROW_ONLY") # QETLAB ignores this flag for vectors.
    inverse = _flag(inverse_permutation, "INV_PERM")
    subsystem_count = length(permutation)
    dims = if dim === nothing
        _equal_default_dimensions(length(x), subsystem_count)
    elseif dim isa Integer
        (dim,)
    elseif dim isa AbstractMatrix
        row_dims, column_dims = _matrix_dimension_rows(dim, subsystem_count)
        prod(row_dims) == length(x) ? row_dims : column_dims
    else
        _dimension_tuple(dim)
    end
    plan = SubsystemPermutationPlan(dims, permutation; inverse=inverse)
    return permute_subsystems(x, plan)
end

function PermuteSystems(
    x::AbstractMatrix, permutation, dim=nothing, row_only=0, inverse_permutation=0
)
    if min(size(x)...) == 1
        oriented_dim = if dim isa AbstractMatrix
            row_dims, column_dims = _matrix_dimension_rows(dim, length(permutation))
            size(x, 1) == 1 ? column_dims : row_dims
        else
            dim
        end
        result = PermuteSystems(
            _state_vector(x), permutation, oriented_dim, row_only, inverse_permutation
        )
        return _restore_state_orientation(result, x)
    end

    rows_only = _flag(row_only, "ROW_ONLY")
    inverse = _flag(inverse_permutation, "INV_PERM")
    subsystem_count = length(permutation)

    row_dims, column_dims = if dim === nothing
        (
            _equal_default_dimensions(size(x, 1), subsystem_count),
            _equal_default_dimensions(size(x, 2), subsystem_count),
        )
    elseif dim isa AbstractMatrix
        _matrix_dimension_rows(dim, subsystem_count)
    else
        dims = _dimension_tuple(dim)
        (dims, dims)
    end

    row_plan = SubsystemPermutationPlan(row_dims, permutation; inverse=inverse)
    if rows_only
        return permute_subsystems(x, row_plan; rows_only=true)
    end
    column_plan = SubsystemPermutationPlan(column_dims, permutation; inverse=inverse)
    return permute_subsystems(x, row_plan, column_plan)
end

"""
    PermutationOperator(DIM, PERM, INV_PERM=0, SP=0)

QETLAB-compatible subsystem permutation operator.  QETLAB's default dense
output is retained; set `SP=1` for the sparse representation.
"""
function PermutationOperator(dim, permutation, inverse_permutation=0, sparse_output=0)
    inverse = _flag(inverse_permutation, "INV_PERM")
    use_sparse = _flag(sparse_output, "SP")
    dims = if dim isa Integer
        local_dimension = _positive_dimension(dim, "DIM")
        ntuple(_ -> local_dimension, maximum(permutation))
    else
        _dimension_tuple(dim)
    end
    return permutation_operator(
        dims, permutation; inverse=inverse, sparse_output=use_sparse
    )
end

"""
    Swap(X, SYS=(1,2), DIM=nothing, ROW_ONLY=0)

Compatibility wrapper for swapping two subsystem positions.
"""
function Swap(
    x::Union{AbstractVector,AbstractMatrix}, systems=(1, 2), dim=nothing, row_only=0
)
    system_tuple = Tuple(systems)
    length(system_tuple) == 2 ||
        throw(ArgumentError("SYS must contain exactly two subsystem indices"))
    system_tuple[1] != system_tuple[2] ||
        throw(ArgumentError("SYS must contain two distinct subsystem indices"))

    dims = if dim === nothing
        if x isa AbstractVector
            local_dimension = round(Int, sqrt(length(x)))
            (local_dimension, local_dimension)
        else
            row_local = round(Int, sqrt(size(x, 1)))
            column_local = round(Int, sqrt(size(x, 2)))
            [row_local row_local; column_local column_local]
        end
    elseif dim isa Integer
        if x isa AbstractVector
            _expand_scalar_dimension(dim, length(x), "DIM")
        else
            row_dims = _expand_scalar_dimension(dim, size(x, 1), "DIM")
            column_dims = _expand_scalar_dimension(dim, size(x, 2), "DIM")
            [row_dims[1] row_dims[2]; column_dims[1] column_dims[2]]
        end
    else
        dim
    end

    subsystem_count = dims isa AbstractMatrix ? size(dims, 2) : length(dims)
    permutation = collect(1:subsystem_count)
    first, second = Int(system_tuple[1]), Int(system_tuple[2])
    1 <= first <= subsystem_count && 1 <= second <= subsystem_count ||
        throw(ArgumentError("SYS entries must be between 1 and $subsystem_count"))
    permutation[first], permutation[second] = permutation[second], permutation[first]
    return PermuteSystems(x, Tuple(permutation), dims, row_only, 0)
end

"""
    SwapOperator(DIM, SP=0)

QETLAB-compatible two-party swap operator.
"""
function SwapOperator(dim, sparse_output=0)
    use_sparse = _flag(sparse_output, "SP")
    dims = dim isa Integer ? (dim, dim) : _dimension_tuple(dim)
    length(dims) == 2 ||
        throw(DimensionMismatch("SwapOperator DIM must contain two subsystems"))
    return swap_operator(dims; sparse_output=use_sparse)
end

"""
    PartialTrace(X, SYS=2, DIM=nothing, MODE=-1)

QETLAB-compatible partial trace.  `MODE` is validated but the Julia kernel
selects its sparse or dense implementation by dispatch.
"""
function PartialTrace(x, systems=2, dim=nothing, mode=-1)
    mode in (-1, 0, 1) || throw(ArgumentError("MODE must be -1, 0, or 1; got $mode"))
    pure_state = x isa AbstractVector || (x isa AbstractMatrix && min(size(x)...) == 1)
    total = pure_state ? length(x) : size(x, 1)
    dims = if dim === nothing
        _expand_scalar_dimension(round(Int, sqrt(total)), total, "DIM")
    elseif dim isa Integer
        _expand_scalar_dimension(dim, total, "DIM")
    elseif dim isa AbstractMatrix
        throw(ArgumentError("PartialTrace requires one dimension vector, not a matrix"))
    else
        _dimension_tuple(dim)
    end
    input = x isa AbstractMatrix && pure_state ? _state_vector(x) : x
    return partial_trace(input, dims; trace_out=systems)
end

"""
    PartialTranspose(X, SYS=2, DIM=nothing)

QETLAB-compatible partial transpose, including rectangular row/column
dimension matrices.
"""
function PartialTranspose(matrix::AbstractMatrix, systems=2, dim=nothing)
    row_dims, column_dims = if dim === nothing
        row_local = round(Int, sqrt(size(matrix, 1)))
        column_local = round(Int, sqrt(size(matrix, 2)))
        ((row_local, row_local), (column_local, column_local))
    elseif dim isa Integer
        size(matrix, 1) == size(matrix, 2) || throw(
            ArgumentError("scalar DIM for PartialTranspose requires a square matrix"),
        )
        dims = _expand_scalar_dimension(dim, size(matrix, 1), "DIM")
        (dims, dims)
    elseif dim isa AbstractMatrix
        _matrix_dimension_rows(dim)
    else
        size(matrix, 1) == size(matrix, 2) || throw(
            ArgumentError("vector DIM for PartialTranspose requires a square matrix"),
        )
        dims = _dimension_tuple(dim)
        (dims, dims)
    end
    return partial_transpose(matrix, row_dims, column_dims; systems=systems)
end

"""
    Realignment(X, DIM=nothing)

QETLAB realignment `|ij><kl| ↦ |ik><jl|`.
"""
function Realignment(matrix::AbstractMatrix, dim=nothing)
    row_dims, column_dims = if dim === nothing
        row_local = round(Int, sqrt(size(matrix, 1)))
        column_local = round(Int, sqrt(size(matrix, 2)))
        ((row_local, row_local), (column_local, column_local))
    elseif dim isa Integer
        size(matrix, 1) == size(matrix, 2) ||
            throw(ArgumentError("scalar DIM for Realignment requires a square matrix"))
        dims = _expand_scalar_dimension(dim, size(matrix, 1), "DIM")
        (dims, dims)
    elseif dim isa AbstractMatrix
        _matrix_dimension_rows(dim, 2)
    else
        size(matrix, 1) == size(matrix, 2) ||
            throw(ArgumentError("vector DIM for Realignment requires a square matrix"))
        dims = _dimension_tuple(dim)
        length(dims) == 2 ||
            throw(DimensionMismatch("Realignment DIM must contain two subsystems"))
        (dims, dims)
    end
    return realign(matrix, row_dims, column_dims; systems=(1,))
end

"""
    InverseRealignment(X, DIM)

Julia migration convenience that inverts `Realignment`; `DIM` describes the
original operator.
"""
function InverseRealignment(matrix::AbstractMatrix, dim)
    row_dims, column_dims = if dim isa AbstractMatrix
        _matrix_dimension_rows(dim, 2)
    else
        dims = _dimension_tuple(dim)
        length(dims) == 2 || throw(
            DimensionMismatch("InverseRealignment DIM must contain two subsystems")
        )
        (dims, dims)
    end
    return inverse_realign(matrix, row_dims, column_dims; systems=(1,))
end

function _projection_mode(mode)
    mode in (-1, 0, 1) || throw(ArgumentError("MODE must be -1, 0, or 1; got $mode"))
    return nothing
end

"""QETLAB-compatible symmetric projection or partial isometry."""
function SymmetricProjection(dim, copies=2, partial=0, mode=-1)
    _projection_mode(mode)
    return if _flag(partial, "PARTIAL")
        symmetric_subspace_basis(dim, copies; sparse_output=true)
    else
        symmetric_projector(dim, copies; sparse_output=true)
    end
end

"""QETLAB-compatible antisymmetric projection or partial isometry."""
function AntisymmetricProjection(dim, copies=2, partial=0, mode=-1)
    _projection_mode(mode)
    return if _flag(partial, "PARTIAL")
        antisymmetric_subspace_basis(dim, copies; sparse_output=true)
    else
        antisymmetric_projector(dim, copies; sparse_output=true)
    end
end

"""Compatibility spelling for [`basis_to_linear`](@ref)."""
BasisToLinear(indices, dims) = basis_to_linear(indices, dims)

"""Compatibility spelling for [`linear_to_basis`](@ref)."""
LinearToBasis(index, dims) = linear_to_basis(index, dims)

"""QETLAB-compatible Pauli constructor."""
Pauli(index, sparse_output=1) = pauli(index; sparse_output=_flag(sparse_output, "SP"))

"""QETLAB-compatible generalized Pauli constructor."""
function GenPauli(first, second, dim, sparse_output=0)
    return generalized_pauli(first, second, dim; sparse_output=_flag(sparse_output, "SP"))
end

"""QETLAB-compatible three-dimensional Gell-Mann constructor."""
function GellMann(index, sparse_output=0)
    return gell_mann(index; sparse_output=_flag(sparse_output, "SP"))
end

"""QETLAB-compatible generalized Gell-Mann constructor."""
function GenGellMann(first, second, dim, sparse_output=0)
    return generalized_gell_mann(
        first, second, dim; sparse_output=_flag(sparse_output, "SP")
    )
end

"""QETLAB-compatible quantum Fourier matrix."""
FourierMatrix(dim) = fourier_matrix(dim)

"""QETLAB-compatible maximally entangled vector."""
function MaxEntangled(dim, sparse_output=0, normalized=1)
    return maximally_entangled(
        dim; sparse_output=_flag(sparse_output, "SP"), normalized=_flag(normalized, "NRML")
    )
end

"""QETLAB-compatible Bell-state constructor."""
function Bell(index=0, sparse_output=0, normalized=1)
    index isa Integer && !(index isa Bool) ||
        throw(ArgumentError("IND must be an integer; got $(repr(index))"))
    return bell_state(
        mod(Int(index), 4);
        sparse_output=_flag(sparse_output, "SP"),
        normalized=_flag(normalized, "NRML"),
    )
end

"""QETLAB-compatible generalized GHZ-state constructor."""
function GHZState(dim, parties, coefficients=nothing)
    return ghz_state(dim, parties; coefficients=coefficients, sparse_output=true)
end

"""QETLAB-compatible W-state constructor."""
function WState(parties, coefficients=nothing)
    return w_state(parties; coefficients=coefficients, sparse_output=true)
end

"""QETLAB-compatible Dicke-state constructor."""
function DickeState(parties, excitations=1, normalized=1)
    return dicke_state(
        parties, excitations; normalized=_flag(normalized, "NRML"), sparse_output=true
    )
end

"""QETLAB-compatible isotropic-state constructor."""
IsotropicState(dim, alpha) = isotropic_state(dim, alpha; sparse_output=true)

"""
QETLAB-compatible bipartite Werner-state constructor.

The upstream multipartite vector-parameter form is not exposed until its
normalization and permutation semantics receive a separate verification.
"""
function WernerState(dim, alpha)
    alpha isa Real ||
        throw(ArgumentError("only the verified scalar bipartite ALPHA form is implemented"))
    return werner_state(dim, alpha; sparse_output=true)
end

"""QETLAB-compatible Horodecki-state constructor."""
HorodeckiState(a, dims=(3, 3)) = horodecki_state(a; dims=dims)

"""QETLAB-compatible Gisin-state constructor."""
GisinState(lambda, theta) = gisin_state(lambda, theta)

"""QETLAB-compatible Breuer-state constructor."""
BreuerState(dim, lambda) = breuer_state(dim, lambda; sparse_output=true)

"""
QETLAB-compatible matrix of unnormalized Brauer states.

The Julia-only `T`, `max_matchings`, and `max_nonzeros` keywords are forwarded
to [`brauer_states`](@ref); the resource guards retain the native defaults.
"""
function BrauerStates(
    dim, pairs; T::Type{<:Number}=Float64, max_matchings=100_000, max_nonzeros=1_000_000
)
    return brauer_states(
        dim, pairs; T=T, max_matchings=max_matchings, max_nonzeros=max_nonzeros
    )
end

"""QETLAB-compatible chessboard-state constructor."""
function ChessboardState(a, b, c, d, m, n, s=nothing, t=nothing)
    return chessboard_state(a, b, c, d, m, n; s=s, t=t)
end

"""
QETLAB-compatible random probability vector with mandatory explicit RNG.
"""
RandomProbabilities(rng::AbstractRNG, n) = random_probabilities(rng, n)

"""
QETLAB-compatible random state vector with mandatory explicit RNG.
"""
function RandomStateVector(rng::AbstractRNG, dim, real_output=0, schmidt_rank=0)
    rank_value = _nonnegative_dimension(schmidt_rank, "K")
    return random_state_vector(
        rng,
        dim;
        real=_flag(real_output, "RE"),
        schmidt_rank=rank_value == 0 ? nothing : rank_value,
    )
end

"""
QETLAB-compatible random density matrix with mandatory explicit RNG.
"""
function RandomDensityMatrix(
    rng::AbstractRNG, dim, real_output=0, rank=dim, distribution="haar"
)
    return random_density_matrix(
        rng, dim; real=_flag(real_output, "RE"), rank=rank, distribution=distribution
    )
end

"""
QETLAB-compatible Haar unitary with mandatory explicit RNG.
"""
function RandomUnitary(rng::AbstractRNG, dim, real_output=0)
    return random_unitary(rng, dim; real=_flag(real_output, "RE"))
end

"""
QETLAB-compatible random graph with mandatory explicit RNG.
"""
function RandomGraph(rng::AbstractRNG, n, edge_probability=0.5)
    return Float64.(random_graph(rng, n; edge_probability=edge_probability))
end

"""
QETLAB-compatible random POVM with mandatory explicit RNG.
"""
function RandomPOVM(rng::AbstractRNG, dim, outcomes, real_output=0)
    return random_povm(rng, dim, outcomes; real=_flag(real_output, "RE"))
end

function _compat_map_dimensions(dim)
    dimensions = if dim isa Integer
        value = _positive_dimension(dim, "DIM")
        (value, value)
    elseif dim isa AbstractMatrix
        Base.require_one_based_indexing(dim)
        size(dim) == (2, 2) || throw(
            DimensionMismatch("a map DIM matrix must have size (2, 2); got $(size(dim))"),
        )
        row_dimensions = Tuple(dim[1, :])
        column_dimensions = Tuple(dim[2, :])
        row_dimensions == column_dimensions || throw(
            ArgumentError(
                "rectangular operator-space maps are not supported by this " *
                "compatibility slice; DIM row and column dimensions must agree",
            ),
        )
        row_dimensions
    elseif dim isa Tuple || dim isa AbstractVector
        Tuple(dim)
    else
        throw(
            ArgumentError(
                "DIM must be an integer, a two-entry tuple/vector, or a 2-by-2 matrix"
            ),
        )
    end
    length(dimensions) == 2 || throw(
        DimensionMismatch(
            "map DIM must contain input and output dimensions; got " *
            "$(length(dimensions)) entries",
        ),
    )
    return (
        _positive_dimension(dimensions[1], "DIM[1]"),
        _positive_dimension(dimensions[2], "DIM[2]"),
    )
end

function _validate_compat_map_dimensions(map, dim)
    dim === nothing && return map
    input_dim, output_dim = _compat_map_dimensions(dim)
    input_dimension(map) == input_dim || throw(
        DimensionMismatch(
            "DIM input dimension $input_dim does not match the map input " *
            "dimension $(input_dimension(map))",
        ),
    )
    output_dimension(map) == output_dim || throw(
        DimensionMismatch(
            "DIM output dimension $output_dim does not match the map output " *
            "dimension $(output_dimension(map))",
        ),
    )
    return map
end

function _compat_choi_dimensions(matrix, dim, input_hint)
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("a Choi matrix must be square; got size $(size(matrix))"))
    total_dimension = size(matrix, 1)
    total_dimension > 0 ||
        throw(ArgumentError("a Choi matrix must have nonzero dimensions"))
    if dim !== nothing
        input_dim, output_dim = _compat_map_dimensions(dim)
    elseif input_hint !== nothing
        input_dim = _positive_dimension(input_hint, "input dimension")
        rem(total_dimension, input_dim) == 0 || throw(
            DimensionMismatch(
                "Choi dimension $total_dimension is not divisible by inferred " *
                "input dimension $input_dim",
            ),
        )
        output_dim = div(total_dimension, input_dim)
    else
        input_dim = isqrt(total_dimension)
        input_dim^2 == total_dimension || throw(
            ArgumentError(
                "cannot infer unequal input/output dimensions from Choi size " *
                "$(size(matrix)); provide DIM=(input, output)",
            ),
        )
        output_dim = input_dim
    end
    input_dim * output_dim == total_dimension || throw(
        DimensionMismatch(
            "DIM=($input_dim, $output_dim) is inconsistent with Choi size " *
            "$(size(matrix))",
        ),
    )
    return input_dim, output_dim
end

function _is_cp_kraus_collection(value)
    (value isa Tuple || value isa AbstractVector) || return false
    isempty(value) && return false
    return all(operator -> operator isa AbstractMatrix && eltype(operator) <: Number, value)
end

function _compat_map(phi, dim=nothing; input_hint=nothing)
    if phi isa AbstractMapRepresentation
        return _validate_compat_map_dimensions(phi, dim)
    elseif phi isa AbstractMatrix && eltype(phi) <: Number
        input_dim, output_dim = _compat_choi_dimensions(phi, dim, input_hint)
        return ChoiRepresentation(phi, input_dim, output_dim)
    elseif _is_cp_kraus_collection(phi)
        return _validate_compat_map_dimensions(KrausRepresentation(phi), dim)
    end
    return throw(
        ArgumentError(
            "PHI must be a numeric Choi matrix, a nonempty vector/tuple of " *
            "completely-positive Kraus matrices, or an AbstractMapRepresentation; " *
            "QETLAB two-column left/right Kraus cells are not yet supported",
        ),
    )
end

"""
    ApplyMap(X, PHI)

QETLAB-compatible application of a map supplied as a numeric Choi matrix or
a vector/tuple of completely-positive Kraus matrices. Choi input/output
dimensions are inferred from `size(X)` when they are unequal. QETLAB's
two-column left/right Kraus-cell representation is intentionally rejected
until a native two-sided representation is available.
"""
function ApplyMap(input::AbstractMatrix, phi)
    size(input, 1) == size(input, 2) || throw(
        DimensionMismatch(
            "ApplyMap currently requires a square input matrix; got $(size(input))"
        ),
    )
    map = _compat_map(phi; input_hint=size(input, 1))
    return apply_channel(input, map)
end

"""
    ChoiMatrix(PHI, SYS=2)

Return the raw Choi matrix of `PHI`. Numeric Choi matrices are copied exactly,
matching QETLAB's early-return behavior. For Kraus input, `SYS=2` uses the
package's input-first convention and `SYS=1` swaps the input/output tensor
factors.
"""
function ChoiMatrix(phi, system=2)
    if phi isa AbstractMatrix && eltype(phi) <: Number
        return copy(phi)
    end
    system isa Integer && !(system isa Bool) && system in (1, 2) ||
        throw(ArgumentError("SYS must be 1 or 2; got $(repr(system))"))
    map = _compat_map(phi)
    matrix = choi_matrix(map)
    system == 2 && return matrix
    return permute_subsystems(
        matrix, (input_dimension(map), output_dimension(map)); permutation=(2, 1)
    )
end

"""
    KrausOperators(PHI, DIM=nothing)

Return canonical completely-positive Kraus matrices. A numeric `PHI` is
interpreted as a Choi matrix; `DIM=(input, output)` is required for unequal
dimensions. General Hermiticity-preserving/non-CP QETLAB two-sided outputs are
not fabricated: such inputs raise `DomainError`.
"""
function KrausOperators(phi, dim=nothing)
    return kraus_operators(_compat_map(phi, dim))
end

"""
    ComplementaryMap(PHI, DIM=nothing)

Return a complementary map, preserving QETLAB's raw representation kind:
Kraus collections produce a vector of matrices and numeric Choi matrices
produce a numeric Choi matrix. The input must be completely positive.
"""
function ComplementaryMap(phi, dim=nothing)
    map = _compat_map(phi, dim)
    complement = complementary_channel(map)
    if phi isa AbstractMapRepresentation
        return complement
    elseif _is_cp_kraus_collection(phi)
        return kraus_operators(complement)
    end
    return choi_matrix(complement)
end

"""
    DualMap(PHI, DIM=nothing)

QETLAB-compatible Hilbert--Schmidt dual. Raw Kraus input returns raw adjoint
Kraus matrices, raw Choi input returns a raw Choi matrix, and project-native
representations retain their representation kind.
"""
function DualMap(phi, dim=nothing)
    map = _compat_map(phi, dim)
    dual = dual_channel(map)
    if phi isa AbstractMapRepresentation
        return dual
    elseif _is_cp_kraus_collection(phi)
        return kraus_operators(dual)
    end
    return choi_matrix(dual)
end

function _compat_operator_dimensions(input::AbstractMatrix, dim)
    size(input, 1) == size(input, 2) || throw(
        DimensionMismatch(
            "PartialMap currently requires a square input matrix; got $(size(input))"
        ),
    )
    total_dimension = size(input, 1)
    if dim === nothing
        local_dimension = isqrt(total_dimension)
        local_dimension^2 == total_dimension || throw(
            ArgumentError(
                "cannot infer two equal subsystem dimensions from input size " *
                "$(size(input)); provide DIM",
            ),
        )
        return (local_dimension, local_dimension)
    elseif dim isa Integer
        return _expand_scalar_dimension(dim, total_dimension, "DIM")
    elseif dim isa AbstractMatrix
        Base.require_one_based_indexing(dim)
        size(dim, 1) == 2 || throw(
            DimensionMismatch(
                "a PartialMap DIM matrix must have two rows; got $(size(dim))"
            ),
        )
        row_dimensions = Tuple(dim[1, :])
        column_dimensions = Tuple(dim[2, :])
        row_dimensions == column_dimensions || throw(
            ArgumentError(
                "rectangular row/column subsystem dimensions are not supported " *
                "by this PartialMap compatibility slice",
            ),
        )
        return row_dimensions
    end
    return _dimension_tuple(dim)
end

"""
    PartialMap(X, PHI, SYS=2, DIM=nothing)

Apply `PHI` to subsystem `SYS` using QETLAB argument order. Square operator
spaces and vector dimensions are supported. QETLAB's rectangular row/column
`DIM` form is rejected explicitly rather than silently applying the wrong
convention.
"""
function PartialMap(input::AbstractMatrix, phi, system=2, dim=nothing)
    system isa Integer && !(system isa Bool) ||
        throw(ArgumentError("SYS must be an integer; got $(repr(system))"))
    dimensions = _compat_operator_dimensions(input, dim)
    1 <= system <= length(dimensions) ||
        throw(ArgumentError("SYS must be between 1 and $(length(dimensions)); got $system"))
    map = _compat_map(phi; input_hint=dimensions[system])
    return partial_map(input, map, system, dimensions)
end

function _compat_finite_real(value, name::AbstractString)
    value isa Real ||
        throw(ArgumentError("$name must be a real number; got $(repr(value))"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    return value
end

function _compat_unnormalized_entangled(dimension::Int, ::Type{T}) where {T}
    total = Base.checked_mul(dimension, dimension)
    indices = 1:(dimension + 1):total
    return sparsevec(indices, fill(one(T), dimension), total)
end

"""
    DepolarizingChannel(DIM, P=0)

Return QETLAB's raw Choi matrix
`(1-P)I/DIM + P|Ω><Ω|`. Unlike the physical-channel native constructor, this
compatibility wrapper preserves QETLAB's finite out-of-CP-range parameter
semantics.
"""
function DepolarizingChannel(dim, p=0)
    dimension = _positive_dimension(dim, "DIM")
    parameter = _compat_finite_real(p, "P")
    scalar_type = promote_type(
        typeof((one(parameter) - parameter) / dimension), typeof(parameter)
    )
    omega = _compat_unnormalized_entangled(dimension, scalar_type)
    total = Base.checked_mul(dimension, dimension)
    diagonal_weight = (one(parameter) - parameter) / dimension
    return spdiagm(0 => fill(convert(scalar_type, diagonal_weight), total)) +
           parameter * (omega * adjoint(omega))
end

"""
    DephasingChannel(DIM, P=0)

Return QETLAB's raw dephasing-channel Choi matrix. Finite values outside the
completely-positive range remain representable for migration compatibility.
"""
function DephasingChannel(dim, p=0)
    dimension = _positive_dimension(dim, "DIM")
    parameter = _compat_finite_real(p, "P")
    scalar_type = promote_type(typeof(one(parameter) - parameter), typeof(parameter))
    omega = _compat_unnormalized_entangled(dimension, scalar_type)
    coherent = omega * adjoint(omega)
    return (one(parameter) - parameter) * spdiagm(0 => diag(coherent)) +
           parameter * coherent
end

function _compat_power(base::Int, exponent::Int, name::AbstractString)
    result = 1
    for _ in 1:exponent
        result = try
            Base.checked_mul(result, base)
        catch error
            error isa OverflowError || rethrow()
            throw(ArgumentError("$name exceeds typemax(Int)"))
        end
    end
    return result
end

"""
    PauliChannel(P)
    PauliChannel(rng, Q)

Return a raw sparse Choi matrix. A probability array follows QETLAB ordering.
The scalar random-channel form requires an explicit leading RNG so the
caller's global random stream is never mutated.
"""
PauliChannel(probabilities::AbstractArray) = choi_matrix(pauli_channel(probabilities))

function PauliChannel(rng::AbstractRNG, qubits)
    qubit_count = _positive_dimension(qubits, "Q")
    probability_count = _compat_power(4, qubit_count, "4^Q")
    probabilities = random_probabilities(rng, probability_count)
    return choi_matrix(pauli_channel(probabilities))
end

"""Return QETLAB's generalized Choi-map matrix."""
ChoiMap(a=1, b=1, c=0) = choi_matrix(choi_map(a, b, c))

"""Return QETLAB's reduction-map matrix."""
ReductionMap(dim, k=1) = choi_matrix(reduction_map(dim, k))

"""
    TraceNorm(X; allow_densify=false)

QETLAB-compatible trace norm. Sparse input is never densified implicitly;
pass `allow_densify=true` after reviewing the allocation.
"""
function TraceNorm(matrix; allow_densify::Bool=false)
    return trace_norm(matrix; allow_densify=allow_densify)
end

"""
    SchattenNorm(X, P; allow_densify=false)

QETLAB-compatible Schatten `P`-norm with an explicit sparse-densification
gate.
"""
function SchattenNorm(matrix, p; allow_densify::Bool=false)
    return schatten_norm(matrix, p; allow_densify=allow_densify)
end

"""
    KyFanNorm(X, K; allow_densify=false)

QETLAB-compatible Ky Fan `K`-norm with an explicit sparse-densification gate.
"""
function KyFanNorm(matrix, k; allow_densify::Bool=false)
    return ky_fan_norm(matrix, k; allow_densify=allow_densify)
end

"""
    Purity(RHO)

Return `real(tr(RHO^2))`, matching QETLAB's unchecked scalar operation. The
Julia-native [`QuantumEntanglementTools.purity`](@ref) should be preferred when density-matrix
validation is required.
"""
function Purity(rho::AbstractMatrix{<:Number})
    Base.require_one_based_indexing(rho)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("RHO must be square; got size $(size(rho))"))
    all(isfinite, rho) || throw(ArgumentError("RHO must contain only finite entries"))
    return real(tr(rho * rho))
end

"""
    Entropy(RHO, BASE=2, ALPHA=1; allow_densify=false)

Compatibility entry point for the verified von Neumann (`ALPHA == 1`) branch
of QETLAB `Entropy`. Rényi orders are rejected explicitly until the native
entropy API covers and validates them.
"""
function Entropy(rho::AbstractMatrix{<:Number}, base=2, alpha=1; allow_densify::Bool=false)
    alpha isa Real && !(alpha isa Bool) ||
        throw(ArgumentError("ALPHA must be a real number"))
    alpha == one(alpha) ||
        throw(ArgumentError("only the verified von Neumann ALPHA=1 branch is implemented"))
    return von_neumann_entropy(rho; base=base, allow_densify=allow_densify)
end

"""
    Fidelity(RHO, SIGMA; allow_densify=false)

Return QETLAB's unsquared Uhlmann root fidelity. Inputs receive the stricter
Julia-native density-matrix validation and are never repaired or normalized.
"""
function Fidelity(rho, sigma; allow_densify::Bool=false)
    return fidelity(rho, sigma; allow_densify=allow_densify)
end

function _compat_bipartite_dimensions(state, dim, name::AbstractString)
    total_dimension = if state isa AbstractVector
        length(state)
    elseif state isa AbstractMatrix
        size(state, 1) == size(state, 2) || throw(
            DimensionMismatch("$name requires a square matrix; got size $(size(state))"),
        )
        size(state, 1)
    else
        throw(ArgumentError("$name input must be a state vector or square matrix"))
    end
    dimensions = if dim === nothing
        local_dimension = isqrt(total_dimension)
        local_dimension^2 == total_dimension || throw(
            ArgumentError(
                "cannot infer two equal subsystem dimensions from total " *
                "dimension $total_dimension; provide DIM",
            ),
        )
        (local_dimension, local_dimension)
    elseif dim isa Integer
        _expand_scalar_dimension(dim, total_dimension, "DIM")
    else
        _dimension_tuple(dim)
    end
    length(dimensions) == 2 ||
        throw(DimensionMismatch("$name DIM must contain exactly two subsystem dimensions"))
    checked = (
        _positive_dimension(dimensions[1], "DIM[1]"),
        _positive_dimension(dimensions[2], "DIM[2]"),
    )
    checked[1] * checked[2] == total_dimension || throw(
        DimensionMismatch(
            "prod(DIM)=$(checked[1] * checked[2]) does not match state " *
            "dimension $total_dimension",
        ),
    )
    return checked
end

"""
    Negativity(RHO, DIM=nothing; allow_densify=false)

QETLAB-compatible bipartite negativity. The Julia implementation validates
normalization and positivity and requires explicit permission before
densifying sparse density matrices.
"""
function Negativity(state, dim=nothing; allow_densify::Bool=false)
    dimensions = _compat_bipartite_dimensions(state, dim, "Negativity")
    return negativity(state, dimensions; systems=(2,), allow_densify=allow_densify)
end

"""
    SchmidtDecomposition(VEC, DIM=nothing, K=0; allow_densify=false)

Return a named tuple `(coefficients, left_vectors, right_vectors)`. `K=0`
keeps QETLAB's numerically nonzero terms, `K=-1` keeps the full thin
decomposition, and positive `K` keeps that many leading terms.
"""
function SchmidtDecomposition(
    state::AbstractVector{<:Number}, dim=nothing, k=0; allow_densify::Bool=false
)
    dimensions = _compat_bipartite_dimensions(state, dim, "SchmidtDecomposition")
    k isa Integer && !(k isa Bool) ||
        throw(ArgumentError("K must be -1, 0, or a positive integer"))
    k >= -1 || throw(ArgumentError("K must be -1, 0, or a positive integer"))
    decomposition = schmidt_decomposition(state, dimensions; allow_densify=allow_densify)
    coefficient_count = length(decomposition.coefficients)
    retained = if k == -1
        coefficient_count
    elseif k == 0
        if coefficient_count == 0
            0
        else
            first_coefficient = first(decomposition.coefficients)
            threshold = maximum(dimensions) * eps(first_coefficient)
            count(value -> value > threshold, decomposition.coefficients)
        end
    else
        k <= coefficient_count ||
            throw(ArgumentError("K=$k exceeds min(DIM)=$coefficient_count"))
        Int(k)
    end
    indices = 1:retained
    return (
        coefficients=decomposition.coefficients[indices],
        left_vectors=decomposition.left_vectors[:, indices],
        right_vectors=decomposition.right_vectors[:, indices],
    )
end

"""
    SchmidtRank(VEC, DIM=nothing, TOL=nothing; allow_densify=false)

Return QETLAB's tolerance-defined bipartite Schmidt rank. Supplying `TOL`
uses it as a strict absolute singular-value threshold.
"""
function SchmidtRank(
    state::AbstractVector{<:Number},
    dim=nothing,
    tolerance=nothing;
    allow_densify::Bool=false,
)
    dimensions = _compat_bipartite_dimensions(state, dim, "SchmidtRank")
    if tolerance === nothing
        _tolerance = sqrt(length(state)) * eps(norm(state))
    else
        _tolerance = _compat_finite_real(tolerance, "TOL")
        _tolerance >= 0 || throw(ArgumentError("TOL must be nonnegative"))
    end
    return schmidt_rank(
        state, dimensions; atol=_tolerance, rtol=0, allow_densify=allow_densify
    )
end

function _compat_product_vector_dimensions(
    vector::AbstractVector, dim, name::AbstractString
)
    total_dimension = length(vector)
    dimensions = if dim === nothing
        local_dimension = isqrt(total_dimension)
        local_dimension^2 == total_dimension || throw(
            ArgumentError(
                "cannot infer two equal subsystem dimensions from vector " *
                "length $total_dimension; provide DIM",
            ),
        )
        (local_dimension, local_dimension)
    elseif dim isa Integer
        _expand_scalar_dimension(dim, total_dimension, "DIM")
    else
        supplied = _dimension_tuple(dim)
        if length(supplied) == 1
            _expand_scalar_dimension(first(supplied), total_dimension, "DIM")
        else
            supplied
        end
    end
    length(dimensions) >= 2 ||
        throw(DimensionMismatch("$name DIM must contain at least two subsystems"))
    checked = ntuple(
        index -> _positive_dimension(dimensions[index], "DIM[$index]"), length(dimensions)
    )
    prod(checked) == total_dimension || throw(
        DimensionMismatch(
            "prod(DIM)=$(prod(checked)) does not match vector length " * "$total_dimension",
        ),
    )
    return checked
end

function _compat_product_operator_dimensions(
    operator::AbstractMatrix, dim, name::AbstractString; bipartite::Bool
)
    row_dimensions, column_dimensions = if dim === nothing
        row_local = isqrt(size(operator, 1))
        column_local = isqrt(size(operator, 2))
        row_local^2 == size(operator, 1) || throw(
            ArgumentError(
                "cannot infer two equal row subsystem dimensions from size " *
                "$(size(operator)); provide DIM",
            ),
        )
        column_local^2 == size(operator, 2) || throw(
            ArgumentError(
                "cannot infer two equal column subsystem dimensions from size " *
                "$(size(operator)); provide DIM",
            ),
        )
        ((row_local, row_local), (column_local, column_local))
    elseif dim isa Integer
        size(operator, 1) == size(operator, 2) ||
            throw(ArgumentError("scalar DIM requires a square operator"))
        dimensions = _expand_scalar_dimension(dim, size(operator, 1), "DIM")
        (dimensions, dimensions)
    elseif dim isa AbstractMatrix
        _matrix_dimension_rows(dim)
    else
        size(operator, 1) == size(operator, 2) || throw(
            ArgumentError(
                "vector DIM requires a square operator; use a two-row DIM " *
                "matrix for rectangular row/column dimensions",
            ),
        )
        supplied = _dimension_tuple(dim)
        dimensions = if length(supplied) == 1
            _expand_scalar_dimension(first(supplied), size(operator, 1), "DIM")
        else
            supplied
        end
        (dimensions, dimensions)
    end
    length(row_dimensions) == length(column_dimensions) || throw(
        DimensionMismatch(
            "DIM row and column layouts must contain the same number of subsystems"
        ),
    )
    length(row_dimensions) >= 2 ||
        throw(DimensionMismatch("$name DIM must contain at least two subsystems"))
    if bipartite && length(row_dimensions) != 2
        throw(DimensionMismatch("$name DIM must contain exactly two subsystems"))
    end
    checked_rows = ntuple(
        index -> _positive_dimension(row_dimensions[index], "DIM[1,$index]"),
        length(row_dimensions),
    )
    checked_columns = ntuple(
        index -> _positive_dimension(column_dimensions[index], "DIM[2,$index]"),
        length(column_dimensions),
    )
    prod(checked_rows) == size(operator, 1) || throw(
        DimensionMismatch(
            "product of DIM row dimensions $(prod(checked_rows)) does not " *
            "match operator row count $(size(operator, 1))",
        ),
    )
    prod(checked_columns) == size(operator, 2) || throw(
        DimensionMismatch(
            "product of DIM column dimensions $(prod(checked_columns)) does " *
            "not match operator column count $(size(operator, 2))",
        ),
    )
    return checked_rows, checked_columns
end

"""
    OperatorSchmidtDecomposition(X, DIM=nothing, K=0;
                                 allow_densify=false)

Return a named tuple `(coefficients, left_factors, right_factors)`. `K=0`
keeps QETLAB's numerically nonzero terms, `K=-1` keeps the full thin
decomposition, and positive `K` keeps that many leading terms. A two-row
`DIM` matrix supplies independent local row and column dimensions.

The native SVD does not reproduce QETLAB's faulty Hermitian-factor repair
branch. Factors reconstruct the operator and are Frobenius-orthonormal, but
individual factors are not promised Hermitian.
"""
function OperatorSchmidtDecomposition(
    operator::AbstractMatrix{<:Number}, dim=nothing, k=0; allow_densify::Bool=false
)
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        operator, dim, "OperatorSchmidtDecomposition"; bipartite=true
    )
    k isa Integer && !(k isa Bool) ||
        throw(ArgumentError("K must be -1, 0, or a positive integer"))
    k >= -1 || throw(ArgumentError("K must be -1, 0, or a positive integer"))
    decomposition = operator_schmidt_decomposition(
        operator, row_dimensions, column_dimensions; allow_densify=allow_densify
    )
    coefficient_count = length(decomposition.coefficients)
    retained = if k == -1
        coefficient_count
    elseif k == 0
        if coefficient_count == 0
            0
        else
            first_coefficient = first(decomposition.coefficients)
            local_operator_dimensions = (
                row_dimensions[1] * column_dimensions[1],
                row_dimensions[2] * column_dimensions[2],
            )
            threshold = maximum(local_operator_dimensions) * eps(first_coefficient)
            count(value -> value > threshold, decomposition.coefficients)
        end
    else
        k <= coefficient_count || throw(
            ArgumentError(
                "K=$k exceeds the full decomposition length " * "$coefficient_count"
            ),
        )
        Int(k)
    end
    indices = 1:retained
    return (
        coefficients=decomposition.coefficients[indices],
        left_factors=decomposition.left_factors[indices],
        right_factors=decomposition.right_factors[indices],
    )
end

"""
    OperatorSchmidtRank(X, DIM=nothing; allow_densify=false)

Return QETLAB's default-tolerance numerical operator Schmidt rank. A two-row
`DIM` matrix supplies independent rectangular local row and column dimensions.
"""
function OperatorSchmidtRank(
    operator::AbstractMatrix{<:Number}, dim=nothing; allow_densify::Bool=false
)
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        operator, dim, "OperatorSchmidtRank"; bipartite=true
    )
    tolerance = sqrt(length(operator)) * eps(norm(operator))
    return operator_schmidt_rank(
        operator,
        row_dimensions,
        column_dimensions;
        atol=tolerance,
        rtol=0,
        allow_densify=allow_densify,
    )
end

"""
    IsProductVector(VEC, DIM=nothing; atol=nothing, rtol=nothing,
                    allow_densify=false) -> ProductAnalysisResult

QETLAB-compatible argument handling with a structured numerical result.
Clear product/nonproduct classifications map to `:within_tolerance` and
`:outside_tolerance`; a numerical equality boundary remains `:boundary`
instead of being collapsed to a Boolean. Product factors are available in
`result.factors`.
"""
function IsProductVector(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dim=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    vector = state isa AbstractMatrix ? _state_vector(state) : state
    dimensions = _compat_product_vector_dimensions(vector, dim, "IsProductVector")
    return is_product_vector(
        vector, dimensions; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

"""
    IsProductOperator(X, DIM=nothing; atol=nothing, rtol=nothing,
                      allow_densify=false) -> ProductAnalysisResult

Analyze an elementary tensor using QETLAB dimension forms and return the
native structured tolerance result. A two-row `DIM` matrix supports
independently rectangular local operator dimensions.
"""
function IsProductOperator(
    operator::AbstractMatrix{<:Number},
    dim=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        operator, dim, "IsProductOperator"; bipartite=false
    )
    return is_product_operator(
        operator,
        row_dimensions,
        column_dimensions;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
    )
end

"""
    Concurrence(RHO; allow_densify=false)

QETLAB-compatible two-qubit concurrence for a normalized pure vector or
density matrix.
"""
function Concurrence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}};
    allow_densify::Bool=false,
)
    if state isa AbstractVector
        return concurrence(state)
    end
    return concurrence(state; allow_densify=allow_densify)
end

"""
    EntFormation(RHO, DIM=nothing; atol=nothing, rtol=nothing,
                 allow_densify=false)

Return base-two entanglement of formation for normalized bipartite pure
vectors and validated two-qubit density matrices. Row and column pure vectors
are accepted. Unlike QETLAB, a higher-dimensional rank-one density matrix is
not silently converted to a vector; pass its state vector explicitly.
"""
function EntFormation(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dim=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    input =
        state isa AbstractMatrix && min(size(state)...) == 1 ? _state_vector(state) : state
    dimensions = _compat_bipartite_dimensions(input, dim, "EntFormation")
    return entanglement_of_formation(
        input, dimensions; base=2, atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

function _compat_default_ball_dimensions(total_dimension::Int)
    total_dimension > 1 ||
        throw(ArgumentError("InSeparableBall requires dimension greater than one"))
    first_dimension = isqrt(total_dimension)
    while first_dimension > 1 && rem(total_dimension, first_dimension) != 0
        first_dimension -= 1
    end
    return (first_dimension, div(total_dimension, first_dimension))
end

"""
    InSeparableBall(X, DIM=nothing; atol=nothing, rtol=nothing,
                    allow_densify=false) -> SeparableBallResult

Apply QETLAB's trace-normalized Gurvits--Barnum ball check while returning the
native structured sufficient-certificate result. `:outside_ball` means only
that this sufficient condition failed, never that the input is entangled.
Unlike QETLAB, positivity and finite-input validation are enforced. `DIM` is a
Julia migration extension; when omitted, a balanced bipartite factorization
of the total dimension is selected.
"""
function InSeparableBall(
    input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dim=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    Base.require_one_based_indexing(input)
    total_dimension = if input isa AbstractVector
        length(input)
    else
        size(input, 1) == size(input, 2) || throw(
            DimensionMismatch(
                "X must be square or an eigenvalue vector; got size " * "$(size(input))",
            ),
        )
        size(input, 1)
    end
    dimensions = if dim === nothing
        _compat_default_ball_dimensions(total_dimension)
    else
        _compat_bipartite_dimensions(input, dim, "InSeparableBall")
    end
    all(isfinite, input) || throw(ArgumentError("X must contain only finite entries"))
    trace_value = input isa AbstractVector ? sum(input) : tr(input)
    isreal(trace_value) || throw(ArgumentError("X must have a real trace"))
    real_trace = real(trace_value)
    real_trace > zero(real_trace) ||
        throw(DomainError(real_trace, "X must have positive trace"))
    normalized = input / real_trace
    return in_separable_ball(
        normalized, dimensions; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

function _compat_majorization_values(
    input::Union{AbstractVector,AbstractMatrix}, name::AbstractString; allow_densify::Bool
)
    Base.require_one_based_indexing(input)
    if issparse(input) && !allow_densify
        throw(
            ArgumentError(
                "$name is sparse; pass allow_densify=true to permit the " *
                "QETLAB-compatible dense comparison",
            ),
        )
    end
    vector_input =
        input isa AbstractVector || (input isa AbstractMatrix && 1 in size(input))
    values = if vector_input
        collect(vec(input))
    else
        for (index, value) in pairs(input)
            value isa Number || throw(
                ArgumentError("$name entry $index must be numeric; got $(repr(value))"),
            )
            isfinite(value) || throw(
                ArgumentError("$name entry $index must be finite; got $(repr(value))"),
            )
        end
        isempty(input) && return Float64[]
        try
            collect(svdvals(Matrix(input)))
        catch err
            err isa MethodError || rethrow()
            throw(
                ArgumentError(
                    "$name has an element type unsupported by dense singular-value decomposition",
                ),
            )
        end
    end
    for (index, value) in pairs(values)
        value isa Real || throw(
            ArgumentError("$name vector entry $index must be real; got $(repr(value))")
        )
        isfinite(value) || throw(
            ArgumentError("$name vector entry $index must be finite; got $(repr(value))"),
        )
    end
    sort!(values; rev=true)
    return values
end

_compat_majorization_widen(value::Integer) = BigInt(value)
_compat_majorization_widen(value::Rational) = Rational{BigInt}(value)
_compat_majorization_widen(value) = value

function _compat_majorization_zero(values)
    isempty(values) && return 0
    return zero(first(values))
end

function _compat_majorization_tolerance(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a nonnegative real number, not Bool"))
    checked = _compat_finite_real(value, name)
    checked >= 0 || throw(ArgumentError("$name must be nonnegative; got $checked"))
    return checked
end

"""
    Majorizes(A, B; atol=0, rtol=eps(Float64)^(3/4),
               allow_densify=false) -> Bool

Preserve pinned QETLAB's weak-majorization contract. Row and column matrices
are treated as vectors; other matrices are replaced by singular values. Each
value list is sorted before shorter input is padded with trailing zeros, and
all prefix sums are compared without requiring equal totals. The one-sided
comparison allowance is `atol + rtol * norm(A_values)`, matching QETLAB's
default when the keywords are omitted while making the tolerance explicit.

Sparse inputs require `allow_densify=true`. Vector entries must be finite real
numbers; general matrix entries may be finite real or complex numbers. Vector
work costs `O(n log n)` and matrix work additionally performs a dense SVD. Use
native `majorizes` for standard strong majorization with symmetric tolerances,
pad-before-sort semantics, and matrix treatment for every matrix shape.
"""
function Majorizes(
    first::Union{AbstractVector,AbstractMatrix},
    second::Union{AbstractVector,AbstractMatrix};
    atol=0,
    rtol=eps(Float64)^(3 / 4),
    allow_densify::Bool=false,
)
    checked_atol = _compat_majorization_tolerance(atol, "atol")
    checked_rtol = _compat_majorization_tolerance(rtol, "rtol")
    first_values = _compat_majorization_values(first, "A"; allow_densify=allow_densify)
    second_values = _compat_majorization_values(second, "B"; allow_densify=allow_densify)
    tolerance = checked_atol + checked_rtol * norm(first_values)

    common_length = max(length(first_values), length(second_values))
    while length(first_values) < common_length
        push!(first_values, _compat_majorization_zero(first_values))
    end
    while length(second_values) < common_length
        push!(second_values, _compat_majorization_zero(second_values))
    end

    first_prefix = BigInt(0)
    second_prefix = BigInt(0)
    for position in 1:common_length
        first_prefix += _compat_majorization_widen(first_values[position])
        second_prefix += _compat_majorization_widen(second_values[position])
        first_prefix + tolerance < second_prefix && return false
    end
    return true
end

function _compat_symmetric_polynomial_values(values::AbstractVector)
    Base.require_one_based_indexing(values)
    return values
end

function _compat_symmetric_polynomial_values(values::AbstractMatrix)
    Base.require_one_based_indexing(values)
    1 in size(values) || throw(
        DimensionMismatch(
            "X must be a vector or a row/column matrix; got size $(size(values))"
        ),
    )
    return vec(values)
end

"""
    ElemSymPoly(X, K)

Evaluate the `K`th elementary symmetric polynomial. Julia vectors and
MATLAB-shaped row or column matrices are accepted. The wrapper delegates to
`elementary_symmetric_polynomial`, including its `e₀ = 1` convention, sparse
vector support, widened exact integer/rational arithmetic, checked narrowing,
and explicit invalid-order errors. Complexity is `O(length(X) * K)` time and
`O(K)` workspace.

This preserves the pinned QETLAB entry-point spelling while using the native
exact-arithmetic implementation.
"""
function ElemSymPoly(values::Union{AbstractVector,AbstractMatrix}, order)
    return elementary_symmetric_polynomial(
        _compat_symmetric_polynomial_values(values), order
    )
end

"""
    CompoundMatrix(A, R; sparse_output=issparse(A))

Return the `R`th multiplicative compound. Order zero is a `1×1` identity.
When `R > min(size(A)...)`, this compatibility entry point preserves pinned
QETLAB's `0×0` result. Other orders delegate to `compound_matrix`, including
its widened exact minor arithmetic and checked narrowing. Sparse inputs remain
sparse by default; `sparse_output` selects the representation explicitly. Use
the native function to retain a mathematically informative zero-by-nonzero
shape when the order exceeds only one dimension.
"""
function CompoundMatrix(matrix::AbstractMatrix, order; sparse_output::Bool=issparse(matrix))
    checked_order = _nonnegative_dimension(order, "R")
    checked_order > min(size(matrix)...) && return Matrix{Float64}(undef, 0, 0)
    return compound_matrix(matrix, checked_order; sparse_output=sparse_output)
end

"""
    AdditiveCompoundMatrix(A, R; sparse_output=issparse(A))

Return the `R`th additive compound. `A` must be square. Order zero reproduces
the reviewed pinned-QETLAB dependency-path error; use native
`additive_compound_matrix(A, 0)` for the mathematically defined `1×1` zero.
Positive orders delegate to the native widened exact arithmetic, and an order
above the dimension returns `0×0`. Sparse inputs remain sparse unless
`sparse_output` requests otherwise.
"""
function AdditiveCompoundMatrix(
    matrix::AbstractMatrix, order; sparse_output::Bool=issparse(matrix)
)
    size(matrix, 1) == size(matrix, 2) || throw(
        DimensionMismatch(
            "AdditiveCompoundMatrix requires a square matrix; got size $(size(matrix))"
        ),
    )
    checked_order = _nonnegative_dimension(order, "R")
    checked_order == 0 && throw(
        ArgumentError(
            "pinned QETLAB AdditiveCompoundMatrix does not define R=0; " *
            "use native additive_compound_matrix(A, 0) for the 1×1 zero",
        ),
    )
    return additive_compound_matrix(matrix, checked_order; sparse_output=sparse_output)
end

function _compat_matrix_predicate_exact(matrix)
    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Number || return false
    real_type = typeof(real(zero(value_type)))
    return real_type <: Integer || real_type <: Rational
end

function _compat_matrix_predicate_tolerance(matrix, tolerance, default_tolerance)
    checked = if tolerance === nothing
        _compat_matrix_predicate_exact(matrix) ? 0 : default_tolerance()
    else
        _compat_finite_real(tolerance, "TOL")
    end
    checked >= 0 || throw(ArgumentError("TOL must be nonnegative"))
    isfinite(checked) || throw(ArgumentError("TOL must be finite; got $(repr(checked))"))
    return checked
end

function _compat_all_minor_default_tolerance(matrix)
    scale = norm(matrix)
    return max(size(matrix)...) * eps(scale)
end

"""
    IsPSD(X, TOL=nothing; allow_densify=false) -> MatrixPredicateResult

Compatibility spelling for QETLAB's positive-semidefinite predicate, with its
`eps(Float64)^(3/4)` default absolute tolerance for floating input. The result
is deliberately structured and three-valued: a numerical boundary is never
collapsed to QETLAB's Boolean. Unlike pinned `IsPSD.m`, this wrapper does not
replace `X` by its Hermitian part, and CVX objects are outside the
dependency-free numeric API. Exact integer and rational matrices use exact
zero tolerance.
"""
function IsPSD(
    matrix::AbstractMatrix{<:Number}, tolerance=nothing; allow_densify::Bool=false
)
    checked_tolerance = _compat_matrix_predicate_tolerance(
        matrix, tolerance, () -> eps(Float64)^(3 / 4)
    )
    return is_positive_semidefinite(
        matrix; atol=checked_tolerance, rtol=0, allow_densify=allow_densify
    )
end

"""
    IsLocallyPSD(X, K; atol=nothing, allow_densify=false,
                 max_submatrices=100_000) -> MatrixPredicateResult

Compatibility spelling for the pinned `K`-local PSD scan. Floating input uses
the absolute tolerance inherited from QETLAB `IsPSD` when `atol` is omitted;
exact input is checked exactly. The result preserves a principal-index witness
and an `unknown` boundary rather than returning a Boolean. No principal
submatrix is silently symmetrized.
"""
function IsLocallyPSD(
    matrix::AbstractMatrix{<:Number},
    order;
    atol=nothing,
    allow_densify::Bool=false,
    max_submatrices=100_000,
)
    checked_tolerance = _compat_matrix_predicate_tolerance(
        matrix, atol, () -> eps(Float64)^(3 / 4)
    )
    return is_locally_positive_semidefinite(
        matrix,
        order;
        atol=checked_tolerance,
        rtol=0,
        allow_densify=allow_densify,
        max_submatrices=max_submatrices,
    )
end

"""
    IsTotallyPositive(
        X, SUB_SIZES=nothing, TOL=nothing;
        allow_densify=false, max_minors=100_000
    ) -> MatrixPredicateResult

Compatibility spelling and positional argument order for QETLAB's all-minor
positive-determinant scan. For floating input, omitted `TOL` is the pinned
single absolute determinant tolerance
`max(size(X)...) * eps(norm(X))`; exact input uses exact comparisons.
Strictly singular minors violate total positivity, and nonzero determinants
inside the tolerance band are `unknown`. This intentionally differs from the
pinned Boolean routine, which accepts zero and small negative determinants.
"""
function IsTotallyPositive(
    matrix::AbstractMatrix{<:Number},
    sub_sizes=nothing,
    tolerance=nothing;
    allow_densify::Bool=false,
    max_minors=100_000,
)
    checked_tolerance = _compat_matrix_predicate_tolerance(
        matrix, tolerance, () -> _compat_all_minor_default_tolerance(matrix)
    )
    return is_totally_positive(
        matrix;
        orders=sub_sizes,
        atol=checked_tolerance,
        rtol=0,
        allow_densify=allow_densify,
        max_minors=max_minors,
    )
end

"""
    IsTotallyNonsingular(
        X, SUB_SIZES=nothing, TOL=nothing;
        allow_densify=false, max_minors=100_000
    ) -> MatrixPredicateResult

Compatibility spelling and positional argument order for QETLAB's all-minor
nonsingularity scan. The pinned default single determinant tolerance is used
for floating input and exact input is decided exactly. A represented singular
minor is a violation; a nonzero determinant inside tolerance is `unknown`
rather than an unsafe Boolean. The native combinatorial guard and explicit
sparse-densification policy remain in force.
"""
function IsTotallyNonsingular(
    matrix::AbstractMatrix{<:Number},
    sub_sizes=nothing,
    tolerance=nothing;
    allow_densify::Bool=false,
    max_minors=100_000,
)
    checked_tolerance = _compat_matrix_predicate_tolerance(
        matrix, tolerance, () -> _compat_all_minor_default_tolerance(matrix)
    )
    return is_totally_nonsingular(
        matrix;
        orders=sub_sizes,
        atol=checked_tolerance,
        rtol=0,
        allow_densify=allow_densify,
        max_minors=max_minors,
    )
end

function _compat_coherence_state(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}}
)
    if state isa AbstractMatrix && one(size(state, 1)) in size(state)
        Base.require_one_based_indexing(state)
        return vec(state)
    end
    return state
end

"""
    L1NormCoherence(RHO; allow_densify=false)

Return QETLAB's `l1`-norm of coherence in the computational basis. Pure row
and column vectors are accepted, while the Julia-native validation,
non-normalization, and sparse-densification policy remains in force.
"""
function L1NormCoherence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}};
    allow_densify::Bool=false,
)
    return l1_coherence(_compat_coherence_state(state); allow_densify=allow_densify)
end

"""
    RelEntCoherence(RHO; base=2, allow_densify=false)

Return the relative entropy of coherence. `base=2` preserves QETLAB's entropy
default; the keyword is exposed explicitly for migration code that needs a
different convention.
"""
function RelEntCoherence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}};
    base=2,
    allow_densify::Bool=false,
)
    return relative_entropy_coherence(
        _compat_coherence_state(state); base=base, allow_densify=allow_densify
    )
end

"""
    CoherenceRank(V, TOL=1e-10, BASIS=nothing; allow_densify=false)

Return the number of basis coefficients whose magnitudes exceed `TOL`.
This follows QETLAB's documented mathematical definition. The pinned
`CoherenceRank.m` implementation instead counts coefficients at or below the
tolerance; that reviewed upstream bug is intentionally not reproduced.
"""
function CoherenceRank(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    tolerance=1e-10,
    basis=nothing;
    allow_densify::Bool=false,
)
    vector = _compat_coherence_state(state)
    vector isa AbstractVector || throw(
        DimensionMismatch("V must be a row or column vector; got size $(size(state))")
    )
    checked_tolerance = _compat_finite_real(tolerance, "TOL")
    checked_tolerance >= 0 || throw(ArgumentError("TOL must be nonnegative"))
    return coherence_rank(
        vector; basis=basis, atol=checked_tolerance, rtol=0, allow_densify=allow_densify
    )
end

"""
    IsPPT(X, SYS=2, DIM=nothing, TOL=sqrt(eps(Float64));
          allow_densify=false) -> CriterionResult

Apply QETLAB's PPT test to a finite Hermitian matrix without requiring unit
trace. The compatibility result remains deliberately tri-state:
`CriterionEntanglementDetected`, `CriterionSatisfied`, or
`CriterionUnknown`. Values within `TOL` of the PSD boundary are never
collapsed to a Boolean.
"""
function IsPPT(
    input::AbstractMatrix{<:Number},
    systems=2,
    dim=nothing,
    tolerance=sqrt(eps(Float64));
    allow_densify::Bool=false,
)
    dimensions = _compat_bipartite_dimensions(input, dim, "IsPPT")
    checked_tolerance = _compat_finite_real(tolerance, "TOL")
    checked_tolerance >= 0 || throw(ArgumentError("TOL must be nonnegative"))
    all(isfinite, input) || throw(ArgumentError("X must contain only finite entries"))
    issparse(input) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "IsPPT requires a dense eigendecomposition; pass " *
                "allow_densify=true to permit converting this sparse matrix",
            ),
        )
    eltype(input) <: Union{Float32,Float64,ComplexF32,ComplexF64} || throw(
        ArgumentError("IsPPT requires a BLAS floating element type; got $(eltype(input))"),
    )
    dense = Matrix(input)
    hermiticity_residual = maximum(
        abs, dense - adjoint(dense); init=zero(typeof(real(zero(eltype(dense)))))
    )
    hermiticity_residual <= checked_tolerance || throw(
        ArgumentError(
            "X is not Hermitian within TOL=$checked_tolerance; maximum " *
            "residual is $hermiticity_residual",
        ),
    )
    work_input = (dense + adjoint(dense)) / 2
    transposed = partial_transpose(work_input, dimensions; systems=systems)
    decomposition = eigen(Hermitian((transposed + adjoint(transposed)) / 2))
    index = argmin(decomposition.values)
    value = decomposition.values[index]
    status = if value < -checked_tolerance
        CriterionEntanglementDetected
    elseif value > checked_tolerance
        CriterionSatisfied
    else
        CriterionUnknown
    end
    witness = status === CriterionSatisfied ? nothing : decomposition.vectors[:, index]
    message = if status === CriterionEntanglementDetected
        "PPT violation certifies entanglement"
    elseif status === CriterionSatisfied
        "PPT necessary condition holds with a tolerance margin; this is not " *
        "a separability certificate"
    else
        "PPT value lies within the numerical tolerance boundary"
    end
    return CriterionResult(
        :ppt, status, value, zero(value), checked_tolerance, witness, message
    )
end

end # module MATLABCompat
