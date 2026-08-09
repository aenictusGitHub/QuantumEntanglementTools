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
    AbstractOptimizationBackend,
    AbsPPTEnumerationConstraintLimit,
    AbsPPTEnumerationEarlyViolation,
    AffineScalar,
    BellScenario,
    CollinsGisinBehavior,
    ChoiRepresentation,
    ComplexAffineMatrix,
    CriterionEntanglementDetected,
    CriterionResult,
    CriterionSatisfied,
    CriterionUnknown,
    KrausRepresentation,
    HermitianAffineMatrix,
    NoOptimizationBackend,
    OperatorSpace,
    OperatorSumRepresentation,
    OptimizationLimits,
    SubsystemPermutationPlan,
    SuperoperatorRepresentation,
    _choi_hermiticity_result,
    _complete_positivity_result_from_choi,
    additive_compound_matrix,
    abs_ppt_constraints,
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
    canonical_map_decomposition,
    cb_norm,
    bcs_game_lower_bound,
    bcs_game_value,
    bell_inequality_bound,
    bell_inequality_qubit_bound,
    channel_distinguishability,
    complementary_channel,
    commutant,
    compound_matrix,
    concurrence,
    copositivity_criterion,
    clique_number_bounds,
    copositive_polynomial,
    dicke_state,
    diamond_norm,
    dual_channel,
    elementary_symmetric_polynomial,
    entanglement_of_formation,
    entangled_subspace,
    filter_normal_form,
    fidelity,
    fourier_matrix,
    generalized_gell_mann,
    generalized_pauli,
    gell_mann,
    ghz_state,
    gisin_state,
    horodecki_state,
    in_separable_ball,
    induced_matrix_norm,
    induced_schatten_lower_bound,
    input_size,
    is_completely_positive,
    is_absolutely_k_incoherent,
    is_block_positive,
    is_entangling_gate,
    is_locally_positive_semidefinite,
    is_hermiticity_preserving,
    is_positive_semidefinite,
    is_k_incoherent,
    is_separable,
    input_dimension,
    inverse_realign,
    isotropic_state,
    is_product_operator,
    is_product_vector,
    is_upb,
    is_abs_ppt,
    is_totally_nonsingular,
    is_totally_positive,
    ky_fan_norm,
    kraus_operators,
    kronecker_sum,
    l1_coherence,
    linear_to_basis,
    local_distinguishability,
    majorizes,
    matsumoto_fidelity,
    matsumoto_fidelity_model,
    maximum_output_fidelity,
    minimum_upb_size,
    upb,
    upb_sep_distinguishable,
    maximally_entangled,
    negativity,
    nonlocal_game_lower_bound,
    npa_membership,
    operator_schmidt_decomposition,
    operator_schmidt_rank,
    operator_sinkhorn,
    operator_sum_factors,
    operator_space,
    output_dimension,
    output_size,
    pauli,
    parallel_repetition,
    partial_trace,
    partial_transpose,
    partial_map,
    pauli_channel,
    polynomial_as_matrix,
    polynomial_bounds,
    polynomial_sos_bounds,
    positive_semidefinite_constraint,
    permutation_operator,
    permute_subsystems,
    random_density_matrix,
    random_graph,
    random_povm,
    random_ppt_state,
    random_probabilities,
    random_superoperator,
    random_state_vector,
    random_unitary,
    realign,
    reduction_map,
    renyi_entropy,
    relative_entropy_coherence,
    robustness_coherence,
    pure_k_coherence_robustness,
    schatten_norm,
    coherence_rank,
    schmidt_decomposition,
    schmidt_k_norm,
    schmidt_rank,
    sk_operator_norm,
    state_distinguishability,
    swap_operator,
    symmetric_extension,
    symmetric_inner_extension,
    symmetric_projector,
    symmetric_subspace_basis,
    tensor_power,
    tensor_product,
    tensor_sum,
    top_k_p_norm,
    top_k_p_norm_epigraph,
    top_k_p_norm_dual_epigraph,
    top_k_p_norm_dual,
    trace_distance_coherence,
    trace_norm,
    twirl,
    generalized_robustness_k_coherence,
    von_neumann_entropy,
    werner_state,
    w_state,
    xor_game_value
using Random: AbstractRNG
using LinearAlgebra: Hermitian, diag, eigen, ishermitian, norm, svdvals, tr
using SparseArrays: AbstractSparseVector, findnz, issparse, sparse, spdiagm, sparsevec

export Tensor,
    TensorSum,
    KroneckerSum,
    ParallelRepetition,
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
    RandomPPTState,
    RandomSuperoperator,
    ApplyMap,
    ChoiMatrix,
    KrausOperators,
    ComplementaryMap,
    DualMap,
    PartialMap,
    IsCP,
    IsHermPreserving,
    DepolarizingChannel,
    DephasingChannel,
    PauliChannel,
    ChoiMap,
    ReductionMap,
    TraceNorm,
    SchattenNorm,
    KyFanNorm,
    kpNorm,
    kpNormDual,
    SkOperatorNorm,
    IsBlockPositive,
    Distinguishability,
    LocalDistinguishability,
    IsSeparable,
    UPBSepDistinguishable,
    DiamondNorm,
    CBNorm,
    ChannelDistinguishability,
    MaximumOutputFidelity,
    InducedMatrixNorm,
    InducedSchattenNorm,
    Twirl,
    Purity,
    Entropy,
    Fidelity,
    MatsumotoFidelity,
    Negativity,
    SchmidtDecomposition,
    SkVectorNorm,
    SchmidtRank,
    OperatorSchmidtDecomposition,
    OperatorSchmidtRank,
    OperatorSinkhorn,
    FilterNormalForm,
    IsProductVector,
    IsProductOperator,
    IsEntanglingGate,
    IsUPB,
    MinUPBSize,
    UPB,
    Concurrence,
    EntFormation,
    EntangledSubspace,
    InSeparableBall,
    AbsPPTConstraints,
    IsAbsPPT,
    SymmetricExtension,
    SymmetricInnerExtension,
    IsPPT,
    L1NormCoherence,
    RelEntCoherence,
    CoherenceRank,
    RobkCohValue,
    IskIncoherent,
    IsAbskIncoh,
    RobustnessCoherence,
    TraceDistanceCoherence,
    GenRobustnesskCoherence,
    Majorizes,
    ElemSymPoly,
    CompoundMatrix,
    AdditiveCompoundMatrix,
    Commutant,
    CopositivePolynomial,
    PolynomialAsMatrix,
    PolynomialOptimize,
    PolynomialSOS,
    IsCopositive,
    CliqueNumber,
    NPAHierarchy,
    NonlocalGameLB,
    XORGameValue,
    BellInequalityMax,
    BellInequalityMaxQubits,
    BCSGameLB,
    BCSGameValue,
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
    ParallelRepetition(V, REPT;
                       max_entries=10_000_000,
                       max_work=100_000_000)

QETLAB-compatible entry point for parallel repetition of a nonlocal-game
coefficient tensor in full-probability notation. `REPT` must be a positive
integer. The native allocation and scalar-work guards remain explicit; pass
`nothing` for either guard only after independently bounding the requested
output.

The pinned MATLAB routine accidentally returns an unchanged input when
`REPT <= 0`. This wrapper deliberately rejects that invalid copy count and
delegates all valid work to [`parallel_repetition`](@ref).
"""
function ParallelRepetition(
    game::AbstractArray{<:Number,4},
    repetitions;
    max_entries=10_000_000,
    max_work=100_000_000,
)
    copy_count = _positive_dimension(repetitions, "REPT")
    return parallel_repetition(game, copy_count; max_entries=max_entries, max_work=max_work)
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

"""
    SymmetricProjection(
        DIM, P=2, PARTIAL=0, MODE=-1;
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

QETLAB-compatible symmetric projection or partial isometry. The Julia-only
resource keywords are forwarded to the native constructor; compatibility
output remains sparse, so `max_dense_entries` is validated but inactive.
"""
function SymmetricProjection(
    dim,
    copies=2,
    partial=0,
    mode=-1;
    max_columns=100_000,
    max_nonzeros=5_000_000,
    max_dense_entries=10_000_000,
    max_work=100_000_000,
)
    _projection_mode(mode)
    return if _flag(partial, "PARTIAL")
        symmetric_subspace_basis(
            dim,
            copies;
            sparse_output=true,
            max_columns=max_columns,
            max_nonzeros=max_nonzeros,
            max_dense_entries=max_dense_entries,
            max_work=max_work,
        )
    else
        symmetric_projector(
            dim,
            copies;
            sparse_output=true,
            max_columns=max_columns,
            max_nonzeros=max_nonzeros,
            max_dense_entries=max_dense_entries,
            max_work=max_work,
        )
    end
end

"""
    AntisymmetricProjection(
        DIM, P=2, PARTIAL=0, MODE=-1;
        max_columns=100_000, max_nonzeros=5_000_000,
        max_dense_entries=10_000_000, max_work=100_000_000,
    )

QETLAB-compatible antisymmetric projection or partial isometry, with the same
Julia-only resource keywords as [`SymmetricProjection`](@ref).
"""
function AntisymmetricProjection(
    dim,
    copies=2,
    partial=0,
    mode=-1;
    max_columns=100_000,
    max_nonzeros=5_000_000,
    max_dense_entries=10_000_000,
    max_work=100_000_000,
)
    _projection_mode(mode)
    return if _flag(partial, "PARTIAL")
        antisymmetric_subspace_basis(
            dim,
            copies;
            sparse_output=true,
            max_columns=max_columns,
            max_nonzeros=max_nonzeros,
            max_dense_entries=max_dense_entries,
            max_work=max_work,
        )
    else
        antisymmetric_projector(
            dim,
            copies;
            sparse_output=true,
            max_columns=max_columns,
            max_nonzeros=max_nonzeros,
            max_dense_entries=max_dense_entries,
            max_work=max_work,
        )
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

"""
    DickeState(
        P, E=1, NRML=1;
        max_nonzeros=1_000_000,
        max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

QETLAB-compatible sparse Dicke-state constructor. The Julia-only resource
keywords retain the native preflight defaults; `max_dense_entries` is
validated but inactive for compatibility's sparse output.
"""
function DickeState(
    parties,
    excitations=1,
    normalized=1;
    max_nonzeros=1_000_000,
    max_dense_entries=10_000_000,
    max_work=100_000_000,
)
    return dicke_state(
        parties,
        excitations;
        normalized=_flag(normalized, "NRML"),
        sparse_output=true,
        max_nonzeros=max_nonzeros,
        max_dense_entries=max_dense_entries,
        max_work=max_work,
    )
end

"""QETLAB-compatible isotropic-state constructor."""
IsotropicState(dim, alpha) = isotropic_state(dim, alpha; sparse_output=true)

"""
    WernerState(DIM, ALPHA; kwargs...)

QETLAB-compatible Werner-state constructor. Scalar `ALPHA` selects the
bipartite family. A vector with `p! - 1` entries selects the multipartite
lexicographic-permutation family and forwards the native physicality and
resource checks. The multipartite implementation intentionally corrects the
pinned loop-overwrite defect, which retained only its last parameter.
"""
function WernerState(dim, alpha; kwargs...)
    return werner_state(dim, alpha; sparse_output=true, kwargs...)
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
    EntangledSubspace(DIM, LOCALDIM, R=1;
                       max_nonzeros=1_000_000,
                       max_work=5_000_000)

Return QETLAB's sparse diagonal-Vandermonde basis for an `R`-entangled
bipartite subspace. A scalar `LOCALDIM` selects equal local dimensions.
Compatibility output uses `Float64`, while the native
[`entangled_subspace`](@ref) API additionally permits an explicit coefficient
type. Positive dimensions, `0 <= R < min(LOCALDIM)`, and the sharp maximal
subspace dimension are validated before construction; no columns are
normalized.
"""
function EntangledSubspace(
    subspace_dimension, local_dims, r=1; max_nonzeros=1_000_000, max_work=5_000_000
)
    return entangled_subspace(
        subspace_dimension,
        local_dims;
        r=r,
        coefficient_type=Float64,
        max_nonzeros=max_nonzeros,
        max_work=max_work,
    )
end

"""
    IsEntanglingGate(U, DIM=nothing; kwargs...)

Compatibility spelling for [`is_entangling_gate`](@ref). `DIM` keeps QETLAB's
scalar, vector, or two-row output/input layout forms. The result is structured
rather than Boolean: inspect `status`, the local-factor/permutation
certificate for `:not_entangling`, or the sparse product witness for
`:entangling`. Numerical boundaries remain `:unknown`.
"""
function IsEntanglingGate(gate, dims=nothing; kwargs...)
    return is_entangling_gate(gate, dims; kwargs...)
end

"""
    IsUPB(U, V, ...;
          return_witness=false,
          structured=false,
          normalization=:allow,
          kwargs...)

Compatibility spelling for QETLAB's local-factor `IsUPB` entry point. The
wrapper explicitly accepts arbitrary nonzero column scaling by default and
delegates every mathematical decision to [`is_upb`](@ref).

The one-output form returns a `Bool` for conclusive results. Set
`return_witness=true` to receive `(boolean, witness_factors)`; the witness is
`nothing` when no extension witness applies. Set `structured=true` to receive
the complete native [`QuantumEntanglementTools.UPBAnalysisResult`](@ref)
instead. A native `:unknown` result raises `DomainError` in the Boolean forms
rather than becoming a false mathematical answer.

Unlike the pinned routine, the wrapper checks mutual orthogonality and
incompleteness and uses Hermitian orthogonality for complex witnesses.
"""
function IsUPB(
    first_factor,
    second_factor,
    remaining_factors...;
    return_witness=false,
    structured=false,
    normalization=:allow,
    kwargs...,
)
    include_witness = _flag(return_witness, "return_witness")
    return_structured = _flag(structured, "structured")
    include_witness &&
        return_structured &&
        throw(ArgumentError("return_witness and structured cannot both be true"))
    result = is_upb(
        first_factor, second_factor, remaining_factors...; normalization, kwargs...
    )
    return_structured && return result
    result.status === :unknown && throw(
        DomainError(
            result,
            "IsUPB reached a numerical or resource boundary; inspect the " *
            "native structured result",
        ),
    )
    return include_witness ? (result.is_upb, result.witness_factors) : result.is_upb
end

"""
    MinUPBSize(DIM, VERBOSE=1)

Compatibility spelling for QETLAB's theorem-table lookup. A known case returns
the exact integer and, when `VERBOSE` is true, prints its recorded primary
reference. An unresolved case raises `DomainError` carrying the native
[`QuantumEntanglementTools.MinimumUPBSizeResult`](@ref); use
[`QuantumEntanglementTools.minimum_upb_size`](@ref) to receive
`:unknown` without an exception.
"""
function MinUPBSize(dims, verbose=1)
    emit_reference = _flag(verbose, "VERBOSE")
    result = minimum_upb_size(dims)
    result.status === :known || throw(
        DomainError(
            result,
            "the reviewed pinned theorem table does not determine this exact " *
            "UPB minimum; inspect the native structured result",
        ),
    )
    emit_reference && println(result.reference)
    return result.size
end

function _upb_compat_output(construction, output)
    output isa Symbol ||
        throw(ArgumentError("output must be :global, :local, or :structured"))
    output === :global && return copy(construction.global_vectors)
    output === :local && return map(copy, construction.local_factors)
    output === :structured && return construction
    return throw(ArgumentError("output must be :global, :local, or :structured"))
end

function _upb_compat_dimensions(input)
    input isa Integer && return (input,)
    input isa Tuple ||
        input isa AbstractVector ||
        throw(ArgumentError("numeric UPB dimensions must be a scalar, tuple, or vector"))
    return input
end

function _upb_compat_construct(rng, input, arguments; output, kwargs)
    if input isa Integer || input isa Tuple || input isa AbstractVector
        length(arguments) <= 1 ||
            throw(ArgumentError("UPB(DIM, VERBOSE) accepts at most one positional flag"))
        emit_reference = _flag(isempty(arguments) ? 1 : arguments[1], "VERBOSE")
        dimensions = _upb_compat_dimensions(input)
        construction =
            rng === nothing ? upb(dimensions; kwargs...) : upb(rng, dimensions; kwargs...)
        emit_reference && println(construction.reference)
        return _upb_compat_output(construction, output)
    end

    construction = if rng === nothing
        upb(input, arguments...; kwargs...)
    else
        upb(rng, input, arguments...; kwargs...)
    end
    return _upb_compat_output(construction, output)
end

"""
    UPB(NAME, family_arguments...; output=:global, kwargs...)
    UPB(DIM, VERBOSE=1; output=:global, kwargs...)
    UPB(rng::AbstractRNG, NAME_OR_DIM, arguments...; output=:global, kwargs...)

Compatibility spelling for QETLAB's UPB catalog. `output=:global` returns the
matrix whose columns are global product states. Because Julia has no MATLAB
`nargout`, use `output=:local` for a tuple containing one local-factor matrix
per party or `output=:structured` for the complete
[`QuantumEntanglementTools.UPBConstruction`](@ref).

A numeric `DIM` may be a scalar, tuple, or vector; the optional positional
`VERBOSE` flag is validated and prints the selected primary reference when
true. Named-family parameters retain their pinned positional order.
Randomized constructions require the leading explicit `rng`; the no-RNG form
raises instead of consuming Julia's global random stream. Entry, work,
full-spark-minor, retry, precision, and tolerance keywords are forwarded to
[`QuantumEntanglementTools.upb`](@ref).

The compatibility surface deliberately retains the native corrections for the
extendible pinned `GenTiles1(2)` result and the nonorthogonal pinned
`John2^4k` reshape-order branch.
"""
function UPB(input, arguments...; output=:global, kwargs...)
    return _upb_compat_construct(nothing, input, arguments; output, kwargs)
end

function UPB(rng::AbstractRNG, input, arguments...; output=:global, kwargs...)
    return _upb_compat_construct(rng, input, arguments; output, kwargs)
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

function _random_superoperator_compat_dimensions(dim)
    if dim isa Integer
        dimension = _positive_dimension(dim, "DIM")
        return dimension, dimension
    elseif dim isa Tuple || dim isa AbstractVector
        dim isa AbstractVector && Base.require_one_based_indexing(dim)
        length(dim) == 2 ||
            throw(DimensionMismatch("DIM must contain exactly two dimensions"))
        return (
            _positive_dimension(dim[1], "DIM[1]"), _positive_dimension(dim[2], "DIM[2]")
        )
    end
    return throw(
        ArgumentError("DIM must be a positive integer or a two-entry tuple/vector")
    )
end

"""
    RandomSuperoperator(
        rng,
        DIM,
        TP=1,
        UN=0,
        RE=0,
        KR=prod(DIM);
        diagnostics=false,
        allow_proportional_unital=false,
        T=Float64,
        atol=0,
        rtol=nothing,
        max_attempts=8,
        max_iterations=1000,
        max_condition_number=nothing,
        max_dimension=4096,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    )

QETLAB-compatible random completely positive map with a mandatory explicit
RNG. The default return is the raw Choi matrix expected by QETLAB callers.
Set `diagnostics=true` to receive the native
[`QuantumEntanglementTools.RandomSuperoperatorResult`](@ref), including
bounded-convergence and marginal residual evidence.

Strict trace preservation plus unitality is impossible for unequal input and
output dimensions. The pinned routine labels its proportional-output branch
as unital; this wrapper rejects that request unless
`allow_proportional_unital=true`, in which case the result explicitly records
`Φ(I_in) ≈ (d_in/d_out)I_out` and does not claim unitality. A failed bounded
construction raises `DomainError` in raw-output mode rather than returning an
uncertified Choi matrix.
"""
function RandomSuperoperator(
    rng::AbstractRNG,
    dim,
    trace_preserving=1,
    unital=0,
    real_output=0,
    kraus_rank=nothing;
    diagnostics=false,
    allow_proportional_unital=false,
    T=Float64,
    atol=0,
    rtol=nothing,
    max_attempts=8,
    max_iterations=1_000,
    max_condition_number=nothing,
    max_dimension=4_096,
    max_entries=10_000_000,
    max_work=1_000_000_000,
)
    input_dimension, output_dimension = _random_superoperator_compat_dimensions(dim)
    requested_trace_preserving = _flag(trace_preserving, "TP")
    requested_unital = _flag(unital, "UN")
    requested_real = _flag(real_output, "RE")
    diagnostics isa Bool ||
        throw(ArgumentError("diagnostics must be Bool; got $(repr(diagnostics))"))
    allow_proportional_unital isa Bool || throw(
        ArgumentError(
            "allow_proportional_unital must be Bool; got " *
            repr(allow_proportional_unital),
        ),
    )

    unequal_balancing =
        requested_trace_preserving &&
        requested_unital &&
        input_dimension != output_dimension
    allow_proportional_unital &&
        !unequal_balancing &&
        throw(
            ArgumentError(
                "allow_proportional_unital=true applies only to an unequal-dimensional " *
                "TP=1, UN=1 request",
            ),
        )
    unequal_balancing &&
        !allow_proportional_unital &&
        throw(
            ArgumentError(
                "TP=1 and UN=1 are incompatible for unequal dimensions; set " *
                "allow_proportional_unital=true only to request the pinned routine's " *
                "corrected proportional-output branch",
            ),
        )

    result = random_superoperator(
        rng,
        (input_dimension, output_dimension);
        trace_preserving=requested_trace_preserving,
        unital=requested_unital && !unequal_balancing,
        proportional_unital=unequal_balancing,
        real=requested_real,
        kraus_rank=kraus_rank,
        representation=:choi,
        T,
        atol,
        rtol,
        max_attempts,
        max_iterations,
        max_condition_number,
        max_dimension,
        max_entries,
        max_work,
    )
    diagnostics && return result
    result.succeeded || throw(
        DomainError(
            result,
            "the bounded random-superoperator construction did not establish its " *
            "requested marginal guarantees; request diagnostics=true to inspect it",
        ),
    )
    return choi_matrix(result.representation)
end

function _compat_map_dimensions(dim)
    if dim isa Integer
        value = _positive_dimension(dim, "DIM")
        return OperatorSpace(value, value, value, value)
    elseif dim isa AbstractMatrix
        Base.require_one_based_indexing(dim)
        size(dim) == (2, 2) || throw(
            DimensionMismatch("a map DIM matrix must have size (2, 2); got $(size(dim))"),
        )
        return OperatorSpace(
            _positive_dimension(dim[1, 1], "DIM[1,1]"),
            _positive_dimension(dim[2, 1], "DIM[2,1]"),
            _positive_dimension(dim[1, 2], "DIM[1,2]"),
            _positive_dimension(dim[2, 2], "DIM[2,2]"),
        )
    elseif dim isa Tuple || dim isa AbstractVector
        dim isa AbstractVector && Base.require_one_based_indexing(dim)
        length(dim) == 2 || throw(
            DimensionMismatch(
                "map DIM must contain input and output dimensions; got " *
                "$(length(dim)) entries",
            ),
        )
        input_dimension = _positive_dimension(dim[1], "DIM[1]")
        output_dimension = _positive_dimension(dim[2], "DIM[2]")
        return OperatorSpace(
            input_dimension, input_dimension, output_dimension, output_dimension
        )
    end
    return throw(
        ArgumentError(
            "DIM must be an integer, a two-entry tuple/vector, or a 2-by-2 matrix"
        ),
    )
end

function _validate_compat_map_dimensions(map, dim)
    dim === nothing && return map
    expected = _compat_map_dimensions(dim)
    actual = operator_space(map)
    actual == expected || throw(
        DimensionMismatch(
            "DIM describes operator space $(input_size(expected)) → " *
            "$(output_size(expected)), but PHI describes $(input_size(actual)) → " *
            "$(output_size(actual))",
        ),
    )
    return map
end

function _compat_input_size(input_hint)
    if input_hint isa Integer
        dimension = _positive_dimension(input_hint, "input dimension")
        return (dimension, dimension)
    elseif input_hint isa Tuple || input_hint isa AbstractVector
        input_hint isa AbstractVector && Base.require_one_based_indexing(input_hint)
        length(input_hint) == 2 ||
            throw(DimensionMismatch("input size must contain exactly two dimensions"))
        return (
            _positive_dimension(input_hint[1], "input row dimension"),
            _positive_dimension(input_hint[2], "input column dimension"),
        )
    end
    return throw(ArgumentError("input size must be an integer or a two-entry shape"))
end

function _compat_choi_space(matrix, dim, input_hint)
    size(matrix, 1) > 0 && size(matrix, 2) > 0 ||
        throw(ArgumentError("a Choi matrix must have nonzero dimensions"))
    space = if dim !== nothing
        _compat_map_dimensions(dim)
    elseif input_hint !== nothing
        input_rows, input_columns = _compat_input_size(input_hint)
        rem(size(matrix, 1), input_rows) == 0 || throw(
            DimensionMismatch(
                "Choi row count $(size(matrix, 1)) is not divisible by inferred " *
                "input row dimension $input_rows",
            ),
        )
        rem(size(matrix, 2), input_columns) == 0 || throw(
            DimensionMismatch(
                "Choi column count $(size(matrix, 2)) is not divisible by inferred " *
                "input column dimension $input_columns",
            ),
        )
        OperatorSpace(
            input_rows,
            input_columns,
            div(size(matrix, 1), input_rows),
            div(size(matrix, 2), input_columns),
        )
    else
        row_dimension = isqrt(size(matrix, 1))
        column_dimension = isqrt(size(matrix, 2))
        row_dimension > 0 &&
        column_dimension > 0 &&
        row_dimension^2 == size(matrix, 1) &&
        column_dimension^2 == size(matrix, 2) || throw(
            ArgumentError(
                "cannot infer operator-space dimensions from Choi size " *
                "$(size(matrix)); provide DIM=[input_rows output_rows; " *
                "input_columns output_columns]",
            ),
        )
        OperatorSpace(row_dimension, column_dimension, row_dimension, column_dimension)
    end
    expected = (
        space.input_rows * space.output_rows, space.input_columns * space.output_columns
    )
    size(matrix) == expected || throw(
        DimensionMismatch(
            "DIM describes Choi size $expected, inconsistent with actual size " *
            "$(size(matrix))",
        ),
    )
    return space
end

function _is_cp_kraus_collection(value)
    (value isa Tuple || value isa AbstractVector) || return false
    isempty(value) && return false
    return all(operator -> operator isa AbstractMatrix && eltype(operator) <: Number, value)
end

function _is_factor_cell_matrix(value)
    value isa AbstractMatrix || return false
    eltype(value) <: Number && return false
    isempty(value) && return false
    return all(operator -> operator isa AbstractMatrix, value)
end

function _compat_factor_map(phi; allow_row_cp::Bool=true)
    if _is_cp_kraus_collection(phi)
        return KrausRepresentation(phi)
    end
    _is_factor_cell_matrix(phi) || throw(
        ArgumentError(
            "a QETLAB-style factor cell must be a nonempty matrix whose entries " *
            "are numeric matrices",
        ),
    )
    Base.require_one_based_indexing(phi)
    factor_rows, factor_columns = size(phi)
    if factor_columns == 1
        return KrausRepresentation([phi[index, 1] for index in 1:factor_rows])
    elseif factor_columns == 2
        return OperatorSumRepresentation(
            [phi[index, 1] for index in 1:factor_rows],
            [phi[index, 2] for index in 1:factor_rows],
        )
    elseif factor_rows == 1 && factor_columns > 2 && allow_row_cp
        return KrausRepresentation([phi[1, index] for index in 1:factor_columns])
    end
    return throw(
        DimensionMismatch(
            "PHI factor cells must have one or two columns" *
            (allow_row_cp ? ", or be a one-row CP collection" : ""),
        ),
    )
end

function _compat_map(phi, dim=nothing; input_hint=nothing, allow_row_cp::Bool=true)
    if phi isa AbstractMapRepresentation
        return _validate_compat_map_dimensions(phi, dim)
    elseif phi isa AbstractMatrix && eltype(phi) <: Number
        space = _compat_choi_space(phi, dim, input_hint)
        return ChoiRepresentation(phi, space)
    elseif _is_cp_kraus_collection(phi) || _is_factor_cell_matrix(phi)
        return _validate_compat_map_dimensions(
            _compat_factor_map(phi; allow_row_cp=allow_row_cp), dim
        )
    end
    return throw(
        ArgumentError(
            "PHI must be a numeric Choi matrix, a nonempty CP vector/tuple, a " *
            "QETLAB-style one/two-column matrix of factors, or an " *
            "AbstractMapRepresentation",
        ),
    )
end

"""
    ApplyMap(X, PHI)

QETLAB-compatible application of a map supplied as a numeric Choi matrix,
a vector/tuple of completely-positive Kraus matrices, or a matrix of factor
cells. One-column factor cells are CP; two-column cells represent
`sum(Aᵢ * X * Bᵢ')`. A one-row collection with more than two entries follows
QETLAB's CP branch. Rectangular input/output matrix spaces are inferred from
`size(X)` and the Choi or factor dimensions.
"""
function ApplyMap(input::AbstractMatrix, phi)
    map = _compat_map(phi; input_hint=size(input))
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
    map = _compat_map(phi; allow_row_cp=false)
    matrix = choi_matrix(map)
    system == 2 && return matrix
    input_rows, input_columns = input_size(map)
    output_rows, output_columns = output_size(map)
    row_plan = SubsystemPermutationPlan((input_rows, output_rows), (2, 1))
    column_plan = SubsystemPermutationPlan((input_columns, output_columns), (2, 1))
    return permute_subsystems(matrix, row_plan, column_plan)
end

"""
    KrausOperators(PHI, DIM=nothing; atol=0, rtol=sqrt(eps(Float64)),
                   allow_densify=false, diagnostics=false)

Compute the pinned canonical factor branches. A completely-positive map
returns a vector of canonical Kraus matrices. A Hermiticity-preserving map
returns a two-column factor matrix with equal positive pairs first and signed
negative pairs second. Every other map returns a two-column SVD factorization
of the unmodified Choi matrix.

Set `diagnostics=true` to return the native
`CanonicalMapDecompositionResult`; otherwise the raw vector/two-column shape
follows QETLAB. A numerical complete-positivity boundary uses the paired
branch and is never silently projected onto a completely-positive map.
Sparse Choi spectral work requires `allow_densify=true`. Full rectangular
`DIM=[input_rows output_rows; input_columns output_columns]` is supported.
"""
function _compat_operator_sum_cells(map::OperatorSumRepresentation)
    factors = operator_sum_factors(map)
    cells = Matrix{AbstractMatrix}(undef, length(factors.left), 2)
    for index in eachindex(factors.left, factors.right)
        cells[index, 1] = factors.left[index]
        cells[index, 2] = factors.right[index]
    end
    return cells
end

function KrausOperators(
    phi,
    dim=nothing;
    atol=0,
    rtol=sqrt(eps(Float64)),
    allow_densify::Bool=false,
    diagnostics::Bool=false,
)
    map = _compat_map(phi, dim; allow_row_cp=false)
    result = canonical_map_decomposition(
        map; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    diagnostics && return result
    factors = operator_sum_factors(result.representation)
    if result.classification === :completely_positive
        return factors.left
    end
    return _compat_operator_sum_cells(result.representation)
end

"""
    ComplementaryMap(PHI, DIM=nothing; atol=0, rtol=sqrt(eps(Float64)),
                     allow_densify=false)

Return the complementary map while preserving a supplied raw dilation.
One-column Kraus data produce one-column complementary factors; two-column
data produce paired factors. A paired complement is defined when the original
output operator space is square, including rectangular input operator spaces,
and is rejected when output row and column dimensions differ.

Numeric Choi input is canonically factorized and returns a numeric Choi
matrix. That spectral path requires explicit sparse densification. Project
representations retain their representation kind.
"""
function ComplementaryMap(
    phi, dim=nothing; atol=0, rtol=sqrt(eps(Float64)), allow_densify::Bool=false
)
    map = _compat_map(phi, dim; allow_row_cp=false)
    complement = complementary_channel(
        map; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    if phi isa AbstractMapRepresentation
        return complement
    elseif _is_cp_kraus_collection(phi)
        return kraus_operators(complement)
    elseif _is_factor_cell_matrix(phi)
        if complement isa KrausRepresentation
            operators = kraus_operators(complement)
            cells = Matrix{AbstractMatrix}(undef, length(operators), 1)
            for index in eachindex(operators)
                cells[index, 1] = operators[index]
            end
            return cells
        end
        return _compat_operator_sum_cells(complement)
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
    if _is_cp_kraus_collection(phi)
        return [copy(adjoint(operator)) for operator in phi]
    elseif _is_factor_cell_matrix(phi)
        return Base.map(operator -> copy(adjoint(operator)), phi)
    end
    map = _compat_map(phi, dim)
    dual = dual_channel(map)
    if phi isa AbstractMapRepresentation
        return dual
    end
    return choi_matrix(dual)
end

function _compat_dimension_product(dimensions, name::AbstractString)
    result = 1
    for dimension in dimensions
        result = try
            Base.checked_mul(result, dimension)
        catch error
            error isa OverflowError || rethrow()
            throw(ArgumentError("$name product exceeds typemax(Int)"))
        end
    end
    return result
end

function _compat_operator_dimensions(input::AbstractMatrix, dim)
    Base.require_one_based_indexing(input)
    row_count, column_count = size(input)
    row_dimensions, column_dimensions = if dim === nothing
        row_dimension = isqrt(row_count)
        column_dimension = isqrt(column_count)
        row_dimension^2 == row_count && column_dimension^2 == column_count || throw(
            ArgumentError(
                "cannot infer two equal row and column subsystem dimensions " *
                "from input size $(size(input)); provide DIM",
            ),
        )
        ((row_dimension, row_dimension), (column_dimension, column_dimension))
    elseif dim isa Integer
        row_count == column_count || throw(
            DimensionMismatch(
                "scalar DIM requires a square input matrix; got $(size(input))"
            ),
        )
        dimensions = _expand_scalar_dimension(dim, row_count, "DIM")
        (dimensions, dimensions)
    elseif dim isa AbstractMatrix
        Base.require_one_based_indexing(dim)
        size(dim, 1) == 2 || throw(
            DimensionMismatch(
                "a PartialMap DIM matrix must have two rows; got $(size(dim))"
            ),
        )
        size(dim, 2) > 0 || throw(
            ArgumentError("a PartialMap DIM matrix must have at least one column")
        )
        (
            Tuple(
                _positive_dimension(dim[1, index], "DIM[1,$index]") for
                index in axes(dim, 2)
            ),
            Tuple(
                _positive_dimension(dim[2, index], "DIM[2,$index]") for
                index in axes(dim, 2)
            ),
        )
    else
        row_count == column_count || throw(
            DimensionMismatch(
                "vector DIM requires a square input matrix; use a two-row DIM " *
                "matrix for rectangular input",
            ),
        )
        dimensions = _dimension_tuple(dim)
        (dimensions, dimensions)
    end
    length(row_dimensions) == length(column_dimensions) || throw(
        DimensionMismatch(
            "row and column DIM layouts must describe the same number of subsystems"
        ),
    )
    row_product = _compat_dimension_product(row_dimensions, "row DIM")
    column_product = _compat_dimension_product(column_dimensions, "column DIM")
    row_product == row_count || throw(
        DimensionMismatch(
            "row DIM product $row_product does not match input row count $row_count"
        ),
    )
    column_product == column_count || throw(
        DimensionMismatch(
            "column DIM product $column_product does not match input column count $column_count",
        ),
    )
    return row_dimensions, column_dimensions
end

"""
    PartialMap(X, PHI, SYS=2, DIM=nothing)

Apply `PHI` to subsystem `SYS` using QETLAB argument order. A vector `DIM`
describes common row and column subsystem dimensions. A two-row `DIM` matrix
describes independent row and column dimensions, so rectangular multipartite
operators and rectangular local maps are supported without constructing a
global Kronecker superoperator.
"""
function PartialMap(input::AbstractMatrix, phi, system=2, dim=nothing)
    system isa Integer && !(system isa Bool) ||
        throw(ArgumentError("SYS must be an integer; got $(repr(system))"))
    row_dimensions, column_dimensions = _compat_operator_dimensions(input, dim)
    1 <= system <= length(row_dimensions) || throw(
        ArgumentError("SYS must be between 1 and $(length(row_dimensions)); got $system"),
    )
    map = _compat_map(phi; input_hint=(row_dimensions[system], column_dimensions[system]))
    return partial_map(input, map, system, row_dimensions, column_dimensions)
end

function _compat_choi_hermiticity_result(
    matrix::AbstractMatrix, tolerance; nonsquare_status::Symbol
)
    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Number || throw(
        ArgumentError(
            "IsHermPreserving requires a concrete numeric Choi element type; got $value_type",
        ),
    )
    real_type = typeof(real(zero(value_type)))
    exact = real_type <: Integer || real_type <: Rational
    absolute_tolerance = exact ? zero(real_type) : tolerance
    return _choi_hermiticity_result(
        matrix;
        atol=absolute_tolerance,
        rtol=zero(absolute_tolerance),
        nonsquare_status=nonsquare_status,
    )
end

function _compat_predicate_tolerance(value, name::AbstractString)
    value isa Bool &&
        throw(ArgumentError("$name must be a finite nonnegative real number, not Bool"))
    checked = _compat_finite_real(value, name)
    checked >= zero(checked) ||
        throw(ArgumentError("$name must be nonnegative; got $(repr(value))"))
    return checked
end

function _compat_choi_complete_positivity_result(
    matrix::AbstractMatrix, tolerance; allow_densify::Bool
)
    value_type = eltype(matrix)
    isconcretetype(value_type) && value_type <: Number || throw(
        ArgumentError(
            "IsCP requires a concrete numeric Choi element type; got $value_type"
        ),
    )
    real_type = typeof(real(zero(value_type)))
    exact = real_type <: Integer || real_type <: Rational
    absolute_tolerance = exact ? zero(real_type) : tolerance
    return _complete_positivity_result_from_choi(
        matrix;
        atol=absolute_tolerance,
        rtol=zero(absolute_tolerance),
        allow_densify=allow_densify,
        nonsquare_status=:violated,
    )
end

"""
    IsCP(PHI, TOL=eps(Float64)^(3/4); allow_densify=false)

QETLAB-compatible complete-positivity diagnostic with a structured
`MatrixPredicateResult`. One-column Kraus input is satisfied by construction.
Other raw factor data are converted to the unmodified Choi matrix and tested
for Hermiticity and positive semidefiniteness. A robust violation is
`MatrixPredicateViolated`; a nonzero tolerance-boundary defect is
`MatrixPredicateUnknown`, never a Boolean negative or a repaired Choi matrix.

Sparse spectral work requires `allow_densify=true`. Exact integer and rational
input is decided exactly, independent of the floating QETLAB default
tolerance.
"""
function IsCP(phi, tolerance=eps(Float64)^(3 / 4); allow_densify::Bool=false)
    checked_tolerance = _compat_predicate_tolerance(tolerance, "TOL")
    if phi isa AbstractMapRepresentation
        matrix = choi_matrix(phi)
        real_type = typeof(real(zero(eltype(matrix))))
        exact = real_type <: Integer || real_type <: Rational
        absolute_tolerance = exact ? zero(real_type) : checked_tolerance
        return is_completely_positive(
            phi;
            atol=absolute_tolerance,
            rtol=zero(absolute_tolerance),
            allow_densify=allow_densify,
        )
    elseif phi isa AbstractMatrix && eltype(phi) <: Number
        return _compat_choi_complete_positivity_result(
            phi, checked_tolerance; allow_densify=allow_densify
        )
    end

    map = _compat_map(phi; allow_row_cp=false)
    map isa KrausRepresentation && return is_completely_positive(map; atol=0, rtol=0)
    return _compat_choi_complete_positivity_result(
        choi_matrix(map), checked_tolerance; allow_densify=allow_densify
    )
end

"""
    IsHermPreserving(PHI, TOL=eps(Float64)^(3/4))

QETLAB-compatible Choi-Hermiticity diagnostic. The return value is a
`MatrixPredicateResult`: exact Hermiticity is `MatrixPredicateSatisfied`, a
defect larger than `TOL` is `MatrixPredicateViolated`, and a nonzero defect
inside the tolerance boundary is deliberately `MatrixPredicateUnknown`
instead of being reported as a mathematical proof. Numeric nonsquare Choi
matrices reproduce QETLAB's negative result.
"""
function IsHermPreserving(phi, tolerance=eps(Float64)^(3 / 4))
    checked_tolerance = _compat_predicate_tolerance(tolerance, "TOL")
    if phi isa AbstractMapRepresentation
        matrix = choi_matrix(phi)
        real_type = typeof(real(zero(eltype(matrix))))
        exact = real_type <: Integer || real_type <: Rational
        if exact
            return is_hermiticity_preserving(
                phi; atol=zero(real_type), rtol=zero(real_type)
            )
        end
        return is_hermiticity_preserving(
            phi; atol=checked_tolerance, rtol=zero(checked_tolerance)
        )
    elseif phi isa AbstractMatrix && eltype(phi) <: Number
        return _compat_choi_hermiticity_result(
            phi, checked_tolerance; nonsquare_status=:violated
        )
    end
    map = _compat_map(phi; allow_row_cp=false)
    return _compat_choi_hermiticity_result(
        choi_matrix(map), checked_tolerance; nonsquare_status=:violated
    )
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
    Twirl(
        X, TYPE="werner", P=2;
        sparse_output=true,
        allow_densify=false,
        max_basis_size=256,
        max_nonzeros=5_000_000,
        max_dense_entries=1_000_000,
        max_work=1_000_000_000,
    )

QETLAB-compatible positional wrapper for [`twirl`](@ref). `TYPE` is a
case-insensitive string or symbol naming `"werner"`, `"isotropic"`, `"real"`,
or `"pauli"`. The wrapper preserves QETLAB's sparse result storage by default,
while forwarding the native allocation and work guards.

Unlike the pinned routine, isotropic and Pauli twirls require exactly `P=2`;
dimension roots are exact, and malformed inputs are rejected before any
spanning-family construction.
"""
function Twirl(
    input,
    type="werner",
    copies=2;
    sparse_output::Bool=true,
    allow_densify::Bool=false,
    max_basis_size=256,
    max_nonzeros=5_000_000,
    max_dense_entries=1_000_000,
    max_work=1_000_000_000,
)
    (type isa AbstractString || type isa Symbol) || throw(
        ArgumentError(
            "TYPE must be one of \"werner\", \"isotropic\", \"real\", or \"pauli\""
        ),
    )
    kind = Symbol(lowercase(String(type)))
    return twirl(
        input;
        kind,
        copies,
        sparse_output,
        allow_densify,
        max_basis_size,
        max_nonzeros,
        max_dense_entries,
        max_work,
    )
end

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

function _compat_vectorize_kp_input(matrix::AbstractMatrix{<:Number})
    Base.require_one_based_indexing(matrix)
    if issparse(matrix)
        rows, columns, values = findnz(sparse(matrix))
        linear_indices = rows .+ (columns .- 1) .* size(matrix, 1)
        return sparsevec(linear_indices, values, length(matrix))
    end
    return vec(matrix)
end

function _compat_kp_input(input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}})
    if input isa AbstractMatrix && min(size(input)...) == 1
        return _compat_vectorize_kp_input(input)
    end
    return input
end

"""
    kpNorm(X, K, P; allow_densify=false)

Numeric-array compatibility entry point for QETLAB's `(K,P)` norm. Julia
vectors and one-row or one-column matrices use the magnitudes of their `K`
largest entries; other matrices use their `K` largest singular values.
`K` is clipped to the available spectrum. The pinned QETLAB numeric vector
code sorts signed entries despite documenting a magnitude norm; this wrapper
uses magnitudes and therefore deliberately corrects negative-vector results.
CVX/model expressions are outside the dependency-free compatibility surface.
"""
function kpNorm(
    input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k,
    p;
    allow_densify::Bool=false,
)
    compatible_input = _compat_kp_input(input)
    return top_k_p_norm(compatible_input, k, p; allow_densify=allow_densify)
end

"""
    kpNorm(X::ComplexAffineMatrix, K, P; limits=OptimizationLimits())

Compatibility mapping for the pinned CVX-expression branch. The returned
[`QuantumEntanglementTools.TopKPNormEpigraph`](@ref) is solver-neutral model
data; it does not open a nested optimizer or claim a numerical norm value. Materialize it explicitly
with `add_top_k_p_norm_epigraph!` after loading the optional JuMP extension.
"""
function kpNorm(
    input::ComplexAffineMatrix, k, p; limits::OptimizationLimits=OptimizationLimits()
)
    return top_k_p_norm_epigraph(input, k, p; limits=limits)
end

"""
    kpNormDual(X, K, P; allow_densify=false)

Numeric-array compatibility entry point for the dual of [`kpNorm`](@ref).
Vector/matrix classification, magnitude ordering, clipping, and explicit
sparse-matrix densification follow `kpNorm`; CVX/model expressions are not
accepted.
"""
function kpNormDual(
    input::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k,
    p;
    allow_densify::Bool=false,
)
    compatible_input = _compat_kp_input(input)
    return top_k_p_norm_dual(compatible_input, k, p; allow_densify=allow_densify)
end

"""
    kpNormDual(X::ComplexAffineMatrix, K, P; limits=OptimizationLimits())

Compatibility mapping for the pinned CVX-expression branch. The result is an
exact solver-neutral epigraph atom; it does not start a nested optimizer or
claim a numerical norm value. Materialize it only after explicitly loading
the optional JuMP extension.
"""
function kpNormDual(
    input::ComplexAffineMatrix, k, p; limits::OptimizationLimits=OptimizationLimits()
)
    return top_k_p_norm_dual_epigraph(input, k, p; limits=limits)
end

function _compat_sk_tolerance(tolerance)
    tolerance === nothing && return (atol=nothing, rtol=nothing)
    tolerance isa Real && !(tolerance isa Bool) && isfinite(tolerance) ||
        throw(ArgumentError("TOL must be a finite nonnegative real number"))
    tolerance >= zero(tolerance) ||
        throw(ArgumentError("TOL must be a finite nonnegative real number"))
    return (atol=tolerance, rtol=zero(tolerance))
end

"""
    SkOperatorNorm(
        rng, X, K=1, DIM=nothing, STR=2, TARGET=-1,
        TOL=eps(Float64)^(3/8); structured=true, backend=NoOptimizationBackend(),
        kwargs...
    )

Compatibility spelling for QETLAB's randomized S(`k`) norm bounds. The
mandatory leading RNG replaces the pinned routine's implicit global random
stream. The default returns the status-rich native result. With
`structured=false`, the four QETLAB output positions are returned as
`(lb, lwit, ub, uwit)` without relabeling coincident numerical bounds as an
exact value.
"""
function SkOperatorNorm(
    rng::AbstractRNG,
    input::AbstractMatrix{<:Number},
    k=1,
    dim=nothing,
    strength=2,
    target=-1,
    tolerance=eps(Float64)^(3 / 8);
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    tolerances = _compat_sk_tolerance(tolerance)
    result = sk_operator_norm(
        rng,
        input;
        k=k,
        dims=dim,
        strength=strength,
        target=target,
        atol=tolerances.atol,
        rtol=tolerances.rtol,
        backend=backend,
        kwargs...,
    )
    structured && return result
    return (
        lb=result.lower_bound,
        lwit=result.lower_witness,
        ub=result.upper_bound,
        uwit=result.upper_witness,
    )
end

"""
    IsBlockPositive(
        rng, X, K=1, DIM=nothing, STR=2, TOL=eps(Float64)^(3/8);
        structured=true, backend=NoOptimizationBackend(), kwargs...
    )

Compatibility spelling for the tri-state block-positivity criterion. A
mandatory explicit RNG replaces transitive global randomness. Structured
results are returned by default. With `structured=false`, the two QETLAB
output positions are `(ibp, wit)`, where `ibp` is `1`, `0`, or `-1` for
certified true, certified false, or inconclusive respectively; `wit` is a
validated negative Schmidt-rank witness vector when available.
"""
function IsBlockPositive(
    rng::AbstractRNG,
    input::AbstractMatrix{<:Number},
    k=1,
    dim=nothing,
    strength=2,
    tolerance=eps(Float64)^(3 / 8);
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    tolerances = _compat_sk_tolerance(tolerance)
    result = is_block_positive(
        rng,
        input;
        k=k,
        dims=dim,
        strength=strength,
        atol=tolerances.atol,
        rtol=tolerances.rtol,
        backend=backend,
        kwargs...,
    )
    structured && return result
    answer = result.verdict === nothing ? -1 : (result.verdict ? 1 : 0)
    witness = result.witness === nothing ? nothing : result.witness.vector
    return (ibp=answer, wit=witness)
end

function _compat_induced_order(order, name::AbstractString)
    if order isa AbstractString
        lowercase(strip(order)) == "fro" ||
            throw(ArgumentError("$name string input must be \"fro\" (case-insensitive)"))
        return 2
    end
    return order
end

function _compat_induced_initial_vector(initial_vector)
    initial_vector === nothing && return nothing
    initial_vector isa Number && return nothing
    if initial_vector isa AbstractVector{<:Number}
        return initial_vector
    elseif initial_vector isa AbstractMatrix{<:Number}
        min(size(initial_vector)...) == 1 ||
            throw(DimensionMismatch("V0 must be a vector, row matrix, or column matrix"))
        return _compat_vectorize_kp_input(initial_vector)
    end
    return throw(
        ArgumentError("V0 must be numeric vector data or a scalar random-start sentinel")
    )
end

"""
    InducedMatrixNorm(
        rng, X, P, Q=P, TOL=nothing, V0=nothing;
        max_iterations=1000, max_work=1_000_000_000, allow_densify=false
    ) -> InducedMatrixNormResult

Numeric-array compatibility entry point for QETLAB's induced `P -> Q` norm
routine. `P` and `Q` accept numeric orders or case-insensitive `"fro"` for
order two. A scalar `V0` retains the pinned random-start sentinel convention;
row and column matrices are vectorized.

The mandatory leading `rng::AbstractRNG` replaces QETLAB's global `randn`
calls. An invalid non-scalar `V0` raises instead of silently warning and
switching to randomness. The structured result preserves the pinned optional
right-vector output and, critically, distinguishes exact closed-form values
from iterative lower bounds. Iteration convergence never implies exactness.
"""
function InducedMatrixNorm(
    rng::AbstractRNG,
    matrix::AbstractMatrix{<:Number},
    p,
    q=p,
    tolerance=nothing,
    initial_vector=nothing;
    max_iterations=1_000,
    max_work=1_000_000_000,
    allow_densify::Bool=false,
)
    checked_p = _compat_induced_order(p, "P")
    checked_q = _compat_induced_order(q, "Q")
    checked_initial_vector = _compat_induced_initial_vector(initial_vector)
    return induced_matrix_norm(
        rng,
        matrix,
        checked_p;
        q=checked_q,
        tolerance=tolerance,
        initial_vector=checked_initial_vector,
        max_iterations=max_iterations,
        max_work=max_work,
        allow_densify=allow_densify,
    )
end

function _compat_induced_initial_matrix(initial_matrix)
    initial_matrix === nothing && return nothing
    initial_matrix isa Number && return nothing
    initial_matrix isa AbstractMatrix{<:Number} || throw(
        ArgumentError("X0 must be a numeric matrix or a scalar random-start sentinel")
    )
    return initial_matrix
end

"""
    InducedSchattenNorm(
        rng, PHI, P, Q=P, DIM=nothing, TOL=nothing, X0=nothing;
        max_iterations=1000, max_work=1_000_000_000,
        max_dense_entries=1_000_000, allow_densify=false
    ) -> InducedSchattenNormResult

Compatibility entry point for QETLAB's numeric induced Schatten lower bound.
`PHI` may be a package map representation, numeric Choi matrix, Kraus
collection, or one/two-column factor cell matrix. `DIM` is validated against
the resulting square input/output algebras. Orders accept numbers or
case-insensitive `"fro"`.

The mandatory leading `rng` replaces the pinned global `randn`. A scalar `X0`
retains QETLAB's random-start sentinel; a malformed matrix raises instead of
silently switching to randomness. The returned structured result preserves
the optional witness while distinguishing the exact `2 -> 2` branch from
every iterative lower bound.
"""
function InducedSchattenNorm(
    rng::AbstractRNG,
    phi,
    p,
    q=p,
    dim=nothing,
    tolerance=nothing,
    initial_matrix=nothing;
    max_iterations=1_000,
    max_work=1_000_000_000,
    max_dense_entries=1_000_000,
    allow_densify::Bool=false,
)
    checked_p = _compat_induced_order(p, "P")
    checked_q = _compat_induced_order(q, "Q")
    map = _compat_map(phi, dim; allow_row_cp=false)
    checked_initial = _compat_induced_initial_matrix(initial_matrix)
    return induced_schatten_lower_bound(
        rng,
        map,
        checked_p;
        q=checked_q,
        initial_matrix=checked_initial,
        tolerance=tolerance,
        max_iterations=max_iterations,
        max_work=max_work,
        max_dense_entries=max_dense_entries,
        allow_densify=allow_densify,
    )
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

Numeric-array compatibility entry point for QETLAB von Neumann and Rényi
entropy. `ALPHA` may be any nonnegative finite real order or `Inf`; the input
receives strict density-matrix validation and is never normalized or clipped.
"""
function Entropy(rho::AbstractMatrix{<:Number}, base=2, alpha=1; allow_densify::Bool=false)
    return renyi_entropy(rho, alpha; base=base, allow_densify=allow_densify)
end

"""
    Fidelity(RHO, SIGMA; allow_densify=false)

Return QETLAB's unsquared Uhlmann root fidelity. Inputs receive the stricter
Julia-native density-matrix validation and are never repaired or normalized.
"""
function Fidelity(rho, sigma; allow_densify::Bool=false)
    return fidelity(rho, sigma; allow_densify=allow_densify)
end

"""
    MatsumotoFidelity(
        RHO, SIGMA; atol=nothing, rtol=nothing, allow_densify=false,
        support_boundary_policy=:reject
    )

Numeric-array compatibility entry point for `tr(RHO # SIGMA)`, where `#` is
the support-aware matrix geometric mean. Unlike the pinned fast branch, this
wrapper never adds an identity regularizer or forms an explicit inverse.
Singular support boundaries are rejected by default and can be projected only
with the named opt-in policy. The pinned CVX-expression branch is outside this
numeric compatibility wrapper.
"""
function MatsumotoFidelity(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    support_boundary_policy::Symbol=:reject,
)
    return matsumoto_fidelity(
        rho,
        sigma;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        support_boundary_policy=support_boundary_policy,
    )
end

"""
    MatsumotoFidelity(RHO::HermitianAffineMatrix,
                      SIGMA::HermitianAffineMatrix;
                      limits=OptimizationLimits())

Compatibility mapping for the pinned CVX-expression branch. The result is a
composable package-owned semidefinite lift, not a solver value. Positivity,
trace, and other constraints on symbolic state variables remain the caller's
explicit responsibility.
"""
function MatsumotoFidelity(
    rho::HermitianAffineMatrix,
    sigma::HermitianAffineMatrix;
    limits::OptimizationLimits=OptimizationLimits(),
)
    return matsumoto_fidelity_model(rho, sigma; limits=limits)
end

"""
    Distinguishability(X, P=nothing; structured=true, backend=NoOptimizationBackend(), ...)

Compatibility entry point for minimum-error state discrimination. The default
returns the status-rich native result. Set `structured=false` only when a
conclusive probability and residual-checked POVM are available; that form
returns `(dist, meas)` in QETLAB output order and otherwise throws a
`DomainError` containing the native result.

Unlike the pinned routine, states and priors are never silently normalized.
"""
function Distinguishability(
    states,
    priors=nothing;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    result = state_distinguishability(states; priors=priors, backend=backend, kwargs...)
    structured && return result
    result.success_probability === nothing && throw(
        DomainError(
            result,
            "Distinguishability did not establish a single success probability; " *
            "request structured=true to inspect bounds and solver status",
        ),
    )
    result.measurement === nothing && throw(
        DomainError(
            result,
            "Distinguishability has no residual-checked measurement; request " *
            "structured=true to inspect the result",
        ),
    )
    return (dist=result.success_probability, meas=result.measurement)
end

function _compat_channel_scalar(result, field::Symbol, name::AbstractString)
    value = getproperty(result, field)
    value === nothing && throw(
        DomainError(
            result,
            "$name did not establish a single value; request structured=true " *
            "to inspect bounds, statuses, residuals, and solver evidence",
        ),
    )
    return value
end

"""
    DiamondNorm(PHI, DIM=nothing; structured=true,
                backend=NoOptimizationBackend(), kwargs...)

Compatibility entry point for the completely bounded trace norm. `PHI` may
be a native map, numeric Choi matrix, Kraus collection, or one/two-column
factor cell. `DIM` follows the reviewed map-dimension compatibility contract.
The default preserves the native status-rich result; `structured=false`
returns QETLAB's scalar shape only when an analytic theorem or matching
residual-checked bounds establish a value.
"""
function DiamondNorm(
    phi,
    dim=nothing;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    map = _compat_map(phi, dim; allow_row_cp=false)
    result = diamond_norm(map; backend=backend, kwargs...)
    return structured ? result : _compat_channel_scalar(result, :value, "DiamondNorm")
end

"""
    CBNorm(PHI, DIM=nothing; structured=true,
           backend=NoOptimizationBackend(), kwargs...)

Compatibility entry point for the completely bounded operator norm. The
default keeps the native adjoint-reduction evidence and solver status;
`structured=false` returns a scalar only for a conclusive result.
"""
function CBNorm(
    phi,
    dim=nothing;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    map = _compat_map(phi, dim; allow_row_cp=false)
    result = cb_norm(map; backend=backend, kwargs...)
    return structured ? result : _compat_channel_scalar(result, :value, "CBNorm")
end

"""
    ChannelDistinguishability(
        PHI, PSI, P=nothing, DIM=nothing;
        structured=true, backend=NoOptimizationBackend(), kwargs...
    )

Compatibility spelling for one-use channel discrimination. Priors retain
QETLAB's third positional slot but are validated, never normalized, and must
be nonnegative. The native implementation applies the full channel
Holevo--Helstrom affine conversion, correcting the pinned routine's
non-deterministic probability defect. A legacy scalar is available only with
`structured=false` and a conclusive result.
"""
function ChannelDistinguishability(
    phi,
    psi,
    priors=nothing,
    dim=nothing;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    first_map = _compat_map(phi, dim; allow_row_cp=false)
    second_map = _compat_map(psi, dim; allow_row_cp=false)
    result = channel_distinguishability(
        first_map, second_map; priors=priors, backend=backend, kwargs...
    )
    return if structured
        result
    else
        _compat_channel_scalar(result, :success_probability, "ChannelDistinguishability")
    end
end

"""
    MaximumOutputFidelity(
        PHI, PSI; structured=true,
        backend=NoOptimizationBackend(), kwargs...
    )

Compatibility spelling for maximum output root fidelity. The direct native
SDP is invariant under Kraus representation and does not reproduce the
pinned unequal-Kraus-rank truncation defect. The structured result is the
default; `structured=false` returns the single QETLAB-shaped scalar only when
a value is established.
"""
function MaximumOutputFidelity(
    phi,
    psi;
    structured::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    first_map = _compat_map(phi; allow_row_cp=false)
    second_map = _compat_map(psi; allow_row_cp=false)
    result = maximum_output_fidelity(first_map, second_map; backend=backend, kwargs...)
    return if structured
        result
    else
        _compat_channel_scalar(result, :value, "MaximumOutputFidelity")
    end
end

function _compat_rounded_sqrt_dimension(total_dimension::Int)
    lower_root = isqrt(total_dimension)
    lower_distance = total_dimension - lower_root^2
    next_root_gap = 2 * lower_root + 1
    return 2 * lower_distance < next_root_gap ? lower_root : lower_root + 1
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
    SkVectorNorm(VEC, K=1, DIM=nothing; allow_densify=false)

Return the Euclidean norm of the `K` largest Schmidt coefficients of a
bipartite vector. A scalar `DIM` denotes the first subsystem dimension; an
omitted `DIM` uses QETLAB's nearest-integer `sqrt(length(VEC))` default and
then requires that value to divide the vector length. This can infer a
rectangular bipartition. `K` is clipped to the smaller subsystem dimension.
Sparse input is evaluated without densification when `K` covers the complete
Schmidt spectrum; a truncated computation requires `allow_densify=true`.
"""
function SkVectorNorm(
    state::AbstractVector{<:Number}, k=1, dim=nothing; allow_densify::Bool=false
)
    checked_dim = dim === nothing ? _compat_rounded_sqrt_dimension(length(state)) : dim
    dimensions = _compat_bipartite_dimensions(state, checked_dim, "SkVectorNorm")
    return schmidt_k_norm(state, dimensions, k; allow_densify=allow_densify)
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
                                 allow_densify=false,
                                 hermitian_factors=nothing)

Return a named tuple `(coefficients, left_factors, right_factors)`. `K=0`
keeps QETLAB's numerically nonzero terms, `K=-1` keeps the full thin
decomposition, and positive `K` keeps that many leading terms. A two-row
`DIM` matrix supplies independent local row and column dimensions. The named
tuple has exactly three fields and can be destructured positionally as
`s, U, V = OperatorSchmidtDecomposition(...)`.

By default, an exactly Hermitian operator with locally square dimensions uses
the native real Hermitian-basis SVD, so every returned factor is Hermitian.
This corrects the pinned branch's linear-indexing failure for unequal local
dimensions and reapplies `K` after the repair instead of accidentally ignoring
it. A nonzero Hermiticity residual is never projected away; pass
`hermitian_factors=true` to require the Hermitian convention or `false` to
disable it explicitly. Locally rectangular factor spaces cannot contain
Hermitian factors and use the general convention in automatic mode.
"""
function OperatorSchmidtDecomposition(
    operator::AbstractMatrix{<:Number},
    dim=nothing,
    k=0;
    allow_densify::Bool=false,
    hermitian_factors::Union{Nothing,Bool}=nothing,
)
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        operator, dim, "OperatorSchmidtDecomposition"; bipartite=true
    )
    k isa Integer && !(k isa Bool) ||
        throw(ArgumentError("K must be -1, 0, or a positive integer"))
    k >= -1 || throw(ArgumentError("K must be -1, 0, or a positive integer"))
    checked_k = try
        Int(k)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("K=$k cannot be represented as Int"))
    end
    use_hermitian_factors = if hermitian_factors === nothing
        row_dimensions == column_dimensions && ishermitian(operator)
    else
        hermitian_factors
    end
    decomposition = operator_schmidt_decomposition(
        operator,
        row_dimensions,
        column_dimensions;
        allow_densify=allow_densify,
        hermitian_factors=use_hermitian_factors,
    )
    coefficient_count = length(decomposition.coefficients)
    retained = if checked_k == -1
        coefficient_count
    elseif checked_k == 0
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
        checked_k <= coefficient_count || throw(
            ArgumentError(
                "K=$checked_k exceeds the full decomposition length " *
                "$coefficient_count",
            ),
        )
        checked_k
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
    OperatorSinkhorn(
        RHO, DIM=nothing, TOL=nothing;
        max_iterations=1000,
        allow_densify=false,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    )

Return exactly two named fields `(sigma, filters)`, which can also be
destructured positionally as `sigma, F = OperatorSinkhorn(...)`. `DIM`
retains QETLAB's default, scalar, and multipartite vector forms. `TOL`
defaults to `sqrt(eps(R))`, is used as the absolute aggregate marginal
residual threshold, and sets the pinned conditioning heuristic `1/TOL`.

Only a checked `:converged` native result is returned. Singular marginals,
ill-conditioning, an iteration or work limit, and numerical failures raise a
`DomainError` containing the structured native result. This wrapper never
changes global warning state and rejects `TOL <= 0` because the iteration is
always bounded.
"""
function OperatorSinkhorn(
    rho::AbstractMatrix{<:Number},
    dim=nothing,
    tolerance=nothing;
    max_iterations=1_000,
    allow_densify::Bool=false,
    max_entries=10_000_000,
    max_work=1_000_000_000,
)
    dim isa AbstractMatrix && throw(
        ArgumentError(
            "OperatorSinkhorn DIM must be a scalar or a subsystem vector, " *
            "not a two-row operator dimension matrix",
        ),
    )
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        rho, dim, "OperatorSinkhorn"; bipartite=false
    )
    row_dimensions == column_dimensions ||
        throw(DimensionMismatch("OperatorSinkhorn requires square local dimensions"))
    real_type = typeof(real(zero(eltype(rho))))
    real_type <: Union{Float32,Float64} || throw(
        ArgumentError(
            "OperatorSinkhorn requires Float32, Float64, ComplexF32, or " *
            "ComplexF64 input and never changes precision implicitly",
        ),
    )
    checked_tolerance = if tolerance === nothing
        sqrt(eps(real_type))
    else
        _compat_finite_real(tolerance, "TOL")
    end
    checked_tolerance > zero(checked_tolerance) ||
        throw(ArgumentError("TOL must be strictly positive for a bounded iteration"))
    result = operator_sinkhorn(
        rho,
        row_dimensions;
        atol=checked_tolerance,
        rtol=0,
        max_iterations=max_iterations,
        max_condition_number=one(checked_tolerance) / checked_tolerance,
        allow_densify=allow_densify,
        max_entries=max_entries,
        max_work=max_work,
    )
    result.converged || throw(
        DomainError(
            result,
            "OperatorSinkhorn ended with status $(result.status); inspect the " *
            "native result for residual and conditioning diagnostics",
        ),
    )
    return (sigma=result.scaled_operator, filters=result.local_filters)
end

"""
    FilterNormalForm(
        RHO, DIM=nothing, TOL=nothing;
        max_iterations=1000,
        allow_densify=false,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    )

Return exactly five named fields `(xi, GA, GB, FA, FB)`, which can be
destructured positionally in QETLAB output order. `DIM` retains the pinned
default, scalar, and two-element vector forms. A multipartite vector or
two-row operator-dimension matrix is rejected because the pinned
post-processing is bipartite.

`TOL` defaults to `sqrt(eps(R))` and is forwarded as the absolute aggregate
Sinkhorn residual threshold with condition limit `1/TOL`. This intentionally
corrects the pinned `FilterNormalForm.m`, which parses `TOL` but accidentally
omits it from its `OperatorSinkhorn` call. The native result retains the full
thin decomposition; this wrapper returns only terms above the reviewed
QETLAB default operator-Schmidt threshold.

Only a checked `:converged` native result is returned. Nonconvergence,
singular marginals, ill-conditioning, work exhaustion, or a failed
reconstruction check raises a `DomainError` containing the structured native
result.
"""
function FilterNormalForm(
    rho::AbstractMatrix{<:Number},
    dim=nothing,
    tolerance=nothing;
    max_iterations=1_000,
    allow_densify::Bool=false,
    max_entries=10_000_000,
    max_work=1_000_000_000,
)
    dim isa AbstractMatrix && throw(
        ArgumentError(
            "FilterNormalForm DIM must be a scalar or a two-element subsystem " *
            "vector, not a two-row operator dimension matrix",
        ),
    )
    row_dimensions, column_dimensions = _compat_product_operator_dimensions(
        rho, dim, "FilterNormalForm"; bipartite=true
    )
    row_dimensions == column_dimensions ||
        throw(DimensionMismatch("FilterNormalForm requires square local dimensions"))
    real_type = typeof(real(zero(eltype(rho))))
    real_type <: Union{Float32,Float64} || throw(
        ArgumentError(
            "FilterNormalForm requires Float32, Float64, ComplexF32, or " *
            "ComplexF64 input and never changes precision implicitly",
        ),
    )
    checked_tolerance = if tolerance === nothing
        sqrt(eps(real_type))
    else
        _compat_finite_real(tolerance, "TOL")
    end
    checked_tolerance > zero(checked_tolerance) ||
        throw(ArgumentError("TOL must be strictly positive for a bounded iteration"))

    result = filter_normal_form(
        rho,
        row_dimensions;
        balance_atol=checked_tolerance,
        balance_rtol=0,
        max_iterations=max_iterations,
        max_condition_number=one(checked_tolerance) / checked_tolerance,
        allow_densify=allow_densify,
        max_entries=max_entries,
        max_work=max_work,
    )
    result.converged || throw(
        DomainError(
            result,
            "FilterNormalForm ended with status $(result.status); inspect the " *
            "native result for rank, residual, conditioning, and work diagnostics",
        ),
    )
    retained = result.coefficient_numerical_rank
    indices = 1:retained
    return (
        xi=result.coefficients[indices],
        GA=result.left_operators[indices],
        GB=result.right_operators[indices],
        FA=result.left_filter,
        FB=result.right_filter,
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
                 allow_densify=false, psd_boundary_policy=:reject,
                 rank_boundary_policy=:reject,
                 range_boundary_policy=:reject)

Return base-two entanglement of formation for normalized bipartite pure
vectors, rank-one density matrices in arbitrary bipartite dimensions, and
validated two-qubit density matrices. Row and column pure vectors are
accepted. An omitted `DIM` uses QETLAB's nearest-integer square-root default
and can therefore infer a rectangular bipartition.

Unlike the pinned implementation's implicit numerical-rank projection,
boundary PSD, rank-one, and concurrence-range projections require the
corresponding explicit `:project` policy. The default `:reject` policies do
not repair boundary data. The corrected zero-concurrence limit is zero rather
than the pinned implementation's `NaN`.
"""
function EntFormation(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dim=nothing;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    psd_boundary_policy::Symbol=:reject,
    rank_boundary_policy::Symbol=:reject,
    range_boundary_policy::Symbol=:reject,
)
    input =
        state isa AbstractMatrix && min(size(state)...) == 1 ? _state_vector(state) : state
    total_dimension = input isa AbstractVector ? length(input) : size(input, 1)
    checked_dim = dim === nothing ? _compat_rounded_sqrt_dimension(total_dimension) : dim
    dimensions = _compat_bipartite_dimensions(input, checked_dim, "EntFormation")
    if input isa AbstractVector
        psd_boundary_policy === :reject ||
            throw(ArgumentError("psd_boundary_policy is meaningful only for matrix input"))
        rank_boundary_policy === :reject ||
            throw(ArgumentError("rank_boundary_policy is meaningful only for matrix input"))
        range_boundary_policy === :reject || throw(
            ArgumentError("range_boundary_policy is meaningful only for matrix input")
        )
        return entanglement_of_formation(
            input, dimensions; base=2, atol=atol, rtol=rtol, allow_densify=allow_densify
        )
    end
    return entanglement_of_formation(
        input,
        dimensions;
        base=2,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        psd_boundary_policy=psd_boundary_policy,
        rank_boundary_policy=rank_boundary_policy,
        range_boundary_policy=range_boundary_policy,
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

function _compat_majorization_accumulator_type(first_values, second_values)
    values = collect(Iterators.flatten((first_values, second_values)))
    isempty(values) && return BigInt
    exact = all(value -> value isa Integer || value isa Rational, values)
    types = exact ? typeof.(_compat_majorization_widen.(values)) : typeof.(values)
    return foldl(promote_type, types)
end

function _compat_majorization_accumulator_value(value, accumulator_type::Type)
    widened =
        value isa Integer || value isa Rational ? _compat_majorization_widen(value) : value
    return convert(accumulator_type, widened)
end

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

    accumulator_type = _compat_majorization_accumulator_type(first_values, second_values)
    first_prefix = zero(accumulator_type)
    second_prefix = zero(accumulator_type)
    for position in 1:common_length
        first_prefix += _compat_majorization_accumulator_value(
            first_values[position], accumulator_type
        )
        second_prefix += _compat_majorization_accumulator_value(
            second_values[position], accumulator_type
        )
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
    CompoundMatrix(
        A, R;
        sparse_output=issparse(A),
        max_entries=10_000_000,
        max_work=100_000_000,
    )

Return the `R`th multiplicative compound. Order zero is a `1×1` identity.
When `R > min(size(A)...)`, this compatibility entry point preserves pinned
QETLAB's `0×0` result. Other orders delegate to `compound_matrix`, including
its widened exact minor arithmetic and checked narrowing. Sparse inputs remain
sparse by default; `sparse_output` selects the representation explicitly. Use
the native function to retain a mathematically informative zero-by-nonzero
shape when the order exceeds only one dimension.

The Julia-only `max_entries` and `max_work` keywords are forwarded for every
nonempty native construction. The pinned early `0×0` compatibility branch
allocates no combinatorial workspace.
"""
function CompoundMatrix(
    matrix::AbstractMatrix,
    order;
    sparse_output::Bool=issparse(matrix),
    max_entries=10_000_000,
    max_work=100_000_000,
)
    checked_order = _nonnegative_dimension(order, "R")
    checked_order > min(size(matrix)...) && return Matrix{Float64}(undef, 0, 0)
    return compound_matrix(
        matrix,
        checked_order;
        sparse_output=sparse_output,
        max_entries=max_entries,
        max_work=max_work,
    )
end

"""
    AdditiveCompoundMatrix(
        A, R;
        sparse_output=issparse(A),
        max_entries=10_000_000,
        max_work=100_000_000,
    )

Return the `R`th additive compound. `A` must be square. Order zero reproduces
the reviewed pinned-QETLAB dependency-path error; use native
`additive_compound_matrix(A, 0)` for the mathematically defined `1×1` zero.
Positive orders delegate to the native widened exact arithmetic, and an order
above the dimension returns `0×0`. Sparse inputs remain sparse unless
`sparse_output` requests otherwise. The Julia-only `max_entries` and
`max_work` keywords retain the native preflight policy.
"""
function AdditiveCompoundMatrix(
    matrix::AbstractMatrix,
    order;
    sparse_output::Bool=issparse(matrix),
    max_entries=10_000_000,
    max_work=100_000_000,
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
    return additive_compound_matrix(
        matrix,
        checked_order;
        sparse_output=sparse_output,
        max_entries=max_entries,
        max_work=max_work,
    )
end

function _compat_commutant_all_sparse(input)
    input isa AbstractMatrix && return issparse(input)
    input isa Union{Tuple,AbstractVector} || return false
    isempty(input) && return false
    return all(matrix -> matrix isa AbstractMatrix && issparse(matrix), input)
end

"""
    Commutant(A;
              atol=0,
              rtol=nothing,
              allow_densify=false,
              max_entries=10_000_000,
              max_work=1_000_000_000)

QETLAB-compatible spelling for [`commutant`](@ref). `A` may be one square
matrix or a nonempty tuple/vector of equally sized square matrices. The return
value is a vector of Hilbert--Schmidt-orthonormal basis matrices.

This wrapper delegates the numerical calculation, tolerance convention, and
allocation guards to the native API. When every supplied generator is sparse,
`allow_densify=true` is still required for the bounded dense SVD, after which
the wrapper converts each returned basis matrix to sparse storage to preserve
QETLAB's output-storage convention. Such a basis can be structurally dense;
the native dense result is usually more efficient.
"""
function Commutant(
    input;
    atol=0,
    rtol=nothing,
    allow_densify::Bool=false,
    max_entries=10_000_000,
    max_work=1_000_000_000,
)
    basis = commutant(
        input;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_entries=max_entries,
        max_work=max_work,
    )
    return _compat_commutant_all_sparse(input) ? sparse.(basis) : basis
end

function _compat_polynomial_coefficients(coefficients::AbstractVector)
    Base.require_one_based_indexing(coefficients)
    return coefficients
end

function _compat_polynomial_coefficients(coefficients::AbstractMatrix)
    Base.require_one_based_indexing(coefficients)
    1 in size(coefficients) || throw(
        DimensionMismatch(
            "P must be a vector or a one-row/one-column matrix; got " *
            "size $(size(coefficients))",
        ),
    )
    return vec(coefficients)
end

function _compat_polynomial_sense(value)
    (value isa AbstractString || value isa Symbol) ||
        throw(ArgumentError("OPTTYPE must be \"min\" or \"max\"; got $(repr(value))"))
    sense = Symbol(lowercase(String(value)))
    sense in (:min, :max) ||
        throw(ArgumentError("OPTTYPE must be \"min\" or \"max\"; got $(repr(value))"))
    return sense
end

function _compat_polynomial_target(value)
    value === nothing && return nothing
    if value isa AbstractString || value isa Symbol
        lowercase(String(value)) == "none" || throw(
            ArgumentError(
                "TARGET must be \"none\" or a finite real number; got " * "$(repr(value))",
            ),
        )
        return nothing
    end
    value isa Real || throw(
        ArgumentError(
            "TARGET must be \"none\" or a finite real number; got $(repr(value))"
        ),
    )
    return value
end

"""
    CopositivePolynomial(C; dense_output=false, max_terms=100_000)

Return the QETLAB-ordered coefficient vector for
[`copositive_polynomial`](@ref). Sparse coefficients are retained by default.
Set `dense_output=true` to request QETLAB's dense vector explicitly; the
native `max_terms` guard is checked before either representation is created.

Unlike the pinned routine, this wrapper rejects a nonsymmetric matrix instead
of silently replacing it by `(C+C')/2`.
"""
function CopositivePolynomial(
    matrix::AbstractMatrix; dense_output::Bool=false, max_terms=100_000
)
    polynomial = copositive_polynomial(matrix; max_terms=max_terms)
    coefficients = copy(polynomial.coefficients)
    return dense_output ? collect(coefficients) : coefficients
end

"""
    PolynomialAsMatrix(
        P, N, D, K=0;
        sparse_output=true,
        max_terms=100_000,
        max_degree=256,
        max_dimension=10_000,
        max_exponent_entries=2_000_000,
        max_dense_entries=10_000_000,
        max_nonzeros=2_000_000,
        max_work=100_000_000,
    )

Numeric compatibility wrapper for [`polynomial_as_matrix`](@ref), preserving
QETLAB's positional `(P,N,D,K)` order and row/column coefficient vectors.
The result is sparse by default. `sparse_output=false` explicitly requests
guarded dense output. CVX expressions are outside the dependency-free numeric
wrapper.
"""
function PolynomialAsMatrix(
    coefficients::Union{AbstractVector,AbstractMatrix},
    variables,
    half_degree,
    level=0;
    kwargs...,
)
    return polynomial_as_matrix(
        _compat_polynomial_coefficients(coefficients),
        variables,
        half_degree;
        level=level,
        kwargs...,
    )
end

"""
    PolynomialOptimize(
        rng, P, N, D, K, OPTTYPE="max", TARGET="none";
        inner_samples=0,
        allow_densify=false,
        ...
    ) -> PolynomialOptimizationResult

Numeric compatibility wrapper for [`polynomial_bounds`](@ref). The mandatory
leading `rng::AbstractRNG` replaces QETLAB's global random stream, and
`inner_samples` is an exact deterministic count rather than a wall-clock
budget. Set `allow_densify=true` explicitly for the guarded generalized
eigensolver.

The structured result keeps hierarchy outer bounds separate from sampled
feasible inner values. Invalid `OPTTYPE` and `TARGET` values are rejected. The
maximization path corrects the pinned recursion defect by applying the target
in the original, unnegated objective convention.
"""
function PolynomialOptimize(
    rng::AbstractRNG,
    coefficients::Union{AbstractVector,AbstractMatrix},
    variables,
    half_degree,
    level,
    optimization_type="max",
    target="none";
    kwargs...,
)
    return polynomial_bounds(
        rng,
        _compat_polynomial_coefficients(coefficients),
        variables,
        half_degree,
        level;
        sense=_compat_polynomial_sense(optimization_type),
        target=_compat_polynomial_target(target),
        kwargs...,
    )
end

"""
    PolynomialSOS(
        rng, P, N, D, K, OPTTYPE="max", TARGET="none";
        backend=NoOptimizationBackend(), inner_samples=0, structured=true, ...
    )

Compatibility entry point for the pinned SOS hierarchy. The mandatory RNG and
exact `inner_samples` replace global, elapsed-time-dependent sampling.
`structured=true` returns
[`QuantumEntanglementTools.PolynomialSOSResult`](@ref), preserving optimizer
status, outer/inner bound directions, the moment matrix, and samples.

Set `structured=false` only after a usable outer bound exists; it returns
`(ob, ib)` in QETLAB output order and otherwise throws a `DomainError`
containing the structured result.
"""
function PolynomialSOS(
    rng::AbstractRNG,
    coefficients::Union{AbstractVector,AbstractMatrix},
    variables,
    half_degree,
    level,
    optimization_type="max",
    target="none";
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    inner_samples=0,
    structured::Bool=true,
    kwargs...,
)
    result = polynomial_sos_bounds(
        rng,
        _compat_polynomial_coefficients(coefficients),
        variables,
        half_degree,
        level;
        backend=backend,
        sense=_compat_polynomial_sense(optimization_type),
        target=_compat_polynomial_target(target),
        inner_samples=inner_samples,
        kwargs...,
    )
    structured && return result
    result.outer_bound === nothing && throw(
        DomainError(
            result,
            "PolynomialSOS has no usable outer bound; request structured=true " *
            "to inspect optimizer status and inner evidence",
        ),
    )
    return (ob=result.outer_bound, ib=result.inner_bound)
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
    IsPSD(X::HermitianAffineMatrix, TOL=nothing)

Return an owned solver-neutral PSD constraint for the pinned CVX-expression
branch. A numerical tolerance has no meaning for symbolic cone membership and
is rejected. Use the numeric-matrix method for a tri-state predicate.
"""
function IsPSD(matrix::HermitianAffineMatrix, tolerance=nothing)
    tolerance === nothing ||
        throw(ArgumentError("TOL is not accepted for an affine PSD model constraint"))
    return positive_semidefinite_constraint(matrix)
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
    RobkCohValue(V, K; atol=nothing, rtol=nothing) -> (ROB, L)

Return the two positional outputs of QETLAB's pure-state `k`-coherence
formula: robustness `ROB` and theorem branch index `L`. Julia vectors and
one-row or one-column matrices are accepted.

The wrapper delegates to [`pure_k_coherence_robustness`](@ref), so it
intentionally sorts coefficient magnitudes, supports complex phases, validates
normalization and `K`, and never normalizes the input. This corrects the pinned
routine's unsafe assumption that `V` is already a sorted, nonnegative,
normalized coefficient vector.
"""
function RobkCohValue(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k;
    atol=nothing,
    rtol=nothing,
)
    vector = _compat_coherence_state(state)
    vector isa AbstractVector ||
        throw(DimensionMismatch("V must be a vector, row matrix, or column matrix"))
    result = pure_k_coherence_robustness(vector, k; atol=atol, rtol=rtol)
    return result.value, result.branch_index
end

function _compat_coherence_conclusive(result, name::AbstractString)
    result.verdict === nothing && throw(
        DomainError(
            result,
            "$name is inconclusive; request structured=true to inspect theorem, " *
            "boundary, resource, and solver evidence",
        ),
    )
    return result.verdict ? 1 : 0
end

function _compat_coherence_value(result, name::AbstractString)
    result.value === nothing && throw(
        DomainError(
            result,
            "$name has no usable value; request structured=true to inspect the " *
            "optimization status and available evidence",
        ),
    )
    return result.value
end

"""
    IskIncoherent(X, K; structured=true, kwargs...)

Compatibility spelling for the coherence-number criterion. The default keeps
the native theorem/solver status and certificates. Set `structured=false` only
for a conclusive result, which returns QETLAB's `1` or `0`; inconclusive
boundaries and backend failures throw with the structured result attached.
"""
function IskIncoherent(state, k; structured::Bool=true, kwargs...)
    result = is_k_incoherent(_compat_coherence_state(state), k; kwargs...)
    return structured ? result : _compat_coherence_conclusive(result, "IskIncoherent")
end

"""
    IsAbskIncoh(X, K; structured=true, kwargs...)

Compatibility spelling for absolute `k`-incoherence. The default preserves
the native tri-state result rather than collapsing a one-sided theorem or
numerical boundary to a Boolean.
"""
function IsAbskIncoh(state, k; structured::Bool=true, kwargs...)
    result = is_absolutely_k_incoherent(_compat_coherence_state(state), k; kwargs...)
    return structured ? result : _compat_coherence_conclusive(result, "IsAbskIncoh")
end

"""
    RobustnessCoherence(RHO; structured=true, kwargs...)

Return the status-rich native robustness result by default. With
`structured=false`, return the scalar only when an analytic or residual-checked
optimizer branch produced one.
"""
function RobustnessCoherence(state; structured::Bool=true, kwargs...)
    result = robustness_coherence(_compat_coherence_state(state); kwargs...)
    return structured ? result : _compat_coherence_value(result, "RobustnessCoherence")
end

"""
    TraceDistanceCoherence(RHO; structured=true, kwargs...)

The default retains the closest state and optimization evidence. A conclusive
`structured=false` call returns `(TDC, D)` with `D` as the closest state's
diagonal vector, consistently across analytic and solver-backed branches.
"""
function TraceDistanceCoherence(state; structured::Bool=true, kwargs...)
    result = trace_distance_coherence(_compat_coherence_state(state); kwargs...)
    structured && return result
    value = _compat_coherence_value(result, "TraceDistanceCoherence")
    result.free_state === nothing &&
        throw(DomainError(result, "TraceDistanceCoherence has no validated closest state"))
    return (tdc=value, diagonal=real.(diag(result.free_state)))
end

"""
    GenRobustnesskCoherence(RHO, K; structured=true, kwargs...)

The structured default preserves the factor-width decomposition and solver
evidence. A conclusive legacy-shaped call returns `(robk, sig)`, where `sig`
is the normalized noise state or `nothing` when the exact robustness is zero.
"""
function GenRobustnesskCoherence(state, k; structured::Bool=true, kwargs...)
    result = generalized_robustness_k_coherence(
        _compat_coherence_state(state), k; kwargs...
    )
    structured && return result
    value = _compat_coherence_value(result, "GenRobustnesskCoherence")
    return (robk=value, sig=result.noise_state)
end

"""
    IsPPT(X, SYS=2, DIM=nothing, TOL=sqrt(eps(Float64));
          allow_densify=false) -> CriterionResult

Apply QETLAB's PPT test to a finite Hermitian matrix without requiring unit
trace. The compatibility result remains deliberately tri-state:
`CriterionEntanglementDetected`, `CriterionSatisfied`, or
`CriterionUnknown`. Values within `TOL` of the PSD boundary are never
collapsed to a Boolean. A nonzero Hermiticity residual within `TOL` returns
`CriterionUnknown` with residual evidence; a larger residual raises. The input
is never replaced by its Hermitian part.
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
    hermiticity_defect = dense - adjoint(dense)
    hermiticity_residual = maximum(
        abs, hermiticity_defect; init=zero(typeof(real(zero(eltype(dense)))))
    )
    hermiticity_residual <= checked_tolerance || throw(
        ArgumentError(
            "X is not Hermitian within TOL=$checked_tolerance; maximum " *
            "residual is $hermiticity_residual",
        ),
    )
    if !iszero(hermiticity_residual)
        index = argmax(abs.(hermiticity_defect))
        witness = (
            kind=:hermiticity_boundary,
            indices=Tuple(index),
            difference=hermiticity_defect[index],
            residual=hermiticity_residual,
        )
        return CriterionResult(
            :ppt,
            CriterionUnknown,
            hermiticity_residual,
            zero(hermiticity_residual),
            checked_tolerance,
            witness,
            "the Hermiticity residual lies within the numerical tolerance " *
            "boundary; the input was not symmetrized and no PPT conclusion is reported",
        )
    end
    transposed = partial_transpose(dense, dimensions; systems=systems)
    decomposition = eigen(Hermitian(transposed))
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

include("absolute_ppt.jl")
include("symmetric_extensions.jl")
include("copositivity_clique.jl")
include("separability_optimization.jl")
include("nonlocal_games.jl")

end # module MATLABCompat
