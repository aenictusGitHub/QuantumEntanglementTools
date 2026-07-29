# Source-informed independent Julia implementations based on the
# specifications in QETLAB L1NormCoherence.m, RelEntCoherence.m, and
# CoherenceRank.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and named coauthors,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

function _coherence_density(
    state::AbstractVector{<:Number}; atol, rtol, operation::AbstractString
)
    validation = _tierd_validate_pure_state(
        state; atol=atol, rtol=rtol, operation=operation
    )
    return (
        vector=state,
        norm_squared=validation.norm_squared,
        atol=validation.atol,
        rtol=validation.rtol,
    )
end

function _coherence_density(
    state::AbstractMatrix{<:Number};
    atol,
    rtol,
    allow_densify::Bool,
    operation::AbstractString,
    boundary_policy::Symbol=:reject,
)
    return _tierd_density_analysis(
        state;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=boundary_policy,
        operation=operation,
    )
end

"""
    l1_coherence(state; atol=nothing, rtol=nothing, allow_densify=false)

Return the `l1`-norm of coherence in the computational basis: the sum of the
absolute off-diagonal entries of the density operator. `state` may be a
normalized pure vector or a validated density matrix.

For a pure vector the equivalent allocation-free identity
`sum(abs, state)^2 - 1` is used, including for sparse vectors. Matrix inputs
are validated as finite, normalized, Hermitian, and positive semidefinite;
sparse matrices require explicit `allow_densify=true` because validation uses
a full eigendecomposition. Inputs are never normalized, clipped, or mutated.
"""
function l1_coherence(
    state::AbstractVector{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    _ = allow_densify
    analysis = _coherence_density(state; atol=atol, rtol=rtol, operation="l1_coherence")
    return sum(abs, analysis.vector)^2 - analysis.norm_squared
end

function l1_coherence(
    state::AbstractMatrix{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    analysis = _coherence_density(
        state;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="l1_coherence",
        boundary_policy=:record,
    )
    matrix = analysis.matrix
    return sum(abs, matrix) - sum(abs, diag(matrix))
end

function _coherence_entropy(probabilities, base)
    logarithm_base = _tierd_validate_log_base(base)
    all(isfinite, probabilities) ||
        throw(ArgumentError("probabilities must contain only finite values"))
    all(probability -> probability >= zero(probability), probabilities) || throw(
        DomainError(
            minimum(probabilities),
            "entropy probabilities must be nonnegative; values are never clipped",
        ),
    )
    entropy = zero(eltype(probabilities))
    for probability in probabilities
        iszero(probability) && continue
        entropy -= probability * log(probability)
    end
    return entropy / log(logarithm_base)
end

"""
    relative_entropy_coherence(state; base, atol=nothing, rtol=nothing,
                               allow_densify=false)

Return `S(diag(state)) - S(state)` in the computational basis, where `S` is
the von Neumann entropy in the explicitly supplied logarithm `base`.

For a normalized pure vector this reduces to the Shannon entropy of
`abs2.(state)`. Density matrices receive the same strict validation and sparse
densification gate as [`von_neumann_entropy`](@ref). No probability,
eigenvalue, or final result is clipped.
"""
function relative_entropy_coherence(
    state::AbstractVector{<:Number};
    base,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    _ = allow_densify
    analysis = _coherence_density(
        state; atol=atol, rtol=rtol, operation="relative_entropy_coherence"
    )
    return _coherence_entropy(abs2.(analysis.vector), base)
end

function relative_entropy_coherence(
    state::AbstractMatrix{<:Number};
    base,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    analysis = _coherence_density(
        state;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="relative_entropy_coherence",
    )
    diagonal_probabilities = real.(diag(analysis.matrix))
    diagonal_entropy = _coherence_entropy(diagonal_probabilities, base)
    state_entropy = _coherence_entropy(analysis.eigenvalues, base)
    return diagonal_entropy - state_entropy
end

"""
    coherence_rank(state; basis=nothing, atol=nothing, rtol=nothing,
                   allow_densify=false)

Return the numerical coherence rank of a normalized pure vector: the number of
basis coefficients strictly larger than
`atol + rtol * maximum(abs, coefficients)`.

`basis=nothing` selects the computational basis and preserves sparse-vector
storage. A supplied basis must be a finite square unitary matrix whose columns
are the basis vectors; coordinates are computed with `basis \\ state`, never
with an explicit inverse. Sparse coordinate transforms require
`allow_densify=true` because the result can be dense.
"""
function coherence_rank(
    state::AbstractVector{<:Number};
    basis=nothing,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    validation = _tierd_validate_pure_state(
        state; atol=atol, rtol=rtol, operation="coherence_rank"
    )
    coefficients = if basis === nothing
        state
    else
        basis isa AbstractMatrix{<:Number} ||
            throw(ArgumentError("basis must be a numeric matrix or `nothing`"))
        Base.require_one_based_indexing(basis)
        size(basis) == (length(state), length(state)) || throw(
            DimensionMismatch(
                "basis has size $(size(basis)); expected " *
                "($(length(state)), $(length(state)))",
            ),
        )
        all(isfinite, basis) ||
            throw(ArgumentError("basis must contain only finite entries"))
        (issparse(state) || issparse(basis)) &&
            !allow_densify &&
            throw(
                ArgumentError(
                    "a sparse basis-coordinate transform may produce dense output; " *
                    "pass allow_densify=true to permit it",
                ),
            )
        identity_matrix = Matrix{eltype(basis)}(I, length(state), length(state))
        isapprox(
            adjoint(basis) * basis,
            identity_matrix;
            atol=validation.atol,
            rtol=validation.rtol,
        ) || throw(
            ArgumentError("basis columns are not unitary within the requested tolerances"),
        )
        coordinate_basis = issparse(basis) ? Matrix(basis) : basis
        coordinate_state = issparse(state) ? Vector(state) : state
        coordinate_basis \ coordinate_state
    end
    real_type = typeof(real(zero(eltype(coefficients))))
    scale = maximum(abs, coefficients; init=zero(real_type))
    threshold = validation.atol + validation.rtol * scale
    return count(value -> abs(value) > threshold, coefficients)
end
