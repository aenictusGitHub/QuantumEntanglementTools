# Source-informed independent implementation of the public behavior in QETLAB
# IsBlockPositive.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2022 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    BlockPositivityStatus

Certificate-aware outcome of [`is_block_positive`](@ref). An optimizer limit,
floating-point boundary, or unsuccessful randomized search is inconclusive,
never a negative mathematical answer.
"""
@enum BlockPositivityStatus::UInt8 begin
    BlockPositivityCertified
    BlockPositivityViolated
    BlockPositivityUnknown
    BlockPositivityNumericalBoundary
    BlockPositivityResourceLimited
    BlockPositivityBackendFailure
end

"""
    BlockPositivityWitness

Explicit Schmidt-rank-`<= k` vector with a checked negative expectation. Such a
witness certifies failure of `k`-block positivity.
"""
struct BlockPositivityWitness{T<:Real,V<:AbstractVector{Complex{T}}}
    vector::V
    expectation::T
    schmidt_rank::Int
    kind::Symbol
    validation_residual::T
    tolerance::T
    validated::Bool
end

"""
    BlockPositivityResult

Tri-state, certificate-aware `k`-block-positivity result.

`verdict` is populated only by a checked spectral certificate, an explicit
negative Schmidt-rank witness, or a numerically validated S(`k`) upper
relaxation. `sk_norm_result` retains the full lower/upper-bound calculation and
backend metadata. A positive relaxation certificate is floating-point evidence
at the reported tolerance, not an exact symbolic proof.
"""
struct BlockPositivityResult{T<:Real,W,U,S}
    status::BlockPositivityStatus
    verdict::Union{Nothing,Bool}
    certificate_kind::Union{Nothing,Symbol}
    dimensions::NTuple{2,Int}
    k::Int
    witness::W
    upper_witness::U
    sk_norm_result::S
    minimum_eigenvalue::T
    operator_norm::T
    tolerance::T
    certified::Bool
    exact::Bool
    warnings::Tuple{Vararg{String}}
    message::String
end

function Base.show(io::IO, result::BlockPositivityResult)
    return print(
        io,
        "BlockPositivityResult(status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", k=",
        result.k,
        ", certified=",
        result.certified,
        ")",
    )
end

function _block_positive_witness(
    operator, vector, dimensions::NTuple{2,Int}, k::Int, tolerance, kind::Symbol
)
    T = _sknorm_real_type(eltype(operator))
    candidate = Complex{T}.(vector)
    candidate ./= norm(candidate)
    raw = dot(candidate, operator * candidate)
    expectation = convert(T, real(raw))
    imaginary_residual = convert(T, abs(imag(raw)))
    rank_value = schmidt_rank(candidate, dimensions; atol=zero(T), rtol=sqrt(eps(T)))
    validated =
        rank_value <= k && expectation < -tolerance && imaginary_residual <= tolerance
    return BlockPositivityWitness(
        candidate,
        expectation,
        rank_value,
        kind,
        imaginary_residual,
        convert(T, tolerance),
        validated,
    )
end

function _block_result(
    status,
    verdict,
    certificate_kind,
    dimensions,
    k,
    witness,
    upper_witness,
    sk_result,
    minimum_eigenvalue,
    operator_norm,
    tolerance,
    certified,
    exact,
    warnings,
    message,
)
    T = promote_type(typeof(minimum_eigenvalue), typeof(operator_norm), typeof(tolerance))
    return BlockPositivityResult{T,typeof(witness),typeof(upper_witness),typeof(sk_result)}(
        status,
        verdict,
        certificate_kind,
        dimensions,
        k,
        witness,
        upper_witness,
        sk_result,
        convert(T, minimum_eigenvalue),
        convert(T, operator_norm),
        convert(T, tolerance),
        certified,
        exact,
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

@doc raw"""
    is_block_positive(rng::AbstractRNG, operator; kwargs...)

Determine whether an exactly Hermitian bipartite operator is `k`-block
positive, returning a [`BlockPositivityResult`](@ref).

Positive definiteness is an immediate positive certificate. A negative
eigenvector or computational-basis vector certifies violation only when its
checked Schmidt rank is at most `k`. Otherwise the function applies the
identity

```math
X\text{ is }k\text{-block positive}
\quad\Longleftrightarrow\quad
\|\,\|X\|I-X\,\|_{S(k)} \le \|X\|
```

through [`sk_operator_norm`](@ref). A lower-bound witness can certify a
violation; a numerically validated SDP upper relaxation can certify positivity
away from the tolerance boundary. All other outcomes are `unknown`.

The randomized route requires the explicit `rng`; direct certificates and
preflight exits consume no random values. The input is never symmetrized,
normalized, clipped, or projected. Sparse input and dense/model work use the
same explicit limits as [`sk_operator_norm`](@ref).
"""
function is_block_positive(
    rng::AbstractRNG,
    operator::AbstractMatrix{<:Number};
    k=1,
    dims=nothing,
    strength=2,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_work=1_000_000_000,
    max_restarts=625,
    max_iterations=100,
    max_hierarchy_level=3,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_k = _sknorm_positive_integer(k, "k")
    prepared = _sknorm_prepare_operator(
        operator,
        dims;
        atol,
        rtol,
        allow_densify,
        max_dense_entries,
        max_work,
        require_hermitian=true,
    )
    dimensions = prepared.dimensions
    checked_k = min(checked_k, min(dimensions...))
    matrix = prepared.matrix
    T = prepared.real_type
    tolerance = prepared.tolerance
    operator_norm = prepared.operator_norm
    eigendecomposition = eigen(Hermitian(matrix))
    minimum_index = argmin(eigendecomposition.values)
    minimum_eigenvalue = convert(T, eigendecomposition.values[minimum_index])
    minimum_vector = eigendecomposition.vectors[:, minimum_index]
    warnings = String[]

    if all(iszero, matrix)
        return _block_result(
            BlockPositivityCertified,
            true,
            :zero_operator_exact,
            dimensions,
            checked_k,
            nothing,
            nothing,
            nothing,
            minimum_eigenvalue,
            operator_norm,
            tolerance,
            true,
            true,
            warnings,
            "the zero operator is exactly k-block positive",
        )
    end

    if minimum_eigenvalue > tolerance
        return _block_result(
            BlockPositivityCertified,
            true,
            :positive_definite_spectral_certificate,
            dimensions,
            checked_k,
            nothing,
            nothing,
            nothing,
            minimum_eigenvalue,
            operator_norm,
            tolerance,
            true,
            false,
            warnings,
            "positive definiteness implies k-block positivity",
        )
    end

    if minimum_eigenvalue < -tolerance
        eigen_witness = _block_positive_witness(
            matrix, minimum_vector, dimensions, checked_k, tolerance, :negative_eigenvector
        )
        if eigen_witness.validated
            return _block_result(
                BlockPositivityViolated,
                false,
                :negative_schmidt_rank_witness,
                dimensions,
                checked_k,
                eigen_witness,
                nothing,
                nothing,
                minimum_eigenvalue,
                operator_norm,
                tolerance,
                true,
                checked_k >= min(dimensions...),
                warnings,
                "a negative eigenvector has Schmidt rank at most k",
            )
        end
        diagonal_index = argmin(real.(diag(matrix)))
        if real(matrix[diagonal_index, diagonal_index]) < -tolerance
            basis_vector = zeros(Complex{T}, size(matrix, 1))
            basis_vector[diagonal_index] = one(Complex{T})
            basis_witness = _block_positive_witness(
                matrix,
                basis_vector,
                dimensions,
                checked_k,
                tolerance,
                :negative_computational_basis_expectation,
            )
            if basis_witness.validated
                return _block_result(
                    BlockPositivityViolated,
                    false,
                    :negative_product_vector_witness,
                    dimensions,
                    checked_k,
                    basis_witness,
                    nothing,
                    nothing,
                    minimum_eigenvalue,
                    operator_norm,
                    tolerance,
                    true,
                    false,
                    warnings,
                    "a computational product vector has negative expectation",
                )
            end
        end
    end

    if checked_k >= min(dimensions...)
        return _block_result(
            BlockPositivityNumericalBoundary,
            nothing,
            nothing,
            dimensions,
            checked_k,
            nothing,
            nothing,
            nothing,
            minimum_eigenvalue,
            operator_norm,
            tolerance,
            false,
            false,
            (
                "the PSD/block-positive decision lies inside the requested " *
                "eigenvalue tolerance band",
            ),
            "full-rank block positivity equals PSD, but the spectral result is " *
            "on the numerical boundary",
        )
    end

    identity_matrix = Matrix{eltype(matrix)}(I, size(matrix, 1), size(matrix, 1))
    shifted = operator_norm * identity_matrix - matrix
    sk_result = sk_operator_norm(
        rng,
        shifted;
        k=checked_k,
        dims=dimensions,
        strength,
        target=operator_norm,
        atol=zero(T),
        rtol=tolerance / max(one(T), operator_norm),
        allow_densify=true,
        max_dense_entries,
        max_work,
        max_restarts,
        max_iterations,
        max_hierarchy_level,
        backend,
        limits,
    )

    lower_witness = sk_result.lower_witness
    if lower_witness !== nothing &&
        lower_witness.quadratic &&
        lower_witness.validated &&
        lower_witness.value >= operator_norm + tolerance
        witness = _block_positive_witness(
            matrix,
            lower_witness.left,
            dimensions,
            checked_k,
            tolerance,
            :translated_sk_lower_witness,
        )
        if witness.validated
            return _block_result(
                BlockPositivityViolated,
                false,
                :negative_schmidt_rank_witness,
                dimensions,
                checked_k,
                witness,
                sk_result.upper_witness,
                sk_result,
                minimum_eigenvalue,
                operator_norm,
                tolerance,
                true,
                false,
                warnings,
                "an explicit S(k) lower-bound witness translates to a negative " *
                "expectation for the input operator",
            )
        end
    end

    upper_witness = sk_result.upper_witness
    if upper_witness !== nothing &&
        upper_witness.numerically_validated &&
        sk_result.upper_bound <= operator_norm - tolerance
        return _block_result(
            BlockPositivityCertified,
            true,
            :sk_sdp_upper_relaxation,
            dimensions,
            checked_k,
            nothing,
            upper_witness,
            sk_result,
            minimum_eigenvalue,
            operator_norm,
            tolerance,
            true,
            false,
            warnings,
            "a checked S(k) upper relaxation lies strictly below the block-" *
            "positivity threshold",
        )
    elseif sk_result.exact && sk_result.upper_bound <= operator_norm
        return _block_result(
            BlockPositivityCertified,
            true,
            :exact_sk_norm,
            dimensions,
            checked_k,
            nothing,
            upper_witness,
            sk_result,
            minimum_eigenvalue,
            operator_norm,
            tolerance,
            true,
            true,
            warnings,
            "an exact S(k)-norm branch meets the block-positivity threshold",
        )
    end

    status = if sk_result.status === SKOperatorNormResourceLimited
        BlockPositivityResourceLimited
    elseif sk_result.optimization_result !== nothing &&
        sk_result.optimization_result.status in (
        OptimizationMalformedBackend,
        OptimizationNumericalFailure,
        OptimizationInconsistent,
    )
        BlockPositivityBackendFailure
    elseif abs(minimum_eigenvalue) <= tolerance ||
        sk_result.upper_bound - sk_result.lower_bound <= tolerance
        BlockPositivityNumericalBoundary
    else
        BlockPositivityUnknown
    end
    append!(warnings, sk_result.warnings)
    push!(
        warnings,
        "necessary bounds or an unsuccessful heuristic search are not a " *
        "block-positivity certificate",
    )
    return _block_result(
        status,
        nothing,
        nothing,
        dimensions,
        checked_k,
        nothing,
        upper_witness,
        sk_result,
        minimum_eigenvalue,
        operator_norm,
        tolerance,
        false,
        false,
        warnings,
        "the available bounds overlap the block-positivity threshold",
    )
end
