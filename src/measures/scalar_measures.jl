# Source-informed independent Julia implementations based on the specifications
# in QETLAB TraceNorm.m, SchattenNorm.m, KyFanNorm.m, kpNorm.m, kpNormDual.m,
# InducedMatrixNorm.m, Purity.m, Entropy.m, Fidelity.m, MatsumotoFidelity.m,
# Negativity.m, SchmidtDecomposition.m, SkVectorNorm.m, SchmidtRank.m, and
# Concurrence.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and named coauthors,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

using Random: AbstractRNG, randn

"""
    SchmidtDecompositionResult

Thin Schmidt decomposition returned by [`schmidt_decomposition`](@ref).
If `result = schmidt_decomposition(ψ, (dA, dB))`, then

```julia
ψ ≈ sum(result.coefficients[k] *
        kron(result.left_vectors[:, k], result.right_vectors[:, k])
        for k in eachindex(result.coefficients))
```

The coefficient vector contains all `min(dA, dB)` singular values, including
numerical or exact zeros.  Use [`schmidt_rank`](@ref) to apply an explicit
rank tolerance.
"""
struct SchmidtDecompositionResult{C<:AbstractVector,L<:AbstractMatrix,R<:AbstractMatrix}
    coefficients::C
    left_vectors::L
    right_vectors::R
end

"""
    InducedMatrixNormResult

Structured result returned by [`induced_matrix_norm`](@ref).

`bound_kind` is `:exact` only for a proved closed-form branch (or the zero
operator); otherwise it is `:lower_bound`, including when the bounded
iteration converges. `status` is one of `:exact`, `:exact_zero`,
`:converged_lower_bound`, `:iteration_limit`, or `:work_limit`.

`witness` is a normalized right-multiplication vector. Its independently
recorded `witness_input_norm` and `witness_output_norm` obey
`witness_output_norm / witness_input_norm <= ||matrix||_{p->q}`.
`normalization_residual` and `value_residual` make that lower-bound evidence
auditable. `iteration_residual` is the absolute change in the objective at
the final completed alternating update; it is `nothing` when no update ran.
Convergence of the iteration never changes a lower bound into an exact value.
"""
struct InducedMatrixNormResult{R<:AbstractFloat,P<:Real,Q<:Real,W<:AbstractVector}
    value::R
    p::P
    q::Q
    bound_kind::Symbol
    status::Symbol
    exact::Bool
    converged::Bool
    witness::W
    witness_input_norm::R
    witness_output_norm::R
    normalization_residual::R
    value_residual::R
    iteration_residual::Union{Nothing,R}
    tolerance::R
    iterations::Int
    work_used::BigInt
    max_iterations::Int
    max_work::Union{Nothing,BigInt}
    message::String
end

function Base.show(io::IO, result::InducedMatrixNormResult)
    return print(
        io,
        "InducedMatrixNormResult(value=",
        result.value,
        ", p=",
        result.p,
        ", q=",
        result.q,
        ", bound_kind=",
        result.bound_kind,
        ", status=",
        result.status,
        ", iterations=",
        result.iterations,
        ")",
    )
end

function _tierd_validate_tolerance(value, name::AbstractString)
    value === nothing && return nothing
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative real number or `nothing`"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= zero(value) || throw(ArgumentError("$name must be nonnegative; got $value"))
    return value
end

function _tierd_tolerances(::Type{R}, atol, rtol) where {R<:AbstractFloat}
    absolute = _tierd_validate_tolerance(atol, "atol")
    relative = _tierd_validate_tolerance(rtol, "rtol")
    return (
        absolute === nothing ? zero(R) : absolute,
        relative === nothing ? sqrt(eps(R)) : relative,
    )
end

function _tierd_real_type(::Type{T}) where {T<:LinearAlgebra.BlasFloat}
    return typeof(real(zero(T)))
end

function _tierd_require_spectral_eltype(array, operation::AbstractString)
    eltype(array) <: LinearAlgebra.BlasFloat || throw(
        ArgumentError(
            "$operation requires Float32, Float64, ComplexF32, or ComplexF64 " *
            "input in the dependency-free core; got eltype $(eltype(array)). " *
            "No implicit precision-changing conversion is performed.",
        ),
    )
    return _tierd_real_type(eltype(array))
end

function _tierd_require_pure_state_eltype(array, operation::AbstractString)
    real_type = typeof(real(zero(eltype(array))))
    real_type <: AbstractFloat || throw(
        ArgumentError(
            "$operation requires a real or complex floating-point pure-state " *
            "element type; got eltype $(eltype(array)). No implicit " *
            "precision-changing conversion is performed.",
        ),
    )
    return real_type
end

function _tierd_require_finite(array, name::AbstractString)
    Base.require_one_based_indexing(array)
    all(isfinite, array) || throw(ArgumentError("$name must contain only finite entries"))
    return nothing
end

function _tierd_dense_matrix(
    matrix::AbstractMatrix{<:Number}; allow_densify::Bool, operation::AbstractString
)
    if issparse(matrix) && !allow_densify
        throw(
            ArgumentError(
                "$operation requires a dense spectral factorization; pass " *
                "`allow_densify=true` to explicitly permit converting this " *
                "$(size(matrix)) sparse matrix to dense storage",
            ),
        )
    end
    _tierd_require_finite(matrix, "matrix")
    real_type = _tierd_require_spectral_eltype(matrix, operation)
    return Matrix(matrix), real_type
end

function _tierd_dense_vector(
    vector::AbstractVector{<:Number}; allow_densify::Bool, operation::AbstractString
)
    if issparse(vector) && !allow_densify
        throw(
            ArgumentError(
                "$operation requires a dense spectral factorization; pass " *
                "`allow_densify=true` to explicitly permit converting this " *
                "length-$(length(vector)) sparse vector to dense storage",
            ),
        )
    end
    _tierd_require_finite(vector, "vector")
    _tierd_require_spectral_eltype(vector, operation)
    return Vector(vector)
end

function _tierd_scale(array, ::Type{R}) where {R<:AbstractFloat}
    entry_scale = maximum(abs, array; init=zero(R))
    return max(one(R), entry_scale)
end

_tierd_threshold(scale, atol, rtol) = atol + rtol * max(one(scale), scale)

function _tierd_singular_values(
    matrix::AbstractMatrix{<:Number}; allow_densify::Bool, operation::AbstractString
)
    dense, _ = _tierd_dense_matrix(matrix; allow_densify=allow_densify, operation=operation)
    isempty(dense) && return zeros(_tierd_real_type(eltype(dense)), 0)
    return svdvals(dense)
end

function _tierd_validate_log_base(base)
    base isa Real && !(base isa Bool) ||
        throw(ArgumentError("base must be a finite real number greater than 1"))
    isfinite(base) && base > one(base) ||
        throw(ArgumentError("base must be finite and greater than 1; got $(repr(base))"))
    return base
end

function _tierd_density_analysis(
    rho::AbstractMatrix{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    boundary_policy::Symbol=:reject,
    operation::AbstractString="this operation",
)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("rho must be square; got size $(size(rho))"))
    !isempty(rho) || throw(ArgumentError("rho must have positive dimension"))
    dense, real_type = _tierd_dense_matrix(
        rho; allow_densify=allow_densify, operation=operation
    )
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    scale = _tierd_scale(dense, real_type)
    validation_tolerance = _tierd_threshold(scale, absolute, relative)

    hermiticity_residual = maximum(abs, dense - adjoint(dense); init=zero(real_type))
    hermiticity_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho is not Hermitian within atol=$absolute and rtol=$relative; " *
            "maximum residual is $hermiticity_residual",
        ),
    )

    trace_value = tr(dense)
    trace_imaginary_residual = abs(imag(trace_value))
    trace_imaginary_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho has a non-real trace outside tolerance; trace(rho)=$trace_value"
        ),
    )
    real_trace = real(trace_value)
    normalization_residual = abs(real_trace - one(real_trace))
    normalization_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho is not normalized within atol=$absolute and rtol=$relative; " *
            "trace(rho)=$trace_value. The input is never normalized implicitly.",
        ),
    )

    # Forming the Hermitian part after the explicit residual check avoids
    # handing a nearly-Hermitian matrix to a one-triangle eigensolver.  This is
    # a documented numerical work matrix; the caller's input is not mutated.
    work_matrix = (dense + adjoint(dense)) / 2
    decomposition = eigen(Hermitian(work_matrix))
    eigenvalues = decomposition.values
    spectral_scale = max(scale, maximum(abs, eigenvalues; init=zero(real_type)))
    spectral_tolerance = _tierd_threshold(spectral_scale, absolute, relative)
    minimum_eigenvalue = minimum(eigenvalues)
    minimum_eigenvalue < -spectral_tolerance && throw(
        DomainError(
            minimum_eigenvalue,
            "rho is not positive semidefinite within atol=$absolute and " *
            "rtol=$relative",
        ),
    )

    spectral_boundary_uncertain = minimum_eigenvalue < zero(minimum_eigenvalue)
    if spectral_boundary_uncertain && boundary_policy === :reject
        throw(
            DomainError(
                minimum_eigenvalue,
                "rho has a small negative eigenvalue inside the validation " *
                "tolerance. Refusing to clip it to zero; tighten the input or " *
                "use a status-returning criterion.",
            ),
        )
    elseif boundary_policy !== :reject && boundary_policy !== :record
        throw(
            ArgumentError(
                "internal boundary policy must be :reject or :record; got $boundary_policy"
            ),
        )
    end
    structural_boundary_uncertain =
        spectral_boundary_uncertain ||
        !iszero(hermiticity_residual) ||
        !iszero(trace_imaginary_residual)
    boundary_uncertain = structural_boundary_uncertain || !iszero(normalization_residual)

    return (
        matrix=Matrix(work_matrix),
        decomposition=decomposition,
        eigenvalues=eigenvalues,
        trace_value=real_trace,
        atol=absolute,
        rtol=relative,
        tolerance=spectral_tolerance,
        boundary_uncertain=boundary_uncertain,
        structural_boundary_uncertain=structural_boundary_uncertain,
        spectral_boundary_uncertain=spectral_boundary_uncertain,
        minimum_eigenvalue=minimum_eigenvalue,
        hermiticity_residual=hermiticity_residual,
        trace_imaginary_residual=trace_imaginary_residual,
        normalization_residual=normalization_residual,
    )
end

function _tierd_validate_pure_state(
    psi::AbstractVector{<:Number}; atol=nothing, rtol=nothing, operation::AbstractString
)
    !isempty(psi) || throw(ArgumentError("psi must have positive length"))
    _tierd_require_finite(psi, "psi")
    real_type = _tierd_require_pure_state_eltype(psi, operation)
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    norm_squared = real(dot(psi, psi))
    tolerance = _tierd_threshold(max(one(real_type), abs(norm_squared)), absolute, relative)
    abs(norm_squared - one(norm_squared)) <= tolerance || throw(
        ArgumentError(
            "psi is not normalized within atol=$absolute and rtol=$relative; " *
            "squared norm is $norm_squared. The input is never normalized implicitly.",
        ),
    )
    return (norm_squared=norm_squared, atol=absolute, rtol=relative)
end

function _tierd_bipartite_layout(dims, total_dimension::Int)
    layout = _as_layout(dims)
    length(layout) == 2 || throw(
        ArgumentError("dims must describe exactly two subsystems; got $(layout.dims)")
    )
    layout.total_dimension == total_dimension || throw(
        DimensionMismatch(
            "prod(dims)=$(layout.total_dimension) does not match state " *
            "dimension $total_dimension",
        ),
    )
    return layout
end

"""
    trace_norm(matrix; allow_densify=false)

Return the Schatten 1-norm, the sum of all singular values.

The dependency-free implementation uses a dense full singular-value
decomposition.  Sparse inputs are therefore rejected unless
`allow_densify=true` is supplied explicitly.  Entries must be finite and the
element type must be a BLAS floating type; no precision-changing conversion is
performed.  For an `m x n` matrix the full SVD costs
`O(m*n*min(m,n))` time and `O(m*n)` dense workspace.
"""
function trace_norm(matrix::AbstractMatrix{<:Number}; allow_densify::Bool=false)
    return sum(
        _tierd_singular_values(matrix; allow_densify=allow_densify, operation="trace_norm")
    )
end

"""
    schatten_norm(matrix, p; allow_densify=false)

Return the Schatten `p`-norm for real `p >= 1`, including `p == Inf`.
Sparse inputs require explicit `allow_densify=true`.  This implementation
computes the full SVD with the complexity described for [`trace_norm`](@ref).
"""
function schatten_norm(matrix::AbstractMatrix{<:Number}, p; allow_densify::Bool=false)
    p isa Real && !(p isa Bool) ||
        throw(ArgumentError("p must be a real number in [1, Inf]"))
    !isnan(p) && p >= one(p) ||
        throw(ArgumentError("p must lie in [1, Inf]; got $(repr(p))"))
    singular_values = _tierd_singular_values(
        matrix; allow_densify=allow_densify, operation="schatten_norm"
    )
    isempty(singular_values) && return zero(eltype(singular_values))
    isinf(p) && return maximum(singular_values)
    isfinite(p) || throw(ArgumentError("p must be finite or Inf; got $(repr(p))"))
    real_type = eltype(singular_values)
    exponent = try
        convert(real_type, p)
    catch error
        error isa InexactError || rethrow()
        throw(
            ArgumentError(
                "p=$(repr(p)) cannot be represented in the matrix's $real_type precision",
            ),
        )
    end
    isfinite(exponent) || throw(
        ArgumentError(
            "p=$(repr(p)) cannot be represented finitely in the matrix's $real_type precision",
        ),
    )
    return sum(value -> value^exponent, singular_values)^(inv(exponent))
end

"""
    ky_fan_norm(matrix, k; allow_densify=false)

Return the sum of the `k` largest singular values.  `k` must be an integer in
`1:min(size(matrix)...)`.  Sparse inputs require explicit
`allow_densify=true`.  This dependency-free implementation computes all
singular values rather than an iterative truncated SVD.
"""
function ky_fan_norm(matrix::AbstractMatrix{<:Number}, k; allow_densify::Bool=false)
    k isa Integer && !(k isa Bool) || throw(ArgumentError("k must be a positive integer"))
    maximum_k = min(size(matrix)...)
    1 <= k <= maximum_k || throw(
        ArgumentError(
            "k=$k is outside the valid range 1:$maximum_k for size $(size(matrix))"
        ),
    )
    singular_values = _tierd_singular_values(
        matrix; allow_densify=allow_densify, operation="ky_fan_norm"
    )
    return sum(@view singular_values[1:Int(k)])
end

function _tierd_validate_top_k(k, maximum_k::Int)
    k isa Integer && !(k isa Bool) || throw(ArgumentError("k must be a positive integer"))
    k >= one(k) || throw(ArgumentError("k must be positive; got $k"))
    return k >= maximum_k ? maximum_k : Int(k)
end

function _tierd_validate_p_order(p)
    p isa Real && !(p isa Bool) ||
        throw(ArgumentError("p must be a real number in [1, Inf]"))
    !isnan(p) && p >= one(p) ||
        throw(ArgumentError("p must lie in [1, Inf]; got $(repr(p))"))
    (isfinite(p) || p == Inf) ||
        throw(ArgumentError("p must be finite or positive Inf; got $(repr(p))"))
    return p
end

function _tierd_sorted_vector_magnitudes(vector::AbstractVector{<:Number})
    _tierd_require_finite(vector, "vector")
    values = if vector isa AbstractSparseVector
        map(value -> float(abs(value)), nonzeros(vector))
    else
        map(value -> float(abs(value)), vector)
    end
    sort!(values; rev=true)
    return values, length(vector)
end

function _tierd_sorted_matrix_singular_values(
    matrix::AbstractMatrix{<:Number}; allow_densify::Bool, operation::AbstractString
)
    if matrix isa Diagonal
        _tierd_require_finite(matrix, "matrix")
        values = [float(abs(value)) for value in diag(matrix)]
        sort!(values; rev=true)
        return values, length(values)
    end
    values = _tierd_singular_values(
        matrix; allow_densify=allow_densify, operation=operation
    )
    return values, length(values)
end

function _tierd_top_k_p_from_spectrum(values, ambient_length::Int, k, p)
    effective_k = _tierd_validate_top_k(k, ambient_length)
    selected_count = min(effective_k, length(values))
    return norm(@view(values[1:selected_count]), p)
end

function _tierd_spectrum_value(values, index::Int)
    return index <= length(values) ? values[index] : zero(eltype(values))
end

function _tierd_plateau_norm(values, prefix_length::Int, average, copies::Int, q)
    explicit_prefix = min(prefix_length, length(values))
    scale = maximum(@view(values[1:explicit_prefix]); init=zero(eltype(values)))
    scale = max(scale, average)
    iszero(scale) && return zero(scale)
    isinf(q) && return scale

    scaled_sum = zero((scale / scale)^q)
    for index in 1:explicit_prefix
        scaled_sum += (values[index] / scale)^q
    end
    scaled_sum += copies * (average / scale)^q
    return scale * scaled_sum^(inv(q))
end

function _tierd_top_k_p_dual_from_spectrum(values, ambient_length::Int, k, p)
    effective_k = _tierd_validate_top_k(k, ambient_length)
    iszero(ambient_length) && return norm(values, p)
    isinf(p) && return norm(values, 1)

    r = min(effective_k - 1, length(values))
    tail_sum = if effective_k <= length(values)
        sum(@view values[effective_k:length(values)])
    else
        zero(eltype(values))
    end
    average = tail_sum / (effective_k - r)
    while r > 0
        _tierd_spectrum_value(values, r) > average && break
        tail_sum += _tierd_spectrum_value(values, r)
        r -= 1
        average = tail_sum / (effective_k - r)
    end

    if p == one(p)
        return max(_tierd_spectrum_value(values, 1), average)
    end
    q = p / (p - one(p))
    return _tierd_plateau_norm(values, r, average, effective_k - r, q)
end

"""
    top_k_p_norm(vector, k, p; allow_densify=false)
    top_k_p_norm(matrix, k, p; allow_densify=false)

Return the `p`-norm of the `k` largest magnitudes of a vector or the `k`
largest singular values of a matrix. `k` must be positive and is clipped to
the available vector length or singular-value count, matching the reviewed
numeric behavior of QETLAB `kpNorm`. `p` must lie in `[1, Inf]`.

Sparse vectors are scanned through their stored values and are never
densified; implicit zeros do not affect the result. General sparse matrices
require `allow_densify=true` because the dependency-free core computes a full
SVD. `Diagonal` matrices use their diagonal entries directly, including
`BigFloat` entries, without forming a dense matrix. Other matrix
factorizations require a BLAS floating element type. No CVX or symbolic-model
expression is accepted.

The pinned QETLAB numeric vector branch sorts signed entries rather than their
magnitudes. This implementation follows the documented norm definition, so
negative and complex vector entries are ordered by `abs`.

# Examples

```jldoctest
julia> top_k_p_norm([-5.0, 4.0, 3.0], 2, 1)
9.0
```

Dense-vector sorting costs `O(n log n)` time and `O(n)` workspace; a sparse
vector uses `O(nnz log nnz)` time and `O(nnz)` workspace. General matrix
complexity is dominated by the full SVD described for [`trace_norm`](@ref).
"""
function top_k_p_norm(vector::AbstractVector{<:Number}, k, p; allow_densify::Bool=false)
    checked_p = _tierd_validate_p_order(p)
    values, ambient_length = _tierd_sorted_vector_magnitudes(vector)
    return _tierd_top_k_p_from_spectrum(values, ambient_length, k, checked_p)
end

function top_k_p_norm(matrix::AbstractMatrix{<:Number}, k, p; allow_densify::Bool=false)
    checked_p = _tierd_validate_p_order(p)
    values, ambient_length = _tierd_sorted_matrix_singular_values(
        matrix; allow_densify=allow_densify, operation="top_k_p_norm"
    )
    return _tierd_top_k_p_from_spectrum(values, ambient_length, k, checked_p)
end

"""
    top_k_p_norm_dual(vector, k, p; allow_densify=false)
    top_k_p_norm_dual(matrix, k, p; allow_densify=false)

Return the dual norm of [`top_k_p_norm`](@ref). The input spectrum is sorted
in descending order. For finite `p`, the implementation uses the verified
QETLAB plateau formula; `p == 1` yields
`max(maximum(spectrum), sum(spectrum) / k)` at its selected plateau and
`p == Inf` yields the full vector `1`-norm or matrix trace norm. `k` clipping,
vector magnitude ordering, sparse-storage behavior, supported element types,
and the exclusion of symbolic model expressions are the same as for
[`top_k_p_norm`](@ref).

# Examples

```jldoctest
julia> isapprox(top_k_p_norm_dual([3.0, 2.0, 1.0], 2, 2), sqrt(18))
true
```

After obtaining the sorted spectrum, the plateau search is linear in the
stored spectrum length; runs of implicit sparse zeros are skipped.
"""
function top_k_p_norm_dual(
    vector::AbstractVector{<:Number}, k, p; allow_densify::Bool=false
)
    checked_p = _tierd_validate_p_order(p)
    values, ambient_length = _tierd_sorted_vector_magnitudes(vector)
    return _tierd_top_k_p_dual_from_spectrum(values, ambient_length, k, checked_p)
end

function top_k_p_norm_dual(
    matrix::AbstractMatrix{<:Number}, k, p; allow_densify::Bool=false
)
    checked_p = _tierd_validate_p_order(p)
    values, ambient_length = _tierd_sorted_matrix_singular_values(
        matrix; allow_densify=allow_densify, operation="top_k_p_norm_dual"
    )
    return _tierd_top_k_p_dual_from_spectrum(values, ambient_length, k, checked_p)
end

const _TIERD_INDUCED_DEFAULT_MAX_ITERATIONS = 1_000
const _TIERD_INDUCED_DEFAULT_MAX_WORK = 1_000_000_000

function _tierd_induced_order(
    order, name::AbstractString, ::Type{R}
) where {R<:AbstractFloat}
    order isa Real && !(order isa Bool) ||
        throw(ArgumentError("$name must be a real number in [1, Inf]"))
    (isfinite(order) || order == Inf) && order >= one(order) ||
        throw(ArgumentError("$name must lie in [1, Inf]; got $(repr(order))"))
    checked = try
        convert(R, order)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$(repr(order)) cannot be represented as $R"))
    end
    isfinite(order) &&
        !isfinite(checked) &&
        throw(ArgumentError("$name=$(repr(order)) overflows the matrix real type $R"))
    return checked
end

function _tierd_induced_tolerance(value, ::Type{R}) where {R<:AbstractFloat}
    value === nothing && return sqrt(eps(R))
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("tolerance must be a finite nonnegative real number"))
    isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("tolerance must be finite and nonnegative; got $(repr(value))"))
    checked = try
        convert(R, value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("tolerance=$(repr(value)) cannot be represented as $R"))
    end
    isfinite(checked) ||
        throw(ArgumentError("tolerance=$(repr(value)) overflows the matrix real type $R"))
    return checked
end

function _tierd_induced_iterations(value)
    value isa Bool &&
        throw(ArgumentError("max_iterations must be a nonnegative integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("max_iterations must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("max_iterations must be nonnegative; got $value"))
    return try
        Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("max_iterations=$value cannot be represented as Int"))
    end
end

function _tierd_induced_work_limit(value)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("max_work must be a positive integer or nothing, not Bool"))
    value isa Integer ||
        throw(ArgumentError("max_work must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("max_work must be positive; got $value"))
    return BigInt(value)
end

function _tierd_induced_require_work(required::BigInt, limit, context::AbstractString)
    limit !== nothing &&
        required > limit &&
        throw(
            ArgumentError(
                "$context requires $required deterministic work units, exceeding " *
                "max_work=$limit",
            ),
        )
    return required
end

function _tierd_induced_argmax_abs(vector::AbstractVector)
    best_index = firstindex(vector)
    best_value = abs(vector[best_index])
    for index in Iterators.drop(eachindex(vector), 1)
        candidate = abs(vector[index])
        if candidate > best_value
            best_index = index
            best_value = candidate
        end
    end
    return best_index, best_value
end

_tierd_induced_phase(value) = iszero(value) ? one(value) : value / abs(value)

function _tierd_induced_left_dual(vector::AbstractVector, q)
    index, scale = _tierd_induced_argmax_abs(vector)
    if q == one(q)
        return map(_tierd_induced_phase, vector)
    elseif isinf(q)
        dual = zeros(eltype(vector), length(vector))
        dual[index] = _tierd_induced_phase(vector[index])
        return dual
    elseif iszero(scale)
        dual = zeros(eltype(vector), length(vector))
        dual[firstindex(dual)] = one(eltype(dual))
        return dual
    end

    magnitudes = map(value -> (abs(value) / scale)^(q - one(q)), vector)
    dual_order = q / (q - one(q))
    denominator = norm(magnitudes, dual_order)
    return map(_tierd_induced_phase, vector) .* magnitudes ./ denominator
end

function _tierd_induced_right_dual(vector::AbstractVector, p)
    index, scale = _tierd_induced_argmax_abs(vector)
    if p == one(p)
        dual = zeros(eltype(vector), length(vector))
        dual[index] = _tierd_induced_phase(vector[index])
        return dual
    elseif isinf(p)
        return map(_tierd_induced_phase, vector)
    elseif iszero(scale)
        dual = zeros(eltype(vector), length(vector))
        dual[firstindex(dual)] = one(eltype(dual))
        return dual
    end

    magnitudes = map(value -> (abs(value) / scale)^(inv(p - one(p))), vector)
    denominator = norm(magnitudes, p)
    return map(_tierd_induced_phase, vector) .* magnitudes ./ denominator
end

function _tierd_induced_start_vector(
    rng::AbstractRNG, matrix::AbstractMatrix, p, initial_vector, ::Type{R}
) where {R<:AbstractFloat}
    columns = size(matrix, 2)
    matrix_type = eltype(matrix)
    vector = if initial_vector === nothing
        if matrix_type <: Real
            randn(rng, R, columns)
        else
            complex.(randn(rng, R, columns), randn(rng, R, columns))
        end
    else
        initial_vector isa AbstractVector{<:Number} ||
            throw(ArgumentError("initial_vector must be a numeric vector or nothing"))
        Base.require_one_based_indexing(initial_vector)
        length(initial_vector) == columns || throw(
            DimensionMismatch(
                "initial_vector must have length $columns; got $(length(initial_vector))",
            ),
        )
        all(isfinite, initial_vector) ||
            throw(ArgumentError("initial_vector must contain only finite entries"))
        initial_real_type = typeof(real(zero(eltype(initial_vector))))
        initial_real_type <: AbstractFloat &&
            initial_real_type !== R &&
            throw(
                ArgumentError(
                    "initial_vector has real precision $initial_real_type but the " *
                    "matrix uses $R; no implicit precision-changing conversion is performed",
                ),
            )
        working_type =
            if matrix_type <: Real && any(value -> !isreal(value), initial_vector)
                Complex{R}
            else
                matrix_type
            end
        try
            working_type.(initial_vector)
        catch err
            err isa InexactError || rethrow()
            throw(
                ArgumentError(
                    "initial_vector entries cannot be represented in working type " *
                    "$working_type",
                ),
            )
        end
    end

    input_norm = norm(vector, p)
    isfinite(input_norm) || throw(ArgumentError("initial_vector has a non-finite p-norm"))
    if iszero(input_norm)
        initial_vector === nothing || throw(
            DomainError(input_norm, "initial_vector must have strictly positive p-norm")
        )
        vector[firstindex(vector)] = one(eltype(vector))
        input_norm = one(R)
    end
    return vector ./ input_norm
end

function _tierd_induced_result(
    value,
    p,
    q,
    bound_kind::Symbol,
    status::Symbol,
    witness::AbstractVector,
    witness_output_norm,
    iteration_residual,
    tolerance,
    iterations::Int,
    work_used::BigInt,
    max_iterations::Int,
    max_work,
    message::String,
)
    real_type = typeof(value)
    witness_input_norm = convert(real_type, norm(witness, p))
    checked_output_norm = convert(real_type, witness_output_norm)
    checked_value = convert(real_type, value)
    return InducedMatrixNormResult(
        checked_value,
        p,
        q,
        bound_kind,
        status,
        bound_kind === :exact,
        status in (:exact, :exact_zero, :converged_lower_bound),
        witness,
        witness_input_norm,
        checked_output_norm,
        abs(witness_input_norm - one(real_type)),
        abs(checked_output_norm - checked_value),
        iteration_residual,
        tolerance,
        iterations,
        work_used,
        max_iterations,
        max_work,
        message,
    )
end

function _tierd_induced_exact_one(matrix, p, q, tolerance, max_iterations, max_work)
    rows, columns = size(matrix)
    entry_work = BigInt(rows) * columns
    work_used = _tierd_induced_require_work(
        2 * entry_work + rows + columns, max_work, "the exact 1-to-1 induced norm branch"
    )
    column_sums = vec(sum(abs, matrix; dims=1))
    value, index = findmax(column_sums)
    witness = zeros(eltype(matrix), columns)
    witness[index] = one(eltype(matrix))
    output_norm = norm(matrix * witness, q)
    return _tierd_induced_result(
        value,
        p,
        q,
        :exact,
        :exact,
        witness,
        output_norm,
        nothing,
        tolerance,
        0,
        work_used,
        max_iterations,
        max_work,
        "the maximum absolute column sum is the exact 1-to-1 induced norm",
    )
end

function _tierd_induced_exact_two(
    matrix, p, q, tolerance, max_iterations, max_work, allow_densify::Bool
)
    rows, columns = size(matrix)
    issparse(matrix) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "the exact 2-to-2 branch requires a dense SVD; pass " *
                "`allow_densify=true` to permit converting this $(size(matrix)) " *
                "sparse matrix to dense storage",
            ),
        )
    entry_work = BigInt(rows) * columns
    work_used = _tierd_induced_require_work(
        entry_work * max(1, min(rows, columns)) + entry_work + rows + columns,
        max_work,
        "the exact 2-to-2 induced norm branch",
    )
    decomposition = svd(Matrix(matrix); full=false)
    value = first(decomposition.S)
    witness = Vector(adjoint(decomposition.Vt)[:, 1])
    output_norm = norm(matrix * witness, q)
    return _tierd_induced_result(
        value,
        p,
        q,
        :exact,
        :exact,
        witness,
        output_norm,
        nothing,
        tolerance,
        0,
        work_used,
        max_iterations,
        max_work,
        "the largest singular value is the exact 2-to-2 induced norm",
    )
end

function _tierd_induced_exact_infinity(matrix, p, q, tolerance, max_iterations, max_work)
    rows, columns = size(matrix)
    entry_work = BigInt(rows) * columns
    work_used = _tierd_induced_require_work(
        2 * entry_work + rows + columns,
        max_work,
        "the exact infinity-to-infinity induced norm branch",
    )
    row_sums = vec(sum(abs, matrix; dims=2))
    value, row = findmax(row_sums)
    witness = Vector{eltype(matrix)}(undef, columns)
    for column in axes(matrix, 2)
        entry = matrix[row, column]
        witness[column] = iszero(entry) ? one(entry) : conj(entry) / abs(entry)
    end
    output_norm = norm(matrix * witness, q)
    return _tierd_induced_result(
        value,
        p,
        q,
        :exact,
        :exact,
        witness,
        output_norm,
        nothing,
        tolerance,
        0,
        work_used,
        max_iterations,
        max_work,
        "the maximum absolute row sum is the exact infinity-to-infinity induced norm",
    )
end

"""
    induced_matrix_norm(
        rng, matrix, p; q=p, tolerance=nothing, initial_vector=nothing,
        max_iterations=1000, max_work=1_000_000_000, allow_densify=false
    ) -> InducedMatrixNormResult

Compute an exact induced matrix norm in the proved `1 -> 1`, `2 -> 2`, and
`Inf -> Inf` cases. All other orders use the pinned alternating
Hölder-equality iteration and return a witnessed lower bound. An explicit
`rng::AbstractRNG` is mandatory. Exact branches and a supplied
`initial_vector` do not consume it, and the global random stream is never
read or mutated.

The matrix must be nonempty, finite, and have element type `Float32`,
`Float64`, `ComplexF32`, or `ComplexF64`. Orders lie in `[1, Inf]`.
`tolerance` defaults to `sqrt(eps(R))` and tests the absolute objective change
between completed updates. `max_iterations` and `max_work` are deterministic
budgets; no wall clock is consulted. A work budget too small to complete an
exact branch or initialize lower-bound evidence raises before calculation.
After lower-bound initialization, reaching either budget preserves the best
witness found so far and records `:iteration_limit` or `:work_limit`.

General sparse matrices stay sparse throughout matrix-vector products. The
exact `1 -> 1` and `Inf -> Inf` branches also preserve sparse storage. The
exact `2 -> 2` branch requires a full SVD and therefore requires
`allow_densify=true` for sparse input. No element-type conversion or
normalization of `matrix` is performed.

For lower-bound results, `converged` means only that the alternating
objective change met `tolerance`; it is not a certificate of the global
maximum. Inspect `bound_kind`, `status`, `witness`, the witness norms and
residuals, `iterations`, and `work_used`.

# Examples

```jldoctest
julia> using Random

julia> result = induced_matrix_norm(
           MersenneTwister(1),
           [1.0 -2.0; 3.0 4.0; -5.0 1.0],
           1,
       );

julia> (result.value, result.exact, result.status, result.witness)
(9.0, true, :exact, [1.0, 0.0])
```
"""
function induced_matrix_norm(
    rng::AbstractRNG,
    matrix::AbstractMatrix{<:Number},
    p;
    q=p,
    tolerance=nothing,
    initial_vector=nothing,
    max_iterations=_TIERD_INDUCED_DEFAULT_MAX_ITERATIONS,
    max_work=_TIERD_INDUCED_DEFAULT_MAX_WORK,
    allow_densify::Bool=false,
)
    Base.require_one_based_indexing(matrix)
    rows, columns = size(matrix)
    rows > 0 && columns > 0 ||
        throw(ArgumentError("matrix must have positive row and column dimensions"))
    _tierd_require_finite(matrix, "matrix")
    real_type = _tierd_require_spectral_eltype(matrix, "induced_matrix_norm")
    checked_p = _tierd_induced_order(p, "p", real_type)
    checked_q = _tierd_induced_order(q, "q", real_type)
    checked_tolerance = _tierd_induced_tolerance(tolerance, real_type)
    checked_iterations = _tierd_induced_iterations(max_iterations)
    checked_work = _tierd_induced_work_limit(max_work)

    if checked_p == one(checked_p) && checked_q == one(checked_q)
        return _tierd_induced_exact_one(
            matrix,
            checked_p,
            checked_q,
            checked_tolerance,
            checked_iterations,
            checked_work,
        )
    elseif checked_p == real_type(2) && checked_q == real_type(2)
        return _tierd_induced_exact_two(
            matrix,
            checked_p,
            checked_q,
            checked_tolerance,
            checked_iterations,
            checked_work,
            allow_densify,
        )
    elseif isinf(checked_p) && isinf(checked_q)
        return _tierd_induced_exact_infinity(
            matrix,
            checked_p,
            checked_q,
            checked_tolerance,
            checked_iterations,
            checked_work,
        )
    end

    entry_work = BigInt(rows) * columns
    initialization_work = entry_work + rows + columns
    _tierd_induced_require_work(
        initialization_work, checked_work, "the lower-bound initialization"
    )

    if all(iszero, matrix)
        witness = zeros(eltype(matrix), columns)
        witness[firstindex(witness)] = one(eltype(witness))
        return _tierd_induced_result(
            zero(real_type),
            checked_p,
            checked_q,
            :exact,
            :exact_zero,
            witness,
            zero(real_type),
            nothing,
            checked_tolerance,
            0,
            initialization_work,
            checked_iterations,
            checked_work,
            "the zero operator has exact induced norm zero for every p and q",
        )
    end

    current_witness = _tierd_induced_start_vector(
        rng, matrix, checked_p, initial_vector, real_type
    )
    current_output = matrix * current_witness
    current_value = convert(real_type, norm(current_output, checked_q))
    best_witness = copy(current_witness)
    best_value = current_value
    best_output_norm = current_value
    work_used = initialization_work
    iteration_work = 2 * entry_work + 3 * (rows + columns)
    last_residual = nothing
    completed_iterations = 0

    while completed_iterations < checked_iterations
        if checked_work !== nothing && work_used + iteration_work > checked_work
            return _tierd_induced_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :work_limit,
                best_witness,
                best_output_norm,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                "the deterministic work budget was reached; the witnessed " *
                "value remains a lower bound, not an exact norm",
            )
        end

        left_witness = _tierd_induced_left_dual(current_output, checked_q)
        right_gradient = adjoint(matrix) * left_witness
        new_witness = _tierd_induced_right_dual(right_gradient, checked_p)
        new_output = matrix * new_witness
        new_value = convert(real_type, norm(new_output, checked_q))
        isfinite(new_value) || throw(
            DomainError(
                new_value, "the induced-norm iteration produced a non-finite objective"
            ),
        )

        completed_iterations += 1
        work_used += iteration_work
        last_residual = abs(new_value - current_value)
        if new_value > best_value
            best_witness = copy(new_witness)
            best_value = new_value
            best_output_norm = new_value
        end
        current_witness = new_witness
        current_output = new_output
        current_value = new_value

        if last_residual <= checked_tolerance
            return _tierd_induced_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :converged_lower_bound,
                best_witness,
                best_output_norm,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                "the alternating objective change met tolerance; the witnessed " *
                "value is still only a lower bound",
            )
        end
    end

    return _tierd_induced_result(
        best_value,
        checked_p,
        checked_q,
        :lower_bound,
        :iteration_limit,
        best_witness,
        best_output_norm,
        last_residual,
        checked_tolerance,
        completed_iterations,
        work_used,
        checked_iterations,
        checked_work,
        "the iteration budget was reached; the witnessed value remains a " *
        "lower bound, not an exact norm",
    )
end

"""
    purity(rho; atol=nothing, rtol=nothing, allow_densify=false)

Return `tr(rho^2)` for a validated density matrix.

`atol=nothing` means zero absolute tolerance and `rtol=nothing` means
`sqrt(eps(R))`, where `R` is the real floating type of `rho`.  The input must
be finite, square, normalized, approximately Hermitian, and positive
semidefinite.  It is never normalized, clipped, or mutated.  After the
Hermiticity residual passes the requested tolerance, the explicitly formed
Hermitian part is used as the numerical work matrix.  The dense Hermitian
eigendecomposition costs `O(n^3)` time and `O(n^2)` workspace.
"""
function purity(
    rho::AbstractMatrix{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    analysis = _tierd_density_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="purity"
    )
    return sum(abs2, analysis.eigenvalues)
end

function _tierd_validate_entropy_order(alpha)
    alpha isa Real && !(alpha isa Bool) ||
        throw(ArgumentError("alpha must be a real number in [0, Inf]"))
    !isnan(alpha) && alpha >= zero(alpha) ||
        throw(ArgumentError("alpha must lie in [0, Inf]; got $(repr(alpha))"))
    (isfinite(alpha) || alpha == Inf) ||
        throw(ArgumentError("alpha must be finite or positive Inf; got $(repr(alpha))"))
    return alpha
end

function _tierd_diagonal_density_eigenvalues(
    rho::Diagonal{<:Number}; atol, rtol, operation::AbstractString
)
    !isempty(rho) || throw(ArgumentError("rho must have positive dimension"))
    _tierd_require_finite(rho, "rho")
    real_type = typeof(real(zero(eltype(rho))))
    real_type <: AbstractFloat || throw(
        ArgumentError(
            "$operation requires a real or complex floating-point element " *
            "type; got eltype $(eltype(rho)). No implicit precision-changing " *
            "conversion is performed.",
        ),
    )
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    diagonal_values = diag(rho)
    scale = max(one(real_type), maximum(abs, diagonal_values; init=zero(real_type)))
    validation_tolerance = _tierd_threshold(scale, absolute, relative)

    hermiticity_residual = maximum(
        value -> abs(value - conj(value)), diagonal_values; init=zero(real_type)
    )
    hermiticity_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho is not Hermitian within atol=$absolute and rtol=$relative; " *
            "maximum residual is $hermiticity_residual",
        ),
    )

    trace_value = sum(diagonal_values)
    trace_imaginary_residual = abs(imag(trace_value))
    trace_imaginary_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho has a non-real trace outside tolerance; trace(rho)=$trace_value"
        ),
    )
    real_trace = real(trace_value)
    normalization_residual = abs(real_trace - one(real_trace))
    normalization_residual <= validation_tolerance || throw(
        ArgumentError(
            "rho is not normalized within atol=$absolute and rtol=$relative; " *
            "trace(rho)=$trace_value. The input is never normalized implicitly.",
        ),
    )

    eigenvalues = real.(diagonal_values)
    spectral_scale = max(scale, maximum(abs, eigenvalues; init=zero(real_type)))
    spectral_tolerance = _tierd_threshold(spectral_scale, absolute, relative)
    minimum_eigenvalue = minimum(eigenvalues)
    minimum_eigenvalue < -spectral_tolerance && throw(
        DomainError(
            minimum_eigenvalue,
            "rho is not positive semidefinite within atol=$absolute and " *
            "rtol=$relative",
        ),
    )
    minimum_eigenvalue < zero(minimum_eigenvalue) && throw(
        DomainError(
            minimum_eigenvalue,
            "rho has a small negative eigenvalue inside the validation " *
            "tolerance. Refusing to clip it to zero.",
        ),
    )
    return eigenvalues
end

function _tierd_entropy_eigenvalues(
    rho::AbstractMatrix{<:Number};
    atol,
    rtol,
    allow_densify::Bool,
    operation::AbstractString,
)
    if rho isa Diagonal
        return _tierd_diagonal_density_eigenvalues(
            rho; atol=atol, rtol=rtol, operation=operation
        )
    end
    analysis = _tierd_density_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation=operation
    )
    return analysis.eigenvalues
end

function _tierd_entropy_from_eigenvalues(eigenvalues, base, alpha)
    working_type = promote_type(
        eltype(eigenvalues), typeof(float(base)), typeof(float(alpha))
    )
    checked_base = convert(working_type, base)
    checked_alpha = convert(working_type, alpha)
    logarithm_base = log(checked_base)
    positive = working_type[
        convert(working_type, value) for value in eigenvalues if value > zero(value)
    ]
    isempty(positive) && throw(
        DomainError(
            eigenvalues,
            "rho has empty positive spectral support and no entropy is defined",
        ),
    )

    if iszero(checked_alpha)
        return log(convert(working_type, length(positive))) / logarithm_base
    elseif checked_alpha == one(checked_alpha)
        numerator = zero(working_type)
        for value in positive
            numerator -= value * log(value)
        end
        return numerator / logarithm_base
    elseif isinf(checked_alpha)
        return -log(maximum(positive)) / logarithm_base
    end

    maximum_log = maximum(log, positive)
    denominator = one(checked_alpha) - checked_alpha
    correction = sum(value -> exp(checked_alpha * (log(value) - maximum_log)), positive)
    natural_entropy =
        (checked_alpha / denominator) * maximum_log + log(correction) / denominator
    return natural_entropy / logarithm_base
end

"""
    renyi_entropy(rho, alpha; base, atol=nothing, rtol=nothing,
                  allow_densify=false)

Return the Rényi entropy of order `alpha` for a validated density matrix:

```math
H_\\alpha(\\rho) =
\\frac{\\log\\!\\left(\\sum_i \\lambda_i^\\alpha\\right)}
     {(1-\\alpha)\\log(\\mathrm{base})}.
```

The continuous endpoint definitions are used exactly: `alpha == 0` is the
logarithm of the positive spectral support, `alpha == 1` is the von Neumann
entropy, and `alpha == Inf` is `-log(maximum(eigenvalues)) / log(base)`.
Finite orders must be nonnegative, and `base` is a required keyword greater
than one. Zero eigenvalues never enter a logarithm.

The input must be finite, square, Hermitian, positive semidefinite, and
normalized within the requested tolerances. It is never normalized or
clipped. After the documented Hermiticity-residual check, the Hermitian work
spectrum is used without mutating the caller's matrix. General
sparse inputs require `allow_densify=true`; the
dependency-free eigensolver otherwise supports BLAS floating element types.
`Diagonal` states are evaluated directly in `O(n)` time and workspace and
therefore preserve generic floating types such as `BigFloat` without
densification. A general dense state uses `O(n^3)` time and `O(n^2)`
workspace.

# Examples

```jldoctest
julia> using LinearAlgebra

julia> renyi_entropy(Diagonal([0.5, 0.25, 0.25]), 2; base=2)
1.4150374992788437
```
"""
function renyi_entropy(
    rho::AbstractMatrix{<:Number},
    alpha;
    base,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    checked_base = _tierd_validate_log_base(base)
    checked_alpha = _tierd_validate_entropy_order(alpha)
    eigenvalues = _tierd_entropy_eigenvalues(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="renyi_entropy"
    )
    return _tierd_entropy_from_eigenvalues(eigenvalues, checked_base, checked_alpha)
end

"""
    von_neumann_entropy(rho; base, atol=nothing, rtol=nothing,
                        allow_densify=false)

Return `-tr(rho * log(rho)) / log(base)` for a validated density matrix. This
is [`renyi_entropy(rho, 1; base=base)`](@ref renyi_entropy).
The logarithm base is a required keyword and must be finite and greater than
one.  Exact zero eigenvalues contribute zero.  Negative
eigenvalues, even if small enough to pass a coarse PSD tolerance, are not
clipped and cause an error.

The dense Hermitian eigendecomposition costs `O(n^3)` time and `O(n^2)`
workspace.
"""
function von_neumann_entropy(
    rho::AbstractMatrix{<:Number};
    base,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    return renyi_entropy(
        rho, 1; base=base, atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

function _tierd_psd_sqrt(analysis; boundary_policy::Symbol=:reject)
    boundary_policy === :reject ||
        boundary_policy === :project ||
        throw(
            ArgumentError(
                "PSD square-root boundary policy must be :reject or :project; " *
                "got $boundary_policy",
            ),
        )
    if analysis.spectral_boundary_uncertain && boundary_policy === :reject
        throw(
            DomainError(
                analysis.minimum_eigenvalue,
                "the PSD square root encountered a negative eigenvalue inside " *
                "the validation tolerance; use an explicitly documented " *
                "boundary projection only in an operation that supports it",
            ),
        )
    end
    eigenvalues = if boundary_policy === :project
        map(value -> value < zero(value) ? zero(value) : value, analysis.eigenvalues)
    else
        analysis.eigenvalues
    end
    square_roots = sqrt.(eigenvalues)
    vectors = analysis.decomposition.vectors
    return vectors * Diagonal(square_roots) * adjoint(vectors)
end

"""
    fidelity(rho, sigma; squared=false, atol=nothing, rtol=nothing,
             allow_densify=false)

Return the Uhlmann root fidelity
`||sqrt(rho) * sqrt(sigma)||_1` by default, matching QETLAB's `Fidelity`
convention.  Set `squared=true` for its square.  Inputs are independently
validated density matrices and are never normalized, clipped, or mutated.
The result is not forcibly clipped to `[0, 1]`.  Inputs that pass approximate
Hermiticity validation use their explicitly formed Hermitian work matrices.
The eigendecompositions and SVD cost `O(n^3)` time and `O(n^2)` workspace.
"""
function fidelity(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    squared::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    size(rho) == size(sigma) || throw(
        DimensionMismatch(
            "rho and sigma must have the same size; got $(size(rho)) and $(size(sigma))"
        ),
    )
    rho_analysis = _tierd_density_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="fidelity"
    )
    sigma_analysis = _tierd_density_analysis(
        sigma; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="fidelity"
    )
    root_product = _tierd_psd_sqrt(rho_analysis) * _tierd_psd_sqrt(sigma_analysis)
    root_fidelity = sum(svdvals(root_product))
    return squared ? root_fidelity^2 : root_fidelity
end

function _tierd_matsumoto_support_policy(policy::Symbol)
    policy === :reject ||
        policy === :project ||
        throw(
            ArgumentError(
                "support_boundary_policy must be :reject or :project; got $policy"
            ),
        )
    return policy
end

function _tierd_exact_rank_one_projector_basis(analysis)
    matrix = analysis.matrix
    exact_projector =
        iszero(analysis.hermiticity_residual) &&
        iszero(analysis.trace_imaginary_residual) &&
        iszero(analysis.normalization_residual) &&
        all(iszero, matrix * matrix - matrix)
    exact_projector || return nothing
    diagonal_values = real.(diag(matrix))
    pivot = argmax(diagonal_values)
    pivot_value = diagonal_values[pivot]
    pivot_value > zero(pivot_value) || return nothing
    vector = copy(@view(matrix[:, pivot])) / sqrt(pivot_value)
    real_type = typeof(real(zero(eltype(matrix))))
    return (basis=reshape(vector, :, 1), eigenvalues=real_type[one(real_type)])
end

function _tierd_matsumoto_support(analysis, policy::Symbol, name::AbstractString)
    exact_projector = _tierd_exact_rank_one_projector_basis(analysis)
    exact_projector === nothing || return merge(exact_projector, (projected=false,))

    eigenvalues = analysis.eigenvalues
    negative = findall(value -> value < zero(value), eigenvalues)
    isempty(negative) || throw(
        DomainError(
            eigenvalues[negative],
            "$name has negative eigenvalues inside the density-validation " *
            "tolerance. Matsumoto fidelity does not repair PSD boundaries.",
        ),
    )
    positive_boundary = findall(
        value -> value > zero(value) && value <= analysis.tolerance, eigenvalues
    )
    if !isempty(positive_boundary) && policy === :reject
        throw(
            DomainError(
                eigenvalues[positive_boundary],
                "$name has positive eigenvalues inside the numerical support " *
                "boundary (tolerance $(analysis.tolerance)); " *
                "support_boundary_policy=:reject refuses to drop them",
            ),
        )
    end
    positive = findall(value -> value > analysis.tolerance, eigenvalues)
    isempty(positive) && throw(
        DomainError(
            eigenvalues,
            "$name has no positive spectral support outside the numerical boundary",
        ),
    )
    return (
        basis=analysis.decomposition.vectors[:, positive],
        eigenvalues=eigenvalues[positive],
        projected=(!isempty(positive_boundary)),
    )
end

function _tierd_matsumoto_intersection_bases(
    left_basis::AbstractMatrix, right_basis::AbstractMatrix, tolerance, policy::Symbol
)
    dimension = size(left_basis, 1)
    size(right_basis, 1) == dimension ||
        throw(DimensionMismatch("support bases must have the same ambient dimension"))
    left_rank = size(left_basis, 2)
    right_rank = size(right_basis, 2)
    if left_rank == dimension
        return right_basis, right_basis
    elseif right_rank == dimension
        return left_basis, left_basis
    end

    left_projector = left_basis * adjoint(left_basis)
    right_projector = right_basis * adjoint(right_basis)
    if left_projector == right_projector
        return left_basis, left_basis
    elseif right_projector * left_basis == left_basis
        return left_basis, left_basis
    elseif left_projector * right_basis == right_basis
        return right_basis, right_basis
    end

    principal = svd(adjoint(left_basis) * right_basis; full=false)
    singular_values = principal.S
    real_type = eltype(singular_values)
    roundoff_tolerance =
        8 * dimension * eps(real_type) * max(one(real_type), maximum(singular_values))
    intersection_tolerance = max(convert(real_type, tolerance), roundoff_tolerance)
    any(value -> value > one(value) + intersection_tolerance, singular_values) &&
        error("principal support cosine exceeded one outside derived roundoff tolerance")

    guaranteed_dimension = max(left_rank + right_rank - dimension, 0)
    selected = falses(length(singular_values))
    for index in 1:guaranteed_dimension
        distance = abs(one(singular_values[index]) - singular_values[index])
        distance <= intersection_tolerance || error(
            "support ranks imply an intersection of dimension " *
            "$guaranteed_dimension, but its principal cosine residual $distance " *
            "exceeds tolerance $intersection_tolerance",
        )
        selected[index] = true
    end

    boundary_indices = Int[]
    for index in (guaranteed_dimension + 1):length(singular_values)
        value = singular_values[index]
        distance = abs(one(value) - value)
        if iszero(distance)
            selected[index] = true
        elseif distance <= intersection_tolerance
            push!(boundary_indices, index)
        end
    end
    if !isempty(boundary_indices)
        policy === :project || throw(
            DomainError(
                singular_values[boundary_indices],
                "the two spectral supports have near-unit principal cosines " *
                "inside tolerance $intersection_tolerance; the exact " *
                "intersection dimension is numerically ambiguous",
            ),
        )
        selected[boundary_indices] .= true
    end

    indices = findall(selected)
    isempty(indices) && return (
        zeros(eltype(left_basis), dimension, 0),
        zeros(eltype(right_basis), dimension, 0),
    )
    left_intersection = left_basis * principal.U[:, indices]
    right_intersection = right_basis * principal.V[:, indices]
    return left_intersection, right_intersection
end

function _tierd_matsumoto_shorted_coordinates(
    support_basis::AbstractMatrix,
    support_eigenvalues::AbstractVector,
    intersection_basis::AbstractMatrix,
)
    coordinates = adjoint(support_basis) * intersection_basis
    inverse_square_roots = map(value -> one(value) / sqrt(value), support_eigenvalues)
    weighted = Diagonal(inverse_square_roots) * coordinates
    compressed_support_inverse = adjoint(weighted) * weighted
    factorization = cholesky(Hermitian(compressed_support_inverse); check=true)
    identity_matrix = Matrix{eltype(compressed_support_inverse)}(
        I, size(compressed_support_inverse)...
    )
    shorted = factorization \ identity_matrix
    return (shorted + adjoint(shorted)) / 2
end

function _tierd_positive_definite_geometric_mean_trace(
    left::AbstractMatrix, right::AbstractMatrix
)
    factorization = cholesky(Hermitian(left); check=true)
    factor = factorization.L
    relative = (factor \ right) / adjoint(factor)
    relative_work = (relative + adjoint(relative)) / 2
    decomposition = eigen(Hermitian(relative_work))
    minimum_eigenvalue = minimum(decomposition.values)
    minimum_eigenvalue > zero(minimum_eigenvalue) || throw(
        DomainError(
            minimum_eigenvalue,
            "the support-restricted relative operator is not positive definite; " *
            "no eigenvalue is clipped or regularized",
        ),
    )
    relative_root =
        decomposition.vectors *
        Diagonal(sqrt.(decomposition.values)) *
        adjoint(decomposition.vectors)
    geometric_mean = factor * relative_root * adjoint(factor)
    trace_value = tr(geometric_mean)
    real_type = typeof(real(zero(eltype(geometric_mean))))
    imaginary_residual = abs(imag(trace_value))
    scale = max(one(real_type), abs(real(trace_value)))
    tolerance = 32 * size(left, 1) * eps(real_type) * scale
    imaginary_residual <= tolerance || error(
        "matrix geometric-mean trace has imaginary residual " *
        "$imaginary_residual above derived roundoff tolerance $tolerance",
    )
    result = real(trace_value)
    result >= zero(result) || throw(
        DomainError(
            result,
            "computed Matsumoto fidelity is negative; no range clipping is performed",
        ),
    )
    return result
end

"""
    matsumoto_fidelity(
        rho, sigma; atol=nothing, rtol=nothing, allow_densify=false,
        support_boundary_policy=:reject
    )

Return the Matsumoto fidelity `tr(rho # sigma)`, where `#` is the
Kubo--Ando matrix geometric mean. This is distinct from [`fidelity`](@ref),
which returns Uhlmann root fidelity.

Positive-definite inputs use Cholesky solves and one Hermitian
eigendecomposition; no explicit matrix inverse or additive regularizer is
formed. For singular inputs, the implementation finds the intersection of
the proved spectral supports, constructs each state's shorted operator on
that intersection through support-restricted solves, and evaluates their
positive-definite geometric mean. Consequently, distinct pure states have
zero Matsumoto fidelity even when their overlap is nonzero.

Inputs must be finite, square, exactly Hermitian, positive semidefinite, and
normalized within `atol` and `rtol`; they are never normalized or
symmetrized. Small negative eigenvalues are rejected even when they lie inside
the density-validation tolerance. A positive eigenvalue in `(0, tolerance]`
or near-coincident singular supports make the exact support numerically
ambiguous. `support_boundary_policy=:reject` reports that boundary.
`:project` explicitly drops the bounded positive spectral tail and permits
canonical identification of near-coincident principal support directions.
No support is dropped silently.

General sparse matrices require `allow_densify=true`. The dense spectral path
supports standard-library BLAS floating types without precision conversion
and costs `O(n^3)` time and `O(n^2)` workspace. Two `Diagonal`
density matrices use the commuting formula `sum(sqrt.(p) .* sqrt.(q))`
directly in `O(n)` time, preserving generic floating types such as
`BigFloat`.

# Examples

```jldoctest
julia> using LinearAlgebra

julia> matsumoto_fidelity(
           Diagonal([0.7, 0.2, 0.1]),
           Diagonal([0.4, 0.5, 0.1]),
       )
0.945378028229756
```
"""
function matsumoto_fidelity(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    support_boundary_policy::Symbol=:reject,
)
    size(rho) == size(sigma) || throw(
        DimensionMismatch(
            "rho and sigma must have the same size; got $(size(rho)) and $(size(sigma))"
        ),
    )
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("rho and sigma must be square; got size $(size(rho))"))
    policy = _tierd_matsumoto_support_policy(support_boundary_policy)
    ishermitian(rho) || throw(
        ArgumentError(
            "rho must be exactly Hermitian; Matsumoto fidelity never symmetrizes input"
        ),
    )
    ishermitian(sigma) || throw(
        ArgumentError(
            "sigma must be exactly Hermitian; Matsumoto fidelity never symmetrizes input",
        ),
    )
    rho_analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="matsumoto_fidelity",
    )
    sigma_analysis = _tierd_density_analysis(
        sigma;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="matsumoto_fidelity",
    )
    rho_support = _tierd_matsumoto_support(rho_analysis, policy, "rho")
    sigma_support = _tierd_matsumoto_support(sigma_analysis, policy, "sigma")
    dimension = size(rho, 1)

    if size(rho_support.basis, 2) == dimension && size(sigma_support.basis, 2) == dimension
        return _tierd_positive_definite_geometric_mean_trace(
            rho_analysis.matrix, sigma_analysis.matrix
        )
    end

    support_tolerance = max(rho_analysis.tolerance, sigma_analysis.tolerance)
    rho_intersection, sigma_intersection = _tierd_matsumoto_intersection_bases(
        rho_support.basis, sigma_support.basis, support_tolerance, policy
    )
    intersection_dimension = size(rho_intersection, 2)
    real_type = promote_type(
        eltype(rho_support.eigenvalues), eltype(sigma_support.eigenvalues)
    )
    iszero(intersection_dimension) && return zero(real_type)
    rho_shorted = _tierd_matsumoto_shorted_coordinates(
        rho_support.basis, rho_support.eigenvalues, rho_intersection
    )
    sigma_shorted = _tierd_matsumoto_shorted_coordinates(
        sigma_support.basis, sigma_support.eigenvalues, sigma_intersection
    )
    return _tierd_positive_definite_geometric_mean_trace(rho_shorted, sigma_shorted)
end

function matsumoto_fidelity(
    rho::Diagonal{<:Number},
    sigma::Diagonal{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    support_boundary_policy::Symbol=:reject,
)
    size(rho) == size(sigma) || throw(
        DimensionMismatch(
            "rho and sigma must have the same size; got $(size(rho)) and $(size(sigma))"
        ),
    )
    _tierd_matsumoto_support_policy(support_boundary_policy)
    rho_values = _tierd_diagonal_density_eigenvalues(
        rho; atol=atol, rtol=rtol, operation="matsumoto_fidelity"
    )
    sigma_values = _tierd_diagonal_density_eigenvalues(
        sigma; atol=atol, rtol=rtol, operation="matsumoto_fidelity"
    )
    working_type = promote_type(eltype(rho_values), eltype(sigma_values))
    result = zero(working_type)
    for (rho_value, sigma_value) in zip(rho_values, sigma_values)
        result +=
            sqrt(convert(working_type, rho_value)) *
            sqrt(convert(working_type, sigma_value))
    end
    return result
end

"""
    trace_distance(rho, sigma; atol=nothing, rtol=nothing,
                   allow_densify=false)

Return `trace_norm(rho - sigma) / 2` for two validated density matrices.
No normalization or range clipping is performed.  The validation
eigendecompositions and trace-norm SVD cost `O(n^3)` time and `O(n^2)`
workspace.
"""
function trace_distance(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    size(rho) == size(sigma) || throw(
        DimensionMismatch(
            "rho and sigma must have the same size; got $(size(rho)) and $(size(sigma))"
        ),
    )
    rho_analysis = _tierd_density_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="trace_distance"
    )
    sigma_analysis = _tierd_density_analysis(
        sigma; atol=atol, rtol=rtol, allow_densify=allow_densify, operation="trace_distance"
    )
    return sum(svdvals(rho_analysis.matrix - sigma_analysis.matrix)) / 2
end

function _tierd_density_from_state(
    state::AbstractVector{<:Number}; atol, rtol, operation::AbstractString
)
    _tierd_validate_pure_state(state; atol=atol, rtol=rtol, operation=operation)
    return state * adjoint(state)
end

function _tierd_density_from_state(
    state::AbstractMatrix{<:Number}; atol, rtol, operation::AbstractString
)
    return state
end

function _tierd_negativity_spectrum(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dims;
    systems,
    atol,
    rtol,
    allow_densify::Bool,
    operation::AbstractString,
)
    rho = _tierd_density_from_state(state; atol=atol, rtol=rtol, operation=operation)
    analysis = _tierd_density_analysis(
        rho; atol=atol, rtol=rtol, allow_densify=allow_densify, operation=operation
    )
    layout = _as_layout(dims)
    length(layout) >= 2 ||
        throw(ArgumentError("$operation requires at least two subsystems"))
    layout.total_dimension == size(analysis.matrix, 1) || throw(
        DimensionMismatch(
            "prod(dims)=$(layout.total_dimension) does not match state " *
            "dimension $(size(analysis.matrix, 1))",
        ),
    )
    selected = _normalize_systems(
        systems,
        length(layout);
        name="systems",
        allow_empty=false,
        proper=true,
        sort_result=true,
    )
    transpose_plan = PartialTransposePlan(layout, selected)
    transposed = partial_transpose(analysis.matrix, transpose_plan)
    hermitian_transpose = Hermitian((transposed + adjoint(transposed)) / 2)
    return eigvals(hermitian_transpose)
end

"""
    negativity(state, dims; systems=(2,), atol=nothing, rtol=nothing,
               allow_densify=false)

Return the sum of absolute negative eigenvalues of the selected partial
transpose.  `state` may be a normalized pure vector or a density matrix.
`dims` fixes subsystem order and `systems` defines one nonempty, proper side
of the bipartition.  The state is validated but never normalized or clipped.
The dense eigendecompositions cost `O(n^3)` time and `O(n^2)` workspace.
"""
function negativity(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dims;
    systems=(2,),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    eigenvalues = _tierd_negativity_spectrum(
        state,
        dims;
        systems=systems,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="negativity",
    )
    return sum(value -> value < zero(value) ? -value : zero(value), eigenvalues)
end

"""
    logarithmic_negativity(state, dims; systems=(2,), base,
                           atol=nothing, rtol=nothing,
                           allow_densify=false)

Return `log(sum(abs, eigvals(partial_transpose(state)))) / log(base)`.
The logarithm base is required explicitly and must be greater than one.  The
state is validated but no eigenvalue or result is clipped.  Complexity is
dominated by dense `O(n^3)` eigendecompositions.
"""
function logarithmic_negativity(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    dims;
    systems=(2,),
    base,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    checked_base = _tierd_validate_log_base(base)
    eigenvalues = _tierd_negativity_spectrum(
        state,
        dims;
        systems=systems,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="logarithmic_negativity",
    )
    return log(sum(abs, eigenvalues)) / log(checked_base)
end

"""
    schmidt_decomposition(psi, dims; allow_densify=false)

Compute the thin Schmidt decomposition of a finite bipartite vector.
`dims == (dA, dB)` uses subsystem 1 as the slowest-varying tensor factor.
The vector need not be normalized.  All `min(dA, dB)` coefficients are
returned in descending order.  Sparse input requires explicit
`allow_densify=true`.  The full SVD costs
`O(dA*dB*min(dA,dB))` time and `O(dA*dB)` workspace.
"""
function schmidt_decomposition(
    psi::AbstractVector{<:Number}, dims; allow_densify::Bool=false
)
    layout = _tierd_bipartite_layout(dims, length(psi))
    dense = _tierd_dense_vector(
        psi; allow_densify=allow_densify, operation="schmidt_decomposition"
    )
    coefficient_matrix = reshape(dense, layout.dims[2], layout.dims[1])
    factorization = svd(coefficient_matrix; full=false)
    return SchmidtDecompositionResult(
        factorization.S, conj(factorization.V), factorization.U
    )
end

"""
    schmidt_coefficients(psi, dims; allow_densify=false)

Return all Schmidt coefficients, including zeros, in descending order.
"""
function schmidt_coefficients(
    psi::AbstractVector{<:Number}, dims; allow_densify::Bool=false
)
    return schmidt_decomposition(psi, dims; allow_densify=allow_densify).coefficients
end

"""
    schmidt_k_norm(psi, dims, k=1; allow_densify=false)

Return the Euclidean norm of the `k` largest Schmidt coefficients of a finite
bipartite vector. `dims` describes exactly two subsystems in the same
slowest-factor-first order as [`schmidt_decomposition`](@ref). `k` must be
positive and is clipped to `min(dims...)`, matching the reviewed numeric
behavior of QETLAB `SkVectorNorm`. The vector need not be normalized.

When `k >= min(dims...)`, the result is `norm(psi)` and sparse or generic
floating-point vectors, including `BigFloat`, are handled without
densification. A smaller `k` requires a full SVD and therefore requires
`allow_densify=true` for sparse input and a BLAS floating element type in the
dependency-free core. The SVD path costs
`O(dA*dB*min(dA,dB))` time and `O(dA*dB)` workspace.

# Examples

```jldoctest
julia> schmidt_k_norm([3.0, 0, 0, 0, 0, 4], (2, 3), 1)
4.0
```
"""
function schmidt_k_norm(psi::AbstractVector{<:Number}, dims, k=1; allow_densify::Bool=false)
    layout = _tierd_bipartite_layout(dims, length(psi))
    _tierd_require_finite(psi, "psi")
    maximum_k = min(layout.dims...)
    effective_k = _tierd_validate_top_k(k, maximum_k)
    effective_k == maximum_k && return norm(psi)
    coefficients = schmidt_coefficients(psi, layout; allow_densify=allow_densify)
    return norm(@view coefficients[1:effective_k])
end

"""
    schmidt_rank(psi, dims; atol=nothing, rtol=nothing,
                 allow_densify=false)

Count Schmidt coefficients strictly larger than
`atol + rtol * maximum(coefficients)`.  Defaults are zero absolute tolerance
and `sqrt(eps(R))` relative tolerance.  This function returns the
tolerance-defined numerical rank; it does not claim an exact algebraic rank.
"""
function schmidt_rank(
    psi::AbstractVector{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    coefficients = schmidt_coefficients(psi, dims; allow_densify=allow_densify)
    real_type = eltype(coefficients)
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    scale = maximum(coefficients; init=zero(real_type))
    threshold = absolute + relative * scale
    return count(value -> value > threshold, coefficients)
end

"""
    concurrence(psi; atol=nothing, rtol=nothing)
    concurrence(
        rho; atol=nothing, rtol=nothing, allow_densify=false,
        psd_boundary_policy=:project, range_boundary_policy=:project
    )

Return Wootters concurrence for a normalized two-qubit pure vector or
`4 x 4` density matrix.  Pure vectors use `2abs(ad-bc)`.  Mixed states use
the singular-value form of the spin-flip construction.  The final
`max(s1-s2-s3-s4, 0)` is part of the mathematical definition; density
eigenvalues and returned values are otherwise never clipped. Rank-deficient
density matrices can acquire negative eigenvalues inside the requested
validation tolerance from roundoff. The explicit default
`psd_boundary_policy=:project` replaces only those already-bounded negative
eigenvalues by zero for the PSD square root; use `:reject` to refuse this
projection. Eigenvalues below the validation tolerance always raise.
Likewise, `range_boundary_policy=:project` maps a computed value in
`(1, 1 + tolerance]` to the exact upper endpoint; `:reject` refuses that
bounded repair, and larger values always raise.

Sparse pure vectors are evaluated directly without densification.  The
`allow_densify` keyword applies to the density-matrix method.  The pure-state
formula is constant-time after validation; the mixed-state method costs
`O(n^3)` time and `O(n^2)` workspace (with `n == 4` here).
"""
function concurrence(state::AbstractVector{<:Number}; atol=nothing, rtol=nothing)
    length(state) == 4 || throw(
        DimensionMismatch(
            "two-qubit concurrence requires a length-4 vector; got length $(length(state))",
        ),
    )
    _tierd_validate_pure_state(state; atol=atol, rtol=rtol, operation="concurrence")
    return 2 * abs(state[1] * state[4] - state[2] * state[3])
end

function concurrence(
    rho::AbstractMatrix{<:Number};
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    psd_boundary_policy::Symbol=:project,
    range_boundary_policy::Symbol=:project,
)
    size(rho) == (4, 4) || throw(
        DimensionMismatch(
            "two-qubit concurrence requires a 4 x 4 density matrix; got $(size(rho))"
        ),
    )
    analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="concurrence",
    )
    real_type = eltype(analysis.eigenvalues)
    imaginary_unit = complex(zero(real_type), one(real_type))
    pauli_y = Complex{real_type}[
        zero(real_type) -imaginary_unit
        imaginary_unit zero(real_type)
    ]
    spin_flip = kron(pauli_y, pauli_y)
    square_root = _tierd_psd_sqrt(analysis; boundary_policy=psd_boundary_policy)
    spin_amplitude = square_root * spin_flip * conj(square_root)
    singular_values = sort(svdvals(spin_amplitude); rev=true)
    candidate =
        singular_values[1] - singular_values[2] - singular_values[3] - singular_values[4]
    result = max(candidate, zero(candidate))
    if result > one(result)
        range_tolerance = _tierd_threshold(one(result), analysis.atol, analysis.rtol)
        result <= one(result) + range_tolerance || throw(
            DomainError(
                result,
                "computed concurrence exceeds one by more than the requested " *
                "numerical tolerance",
            ),
        )
        range_boundary_policy === :project ||
            range_boundary_policy === :reject ||
            throw(
                ArgumentError(
                    "concurrence range boundary policy must be :project or " *
                    ":reject; got $range_boundary_policy",
                ),
            )
        range_boundary_policy === :reject && throw(
            DomainError(
                result,
                "computed concurrence lies just above one inside the numerical " *
                "tolerance and range_boundary_policy=:reject refuses projection",
            ),
        )
        return one(result)
    end
    range_boundary_policy === :project ||
        range_boundary_policy === :reject ||
        throw(
            ArgumentError(
                "concurrence range boundary policy must be :project or :reject; " *
                "got $range_boundary_policy",
            ),
        )
    return result
end
