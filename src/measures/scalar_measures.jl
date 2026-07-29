# Source-informed independent Julia implementations based on the specifications
# in QETLAB TraceNorm.m, SchattenNorm.m, KyFanNorm.m, Purity.m, Entropy.m,
# Fidelity.m, Negativity.m, SchmidtDecomposition.m, SchmidtRank.m, and
# Concurrence.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and named coauthors,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

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

"""
    von_neumann_entropy(rho; base, atol=nothing, rtol=nothing,
                        allow_densify=false)

Return `-tr(rho * log(rho)) / log(base)` for a validated density matrix.
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
    checked_base = _tierd_validate_log_base(base)
    analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        operation="von_neumann_entropy",
    )
    numerator = zero(eltype(analysis.eigenvalues))
    for value in analysis.eigenvalues
        iszero(value) && continue
        numerator -= value * log(value)
    end
    return numerator / log(checked_base)
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
