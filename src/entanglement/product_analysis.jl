# Source-informed independent Julia implementations based on the specifications
# in QETLAB OperatorSchmidtDecomposition.m, OperatorSchmidtRank.m,
# IsProductVector.m, IsProductOperator.m, EntFormation.m, and
# InSeparableBall.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    OperatorSchmidtDecompositionResult

Operator Schmidt decomposition returned by
[`operator_schmidt_decomposition`](@ref). If `result` decomposes `X`, then

```julia
X ≈ tensor_sum(
    result.left_factors,
    result.right_factors;
    weights = result.coefficients,
)
```

The factors are Frobenius-orthonormal within each side. Singular-vector phases
are not unique, and Hermitian input does not imply that each returned factor is
Hermitian. `row_dims` and `column_dims` record the two local row and column
dimensions used for the decomposition.
"""
struct OperatorSchmidtDecompositionResult{
    C<:AbstractVector,L<:AbstractVector,R<:AbstractVector,RD<:Tuple,CD<:Tuple
}
    coefficients::C
    left_factors::L
    right_factors::R
    row_dims::RD
    column_dims::CD
end

"""
    ProductAnalysisResult

Tolerance-aware result returned by [`is_product_vector`](@ref) and
[`is_product_operator`](@ref).

`status` is one of:

- `:within_tolerance`: every recursive Schmidt-tail residual and the final
  product-approximation residual are below their thresholds;
- `:outside_tolerance`: at least one residual exceeds its threshold;
- `:boundary`: no residual exceeds its threshold and at least one positive
  residual equals it.

This is a numerical classification, not an exact symbolic proof. `factors`
contains the sequential leading-Schmidt product approximation in subsystem
order, even when the status is not `:within_tolerance`.
`reconstruction_residual` is its Euclidean or Frobenius residual in the
original input's units, and `approximation_threshold` is the corresponding
global threshold. `cut_residuals` are full Schmidt-tail norms scaled back into
the original input's units; `cut_thresholds` expose their numerical decisions.
"""
struct ProductAnalysisResult{F,R,AT,CR<:AbstractVector,CT<:AbstractVector}
    status::Symbol
    factors::F
    reconstruction_residual::R
    approximation_threshold::AT
    cut_residuals::CR
    cut_thresholds::CT
    message::String
end

"""
    SeparableBallResult

Structured Gurvits--Barnum separable-ball result returned by
[`in_separable_ball`](@ref).

`status == :separable_certified` is a sufficient separability certificate.
`:outside_ball` means only that this particular sufficient condition did not
certify the input; it is never an entanglement conclusion. `:unknown` denotes
a numerical boundary in either input PSD validation or the ball inequality.
The reported `purity`, exact ball `boundary`, and numerical `tolerance` make
the decision auditable.
"""
struct SeparableBallResult{P,B,T}
    status::Symbol
    purity::P
    boundary::B
    tolerance::T
    message::String
end

function Base.show(io::IO, result::OperatorSchmidtDecompositionResult)
    return print(
        io,
        "OperatorSchmidtDecompositionResult(",
        length(result.coefficients),
        " terms, row_dims=",
        result.row_dims,
        ", column_dims=",
        result.column_dims,
        ")",
    )
end

function Base.show(io::IO, result::ProductAnalysisResult)
    return print(
        io,
        "ProductAnalysisResult(status=",
        result.status,
        ", residual=",
        result.reconstruction_residual,
        ", threshold=",
        result.approximation_threshold,
        ", cuts=",
        length(result.cut_residuals),
        ")",
    )
end

function Base.show(io::IO, result::SeparableBallResult)
    return print(
        io,
        "SeparableBallResult(status=",
        result.status,
        ", purity=",
        result.purity,
        ", boundary=",
        result.boundary,
        ", tolerance=",
        result.tolerance,
        ")",
    )
end

function _tiere_operator_layouts(
    operator::AbstractMatrix,
    row_dims,
    column_dims;
    bipartite::Bool,
    operation::AbstractString,
)
    row_layout = _as_layout(row_dims)
    column_layout = _as_layout(column_dims)
    length(row_layout) == length(column_layout) || throw(
        DimensionMismatch(
            "row dims has $(length(row_layout)) subsystems but column dims " *
            "has $(length(column_layout))",
        ),
    )
    minimum_count = 2
    length(row_layout) >= minimum_count ||
        throw(ArgumentError("$operation requires at least two subsystems"))
    if bipartite && length(row_layout) != 2
        throw(
            ArgumentError(
                "$operation requires exactly two subsystems; got " * "$(row_layout.dims)"
            ),
        )
    end
    _validate_matrix_dimensions(operator, row_layout, column_layout)
    return row_layout, column_layout
end

function _tiere_operator_schmidt_state(
    operator::AbstractMatrix{<:Number},
    row_dims,
    column_dims;
    allow_densify::Bool,
    operation::AbstractString,
)
    row_layout, column_layout = _tiere_operator_layouts(
        operator, row_dims, column_dims; bipartite=true, operation=operation
    )
    dense, _ = _tierd_dense_matrix(
        operator; allow_densify=allow_densify, operation=operation
    )
    plan = RealignmentPlan(row_layout, column_layout; systems=(1,))
    realigned = realign(dense, plan)
    # `realign` orders each local matrix in row-major pair order. Transposing
    # without conjugation before `vec` converts that matrix into the
    # slowest-subsystem-first state convention used by `schmidt_decomposition`.
    state = vec(permutedims(realigned, (2, 1)))
    return (
        state=state,
        state_dims=plan.output_size,
        row_layout=row_layout,
        column_layout=column_layout,
    )
end

function _tiere_unrowmajor(
    vector::AbstractVector, row_dimension::Int, column_dimension::Int
)
    return permutedims(reshape(copy(vector), column_dimension, row_dimension), (2, 1))
end

"""
    operator_schmidt_decomposition(operator, dims; allow_densify=false)
    operator_schmidt_decomposition(
        operator, row_dims, column_dims; allow_densify=false
    )

Compute the full operator Schmidt decomposition of a bipartite operator.
`row_dims == (rA, rB)` and `column_dims == (cA, cB)` describe an
`(rA*rB) x (cA*cB)` matrix, so independently rectangular local operator
spaces are supported. The one-layout method is the square-operator
convenience form.

The implementation reuses [`realign`](@ref) and
[`schmidt_decomposition`](@ref). It returns all
`min(rA*cA, rB*cB)` coefficients, including numerical zeros, in descending
order. The full dense SVD costs
`O(rA*cA*rB*cB*min(rA*cA,rB*cB))` time and dense workspace proportional to
the operator size. Sparse input is rejected unless `allow_densify=true`.
Only BLAS floating element types are supported; no precision-changing
conversion, normalization, clipping, or Hermitian-factor repair is performed.
"""
function operator_schmidt_decomposition(
    operator::AbstractMatrix{<:Number}, dims; allow_densify::Bool=false
)
    return operator_schmidt_decomposition(operator, dims, dims; allow_densify=allow_densify)
end

function operator_schmidt_decomposition(
    operator::AbstractMatrix{<:Number}, row_dims, column_dims; allow_densify::Bool=false
)
    data = _tiere_operator_schmidt_state(
        operator,
        row_dims,
        column_dims;
        allow_densify=allow_densify,
        operation="operator_schmidt_decomposition",
    )
    decomposition = schmidt_decomposition(data.state, data.state_dims; allow_densify=false)
    term_count = length(decomposition.coefficients)
    left_factors = [
        _tiere_unrowmajor(
            decomposition.left_vectors[:, index], data.row_layout[1], data.column_layout[1]
        ) for index in 1:term_count
    ]
    right_factors = [
        _tiere_unrowmajor(
            decomposition.right_vectors[:, index], data.row_layout[2], data.column_layout[2]
        ) for index in 1:term_count
    ]
    return OperatorSchmidtDecompositionResult(
        decomposition.coefficients,
        left_factors,
        right_factors,
        data.row_layout.dims,
        data.column_layout.dims,
    )
end

"""
    operator_schmidt_coefficients(operator, dims; allow_densify=false)
    operator_schmidt_coefficients(
        operator, row_dims, column_dims; allow_densify=false
    )

Return all operator Schmidt coefficients, including numerical zeros, in
descending order. Dimension, sparse-densification, element-type, and
complexity behavior match [`operator_schmidt_decomposition`](@ref).
"""
function operator_schmidt_coefficients(
    operator::AbstractMatrix{<:Number}, dims; allow_densify::Bool=false
)
    return operator_schmidt_coefficients(operator, dims, dims; allow_densify=allow_densify)
end

function operator_schmidt_coefficients(
    operator::AbstractMatrix{<:Number}, row_dims, column_dims; allow_densify::Bool=false
)
    data = _tiere_operator_schmidt_state(
        operator,
        row_dims,
        column_dims;
        allow_densify=allow_densify,
        operation="operator_schmidt_coefficients",
    )
    return schmidt_coefficients(data.state, data.state_dims; allow_densify=false)
end

"""
    operator_schmidt_rank(
        operator, dims; atol=nothing, rtol=nothing, allow_densify=false
    )
    operator_schmidt_rank(
        operator, row_dims, column_dims;
        atol=nothing, rtol=nothing, allow_densify=false
    )

Return the tolerance-defined numerical operator Schmidt rank. Coefficients are
counted when they are strictly larger than
`atol + rtol * maximum(coefficients)`. Defaults are zero absolute tolerance
and `sqrt(eps(R))` relative tolerance. This is not an exact symbolic rank.
"""
function operator_schmidt_rank(
    operator::AbstractMatrix{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    return operator_schmidt_rank(
        operator, dims, dims; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

function operator_schmidt_rank(
    operator::AbstractMatrix{<:Number},
    row_dims,
    column_dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    data = _tiere_operator_schmidt_state(
        operator,
        row_dims,
        column_dims;
        allow_densify=allow_densify,
        operation="operator_schmidt_rank",
    )
    return schmidt_rank(
        data.state, data.state_dims; atol=atol, rtol=rtol, allow_densify=false
    )
end

function _tiere_update_product_status(status::Symbol, residual, threshold)
    residual > threshold && return :outside_tolerance
    if status !== :outside_tolerance && !iszero(residual) && residual == threshold
        return :boundary
    end
    return status
end

function _tiere_product_message(status::Symbol, noun::AbstractString)
    status === :within_tolerance &&
        return "$noun is product-approximable within every requested residual tolerance"
    status === :outside_tolerance &&
        return "$noun has a Schmidt-tail or reconstruction residual above tolerance"
    return "$noun lies on a requested residual-tolerance boundary"
end

function _tiere_factor_tuple(factors::Vector{T}, ::SubsystemLayout{N}) where {T,N}
    return ntuple(index -> factors[index], Val(N))
end

"""
    is_product_vector(
        vector, dims; atol=nothing, rtol=nothing, allow_densify=false
    ) -> ProductAnalysisResult

Analyze whether a nonzero vector is a multipartite elementary tensor. `dims`
contains at least two subsystem dimensions in tensor-product order. The method
recursively reuses [`schmidt_decomposition`](@ref), comparing the full
nonleading Schmidt-tail norm of every first-subsystem-versus-remainder cut
with `atol + rtol*s₁`. Residuals and leading coefficients are scaled back to
the original vector's units at later recursive cuts. The final reconstructed
product is also checked against `atol + rtol*norm(vector)`. Defaults are zero
absolute tolerance and `sqrt(eps(R))` relative tolerance.

The vector need not be normalized, but the zero vector is rejected because it
has no meaningful state-like product decomposition. Sparse input requires
`allow_densify=true`; only BLAS floating element types are supported. No input
is normalized, clipped, or mutated. The sequence of dense SVDs is polynomial
in the vector length and stores dense intermediate factors.
"""
function is_product_vector(
    vector::AbstractVector{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    layout = _as_layout(dims)
    length(layout) >= 2 ||
        throw(ArgumentError("is_product_vector requires at least two subsystems"))
    _validate_vector_dimension(vector, layout)
    dense = _tierd_dense_vector(
        vector; allow_densify=allow_densify, operation="is_product_vector"
    )
    input_norm = norm(dense)
    iszero(input_norm) &&
        throw(DomainError(input_norm, "is_product_vector requires a nonzero vector"))
    real_type = _tierd_real_type(eltype(dense))
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    threshold_type = promote_type(real_type, typeof(absolute), typeof(relative))
    cut_residuals = Vector{threshold_type}()
    cut_thresholds = Vector{threshold_type}()
    factors = Vector{Vector{eltype(dense)}}()
    current = dense
    leading_scale = one(real_type)
    status = :within_tolerance

    for subsystem in 1:(length(layout) - 1)
        remainder_dimension = _checked_product(
            layout.dims[(subsystem + 1):end], "remaining product-vector dimensions"
        )
        decomposition = schmidt_decomposition(
            current, (layout[subsystem], remainder_dimension); allow_densify=false
        )
        coefficients = decomposition.coefficients
        leading = coefficients[1]
        residual =
            leading_scale *
            (length(coefficients) == 1 ? zero(leading) : norm(@view coefficients[2:end]))
        threshold = absolute + relative * (leading_scale * leading)
        push!(cut_residuals, residual)
        push!(cut_thresholds, threshold)
        status = _tiere_update_product_status(status, residual, threshold)
        push!(factors, Vector(decomposition.left_vectors[:, 1]))
        current = Vector(decomposition.right_vectors[:, 1])
        leading_scale *= leading
    end
    push!(factors, current)
    factors[1] .*= leading_scale
    reconstruction = tensor_product(factors...)
    reconstruction_residual = norm(dense - reconstruction)
    approximation_threshold = absolute + relative * input_norm
    status = _tiere_update_product_status(
        status, reconstruction_residual, approximation_threshold
    )
    return ProductAnalysisResult(
        status,
        _tiere_factor_tuple(factors, layout),
        reconstruction_residual,
        approximation_threshold,
        cut_residuals,
        cut_thresholds,
        _tiere_product_message(status, "vector"),
    )
end

"""
    is_product_operator(
        operator, dims; atol=nothing, rtol=nothing, allow_densify=false
    )
    is_product_operator(
        operator, row_dims, column_dims;
        atol=nothing, rtol=nothing, allow_densify=false
    ) -> ProductAnalysisResult

Analyze whether a nonzero multipartite operator is an elementary tensor.
Separate row and column layouts support independently rectangular local
operators. The method recursively applies
[`operator_schmidt_decomposition`](@ref) to the first local operator versus
the remainder. Every cut records the full nonleading Schmidt-tail norm in the
original operator's units, and the reconstructed product is checked against
`atol + rtol*norm(operator)`.

The zero operator is rejected because its factorization is nonunique and
cannot carry meaningful factors. Sparse input requires explicit densification,
and only BLAS floating element types are supported. The input is never
normalized, clipped, symmetrized, or mutated.
"""
function is_product_operator(
    operator::AbstractMatrix{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    return is_product_operator(
        operator, dims, dims; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
end

function is_product_operator(
    operator::AbstractMatrix{<:Number},
    row_dims,
    column_dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    row_layout, column_layout = _tiere_operator_layouts(
        operator, row_dims, column_dims; bipartite=false, operation="is_product_operator"
    )
    dense, real_type = _tierd_dense_matrix(
        operator; allow_densify=allow_densify, operation="is_product_operator"
    )
    input_norm = norm(dense)
    iszero(input_norm) &&
        throw(DomainError(input_norm, "is_product_operator requires a nonzero operator"))
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    threshold_type = promote_type(real_type, typeof(absolute), typeof(relative))
    cut_residuals = Vector{threshold_type}()
    cut_thresholds = Vector{threshold_type}()
    factors = Vector{Matrix{eltype(dense)}}()
    current = dense
    leading_scale = one(real_type)
    status = :within_tolerance

    for subsystem in 1:(length(row_layout) - 1)
        remaining_rows = _checked_product(
            row_layout.dims[(subsystem + 1):end],
            "remaining product-operator row dimensions",
        )
        remaining_columns = _checked_product(
            column_layout.dims[(subsystem + 1):end],
            "remaining product-operator column dimensions",
        )
        decomposition = operator_schmidt_decomposition(
            current,
            (row_layout[subsystem], remaining_rows),
            (column_layout[subsystem], remaining_columns);
            allow_densify=false,
        )
        coefficients = decomposition.coefficients
        leading = coefficients[1]
        residual =
            leading_scale *
            (length(coefficients) == 1 ? zero(leading) : norm(@view coefficients[2:end]))
        threshold = absolute + relative * (leading_scale * leading)
        push!(cut_residuals, residual)
        push!(cut_thresholds, threshold)
        status = _tiere_update_product_status(status, residual, threshold)
        push!(factors, copy(decomposition.left_factors[1]))
        current = copy(decomposition.right_factors[1])
        leading_scale *= leading
    end
    push!(factors, current)
    factors[1] .*= leading_scale
    reconstruction = tensor_product(factors...)
    reconstruction_residual = norm(dense - reconstruction)
    approximation_threshold = absolute + relative * input_norm
    status = _tiere_update_product_status(
        status, reconstruction_residual, approximation_threshold
    )
    return ProductAnalysisResult(
        status,
        _tiere_factor_tuple(factors, row_layout),
        reconstruction_residual,
        approximation_threshold,
        cut_residuals,
        cut_thresholds,
        _tiere_product_message(status, "operator"),
    )
end

function _tiere_entropy_from_probabilities(probabilities, base)
    real_type = eltype(probabilities)
    converted_base = try
        convert(real_type, base)
    catch error
        error isa InexactError || rethrow()
        throw(
            ArgumentError(
                "base=$(repr(base)) cannot be represented in the state's " *
                "$real_type precision",
            ),
        )
    end
    isfinite(converted_base) && converted_base > one(converted_base) || throw(
        ArgumentError(
            "base=$(repr(base)) must remain finite and greater than one in " *
            "the state's $real_type precision",
        ),
    )
    result = zero(real_type)
    for probability in probabilities
        probability >= zero(probability) || throw(
            DomainError(
                probability,
                "entropy probabilities must be nonnegative; no clipping is performed",
            ),
        )
        iszero(probability) && continue
        result -= probability * log(probability)
    end
    return result / log(converted_base)
end

"""
    entanglement_of_formation(
        state, dims; base=2, atol=nothing, rtol=nothing,
        allow_densify=false, psd_boundary_policy=:project,
        range_boundary_policy=:project
    )

Return the closed-form bipartite entanglement of formation in the domains where
it is known exactly:

- a validated normalized pure vector in arbitrary bipartite dimensions, using
  the entropy of its squared Schmidt coefficients;
- a validated `4 x 4` two-qubit density matrix with `dims == (2, 2)`, using
  Wootters concurrence and binary entropy.

`base` defaults to `2` and must be finite and greater than one. Sparse spectral
work requires `allow_densify=true`, and the dependency-free implementation
supports only BLAS floating element types. Mixed states outside `2 x 2` local
dimensions are rejected, including rank-one matrices: callers must pass pure
states as vectors. For the mixed-state formula, `psd_boundary_policy=:project`
explicitly projects only negative eigenvalues already proven to lie inside
the requested PSD validation tolerance to zero; use `:reject` to refuse that
roundoff-boundary repair. Eigenvalues below the tolerance always raise.
`range_boundary_policy` analogously controls a bounded projection of a
mixed-state concurrence just above its exact upper endpoint.
Inputs and intermediate probabilities are never normalized. Complexity is
dominated by dense SVD/eigendecomposition work.
"""
function entanglement_of_formation(
    state::AbstractVector{<:Number},
    dims;
    base=2,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    checked_base = _tierd_validate_log_base(base)
    _tierd_bipartite_layout(dims, length(state))
    _tierd_validate_pure_state(
        state; atol=atol, rtol=rtol, operation="entanglement_of_formation"
    )
    coefficients = schmidt_coefficients(state, dims; allow_densify=allow_densify)
    probabilities = abs2.(coefficients)
    return _tiere_entropy_from_probabilities(probabilities, checked_base)
end

function entanglement_of_formation(
    rho::AbstractMatrix{<:Number},
    dims;
    base=2,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    psd_boundary_policy::Symbol=:project,
    range_boundary_policy::Symbol=:project,
)
    checked_base = _tierd_validate_log_base(base)
    layout = _tierd_bipartite_layout(dims, size(rho, 1))
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("rho must be square; got size $(size(rho))"))
    layout.dims == (2, 2) || throw(
        ArgumentError(
            "mixed-state entanglement_of_formation is implemented only " *
            "for dims=(2, 2); got $(layout.dims)",
        ),
    )
    concurrence_value = concurrence(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        psd_boundary_policy=psd_boundary_policy,
        range_boundary_policy=range_boundary_policy,
    )
    radicand = one(concurrence_value) - concurrence_value^2
    radicand >= zero(radicand) || throw(
        DomainError(
            concurrence_value,
            "computed concurrence exceeds one; refusing to clip the " *
            "Wootters formula radicand",
        ),
    )
    first_probability = (one(concurrence_value) + sqrt(radicand)) / 2
    probabilities = [first_probability, one(first_probability) - first_probability]
    return _tiere_entropy_from_probabilities(probabilities, checked_base)
end

function _tiere_separable_ball_result(
    purity_value, dimension::Int, absolute, relative; input_boundary_uncertain::Bool
)
    dimension > 1 || throw(
        ArgumentError(
            "the Gurvits--Barnum ball formula requires total dimension greater than one"
        ),
    )
    boundary = one(purity_value) / (dimension - 1)
    scale = max(one(purity_value), abs(purity_value), abs(boundary))
    tolerance = _tierd_threshold(scale, absolute, relative)
    status = if input_boundary_uncertain
        :unknown
    elseif purity_value <= boundary - tolerance
        :separable_certified
    elseif purity_value > boundary + tolerance
        :outside_ball
    else
        :unknown
    end
    message = if status === :separable_certified
        "the Gurvits--Barnum purity bound certifies bipartite separability"
    elseif status === :outside_ball
        "the state is outside this sufficient separable ball; no entanglement conclusion follows"
    elseif input_boundary_uncertain
        "input validation has a nonzero residual inside the requested numerical tolerance, so no certificate is reported"
    else
        "the state lies inside the numerical tolerance band around the separable-ball boundary"
    end
    return SeparableBallResult(status, purity_value, boundary, tolerance, message)
end

function _tiere_eigenvalue_analysis(eigenvalues::AbstractVector{<:Number}; atol, rtol)
    !isempty(eigenvalues) || throw(ArgumentError("eigenvalues must have positive length"))
    eltype(eigenvalues) <: AbstractFloat || throw(
        ArgumentError(
            "the eigenvalue-vector path requires a concrete floating " *
            "element type; got $(eltype(eigenvalues)). No implicit " *
            "precision-changing conversion is performed.",
        ),
    )
    isconcretetype(eltype(eigenvalues)) || throw(
        ArgumentError(
            "eigenvalues must have a concrete floating element type; got " *
            "$(eltype(eigenvalues))",
        ),
    )
    all(isfinite, eigenvalues) ||
        throw(ArgumentError("eigenvalues must contain only finite entries"))
    real_type = eltype(eigenvalues)
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    trace_value = sum(eigenvalues)
    trace_tolerance = _tierd_threshold(abs(trace_value), absolute, relative)
    normalization_residual = abs(trace_value - one(trace_value))
    normalization_residual <= trace_tolerance || throw(
        ArgumentError(
            "eigenvalues are not normalized within atol=$absolute and " *
            "rtol=$relative; their sum is $trace_value. The input is " *
            "never normalized implicitly.",
        ),
    )
    entry_scale = maximum(abs, eigenvalues; init=zero(real_type))
    positivity_tolerance = _tierd_threshold(entry_scale, absolute, relative)
    minimum_value = minimum(eigenvalues)
    minimum_value < -positivity_tolerance && throw(
        DomainError(
            minimum_value,
            "eigenvalues are not nonnegative within the requested tolerances",
        ),
    )
    return (
        purity=sum(abs2, eigenvalues),
        atol=absolute,
        rtol=relative,
        boundary_uncertain=minimum_value < zero(minimum_value) ||
                           !iszero(normalization_residual),
        minimum_value=minimum_value,
        normalization_residual=normalization_residual,
    )
end

"""
    in_separable_ball(
        rho_or_eigenvalues, dims;
        atol=nothing, rtol=nothing, allow_densify=false
    ) -> SeparableBallResult

Apply the Gurvits--Barnum sufficient separability condition to a validated
bipartite density matrix or to its supplied eigenvalue vector. For total
dimension `D > 1`, a trace-one positive state is certified when

```math
\\operatorname{tr}(\\rho^2) \\leq \\frac{1}{D-1}.
```

A general matrix input is fully density-validated and requires a dense
Hermitian eigendecomposition; sparse matrices therefore require
`allow_densify=true` and only BLAS floating element types are supported.
`Diagonal` inputs take a structure-aware `O(D)` path without densification.
The eigenvalue-vector path likewise performs no densification and supports
concrete floating types, including `BigFloat`.

Unlike the pinned QETLAB entry point, this API never silently divides by the
trace. Inputs must already be normalized and are never clipped or repaired.
Failure of this sufficient test returns `:outside_ball` or `:unknown`, never an
entanglement result. Matrix validation costs `O(D^3)`; the trusted-eigenvalue
path costs `O(D)`.
"""
function in_separable_ball(
    rho::Diagonal{<:Number}, dims; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    layout = _tierd_bipartite_layout(dims, size(rho, 1))
    diagonal = diag(rho)
    _tierd_require_finite(diagonal, "diagonal entries")
    real_type = typeof(real(zero(eltype(diagonal))))
    real_type <: AbstractFloat || throw(
        ArgumentError(
            "the structure-aware Diagonal path requires real or complex " *
            "floating entries; got eltype $(eltype(diagonal))",
        ),
    )
    absolute, relative = _tierd_tolerances(real_type, atol, rtol)
    scale = maximum(abs, diagonal; init=zero(real_type))
    validation_tolerance = _tierd_threshold(scale, absolute, relative)
    imaginary_residual = maximum(abs ∘ imag, diagonal; init=zero(real_type))
    imaginary_residual <= validation_tolerance || throw(
        ArgumentError(
            "a diagonal density matrix must have real diagonal entries; " *
            "maximum imaginary residual is $imaginary_residual",
        ),
    )
    real_diagonal = real.(diagonal)
    analysis = _tiere_eigenvalue_analysis(real_diagonal; atol=absolute, rtol=relative)
    return _tiere_separable_ball_result(
        analysis.purity,
        layout.total_dimension,
        analysis.atol,
        analysis.rtol;
        input_boundary_uncertain=analysis.boundary_uncertain || !iszero(imaginary_residual),
    )
end

function in_separable_ball(
    rho::AbstractMatrix{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("rho must be square; got size $(size(rho))"))
    layout = _tierd_bipartite_layout(dims, size(rho, 1))
    analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="in_separable_ball",
    )
    purity_value = sum(abs2, analysis.eigenvalues)
    return _tiere_separable_ball_result(
        purity_value,
        layout.total_dimension,
        analysis.atol,
        analysis.rtol;
        input_boundary_uncertain=analysis.boundary_uncertain,
    )
end

function in_separable_ball(
    eigenvalues::AbstractVector{<:Number},
    dims;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
)
    allow_densify && throw(
        ArgumentError("allow_densify is not meaningful for the eigenvalue-vector path")
    )
    layout = _tierd_bipartite_layout(dims, length(eigenvalues))
    analysis = _tiere_eigenvalue_analysis(eigenvalues; atol=atol, rtol=rtol)
    return _tiere_separable_ball_result(
        analysis.purity,
        layout.total_dimension,
        analysis.atol,
        analysis.rtol;
        input_boundary_uncertain=analysis.boundary_uncertain,
    )
end
