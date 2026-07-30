# Source-informed independent implementation of the public behavior in QETLAB
# SkOperatorNorm.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The native contract intentionally separates witness-backed lower bounds,
# analytic upper bounds, numerically validated SDP relaxations, and exact
# branches. Randomized searches require an explicit RNG and bounded work.

using Random

"""
    SKOperatorNormStatus

Outcome classification for [`sk_operator_norm`](@ref). Only
`SKOperatorNormExact` claims an analytic exact value. A closed numerical gap,
an optimizer result, or a randomized search remains a bounds result.
"""
@enum SKOperatorNormStatus::UInt8 begin
    SKOperatorNormExact
    SKOperatorNormBounds
    SKOperatorNormTargetAbove
    SKOperatorNormTargetBelow
    SKOperatorNormResourceLimited
    SKOperatorNormNumericalBoundary
end

"""
    SKOperatorNormLowerWitness

Explicit Schmidt-rank-constrained vectors verifying a lower bound. The value is
`abs(left' * operator * right)`; `quadratic=true` means the two vectors are
identical and the value is a positive-semidefinite quadratic expectation.
"""
struct SKOperatorNormLowerWitness{T<:Real,V<:AbstractVector{Complex{T}}}
    left::V
    right::V
    value::T
    left_schmidt_rank::Int
    right_schmidt_rank::Int
    quadratic::Bool
    kind::Symbol
    validation_residual::T
    validated::Bool
end

"""
    SKOperatorNormUpperWitness

Auditable upper-relaxation evidence. `numerically_validated=true` means the
package checked the solver's cone and stationarity residuals at `tolerance`;
it does not turn floating-point evidence into an exact symbolic certificate.
"""
struct SKOperatorNormUpperWitness{T<:Real,D}
    bound::T
    formulation::Symbol
    hierarchy_level::Int
    dual_data::D
    residual::Union{Nothing,T}
    tolerance::T
    numerically_validated::Bool
    exact::Bool
end

"""
    SKOperatorNormProblem

Package-owned SDP relaxation for the positive-semidefinite S(`k`) operator
norm. Level one is the PPT relaxation for `k == 1` and the reduction-type
relaxation for `k > 1`. Higher levels use a bosonic-compressed PPT symmetric
extension and are available only for `k == 1`.
"""
struct SKOperatorNormProblem{T<:Real,O,P,S,E,B,M}
    operator::O
    dimensions::NTuple{2,Int}
    k::Int
    hierarchy_level::Int
    formulation::Symbol
    program::P
    state_view::S
    extension_view::E
    bosonic_basis::B
    ambient_dimension::Int
    variable_matrix_dimension::Int
    tolerance::T
    metadata::M
end

"""
    SKOperatorNormResult

Status-rich bounds for the S(`k`) operator norm.

`lower_bound` and `upper_bound` delimit the value. A lower witness is present
only when explicit Schmidt-rank-constrained vectors attain
`witness_lower_bound`; theorem-only lower bounds are labeled separately.
`relaxation_upper_bound` is a solver-derived numerical bound and is never
silently presented as exact. `optimization_result` retains all backend
termination and residual metadata.
"""
struct SKOperatorNormResult{T<:Real,L,U,P,O}
    status::SKOperatorNormStatus
    lower_bound::T
    upper_bound::T
    lower_bound_kind::Symbol
    upper_bound_kind::Symbol
    witness_lower_bound::T
    lower_witness::L
    upper_witness::U
    relaxation_upper_bound::Union{Nothing,T}
    exact::Bool
    dimensions::NTuple{2,Int}
    k::Int
    hierarchy_level::Int
    target::Union{Nothing,T}
    target_decision::Union{Nothing,Symbol}
    random_restarts::Int
    random_iterations::Int
    work_used::BigInt
    problem::P
    optimization_result::O
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

function Base.show(io::IO, result::SKOperatorNormResult)
    return print(
        io,
        "SKOperatorNormResult(status=",
        result.status,
        ", bounds=(",
        result.lower_bound,
        ", ",
        result.upper_bound,
        "), k=",
        result.k,
        ", exact=",
        result.exact,
        ")",
    )
end

_sknorm_real_type(::Type{T}) where {T<:Real} = T
_sknorm_real_type(::Type{Complex{T}}) where {T<:Real} = T

function _sknorm_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _sknorm_nonnegative_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _sknorm_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be a positive integer or nothing"))
    return BigInt(value)
end

function _sknorm_dimensions(dims, dimension::Int)
    if dims === nothing
        root = isqrt(dimension)
        root * root == dimension || throw(
            ArgumentError(
                "dims must be supplied when operator dimension $dimension is not " *
                "a perfect square",
            ),
        )
        return (root, root)
    elseif dims isa Integer
        first_dimension = _sknorm_positive_integer(dims, "dims")
        dimension % first_dimension == 0 || throw(
            DimensionMismatch(
                "scalar dims=$first_dimension does not divide operator dimension " *
                "$dimension",
            ),
        )
        return (first_dimension, dimension ÷ first_dimension)
    elseif dims isa Union{Tuple,AbstractVector}
        dims isa AbstractVector && Base.require_one_based_indexing(dims)
        length(dims) == 2 ||
            throw(ArgumentError("dims must contain exactly two local dimensions"))
        checked = (
            _sknorm_positive_integer(dims[1], "dims[1]"),
            _sknorm_positive_integer(dims[2], "dims[2]"),
        )
        BigInt(checked[1]) * checked[2] == dimension || throw(
            DimensionMismatch(
                "prod(dims)=$(BigInt(checked[1]) * checked[2]) does not match " *
                "operator dimension $dimension",
            ),
        )
        return checked
    end
    return throw(
        ArgumentError(
            "dims must be nothing, a positive integer, or a two-entry tuple/vector"
        ),
    )
end

function _sknorm_tolerance(::Type{R}, scale; atol, rtol) where {R<:AbstractFloat}
    absolute = atol === nothing ? zero(R) : atol
    relative = rtol === nothing ? sqrt(eps(R)) : rtol
    for (name, value) in (("atol", absolute), ("rtol", relative))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    a, r, s = promote(absolute, relative, scale)
    return convert(R, a + r * max(one(s), s))
end

function _sknorm_prepare_operator(
    operator::AbstractMatrix{<:Number},
    dims;
    atol,
    rtol,
    allow_densify::Bool,
    max_dense_entries,
    max_work,
    require_hermitian::Bool=false,
    require_psd::Bool=false,
)
    Base.require_one_based_indexing(operator)
    size(operator, 1) == size(operator, 2) ||
        throw(DimensionMismatch("S(k) operator norms require a square matrix"))
    dimension = size(operator, 1)
    dimension > 0 || throw(ArgumentError("operator must not be empty"))
    dimensions = _sknorm_dimensions(dims, dimension)
    for (index, value) in pairs(operator)
        value isa Number && !(value isa Bool) ||
            throw(ArgumentError("operator entry $index must be numeric and non-Boolean"))
        isfinite(value) || throw(ArgumentError("operator entry $index must be finite"))
    end
    require_hermitian &&
        !ishermitian(operator) &&
        throw(
            ArgumentError(
                "the requested formulation requires an exactly Hermitian operator; " *
                "no symmetrization is applied",
            ),
        )
    dense_limit = _sknorm_limit(max_dense_entries, "max_dense_entries")
    needed_entries = BigInt(dimension)^2
    dense_limit !== nothing &&
        needed_entries > dense_limit &&
        throw(
            ArgumentError(
                "the spectral path needs $needed_entries dense entries, exceeding " *
                "max_dense_entries=$dense_limit",
            ),
        )
    issparse(operator) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "the spectral path requires densification of sparse input; pass " *
                "allow_densify=true after reviewing the $(size(operator)) allocation",
            ),
        )
    work_limit = _sknorm_limit(max_work, "max_work")
    spectral_work = BigInt(dimension)^3
    work_limit !== nothing &&
        spectral_work > work_limit &&
        throw(
            ArgumentError(
                "spectral validation needs $spectral_work work units, exceeding " *
                "max_work=$work_limit",
            ),
        )
    dense = Matrix(operator)
    R = _sknorm_real_type(eltype(dense))
    R <: Union{Float32,Float64} || throw(
        ArgumentError(
            "S(k) spectral and optimizer paths require Float32 or Float64 real " *
            "components; got $R",
        ),
    )
    singular_values = svdvals(dense)
    operator_norm = convert(R, maximum(singular_values; init=zero(R)))
    tolerance = _sknorm_tolerance(R, operator_norm; atol, rtol)
    hermitian = ishermitian(dense)
    eigenvalues = hermitian ? real.(eigvals(Hermitian(dense))) : R[]
    minimum_eigenvalue = hermitian ? minimum(eigenvalues) : nothing
    psd = hermitian && minimum_eigenvalue >= zero(R)
    psd_boundary =
        hermitian && minimum_eigenvalue < zero(R) && minimum_eigenvalue >= -tolerance
    require_psd &&
        !psd &&
        throw(
            ArgumentError(
                if psd_boundary
                    "the operator lies on the numerical PSD boundary; no projection is applied"
                else
                    "the requested formulation requires a positive semidefinite operator"
                end,
            ),
        )
    return (
        matrix=dense,
        dimensions=dimensions,
        real_type=R,
        singular_values=singular_values,
        operator_norm=operator_norm,
        tolerance=tolerance,
        hermitian=hermitian,
        eigenvalues=eigenvalues,
        minimum_eigenvalue=minimum_eigenvalue,
        psd=psd,
        psd_boundary=psd_boundary,
        spectral_work=spectral_work,
        work_limit=work_limit,
    )
end

function _sknorm_affine_combination(
    left::HermitianAffineMatrix{T},
    right::HermitianAffineMatrix{T},
    left_weight,
    right_weight,
    name::Symbol,
) where {T<:Real}
    left.dimension == right.dimension ||
        throw(DimensionMismatch("affine matrix dimensions do not match"))
    left.variable_count == right.variable_count ||
        throw(DimensionMismatch("affine variable counts do not match"))
    left_terms = Dict(term.variable => term.coefficient for term in left.terms)
    right_terms = Dict(term.variable => term.coefficient for term in right.terms)
    zero_matrix = spzeros(Complex{T}, left.dimension, left.dimension)
    variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    for variable in sort!(collect(union(keys(left_terms), keys(right_terms))))
        coefficient =
            left_weight * get(left_terms, variable, zero_matrix) +
            right_weight * get(right_terms, variable, zero_matrix)
        nnz(sparse(coefficient)) == 0 && continue
        push!(variables, variable)
        push!(coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name,
        left_weight * left.constant + right_weight * right.constant,
        variables,
        coefficients,
        left.variable_count,
    )
end

function _sknorm_affine_congruence(
    matrix::HermitianAffineMatrix, factor::AbstractMatrix{<:Number}, name::Symbol
)
    Base.require_one_based_indexing(factor)
    size(factor, 2) == matrix.dimension ||
        throw(DimensionMismatch("congruence factor has the wrong column dimension"))
    variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    for term in matrix.terms
        coefficient = factor * term.coefficient * adjoint(factor)
        nnz(sparse(coefficient)) == 0 && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name,
        factor * matrix.constant * adjoint(factor),
        variables,
        coefficients,
        matrix.variable_count,
    )
end

function _sknorm_trace_pairing_affine(
    operator::AbstractMatrix{<:Number}, matrix::HermitianAffineMatrix{T}
) where {T<:Real}
    variables = Int[]
    coefficients = T[]
    for term in matrix.terms
        value = convert(T, real(tr(operator * term.coefficient)))
        iszero(value) && continue
        push!(variables, term.variable)
        push!(coefficients, value)
    end
    constant = convert(T, real(tr(operator * matrix.constant)))
    return AffineScalar(constant, variables, coefficients, matrix.variable_count)
end

function _sknorm_problem_preflight(
    dimensions::NTuple{2,Int}, k::Int, level::Int, limits::OptimizationLimits
)
    dimension_a, dimension_b = dimensions
    if level > 1 && k != 1
        throw(
            ArgumentError(
                "higher symmetric-extension levels are defined only for k == 1; " *
                "use hierarchy_level=1 for k > 1",
            ),
        )
    end
    ambient = if level == 1
        BigInt(dimension_a * dimension_b)
    else
        BigInt(dimension_a) * BigInt(dimension_b)^level
    end
    variable_dimension = if level == 1
        ambient
    else
        BigInt(dimension_a) * binomial(BigInt(dimension_b + level - 1), BigInt(level))
    end
    variable_count = variable_dimension^2
    psd_blocks = BigInt(2)
    largest_psd = max(ambient, variable_dimension)
    real_entries =
        (2variable_dimension) * (2variable_dimension + 1) ÷ 2 +
        (2ambient) * (2ambient + 1) ÷ 2
    reasons = String[]
    variable_count > limits.max_variables && push!(
        reasons,
        "model needs $variable_count variables, exceeding " *
        "max_variables=$(limits.max_variables)",
    )
    psd_blocks > limits.max_psd_blocks && push!(
        reasons,
        "model needs $psd_blocks PSD blocks, exceeding " *
        "max_psd_blocks=$(limits.max_psd_blocks)",
    )
    largest_psd > limits.max_psd_dimension && push!(
        reasons,
        "largest PSD block has dimension $largest_psd, exceeding " *
        "max_psd_dimension=$(limits.max_psd_dimension)",
    )
    real_entries > limits.max_model_entries && push!(
        reasons,
        "real-block PSD scalarization needs $real_entries entries, exceeding " *
        "max_model_entries=$(limits.max_model_entries)",
    )
    return (
        allowed=isempty(reasons),
        reason=join(reasons, "; "),
        ambient=ambient,
        variable_dimension=variable_dimension,
        variable_count=variable_count,
        psd_blocks=psd_blocks,
        real_entries=real_entries,
    )
end

"""
    sk_operator_norm_problem(operator; k=1, dims=nothing,
                             hierarchy_level=1, kwargs...)

Build a package-owned SDP upper relaxation for a positive-semidefinite
bipartite operator. The input must be exactly Hermitian; entries are never
symmetrized, clipped, normalized, or projected.

Level one imposes PPT when `k == 1`, or
`k * Tr_B(rho) ⊗ I_B - rho >= 0` when `k > 1`. Higher levels use a
bosonic-compressed extension with the pinned PPT cut and are restricted to
`k == 1`. All combinatorial and dense real-block sizes are checked before the
model is allocated.
"""
function sk_operator_norm_problem(
    operator::AbstractMatrix{<:Number};
    k=1,
    dims=nothing,
    hierarchy_level=1,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_work=1_000_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_k = _sknorm_positive_integer(k, "k")
    level = _sknorm_positive_integer(hierarchy_level, "hierarchy_level")
    prepared = _sknorm_prepare_operator(
        operator,
        dims;
        atol,
        rtol,
        allow_densify,
        max_dense_entries,
        max_work,
        require_hermitian=true,
        require_psd=true,
    )
    checked_k = min(checked_k, min(prepared.dimensions...))
    preflight = _sknorm_problem_preflight(prepared.dimensions, checked_k, level, limits)
    preflight.allowed || throw(ArgumentError(preflight.reason))
    variable_dimension = Int(preflight.variable_dimension)
    variable_count = Int(preflight.variable_count)
    R = prepared.real_type
    variable = hermitian_variable(
        level == 1 ? :rho : :bosonic_extension,
        variable_dimension;
        variable_count,
        coefficient_type=R,
    )

    basis = nothing
    extension = variable
    state_view = variable
    full_dimensions = prepared.dimensions
    if level > 1
        basis_b = symmetric_subspace_basis(
            prepared.dimensions[2], level; T=R, sparse_output=true
        )
        identity_a = sparse(
            1:(prepared.dimensions[1]),
            1:(prepared.dimensions[1]),
            fill(one(R), prepared.dimensions[1]),
            prepared.dimensions[1],
            prepared.dimensions[1],
        )
        basis = kron(identity_a, basis_b)
        extension = _sknorm_affine_congruence(variable, basis, :extension)
        full_dimensions = (
            prepared.dimensions[1], ntuple(_ -> prepared.dimensions[2], level)...
        )
        state_view = partial_trace_affine(
            extension, full_dimensions; trace_out=Tuple(3:(level + 1)), name=:rho_ab
        )
    end

    psd_constraints = HermitianAffineMatrix[variable]
    formulation = if level == 1 && checked_k == 1
        push!(
            psd_constraints,
            partial_transpose_affine(
                state_view,
                prepared.dimensions;
                systems=(2,),
                name=:rho_partial_transpose,
            ),
        )
        :ppt_relaxation
    elseif level == 1
        reduced = partial_trace_affine(
            state_view, prepared.dimensions; trace_out=(2,), name=:rho_a
        )
        identity_b = sparse(
            1:(prepared.dimensions[2]),
            1:(prepared.dimensions[2]),
            fill(one(R), prepared.dimensions[2]),
            prepared.dimensions[2],
            prepared.dimensions[2],
        )
        lifted = tensor_affine(reduced; right=identity_b, name=:rho_a_tensor_identity)
        push!(
            psd_constraints,
            _sknorm_affine_combination(
                lifted, state_view, checked_k, -one(R), :schmidt_number_reduction
            ),
        )
        :schmidt_number_reduction_relaxation
    else
        last_system = cld(level, 2) + 1
        push!(
            psd_constraints,
            partial_transpose_affine(
                extension,
                full_dimensions;
                systems=Tuple(1:last_system),
                name=:extension_partial_transpose,
            ),
        )
        :bosonic_ppt_symmetric_extension
    end
    trace_interval = AffineInterval(trace_affine(variable), nothing, one(R), :trace_le_one)
    objective = _sknorm_trace_pairing_affine(prepared.matrix, state_view)
    views = if level == 1
        HermitianAffineMatrix[state_view]
    else
        HermitianAffineMatrix[state_view, extension, variable]
    end
    metadata = (
        formulation,
        dimensions=prepared.dimensions,
        k=checked_k,
        hierarchy_level=level,
        ambient_dimension=Int(preflight.ambient),
        variable_matrix_dimension=variable_dimension,
        relaxation_only=true,
        exact_low_dimension_theorem=(
            checked_k == 1 &&
            min(prepared.dimensions...) == 2 &&
            max(prepared.dimensions...) <= 3 &&
            level == 1
        ),
        subsystem_order=:A_then_B_copies,
        preflight,
    )
    program = SemidefiniteProgram(
        :sk_operator_norm_upper_relaxation,
        :maximize,
        variable_count,
        objective;
        intervals=[trace_interval],
        psd_constraints,
        primal_views=views,
        known_feasible_point=zeros(R, variable_count),
        limits,
        metadata,
    )
    return SKOperatorNormProblem(
        copy(prepared.matrix),
        prepared.dimensions,
        checked_k,
        level,
        formulation,
        program,
        state_view,
        extension,
        basis,
        Int(preflight.ambient),
        variable_dimension,
        prepared.tolerance,
        metadata,
    )
end

function _sknorm_schmidt_projection(
    vector::AbstractVector{<:Number}, dimensions::NTuple{2,Int}, k::Int
)
    coefficient_matrix = reshape(vector, dimensions[2], dimensions[1])
    factorization = svd(coefficient_matrix; full=false)
    keep = min(k, length(factorization.S))
    projected =
        factorization.U[:, 1:keep] *
        Diagonal(factorization.S[1:keep]) *
        factorization.Vt[1:keep, :]
    flattened = Complex.(vec(projected))
    magnitude = norm(flattened)
    if iszero(magnitude)
        flattened .= zero(eltype(flattened))
        flattened[1] = one(eltype(flattened))
    else
        flattened ./= magnitude
    end
    rank_value = count(value -> !iszero(value), svdvals(projected))
    return flattened, min(rank_value, k)
end

function _sknorm_witness(operator, left, right, dimensions, k, kind; quadratic::Bool=false)
    T = _sknorm_real_type(eltype(operator))
    left_vector = Complex{T}.(left)
    right_vector = Complex{T}.(right)
    left_vector ./= norm(left_vector)
    right_vector ./= norm(right_vector)
    value = convert(T, abs(dot(left_vector, operator * right_vector)))
    left_rank = schmidt_rank(left_vector, dimensions; atol=zero(T), rtol=sqrt(eps(T)))
    right_rank = schmidt_rank(right_vector, dimensions; atol=zero(T), rtol=sqrt(eps(T)))
    residual = convert(T, abs(value - abs(dot(left_vector, operator * right_vector))))
    validated =
        left_rank <= k && right_rank <= k && residual <= sqrt(eps(T)) * max(one(T), value)
    return SKOperatorNormLowerWitness(
        left_vector,
        right_vector,
        value,
        left_rank,
        right_rank,
        quadratic,
        kind,
        residual,
        validated,
    )
end

function _sknorm_computational_witness(operator, dimensions, k, psd::Bool)
    dimension = size(operator, 1)
    if psd
        index = argmax(real.(diag(operator)))
        vector = zeros(Complex{_sknorm_real_type(eltype(operator))}, dimension)
        vector[index] = one(eltype(vector))
        return _sknorm_witness(
            operator, vector, vector, dimensions, k, :computational_basis; quadratic=true
        )
    end
    row, column = Tuple(argmax(abs.(operator)))
    left = zeros(Complex{_sknorm_real_type(eltype(operator))}, dimension)
    right = zeros(eltype(left), dimension)
    left[row] = one(eltype(left))
    right[column] = one(eltype(right))
    return _sknorm_witness(operator, left, right, dimensions, k, :computational_basis)
end

function _sknorm_projected_search(
    rng::AbstractRNG,
    operator,
    dimensions,
    k,
    restarts,
    max_iterations,
    tolerance,
    psd::Bool,
    initial_witness,
)
    T = _sknorm_real_type(eltype(operator))
    dimension = size(operator, 1)
    best = initial_witness
    iterations_used = 0
    for _ in 1:restarts
        start = randn(rng, Complex{T}, dimension)
        right, _ = _sknorm_schmidt_projection(start, dimensions, k)
        left = right
        previous = -one(T)
        for _ in 1:max_iterations
            iterations_used += 1
            if psd
                candidate, _ = _sknorm_schmidt_projection(operator * right, dimensions, k)
                left = candidate
                right = candidate
                value = convert(T, real(dot(right, operator * right)))
            else
                left, _ = _sknorm_schmidt_projection(operator * right, dimensions, k)
                right, _ = _sknorm_schmidt_projection(
                    adjoint(operator) * left, dimensions, k
                )
                value = convert(T, abs(dot(left, operator * right)))
            end
            if value > best.value
                best = _sknorm_witness(
                    operator,
                    left,
                    right,
                    dimensions,
                    k,
                    :bounded_projected_iteration;
                    quadratic=psd,
                )
            end
            abs(value - previous) <= tolerance * max(one(T), abs(value)) && break
            previous = value
        end
    end
    return best, iterations_used
end

function _sknorm_target_decision(lower, upper, target, tolerance)
    target === nothing && return nothing
    lower >= target + tolerance && return :above
    upper <= target - tolerance && return :below
    return nothing
end

function _sknorm_result(
    status,
    lower,
    upper,
    lower_kind,
    upper_kind,
    witness,
    upper_witness,
    relaxation_upper,
    exact,
    dimensions,
    k,
    level,
    target,
    decision,
    restarts,
    iterations,
    work,
    problem,
    optimization_result,
    tolerance,
    warnings,
    message,
)
    T = promote_type(
        typeof(lower),
        typeof(upper),
        typeof(tolerance),
        target === nothing ? typeof(lower) : typeof(target),
    )
    return SKOperatorNormResult{
        T,typeof(witness),typeof(upper_witness),typeof(problem),typeof(optimization_result)
    }(
        status,
        convert(T, lower),
        convert(T, upper),
        lower_kind,
        upper_kind,
        convert(T, witness === nothing ? zero(T) : witness.value),
        witness,
        upper_witness,
        relaxation_upper === nothing ? nothing : convert(T, relaxation_upper),
        exact,
        dimensions,
        k,
        level,
        target === nothing ? nothing : convert(T, target),
        decision,
        restarts,
        iterations,
        work,
        problem,
        optimization_result,
        convert(T, tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

"""
    sk_operator_norm(rng::AbstractRNG, operator; kwargs...)

Compute auditable lower and upper bounds on the bipartite S(`k`) operator
norm. Randomized projected iterations require the explicit `rng`; the
function never touches Julia's global random stream.

Exact analytic branches cover the zero operator, `k >= min(dims)`, and exact
rank-one matrices. General lower bounds come from theorem estimates and
explicit Schmidt-rank-constrained witnesses. Positive-semidefinite inputs may
also use a bounded projected search and an optional SDP upper relaxation.
Optimizer outcomes, numerical relaxation bounds, and exact values remain
distinct in the returned [`SKOperatorNormResult`](@ref).

`strength=0` disables random and SDP work. For `k == 1`, higher strengths use
the matching bosonic PPT hierarchy level (capped by `max_hierarchy_level`);
for `k > 1`, the reviewed level-one reduction relaxation is used. Sparse input
requires explicit `allow_densify=true`. All dense, work, restart, iteration,
hierarchy, and solver-model budgets are checked before the relevant
allocation or random draw.
"""
function sk_operator_norm(
    rng::AbstractRNG,
    operator::AbstractMatrix{<:Number};
    k=1,
    dims=nothing,
    strength=2,
    target=nothing,
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
    checked_strength = _sknorm_nonnegative_integer(strength, "strength")
    restart_cap = _sknorm_nonnegative_integer(max_restarts, "max_restarts")
    iteration_cap = _sknorm_nonnegative_integer(max_iterations, "max_iterations")
    hierarchy_cap = _sknorm_positive_integer(max_hierarchy_level, "max_hierarchy_level")
    prepared = _sknorm_prepare_operator(
        operator, dims; atol, rtol, allow_densify, max_dense_entries, max_work
    )
    dimensions = prepared.dimensions
    checked_k = min(checked_k, min(dimensions...))
    T = prepared.real_type
    tolerance = prepared.tolerance
    target_value = if target === nothing || target == -1
        nothing
    else
        target isa Real && !(target isa Bool) && isfinite(target) && target >= zero(target) || throw(
            ArgumentError(
                "target must be nothing, -1, or a finite nonnegative real number"
            ),
        )
        convert(T, target)
    end
    matrix = prepared.matrix
    operator_norm = prepared.operator_norm
    warnings = String[]
    work_used = prepared.spectral_work

    if iszero(operator_norm)
        witness = _sknorm_computational_witness(matrix, dimensions, checked_k, true)
        return _sknorm_result(
            SKOperatorNormExact,
            zero(T),
            zero(T),
            :zero_operator,
            :zero_operator,
            witness,
            SKOperatorNormUpperWitness(
                zero(T), :zero_operator, 0, nothing, zero(T), tolerance, true, true
            ),
            nothing,
            true,
            dimensions,
            checked_k,
            0,
            target_value,
            _sknorm_target_decision(zero(T), zero(T), target_value, tolerance),
            0,
            0,
            work_used,
            nothing,
            nothing,
            tolerance,
            warnings,
            "the zero operator has exact S(k) norm zero",
        )
    end

    factorization = svd(matrix; full=false)
    if checked_k >= min(dimensions...)
        witness = _sknorm_witness(
            matrix,
            factorization.U[:, 1],
            adjoint(factorization.Vt)[:, 1],
            dimensions,
            checked_k,
            :top_singular_vectors,
        )
        exact_value = convert(T, operator_norm)
        decision = _sknorm_target_decision(
            exact_value, exact_value, target_value, tolerance
        )
        return _sknorm_result(
            SKOperatorNormExact,
            exact_value,
            exact_value,
            :full_schmidt_rank_exact,
            :full_schmidt_rank_exact,
            witness,
            SKOperatorNormUpperWitness(
                exact_value, :operator_norm, 0, nothing, zero(T), tolerance, true, true
            ),
            nothing,
            true,
            dimensions,
            checked_k,
            0,
            target_value,
            decision,
            0,
            0,
            work_used,
            nothing,
            nothing,
            tolerance,
            warnings,
            "k reaches the smaller local dimension, so the S(k) norm equals " *
            "the operator norm",
        )
    end

    exact_rank = count(!iszero, factorization.S)
    if exact_rank == 1
        left, _ = _sknorm_schmidt_projection(factorization.U[:, 1], dimensions, checked_k)
        right, _ = _sknorm_schmidt_projection(
            adjoint(factorization.Vt)[:, 1], dimensions, checked_k
        )
        witness = _sknorm_witness(
            matrix, left, right, dimensions, checked_k, :rank_one_exact
        )
        exact_value = witness.value
        decision = _sknorm_target_decision(
            exact_value, exact_value, target_value, tolerance
        )
        return _sknorm_result(
            SKOperatorNormExact,
            exact_value,
            exact_value,
            :rank_one_exact,
            :rank_one_exact,
            witness,
            SKOperatorNormUpperWitness(
                exact_value,
                :rank_one_schmidt_formula,
                0,
                nothing,
                zero(T),
                tolerance,
                true,
                true,
            ),
            nothing,
            true,
            dimensions,
            checked_k,
            0,
            target_value,
            decision,
            0,
            0,
            work_used,
            nothing,
            nothing,
            tolerance,
            warnings,
            "the exact rank-one Schmidt formula determines the norm",
        )
    elseif count(value -> value > tolerance, factorization.S) == 1
        push!(
            warnings,
            "the operator is numerically rank one but not exactly rank one; the " *
            "rank-one shortcut was not used",
        )
    end

    lower = convert(T, checked_k / min(dimensions...) * operator_norm)
    lower_kind = :analytic_norm_equivalence
    upper = convert(T, operator_norm)
    upper_kind = :operator_norm
    witness = _sknorm_computational_witness(matrix, dimensions, checked_k, prepared.psd)
    if witness.value > lower
        lower = witness.value
        lower_kind = :explicit_witness
    end

    if prepared.psd
        eigendecomposition = eigen(Hermitian(matrix))
        top_vector = eigendecomposition.vectors[:, argmax(eigendecomposition.values)]
        projected, _ = _sknorm_schmidt_projection(top_vector, dimensions, checked_k)
        candidate = _sknorm_witness(
            matrix,
            projected,
            projected,
            dimensions,
            checked_k,
            :projected_top_eigenvector;
            quadratic=true,
        )
        if candidate.value > witness.value
            witness = candidate
        end
        if candidate.value > lower
            lower = candidate.value
            lower_kind = :explicit_witness
        end
    elseif prepared.psd_boundary
        push!(
            warnings,
            "the operator lies inside the numerical PSD boundary band; PSD-only " *
            "random and SDP routes were not used",
        )
    end

    decision = _sknorm_target_decision(lower, upper, target_value, tolerance)
    if decision !== nothing
        status = decision === :above ? SKOperatorNormTargetAbove : SKOperatorNormTargetBelow
        return _sknorm_result(
            status,
            lower,
            upper,
            lower_kind,
            upper_kind,
            witness,
            SKOperatorNormUpperWitness(
                upper, :operator_norm, 0, nothing, zero(T), tolerance, true, false
            ),
            nothing,
            false,
            dimensions,
            checked_k,
            0,
            target_value,
            decision,
            0,
            0,
            work_used,
            nothing,
            nothing,
            tolerance,
            warnings,
            "the analytic bounds meet the requested target before random or " *
            "optimizer work",
        )
    end

    requested_restarts = if checked_strength == 0 || restart_cap == 0 || iteration_cap == 0
        0
    else
        requested = min(BigInt(5)^checked_strength, BigInt(restart_cap))
        Int(requested)
    end
    random_iterations = 0
    actual_restarts = requested_restarts
    if requested_restarts > 0
        random_work =
            BigInt(requested_restarts) *
            BigInt(iteration_cap) *
            (BigInt(size(matrix, 1))^2 + BigInt(size(matrix, 1)) * min(dimensions...))
        if prepared.work_limit !== nothing && work_used + random_work > prepared.work_limit
            actual_restarts = 0
            push!(
                warnings,
                "the randomized lower-bound search was skipped because its " *
                "$random_work work units exceed max_work=$(prepared.work_limit)",
            )
        elseif iteration_cap > 0
            witness, random_iterations = _sknorm_projected_search(
                rng,
                matrix,
                dimensions,
                checked_k,
                requested_restarts,
                iteration_cap,
                tolerance,
                prepared.psd,
                witness,
            )
            work_per_iteration =
                BigInt(size(matrix, 1))^2 + BigInt(size(matrix, 1)) * min(dimensions...)
            work_used += BigInt(random_iterations) * work_per_iteration
            if witness.value > lower
                lower = witness.value
                lower_kind = :explicit_randomized_witness
            end
        end
    end

    problem = nothing
    optimization_result = nothing
    relaxation_upper = nothing
    upper_witness = SKOperatorNormUpperWitness(
        upper, :operator_norm, 0, nothing, zero(T), tolerance, true, false
    )
    hierarchy_level = 0
    if checked_strength > 0 && prepared.psd
        requested_level = checked_k == 1 ? min(checked_strength, hierarchy_cap) : 1
        for level in requested_level:-1:1
            built = try
                sk_operator_norm_problem(
                    matrix;
                    k=checked_k,
                    dims=dimensions,
                    hierarchy_level=level,
                    atol=zero(T),
                    rtol=tolerance / max(one(T), operator_norm),
                    allow_densify=true,
                    max_dense_entries,
                    max_work,
                    limits,
                )
            catch error
                if error isa ArgumentError && (
                    occursin("exceeding", sprint(showerror, error)) ||
                    occursin("needs", sprint(showerror, error))
                )
                    push!(
                        warnings,
                        "hierarchy level $level was skipped by preflight: " *
                        sprint(showerror, error),
                    )
                    nothing
                else
                    rethrow()
                end
            end
            built === nothing && continue
            problem = built
            hierarchy_level = level
            optimization_result = solve_optimization(built.program, backend)
            if optimization_result.status === OptimizationOptimal &&
                optimization_result.dual !== nothing &&
                optimization_result.dual_residual !== nothing &&
                optimization_result.dual_residual <= tolerance
                raw_bound = if optimization_result.objective_bound !== nothing
                    optimization_result.objective_bound
                else
                    optimization_result.dual_objective_value
                end
                if raw_bound !== nothing
                    safety = tolerance * (one(T) + sqrt(T(built.program.variable_count)))
                    candidate_upper = convert(T, raw_bound + safety)
                    if candidate_upper + tolerance >= lower
                        relaxation_upper = candidate_upper
                        upper = min(upper, candidate_upper)
                        upper_kind = :numerically_validated_sdp_relaxation
                        upper_witness = SKOperatorNormUpperWitness(
                            candidate_upper,
                            built.formulation,
                            level,
                            optimization_result.dual,
                            convert(T, optimization_result.dual_residual),
                            tolerance,
                            true,
                            false,
                        )
                    else
                        push!(
                            warnings,
                            "the numerical SDP upper bound is below a validated lower " *
                            "witness and was discarded",
                        )
                    end
                end
            elseif optimization_result.status in (
                OptimizationLimit,
                OptimizationNumericalFailure,
                OptimizationInconsistent,
                OptimizationMalformedBackend,
            )
                push!(
                    warnings,
                    "the optional upper relaxation was inconclusive: " *
                    optimization_result.message,
                )
            end
            break
        end
    end

    decision = _sknorm_target_decision(lower, upper, target_value, tolerance)
    status = if decision === :above
        SKOperatorNormTargetAbove
    elseif decision === :below
        SKOperatorNormTargetBelow
    elseif any(
        warning -> occursin("skipped", warning) || occursin("preflight", warning),
        warnings,
    )
        SKOperatorNormResourceLimited
    elseif prepared.psd_boundary
        SKOperatorNormNumericalBoundary
    else
        SKOperatorNormBounds
    end
    message = if decision === :above
        "a validated lower bound is above the requested target"
    elseif decision === :below
        "a validated upper bound is below the requested target"
    elseif upper - lower <= tolerance
        "the reported numerical bounds close within tolerance; this is not an " *
        "analytic exactness claim"
    else
        "the result contains certified analytic/witness bounds and any separate " *
        "numerical relaxation evidence"
    end
    return _sknorm_result(
        status,
        lower,
        upper,
        lower_kind,
        upper_kind,
        witness,
        upper_witness,
        relaxation_upper,
        false,
        dimensions,
        checked_k,
        hierarchy_level,
        target_value,
        decision,
        actual_restarts,
        random_iterations,
        work_used,
        problem,
        optimization_result,
        tolerance,
        warnings,
        message,
    )
end
