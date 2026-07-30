# Source-informed independent Julia implementation based on the specifications
# in QETLAB SymmetricExtension.m, SymmetricInnerExtension.m,
# RandomPPTState.m, and helpers/jacobi_poly.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# Mathematical formulations were independently checked against:
# A. C. Doherty, P. A. Parrilo, and F. M. Spedalieri,
# https://arxiv.org/abs/quant-ph/0308032;
# M. Navascués, M. Owari, and M. B. Plenio,
# https://arxiv.org/abs/0906.2735;
# J. Chen, Z. Ji, D. Kribs, N. Lütkenhaus, and B. Zeng,
# https://arxiv.org/abs/1310.3530; and
# J. M. Leinaas, J. Myrheim, and P. Ø. Sollid,
# https://arxiv.org/abs/1002.1949.

using Random

"""
    SymmetricExtensionStatus

Certificate-aware outcome of [`symmetric_extension`](@ref) or
[`symmetric_inner_extension`](@ref). A backend failure, resource limit, or
floating-point boundary is always inconclusive.
"""
@enum SymmetricExtensionStatus::UInt8 begin
    SymmetricExtensionExactPresent
    SymmetricExtensionExactAbsent
    SymmetricExtensionAnalyticPresent
    SymmetricExtensionAnalyticAbsent
    SymmetricExtensionSolverPresent
    SymmetricExtensionSolverAbsent
    SymmetricExtensionNumericalBoundary
    SymmetricExtensionResourceLimit
    SymmetricExtensionBackendUnavailable
    SymmetricExtensionBackendFailure
    SymmetricExtensionInvalidCertificate
end

"""
    SymmetricExtensionProblem

Package-owned, solver-independent model for a `k`-copy extension of subsystem
`B`. `program` is a [`SemidefiniteProgram`](@ref), never a JuMP or optimizer
object. `extension_view` is the ambient `A ⊗ B^k` affine matrix reconstructed
from the program coordinates.
"""
struct SymmetricExtensionProblem{T<:Real,S,P,E,B,M}
    state::S
    dimensions::NTuple{2,Int}
    order::Int
    ppt::Bool
    bosonic::Bool
    formulation::Symbol
    program::P
    extension_view::E
    bosonic_basis::B
    marginal_equality_count::Int
    symmetry_equality_count::Int
    ppt_block_count::Int
    ambient_dimension::Int
    variable_matrix_dimension::Int
    tolerance::T
    metadata::M
end

"""
    SymmetricInnerExtensionProblem

Package-owned Navascués--Owari--Plenio inner-hierarchy model. `mixing_parameter`
is the Jacobi-root coefficient used by the PPT hierarchy; it is zero for the
non-PPT hierarchy. The modeled positive matrix is compressed to the bosonic
subspace and `extension_view` is its ambient lift.
"""
struct SymmetricInnerExtensionProblem{T<:Real,S,P,E,B,M}
    state::S
    dimensions::NTuple{2,Int}
    order::Int
    ppt::Bool
    program::P
    extension_view::E
    bosonic_basis::B
    mixing_parameter::T
    marginal_equality_count::Int
    ppt_block_count::Int
    ambient_dimension::Int
    variable_matrix_dimension::Int
    tolerance::T
    metadata::M
end

"""
    SymmetricExtensionWitness

Numerically checked dual separator reconstructed from the marginal-equality
dual coordinates. For an outer symmetric-extension infeasibility certificate,
`entanglement_witness=true` means the separator was checked against the full
dual stationarity/cone residual. For an inner-hierarchy separator this field is
always false: being outside an inner approximation does **not** prove
entanglement.
"""
struct SymmetricExtensionWitness{M,T,D}
    operator::M
    expectation::T
    normalization_residual::T
    dual_stationarity_residual::Union{Nothing,T}
    tolerance::T
    dual_data::D
    separator_validated::Bool
    entanglement_witness::Bool
    warning::String
end

"""
    SymmetricExtensionResult

Status-rich result for the outer or inner symmetric-extension hierarchy.

`verdict` is `true` or `false` only for a validated analytic, primal, or dual
certificate. `optimization_result` retains the backend termination, primal,
dual, objective, gap, residual, version, option, and resource metadata.
`residuals` records PSD, trace, marginal, permutation, bosonic-support, and
requested partial-transpose checks whenever an extension candidate exists.
"""
struct SymmetricExtensionResult{P,E,W,O,R,T}
    status::SymmetricExtensionStatus
    verdict::Union{Nothing,Bool}
    certificate_kind::Union{Nothing,Symbol}
    hierarchy::Symbol
    dimensions::NTuple{2,Int}
    order::Int
    ppt::Bool
    bosonic::Bool
    problem::P
    extension::E
    witness::W
    optimization_result::O
    residuals::R
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

function Base.show(io::IO, result::SymmetricExtensionResult)
    return print(
        io,
        "SymmetricExtensionResult(status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", hierarchy=",
        result.hierarchy,
        ", order=",
        result.order,
        ", ppt=",
        result.ppt,
        ", bosonic=",
        result.bosonic,
        ")",
    )
end

"""
    RandomPPTStatus

Outcome of [`random_ppt_state`](@ref). Only `RandomPPTConstructed` carries a
verified state.
"""
@enum RandomPPTStatus::UInt8 begin
    RandomPPTConstructed
    RandomPPTResourceLimit
    RandomPPTIterationLimit
    RandomPPTNumericalFailure
    RandomPPTVerificationFailed
end

"""
    RandomPPTStateResult

Auditable random-PPT construction result. `state` is populated only after
independent PSD, trace, requested partial-transpose, and rank-bound checks.
`candidate` retains an unverified failed candidate without presenting it as a
density matrix.
"""
struct RandomPPTStateResult{S,C,H,T}
    status::RandomPPTStatus
    state::S
    candidate::C
    dimensions::NTuple{2,Int}
    requested_ranks::NTuple{2,Int}
    numerical_ranks::Union{Nothing,NTuple{2,Int}}
    construction::Symbol
    real_output::Bool
    iterations::Int
    max_iterations::Int
    convergence_history::H
    minimum_eigenvalue::Union{Nothing,T}
    minimum_partial_transpose_eigenvalue::Union{Nothing,T}
    trace_residual::Union{Nothing,T}
    hermiticity_residual::Union{Nothing,T}
    partial_transpose_hermiticity_residual::Union{Nothing,T}
    tolerance::T
    work_used::BigInt
    max_work::Union{Nothing,BigInt}
    max_dense_entries::Union{Nothing,BigInt}
    normalized_by_construction::Bool
    verified::Bool
    message::String
end

function Base.show(io::IO, result::RandomPPTStateResult)
    return print(
        io,
        "RandomPPTStateResult(status=",
        result.status,
        ", dimensions=",
        result.dimensions,
        ", requested_ranks=",
        result.requested_ranks,
        ", construction=",
        result.construction,
        ", verified=",
        result.verified,
        ")",
    )
end

function _symext_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _symext_nonnegative_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _symext_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be a positive integer or nothing"))
    return BigInt(value)
end

function _symext_checked_int(value::BigInt, name::AbstractString)
    value <= typemax(Int) ||
        throw(ArgumentError("$name=$value exceeds the addressable Int range"))
    return Int(value)
end

function _symext_dimensions(dims, dimension::Int)
    if dims === nothing
        root = isqrt(dimension)
        root * root == dimension || throw(
            ArgumentError(
                "dims must be supplied when matrix dimension $dimension is not a perfect square",
            ),
        )
        return (root, root)
    elseif dims isa Integer
        first_dimension = _symext_positive_integer(dims, "dims")
        dimension % first_dimension == 0 || throw(
            DimensionMismatch(
                "scalar dims=$first_dimension does not divide matrix dimension $dimension",
            ),
        )
        return (first_dimension, dimension ÷ first_dimension)
    elseif dims isa Union{Tuple,AbstractVector}
        dims isa AbstractVector && Base.require_one_based_indexing(dims)
        length(dims) == 2 ||
            throw(ArgumentError("dims must contain exactly two local dimensions"))
        checked = (
            _symext_positive_integer(dims[1], "dims[1]"),
            _symext_positive_integer(dims[2], "dims[2]"),
        )
        BigInt(checked[1]) * checked[2] == dimension || throw(
            DimensionMismatch(
                "prod(dims)=$(BigInt(checked[1]) * checked[2]) does not match matrix dimension $dimension",
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

_symext_real_type(::Type{T}) where {T<:Real} = T
_symext_real_type(::Type{Complex{T}}) where {T<:Real} = T

function _symext_tolerance(::Type{R}, scale; atol, rtol) where {R<:Real}
    exact = R <: Union{Integer,Rational}
    absolute = atol === nothing ? zero(R) : atol
    relative = rtol === nothing ? (exact ? zero(R) : sqrt(eps(R))) : rtol
    for (name, value) in (("atol", absolute), ("rtol", relative))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    exact &&
        (!iszero(absolute) || !iszero(relative)) &&
        throw(ArgumentError("exact inputs require atol=rtol=0"))
    a, r, s = promote(absolute, relative, scale)
    return a + r * max(one(s), s)
end

function _symext_prepare_operator(
    state::AbstractMatrix{<:Number},
    dims;
    atol,
    rtol,
    allow_densify::Bool,
    max_dense_entries,
    operation::AbstractString,
)
    Base.require_one_based_indexing(state)
    size(state, 1) == size(state, 2) ||
        throw(DimensionMismatch("$operation requires a square matrix"))
    dimension = size(state, 1)
    dimension > 0 || throw(ArgumentError("$operation requires a nonempty matrix"))
    dimensions = _symext_dimensions(dims, dimension)
    for (index, value) in pairs(state)
        value isa Number && !(value isa Bool) ||
            throw(ArgumentError("$operation entry $index must be numeric and non-Boolean"))
        isfinite(value) || throw(ArgumentError("$operation entry $index must be finite"))
    end
    ishermitian(state) || throw(
        ArgumentError(
            "$operation requires an exactly Hermitian input; no implicit symmetrization is applied",
        ),
    )
    entry_limit = _symext_limit(max_dense_entries, "max_dense_entries")
    needed_entries = BigInt(dimension)^2
    entry_limit !== nothing &&
        needed_entries > entry_limit &&
        throw(
            ArgumentError(
                "$operation needs $needed_entries dense entries, exceeding max_dense_entries=$entry_limit",
            ),
        )
    if issparse(state) && !allow_densify
        throw(
            ArgumentError(
                "$operation requires dense spectral validation; pass allow_densify=true after reviewing the $(size(state)) allocation",
            ),
        )
    end
    dense = Matrix(state)
    R = _symext_real_type(eltype(dense))
    R <: AbstractFloat || throw(
        ArgumentError(
            "$operation optimizer models require Float32, Float64, or compatible floating real components; got $R",
        ),
    )
    scale = maximum(abs, dense; init=one(R))
    tolerance = _symext_tolerance(R, scale; atol, rtol)
    eigenvalues = eigvals(Hermitian(dense))
    minimum_eigenvalue = minimum(real, eigenvalues)
    minimum_eigenvalue < -tolerance && throw(
        ArgumentError(
            "$operation requires a positive semidefinite input; minimum eigenvalue $minimum_eigenvalue is below -$tolerance",
        ),
    )
    return (
        matrix=dense,
        dimensions=dimensions,
        real_type=R,
        tolerance=convert(R, tolerance),
        minimum_eigenvalue=convert(R, minimum_eigenvalue),
        psd_boundary=abs(minimum_eigenvalue) <= tolerance,
        scale=convert(R, scale),
    )
end

function _symext_affine_congruence(
    matrix::HermitianAffineMatrix{T}, factor::AbstractMatrix{<:Number}, name::Symbol
) where {T<:Real}
    Base.require_one_based_indexing(factor)
    size(factor, 2) == matrix.dimension ||
        throw(DimensionMismatch("congruence factor has the wrong column dimension"))
    constant = factor * matrix.constant * adjoint(factor)
    variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    for term in matrix.terms
        coefficient = factor * term.coefficient * adjoint(factor)
        nnz(sparse(coefficient)) == 0 && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name, constant, variables, coefficients, matrix.variable_count
    )
end

function _symext_affine_combination(
    left::HermitianAffineMatrix{T},
    right::HermitianAffineMatrix{T},
    left_weight,
    right_weight,
    name::Symbol,
) where {T<:Real}
    left.dimension == right.dimension ||
        throw(DimensionMismatch("affine matrices must have the same dimension"))
    left.variable_count == right.variable_count ||
        throw(DimensionMismatch("affine matrices must have the same variable count"))
    left_terms = Dict(term.variable => term.coefficient for term in left.terms)
    right_terms = Dict(term.variable => term.coefficient for term in right.terms)
    variables = sort!(collect(union(keys(left_terms), keys(right_terms))))
    kept_variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    zero_matrix = spzeros(Complex{T}, left.dimension, left.dimension)
    for variable in variables
        coefficient =
            left_weight * get(left_terms, variable, zero_matrix) +
            right_weight * get(right_terms, variable, zero_matrix)
        nnz(sparse(coefficient)) == 0 && continue
        push!(kept_variables, variable)
        push!(coefficients, coefficient)
    end
    constant = left_weight * left.constant + right_weight * right.constant
    return HermitianAffineMatrix(
        name, constant, kept_variables, coefficients, left.variable_count
    )
end

function _symext_affine_difference(
    left::HermitianAffineMatrix, right::HermitianAffineMatrix, name::Symbol
)
    return _symext_affine_combination(left, right, 1, -1, name)
end

function _symext_bosonic_lift(
    dimensions::NTuple{2,Int}, order::Int, ::Type{R}
) where {R<:AbstractFloat}
    dimension_a, dimension_b = dimensions
    basis_b = symmetric_subspace_basis(dimension_b, order; T=R, sparse_output=true)
    identity_a = sparse(
        1:dimension_a, 1:dimension_a, fill(one(R), dimension_a), dimension_a, dimension_a
    )
    return kron(identity_a, basis_b)
end

function _symext_preflight(
    dimensions::NTuple{2,Int},
    order::Int,
    bosonic::Bool,
    ppt::Bool,
    limits::OptimizationLimits,
    ;
    inner::Bool=false,
)
    dimension_a, dimension_b = dimensions
    ambient_big = BigInt(dimension_a) * BigInt(dimension_b)^order
    symmetric_dimension_big =
        BigInt(dimension_a) * binomial(BigInt(dimension_b + order - 1), order)
    variable_matrix_big = bosonic ? symmetric_dimension_big : ambient_big
    variables_big = variable_matrix_big^2
    symmetry_equalities_big = bosonic ? BigInt(0) : BigInt(order - 1) * ambient_big^2
    marginal_equalities_big = BigInt(dimension_a * dimension_b)^2
    equalities_big = marginal_equalities_big + symmetry_equalities_big
    ppt_blocks_big = ppt ? BigInt(inner ? 1 : order) : BigInt(0)
    psd_blocks_big = BigInt(1) + ppt_blocks_big
    largest_psd_big = ppt ? max(variable_matrix_big, ambient_big) : variable_matrix_big
    variable_real_block_entries = (2variable_matrix_big) * (2variable_matrix_big + 1) ÷ 2
    ppt_real_block_entries = ppt_blocks_big * (2ambient_big) * (2ambient_big + 1) ÷ 2
    real_block_entries_big = variable_real_block_entries + ppt_real_block_entries
    reasons = String[]
    variables_big > limits.max_variables && push!(
        reasons,
        "model needs $variables_big variables, exceeding max_variables=$(limits.max_variables)",
    )
    equalities_big > limits.max_equalities && push!(
        reasons,
        "model needs $equalities_big equalities, exceeding max_equalities=$(limits.max_equalities)",
    )
    psd_blocks_big > limits.max_psd_blocks && push!(
        reasons,
        "model needs $psd_blocks_big PSD blocks, exceeding max_psd_blocks=$(limits.max_psd_blocks)",
    )
    largest_psd_big > limits.max_psd_dimension && push!(
        reasons,
        "largest PSD block has dimension $largest_psd_big, exceeding max_psd_dimension=$(limits.max_psd_dimension)",
    )
    real_block_entries_big > limits.max_model_entries && push!(
        reasons,
        "real-block PSD scalarization needs $real_block_entries_big triangle entries, exceeding max_model_entries=$(limits.max_model_entries)",
    )
    return (
        allowed=isempty(reasons),
        reason=join(reasons, "; "),
        ambient_dimension=ambient_big,
        symmetric_dimension=symmetric_dimension_big,
        variable_matrix_dimension=variable_matrix_big,
        variable_count=variables_big,
        marginal_equalities=marginal_equalities_big,
        symmetry_equalities=symmetry_equalities_big,
        psd_blocks=psd_blocks_big,
        ppt_blocks=ppt_blocks_big,
        real_block_triangle_entries=real_block_entries_big,
    )
end

function _symext_permuted_affine(
    extension::HermitianAffineMatrix, dimensions::NTuple{2,Int}, order::Int, first_b::Int
)
    full_dimensions = (dimensions[1], ntuple(_ -> dimensions[2], order)...)
    permutation = collect(1:(order + 1))
    left_position = first_b + 1
    right_position = left_position + 1
    permutation[left_position], permutation[right_position] = permutation[right_position],
    permutation[left_position]
    plan = SubsystemPermutationPlan(full_dimensions, Tuple(permutation))
    return _optimization_apply_matrix_map(
        coefficient -> permute_subsystems(coefficient, plan),
        extension,
        Symbol(:extension_swap_, first_b, :_, first_b + 1),
    )
end

function _symext_build_outer_problem(
    prepared, order::Int, ppt::Bool, bosonic::Bool, limits::OptimizationLimits
)
    dimensions = prepared.dimensions
    R = prepared.real_type
    preflight = _symext_preflight(dimensions, order, bosonic, ppt, limits)
    preflight.allowed || return (problem=nothing, preflight=preflight)
    ambient = _symext_checked_int(preflight.ambient_dimension, "ambient_dimension")
    variable_dimension = _symext_checked_int(
        preflight.variable_matrix_dimension, "variable_matrix_dimension"
    )
    variable_count = _symext_checked_int(preflight.variable_count, "variable_count")
    variable = hermitian_variable(
        bosonic ? :bosonic_coordinate : :extension,
        variable_dimension;
        variable_count=variable_count,
        coefficient_type=R,
    )
    basis = bosonic ? _symext_bosonic_lift(dimensions, order, R) : nothing
    extension = if bosonic
        _symext_affine_congruence(variable, basis, :extension)
    else
        variable
    end
    full_dimensions = (dimensions[1], ntuple(_ -> dimensions[2], order)...)
    trace_out = order == 1 ? () : Tuple(3:(order + 1))
    marginal = if isempty(trace_out)
        extension
    else
        partial_trace_affine(
            extension, full_dimensions; trace_out=trace_out, name=:extension_ab_marginal
        )
    end
    equalities = hermitian_equalities(
        marginal; target=prepared.matrix, name_prefix=:extension_marginal
    )
    marginal_equality_count = length(equalities)
    symmetry_equality_count = 0
    if !bosonic
        for first_b in 1:(order - 1)
            permuted = _symext_permuted_affine(extension, dimensions, order, first_b)
            difference = _symext_affine_difference(
                extension, permuted, Symbol(:extension_symmetry_, first_b)
            )
            constraints = hermitian_equalities(
                difference; name_prefix=Symbol(:extension_symmetry_, first_b)
            )
            symmetry_equality_count += length(constraints)
            append!(equalities, constraints)
        end
    end
    psd_constraints = HermitianAffineMatrix[variable]
    ppt_block_count = 0
    if ppt
        for last_system in 2:(order + 1)
            transposed = partial_transpose_affine(
                extension,
                full_dimensions;
                systems=Tuple(2:last_system),
                name=Symbol(:extension_ppt_, last_system - 1),
            )
            push!(psd_constraints, transposed)
            ppt_block_count += 1
        end
    end
    views = HermitianAffineMatrix[extension]
    bosonic && push!(views, variable)
    objective = AffineScalar(zero(R), zeros(R, variable_count))
    metadata = (
        formulation=:dps_symmetric_extension,
        dimensions=dimensions,
        order=order,
        ppt=ppt,
        bosonic=bosonic,
        ambient_dimension=ambient,
        variable_matrix_dimension=variable_dimension,
        marginal_equality_count=marginal_equality_count,
        symmetry_equality_count=symmetry_equality_count,
        ppt_block_count=ppt_block_count,
        subsystem_order=:A_then_B_copies,
        preflight=preflight,
    )
    program = SemidefiniteProgram(
        :symmetric_extension,
        :feasibility,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=psd_constraints,
        primal_views=views,
        limits=limits,
        metadata=metadata,
    )
    problem = SymmetricExtensionProblem{
        R,
        typeof(prepared.matrix),
        typeof(program),
        typeof(extension),
        typeof(basis),
        typeof(metadata),
    }(
        copy(prepared.matrix),
        dimensions,
        order,
        ppt,
        bosonic,
        :dps_symmetric_extension,
        program,
        extension,
        basis,
        marginal_equality_count,
        symmetry_equality_count,
        ppt_block_count,
        ambient,
        variable_dimension,
        prepared.tolerance,
        metadata,
    )
    return (problem=problem, preflight=preflight)
end

# Descending-power coefficients. This private, independently tested recurrence
# supersedes the pinned helpers/jacobi_poly.m without exposing a helper alias.
function _jacobi_polynomial_coefficients(alpha::Real, beta::Real, degree::Integer)
    n = _symext_nonnegative_integer(degree, "degree")
    a, b = promote(alpha, beta)
    if n == 0
        return [one(a + b)]
    elseif n == 1
        two = one(a + b) + one(a + b)
        return [(a + b + two) / two, (a - b) / two]
    end
    previous = _jacobi_polynomial_coefficients(a, b, n - 1)
    earlier = _jacobi_polynomial_coefficients(a, b, n - 2)
    T = promote_type(eltype(previous), eltype(earlier))
    padded_constant = vcat(zero(T), previous)
    padded_linear = vcat(previous, zero(T))
    padded_earlier = vcat(zero(T), zero(T), earlier)
    nT = convert(T, n)
    coefficient_constant = (2nT + a + b - 1) * (a^2 - b^2)
    coefficient_linear = (2nT + a + b - 1) * (2nT + a + b) * (2nT + a + b - 2)
    coefficient_earlier = 2 * (nT + a - 1) * (nT + b - 1) * (2nT + a + b)
    denominator = 2nT * (nT + a + b) * (2nT + a + b - 2)
    iszero(denominator) &&
        throw(ArgumentError("Jacobi recurrence is singular for the supplied parameters"))
    return (
        coefficient_constant * padded_constant + coefficient_linear * padded_linear -
        coefficient_earlier * padded_earlier
    ) / denominator
end

function _jacobi_roots(alpha::R, beta::R, degree::Int) where {R<:AbstractFloat}
    coefficients = R.(_jacobi_polynomial_coefficients(alpha, beta, degree))
    leading = first(coefficients)
    isfinite(leading) && !iszero(leading) ||
        return (roots=nothing, residual=R(Inf), message="invalid leading coefficient")
    n = length(coefficients) - 1
    n == 0 && return (roots=R[], residual=zero(R), message="constant polynomial")
    companion = zeros(R, n, n)
    n > 1 && (companion[2:n, 1:(n - 1)] .= Matrix{R}(I, n - 1, n - 1))
    companion[:, n] .= -reverse(coefficients[2:end]) / leading
    values = eigvals(companion)
    imaginary_residual = maximum(abs ∘ imag, values; init=zero(R))
    root_scale = maximum(abs, values; init=one(R))
    tolerance = 256 * eps(R) * max(one(R), root_scale)
    imaginary_residual <= tolerance || return (
        roots=nothing,
        residual=imaginary_residual,
        message="Jacobi companion roots have non-negligible imaginary parts",
    )
    roots = sort!(R[real(value) for value in values])
    polynomial_residual = zero(R)
    for root in roots
        value = first(coefficients)
        for coefficient in coefficients[2:end]
            value = muladd(value, root, coefficient)
        end
        polynomial_residual = max(polynomial_residual, abs(value))
    end
    return (
        roots=roots,
        residual=max(imaginary_residual, polynomial_residual),
        message="Jacobi roots computed from the reviewed private recurrence",
    )
end

function _symext_inner_mixing_parameter(
    dimension_b::Int, order::Int, ::Type{R}
) where {R<:AbstractFloat}
    dimension_b >= 2 || return (
        value=zero(R),
        residual=zero(R),
        roots=R[],
        message="the one-dimensional B system needs no Jacobi mixing",
    )
    degree = fld(order, 2) + 1
    root_result = _jacobi_roots(
        convert(R, dimension_b - 2), convert(R, mod(order, 2)), degree
    )
    root_result.roots === nothing && return (
        value=nothing,
        residual=root_result.residual,
        roots=nothing,
        message=root_result.message,
    )
    largest_root = maximum(root_result.roots)
    value =
        (one(R) - largest_root) * convert(R, dimension_b) /
        (2 * convert(R, dimension_b - 1))
    isfinite(value) && zero(R) <= value <= one(R) || return (
        value=nothing,
        residual=root_result.residual,
        roots=root_result.roots,
        message="Jacobi mixing parameter lies outside [0,1]",
    )
    return (
        value=value,
        residual=root_result.residual,
        roots=root_result.roots,
        message=root_result.message,
    )
end

function _symext_build_inner_problem(
    prepared, order::Int, ppt::Bool, limits::OptimizationLimits
)
    dimensions = prepared.dimensions
    R = prepared.real_type
    preflight = _symext_preflight(dimensions, order, true, ppt, limits; inner=true)
    preflight.allowed || return (problem=nothing, preflight=preflight, mixing=nothing)
    mixing = if ppt
        _symext_inner_mixing_parameter(dimensions[2], order, R)
    else
        (
            value=zero(R),
            residual=zero(R),
            roots=R[],
            message="non-PPT inner hierarchy uses the rational depolarizing coefficient",
        )
    end
    mixing.value === nothing ||
        isfinite(mixing.value) ||
        return (problem=nothing, preflight=preflight, mixing=mixing)
    mixing.value === nothing && return (problem=nothing, preflight=preflight, mixing=mixing)
    ambient = _symext_checked_int(preflight.ambient_dimension, "ambient_dimension")
    variable_dimension = _symext_checked_int(
        preflight.variable_matrix_dimension, "variable_matrix_dimension"
    )
    variable_count = _symext_checked_int(preflight.variable_count, "variable_count")
    variable = hermitian_variable(
        :inner_bosonic_coordinate,
        variable_dimension;
        variable_count=variable_count,
        coefficient_type=R,
    )
    basis = _symext_bosonic_lift(dimensions, order, R)
    extension = _symext_affine_congruence(variable, basis, :extension)
    full_dimensions = (dimensions[1], ntuple(_ -> dimensions[2], order)...)
    rho_ab = partial_trace_affine(
        extension, full_dimensions; trace_out=Tuple(3:(order + 1)), name=:inner_ab_marginal
    )
    rho_a = partial_trace_affine(
        extension, full_dimensions; trace_out=Tuple(2:(order + 1)), name=:inner_a_marginal
    )
    identity_b = sparse(
        1:dimensions[2],
        1:dimensions[2],
        fill(one(R), dimensions[2]),
        dimensions[2],
        dimensions[2],
    )
    rho_a_identity = tensor_affine(rho_a; right=identity_b, name=:inner_a_identity_b)
    transformed = if ppt
        _symext_affine_combination(
            rho_ab,
            rho_a_identity,
            one(R) - mixing.value,
            mixing.value / dimensions[2],
            :inner_target_map,
        )
    else
        denominator = convert(R, order + dimensions[2])
        _symext_affine_combination(
            rho_ab,
            rho_a_identity,
            convert(R, order) / denominator,
            inv(denominator),
            :inner_target_map,
        )
    end
    equalities = hermitian_equalities(
        transformed; target=prepared.matrix, name_prefix=:inner_marginal
    )
    psd_constraints = HermitianAffineMatrix[variable]
    ppt_block_count = 0
    if ppt
        last_system = cld(order, 2) + 1
        transposed = partial_transpose_affine(
            extension,
            full_dimensions;
            systems=Tuple(1:last_system),
            name=:inner_extension_ppt,
        )
        push!(psd_constraints, transposed)
        ppt_block_count = 1
    end
    objective = AffineScalar(zero(R), zeros(R, variable_count))
    metadata = (
        formulation=:navascues_owari_plenio_inner_extension,
        dimensions=dimensions,
        order=order,
        ppt=ppt,
        bosonic=true,
        mixing_parameter=mixing.value,
        jacobi_roots=mixing.roots,
        jacobi_residual=mixing.residual,
        ambient_dimension=ambient,
        variable_matrix_dimension=variable_dimension,
        marginal_equality_count=length(equalities),
        ppt_block_count=ppt_block_count,
        subsystem_order=:A_then_B_copies,
        preflight=preflight,
    )
    program = SemidefiniteProgram(
        :symmetric_inner_extension,
        :feasibility,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=psd_constraints,
        primal_views=[extension, variable],
        limits=limits,
        metadata=metadata,
    )
    problem = SymmetricInnerExtensionProblem{
        R,
        typeof(prepared.matrix),
        typeof(program),
        typeof(extension),
        typeof(basis),
        typeof(metadata),
    }(
        copy(prepared.matrix),
        dimensions,
        order,
        ppt,
        program,
        extension,
        basis,
        mixing.value,
        length(equalities),
        ppt_block_count,
        ambient,
        variable_dimension,
        prepared.tolerance,
        metadata,
    )
    return (problem=problem, preflight=preflight, mixing=mixing)
end

function _symext_matrix_residual(matrix)
    return maximum(abs, matrix; init=zero(real(eltype(matrix))))
end

function _symext_validation_tolerance(problem, backend)
    if backend isa JuMPBackend
        scale = maximum(abs, problem.state; init=one(problem.tolerance))
        return convert(
            typeof(problem.tolerance), backend.atol + backend.rtol * max(one(scale), scale)
        )
    end
    return problem.tolerance
end

function _symext_validate_extension(
    extension::AbstractMatrix,
    problem::Union{SymmetricExtensionProblem,SymmetricInnerExtensionProblem};
    tolerance,
    inner::Bool,
)
    dimensions = problem.dimensions
    order = problem.order
    full_dimensions = (dimensions[1], ntuple(_ -> dimensions[2], order)...)
    expected_dimension = prod(full_dimensions)
    size(extension) == (expected_dimension, expected_dimension) || return (
        valid=false,
        minimum_eigenvalue=nothing,
        psd_violation=Inf,
        trace_residual=Inf,
        marginal_residual=Inf,
        permutation_residual=Inf,
        bosonic_support_residual=Inf,
        ppt_minimum_eigenvalues=(),
        ppt_violation=Inf,
        hermiticity_residual=Inf,
        message="extension has the wrong matrix dimension",
    )
    dense = Matrix(extension)
    hermiticity_residual = _symext_matrix_residual(dense - adjoint(dense))
    hermiticity_residual <= tolerance || return (
        valid=false,
        minimum_eigenvalue=nothing,
        psd_violation=Inf,
        trace_residual=abs(real(tr(dense)) - real(tr(problem.state))),
        marginal_residual=Inf,
        permutation_residual=Inf,
        bosonic_support_residual=Inf,
        ppt_minimum_eigenvalues=(),
        ppt_violation=Inf,
        hermiticity_residual=hermiticity_residual,
        message="extension is non-Hermitian outside tolerance",
    )
    minimum_eigenvalue = eigmin(Hermitian(dense))
    psd_violation = max(zero(minimum_eigenvalue), -minimum_eigenvalue)
    trace_residual = abs(real(tr(dense)) - real(tr(problem.state)))
    trace_out = order == 1 ? () : Tuple(3:(order + 1))
    marginal = if isempty(trace_out)
        dense
    else
        partial_trace(dense, full_dimensions; trace_out=trace_out)
    end
    if inner
        rho_a = partial_trace(dense, full_dimensions; trace_out=Tuple(2:(order + 1)))
        identity_b = Matrix{eltype(dense)}(I, dimensions[2], dimensions[2])
        transformed = if problem.ppt
            (one(problem.mixing_parameter) - problem.mixing_parameter) * marginal +
            problem.mixing_parameter * tensor_product(rho_a, identity_b) / dimensions[2]
        else
            (order * marginal + tensor_product(rho_a, identity_b)) / (order + dimensions[2])
        end
        marginal_residual = _symext_matrix_residual(transformed - problem.state)
    else
        marginal_residual = _symext_matrix_residual(marginal - problem.state)
    end
    permutation_residual = zero(real(eltype(dense)))
    for first_b in 1:(order - 1)
        permutation = collect(1:(order + 1))
        left_position = first_b + 1
        right_position = left_position + 1
        permutation[left_position], permutation[right_position] = permutation[right_position],
        permutation[left_position]
        permuted = permute_subsystems(
            dense, full_dimensions; permutation=Tuple(permutation)
        )
        permutation_residual = max(
            permutation_residual, _symext_matrix_residual(permuted - dense)
        )
    end
    bosonic = inner || problem.bosonic
    bosonic_support_residual = if bosonic
        basis = problem.bosonic_basis
        projector = basis * adjoint(basis)
        _symext_matrix_residual(projector * dense * projector - dense)
    else
        zero(permutation_residual)
    end
    ppt_minimum_eigenvalues = if problem.ppt
        if inner
            last_system = cld(order, 2) + 1
            transposed = partial_transpose(
                dense, full_dimensions; systems=Tuple(1:last_system)
            )
            (eigmin(Hermitian(Matrix(transposed))),)
        else
            Tuple(
                eigmin(
                    Hermitian(
                        Matrix(
                            partial_transpose(
                                dense, full_dimensions; systems=Tuple(2:last_system)
                            ),
                        ),
                    ),
                ) for last_system in 2:(order + 1)
            )
        end
    else
        ()
    end
    ppt_violation = maximum(
        (max(zero(value), -value) for value in ppt_minimum_eigenvalues);
        init=zero(minimum_eigenvalue),
    )
    valid =
        maximum((
            psd_violation,
            trace_residual,
            marginal_residual,
            permutation_residual,
            bosonic_support_residual,
            ppt_violation,
            hermiticity_residual,
        ),) <= tolerance
    return (
        valid=valid,
        minimum_eigenvalue=minimum_eigenvalue,
        psd_violation=psd_violation,
        trace_residual=trace_residual,
        marginal_residual=marginal_residual,
        permutation_residual=permutation_residual,
        bosonic_support_residual=bosonic_support_residual,
        ppt_minimum_eigenvalues=ppt_minimum_eigenvalues,
        ppt_violation=ppt_violation,
        hermiticity_residual=hermiticity_residual,
        message=if valid
            "extension passed independent primal checks"
        else
            "one or more extension residuals exceed tolerance"
        end,
    )
end

function _symext_hermitian_from_duals(values::AbstractVector{T}, dimension::Int) where {T}
    length(values) == dimension^2 ||
        throw(DimensionMismatch("wrong number of Hermitian equality duals"))
    matrix = zeros(Complex{T}, dimension, dimension)
    position = 1
    for index in 1:dimension
        matrix[index, index] = values[position]
        position += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        value = values[position] / 2
        matrix[row, column] += value
        matrix[column, row] += value
        position += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        value = values[position] / 2
        matrix[row, column] += im * value
        matrix[column, row] -= im * value
        position += 1
    end
    return matrix
end

function _symext_dual_witness(problem, optimization, tolerance; inner::Bool)
    optimization.dual === nothing && return nothing
    count = problem.marginal_equality_count
    length(optimization.dual.equalities) >= count || return nothing
    dimension = prod(problem.dimensions)
    equality_dual = _symext_hermitian_from_duals(
        view(optimization.dual.equalities, 1:count), dimension
    )
    # JuMP/MOI's equality-dual convention is opposite to the separating
    # operator convention used by the DPS cone: the witness is `-W`, as in
    # the pinned routine. The full, unmodified dual tuple remains attached
    # below and its stationarity/cone residual is checked independently.
    raw = -equality_dual
    expectation = real(dot(raw, problem.state))
    expectation < -tolerance || return SymmetricExtensionWitness(
        raw,
        expectation,
        abs(expectation + one(expectation)),
        optimization.dual_residual,
        tolerance,
        optimization.dual,
        false,
        false,
        if inner
            "the dual object is not automatically an entanglement witness"
        else
            "the dual orientation did not produce a robust negative expectation"
        end,
    )
    normalized = raw / (-expectation)
    normalized_expectation = real(dot(normalized, problem.state))
    dual_valid =
        optimization.dual_residual !== nothing &&
        optimization.dual_residual <= tolerance &&
        abs(normalized_expectation + one(normalized_expectation)) <= 8tolerance
    return SymmetricExtensionWitness(
        normalized,
        normalized_expectation,
        abs(normalized_expectation + one(normalized_expectation)),
        optimization.dual_residual,
        tolerance,
        optimization.dual,
        dual_valid,
        !inner && dual_valid,
        if inner
            "a negative dual separator for the inner cone is not automatically an entanglement witness"
        else
            "the outer-hierarchy dual separator contains all separable operators in its nonnegative cone"
        end,
    )
end

function _symext_result(
    status,
    verdict,
    certificate_kind,
    hierarchy,
    dimensions,
    order,
    ppt,
    bosonic;
    problem=nothing,
    extension=nothing,
    witness=nothing,
    optimization_result=nothing,
    residuals=nothing,
    tolerance,
    warnings=(),
    message,
)
    return SymmetricExtensionResult(
        status,
        verdict,
        certificate_kind,
        hierarchy,
        dimensions,
        order,
        ppt,
        bosonic,
        problem,
        extension,
        witness,
        optimization_result,
        residuals,
        tolerance,
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

function _symext_backend_status(optimization)
    optimization.status === OptimizationBackendUnavailable &&
        return SymmetricExtensionBackendUnavailable
    optimization.status === OptimizationLimit && return SymmetricExtensionResourceLimit
    optimization.status in (
        OptimizationMalformedBackend,
        OptimizationNumericalFailure,
        OptimizationUnsupported,
        OptimizationUnbounded,
        OptimizationInconsistent,
        OptimizationUnknown,
    ) && return SymmetricExtensionBackendFailure
    return SymmetricExtensionBackendFailure
end

function _symext_is_model_limit_error(error)
    error isa ArgumentError || return false
    message = sprint(showerror, error)
    return occursin("exceeding max_", message) || occursin("exceeds max_", message)
end

function _symext_solve_problem(problem, backend; inner::Bool)
    optimization = solve_optimization(problem.program, backend)
    tolerance = _symext_validation_tolerance(problem, backend)
    hierarchy = inner ? :inner : :outer
    bosonic = inner ? true : problem.bosonic
    if optimization.status in (OptimizationOptimal, OptimizationFeasible) &&
        optimization.primal !== nothing &&
        hasproperty(optimization.primal.views, :extension)
        extension = Matrix(optimization.primal.views.extension)
        residuals = _symext_validate_extension(
            extension, problem; tolerance=tolerance, inner=inner
        )
        if optimization.status === OptimizationOptimal && residuals.valid
            return _symext_result(
                SymmetricExtensionSolverPresent,
                true,
                :validated_primal_extension,
                hierarchy,
                problem.dimensions,
                problem.order,
                problem.ppt,
                bosonic;
                problem=problem,
                extension=extension,
                optimization_result=optimization,
                residuals=residuals,
                tolerance=tolerance,
                warnings=optimization.warnings,
                message="the optimizer result and independently reconstructed extension satisfy every requested constraint",
            )
        end
        return _symext_result(
            if residuals.valid
                SymmetricExtensionNumericalBoundary
            else
                SymmetricExtensionInvalidCertificate
            end,
            nothing,
            nothing,
            hierarchy,
            problem.dimensions,
            problem.order,
            problem.ppt,
            bosonic;
            problem=problem,
            extension=extension,
            optimization_result=optimization,
            residuals=residuals,
            tolerance=tolerance,
            warnings=(
                optimization.warnings...,
                "a feasible/inaccurate or residual-failing primal candidate is retained but is not promoted to a certificate",
            ),
            message=residuals.message,
        )
    elseif optimization.status === OptimizationInfeasible
        witness = _symext_dual_witness(problem, optimization, tolerance; inner=inner)
        if witness !== nothing && witness.separator_validated
            return _symext_result(
                SymmetricExtensionSolverAbsent,
                false,
                if inner
                    :validated_inner_cone_separator
                else
                    :validated_symmetric_extension_separator
                end,
                hierarchy,
                problem.dimensions,
                problem.order,
                problem.ppt,
                bosonic;
                problem=problem,
                witness=witness,
                optimization_result=optimization,
                tolerance=tolerance,
                warnings=inner ? (witness.warning,) : (),
                message=if inner
                    "a checked dual certificate separates the input from the inner approximation; this alone does not prove entanglement"
                else
                    "a checked dual certificate separates the input from the symmetric-extension cone"
                end,
            )
        end
        return _symext_result(
            SymmetricExtensionInvalidCertificate,
            nothing,
            nothing,
            hierarchy,
            problem.dimensions,
            problem.order,
            problem.ppt,
            bosonic;
            problem=problem,
            witness=witness,
            optimization_result=optimization,
            tolerance=tolerance,
            warnings=(
                "the optimizer reported infeasibility but no independently validated normalized dual separator was available",
            ),
            message="solver infeasibility is retained as evidence but not promoted to a mathematical verdict",
        )
    end
    status = _symext_backend_status(optimization)
    return _symext_result(
        status,
        nothing,
        nothing,
        hierarchy,
        problem.dimensions,
        problem.order,
        problem.ppt,
        bosonic;
        problem=problem,
        optimization_result=optimization,
        tolerance=tolerance,
        warnings=optimization.warnings,
        message=optimization.message,
    )
end

function _symext_ppt_shortcut(prepared, order::Int, bosonic::Bool)
    state = prepared.matrix
    dimensions = prepared.dimensions
    transposed = partial_transpose(state, dimensions; systems=(2,))
    decomposition = eigen(Hermitian(Matrix(transposed)))
    index = argmin(decomposition.values)
    minimum_eigenvalue = decomposition.values[index]
    tolerance = prepared.tolerance
    if minimum_eigenvalue < -tolerance
        vector = decomposition.vectors[:, index]
        raw_witness = partial_transpose(vector * adjoint(vector), dimensions; systems=(2,))
        expectation = real(dot(raw_witness, state))
        witness = if expectation < -tolerance
            normalized = raw_witness / (-expectation)
            SymmetricExtensionWitness(
                normalized,
                real(dot(normalized, state)),
                abs(real(dot(normalized, state)) + 1),
                nothing,
                tolerance,
                nothing,
                true,
                true,
                "the decomposable NPT witness is nonnegative on every PPT extension cone",
            )
        else
            nothing
        end
        return _symext_result(
            SymmetricExtensionAnalyticAbsent,
            false,
            :negative_partial_transpose,
            :outer,
            dimensions,
            order,
            true,
            bosonic;
            witness=witness,
            residuals=(
                partial_transpose_minimum_eigenvalue=minimum_eigenvalue,
                partial_transpose_violation=(-minimum_eigenvalue),
            ),
            tolerance=tolerance,
            message="negative partial transpose excludes every requested PPT symmetric extension",
        )
    elseif minimum_eigenvalue > tolerance
        order > 1 && prod(dimensions) > 6 && return nothing
        one_copy = order == 1
        residuals = if one_copy
            (
                valid=true,
                minimum_eigenvalue=prepared.minimum_eigenvalue,
                psd_violation=max(
                    zero(prepared.minimum_eigenvalue), -prepared.minimum_eigenvalue
                ),
                trace_residual=zero(tolerance),
                marginal_residual=zero(tolerance),
                permutation_residual=zero(tolerance),
                bosonic_support_residual=zero(tolerance),
                ppt_minimum_eigenvalues=(minimum_eigenvalue,),
                ppt_violation=zero(minimum_eigenvalue),
                hermiticity_residual=zero(tolerance),
                message="the input passed the one-copy PSD and PPT checks",
            )
        else
            (
                partial_transpose_minimum_eigenvalue=minimum_eigenvalue,
                partial_transpose_violation=zero(minimum_eigenvalue),
            )
        end
        return _symext_result(
            SymmetricExtensionAnalyticPresent,
            true,
            if one_copy
                :one_copy_ppt_extension
            else
                :low_dimensional_ppt_separability_theorem
            end,
            :outer,
            dimensions,
            order,
            true,
            bosonic;
            extension=one_copy ? copy(state) : nothing,
            residuals=residuals,
            tolerance=tolerance,
            message=if one_copy
                "the input itself is a robustly PPT one-copy extension"
            else
                "for total dimension at most six, robust PPT implies separability and therefore every requested symmetric extension exists"
            end,
        )
    end
    return _symext_result(
        SymmetricExtensionNumericalBoundary,
        nothing,
        nothing,
        :outer,
        dimensions,
        order,
        true,
        bosonic;
        residuals=(
            partial_transpose_minimum_eigenvalue=minimum_eigenvalue,
            partial_transpose_violation=max(zero(minimum_eigenvalue), -minimum_eigenvalue),
        ),
        tolerance=tolerance,
        message=if order == 1
            "the one-copy partial-transpose test lies inside the numerical boundary"
        else
            "the low-dimensional PPT shortcut lies inside the numerical boundary"
        end,
    )
end

function _symext_two_qubit_shortcut(prepared, bosonic::Bool)
    state = prepared.matrix
    dimensions = prepared.dimensions
    rho_b = partial_trace(state, dimensions; trace_out=(1,))
    determinant = det(state)
    determinant_imaginary = abs(imag(determinant))
    determinant_real = real(determinant)
    tolerance = prepared.tolerance
    if determinant_imaginary > tolerance || determinant_real < -tolerance
        return _symext_result(
            SymmetricExtensionNumericalBoundary,
            nothing,
            nothing,
            :outer,
            dimensions,
            2,
            false,
            bosonic;
            residuals=(
                determinant=determinant,
                determinant_imaginary_residual=determinant_imaginary,
            ),
            tolerance=tolerance,
            message="the determinant needed by the two-qubit theorem is numerically inconsistent with a PSD Hermitian state",
        )
    elseif determinant_real < zero(determinant_real)
        return _symext_result(
            SymmetricExtensionNumericalBoundary,
            nothing,
            nothing,
            :outer,
            dimensions,
            2,
            false,
            bosonic;
            residuals=(determinant=determinant,),
            tolerance=tolerance,
            message="the determinant lies just below zero; it is not clipped before applying the analytic theorem",
        )
    end
    left = real(tr(rho_b * rho_b))
    right = real(tr(state * state)) - 4sqrt(determinant_real)
    margin = left - right
    scale = max(one(margin), abs(left), abs(right))
    theorem_tolerance = max(tolerance, sqrt(eps(typeof(margin))) * scale)
    if margin > theorem_tolerance
        return _symext_result(
            SymmetricExtensionAnalyticPresent,
            true,
            :two_qubit_symmetric_extension_theorem,
            :outer,
            dimensions,
            2,
            false,
            bosonic;
            residuals=(left=left, right=right, margin=margin),
            tolerance=theorem_tolerance,
            message="the Chen--Ji--Kribs--Lütkenhaus--Zeng inequality holds with margin",
        )
    elseif margin < -theorem_tolerance
        return _symext_result(
            SymmetricExtensionAnalyticAbsent,
            false,
            :two_qubit_symmetric_extension_theorem,
            :outer,
            dimensions,
            2,
            false,
            bosonic;
            residuals=(left=left, right=right, margin=margin),
            tolerance=theorem_tolerance,
            message="the two-qubit symmetric-extension inequality is violated with margin",
        )
    end
    return _symext_result(
        SymmetricExtensionNumericalBoundary,
        nothing,
        nothing,
        :outer,
        dimensions,
        2,
        false,
        bosonic;
        residuals=(left=left, right=right, margin=margin),
        tolerance=theorem_tolerance,
        message="the two-qubit symmetric-extension inequality lies inside the numerical boundary",
    )
end

"""
    symmetric_extension_problem(state; order=2, dims=nothing, ppt=false,
                                bosonic=false, ...)

Build the solver-independent DPS symmetric-extension SDP. This low-level
builder throws on malformed input or exceeded limits. The high-level
[`symmetric_extension`](@ref) converts work-limit failures into structured
inconclusive results.
"""
function symmetric_extension_problem(
    state::AbstractMatrix{<:Number};
    order=2,
    dims=nothing,
    ppt::Bool=false,
    bosonic::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _symext_positive_integer(order, "order")
    checked_max_order = _symext_positive_integer(max_order, "max_order")
    checked_order <= checked_max_order ||
        throw(ArgumentError("order=$checked_order exceeds max_order=$checked_max_order"))
    prepared = _symext_prepare_operator(
        state,
        dims;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        operation="symmetric_extension_problem",
    )
    prepared.minimum_eigenvalue < zero(prepared.minimum_eigenvalue) && throw(
        ArgumentError(
            "the input minimum eigenvalue $(prepared.minimum_eigenvalue) lies on the negative PSD tolerance boundary; use symmetric_extension to retain an inconclusive structured result",
        ),
    )
    built = _symext_build_outer_problem(prepared, checked_order, ppt, bosonic, limits)
    built.problem === nothing && throw(ArgumentError(built.preflight.reason))
    return built.problem
end

"""
    symmetric_extension(state; order=2, dims=nothing, ppt=false,
                        bosonic=false, backend=NoOptimizationBackend(), ...)

Analyze whether `state` has a permutation-invariant `k=order` extension on
subsystem `B`. An explicit backend is required outside theorem branches.
General models impose the marginal, adjacent-copy permutation symmetry, PSD,
and requested representative PPT constraints directly. Bosonic models use an
orthonormal occupation basis and retain the ambient extension.

The native result never converts backend failure or a tolerance boundary to a
Boolean. `OptimizationLimits`, `max_order`, and `max_dense_entries` are checked
before the corresponding large allocations.
"""
function symmetric_extension(
    state::AbstractMatrix{<:Number};
    order=2,
    dims=nothing,
    ppt::Bool=false,
    bosonic::Bool=false,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    prefer_analytic::Bool=true,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _symext_positive_integer(order, "order")
    checked_max_order = _symext_positive_integer(max_order, "max_order")
    prepared = _symext_prepare_operator(
        state,
        dims;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        operation="symmetric_extension",
    )
    if checked_order > checked_max_order
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :outer,
            prepared.dimensions,
            checked_order,
            ppt,
            bosonic;
            tolerance=prepared.tolerance,
            message="order=$checked_order exceeds max_order=$checked_max_order",
        )
    end
    if prepared.minimum_eigenvalue < zero(prepared.minimum_eigenvalue)
        return _symext_result(
            SymmetricExtensionNumericalBoundary,
            nothing,
            nothing,
            :outer,
            prepared.dimensions,
            checked_order,
            ppt,
            bosonic;
            residuals=(
                minimum_eigenvalue=prepared.minimum_eigenvalue,
                psd_violation=(-prepared.minimum_eigenvalue),
            ),
            tolerance=prepared.tolerance,
            message="the input has a negative eigenvalue inside the PSD tolerance boundary; it is not repaired or promoted to an extension verdict",
        )
    end
    if checked_order == 1
        if ppt
            return _symext_ppt_shortcut(prepared, 1, bosonic)
        end
        residuals = (
            valid=true,
            minimum_eigenvalue=prepared.minimum_eigenvalue,
            psd_violation=max(
                zero(prepared.minimum_eigenvalue), -prepared.minimum_eigenvalue
            ),
            trace_residual=zero(prepared.tolerance),
            marginal_residual=zero(prepared.tolerance),
            permutation_residual=zero(prepared.tolerance),
            bosonic_support_residual=zero(prepared.tolerance),
            ppt_minimum_eigenvalues=(),
            ppt_violation=zero(prepared.tolerance),
            hermiticity_residual=zero(prepared.tolerance),
            message="the one-copy extension is the input itself",
        )
        return _symext_result(
            SymmetricExtensionExactPresent,
            true,
            :one_copy_identity_extension,
            :outer,
            prepared.dimensions,
            1,
            false,
            bosonic;
            extension=copy(prepared.matrix),
            residuals=residuals,
            tolerance=prepared.tolerance,
            message="a PSD input is its own one-copy symmetric extension",
        )
    end
    no_backend = backend isa NoOptimizationBackend
    if prefer_analytic && no_backend && ppt
        shortcut = _symext_ppt_shortcut(prepared, checked_order, bosonic)
        shortcut === nothing || return shortcut
    elseif prefer_analytic &&
        no_backend &&
        checked_order == 2 &&
        !ppt &&
        prepared.dimensions == (2, 2)
        return _symext_two_qubit_shortcut(prepared, bosonic)
    end
    preflight = _symext_preflight(prepared.dimensions, checked_order, bosonic, ppt, limits)
    if !preflight.allowed
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :outer,
            prepared.dimensions,
            checked_order,
            ppt,
            bosonic;
            residuals=(preflight=preflight,),
            tolerance=prepared.tolerance,
            message=preflight.reason,
        )
    end
    built = try
        _symext_build_outer_problem(prepared, checked_order, ppt, bosonic, limits)
    catch error
        _symext_is_model_limit_error(error) || rethrow()
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :outer,
            prepared.dimensions,
            checked_order,
            ppt,
            bosonic;
            tolerance=prepared.tolerance,
            message=sprint(showerror, error),
        )
    end
    return _symext_solve_problem(built.problem, backend; inner=false)
end

"""
    symmetric_inner_extension_problem(state; order=2, dims=nothing,
                                      ppt=false, ...)

Build the solver-independent Navascués--Owari--Plenio inner-hierarchy SDP.
The private Jacobi recurrence is used only for the PPT mixing coefficient.
"""
function symmetric_inner_extension_problem(
    state::AbstractMatrix{<:Number};
    order=2,
    dims=nothing,
    ppt::Bool=false,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _symext_positive_integer(order, "order")
    checked_order >= 2 ||
        throw(ArgumentError("symmetric inner extensions require order >= 2"))
    checked_max_order = _symext_positive_integer(max_order, "max_order")
    checked_order <= checked_max_order ||
        throw(ArgumentError("order=$checked_order exceeds max_order=$checked_max_order"))
    prepared = _symext_prepare_operator(
        state,
        dims;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        operation="symmetric_inner_extension_problem",
    )
    prepared.minimum_eigenvalue < zero(prepared.minimum_eigenvalue) && throw(
        ArgumentError(
            "the input minimum eigenvalue $(prepared.minimum_eigenvalue) lies on the negative PSD tolerance boundary; use symmetric_inner_extension to retain an inconclusive structured result",
        ),
    )
    built = _symext_build_inner_problem(prepared, checked_order, ppt, limits)
    built.mixing !== nothing &&
        built.mixing.value === nothing &&
        throw(ArgumentError(built.mixing.message))
    built.problem === nothing && throw(ArgumentError(built.preflight.reason))
    return built.problem
end

"""
    symmetric_inner_extension(state; order=2, dims=nothing, ppt=false,
                              backend=NoOptimizationBackend(), ...)

Test membership in the inner approximation to the separable cone from
Navascués, Owari, and Plenio. A positive result retains a checked bosonic
extension. A negative dual separator proves only nonmembership in this inner
cone: it is deliberately marked `entanglement_witness=false`.
"""
function symmetric_inner_extension(
    state::AbstractMatrix{<:Number};
    order=2,
    dims=nothing,
    ppt::Bool=false,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_order=8,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _symext_positive_integer(order, "order")
    checked_order >= 2 ||
        throw(ArgumentError("symmetric inner extensions require order >= 2"))
    checked_max_order = _symext_positive_integer(max_order, "max_order")
    prepared = _symext_prepare_operator(
        state,
        dims;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        operation="symmetric_inner_extension",
    )
    if checked_order > checked_max_order
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            tolerance=prepared.tolerance,
            warnings=(
                "a negative inner-hierarchy result would not automatically be an entanglement witness",
            ),
            message="order=$checked_order exceeds max_order=$checked_max_order",
        )
    end
    if prepared.minimum_eigenvalue < zero(prepared.minimum_eigenvalue)
        return _symext_result(
            SymmetricExtensionNumericalBoundary,
            nothing,
            nothing,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            residuals=(
                minimum_eigenvalue=prepared.minimum_eigenvalue,
                psd_violation=(-prepared.minimum_eigenvalue),
            ),
            tolerance=prepared.tolerance,
            warnings=(
                "a negative inner-hierarchy result would not automatically be an entanglement witness",
            ),
            message="the input has a negative eigenvalue inside the PSD tolerance boundary; it is not repaired or promoted to an inner-cone verdict",
        )
    end
    if iszero(norm(prepared.matrix))
        ambient_big =
            BigInt(prepared.dimensions[1]) * BigInt(prepared.dimensions[2])^checked_order
        entry_limit = _symext_limit(max_dense_entries, "max_dense_entries")
        extension_entries = ambient_big^2
        if entry_limit !== nothing && extension_entries > entry_limit
            return _symext_result(
                SymmetricExtensionResourceLimit,
                nothing,
                nothing,
                :inner,
                prepared.dimensions,
                checked_order,
                ppt,
                true;
                residuals=(
                    ambient_dimension=ambient_big, extension_entries=extension_entries
                ),
                tolerance=prepared.tolerance,
                message="the explicit zero extension needs $extension_entries dense entries, exceeding max_dense_entries=$entry_limit",
            )
        end
        ambient = _symext_checked_int(ambient_big, "ambient_dimension")
        extension = zeros(eltype(prepared.matrix), ambient, ambient)
        return _symext_result(
            SymmetricExtensionExactPresent,
            true,
            :zero_cone_element,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            extension=extension,
            residuals=(
                valid=true,
                minimum_eigenvalue=zero(prepared.tolerance),
                psd_violation=zero(prepared.tolerance),
                trace_residual=zero(prepared.tolerance),
                marginal_residual=zero(prepared.tolerance),
                permutation_residual=zero(prepared.tolerance),
                bosonic_support_residual=zero(prepared.tolerance),
                ppt_minimum_eigenvalues=ppt ? (zero(prepared.tolerance),) : (),
                ppt_violation=zero(prepared.tolerance),
                hermiticity_residual=zero(prepared.tolerance),
                message="the zero operator is in every inner cone",
            ),
            tolerance=prepared.tolerance,
            message="the zero operator has the zero bosonic extension",
        )
    end
    preflight = _symext_preflight(
        prepared.dimensions, checked_order, true, ppt, limits; inner=true
    )
    if !preflight.allowed
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            residuals=(preflight=preflight,),
            tolerance=prepared.tolerance,
            warnings=(
                "a negative inner-hierarchy result would not automatically be an entanglement witness",
            ),
            message=preflight.reason,
        )
    end
    built = try
        _symext_build_inner_problem(prepared, checked_order, ppt, limits)
    catch error
        _symext_is_model_limit_error(error) || rethrow()
        return _symext_result(
            SymmetricExtensionResourceLimit,
            nothing,
            nothing,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            tolerance=prepared.tolerance,
            warnings=(
                "a negative inner-hierarchy result would not automatically be an entanglement witness",
            ),
            message=sprint(showerror, error),
        )
    end
    if built.mixing !== nothing && built.mixing.value === nothing
        return _symext_result(
            SymmetricExtensionBackendFailure,
            nothing,
            nothing,
            :inner,
            prepared.dimensions,
            checked_order,
            ppt,
            true;
            residuals=(jacobi=built.mixing,),
            tolerance=prepared.tolerance,
            warnings=(
                "the private Jacobi-root calculation failed before a model was built",
            ),
            message=built.mixing.message,
        )
    end
    return _symext_solve_problem(built.problem, backend; inner=true)
end

function _random_ppt_dimensions(dims)
    if dims isa Integer
        dimension = _symext_positive_integer(dims, "dims")
        return (dimension, dimension)
    elseif dims isa Union{Tuple,AbstractVector}
        dims isa AbstractVector && Base.require_one_based_indexing(dims)
        length(dims) == 2 ||
            throw(ArgumentError("dims must be a scalar or a two-entry tuple/vector"))
        return (
            _symext_positive_integer(dims[1], "dims[1]"),
            _symext_positive_integer(dims[2], "dims[2]"),
        )
    end
    return throw(
        ArgumentError("dims must be a positive integer or a two-entry tuple/vector")
    )
end

function _random_ppt_ranks(ranks, dimension::Int)
    values = if ranks === nothing
        (dimension, dimension)
    elseif ranks isa Integer
        checked = _symext_positive_integer(ranks, "ranks")
        (checked, checked)
    elseif ranks isa Union{Tuple,AbstractVector}
        ranks isa AbstractVector && Base.require_one_based_indexing(ranks)
        length(ranks) == 2 ||
            throw(ArgumentError("ranks must be a scalar or a two-entry tuple/vector"))
        (
            _symext_positive_integer(ranks[1], "ranks[1]"),
            _symext_positive_integer(ranks[2], "ranks[2]"),
        )
    else
        throw(ArgumentError("ranks must be nothing, a scalar, or a two-entry tuple/vector"))
    end
    all(value -> value <= dimension, values) || throw(
        ArgumentError("each requested rank must not exceed total dimension $dimension")
    )
    return values
end

function _random_ppt_type(::Type{R}) where {R}
    R in (Float32, Float64) || throw(ArgumentError("T must be Float32 or Float64"))
    return R
end

function _random_ppt_construction(value, ranks, dimension)
    value isa Symbol || throw(
        ArgumentError(
            "construction must be :auto, :shifted_induced, or :separable_mixture"
        ),
    )
    value in (:auto, :shifted_induced, :separable_mixture) || throw(
        ArgumentError(
            "construction must be :auto, :shifted_induced, or :separable_mixture"
        ),
    )
    if value === :auto
        return minimum(ranks) == dimension ? :shifted_induced : :separable_mixture
    elseif value === :shifted_induced && minimum(ranks) < dimension
        throw(
            ArgumentError(
                ":shifted_induced is full-rank and cannot satisfy requested low-rank bounds; use :separable_mixture",
            ),
        )
    end
    return value
end

function _random_ppt_failure(
    status,
    dimensions,
    ranks,
    construction,
    real_output,
    max_iterations,
    tolerance,
    work_used,
    max_work,
    max_dense_entries,
    message;
    candidate=nothing,
    history=typeof(tolerance)[],
    iterations=0,
)
    return RandomPPTStateResult(
        status,
        nothing,
        candidate,
        dimensions,
        ranks,
        nothing,
        construction,
        real_output,
        iterations,
        max_iterations,
        history,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        tolerance,
        work_used,
        max_work,
        max_dense_entries,
        false,
        false,
        String(message),
    )
end

function _random_ppt_draw_ginibre(
    rng::AbstractRNG, ::Type{R}, rows::Int, columns::Int, real_output::Bool
) where {R<:AbstractFloat}
    real_part = randn(rng, R, rows, columns)
    return real_output ? real_part : complex.(real_part, randn(rng, R, rows, columns))
end

function _random_ppt_shifted_induced(
    rng::AbstractRNG, dimensions, ::Type{R}, real_output, tolerance
) where {R<:AbstractFloat}
    dimension = prod(dimensions)
    ginibre = _random_ppt_draw_ginibre(rng, R, dimension, dimension, real_output)
    all(isfinite, ginibre) || return (
        candidate=nothing, normalized=false, message="the RNG produced nonfinite data"
    )
    base = ginibre * adjoint(ginibre)
    base_trace = real(tr(base))
    isfinite(base_trace) && base_trace > zero(R) || return (
        candidate=nothing,
        normalized=false,
        message="the random Gram matrix has invalid trace",
    )
    base /= base_trace
    transposed = partial_transpose(base, dimensions; systems=(2,))
    minimum_pt = eigmin(Hermitian(Matrix(transposed)))
    margin = max(convert(R, 8) * tolerance, convert(R, 64dimension) * eps(R))
    shift = max(zero(R), margin - minimum_pt)
    candidate = base + shift * Matrix{eltype(base)}(I, dimension, dimension)
    normalization = real(tr(candidate))
    isfinite(normalization) && normalization > zero(R) || return (
        candidate=candidate,
        normalized=false,
        message="the shifted state has invalid trace",
    )
    candidate /= normalization
    return (
        candidate=candidate,
        normalized=true,
        message="a random induced state was shifted by an explicit identity multiple before normalization",
    )
end

function _random_ppt_separable_mixture(
    rng::AbstractRNG, dimensions, ranks, ::Type{R}, real_output
) where {R<:AbstractFloat}
    component_count = min(ranks...)
    weights = randexp(rng, R, component_count)
    weight_sum = sum(weights)
    isfinite(weight_sum) && weight_sum > zero(R) || return (
        candidate=nothing, normalized=false, message="the RNG produced invalid weights"
    )
    weights /= weight_sum
    element_type = real_output ? R : Complex{R}
    dimension = prod(dimensions)
    candidate = zeros(element_type, dimension, dimension)
    for component in 1:component_count
        left = vec(_random_ppt_draw_ginibre(rng, R, dimensions[1], 1, real_output))
        right = vec(_random_ppt_draw_ginibre(rng, R, dimensions[2], 1, real_output))
        left_norm = norm(left)
        right_norm = norm(right)
        isfinite(left_norm) && left_norm > zero(R) || return (
            candidate=candidate,
            normalized=false,
            message="the RNG produced an invalid local factor",
        )
        isfinite(right_norm) && right_norm > zero(R) || return (
            candidate=candidate,
            normalized=false,
            message="the RNG produced an invalid local factor",
        )
        left /= left_norm
        right /= right_norm
        product_vector = tensor_product(left, right)
        candidate .+= weights[component] .* (product_vector * adjoint(product_vector))
    end
    return (
        candidate=candidate,
        normalized=true,
        message="a bounded random convex mixture of product pure states was constructed",
    )
end

function _random_ppt_validate(candidate, dimensions, ranks, tolerance)
    hermiticity_residual = _symext_matrix_residual(candidate - adjoint(candidate))
    transposed = partial_transpose(candidate, dimensions; systems=(2,))
    partial_transpose_hermiticity_residual = _symext_matrix_residual(
        transposed - adjoint(transposed)
    )
    if max(hermiticity_residual, partial_transpose_hermiticity_residual) > tolerance
        return (
            valid=false,
            ranks=nothing,
            minimum_eigenvalue=nothing,
            minimum_partial_transpose_eigenvalue=nothing,
            trace_residual=abs(real(tr(candidate)) - one(tolerance)),
            hermiticity_residual=hermiticity_residual,
            partial_transpose_hermiticity_residual=partial_transpose_hermiticity_residual,
            message="candidate is non-Hermitian outside tolerance",
        )
    end
    values = eigvals(Hermitian(Matrix(candidate)))
    pt_values = eigvals(Hermitian(Matrix(transposed)))
    minimum_eigenvalue = minimum(values)
    minimum_pt = minimum(pt_values)
    trace_residual = abs(real(tr(candidate)) - one(tolerance))
    numerical_ranks = (
        count(value -> value > tolerance, values),
        count(value -> value > tolerance, pt_values),
    )
    valid =
        minimum_eigenvalue >= -tolerance &&
        minimum_pt >= -tolerance &&
        trace_residual <= tolerance &&
        numerical_ranks[1] <= ranks[1] &&
        numerical_ranks[2] <= ranks[2]
    return (
        valid=valid,
        ranks=numerical_ranks,
        minimum_eigenvalue=minimum_eigenvalue,
        minimum_partial_transpose_eigenvalue=minimum_pt,
        trace_residual=trace_residual,
        hermiticity_residual=hermiticity_residual,
        partial_transpose_hermiticity_residual=partial_transpose_hermiticity_residual,
        message=if valid
            "candidate passed PSD, trace, PPT, and rank-bound checks"
        else
            "candidate failed at least one PSD, trace, PPT, or rank-bound check"
        end,
    )
end

"""
    random_ppt_state(rng, dims; ranks=nothing, construction=:auto, ...)

Construct a random bipartite PPT state using only the explicit `rng`.

`:shifted_induced` implements the pinned full-rank idea with an explicit
strict-PPT margin. `:separable_mixture` constructs at most
`min(ranks...)` random product projectors and therefore guarantees both rank
bounds without an unbounded nonconvex loop. `:auto` selects the first method
for full ranks and the second otherwise. The documented upstream contract does
not promise a named distribution or entanglement.

No candidate is returned as `state` until Hermiticity, PSD, trace one, partial
transpose PSD, and both requested rank bounds have been checked. Resource
limits are tested before consuming `rng`.
"""
function random_ppt_state(
    rng::AbstractRNG,
    dims;
    ranks=nothing,
    construction::Symbol=:auto,
    real::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
    atol=nothing,
    rtol=nothing,
    max_iterations=1,
    max_dense_entries=1_000_000,
    max_work=100_000_000,
)
    R = _random_ppt_type(T)
    dimensions = _random_ppt_dimensions(dims)
    dimension_big = BigInt(dimensions[1]) * dimensions[2]
    dimension = _symext_checked_int(dimension_big, "prod(dims)")
    requested_ranks = _random_ppt_ranks(ranks, dimension)
    selected_construction = _random_ppt_construction(
        construction, requested_ranks, dimension
    )
    iterations = _symext_nonnegative_integer(max_iterations, "max_iterations")
    entry_limit = _symext_limit(max_dense_entries, "max_dense_entries")
    work_limit = _symext_limit(max_work, "max_work")
    scale = one(R)
    tolerance = _symext_tolerance(R, scale; atol, rtol)
    component_count = min(requested_ranks...)
    entries_needed = if selected_construction === :shifted_induced
        5dimension_big^2
    else
        3dimension_big^2 +
        dimension_big * component_count +
        BigInt(dimensions[1] + dimensions[2]) * component_count
    end
    work_needed = if selected_construction === :shifted_induced
        6dimension_big^3 + 4dimension_big^2
    else
        2dimension_big^3 +
        BigInt(component_count) *
        (BigInt(dimensions[1])^2 + BigInt(dimensions[2])^2 + dimension_big^2)
    end
    if entry_limit !== nothing && entries_needed > entry_limit
        return _random_ppt_failure(
            RandomPPTResourceLimit,
            dimensions,
            requested_ranks,
            selected_construction,
            real,
            iterations,
            tolerance,
            BigInt(0),
            work_limit,
            entry_limit,
            "construction needs $entries_needed dense entries, exceeding max_dense_entries=$entry_limit",
        )
    elseif work_limit !== nothing && work_needed > work_limit
        return _random_ppt_failure(
            RandomPPTResourceLimit,
            dimensions,
            requested_ranks,
            selected_construction,
            real,
            iterations,
            tolerance,
            BigInt(0),
            work_limit,
            entry_limit,
            "construction needs estimated work $work_needed, exceeding max_work=$work_limit",
        )
    elseif iterations == 0
        return _random_ppt_failure(
            RandomPPTIterationLimit,
            dimensions,
            requested_ranks,
            selected_construction,
            real,
            iterations,
            tolerance,
            BigInt(0),
            work_limit,
            entry_limit,
            "max_iterations=0 prevents the single bounded construction step",
        )
    end
    construction_result = if selected_construction === :shifted_induced
        _random_ppt_shifted_induced(rng, dimensions, R, real, tolerance)
    else
        _random_ppt_separable_mixture(rng, dimensions, requested_ranks, R, real)
    end
    construction_result.candidate === nothing && return _random_ppt_failure(
        RandomPPTNumericalFailure,
        dimensions,
        requested_ranks,
        selected_construction,
        real,
        iterations,
        tolerance,
        work_needed,
        work_limit,
        entry_limit,
        construction_result.message;
        iterations=1,
    )
    validation = _random_ppt_validate(
        construction_result.candidate, dimensions, requested_ranks, tolerance
    )
    history = R[max(
        if validation.minimum_eigenvalue === nothing
            R(Inf)
        else
            max(zero(R), -validation.minimum_eigenvalue)
        end,
        if validation.minimum_partial_transpose_eigenvalue === nothing
            R(Inf)
        else
            max(zero(R), -validation.minimum_partial_transpose_eigenvalue)
        end,
        validation.trace_residual,
    ),]
    if !validation.valid
        return RandomPPTStateResult(
            RandomPPTVerificationFailed,
            nothing,
            construction_result.candidate,
            dimensions,
            requested_ranks,
            validation.ranks,
            selected_construction,
            real,
            1,
            iterations,
            history,
            validation.minimum_eigenvalue,
            validation.minimum_partial_transpose_eigenvalue,
            validation.trace_residual,
            validation.hermiticity_residual,
            validation.partial_transpose_hermiticity_residual,
            tolerance,
            work_needed,
            work_limit,
            entry_limit,
            construction_result.normalized,
            false,
            validation.message,
        )
    end
    return RandomPPTStateResult(
        RandomPPTConstructed,
        copy(construction_result.candidate),
        construction_result.candidate,
        dimensions,
        requested_ranks,
        validation.ranks,
        selected_construction,
        real,
        1,
        iterations,
        history,
        validation.minimum_eigenvalue,
        validation.minimum_partial_transpose_eigenvalue,
        validation.trace_residual,
        validation.hermiticity_residual,
        validation.partial_transpose_hermiticity_residual,
        tolerance,
        work_needed,
        work_limit,
        entry_limit,
        construction_result.normalized,
        true,
        construction_result.message * "; " * validation.message,
    )
end
