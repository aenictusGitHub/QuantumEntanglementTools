# Source-informed, independently structured implementation of the primal SDP
# in QETLAB PolynomialSOS.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
#
# The wall-clock-derived random sampling loop is replaced by an exact caller-
# supplied sample count and explicit RNG. Solver data uses the package-owned
# affine SDP interface.

"""
    PolynomialSOSResult

Status-rich output from [`polynomial_sos_bounds`](@ref).

`outer_bound` is the numerical SOS-relaxation bound when a compatible backend
returns usable primal/dual evidence. `inner_bound` is an attained sampled
value. A numerical optimizer result is never relabeled as an exact theorem
certificate; inspect `optimization_result`, residuals, and the returned moment
matrix.
"""
struct PolynomialSOSResult{T,V,M,P,O,U}
    sense::Symbol
    hierarchy_level::Int
    status::OptimizationStatus
    outer_bound::Union{Nothing,T}
    inner_bound::Union{Nothing,T}
    outer_kind::Symbol
    inner_kind::Symbol
    certified_outer::Bool
    best_point::Union{Nothing,V}
    moment_matrix::Union{Nothing,M}
    samples_requested::Int
    samples_evaluated::Int
    target::U
    target_status::Symbol
    problem::P
    optimization_result::O
    message::String
end

function Base.show(io::IO, result::PolynomialSOSResult)
    return print(
        io,
        "PolynomialSOSResult(sense=",
        result.sense,
        ", level=",
        result.hierarchy_level,
        ", status=",
        result.status,
        ", outer_bound=",
        result.outer_bound,
        ", inner_bound=",
        result.inner_bound,
        ")",
    )
end

function _polynomial_sos_limit(value, name::AbstractString; allow_zero::Bool=false)
    return _polynomial_limit(value, name; allow_zero)
end

function _polynomial_sos_preflight(
    polynomial::HomogeneousPolynomial,
    level::Int,
    limits::OptimizationLimits;
    max_terms,
    max_full_dimension,
    max_projection_entries,
    max_work,
)
    half_degree = div(polynomial.degree, 2)
    copies = _polynomial_checked_sum(half_degree, level, "half_degree + level")
    copies > 0 || throw(ArgumentError("a constant polynomial has no nontrivial SOS model"))
    n = polynomial.variables
    full_dimension = BigInt(n)^copies
    full_limit = _polynomial_sos_limit(max_full_dimension, "max_full_dimension")
    full_limit !== nothing &&
        full_dimension > full_limit &&
        throw(
            ArgumentError(
                "SOS symmetric-tensor dimension $full_dimension exceeds " *
                "max_full_dimension=$full_limit",
            ),
        )
    projection_limit = _polynomial_sos_limit(
        max_projection_entries, "max_projection_entries"
    )
    projection_limit !== nothing &&
        full_dimension > projection_limit &&
        throw(
            ArgumentError(
                "the symmetric basis stores $full_dimension entries, exceeding " *
                "max_projection_entries=$projection_limit",
            ),
        )
    matrix_dimension = _polynomial_monomial_count(n, copies; max_terms=max_terms)
    variable_count = BigInt(matrix_dimension)^2
    variable_count <= limits.max_variables || throw(
        ArgumentError(
            "SOS moment matrix needs $variable_count variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    matrix_dimension <= limits.max_psd_dimension || throw(
        ArgumentError(
            "SOS moment matrix dimension $matrix_dimension exceeds " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    equality_count = full_dimension^2 + 1
    equality_count <= limits.max_equalities || throw(
        ArgumentError(
            "SOS invariance needs $equality_count scalar equalities, exceeding " *
            "max_equalities=$(limits.max_equalities)",
        ),
    )
    work = full_dimension^2 * variable_count
    work_limit = _polynomial_sos_limit(max_work, "max_work")
    work_limit !== nothing &&
        work > work_limit &&
        throw(
            ArgumentError(
                "SOS affine-invariance construction work estimate $work exceeds " *
                "max_work=$work_limit",
            ),
        )
    return (
        copies=copies,
        full_dimension=Int(full_dimension),
        matrix_dimension=matrix_dimension,
        variable_count=Int(variable_count),
        estimated_work=work,
    )
end

function _polynomial_sos_objective(
    polynomial_matrix::AbstractMatrix{<:Real}, moment::HermitianAffineMatrix{T}
) where {T<:Real}
    variables = Int[]
    values = T[]
    for term in moment.terms
        value = real(tr(polynomial_matrix * term.coefficient))
        iszero(value) && continue
        push!(variables, term.variable)
        push!(values, convert(T, value))
    end
    constant = convert(T, real(tr(polynomial_matrix * moment.constant)))
    return AffineScalar(constant, variables, values, moment.variable_count)
end

function _polynomial_sos_apply_basis(
    moment::HermitianAffineMatrix, basis::AbstractMatrix{<:Number}, name::Symbol
)
    basis_adjoint = adjoint(basis)
    return HermitianAffineMatrix(
        name,
        sparse(basis * moment.constant * basis_adjoint),
        [term.variable for term in moment.terms],
        [sparse(basis * term.coefficient * basis_adjoint) for term in moment.terms],
        moment.variable_count,
    )
end

function _polynomial_sos_affine_difference(
    left::HermitianAffineMatrix{T}, right::HermitianAffineMatrix{T}, name::Symbol
) where {T<:Real}
    left.dimension == right.dimension ||
        throw(DimensionMismatch("affine matrix dimensions differ"))
    left.variable_count == right.variable_count ||
        throw(DimensionMismatch("affine coordinate counts differ"))
    coefficients = Dict{Int,SparseMatrixCSC{Complex{T},Int}}()
    for term in left.terms
        coefficients[term.variable] = copy(term.coefficient)
    end
    for term in right.terms
        if haskey(coefficients, term.variable)
            coefficients[term.variable] = coefficients[term.variable] - term.coefficient
        else
            coefficients[term.variable] = -term.coefficient
        end
    end
    variables = sort!([
        variable for (variable, coefficient) in coefficients if nnz(coefficient) > 0
    ],)
    return HermitianAffineMatrix(
        name,
        left.constant - right.constant,
        variables,
        [coefficients[variable] for variable in variables],
        left.variable_count,
    )
end

"""
    polynomial_sos_problem(polynomial; level=0, sense=:max, kwargs...)

Construct the package-owned primal SDP for the QETLAB sum-of-squares
hierarchy. For degree `2d` and hierarchy level `k`, the moment variable acts
on the normalized degree-`d+k` symmetric monomial basis. Its lifted
`n^(d+k)`-dimensional operator is positive semidefinite, trace one, and
invariant under partial transpose of the first tensor factor.

No solver is selected. Exponential tensor dimensions, affine work, model
entries, variables, equalities, and PSD blocks are checked before allocation.
"""
function polynomial_sos_problem(
    polynomial::HomogeneousPolynomial;
    level=0,
    sense::Symbol=:max,
    limits::OptimizationLimits=OptimizationLimits(),
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_full_dimension=512,
    max_projection_entries=100_000,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    _validate_polynomial(polynomial; max_terms)
    eltype(polynomial) <: Real || throw(
        ArgumentError(
            "polynomial_sos_problem requires real coefficients; got " *
            "$(eltype(polynomial))",
        ),
    )
    iseven(polynomial.degree) || throw(
        ArgumentError(
            "polynomial_sos_problem requires an even degree; got " * "$(polynomial.degree)",
        ),
    )
    checked_sense = _polynomial_sense(sense)
    checked_level = _polynomial_nonnegative_integer(level, "level")
    preflight = _polynomial_sos_preflight(
        polynomial,
        checked_level,
        limits;
        max_terms,
        max_full_dimension,
        max_projection_entries,
        max_work,
    )
    solver_type = _polynomial_solver_type(eltype(polynomial))
    converted = HomogeneousPolynomial(
        solver_type.(polynomial.coefficients),
        polynomial.variables,
        polynomial.degree;
        max_terms,
    )
    polynomial_matrix = polynomial_as_matrix(
        converted;
        level=checked_level,
        sparse_output=true,
        max_terms,
        max_degree,
        max_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_work=nothing,
    )
    size(polynomial_matrix, 1) == preflight.matrix_dimension ||
        error("internal polynomial matrix dimension is inconsistent")

    moment = hermitian_variable(
        :polynomial_sos_moment, preflight.matrix_dimension; coefficient_type=solver_type
    )
    basis = symmetric_subspace_basis(
        polynomial.variables, preflight.copies; T=solver_type, sparse_output=true
    )
    size(basis) == (preflight.full_dimension, preflight.matrix_dimension) ||
        error("internal symmetric-basis dimension is inconsistent")
    lifted = _polynomial_sos_apply_basis(moment, basis, :polynomial_sos_lifted_moment)
    transposed = partial_transpose_affine(
        lifted,
        ntuple(_ -> polynomial.variables, preflight.copies);
        systems=(1,),
        name=:polynomial_sos_partial_transpose,
    )
    invariance = _polynomial_sos_affine_difference(
        transposed, lifted, :polynomial_sos_invariance
    )
    equalities = hermitian_equalities(invariance; name_prefix=:polynomial_sos_invariance)
    trace_data = trace_affine(moment)
    trace_indices, trace_values = findnz(trace_data.coefficients)
    push!(
        equalities,
        AffineEquality(
            AffineScalar(
                trace_data.constant - one(solver_type),
                trace_indices,
                trace_values,
                moment.variable_count,
            ),
            :polynomial_sos_unit_trace,
        ),
    )
    objective = _polynomial_sos_objective(polynomial_matrix, moment)
    return SemidefiniteProgram(
        :polynomial_sos,
        checked_sense === :max ? :maximize : :minimize,
        moment.variable_count,
        objective;
        equalities,
        psd_constraints=[moment],
        primal_views=[moment],
        limits,
        metadata=(
            upstream_function="PolynomialSOS",
            formulation=:symmetric_moment_partial_transpose_relaxation,
            hierarchy_level=checked_level,
            polynomial_variables=polynomial.variables,
            polynomial_degree=polynomial.degree,
            tensor_copies=preflight.copies,
            full_tensor_dimension=preflight.full_dimension,
            moment_dimension=preflight.matrix_dimension,
            estimated_work=preflight.estimated_work,
            outer_semantics=checked_sense === :max ? :upper_bound : :lower_bound,
        ),
    )
end

function _polynomial_sos_result_bound(result::OptimizationResult, sense::Symbol)
    result.status in (OptimizationOptimal, OptimizationFeasible) || return nothing
    if result.objective_bound !== nothing
        return result.objective_bound
    end
    result.status === OptimizationOptimal || return nothing
    return result.objective_value
end

function _polynomial_sos_sample(
    rng::AbstractRNG,
    polynomial::HomogeneousPolynomial,
    level::Int,
    sense::Symbol,
    inner_samples::Int;
    max_terms,
    max_degree,
    max_dimension,
    max_exponent_entries,
    max_dense_entries,
    max_nonzeros,
    max_samples,
    max_work,
)
    inner_samples == 0 &&
        return (bound=nothing, point=nothing, evaluated=0, kind=:not_requested)
    sampled = polynomial_bounds(
        rng,
        polynomial;
        level,
        sense,
        target=nothing,
        inner_samples,
        allow_densify=true,
        max_terms,
        max_degree,
        max_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_samples,
        max_work,
    )
    return (
        bound=sampled.inner_bound,
        point=sampled.best_point,
        evaluated=sampled.samples_evaluated,
        kind=sampled.inner_kind,
    )
end

"""
    polynomial_sos_bounds(rng, polynomial; backend=NoOptimizationBackend(),
                          level=0, sense=:max, inner_samples=0, kwargs...)

Solve or report the SOS hierarchy outer bound and optionally evaluate an
exact caller-selected number of random unit-sphere points.

The optimizer is explicit. Missing, limited, inaccurate, inconsistent, and
failed outcomes retain their [`OptimizationStatus`](@ref). A sampled point is
an attained inner value, while a floating conic solve is numerical relaxation
evidence and sets `certified_outer=false`.
"""
function polynomial_sos_bounds(
    rng::AbstractRNG,
    polynomial::HomogeneousPolynomial;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    level=0,
    sense::Symbol=:max,
    target=nothing,
    inner_samples=0,
    limits::OptimizationLimits=OptimizationLimits(),
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    max_degree=_POLYNOMIAL_DEFAULT_MAX_DEGREE,
    max_dimension=_POLYNOMIAL_DEFAULT_MAX_DIMENSION,
    max_exponent_entries=_POLYNOMIAL_DEFAULT_MAX_EXPONENT_ENTRIES,
    max_dense_entries=_POLYNOMIAL_DEFAULT_MAX_DENSE_ENTRIES,
    max_nonzeros=_POLYNOMIAL_DEFAULT_MAX_NONZEROS,
    max_samples=_POLYNOMIAL_DEFAULT_MAX_SAMPLES,
    max_full_dimension=512,
    max_projection_entries=100_000,
    max_work=_POLYNOMIAL_DEFAULT_MAX_WORK,
)
    checked_sense = _polynomial_sense(sense)
    checked_level = _polynomial_nonnegative_integer(level, "level")
    checked_samples = _polynomial_nonnegative_integer(inner_samples, "inner_samples")
    sample_limit = _polynomial_limit(max_samples, "max_samples"; allow_zero=true)
    sample_limit !== nothing &&
        checked_samples > sample_limit &&
        throw(
            ArgumentError(
                "inner_samples=$checked_samples exceeds max_samples=$sample_limit"
            ),
        )
    solver_type = _polynomial_solver_type(eltype(polynomial))
    checked_target = _polynomial_target(target, solver_type)

    if polynomial.degree == 0
        _validate_polynomial(polynomial; max_terms)
        value = convert(solver_type, only(polynomial.coefficients))
        status = OptimizationOptimal
        target_status = _polynomial_target_status(checked_sense, value, checked_target)
        return PolynomialSOSResult{
            solver_type,
            Vector{solver_type},
            Matrix{solver_type},
            Nothing,
            Nothing,
            typeof(checked_target),
        }(
            checked_sense,
            checked_level,
            status,
            value,
            value,
            :exact_constant,
            :exact_constant,
            true,
            zeros(solver_type, polynomial.variables),
            nothing,
            checked_samples,
            0,
            checked_target,
            target_status,
            nothing,
            nothing,
            "a degree-zero homogeneous polynomial is constant",
        )
    end

    problem = polynomial_sos_problem(
        polynomial;
        level=checked_level,
        sense=checked_sense,
        limits,
        max_terms,
        max_degree,
        max_dimension,
        max_exponent_entries,
        max_dense_entries,
        max_nonzeros,
        max_full_dimension,
        max_projection_entries,
        max_work,
    )
    optimization = solve_optimization(problem, backend)
    outer = _polynomial_sos_result_bound(optimization, checked_sense)
    target_status = if outer === nothing
        checked_target === nothing ? :not_requested : :unknown
    else
        _polynomial_target_status(checked_sense, outer, checked_target)
    end
    target_met = target_status in (:outer_proves_at_most, :outer_proves_at_least)
    sample = if target_met
        (bound=nothing, point=nothing, evaluated=0, kind=:skipped_by_outer_target)
    else
        _polynomial_sos_sample(
            rng,
            polynomial,
            checked_level,
            checked_sense,
            checked_samples;
            max_terms,
            max_degree,
            max_dimension,
            max_exponent_entries,
            max_dense_entries,
            max_nonzeros,
            max_samples,
            max_work,
        )
    end
    moment = if optimization.primal === nothing
        nothing
    else
        Matrix(optimization.primal.views.polynomial_sos_moment)
    end
    message = if optimization.status === OptimizationBackendUnavailable
        "no optimization backend was supplied; only requested sampled inner evidence is available"
    elseif optimization.status in (OptimizationOptimal, OptimizationFeasible)
        "the outer value is a numerical SOS relaxation bound; the inner value, when present, is attained by the returned point"
    else
        "the optimizer did not produce a usable outer bound; its structured status and diagnostics are retained"
    end
    return PolynomialSOSResult{
        solver_type,
        Vector{solver_type},
        Matrix{Complex{solver_type}},
        typeof(problem),
        typeof(optimization),
        typeof(checked_target),
    }(
        checked_sense,
        checked_level,
        optimization.status,
        outer,
        sample.bound,
        :sos_relaxation_numerical,
        sample.kind,
        false,
        sample.point,
        moment,
        checked_samples,
        sample.evaluated,
        checked_target,
        target_status,
        problem,
        optimization,
        message,
    )
end

function polynomial_sos_bounds(
    rng::AbstractRNG,
    coefficients::AbstractVector,
    variables,
    half_degree,
    level;
    max_terms=_POLYNOMIAL_DEFAULT_MAX_TERMS,
    kwargs...,
)
    checked_half_degree = _polynomial_nonnegative_integer(half_degree, "half_degree")
    degree = Base.checked_mul(2, checked_half_degree)
    polynomial = HomogeneousPolynomial(coefficients, variables, degree; max_terms)
    return polynomial_sos_bounds(rng, polynomial; level, max_terms, kwargs...)
end
