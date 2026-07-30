# Source-informed independent Julia implementation based on the executable
# contract in QETLAB Distinguishability.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The solver-independent primal/dual formulation follows Section 3.1 of
# John Watrous, The Theory of Quantum Information:
# https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf

"""
    StateDiscriminationStatus

Certificate-aware status returned by [`state_distinguishability`](@ref).
Analytic and exactly orthogonal branches are separate from numerical SDP
outcomes. A missing backend, resource limit, or failed residual check is always
inconclusive.
"""
@enum StateDiscriminationStatus::UInt8 begin
    StateDiscriminationTrivialOptimal
    StateDiscriminationHelstromOptimal
    StateDiscriminationOrthogonalOptimal
    StateDiscriminationSolverOptimal
    StateDiscriminationSolverFeasible
    StateDiscriminationNumericalBoundary
    StateDiscriminationResourceLimit
    StateDiscriminationBackendUnavailable
    StateDiscriminationBackendFailure
    StateDiscriminationInvalidCertificate
end

"""
    StateDiscriminationProblem

Package-owned minimum-error state-discrimination SDP. `program` contains no
JuMP or optimizer objects. The primal measurement effects are the named affine
views in `measurement_views`.
"""
struct StateDiscriminationProblem{T<:AbstractFloat,S,P,V,M}
    states::S
    priors::P
    input_kind::Symbol
    dimension::Int
    state_count::Int
    program::V
    measurement_views::M
    tolerance::T
end

"""
    StateDiscriminationResult

Status-rich minimum-error discrimination result.

`lower_bound` is the success probability of the returned POVM whenever a
validated POVM is available. `upper_bound` is analytic for the exact branches
and otherwise comes only from an available dual/optimizer bound. The scalar
`success_probability` is populated only when the two bounds agree within the
recorded tolerance. `measurement` is a tuple of owned effects; it is never
returned unless positivity and completeness have been checked independently.

`certified=true` is reserved for the one-state/deterministic-prior theorem, the
Helstrom theorem, or an exactly orthogonal support decomposition. A numerical
solver optimum retains its primal/dual diagnostics but is not promoted to an
exact mathematical certificate.
"""
struct StateDiscriminationResult{T<:AbstractFloat,P,Q,M,D,O,R}
    status::StateDiscriminationStatus
    success_probability::Union{Nothing,T}
    lower_bound::Union{Nothing,T}
    upper_bound::Union{Nothing,T}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    dimension::Int
    state_count::Int
    priors::P
    input_kind::Symbol
    problem::Q
    measurement::M
    dual_operator::D
    optimization_result::O
    residuals::R
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

function Base.show(io::IO, result::StateDiscriminationResult)
    return print(
        io,
        "StateDiscriminationResult(status=",
        result.status,
        ", success_probability=",
        result.success_probability,
        ", lower_bound=",
        result.lower_bound,
        ", upper_bound=",
        result.upper_bound,
        ", state_count=",
        result.state_count,
        ")",
    )
end

function _discrimination_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _discrimination_real_type(::Type{T}) where {T<:Number}
    R = typeof(real(zero(T)))
    R in (Float32, Float64) || throw(
        ArgumentError(
            "state discrimination currently requires Float32, Float64, " *
            "ComplexF32, or ComplexF64 inputs; got $T. No implicit " *
            "precision-changing conversion is performed.",
        ),
    )
    return R
end

function _discrimination_common_real_type(states, priors)
    types = Type[]
    for state in states
        push!(types, _discrimination_real_type(eltype(state)))
    end
    if priors !== nothing && eltype(priors) in (Float32, Float64)
        push!(types, eltype(priors))
    end
    T = promote_type(types...)
    T in (Float32, Float64) ||
        throw(ArgumentError("state and prior types must promote to Float32 or Float64"))
    return T
end

function _discrimination_validate_limits(
    dimension::Int, state_count::Int; max_dimension, max_states, max_variables
)
    dimension <= _discrimination_positive_integer(max_dimension, "max_dimension") || throw(
        ArgumentError("state dimension $dimension exceeds max_dimension=$max_dimension")
    )
    state_count <= _discrimination_positive_integer(max_states, "max_states") ||
        throw(ArgumentError("state count $state_count exceeds max_states=$max_states"))
    variable_count = BigInt(state_count) * BigInt(dimension)^2
    variable_count <= _discrimination_positive_integer(max_variables, "max_variables") ||
        throw(
            ArgumentError(
                "the POVM model needs $variable_count real variables, exceeding " *
                "max_variables=$max_variables",
            ),
        )
    variable_count <= typemax(Int) ||
        throw(ArgumentError("the POVM variable count is too large for Int"))
    return Int(variable_count)
end

function _discrimination_priors(priors, state_count::Int, ::Type{T}; atol, rtol) where {T}
    absolute, relative = _tierd_tolerances(T, atol, rtol)
    values = if priors === nothing
        fill(one(T) / convert(T, state_count), state_count)
    else
        priors isa AbstractVector ||
            throw(ArgumentError("priors must be nothing or an AbstractVector"))
        Base.require_one_based_indexing(priors)
        length(priors) == state_count || throw(
            DimensionMismatch("expected $state_count priors, got $(length(priors))")
        )
        Vector{T}(undef, state_count)
    end
    if priors !== nothing
        for index in 1:state_count
            value = priors[index]
            value isa Real && !(value isa Bool) && isfinite(value) ||
                throw(ArgumentError("priors must contain finite real values"))
            value >= zero(value) || throw(DomainError(value, "priors must be nonnegative"))
            values[index] = convert(T, value)
        end
    end
    total = sum(values)
    tolerance = absolute + relative * max(one(T), abs(total))
    abs(total - one(T)) <= tolerance || throw(
        ArgumentError(
            "priors must sum to one within atol=$absolute and rtol=$relative; " *
            "the sum is $total. Priors are never normalized implicitly.",
        ),
    )
    return values, tolerance
end

function _discrimination_prepare_density_states(
    states;
    priors,
    atol,
    rtol,
    allow_densify::Bool,
    max_dimension,
    max_states,
    max_variables,
)
    states isa Union{Tuple,AbstractVector} ||
        throw(ArgumentError("density states must be supplied as a tuple or vector"))
    states isa AbstractArray && Base.require_one_based_indexing(states)
    isempty(states) && throw(ArgumentError("at least one state is required"))
    all(state -> state isa AbstractMatrix{<:Number}, states) ||
        throw(ArgumentError("every density-state entry must be a numeric matrix"))
    dimension = size(first(states), 1)
    size(first(states), 2) == dimension ||
        throw(DimensionMismatch("density states must be square"))
    dimension > 0 || throw(ArgumentError("states must have positive dimension"))
    state_count = length(states)
    _discrimination_validate_limits(
        dimension,
        state_count;
        max_dimension=max_dimension,
        max_states=max_states,
        max_variables=max_variables,
    )
    T = _discrimination_common_real_type(states, priors)
    prior_values, prior_tolerance = _discrimination_priors(
        priors, state_count, T; atol=atol, rtol=rtol
    )
    prepared = Matrix{Complex{T}}[]
    tolerances = T[prior_tolerance]
    boundary = false
    normalization_residual = zero(T)
    for (index, state) in enumerate(states)
        size(state) == (dimension, dimension) || throw(
            DimensionMismatch(
                "state $index has size $(size(state)); expected " *
                "($dimension, $dimension)",
            ),
        )
        ishermitian(state) || throw(
            ArgumentError(
                "state $index must be exactly Hermitian; no symmetrization is applied"
            ),
        )
        analysis = _tierd_density_analysis(
            state;
            atol=atol,
            rtol=rtol,
            allow_densify=allow_densify,
            boundary_policy=:reject,
            operation="state_distinguishability",
        )
        # _tierd_density_analysis checks a residual before constructing a
        # Hermitian work matrix. Requiring exact Hermiticity here ensures that
        # this API never replaces a caller's state by its Hermitian part.
        push!(prepared, Matrix{Complex{T}}(state))
        push!(tolerances, convert(T, analysis.tolerance))
        normalization_residual = max(
            normalization_residual, convert(T, analysis.normalization_residual)
        )
        boundary |= analysis.boundary_uncertain
    end
    return (
        states=prepared,
        priors=prior_values,
        input_kind=:density_matrices,
        dimension=dimension,
        state_count=state_count,
        tolerance=maximum(tolerances),
        input_boundary=boundary,
        normalization_residual=normalization_residual,
    )
end

function _discrimination_prepare_pure_columns(
    columns::AbstractMatrix{<:Number};
    priors,
    atol,
    rtol,
    allow_densify::Bool,
    max_dimension,
    max_states,
    max_variables,
)
    Base.require_one_based_indexing(columns)
    dimension, state_count = size(columns)
    dimension > 0 && state_count > 0 ||
        throw(ArgumentError("the pure-state matrix must have at least one row and column"))
    issparse(columns) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "pure-state column validation and POVM construction densify sparse input; " *
                "pass allow_densify=true",
            ),
        )
    _discrimination_validate_limits(
        dimension,
        state_count;
        max_dimension=max_dimension,
        max_states=max_states,
        max_variables=max_variables,
    )
    T = _discrimination_common_real_type((columns,), priors)
    prior_values, prior_tolerance = _discrimination_priors(
        priors, state_count, T; atol=atol, rtol=rtol
    )
    prepared = Matrix{Complex{T}}[]
    tolerance = prior_tolerance
    normalization_residual = zero(T)
    for index in 1:state_count
        vector = Vector{Complex{T}}(@view columns[:, index])
        validation = _tierd_validate_pure_state(
            vector; atol=atol, rtol=rtol, operation="state_distinguishability"
        )
        tolerance = max(tolerance, convert(T, validation.atol + validation.rtol))
        normalization_residual = max(
            normalization_residual, abs(convert(T, validation.norm_squared) - one(T))
        )
        push!(prepared, vector * adjoint(vector))
    end
    return (
        states=prepared,
        priors=prior_values,
        input_kind=:pure_columns,
        dimension=dimension,
        state_count=state_count,
        tolerance=tolerance,
        input_boundary=(!iszero(normalization_residual)),
        normalization_residual=normalization_residual,
    )
end

function _discrimination_prepare(
    states;
    priors=nothing,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dimension=64,
    max_states=64,
    max_variables=100_000,
)
    if states isa AbstractMatrix{<:Number}
        return _discrimination_prepare_pure_columns(
            states;
            priors=priors,
            atol=atol,
            rtol=rtol,
            allow_densify=allow_densify,
            max_dimension=max_dimension,
            max_states=max_states,
            max_variables=max_variables,
        )
    end
    return _discrimination_prepare_density_states(
        states;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dimension=max_dimension,
        max_states=max_states,
        max_variables=max_variables,
    )
end

function _discrimination_affine_sum(
    matrices::AbstractVector{<:HermitianAffineMatrix{T}}, name::Symbol
) where {T<:Real}
    isempty(matrices) && throw(ArgumentError("at least one affine matrix is required"))
    dimension = first(matrices).dimension
    variable_count = first(matrices).variable_count
    all(
        matrix -> matrix.dimension == dimension && matrix.variable_count == variable_count,
        matrices,
    ) || throw(DimensionMismatch("affine matrices have incompatible dimensions"))
    constant = spzeros(Complex{T}, dimension, dimension)
    variables = Int[]
    coefficients = SparseMatrixCSC{Complex{T},Int}[]
    for matrix in matrices
        constant += matrix.constant
        for term in matrix.terms
            push!(variables, term.variable)
            push!(coefficients, term.coefficient)
        end
    end
    return HermitianAffineMatrix(name, constant, variables, coefficients, variable_count)
end

function _discrimination_objective(prepared, views, variable_count::Int)
    T = eltype(prepared.priors)
    variables = Int[]
    coefficients = T[]
    for (prior, state, view) in zip(prepared.priors, prepared.states, views)
        for term in view.terms
            value = prior * real(dot(state, term.coefficient))
            iszero(value) && continue
            push!(variables, term.variable)
            push!(coefficients, convert(T, value))
        end
    end
    return AffineScalar(zero(T), variables, coefficients, variable_count)
end

function _discrimination_known_povm_coordinates(prepared, variable_count::Int)
    T = eltype(prepared.priors)
    coordinates = zeros(T, variable_count)
    weight = one(T) / convert(T, prepared.state_count)
    block = prepared.dimension^2
    for state_index in 1:prepared.state_count
        first_variable = (state_index - 1) * block + 1
        coordinates[first_variable:(first_variable + prepared.dimension - 1)] .= weight
    end
    return coordinates
end

function _discrimination_build_problem(
    prepared; limits::OptimizationLimits=OptimizationLimits()
)
    dimension = prepared.dimension
    state_count = prepared.state_count
    variable_count = state_count * dimension^2
    T = eltype(prepared.priors)
    views = HermitianAffineMatrix{T}[]
    for index in 1:state_count
        first_variable = (index - 1) * dimension^2 + 1
        push!(
            views,
            hermitian_variable(
                Symbol(:measurement_, index),
                dimension;
                first_variable=first_variable,
                variable_count=variable_count,
                coefficient_type=T,
            ),
        )
    end
    sum_view = _discrimination_affine_sum(views, :measurement_sum)
    equalities = hermitian_equalities(
        sum_view;
        target=Matrix{Complex{T}}(I, dimension, dimension),
        name_prefix=:measurement_complete,
    )
    objective = _discrimination_objective(prepared, views, variable_count)
    feasible = _discrimination_known_povm_coordinates(prepared, variable_count)
    program = SemidefiniteProgram(
        :minimum_error_state_discrimination,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=views,
        primal_views=views,
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:minimum_error_povm_primal,
            state_count=state_count,
            dimension=dimension,
            input_kind=prepared.input_kind,
            source=:watrous_state_discrimination_sdp,
        ),
    )
    return StateDiscriminationProblem(
        prepared.states,
        prepared.priors,
        prepared.input_kind,
        dimension,
        state_count,
        program,
        Tuple(views),
        prepared.tolerance,
    )
end

"""
    state_discrimination_problem(states; priors=nothing, ...)

Build the solver-neutral minimum-error POVM SDP

```math
\\max_{M_j \\succeq 0}\\;\\sum_j p_j\\operatorname{tr}(M_j\\rho_j),
\\qquad \\sum_j M_j=I.
```

`states` is either a tuple/vector of normalized density matrices or a matrix
whose columns are normalized pure states. Inputs and priors are never
normalized, symmetrized, clipped, or otherwise repaired. Sparse spectral
validation requires `allow_densify=true`.
"""
function state_discrimination_problem(
    states;
    priors=nothing,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dimension=64,
    max_states=64,
    max_variables=100_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _discrimination_prepare(
        states;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dimension=max_dimension,
        max_states=max_states,
        max_variables=max_variables,
    )
    return _discrimination_build_problem(prepared; limits=limits)
end

function _discrimination_measurement_residuals(measurement, prepared, tolerance)
    T = eltype(prepared.priors)
    dimension = prepared.dimension
    total = zeros(Complex{T}, dimension, dimension)
    minimum_eigenvalue = convert(T, Inf)
    hermiticity_residual = zero(T)
    for effect in measurement
        size(effect) == (dimension, dimension) || return (
            valid=false,
            success_probability=nothing,
            completeness_residual=convert(T, Inf),
            positivity_violation=convert(T, Inf),
            minimum_eigenvalue=(-convert(T, Inf)),
            hermiticity_residual=convert(T, Inf),
        )
        hermiticity_residual = max(
            hermiticity_residual, maximum(abs, effect - adjoint(effect); init=zero(T))
        )
        values = eigvals(Hermitian(Matrix(effect)))
        minimum_eigenvalue = min(minimum_eigenvalue, minimum(real, values))
        total .+= effect
    end
    completeness_residual = maximum(
        abs, total - Matrix{Complex{T}}(I, dimension, dimension); init=zero(T)
    )
    positivity_violation = max(zero(T), -minimum_eigenvalue)
    success = zero(T)
    for (prior, state, effect) in zip(prepared.priors, prepared.states, measurement)
        success += prior * real(dot(state, effect))
    end
    valid =
        hermiticity_residual <= tolerance &&
        completeness_residual <= tolerance &&
        positivity_violation <= tolerance &&
        success >= -tolerance &&
        success <= one(T) + tolerance
    return (
        valid=valid,
        success_probability=success,
        completeness_residual=completeness_residual,
        positivity_violation=positivity_violation,
        minimum_eigenvalue=minimum_eigenvalue,
        hermiticity_residual=hermiticity_residual,
    )
end

function _discrimination_result(
    status,
    prepared;
    success_probability=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    certified=false,
    certificate_kind=nothing,
    problem=nothing,
    measurement=nothing,
    dual_operator=nothing,
    optimization_result=nothing,
    residuals=NamedTuple(),
    warnings=(),
    message,
)
    T = eltype(prepared.priors)
    return StateDiscriminationResult(
        status,
        success_probability === nothing ? nothing : convert(T, success_probability),
        lower_bound === nothing ? nothing : convert(T, lower_bound),
        upper_bound === nothing ? nothing : convert(T, upper_bound),
        certified,
        certificate_kind,
        prepared.dimension,
        prepared.state_count,
        copy(prepared.priors),
        prepared.input_kind,
        problem,
        measurement,
        dual_operator,
        optimization_result,
        residuals,
        convert(T, prepared.tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

function _discrimination_trivial(prepared)
    T = eltype(prepared.priors)
    deterministic = findfirst(==(one(T)), prepared.priors)
    prepared.state_count == 1 || deterministic !== nothing || return nothing
    chosen = prepared.state_count == 1 ? 1 : deterministic
    measurement = Tuple(
        if index == chosen
            Matrix{Complex{T}}(I, prepared.dimension, prepared.dimension)
        else
            zeros(Complex{T}, prepared.dimension, prepared.dimension)
        end for index in 1:prepared.state_count
    )
    residuals = _discrimination_measurement_residuals(
        measurement, prepared, prepared.tolerance
    )
    value = residuals.success_probability
    return _discrimination_result(
        StateDiscriminationTrivialOptimal,
        prepared;
        success_probability=value,
        lower_bound=value,
        upper_bound=value,
        certified=residuals.valid,
        certificate_kind=:deterministic_prior,
        measurement=measurement,
        residuals=residuals,
        message="a single possible state is distinguished with probability one",
    )
end

function _discrimination_helstrom(prepared)
    prepared.state_count == 2 || return nothing
    T = eltype(prepared.priors)
    difference =
        prepared.priors[1] * prepared.states[1] - prepared.priors[2] * prepared.states[2]
    decomposition = eigen(Hermitian(difference))
    values = decomposition.values
    positive_indices = findall(value -> value >= zero(value), values)
    first_effect = if isempty(positive_indices)
        zeros(Complex{T}, prepared.dimension, prepared.dimension)
    else
        vectors = decomposition.vectors[:, positive_indices]
        Matrix{Complex{T}}(vectors * adjoint(vectors))
    end
    second_effect =
        Matrix{Complex{T}}(I, prepared.dimension, prepared.dimension) - first_effect
    measurement = (first_effect, second_effect)
    weighted_trace = sum(
        prior * real(tr(state)) for (prior, state) in zip(prepared.priors, prepared.states)
    )
    value = (weighted_trace + sum(abs, values)) / convert(T, 2)
    residuals = _discrimination_measurement_residuals(
        measurement, prepared, 8 * prepared.tolerance
    )
    objective_residual = if residuals.success_probability === nothing
        convert(T, Inf)
    else
        abs(residuals.success_probability - value)
    end
    spectral_boundary = any(
        eigenvalue -> !iszero(eigenvalue) && abs(eigenvalue) <= 8 * prepared.tolerance,
        values,
    )
    valid = residuals.valid && objective_residual <= 8 * prepared.tolerance
    warnings = String[]
    spectral_boundary && push!(
        warnings,
        "the Helstrom operator has a nonzero eigenvalue inside the numerical boundary",
    )
    return _discrimination_result(
        if spectral_boundary
            StateDiscriminationNumericalBoundary
        else
            StateDiscriminationHelstromOptimal
        end,
        prepared;
        success_probability=value,
        lower_bound=value,
        upper_bound=value,
        certified=valid && !spectral_boundary,
        certificate_kind=valid && !spectral_boundary ? :helstrom_theorem : nothing,
        measurement=valid ? measurement : nothing,
        residuals=merge(
            residuals,
            (
                helstrom_objective_residual=objective_residual,
                helstrom_minimum_absolute_nonzero_eigenvalue=minimum(
                    (abs(value) for value in values if !iszero(value)); init=convert(T, Inf)
                ),
            ),
        ),
        warnings=warnings,
        message=if spectral_boundary
            "the Helstrom value was computed, but its support projector is on a numerical boundary"
        else
            "the two-state optimum and POVM follow from the Helstrom theorem"
        end,
    )
end

function _discrimination_exact_orthogonal(prepared)
    prepared.state_count <= prepared.dimension || return nothing
    for left in 1:(prepared.state_count - 1)
        for right in (left + 1):prepared.state_count
            iszero(prepared.states[left] * prepared.states[right]) || return nothing
        end
    end
    T = eltype(prepared.priors)
    support_projectors = Matrix{Complex{T}}[]
    for state in prepared.states
        decomposition = eigen(Hermitian(state))
        indices = findall(value -> !iszero(value), decomposition.values)
        projector = if isempty(indices)
            zeros(Complex{T}, prepared.dimension, prepared.dimension)
        else
            vectors = decomposition.vectors[:, indices]
            Matrix{Complex{T}}(vectors * adjoint(vectors))
        end
        push!(support_projectors, projector)
    end
    leftover = Matrix{Complex{T}}(I, prepared.dimension, prepared.dimension)
    for index in 2:prepared.state_count
        leftover .-= support_projectors[index]
    end
    measurement = Tuple(
        index == 1 ? leftover : support_projectors[index] for
        index in 1:prepared.state_count
    )
    residuals = _discrimination_measurement_residuals(
        measurement, prepared, 8 * prepared.tolerance
    )
    residuals.valid || return nothing
    value = residuals.success_probability
    return _discrimination_result(
        StateDiscriminationOrthogonalOptimal,
        prepared;
        success_probability=value,
        lower_bound=value,
        upper_bound=value,
        certified=true,
        certificate_kind=:orthogonal_supports,
        measurement=measurement,
        residuals=residuals,
        message="the state supports are exactly orthogonal and the support POVM succeeds with probability one",
    )
end

function _discrimination_hermitian_from_duals(values, dimension::Int)
    length(values) == dimension^2 ||
        throw(DimensionMismatch("wrong number of measurement-completeness duals"))
    T = eltype(values)
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

function _discrimination_dual_candidate(optimization, prepared, tolerance)
    optimization.dual === nothing && return nothing
    count = prepared.dimension^2
    length(optimization.dual.equalities) >= count || return nothing
    raw = _discrimination_hermitian_from_duals(
        view(optimization.dual.equalities, 1:count), prepared.dimension
    )
    best = nothing
    for candidate in (raw, -raw)
        minimum_slack = convert(eltype(prepared.priors), Inf)
        for (prior, state) in zip(prepared.priors, prepared.states)
            values = eigvals(Hermitian(candidate - prior * state))
            minimum_slack = min(minimum_slack, minimum(real, values))
        end
        trace_value = real(tr(candidate))
        objective_residual = if optimization.dual_objective_value === nothing
            zero(trace_value)
        else
            abs(trace_value - optimization.dual_objective_value)
        end
        score = max(zero(trace_value), -minimum_slack) + objective_residual
        data = (
            operator=Matrix(candidate),
            minimum_slack_eigenvalue=minimum_slack,
            trace=trace_value,
            objective_residual=objective_residual,
            valid=minimum_slack >= -tolerance && objective_residual <= 8tolerance,
            score=score,
        )
        best === nothing || data.score < best.score ? (best = data) : nothing
    end
    return best
end

function _discrimination_solver_result(prepared, problem, optimization, backend)
    tolerance = if backend isa JuMPBackend
        max(
            prepared.tolerance,
            convert(
                eltype(prepared.priors),
                backend.atol +
                backend.rtol * max(
                    one(eltype(prepared.priors)),
                    if optimization.objective_value === nothing
                        zero(eltype(prepared.priors))
                    else
                        abs(optimization.objective_value)
                    end,
                ),
            ),
        )
    else
        prepared.tolerance
    end
    measurement = nothing
    residuals = NamedTuple()
    lower = nothing
    if optimization.primal !== nothing
        names = Tuple(Symbol(:measurement_, index) for index in 1:prepared.state_count)
        candidate = Tuple(
            Matrix(getproperty(optimization.primal.views, name)) for name in names
        )
        checked = _discrimination_measurement_residuals(candidate, prepared, 8tolerance)
        objective_residual = if optimization.objective_value === nothing
            nothing
        else
            abs(checked.success_probability - optimization.objective_value)
        end
        residuals = merge(checked, (solver_objective_residual=objective_residual,))
        if checked.valid &&
            (objective_residual === nothing || objective_residual <= 8tolerance)
            measurement = candidate
            lower = checked.success_probability
        end
    end
    dual = _discrimination_dual_candidate(optimization, prepared, 8tolerance)
    upper = if dual !== nothing && dual.valid
        dual.trace
    elseif optimization.objective_bound !== nothing
        optimization.objective_bound
    else
        nothing
    end
    matched =
        lower !== nothing &&
        upper !== nothing &&
        abs(lower - upper) <= 16tolerance * max(one(tolerance), abs(lower), abs(upper))
    success = matched ? (lower + upper) / 2 : nothing
    if optimization.status === OptimizationBackendUnavailable
        return _discrimination_result(
            StateDiscriminationBackendUnavailable,
            prepared;
            problem=problem,
            optimization_result=optimization,
            residuals=residuals,
            message=optimization.message,
        )
    elseif optimization.status === OptimizationLimit
        return _discrimination_result(
            StateDiscriminationResourceLimit,
            prepared;
            lower_bound=lower,
            upper_bound=upper,
            problem=problem,
            measurement=measurement,
            dual_operator=dual,
            optimization_result=optimization,
            residuals=residuals,
            message="the optimizer stopped at a resource limit; available bounds are retained",
        )
    elseif optimization.status === OptimizationOptimal && measurement !== nothing
        return _discrimination_result(
            StateDiscriminationSolverOptimal,
            prepared;
            success_probability=success,
            lower_bound=lower,
            upper_bound=upper,
            certified=false,
            certificate_kind=nothing,
            problem=problem,
            measurement=measurement,
            dual_operator=dual,
            optimization_result=optimization,
            residuals=residuals,
            warnings=if matched
                ()
            else
                ("the available primal and dual bounds do not coincide within tolerance",)
            end,
            message=if matched
                "the optimizer returned a residual-checked primal-dual optimum"
            else
                "the optimizer reported optimality, but package bounds did not coincide within tolerance"
            end,
        )
    elseif optimization.status === OptimizationFeasible && measurement !== nothing
        return _discrimination_result(
            StateDiscriminationSolverFeasible,
            prepared;
            lower_bound=lower,
            upper_bound=upper,
            problem=problem,
            measurement=measurement,
            dual_operator=dual,
            optimization_result=optimization,
            residuals=residuals,
            message="a residual-checked POVM supplies a lower bound; solver optimality was not established",
        )
    end
    backend_failure = optimization.status in (
        OptimizationMalformedBackend,
        OptimizationNumericalFailure,
        OptimizationUnsupported,
        OptimizationUnknown,
        OptimizationInfeasible,
        OptimizationUnbounded,
    )
    return _discrimination_result(
        if backend_failure
            StateDiscriminationBackendFailure
        else
            StateDiscriminationInvalidCertificate
        end,
        prepared;
        lower_bound=lower,
        upper_bound=upper,
        problem=problem,
        measurement=measurement,
        dual_operator=dual,
        optimization_result=optimization,
        residuals=residuals,
        message="the optimization outcome does not establish a valid optimal discrimination result",
    )
end

"""
    state_distinguishability(states; priors=nothing,
                             backend=NoOptimizationBackend(), ...)

Compute the minimum-error success probability for a validated ensemble.

Two states use the Helstrom theorem without a solver. Exactly orthogonal
supports and deterministic priors also use solver-free certificates. Three or
more nonorthogonal states build a package-owned SDP and require an explicit
optional backend. The returned [`StateDiscriminationResult`](@ref) retains the
POVM, primal/dual bounds, residuals, solver metadata, and inconclusive status.

Unlike the pinned routine, neither density matrices, pure columns, nor priors
are normalized implicitly.
"""
function state_distinguishability(
    states;
    priors=nothing,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dimension=64,
    max_states=64,
    max_variables=100_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _discrimination_prepare(
        states;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dimension=max_dimension,
        max_states=max_states,
        max_variables=max_variables,
    )
    trivial = _discrimination_trivial(prepared)
    trivial === nothing || return trivial
    helstrom = _discrimination_helstrom(prepared)
    helstrom === nothing || return helstrom
    orthogonal = _discrimination_exact_orthogonal(prepared)
    orthogonal === nothing || return orthogonal
    problem = _discrimination_build_problem(prepared; limits=limits)
    optimization = solve_optimization(problem.program, backend)
    return _discrimination_solver_result(prepared, problem, optimization, backend)
end
