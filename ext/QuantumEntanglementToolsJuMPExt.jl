module QuantumEntanglementToolsJuMPExt

using JuMP: JuMP
using LinearAlgebra: LinearAlgebra
using SparseArrays: SparseArrays
import QuantumEntanglementTools as QET

const MOI = JuMP.MOI

function _jump_affine(
    function_data::QET.AffineScalar{T}, variables::AbstractVector{<:JuMP.AbstractJuMPScalar}
) where {T<:Real}
    pairs = Pair{eltype(variables),T}[]
    indices, coefficients = SparseArrays.findnz(function_data.coefficients)
    sizehint!(pairs, length(indices))
    for (index, coefficient) in zip(indices, coefficients)
        push!(pairs, variables[index] => coefficient)
    end
    return JuMP.GenericAffExpr(function_data.constant, pairs)
end

function _jump_real_block(
    matrix::QET.HermitianAffineMatrix{T},
    variables::AbstractVector{<:JuMP.AbstractJuMPScalar},
) where {T<:Real}
    dimension = matrix.dimension
    expressions = Matrix{JuMP.GenericAffExpr{T,eltype(variables)}}(
        undef, 2dimension, 2dimension
    )
    real_constant = QET.real_block_embedding(matrix.constant)
    for column in 1:(2dimension), row in 1:(2dimension)
        expressions[row, column] = JuMP.GenericAffExpr(
            convert(T, real_constant[row, column]), Pair{eltype(variables),T}[]
        )
    end
    for term in matrix.terms
        block = QET.real_block_embedding(term.coefficient)
        variable = variables[term.variable]
        rows, columns, values = SparseArrays.findnz(SparseArrays.sparse(block))
        for (row, column, coefficient) in zip(rows, columns, values)
            JuMP.add_to_expression!(
                expressions[row, column], convert(T, coefficient), variable
            )
        end
    end
    return LinearAlgebra.Symmetric(expressions)
end

function _set_optimizer_options!(model, backend::QET.JuMPBackend)
    backend.silent && JuMP.set_silent(model)
    backend.time_limit_seconds === nothing ||
        JuMP.set_time_limit_sec(model, backend.time_limit_seconds)
    for (name, value) in pairs(backend.optimizer_options)
        JuMP.set_optimizer_attribute(model, String(name), value)
    end
    return nothing
end

function _build_model(
    problem::QET.SemidefiniteProgram{T}, backend::QET.JuMPBackend
) where {T<:Real}
    !isempty(problem.psd_constraints) &&
        !backend.allow_densify &&
        throw(
            ArgumentError(
                "the JuMP real-block SDP translation allocates dense expression " *
                "matrices; construct JuMPBackend with allow_densify=true",
            ),
        )
    model = JuMP.GenericModel{T}(backend.optimizer_factory)
    _set_optimizer_options!(model, backend)
    variables = JuMP.@variable(model, [1:(problem.variable_count)])
    if problem.initial_point !== nothing
        for index in eachindex(variables)
            JuMP.set_start_value(variables[index], problem.initial_point[index])
        end
    end

    objective = _jump_affine(problem.objective, variables)
    if problem.sense === :minimize
        JuMP.set_objective_sense(model, MOI.MIN_SENSE)
        JuMP.set_objective_function(model, objective)
    elseif problem.sense === :maximize
        JuMP.set_objective_sense(model, MOI.MAX_SENSE)
        JuMP.set_objective_function(model, objective)
    else
        JuMP.set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    end

    equality_refs = JuMP.ConstraintRef[]
    for constraint in problem.equalities
        expression = _jump_affine(constraint.function_data, variables)
        push!(equality_refs, JuMP.@constraint(model, expression == zero(T)))
    end

    interval_refs = JuMP.ConstraintRef[]
    for constraint in problem.intervals
        expression = _jump_affine(constraint.function_data, variables)
        reference = if constraint.lower === nothing
            JuMP.@constraint(model, expression <= constraint.upper)
        elseif constraint.upper === nothing
            JuMP.@constraint(model, expression >= constraint.lower)
        else
            JuMP.@constraint(model, constraint.lower <= expression <= constraint.upper)
        end
        push!(interval_refs, reference)
    end

    psd_refs = JuMP.ConstraintRef[]
    for matrix in problem.psd_constraints
        real_block = _jump_real_block(matrix, variables)
        push!(psd_refs, JuMP.@constraint(model, real_block in JuMP.PSDCone()))
    end
    return (; model, variables, equality_refs, interval_refs, psd_refs)
end

_status_symbol(status) = Symbol(lowercase(string(status)))

function _optional_query(query, model)
    return try
        value = query(model)
        value === missing ? nothing : value
    catch
        nothing
    end
end

function _finite_real_or_nothing(value, ::Type{T}) where {T<:Real}
    value === nothing && return nothing
    value isa Real && !(value isa Bool) && isfinite(value) || return nothing
    return try
        convert(T, value)
    catch
        nothing
    end
end

function _metadata(
    problem::QET.SemidefiniteProgram{T},
    backend::QET.JuMPBackend,
    model,
    raw_status::AbstractString,
) where {T<:Real}
    reported_name = try
        String(JuMP.solver_name(model))
    catch
        "unavailable"
    end
    return QET.OptimizerMetadata(
        :JuMP,
        Base.pkgversion(JuMP),
        Base.pkgversion(MOI),
        backend.optimizer_name,
        reported_name,
        backend.optimizer_version,
        backend.optimizer_options,
        problem.limits,
        T,
        isempty(problem.psd_constraints) ? :not_used : :real_block,
        String(raw_status),
    )
end

function _exception_status(error)
    if error isa Union{
        MOI.UnsupportedConstraint,
        MOI.UnsupportedAttribute,
        MOI.AddConstraintNotAllowed,
        MOI.AddVariableNotAllowed,
    }
        return (
            QET.OptimizationUnsupported,
            :unsupported_model,
            "the optimizer does not support the required public MOI operation",
        )
    elseif error isa ArgumentError &&
        occursin("allow_densify=true", sprint(showerror, error))
        return (
            QET.OptimizationUnsupported,
            :densification_not_allowed,
            sprint(showerror, error),
        )
    end
    return (
        QET.OptimizationMalformedBackend,
        :backend_exception,
        "the optimizer factory, options, or solve call failed",
    )
end

function _failure_from_exception(
    problem::QET.SemidefiniteProgram, backend::QET.JuMPBackend, error; model=nothing
)
    status, termination, prefix = _exception_status(error)
    detail = sprint(showerror, error)
    raw_status = "$prefix: $detail"
    metadata = if model === nothing
        QET.OptimizerMetadata(
            :JuMP,
            Base.pkgversion(JuMP),
            Base.pkgversion(MOI),
            backend.optimizer_name,
            "unavailable",
            backend.optimizer_version,
            backend.optimizer_options,
            problem.limits,
            eltype(problem.objective.coefficients),
            isempty(problem.psd_constraints) ? :not_used : :real_block,
            raw_status,
        )
    else
        _metadata(problem, backend, model, raw_status)
    end
    return QET._optimization_failure_result(
        problem,
        status,
        termination,
        raw_status;
        optimizer=metadata,
        warnings=("backend exception type: $(typeof(error))",),
    )
end

function _proposed_status(
    termination, primal_status, dual_status, known_feasible_point_verified::Bool
)
    if termination === MOI.OPTIMAL
        return if primal_status === MOI.FEASIBLE_POINT && dual_status === MOI.FEASIBLE_POINT
            QET.OptimizationOptimal
        else
            QET.OptimizationInconsistent
        end
    elseif termination in
        (MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
        return if primal_status in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
            QET.OptimizationFeasible
        else
            QET.OptimizationInconsistent
        end
    elseif termination in (
        MOI.TIME_LIMIT,
        MOI.ITERATION_LIMIT,
        MOI.NODE_LIMIT,
        MOI.SOLUTION_LIMIT,
        MOI.MEMORY_LIMIT,
        MOI.OBJECTIVE_LIMIT,
        MOI.NORM_LIMIT,
        MOI.OTHER_LIMIT,
    )
        return QET.OptimizationLimit
    elseif termination in
        (MOI.NUMERICAL_ERROR, MOI.SLOW_PROGRESS, MOI.INVALID_MODEL, MOI.OTHER_ERROR)
        return QET.OptimizationNumericalFailure
    elseif termination in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
        return if dual_status === MOI.INFEASIBILITY_CERTIFICATE
            QET.OptimizationInfeasible
        else
            QET.OptimizationUnknown
        end
    elseif termination in (MOI.DUAL_INFEASIBLE, MOI.ALMOST_DUAL_INFEASIBLE)
        if primal_status === MOI.INFEASIBILITY_CERTIFICATE && known_feasible_point_verified
            return QET.OptimizationUnbounded
        end
        return QET.OptimizationUnknown
    end
    return QET.OptimizationUnknown
end

function _extract_primal(problem::QET.SemidefiniteProgram{T}, built) where {T<:Real}
    JuMP.has_values(built.model) || return nothing
    coordinates = T[JuMP.value(variable) for variable in built.variables]
    all(isfinite, coordinates) || return nothing
    names = Tuple(matrix.name for matrix in problem.primal_views)
    values = Tuple(
        QET.evaluate_affine(matrix, coordinates) for matrix in problem.primal_views
    )
    views = NamedTuple{names}(values)
    return QET.OptimizationPrimal(coordinates, views)
end

function _extract_dual(problem::QET.SemidefiniteProgram{T}, built) where {T<:Real}
    JuMP.has_duals(built.model) || return nothing
    equality_duals = T[JuMP.dual(reference) for reference in built.equality_refs]
    interval_duals = T[JuMP.dual(reference) for reference in built.interval_refs]
    real_blocks = Tuple(Matrix{T}(JuMP.dual(reference)) for reference in built.psd_refs)
    hermitian_blocks = Tuple(
        Matrix(QET.hermitian_dual_from_real_block(block)) for block in real_blocks
    )
    return QET.OptimizationDual(
        equality_duals, interval_duals, real_blocks, hermitian_blocks
    )
end

function _objective_gap(objective_value, dual_objective, ::Type{T}) where {T}
    objective_value === nothing && return nothing
    dual_objective === nothing && return nothing
    return convert(T, abs(objective_value - dual_objective))
end

# This generic package-owned SDP translation is also the optional backend path
# used by `abs_ppt_constraints` affine models, `is_abs_ppt` ordering
# realization, `symmetric_extension`, and `symmetric_inner_extension`; those
# APIs independently validate any solver-derived evidence.
function QET._solve_optional_optimization(
    problem::QET.SemidefiniteProgram{T}, backend::QET.JuMPBackend
) where {T<:Real}
    built = try
        _build_model(problem, backend)
    catch error
        return _failure_from_exception(problem, backend, error)
    end
    try
        JuMP.optimize!(built.model)
    catch error
        return _failure_from_exception(problem, backend, error; model=built.model)
    end

    termination = JuMP.termination_status(built.model)
    primal_status = JuMP.primal_status(built.model)
    dual_status = JuMP.dual_status(built.model)
    raw_status = try
        String(JuMP.raw_status(built.model))
    catch
        string(termination)
    end
    metadata = _metadata(problem, backend, built.model, raw_status)
    primal = try
        _extract_primal(problem, built)
    catch error
        return _failure_from_exception(problem, backend, error; model=built.model)
    end
    dual = try
        _extract_dual(problem, built)
    catch
        nothing
    end

    objective_value = _finite_real_or_nothing(
        _optional_query(JuMP.objective_value, built.model), T
    )
    objective_bound = _finite_real_or_nothing(
        _optional_query(JuMP.objective_bound, built.model), T
    )
    dual_objective = _finite_real_or_nothing(
        _optional_query(JuMP.dual_objective_value, built.model), T
    )
    objective_bound === nothing &&
        dual_objective !== nothing &&
        (objective_bound = dual_objective)
    absolute_gap = _objective_gap(objective_value, dual_objective, T)
    relative_gap = _finite_real_or_nothing(
        _optional_query(JuMP.relative_gap, built.model), T
    )
    if relative_gap === nothing && absolute_gap !== nothing
        scale = max(
            one(T),
            objective_value === nothing ? zero(T) : abs(objective_value),
            dual_objective === nothing ? zero(T) : abs(dual_objective),
        )
        relative_gap = absolute_gap / scale
    end
    primal_residual_value = if primal === nothing
        nothing
    else
        try
            convert(
                T,
                QET.primal_residual(
                    problem, primal.coordinates; allow_densify=backend.allow_densify
                ),
            )
        catch
            nothing
        end
    end
    dual_residual_value = if dual === nothing
        nothing
    else
        try
            convert(
                T,
                QET.dual_stationarity_residual(
                    problem,
                    dual.equalities,
                    dual.intervals,
                    collect(dual.psd_real_blocks),
                ),
            )
        catch
            nothing
        end
    end
    tolerance = QET._optimization_residual_tolerance(problem, backend, objective_value)
    known_feasible_point_verified = if problem.known_feasible_point === nothing
        false
    else
        try
            QET.primal_residual(
                problem,
                problem.known_feasible_point;
                allow_densify=backend.allow_densify,
            ) <= tolerance
        catch
            false
        end
    end
    proposed = _proposed_status(
        termination, primal_status, dual_status, known_feasible_point_verified
    )
    proposed in (QET.OptimizationOptimal, QET.OptimizationFeasible) &&
        primal === nothing &&
        (proposed = QET.OptimizationInconsistent)
    proposed === QET.OptimizationOptimal &&
        dual === nothing &&
        (proposed = QET.OptimizationInconsistent)
    status = QET._optimization_consistent_status(
        proposed, primal_residual_value, dual_residual_value, tolerance
    )

    iterations = _optional_query(JuMP.barrier_iterations, built.model)
    iterations = iterations isa Integer && iterations >= 0 ? Int(iterations) : nothing
    solve_time = _optional_query(JuMP.solve_time, built.model)
    solve_time = if solve_time isa Real && isfinite(solve_time) && solve_time >= 0
        Float64(solve_time)
    else
        nothing
    end
    message = if status === QET.OptimizationOptimal
        "the optimizer reported an optimal primal-dual result within package residual tolerances"
    elseif status === QET.OptimizationFeasible
        "the optimizer returned a feasible or inaccurate solution; it is not an exact certificate"
    elseif status === QET.OptimizationInfeasible
        "the optimizer returned an infeasibility certificate"
    elseif status === QET.OptimizationUnbounded
        "the optimizer returned a primal ray and the model carries an independently checked feasible point"
    elseif status === QET.OptimizationLimit
        "the optimizer stopped at a resource limit"
    elseif status === QET.OptimizationNumericalFailure
        "the optimizer reported a numerical or model failure"
    elseif status === QET.OptimizationInconsistent
        "the solver status conflicts with package-owned primal or dual residual checks"
    else
        "the optimizer outcome is inconclusive"
    end
    warning_values = String[]
    proposed === QET.OptimizationOptimal &&
        dual === nothing &&
        push!(warning_values, "the optimizer supplied no dual solution")
    status === QET.OptimizationInconsistent && push!(
        warning_values,
        "at least one residual exceeds atol + rtol * max(1, |objective|)",
    )
    return QET.OptimizationResult{T,typeof(primal),typeof(dual),typeof(metadata)}(
        status,
        _status_symbol(termination),
        _status_symbol(primal_status),
        _status_symbol(dual_status),
        objective_value,
        objective_bound,
        dual_objective,
        absolute_gap,
        relative_gap,
        primal_residual_value,
        dual_residual_value,
        iterations,
        solve_time,
        false,
        nothing,
        primal,
        dual,
        metadata,
        Tuple(warning_values),
        message,
    )
end

include("psd_constraints.jl")
include("top_k_p_norm_epigraph.jl")
include("top_k_p_norm_dual_epigraph.jl")
include("polynomial_sos.jl")
include("copositivity_clique.jl")
include("sk_operator_norm.jl")
include("block_positivity.jl")
include("separability_optimization.jl")
include("nonlocal_games.jl")

# Loading the JuMP extension permits an explicit positional-backend form while
# retaining the solver-free keyword-only method in the core package.
function QET.state_distinguishability(states, backend::QET.JuMPBackend; kwargs...)
    return QET.state_distinguishability(states; backend=backend, kwargs...)
end

function QET.diamond_norm(map, backend::QET.JuMPBackend; kwargs...)
    return QET.diamond_norm(map; backend=backend, kwargs...)
end

function QET.cb_norm(map, backend::QET.JuMPBackend; kwargs...)
    return QET.cb_norm(map; backend=backend, kwargs...)
end

function QET.channel_distinguishability(
    first_map, second_map, backend::QET.JuMPBackend; kwargs...
)
    return QET.channel_distinguishability(first_map, second_map; backend=backend, kwargs...)
end

function QET.maximum_output_fidelity(
    first_map, second_map, backend::QET.JuMPBackend; kwargs...
)
    return QET.maximum_output_fidelity(first_map, second_map; backend=backend, kwargs...)
end

function QET.is_k_incoherent(state, k, backend::QET.JuMPBackend; kwargs...)
    return QET.is_k_incoherent(state, k; backend=backend, kwargs...)
end

function QET.is_absolutely_k_incoherent(state, k, backend::QET.JuMPBackend; kwargs...)
    return QET.is_absolutely_k_incoherent(state, k; backend=backend, kwargs...)
end

function QET.robustness_coherence(state, backend::QET.JuMPBackend; kwargs...)
    return QET.robustness_coherence(state; backend=backend, kwargs...)
end

function QET.trace_distance_coherence(state, backend::QET.JuMPBackend; kwargs...)
    return QET.trace_distance_coherence(state; backend=backend, kwargs...)
end

function QET.generalized_robustness_k_coherence(
    state, k, backend::QET.JuMPBackend; kwargs...
)
    return QET.generalized_robustness_k_coherence(state, k; backend=backend, kwargs...)
end

end # module QuantumEntanglementToolsJuMPExt
