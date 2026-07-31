# Package-owned, non-mutating helpers for interpreting status-rich results.
# These helpers deliberately leave each result's native `status` field intact.

"""
    conclusion(result) -> Symbol

Return a conservative, human-facing conclusion for a status-rich result.

For entanglement results the vocabulary is `:entangled`, `:separable`, and
`:unknown`. Passing a necessary criterion and falling outside a sufficient
separable ball both return `:unknown`: neither event decides separability or
entanglement. For [`OptimizationResult`](@ref), usable terminal solver outcomes
map to `:optimal`, `:feasible`, `:infeasible`, or `:unbounded`; limited,
unavailable, inconsistent, failed, and unknown outcomes map to `:unknown`.

The original, more detailed status remains available in `result.status`.
"""
function conclusion end

conclusion(result::EntanglementReport) = result.status

function conclusion(result::CriterionResult)
    result.status === CriterionEntanglementDetected && return :entangled
    return :unknown
end

function conclusion(result::SeparableBallResult)
    result.status === :separable_certified && return :separable
    return :unknown
end

function conclusion(result::OptimizationResult)
    result.status === OptimizationOptimal && return :optimal
    result.status === OptimizationFeasible && return :feasible
    result.status === OptimizationInfeasible && return :infeasible
    result.status === OptimizationUnbounded && return :unbounded
    return :unknown
end

"""
    is_conclusive(result) -> Bool

Return whether [`conclusion`](@ref) contains a resolved domain outcome rather
than `:unknown`. This is deliberately distinct from [`is_certified`](@ref):
for example, a solver may report a usable terminal outcome without supplying a
separately validated mathematical certificate.

For an [`OptimizationResult`](@ref), `true` means only that the solver reported
an optimal, feasible, infeasible, or unbounded terminal outcome. It does not
promote that solver-domain outcome to a mathematical certificate; inspect
[`is_certified`](@ref), residuals, gaps, warnings, and backend metadata
separately.
"""
function is_conclusive end

is_conclusive(result::EntanglementReport) = conclusion(result) !== :unknown
is_conclusive(result::CriterionResult) = conclusion(result) !== :unknown
is_conclusive(result::SeparableBallResult) = conclusion(result) !== :unknown
is_conclusive(result::OptimizationResult) = conclusion(result) !== :unknown

"""
    is_certified(result) -> Bool

Return whether the result carries a package-recognized mathematical
certificate. A robust necessary-criterion violation certifies entanglement,
and membership in the sufficient separable ball certifies separability.
Passing those one-sided tests is not a certificate.
"""
function is_certified end

is_certified(result::EntanglementReport) = result.certified
function is_certified(result::CriterionResult)
    return result.status === CriterionEntanglementDetected
end
is_certified(result::SeparableBallResult) = result.status === :separable_certified
is_certified(result::OptimizationResult) = result.certified

"""
    explain(result) -> String

Return the result's human-readable reason without changing or reclassifying
the result. Optimization warnings are appended so that important residual or
backend qualifications are not hidden.
"""
function explain end

explain(result::EntanglementReport) = result.message
explain(result::CriterionResult) = result.message
explain(result::SeparableBallResult) = result.message

function explain(result::OptimizationResult)
    isempty(result.warnings) && return result.message
    warning_text = join(result.warnings, "; ")
    isempty(result.message) && return "Warnings: $warning_text"
    return "$(result.message) Warnings: $warning_text"
end

_result_plain_text(value::Nothing) = "none"
_result_plain_text(value::Bool) = value ? "yes" : "no"
_result_plain_text(value::Symbol) = replace(String(value), '_' => ' ')
_result_plain_text(value) = string(value)

function _show_result_field(io::IO, label::AbstractString, value)
    return print(io, "\n  ", label, ": ", _result_plain_text(value))
end

function _show_result_certification(
    io::IO, certified::Bool, certificate_kind::Union{Nothing,Symbol}=nothing
)
    _show_result_field(io, "Certified", certified)
    if certified && certificate_kind !== nothing
        _show_result_field(io, "Certificate", certificate_kind)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", report::EntanglementReport)
    print(io, "EntanglementReport")
    _show_result_field(io, "Conclusion", conclusion(report))
    _show_result_certification(io, is_certified(report), report.certificate_kind)
    _show_result_field(io, "Method", report.method)
    _show_result_field(io, "Backend", report.backend)
    _show_result_field(io, "Explanation", explain(report))
    _show_result_field(io, "Attempts", length(report.attempts))
    for (index, attempt) in pairs(report.attempts)
        print(
            io,
            "\n    ",
            index,
            ". ",
            _result_plain_text(attempt.method),
            " via ",
            _result_plain_text(attempt.backend),
            ": ",
            _result_plain_text(attempt.status),
            " (certified: ",
            _result_plain_text(attempt.certified),
            ")",
        )
        isempty(attempt.message) ||
            print(io, "\n       Reason: ", _result_plain_text(attempt.message))
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", result::CriterionResult)
    print(io, "CriterionResult")
    _show_result_field(io, "Criterion", result.criterion)
    _show_result_field(io, "Conclusion", conclusion(result))
    _show_result_certification(
        io, is_certified(result), is_certified(result) ? :criterion_violation : nothing
    )
    _show_result_field(io, "Criterion status", result.status)
    _show_result_field(io, "Value", result.value)
    _show_result_field(io, "Threshold", result.threshold)
    _show_result_field(io, "Tolerance", result.tolerance)
    _show_result_field(io, "Explanation", explain(result))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", result::SeparableBallResult)
    print(io, "SeparableBallResult")
    _show_result_field(io, "Conclusion", conclusion(result))
    _show_result_certification(
        io,
        is_certified(result),
        is_certified(result) ? :sufficient_separable_ball : nothing,
    )
    _show_result_field(io, "Ball status", result.status)
    _show_result_field(io, "Purity", result.purity)
    _show_result_field(io, "Boundary", result.boundary)
    _show_result_field(io, "Tolerance", result.tolerance)
    _show_result_field(io, "Explanation", explain(result))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", result::OptimizationResult)
    print(io, "OptimizationResult")
    _show_result_field(io, "Conclusion", conclusion(result))
    _show_result_certification(io, is_certified(result), result.certificate_kind)
    _show_result_field(io, "Optimization status", result.status)
    _show_result_field(io, "Termination", result.termination_status)
    _show_result_field(io, "Primal status", result.primal_status)
    _show_result_field(io, "Dual status", result.dual_status)
    result.objective_value === nothing ||
        _show_result_field(io, "Objective value", result.objective_value)
    result.objective_bound === nothing ||
        _show_result_field(io, "Objective bound", result.objective_bound)
    result.absolute_gap === nothing ||
        _show_result_field(io, "Absolute gap", result.absolute_gap)
    result.relative_gap === nothing ||
        _show_result_field(io, "Relative gap", result.relative_gap)
    result.primal_residual === nothing ||
        _show_result_field(io, "Primal residual", result.primal_residual)
    result.dual_residual === nothing ||
        _show_result_field(io, "Dual residual", result.dual_residual)
    _show_result_field(io, "Explanation", result.message)
    if !isempty(result.warnings)
        print(io, "\n  Warnings:")
        for warning in result.warnings
            print(io, "\n    - ", warning)
        end
    end
    return nothing
end
