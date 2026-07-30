# Independently designed Julia-native orchestration layer. The mathematical
# criteria it invokes carry their QETLAB provenance in criteria.jl and
# PROVENANCE.toml.

"""
    AbstractEntanglementMethod

Supertype for package-owned entanglement-analysis method configurations.
Optional backends extend [`detect_entanglement`](@ref) by dispatching on a
concrete subtype without changing core result semantics.
"""
abstract type AbstractEntanglementMethod end

"""
    AbstractEntanglementBackend

Supertype for backend descriptors returned by
[`available_entanglement_backends`](@ref).
"""
abstract type AbstractEntanglementBackend end

"""
    NativeEntanglementBackend

Descriptor for the dependency-free criteria implemented by this package.
"""
struct NativeEntanglementBackend <: AbstractEntanglementBackend end

"""
    EntanglementDetectionBackend

Descriptor for the optional EntanglementDetection.jl adapter. The descriptor
is returned by [`available_entanglement_backends`](@ref) only after
EntanglementDetection.jl has loaded and activated the package extension.
"""
struct EntanglementDetectionBackend <: AbstractEntanglementBackend end

abstract type AbstractEntanglementDetectionMethod <: AbstractEntanglementMethod end

"""
    EntanglementDetectionSearch(;
        timeout_seconds=120,
        max_iteration=10_000,
        epsilon=1e-6,
        callback_iter=10_000,
        atol=0,
        rtol=sqrt(eps(Float64)),
        allow_densify=false,
    )

Configure the optional EntanglementDetection.jl 0.2.2 heuristic search.
Execution always occurs in a fresh child Julia process because that backend
version seeds its default RNG and has options that can redirect stdout or
change the process-wide BLAS thread count. `timeout_seconds` is a finite,
strictly positive wall-clock limit; a timeout, child-process failure, or
backend exception is returned as an uncertified `:unknown`
[`EntanglementReport`](@ref).

`max_iteration`, `epsilon`, and `callback_iter` are passed to the documented
public `EntanglementDetection.entanglement_detection` entry point. The adapter
fixes backend verbosity to zero and does not accept a logfile. Input density
matrices use the package's normal finite/Hermitian/trace/positivity validation.
Sparse input requires the explicit `allow_densify=true` opt-in.

EntanglementDetection.jl's Boolean-or-`nothing` conclusion is retained only as
package-owned candidate evidence. It is never upgraded by this adapter to a
certified entanglement or separability conclusion.
"""
struct EntanglementDetectionSearch{T<:Real} <: AbstractEntanglementDetectionMethod
    timeout_seconds::Float64
    max_iteration::Int
    epsilon::T
    callback_iter::Int
    atol::T
    rtol::T
    allow_densify::Bool
end

function EntanglementDetectionSearch(;
    timeout_seconds=120,
    max_iteration=10_000,
    epsilon=1e-6,
    callback_iter=10_000,
    atol=0,
    rtol=sqrt(eps(Float64)),
    allow_densify::Bool=false,
)
    timeout_seconds isa Real && !(timeout_seconds isa Bool) ||
        throw(ArgumentError("timeout_seconds must be a finite positive real number"))
    isfinite(timeout_seconds) && timeout_seconds > zero(timeout_seconds) ||
        throw(ArgumentError("timeout_seconds must be a finite positive real number"))
    timeout = try
        Float64(timeout_seconds)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("timeout_seconds cannot be represented as Float64"))
    end
    isfinite(timeout) && timeout > 0 || throw(
        ArgumentError(
            "timeout_seconds must remain finite and positive when represented as Float64",
        ),
    )

    iterations = _positive_int(max_iteration, "max_iteration")
    callback = _positive_int(callback_iter, "callback_iter")
    checked_epsilon = _tierd_validate_tolerance(epsilon, "epsilon")
    checked_epsilon === nothing &&
        throw(ArgumentError("epsilon must be a finite nonnegative real number"))
    absolute = _tierd_validate_tolerance(atol, "atol")
    absolute === nothing &&
        throw(ArgumentError("atol must be a finite nonnegative real number"))
    relative = _tierd_validate_tolerance(rtol, "rtol")
    relative === nothing &&
        throw(ArgumentError("rtol must be a finite nonnegative real number"))
    promoted_epsilon, promoted_absolute, promoted_relative = promote(
        checked_epsilon, absolute, relative
    )
    return EntanglementDetectionSearch(
        timeout,
        iterations,
        promoted_epsilon,
        callback,
        promoted_absolute,
        promoted_relative,
        allow_densify,
    )
end

"""
    NativePPT(; systems=(2,), atol=0.0, rtol=sqrt(eps(Float64)),
                allow_densify=false)

Configure the native positive-partial-transpose method. `systems` selects a
nonempty proper subset of subsystem indices. Tolerances are finite,
nonnegative, and passed directly to [`ppt_criterion`](@ref). Sparse spectral
work requires the explicit `allow_densify=true` opt-in.
"""
struct NativePPT{T<:Real,S} <: AbstractEntanglementMethod
    systems::S
    atol::T
    rtol::T
    allow_densify::Bool
end

function NativePPT(;
    systems=(2,), atol=0.0, rtol=sqrt(eps(Float64)), allow_densify::Bool=false
)
    absolute = _tierd_validate_tolerance(atol, "atol")
    relative = _tierd_validate_tolerance(rtol, "rtol")
    absolute === nothing &&
        throw(ArgumentError("atol must be a finite nonnegative real number"))
    relative === nothing &&
        throw(ArgumentError("rtol must be a finite nonnegative real number"))
    promoted_absolute, promoted_relative = promote(absolute, relative)
    return NativePPT(systems, promoted_absolute, promoted_relative, allow_densify)
end

"""
    EntanglementAttempt

One recorded step in an [`EntanglementReport`](@ref). `status` is one of
`:entangled`, `:separable`, or `:unknown`. `certified` is true only when the
attempt supplies a mathematically valid certificate in the stated domain.
`raw_result` is package-owned criterion/decomposition data, never an
undocumented external backend object. Package-certified attempts are produced
only by validated internal constructors; the public constructor is available
for uncertified user annotations.
"""
struct EntanglementAttempt
    method::Symbol
    backend::Symbol
    status::Symbol
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    raw_result::Any
    message::String
    _package_validated::Bool

    function EntanglementAttempt(
        method::Symbol,
        backend::Symbol,
        status::Symbol,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        raw_result,
        message::AbstractString,
    )
        status in (:entangled, :separable, :unknown) || throw(
            ArgumentError("attempt status must be :entangled, :separable, or :unknown")
        )
        certified && throw(
            ArgumentError(
                "package-certified attempts are produced by analysis routines, not by the public constructor",
            ),
        )
        certificate_kind === nothing ||
            throw(ArgumentError("an uncertified attempt cannot name a certificate_kind"))
        return new(
            method,
            backend,
            status,
            false,
            nothing,
            deepcopy(raw_result),
            String(message),
            false,
        )
    end

    function EntanglementAttempt(
        token::_ValidatedConstructorToken,
        method::Symbol,
        backend::Symbol,
        status::Symbol,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        raw_result,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        status in (:entangled, :separable, :unknown) || throw(
            ArgumentError("attempt status must be :entangled, :separable, or :unknown")
        )
        if certified
            status in (:entangled, :separable) ||
                throw(ArgumentError("a certified attempt must be entangled or separable"))
            certificate_kind === nothing &&
                throw(ArgumentError("a certified attempt must name its certificate_kind"))
            raw_result === nothing &&
                throw(ArgumentError("a certified attempt must retain its evidence"))
        else
            certificate_kind === nothing || throw(
                ArgumentError("an uncertified attempt cannot name a certificate_kind")
            )
        end
        return new(
            method,
            backend,
            status,
            certified,
            certificate_kind,
            _entanglement_owned_evidence(raw_result),
            String(message),
            true,
        )
    end
end

_entanglement_owned_evidence(::Nothing) = nothing
_entanglement_owned_evidence(value::_ReadOnlyPlanArray) = _read_only_plan_array(copy(value))
_entanglement_owned_evidence(value::Array{<:Number}) = _read_only_plan_array(value)
function _entanglement_owned_evidence(value::Diagonal)
    return Diagonal(_read_only_plan_array(Vector(diag(value))))
end
function _entanglement_owned_evidence(value::Tuple)
    return map(_entanglement_owned_evidence, value)
end
function _entanglement_owned_evidence(value::NamedTuple{names}) where {names}
    copied = Tuple(_entanglement_owned_evidence(entry) for entry in values(value))
    return NamedTuple{names}(copied)
end
function _entanglement_owned_evidence(value::CriterionResult)
    return CriterionResult(
        value.criterion,
        value.status,
        value.value,
        value.threshold,
        value.tolerance,
        _entanglement_owned_evidence(value.witness),
        value.message,
    )
end
function _entanglement_owned_evidence(value::SchmidtDecompositionResult)
    return SchmidtDecompositionResult(
        _entanglement_owned_evidence(value.coefficients),
        _entanglement_owned_evidence(value.left_vectors),
        _entanglement_owned_evidence(value.right_vectors),
    )
end
_entanglement_owned_evidence(value) = deepcopy(value)

function _validated_entanglement_attempt(
    method::Symbol,
    backend::Symbol,
    status::Symbol,
    certified::Bool,
    certificate_kind::Union{Nothing,Symbol},
    raw_result,
    message::AbstractString,
)
    return EntanglementAttempt(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        method,
        backend,
        status,
        certified,
        certificate_kind,
        raw_result,
        message,
    )
end

function _is_package_validated_entanglement_attempt(attempt::EntanglementAttempt)
    return getfield(attempt, :_package_validated)
end

"""
    EntanglementReport

Structured high-level conclusion. `status` is `:entangled`, `:separable`, or
`:unknown`; `certified` is never true without a named `certificate_kind`.
`attempts` records every method run in order in owned immutable tuple storage.
A necessary criterion that merely passes is represented as `:unknown`, except
where a separately stated theorem makes it sufficient. Package-certified
reports are produced only by analysis routines.
"""
struct EntanglementReport
    status::Symbol
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    method::Symbol
    backend::Symbol
    evidence::Any
    attempts::Tuple{Vararg{EntanglementAttempt}}
    message::String
    _package_validated::Bool

    function EntanglementReport(
        status::Symbol,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        method::Symbol,
        backend::Symbol,
        evidence,
        attempts::Union{Tuple,AbstractVector},
        message::AbstractString,
    )
        status in (:entangled, :separable, :unknown) || throw(
            ArgumentError("report status must be :entangled, :separable, or :unknown")
        )
        all(attempt -> attempt isa EntanglementAttempt, attempts) ||
            throw(ArgumentError("attempts must contain only EntanglementAttempt values"))
        certified && throw(
            ArgumentError(
                "package-certified reports are produced by analysis routines, not by the public constructor",
            ),
        )
        certificate_kind === nothing ||
            throw(ArgumentError("an uncertified report cannot name a certificate_kind"))
        return new(
            status,
            false,
            nothing,
            method,
            backend,
            deepcopy(evidence),
            Tuple(attempts),
            String(message),
            false,
        )
    end

    function EntanglementReport(
        token::_ValidatedConstructorToken,
        status::Symbol,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        method::Symbol,
        backend::Symbol,
        evidence,
        attempts::Union{Tuple,AbstractVector},
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        status in (:entangled, :separable, :unknown) || throw(
            ArgumentError("report status must be :entangled, :separable, or :unknown")
        )
        all(attempt -> attempt isa EntanglementAttempt, attempts) ||
            throw(ArgumentError("attempts must contain only EntanglementAttempt values"))
        owned_attempts = Tuple(attempts)
        if certified
            status in (:entangled, :separable) ||
                throw(ArgumentError("a certified report must be entangled or separable"))
            certificate_kind === nothing &&
                throw(ArgumentError("a certified report must name its certificate_kind"))
            evidence === nothing &&
                throw(ArgumentError("a certified report must retain its evidence"))
            matching_attempt = any(owned_attempts) do attempt
                return _is_package_validated_entanglement_attempt(attempt) &&
                       attempt.certified &&
                       attempt.status === status &&
                       attempt.certificate_kind === certificate_kind &&
                       attempt.method === method &&
                       attempt.backend === backend
            end
            matching_attempt || throw(
                ArgumentError(
                    "a certified report requires a matching package-validated attempt"
                ),
            )
        else
            certificate_kind === nothing ||
                throw(ArgumentError("an uncertified report cannot name a certificate_kind"))
        end
        return new(
            status,
            certified,
            certificate_kind,
            method,
            backend,
            _entanglement_owned_evidence(evidence),
            owned_attempts,
            String(message),
            true,
        )
    end
end

function _validated_entanglement_report(
    status::Symbol,
    certified::Bool,
    certificate_kind::Union{Nothing,Symbol},
    method::Symbol,
    backend::Symbol,
    evidence,
    attempts,
    message::AbstractString,
)
    return EntanglementReport(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        certified,
        certificate_kind,
        method,
        backend,
        evidence,
        attempts,
        message,
    )
end

function _is_package_validated_entanglement_report(report::EntanglementReport)
    getfield(report, :_package_validated) || return false
    report.certified || return report.certificate_kind === nothing
    return any(report.attempts) do attempt
        return _is_package_validated_entanglement_attempt(attempt) &&
               attempt.certified &&
               attempt.status === report.status &&
               attempt.certificate_kind === report.certificate_kind &&
               attempt.method === report.method &&
               attempt.backend === report.backend
    end
end

function Base.show(io::IO, attempt::EntanglementAttempt)
    return print(
        io,
        "EntanglementAttempt(",
        attempt.method,
        ", status=",
        attempt.status,
        ", certified=",
        attempt.certified,
        ")",
    )
end

function Base.show(io::IO, report::EntanglementReport)
    return print(
        io,
        "EntanglementReport(status=",
        report.status,
        ", certified=",
        report.certified,
        ", method=",
        report.method,
        ", attempts=",
        length(report.attempts),
        ")",
    )
end

"""
    backend_capabilities(backend)

Return stable package-owned metadata for an entanglement backend.
"""
function backend_capabilities(::NativeEntanglementBackend)
    return (
        name=:native,
        version=v"0.1.0",
        methods=(:ppt, :realignment, :reduction, :pure_schmidt),
        conclusions=(:entangled, :separable, :unknown),
        side_effect_free=true,
        optional_dependency=false,
    )
end

function _entanglement_detection_extension()
    return Base.get_extension(
        @__MODULE__, :QuantumEntanglementToolsEntanglementDetectionExt
    )
end

function backend_capabilities(::EntanglementDetectionBackend)
    extension = _entanglement_detection_extension()
    return (
        name=:entanglement_detection,
        version=isnothing(extension) ? nothing : extension.backend_version(),
        methods=(:heuristic_search,),
        conclusions=(:unknown,),
        backend_candidates=(:entangled, :separable, :inconclusive),
        side_effect_free=false,
        caller_state_isolated=true,
        isolation=:child_process,
        optional_dependency=true,
        minimum_resolvable_julia=v"1.11.0",
        loaded=(!isnothing(extension)),
        certifies_conclusions=false,
        source_integrity_enforced=false,
        transport_trust=:same_version_local_worker,
        resource_sandboxed=false,
    )
end

"""
    available_entanglement_backends()

Return the backends currently registered in the dependency-free core. Optional
extensions may add explicit backend descriptors but are never loaded merely by
calling this function.
"""
function available_entanglement_backends()
    native = NativeEntanglementBackend()
    isnothing(_entanglement_detection_extension()) && return (native,)
    return (native, EntanglementDetectionBackend())
end

"""
    detect_entanglement(state, dims, method)

Run one explicitly selected entanglement method and return an
[`EntanglementReport`](@ref). Optional packages extend this generic function by
dispatch; absence or failure of a backend must not be mapped to a mathematical
negative result.
"""
function detect_entanglement end

function detect_entanglement(
    ::AbstractMatrix{<:Number}, dims, ::AbstractEntanglementDetectionMethod
)
    _as_layout(dims)
    return throw(
        ArgumentError(
            "EntanglementDetection.jl is not loaded; load it explicitly with " *
            "`using EntanglementDetection` before using EntanglementDetectionSearch",
        ),
    )
end

function _ppt_low_dimension_separability_domain(layout::SubsystemLayout)
    length(layout) == 2 || return false
    first_dimension, second_dimension = layout.dims
    return (first_dimension == 2 && second_dimension in (2, 3)) ||
           (second_dimension == 2 && first_dimension in (2, 3))
end

function _report_from_attempt(attempt::EntanglementAttempt)
    _is_package_validated_entanglement_attempt(attempt) ||
        throw(ArgumentError("reports may only promote package-validated attempts"))
    return _validated_entanglement_report(
        attempt.status,
        attempt.certified,
        attempt.certificate_kind,
        attempt.method,
        attempt.backend,
        attempt.raw_result,
        EntanglementAttempt[attempt],
        attempt.message,
    )
end

function detect_entanglement(rho::AbstractMatrix{<:Number}, dims, method::NativePPT)
    layout = _as_layout(dims)
    result = ppt_criterion(
        rho,
        layout;
        systems=method.systems,
        atol=method.atol,
        rtol=method.rtol,
        allow_densify=method.allow_densify,
    )
    attempt = if result.status === CriterionEntanglementDetected
        _validated_entanglement_attempt(
            :ppt,
            :native,
            :entangled,
            true,
            :negative_partial_transpose_witness,
            result,
            "a negative partial-transpose eigenpair certifies entanglement",
        )
    elseif result.status === CriterionSatisfied &&
        _ppt_low_dimension_separability_domain(layout)
        _validated_entanglement_attempt(
            :ppt,
            :native,
            :separable,
            true,
            :ppt_low_dimension_theorem,
            result,
            "PPT is sufficient for separability in 2×2 and 2×3 bipartite systems",
        )
    else
        _validated_entanglement_attempt(
            :ppt,
            :native,
            :unknown,
            false,
            nothing,
            result,
            if result.status === CriterionSatisfied
                "PPT is only a necessary condition in these dimensions"
            else
                "the PPT eigenvalue lies within the configured numerical boundary"
            end,
        )
    end
    return _report_from_attempt(attempt)
end

function _criterion_attempt(method::Symbol, result::CriterionResult)
    if result.status === CriterionEntanglementDetected
        certificate = if method === :realignment
            :realignment_cross_norm_violation
        elseif method === :reduction
            :reduction_map_witness
        else
            Symbol(method, :_violation)
        end
        return _validated_entanglement_attempt(
            method, :native, :entangled, true, certificate, result, result.message
        )
    end
    return _validated_entanglement_attempt(
        method,
        :native,
        :unknown,
        false,
        nothing,
        result,
        if result.status === CriterionSatisfied
            "$method is a necessary condition here and gives no separability certificate"
        else
            result.message
        end,
    )
end

"""
    analyze_entanglement(state, dims; strategy=:certificates_first, ...)

Run the dependency-free certificate pipeline. For density matrices the
`:certificates_first` strategy tries PPT, realignment/CCNR, and the reduction
criterion in that order, retaining every attempt. It stops on a valid
entanglement certificate or on the exact PPT-separability theorem in `2×2` and
`2×3`; otherwise it returns `:unknown`.

For bipartite pure vectors, a Schmidt coefficient above the configured
numerical threshold beyond the leading coefficient certifies entanglement.
Separability is reported only when the computed trailing coefficients are
exactly zero; a tolerance-defined rank-one result with nonzero trailing
coefficients remains `:unknown`. Inputs are never normalized or repaired.
"""
function analyze_entanglement(
    rho::AbstractMatrix{<:Number},
    dims;
    strategy::Symbol=:certificates_first,
    systems=(2,),
    atol=0.0,
    rtol=sqrt(eps(Float64)),
    allow_densify::Bool=false,
)
    strategy === :certificates_first || throw(
        ArgumentError(
            "the dependency-free core currently supports only " *
            "strategy=:certificates_first",
        ),
    )
    ppt_report = detect_entanglement(
        rho,
        dims,
        NativePPT(; systems=systems, atol=atol, rtol=rtol, allow_densify=allow_densify),
    )
    ppt_report.certified && return ppt_report
    attempts = collect(ppt_report.attempts)

    realignment_result = realignment_criterion(
        rho, dims; systems=(1,), atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    realignment_attempt = _criterion_attempt(:realignment, realignment_result)
    push!(attempts, realignment_attempt)
    if realignment_attempt.certified
        return _validated_entanglement_report(
            :entangled,
            true,
            realignment_attempt.certificate_kind,
            :realignment,
            :native,
            realignment_result,
            attempts,
            realignment_attempt.message,
        )
    end

    layout = _as_layout(dims)
    if length(layout) == 2
        reduction_result = reduction_criterion(
            rho, layout; side=:both, atol=atol, rtol=rtol, allow_densify=allow_densify
        )
        reduction_attempt = _criterion_attempt(:reduction, reduction_result)
        push!(attempts, reduction_attempt)
        if reduction_attempt.certified
            return _validated_entanglement_report(
                :entangled,
                true,
                reduction_attempt.certificate_kind,
                :reduction,
                :native,
                reduction_result,
                attempts,
                reduction_attempt.message,
            )
        end
    end

    return _validated_entanglement_report(
        :unknown,
        false,
        nothing,
        :certificates_first,
        :native,
        nothing,
        attempts,
        "the configured native criteria produced no valid entanglement or separability certificate",
    )
end

function analyze_entanglement(
    psi::AbstractVector{<:Number},
    dims;
    strategy::Symbol=:certificates_first,
    atol=0.0,
    rtol=sqrt(eps(Float64)),
    allow_densify::Bool=false,
)
    strategy === :certificates_first || throw(
        ArgumentError(
            "the dependency-free core currently supports only " *
            "strategy=:certificates_first",
        ),
    )
    layout = _tierd_bipartite_layout(dims, length(psi))
    validation = _tierd_validate_pure_state(
        psi; atol=atol, rtol=rtol, operation="analyze_entanglement"
    )
    decomposition = schmidt_decomposition(psi, layout; allow_densify=allow_densify)
    scale = maximum(
        decomposition.coefficients; init=zero(eltype(decomposition.coefficients))
    )
    threshold = validation.atol + validation.rtol * scale
    rank = count(value -> value > threshold, decomposition.coefficients)
    trailing_coefficients = Iterators.drop(decomposition.coefficients, 1)
    exact_product = rank == 1 && all(iszero, trailing_coefficients)
    status = rank > 1 ? :entangled : (exact_product ? :separable : :unknown)
    certified = status !== :unknown
    certificate = if status === :entangled
        :pure_state_schmidt_rank
    elseif status === :separable
        :pure_product_decomposition
    else
        nothing
    end
    message = if status === :entangled
        "Schmidt rank $rank certifies pure-state entanglement"
    elseif status === :separable
        "exactly zero trailing Schmidt coefficients certify a bipartite pure product state"
    else
        "the tolerance-defined Schmidt rank is $rank, but nonzero trailing coefficients or an oversized tolerance prevent a separability certificate"
    end
    raw_result = (
        decomposition=decomposition,
        rank=rank,
        norm_squared=validation.norm_squared,
        threshold=threshold,
    )
    attempt = _validated_entanglement_attempt(
        :pure_schmidt, :native, status, certified, certificate, raw_result, message
    )
    return _report_from_attempt(attempt)
end
