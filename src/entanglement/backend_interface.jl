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
undocumented external backend object.
"""
struct EntanglementAttempt
    method::Symbol
    backend::Symbol
    status::Symbol
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    raw_result::Any
    message::String

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
        certified &&
            certificate_kind === nothing &&
            throw(ArgumentError("a certified attempt must name its certificate_kind"))
        return new(
            method,
            backend,
            status,
            certified,
            certificate_kind,
            raw_result,
            String(message),
        )
    end
end

"""
    EntanglementReport

Structured high-level conclusion. `status` is `:entangled`, `:separable`, or
`:unknown`; `certified` is never true without a named `certificate_kind`.
`attempts` records every method run in order. A necessary criterion that merely
passes is represented as `:unknown`, except where a separately stated theorem
makes it sufficient.
"""
struct EntanglementReport
    status::Symbol
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    method::Symbol
    backend::Symbol
    evidence::Any
    attempts::Vector{EntanglementAttempt}
    message::String

    function EntanglementReport(
        status::Symbol,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        method::Symbol,
        backend::Symbol,
        evidence,
        attempts::Vector{EntanglementAttempt},
        message::AbstractString,
    )
        status in (:entangled, :separable, :unknown) || throw(
            ArgumentError("report status must be :entangled, :separable, or :unknown")
        )
        certified &&
            certificate_kind === nothing &&
            throw(ArgumentError("a certified report must name its certificate_kind"))
        status === :unknown &&
            certified &&
            throw(ArgumentError("an :unknown report cannot be certified"))
        return new(
            status,
            certified,
            certificate_kind,
            method,
            backend,
            evidence,
            attempts,
            String(message),
        )
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

"""
    available_entanglement_backends()

Return the backends currently registered in the dependency-free core. Optional
extensions may add explicit backend descriptors but are never loaded merely by
calling this function.
"""
available_entanglement_backends() = (NativeEntanglementBackend(),)

"""
    detect_entanglement(state, dims, method)

Run one explicitly selected entanglement method and return an
[`EntanglementReport`](@ref). Optional packages extend this generic function by
dispatch; absence or failure of a backend must not be mapped to a mathematical
negative result.
"""
function detect_entanglement end

function _ppt_low_dimension_separability_domain(layout::SubsystemLayout)
    length(layout) == 2 || return false
    first_dimension, second_dimension = layout.dims
    return (first_dimension == 2 && second_dimension in (2, 3)) ||
           (second_dimension == 2 && first_dimension in (2, 3))
end

function _report_from_attempt(attempt::EntanglementAttempt)
    return EntanglementReport(
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
        EntanglementAttempt(
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
        EntanglementAttempt(
            :ppt,
            :native,
            :separable,
            true,
            :ppt_low_dimension_theorem,
            result,
            "PPT is sufficient for separability in 2×2 and 2×3 bipartite systems",
        )
    else
        EntanglementAttempt(
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
        return EntanglementAttempt(
            method, :native, :entangled, true, certificate, result, result.message
        )
    end
    return EntanglementAttempt(
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
    attempts = copy(ppt_report.attempts)

    realignment_result = realignment_criterion(
        rho, dims; systems=(1,), atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    realignment_attempt = _criterion_attempt(:realignment, realignment_result)
    push!(attempts, realignment_attempt)
    if realignment_attempt.certified
        return EntanglementReport(
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
            return EntanglementReport(
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

    return EntanglementReport(
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
    attempt = EntanglementAttempt(
        :pure_schmidt, :native, status, certified, certificate, raw_result, message
    )
    return _report_from_attempt(attempt)
end
