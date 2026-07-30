# Source-informed independent Julia foundation for the Bell/nonlocal-game
# contracts in QETLAB BellInequalityMax.m, NonlocalGameLB.m,
# XORGameValue.m, BCSGameLB.m, BCSGameValue.m, and helpers/{bcs_to_nonlocal,
# cg2fc,cg2fp,fc2cg,fc2fp,fp2cg,fp2fc}.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
#
# QETLAB: Copyright 2014-2022 Nathaniel Johnston, Mateus Araújo, Vincent
# Russo, BSD-2-Clause. Full terms: licenses/QETLAB-LICENSE.txt.
#
# The types below are package-owned and solver-neutral. They make the axis
# order explicit and never normalize, clip, symmetrize, or repair caller data.

struct _NonlocalReadOnlyArray{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    storage::A
end

function _nonlocal_read_only(array::AbstractArray{T,N}) where {T,N}
    owned = copy(array)
    return _NonlocalReadOnlyArray{T,N,typeof(owned)}(owned)
end

Base.size(array::_NonlocalReadOnlyArray) = size(getfield(array, :storage))
Base.axes(array::_NonlocalReadOnlyArray) = axes(getfield(array, :storage))
Base.IndexStyle(::Type{_NonlocalReadOnlyArray{T,N,A}}) where {T,N,A} = Base.IndexStyle(A)
Base.copy(array::_NonlocalReadOnlyArray) = copy(getfield(array, :storage))
function Base.similar(array::_NonlocalReadOnlyArray, ::Type{T}, dimensions::Dims) where {T}
    return similar(getfield(array, :storage), T, dimensions)
end
function Base.similar(array::_NonlocalReadOnlyArray, ::Type{T}) where {T}
    return similar(getfield(array, :storage), T)
end
function SparseArrays.issparse(array::_NonlocalReadOnlyArray)
    return SparseArrays.issparse(getfield(array, :storage))
end

Base.@propagate_inbounds function Base.getindex(array::_NonlocalReadOnlyArray, indices...)
    return getfield(array, :storage)[indices...]
end

function Base.getproperty(array::_NonlocalReadOnlyArray, name::Symbol)
    name === :storage && return copy(getfield(array, :storage))
    return getfield(array, name)
end

function _nonlocal_read_only_nested(array::AbstractArray)
    owned_elements = map(array) do element
        return element isa AbstractArray ? _nonlocal_read_only(element) : element
    end
    return _nonlocal_read_only(owned_elements)
end

function _nonlocal_owned_evidence(value)
    value isa _NonlocalReadOnlyArray && return value
    value isa AbstractArray && return _nonlocal_read_only(value)
    return value
end

"""
    BellScenario(alice_outputs, bob_outputs, alice_settings, bob_settings)

Validated bipartite Bell scenario. Full-probability arrays always use axes
`(alice output, bob output, alice setting, bob setting)`.
"""
struct BellScenario
    alice_outputs::Int
    bob_outputs::Int
    alice_settings::Int
    bob_settings::Int

    function BellScenario(
        alice_outputs::Int, bob_outputs::Int, alice_settings::Int, bob_settings::Int
    )
        for (value, name) in (
            (alice_outputs, "alice_outputs"),
            (bob_outputs, "bob_outputs"),
            (alice_settings, "alice_settings"),
            (bob_settings, "bob_settings"),
        )
            value > 0 || throw(ArgumentError("$name must be a positive integer"))
        end
        return new(alice_outputs, bob_outputs, alice_settings, bob_settings)
    end
end

function _nonlocal_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _nonlocal_nonnegative_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _nonlocal_optional_limit(value, name::AbstractString; allow_zero::Bool=false)
    value === nothing && return nothing
    return BigInt(
        if allow_zero
            _nonlocal_nonnegative_integer(value, name)
        else
            _nonlocal_positive_integer(value, name)
        end,
    )
end

function _nonlocal_check_budget(required::Integer, limit, name::AbstractString)
    checked = _nonlocal_optional_limit(limit, name)
    checked === nothing && return BigInt(required)
    BigInt(required) <= checked ||
        throw(ArgumentError("$required required units exceed $name=$(Int(checked))"))
    return BigInt(required)
end

function BellScenario(alice_outputs, bob_outputs, alice_settings, bob_settings)
    return BellScenario(
        _nonlocal_positive_integer(alice_outputs, "alice_outputs"),
        _nonlocal_positive_integer(bob_outputs, "bob_outputs"),
        _nonlocal_positive_integer(alice_settings, "alice_settings"),
        _nonlocal_positive_integer(bob_settings, "bob_settings"),
    )
end

function BellScenario(desc::AbstractVector)
    Base.require_one_based_indexing(desc)
    length(desc) == 4 ||
        throw(DimensionMismatch("a Bell scenario description must have four entries"))
    return BellScenario(desc[1], desc[2], desc[3], desc[4])
end

function Base.Tuple(scenario::BellScenario)
    return (
        scenario.alice_outputs,
        scenario.bob_outputs,
        scenario.alice_settings,
        scenario.bob_settings,
    )
end

function _nonlocal_check_numeric_array(array, name::AbstractString)
    array isa AbstractArray{<:Number} ||
        throw(ArgumentError("$name must be a numeric array"))
    Base.require_one_based_indexing(array)
    for value in array
        value isa Bool && throw(ArgumentError("$name must not contain Boolean values"))
        isfinite(value) || throw(ArgumentError("$name must contain only finite values"))
    end
    return nothing
end

function _nonlocal_check_real_array(array, name::AbstractString)
    array isa AbstractArray{<:Real} || throw(ArgumentError("$name must be a real array"))
    _nonlocal_check_numeric_array(array, name)
    return nothing
end

function _nonlocal_tolerance(::Type{T}, atol, rtol, scale=one(T)) where {T<:Real}
    default = T <: AbstractFloat ? sqrt(eps(T)) : zero(T)
    absolute = atol === nothing ? default : atol
    relative = rtol === nothing ? default : rtol
    for (value, name) in ((absolute, "atol"), (relative, "rtol"))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    return convert(
        promote_type(T, typeof(absolute), typeof(relative)),
        absolute + relative * max(one(scale), abs(scale)),
    )
end

"""
    FullProbabilityBehavior(probabilities, scenario; atol=nothing, rtol=nothing)

Owned full behavior `P[a,b,x,y]`. The constructor checks finite real entries,
per-setting normalization, and no-signalling. Values inside a floating
tolerance band are retained exactly as supplied; no clipping occurs.
"""
struct FullProbabilityBehavior{T<:Real,A<:AbstractArray{T,4},R<:Real}
    scenario::BellScenario
    probabilities::A
    tolerance::R
    minimum_probability::T
    normalization_residual::R
    no_signalling_residual::R
    boundary::Bool

    function FullProbabilityBehavior(
        token::_ValidatedConstructorToken,
        scenario::BellScenario,
        probabilities::A,
        tolerance::R,
        minimum_probability::T,
        normalization_residual::R,
        no_signalling_residual::R,
        boundary::Bool,
    ) where {T<:Real,A<:AbstractArray{T,4},R<:Real}
        _require_validated_constructor_token(token)
        probabilities isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("behavior probabilities require read-only storage"))
        size(probabilities) == Tuple(scenario) ||
            throw(DimensionMismatch("behavior probabilities and scenario do not match"))
        tolerance >= zero(tolerance) && isfinite(tolerance) ||
            throw(ArgumentError("behavior tolerance must be finite and nonnegative"))
        minimum_probability == minimum(probabilities) ||
            throw(ArgumentError("behavior minimum-probability diagnostic is inconsistent"))
        normalization_residual >= zero(normalization_residual) &&
        no_signalling_residual >= zero(no_signalling_residual) ||
            throw(ArgumentError("behavior residuals must be nonnegative"))
        return new{T,A,R}(
            scenario,
            probabilities,
            tolerance,
            minimum_probability,
            normalization_residual,
            no_signalling_residual,
            boundary,
        )
    end
end

function _nonlocal_behavior_diagnostics(probabilities, scenario::BellScenario; atol, rtol)
    T = eltype(probabilities)
    scale = isempty(probabilities) ? one(T) : maximum(abs, probabilities)
    tolerance = _nonlocal_tolerance(T, atol, rtol, scale)
    minimum_probability = minimum(probabilities)
    minimum_probability < -tolerance && throw(
        DomainError(
            minimum_probability,
            "probabilities contain a value below -tolerance; no clipping is applied",
        ),
    )
    normalization_residual = zero(tolerance)
    for x in 1:scenario.alice_settings, y in 1:scenario.bob_settings
        total = sum(
            probabilities[a, b, x, y] for
            a in 1:scenario.alice_outputs, b in 1:scenario.bob_outputs
        )
        normalization_residual = max(normalization_residual, abs(total - one(total)))
    end
    normalization_residual <= tolerance || throw(
        ArgumentError(
            "each setting pair must sum to one; maximum residual is " *
            "$normalization_residual and tolerance is $tolerance. Inputs are never normalized.",
        ),
    )
    no_signalling_residual = zero(tolerance)
    for x in 1:scenario.alice_settings, a in 1:scenario.alice_outputs
        reference = sum(probabilities[a, b, x, 1] for b in 1:scenario.bob_outputs)
        for y in 2:scenario.bob_settings
            value = sum(probabilities[a, b, x, y] for b in 1:scenario.bob_outputs)
            no_signalling_residual = max(no_signalling_residual, abs(value - reference))
        end
    end
    for y in 1:scenario.bob_settings, b in 1:scenario.bob_outputs
        reference = sum(probabilities[a, b, 1, y] for a in 1:scenario.alice_outputs)
        for x in 2:scenario.alice_settings
            value = sum(probabilities[a, b, x, y] for a in 1:scenario.alice_outputs)
            no_signalling_residual = max(no_signalling_residual, abs(value - reference))
        end
    end
    no_signalling_residual <= tolerance || throw(
        ArgumentError(
            "the behavior is signalling; maximum marginal residual is " *
            "$no_signalling_residual and tolerance is $tolerance",
        ),
    )
    boundary =
        minimum_probability < zero(minimum_probability) ||
        normalization_residual > zero(normalization_residual) ||
        no_signalling_residual > zero(no_signalling_residual)
    return (
        tolerance=tolerance,
        minimum_probability=minimum_probability,
        normalization_residual=normalization_residual,
        no_signalling_residual=no_signalling_residual,
        boundary=boundary,
    )
end

function FullProbabilityBehavior(
    probabilities::AbstractArray{<:Real,4},
    scenario::BellScenario;
    atol=nothing,
    rtol=nothing,
)
    _nonlocal_check_real_array(probabilities, "probabilities")
    expected = Tuple(scenario)
    size(probabilities) == expected || throw(
        DimensionMismatch(
            "full behavior has size $(size(probabilities)); expected $expected"
        ),
    )
    diagnostics = _nonlocal_behavior_diagnostics(
        probabilities, scenario; atol=atol, rtol=rtol
    )
    owned = _nonlocal_read_only(probabilities)
    return FullProbabilityBehavior(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        scenario,
        owned,
        diagnostics.tolerance,
        diagnostics.minimum_probability,
        diagnostics.normalization_residual,
        diagnostics.no_signalling_residual,
        diagnostics.boundary,
    )
end

"""
    CollinsGisinBehavior(coefficients, scenario; atol=nothing, rtol=nothing)

Owned Collins--Gisin behavior. Row one/column one contain normalization and
marginals; the remaining blocks contain joint probabilities for all but the
last outcome of each party.
"""
struct CollinsGisinBehavior{T<:Real,M<:AbstractMatrix{T},R<:Real}
    scenario::BellScenario
    coefficients::M
    tolerance::R

    function CollinsGisinBehavior(
        token::_ValidatedConstructorToken,
        scenario::BellScenario,
        coefficients::M,
        tolerance::R,
    ) where {T<:Real,M<:AbstractMatrix{T},R<:Real}
        _require_validated_constructor_token(token)
        coefficients isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("Collins-Gisin behavior requires read-only storage"))
        size(coefficients) == _nonlocal_cg_shape(scenario) ||
            throw(DimensionMismatch("Collins-Gisin behavior and scenario do not match"))
        tolerance >= zero(tolerance) && isfinite(tolerance) ||
            throw(ArgumentError("behavior tolerance must be finite and nonnegative"))
        return new{T,M,R}(scenario, coefficients, tolerance)
    end
end

function _nonlocal_cg_shape(scenario::BellScenario)
    return (
        1 + (scenario.alice_outputs - 1) * scenario.alice_settings,
        1 + (scenario.bob_outputs - 1) * scenario.bob_settings,
    )
end

function _nonlocal_aindex(a::Int, x::Int, scenario::BellScenario)
    return 1 + a + (x - 1) * (scenario.alice_outputs - 1)
end
function _nonlocal_bindex(b::Int, y::Int, scenario::BellScenario)
    return 1 + b + (y - 1) * (scenario.bob_outputs - 1)
end

function CollinsGisinBehavior(
    coefficients::AbstractMatrix{<:Real}, scenario::BellScenario; atol=nothing, rtol=nothing
)
    _nonlocal_check_real_array(coefficients, "Collins-Gisin behavior")
    size(coefficients) == _nonlocal_cg_shape(scenario) || throw(
        DimensionMismatch(
            "Collins-Gisin behavior has size $(size(coefficients)); expected " *
            "$(_nonlocal_cg_shape(scenario))",
        ),
    )
    tolerance = _nonlocal_tolerance(
        eltype(coefficients), atol, rtol, maximum(abs, coefficients)
    )
    abs(coefficients[1, 1] - one(eltype(coefficients))) <= tolerance || throw(
        ArgumentError(
            "Collins-Gisin normalization entry must equal one within tolerance; " *
            "it is never repaired",
        ),
    )
    owned = _nonlocal_read_only(coefficients)
    # Conversion followed by diagnostics validates positivity,
    # normalization, and no-signalling while retaining the exact supplied
    # coefficients.
    full = _nonlocal_cg_to_full_raw(owned, scenario)
    _nonlocal_behavior_diagnostics(full, scenario; atol=tolerance, rtol=zero(tolerance))
    return CollinsGisinBehavior(_VALIDATED_CONSTRUCTOR_TOKEN, scenario, owned, tolerance)
end

function _nonlocal_cg_to_full_raw(cg, scenario::BellScenario)
    oa, ob, ma, mb = Tuple(scenario)
    T = eltype(cg)
    full = zeros(T, oa, ob, ma, mb)
    for x in 1:ma, y in 1:mb
        for a in 1:(oa - 1), b in 1:(ob - 1)
            full[a, b, x, y] = cg[
                _nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)
            ]
        end
        for a in 1:(oa - 1)
            full[a, ob, x, y] =
                cg[_nonlocal_aindex(a, x, scenario), 1] - sum(
                    cg[_nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)]
                    for b in 1:(ob - 1)
                )
        end
        for b in 1:(ob - 1)
            full[oa, b, x, y] =
                cg[1, _nonlocal_bindex(b, y, scenario)] - sum(
                    cg[_nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)]
                    for a in 1:(oa - 1)
                )
        end
        full[oa, ob, x, y] =
            cg[1, 1] - sum(cg[_nonlocal_aindex(a, x, scenario), 1] for a in 1:(oa - 1)) -
            sum(cg[1, _nonlocal_bindex(b, y, scenario)] for b in 1:(ob - 1)) + sum(
                cg[_nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)] for
                a in 1:(oa - 1), b in 1:(ob - 1)
            )
    end
    return full
end

function full_probability_behavior(
    behavior::CollinsGisinBehavior; atol=nothing, rtol=nothing
)
    return FullProbabilityBehavior(
        _nonlocal_cg_to_full_raw(behavior.coefficients, behavior.scenario),
        behavior.scenario;
        atol=isnothing(atol) ? behavior.tolerance : atol,
        rtol=isnothing(rtol) ? zero(behavior.tolerance) : rtol,
    )
end

function collins_gisin_behavior(
    behavior::FullProbabilityBehavior; atol=nothing, rtol=nothing
)
    scenario = behavior.scenario
    oa, ob, ma, mb = Tuple(scenario)
    T = eltype(behavior.probabilities)
    cg = zeros(T, _nonlocal_cg_shape(scenario))
    cg[1, 1] = one(T)
    for x in 1:ma, a in 1:(oa - 1)
        cg[_nonlocal_aindex(a, x, scenario), 1] = sum(
            behavior.probabilities[a, b, x, 1] for b in 1:ob
        )
    end
    for y in 1:mb, b in 1:(ob - 1)
        cg[1, _nonlocal_bindex(b, y, scenario)] = sum(
            behavior.probabilities[a, b, 1, y] for a in 1:oa
        )
    end
    for x in 1:ma, y in 1:mb
        for a in 1:(oa - 1), b in 1:(ob - 1)
            cg[_nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)] = behavior.probabilities[
                a, b, x, y
            ]
        end
    end
    return CollinsGisinBehavior(
        cg,
        scenario;
        atol=isnothing(atol) ? behavior.tolerance : atol,
        rtol=isnothing(rtol) ? zero(behavior.tolerance) : rtol,
    )
end

"""
    BellFunctional(coefficients, scenario; notation=:full_probability)

Owned Bell functional. Internally it stores the full-probability coefficients
`C[a,b,x,y]`; `source_notation` records the caller's notation. Collins--Gisin
and full-correlator inputs are converted algebraically without normalizing or
repairing coefficients.
"""
struct BellFunctional{T<:Real,A<:AbstractArray{T,4},S}
    scenario::BellScenario
    coefficients::A
    source_notation::Symbol
    source_coefficients::S

    function BellFunctional(
        token::_ValidatedConstructorToken,
        scenario::BellScenario,
        coefficients::A,
        source_notation::Symbol,
        source_coefficients::S,
    ) where {T<:Real,A<:AbstractArray{T,4},S}
        _require_validated_constructor_token(token)
        coefficients isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("Bell coefficients require read-only storage"))
        source_coefficients isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("source Bell coefficients require read-only storage"))
        size(coefficients) == Tuple(scenario) ||
            throw(DimensionMismatch("Bell coefficients and scenario do not match"))
        source_notation in
        (:full_probability, :fp, :collins_gisin, :cg, :full_correlator, :fc) ||
            throw(ArgumentError("invalid Bell source notation"))
        return new{T,A,S}(scenario, coefficients, source_notation, source_coefficients)
    end
end

_nonlocal_fraction_type(::Type{T}) where {T<:Integer} = Rational{BigInt}
_nonlocal_fraction_type(::Type{<:Rational}) = Rational{BigInt}
_nonlocal_fraction_type(::Type{T}) where {T<:AbstractFloat} = T
_nonlocal_fraction_type(::Type{T}) where {T<:Real} = promote_type(T, Float64)

function _nonlocal_cg_functional_to_full(cg, scenario::BellScenario)
    oa, ob, ma, mb = Tuple(scenario)
    T = _nonlocal_fraction_type(eltype(cg))
    converted = T.(cg)
    full = zeros(T, oa, ob, ma, mb)
    for x in 1:ma, y in 1:mb
        for a in 1:(oa - 1), b in 1:(ob - 1)
            full[a, b, x, y] =
                converted[1, 1] / (ma * mb) +
                converted[_nonlocal_aindex(a, x, scenario), 1] / mb +
                converted[1, _nonlocal_bindex(b, y, scenario)] / ma +
                converted[
                    _nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)
                ]
        end
        for a in 1:(oa - 1)
            full[a, ob, x, y] =
                converted[1, 1] / (ma * mb) +
                converted[_nonlocal_aindex(a, x, scenario), 1] / mb
        end
        for b in 1:(ob - 1)
            full[oa, b, x, y] =
                converted[1, 1] / (ma * mb) +
                converted[1, _nonlocal_bindex(b, y, scenario)] / ma
        end
        full[oa, ob, x, y] = converted[1, 1] / (ma * mb)
    end
    return full
end

function _nonlocal_fc_functional_to_full(fc, scenario::BellScenario)
    scenario.alice_outputs == 2 && scenario.bob_outputs == 2 || throw(
        ArgumentError("full-correlator notation requires exactly two outcomes per party"),
    )
    ma, mb = scenario.alice_settings, scenario.bob_settings
    size(fc) == (ma + 1, mb + 1) || throw(
        DimensionMismatch(
            "full-correlator matrix has size $(size(fc)); expected $((ma + 1, mb + 1))"
        ),
    )
    T = _nonlocal_fraction_type(eltype(fc))
    converted = T.(fc)
    full = zeros(T, 2, 2, ma, mb)
    for x in 1:ma, y in 1:mb
        full[1, 1, x, y] =
            converted[1, 1] / (ma * mb) +
            converted[x + 1, 1] / mb +
            converted[1, y + 1] / ma +
            converted[x + 1, y + 1]
        full[1, 2, x, y] =
            converted[1, 1] / (ma * mb) + converted[x + 1, 1] / mb -
            converted[1, y + 1] / ma - converted[x + 1, y + 1]
        full[2, 1, x, y] =
            converted[1, 1] / (ma * mb) - converted[x + 1, 1] / mb +
            converted[1, y + 1] / ma - converted[x + 1, y + 1]
        full[2, 2, x, y] =
            converted[1, 1] / (ma * mb) - converted[x + 1, 1] / mb -
            converted[1, y + 1] / ma + converted[x + 1, y + 1]
    end
    return full
end

function BellFunctional(
    coefficients,
    scenario::BellScenario;
    notation::Symbol=:full_probability,
    max_entries=1_000_000,
)
    notation in (:full_probability, :fp, :collins_gisin, :cg, :full_correlator, :fc) ||
        throw(
            ArgumentError(
                "notation must be :full_probability/:fp, :collins_gisin/:cg, " *
                "or :full_correlator/:fc",
            ),
        )
    required = prod(BigInt, Tuple(scenario))
    _nonlocal_check_budget(required, max_entries, "max_entries")
    full = if notation in (:full_probability, :fp)
        coefficients isa AbstractArray{<:Real,4} || throw(
            ArgumentError("full-probability coefficients must be a real 4-D array")
        )
        _nonlocal_check_real_array(coefficients, "Bell coefficients")
        size(coefficients) == Tuple(scenario) || throw(
            DimensionMismatch(
                "Bell coefficients have size $(size(coefficients)); expected " *
                "$(Tuple(scenario))",
            ),
        )
        copy(coefficients)
    elseif notation in (:collins_gisin, :cg)
        coefficients isa AbstractMatrix{<:Real} ||
            throw(ArgumentError("Collins-Gisin coefficients must be a real matrix"))
        _nonlocal_check_real_array(coefficients, "Bell coefficients")
        size(coefficients) == _nonlocal_cg_shape(scenario) || throw(
            DimensionMismatch(
                "Collins-Gisin coefficients have size $(size(coefficients)); expected " *
                "$(_nonlocal_cg_shape(scenario))",
            ),
        )
        _nonlocal_cg_functional_to_full(coefficients, scenario)
    else
        coefficients isa AbstractMatrix{<:Real} ||
            throw(ArgumentError("full-correlator coefficients must be a real matrix"))
        _nonlocal_check_real_array(coefficients, "Bell coefficients")
        _nonlocal_fc_functional_to_full(coefficients, scenario)
    end
    return BellFunctional(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        scenario,
        _nonlocal_read_only(full),
        notation,
        _nonlocal_read_only(coefficients),
    )
end

function _nonlocal_full_functional_to_cg(functional::BellFunctional)
    C = functional.coefficients
    scenario = functional.scenario
    oa, ob, ma, mb = Tuple(scenario)
    T = eltype(C)
    cg = zeros(T, _nonlocal_cg_shape(scenario))
    cg[1, 1] = sum(C[oa, ob, x, y] for x in 1:ma, y in 1:mb)
    for a in 1:(oa - 1), x in 1:ma
        cg[_nonlocal_aindex(a, x, scenario), 1] = sum(
            C[a, ob, x, y] - C[oa, ob, x, y] for y in 1:mb
        )
    end
    for b in 1:(ob - 1), y in 1:mb
        cg[1, _nonlocal_bindex(b, y, scenario)] = sum(
            C[oa, b, x, y] - C[oa, ob, x, y] for x in 1:ma
        )
    end
    for a in 1:(oa - 1), b in 1:(ob - 1), x in 1:ma, y in 1:mb
        cg[_nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)] =
            C[a, b, x, y] - C[a, ob, x, y] - C[oa, b, x, y] + C[oa, ob, x, y]
    end
    return cg
end

function evaluate_bell_functional(
    functional::BellFunctional, behavior::FullProbabilityBehavior
)
    functional.scenario == behavior.scenario ||
        throw(DimensionMismatch("functional and behavior scenarios do not match"))
    return sum(functional.coefficients .* behavior.probabilities)
end

function evaluate_bell_functional(
    functional::BellFunctional, behavior::CollinsGisinBehavior
)
    functional.scenario == behavior.scenario ||
        throw(DimensionMismatch("functional and behavior scenarios do not match"))
    return sum(_nonlocal_full_functional_to_cg(functional) .* behavior.coefficients)
end

"""
    NonlocalGame(probabilities, payoff; atol=nothing, rtol=nothing)

Owned nonlocal game. `probabilities[x,y]` is a normalized question
distribution and `payoff[a,b,x,y]` is a finite real reward. Rewards are not
restricted to `[0,1]` and are never rescaled.
"""
struct NonlocalGame{T<:Real,P<:AbstractMatrix{T},R<:Real,V<:AbstractArray{R,4},Q}
    scenario::BellScenario
    probabilities::P
    payoff::V
    tolerance::Q

    function NonlocalGame(
        token::_ValidatedConstructorToken,
        scenario::BellScenario,
        probabilities::P,
        payoff::V,
        tolerance::Q,
    ) where {T<:Real,P<:AbstractMatrix{T},R<:Real,V<:AbstractArray{R,4},Q}
        _require_validated_constructor_token(token)
        probabilities isa _NonlocalReadOnlyArray && payoff isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("nonlocal-game data require read-only storage"))
        size(probabilities) == (scenario.alice_settings, scenario.bob_settings) ||
            throw(DimensionMismatch("question probabilities and scenario do not match"))
        size(payoff) == Tuple(scenario) ||
            throw(DimensionMismatch("payoff and scenario do not match"))
        tolerance isa Real &&
        !(tolerance isa Bool) &&
        isfinite(tolerance) &&
        tolerance >= zero(tolerance) ||
            throw(ArgumentError("game tolerance must be finite and nonnegative"))
        return new{T,P,R,V,Q}(scenario, probabilities, payoff, tolerance)
    end
end

function NonlocalGame(
    probabilities::AbstractMatrix{<:Real},
    payoff::AbstractArray{<:Real,4};
    atol=nothing,
    rtol=nothing,
    max_entries=1_000_000,
)
    _nonlocal_check_real_array(probabilities, "question probabilities")
    _nonlocal_check_real_array(payoff, "payoff")
    ma, mb = size(probabilities)
    oa, ob, payoff_ma, payoff_mb = size(payoff)
    (payoff_ma, payoff_mb) == (ma, mb) || throw(
        DimensionMismatch(
            "payoff setting dimensions $((payoff_ma, payoff_mb)) do not match " *
            "probability size $((ma, mb))",
        ),
    )
    scenario = BellScenario(oa, ob, ma, mb)
    _nonlocal_check_budget(length(payoff), max_entries, "max_entries")
    T = promote_type(eltype(probabilities), eltype(payoff))
    tolerance = _nonlocal_tolerance(T, atol, rtol, max(maximum(abs, probabilities), one(T)))
    minimum_probability = minimum(probabilities)
    minimum_probability < -tolerance && throw(
        DomainError(minimum_probability, "question probabilities must be nonnegative")
    )
    total = sum(probabilities)
    abs(total - one(total)) <= tolerance || throw(
        ArgumentError(
            "question probabilities must sum to one; their sum is $total. " *
            "They are never normalized implicitly.",
        ),
    )
    return NonlocalGame(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        scenario,
        _nonlocal_read_only(probabilities),
        _nonlocal_read_only(payoff),
        tolerance,
    )
end

function BellFunctional(game::NonlocalGame; max_entries=1_000_000)
    _nonlocal_check_budget(length(game.payoff), max_entries, "max_entries")
    coefficients = similar(
        game.payoff, promote_type(eltype(game.probabilities), eltype(game.payoff))
    )
    for index in CartesianIndices(game.payoff)
        a, b, x, y = Tuple(index)
        coefficients[index] = game.probabilities[x, y] * game.payoff[a, b, x, y]
    end
    return BellFunctional(
        coefficients, game.scenario; notation=:full_probability, max_entries=max_entries
    )
end

"""
    BCSGame(constraints; max_entries=...)

Validated binary constraint-system game. Every constraint is an owned `0/1`
array over the same binary variable axes.
"""
struct BCSGame{C}
    constraints::C
    variable_count::Int

    function BCSGame(
        token::_ValidatedConstructorToken, constraints::C, variable_count::Int
    ) where {C}
        _require_validated_constructor_token(token)
        constraints isa Tuple && !isempty(constraints) ||
            throw(ArgumentError("BCS constraints must be a nonempty immutable collection"))
        variable_count > 0 || throw(ArgumentError("BCS variable count must be positive"))
        all(constraint -> constraint isa _NonlocalReadOnlyArray, constraints) ||
            throw(ArgumentError("BCS constraints require read-only storage"))
        expected_shape = ntuple(_ -> 2, variable_count)
        all(constraint -> size(constraint) == expected_shape, constraints) ||
            throw(DimensionMismatch("BCS constraint shapes are inconsistent"))
        return new{C}(constraints, variable_count)
    end
end

function BCSGame(constraints::AbstractVector; max_entries=1_000_000)
    Base.require_one_based_indexing(constraints)
    isempty(constraints) && throw(ArgumentError("at least one BCS constraint is required"))
    first_constraint = first(constraints)
    first_constraint isa AbstractArray ||
        throw(ArgumentError("each BCS constraint must be an array"))
    variable_count = ndims(first_constraint)
    variable_count > 0 || throw(ArgumentError("BCS constraints need at least one variable"))
    expected_shape = ntuple(_ -> 2, variable_count)
    total_entries = BigInt(length(constraints)) * prod(BigInt, expected_shape)
    _nonlocal_check_budget(total_entries, max_entries, "max_entries")
    owned = Any[]
    for (index, constraint) in enumerate(constraints)
        constraint isa AbstractArray ||
            throw(ArgumentError("constraint $index must be an array"))
        Base.require_one_based_indexing(constraint)
        size(constraint) == expected_shape || throw(
            DimensionMismatch(
                "constraint $index has size $(size(constraint)); expected $expected_shape",
            ),
        )
        all(value -> value in (false, true, 0, 1), constraint) || throw(
            ArgumentError("constraint $index must contain only binary zero/one entries")
        )
        push!(owned, _nonlocal_read_only(Int.(constraint)))
    end
    return BCSGame(_VALIDATED_CONSTRUCTOR_TOKEN, Tuple(owned), variable_count)
end

function _nonlocal_constraint_active(constraint, variable::Int)
    axes_except = ntuple(
        index -> index == variable ? (1:1) : axes(constraint, index), ndims(constraint)
    )
    for index in CartesianIndices(axes_except)
        first_index = Tuple(index)
        second_index = ntuple(
            dimension -> dimension == variable ? 2 : first_index[dimension],
            ndims(constraint),
        )
        constraint[first_index...] == constraint[second_index...] || return true
    end
    return false
end

"""
    nonlocal_game(game::BCSGame; max_entries=...)

Convert a BCS game into the explicit nonlocal-game representation used by the
pinned `bcs_to_nonlocal` helper. The stale upstream `update_odometer`
dependency is not reproduced; checked Julia Cartesian indexing supplies the
same finite enumeration.
"""
function nonlocal_game(game::BCSGame; max_entries=1_000_000)
    constraint_count = length(game.constraints)
    variable_count = game.variable_count
    active = falses(constraint_count, variable_count)
    max_active = 0
    for x in 1:constraint_count
        for y in 1:variable_count
            active[x, y] = _nonlocal_constraint_active(game.constraints[x], y)
        end
        count = sum(view(active, x, :))
        count > 0 || throw(
            ArgumentError(
                "constraint $x contains no active variable, so the pinned question " *
                "distribution would divide by zero",
            ),
        )
        max_active = max(max_active, count)
    end
    alice_outputs_big = big(2)^max_active
    alice_outputs_big <= typemax(Int) ||
        throw(ArgumentError("the BCS Alice output count is too large for Int"))
    alice_outputs = Int(alice_outputs_big)
    required = BigInt(alice_outputs) * 2 * constraint_count * variable_count
    _nonlocal_check_budget(required, max_entries, "max_entries")
    probabilities = zeros(Rational{BigInt}, constraint_count, variable_count)
    for x in 1:constraint_count
        count = sum(view(active, x, :))
        for y in 1:variable_count
            active[x, y] || continue
            probabilities[x, y] = Rational{BigInt}(1, constraint_count * count)
        end
    end
    payoff = zeros(Int, alice_outputs, 2, constraint_count, variable_count)
    assignment = ones(Int, variable_count)
    for encoded in 0:(alice_outputs - 1)
        # QETLAB's `dec2bin` enumerates the first active variable as the
        # most-significant bit. Julia's `digits` is least-significant first,
        # so reverse it to preserve the public Alice-output labels.
        bits = reverse(digits(encoded; base=2, pad=max_active))
        for x in 1:constraint_count
            fill!(assignment, 1)
            cursor = 1
            for variable in 1:variable_count
                if active[x, variable]
                    assignment[variable] = bits[cursor] + 1
                    cursor += 1
                end
            end
            satisfies = game.constraints[x][assignment...] == 1
            satisfies || continue
            for y in 1:variable_count
                # QETLAB fills every payoff question, including variables
                # whose question probability is zero for this constraint.
                # Inactive variables retain the helper's default answer one.
                payoff[encoded + 1, assignment[y], x, y] = 1
            end
        end
    end
    return NonlocalGame(probabilities, payoff; atol=0, rtol=0, max_entries=max_entries)
end

function _nonlocal_decode_strategy(index::Integer, base::Int, digits_count::Int)
    strategy = Vector{Int}(undef, digits_count)
    value = BigInt(index)
    for position in digits_count:-1:1
        strategy[position] = Int(mod(value, base)) + 1
        value = div(value, base)
    end
    return strategy
end

function _nonlocal_classical_optimum(
    game::NonlocalGame; max_strategies=1_000_000, max_work=100_000_000
)
    oa, ob, ma, mb = Tuple(game.scenario)
    alice_count = big(oa)^ma
    bob_count = big(ob)^mb
    enumerate_bob = bob_count <= alice_count
    strategy_count = enumerate_bob ? bob_count : alice_count
    _nonlocal_check_budget(strategy_count, max_strategies, "max_strategies")
    work = if enumerate_bob
        strategy_count * ma * oa * mb
    else
        strategy_count * mb * ob * ma
    end
    _nonlocal_check_budget(work, max_work, "max_work")
    T = promote_type(eltype(game.probabilities), eltype(game.payoff))
    best = nothing
    best_alice = Int[]
    best_bob = Int[]
    if enumerate_bob
        for encoded in BigInt(0):(strategy_count - 1)
            bob = _nonlocal_decode_strategy(encoded, ob, mb)
            alice = Vector{Int}(undef, ma)
            total = zero(T)
            for x in 1:ma
                best_x = nothing
                best_a = 1
                for a in 1:oa
                    value = sum(
                        game.probabilities[x, y] * game.payoff[a, bob[y], x, y] for
                        y in 1:mb
                    )
                    if best_x === nothing || value > best_x
                        best_x = value
                        best_a = a
                    end
                end
                alice[x] = best_a
                total += best_x
            end
            if best === nothing || total > best
                best = total
                best_alice = alice
                best_bob = bob
            end
        end
    else
        for encoded in BigInt(0):(strategy_count - 1)
            alice = _nonlocal_decode_strategy(encoded, oa, ma)
            bob = Vector{Int}(undef, mb)
            total = zero(T)
            for y in 1:mb
                best_y = nothing
                best_b = 1
                for b in 1:ob
                    value = sum(
                        game.probabilities[x, y] * game.payoff[alice[x], b, x, y] for
                        x in 1:ma
                    )
                    if best_y === nothing || value > best_y
                        best_y = value
                        best_b = b
                    end
                end
                bob[y] = best_b
                total += best_y
            end
            if best === nothing || total > best
                best = total
                best_alice = alice
                best_bob = bob
            end
        end
    end
    return (
        value=best,
        alice_strategy=best_alice,
        bob_strategy=best_bob,
        enumerated_party=enumerate_bob ? :bob : :alice,
        strategies_evaluated=strategy_count,
        joint_strategy_count=alice_count * bob_count,
        scalar_work=work,
    )
end
