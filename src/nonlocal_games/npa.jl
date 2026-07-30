# Solver-neutral NPA and Bell/XOR optimization models.
#
# QETLAB source contract: NPAHierarchy.m, BellInequalityMax.m, and
# XORGameValue.m at d8589610f00cff106537268dee2e2a1153f3a601.
# Recorded source SHA-256 values:
#   NPAHierarchy.m       ee0e2bdba84ed1ad2ea7ff7060b9a5883e66cd3eeecc9cfae7773a72854b27b9
#   BellInequalityMax.m  fac680c2990af158c5b44a774079de8105aab89ab8f057233aae55ced61a8df0
#   XORGameValue.m       dbbd2b5088559760c140f4c1ef5034f5e3aed2d0699f93608dc73db94e97d84c
#
# QETLAB source is BSD-2-Clause (licenses/QETLAB-LICENSE.txt).
#
# The hierarchy formulation is independently specified from:
# M. Navascués, S. Pironio, and A. Acín, "A convergent hierarchy of
# semidefinite programs characterizing the set of quantum correlations",
# New J. Phys. 10, 073013 (2008), https://arxiv.org/abs/0803.4290.
#
# The XOR SDP follows:
# R. Cleve, P. Høyer, B. Toner, and J. Watrous, "Consequences and limits of
# nonlocal strategies", https://arxiv.org/abs/quant-ph/0404076.

"""
    NPALetter(party, setting, outcome)

One non-last-outcome projector in a bipartite NPA word. `party` is `:alice`
or `:bob`; settings and outcomes are one-based.
"""
struct NPALetter
    party::UInt8
    setting::Int
    outcome::Int

    function NPALetter(
        token::_ValidatedConstructorToken, party::UInt8, setting::Int, outcome::Int
    )
        _require_validated_constructor_token(token)
        party in (0x01, 0x02) || throw(ArgumentError("invalid internal NPA party code"))
        setting > 0 || throw(ArgumentError("NPA settings must be positive"))
        outcome > 0 || throw(ArgumentError("NPA outcomes must be positive"))
        return new(party, setting, outcome)
    end
end

function NPALetter(party::Symbol, setting, outcome)
    code = if party === :alice
        UInt8(1)
    elseif party === :bob
        UInt8(2)
    else
        throw(ArgumentError("party must be :alice or :bob"))
    end
    return NPALetter(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        code,
        _nonlocal_positive_integer(setting, "setting"),
        _nonlocal_positive_integer(outcome, "outcome"),
    )
end

_npa_party(letter::NPALetter) = letter.party == 0x01 ? :alice : :bob

function Base.:(==)(left::NPALetter, right::NPALetter)
    return left.party == right.party &&
           left.setting == right.setting &&
           left.outcome == right.outcome
end
function Base.hash(letter::NPALetter, seed::UInt)
    return hash((letter.party, letter.setting, letter.outcome), seed)
end
function Base.isless(left::NPALetter, right::NPALetter)
    return (left.party, left.setting, left.outcome) <
           (right.party, right.setting, right.outcome)
end

"""
    NPAWord(letters)

Immutable ordered product of NPA projector letters. Use [`npa_word`](@ref) to
apply commuting-party, idempotence, and orthogonality reductions.
"""
struct NPAWord
    letters::Tuple{Vararg{NPALetter}}
end

NPAWord() = NPAWord(())
NPAWord(letters::AbstractVector{NPALetter}) = NPAWord(Tuple(letters))

Base.length(word::NPAWord) = length(word.letters)
Base.isempty(word::NPAWord) = isempty(word.letters)
Base.iterate(word::NPAWord, state...) = iterate(word.letters, state...)
Base.:(==)(left::NPAWord, right::NPAWord) = left.letters == right.letters
Base.hash(word::NPAWord, seed::UInt) = hash(word.letters, seed)

function _npa_word_isless(left::NPAWord, right::NPAWord)
    n = min(length(left), length(right))
    for index in 1:n
        left.letters[index] == right.letters[index] && continue
        return isless(left.letters[index], right.letters[index])
    end
    return length(left) < length(right)
end

function _npa_reduce_party(letters::Vector{NPALetter})
    reduced = NPALetter[]
    for letter in letters
        if !isempty(reduced) && last(reduced).setting == letter.setting
            last(reduced).outcome == letter.outcome || return nothing
            # Same projector twice: P^2 = P.
            continue
        end
        push!(reduced, letter)
    end
    return reduced
end

"""
    npa_word(letters)

Canonicalize a product. Alice and Bob letters commute, equal adjacent
projectors are idempotent, and different outcomes of one measurement are
orthogonal. `nothing` denotes the zero operator.
"""
function npa_word(letters)
    collected = NPALetter[]
    for letter in letters
        letter isa NPALetter ||
            throw(ArgumentError("NPA words may contain only NPALetter values"))
        push!(collected, letter)
    end
    alice = _npa_reduce_party([letter for letter in collected if letter.party == 0x01])
    alice === nothing && return nothing
    bob = _npa_reduce_party([letter for letter in collected if letter.party == 0x02])
    bob === nothing && return nothing
    return NPAWord(Tuple(vcat(alice, bob)))
end

function _npa_adjoint(word::NPAWord)
    return npa_word(reverse(word.letters))
end

function _npa_product(left::NPAWord, right::NPAWord)
    return npa_word((left.letters..., right.letters...))
end

function _npa_moment_word(left::NPAWord, right::NPAWord)
    adjoint_left = _npa_adjoint(left)
    adjoint_left === nothing && return nothing
    return npa_word((adjoint_left.letters..., right.letters...))
end

struct NPALevel
    label::String
    base_level::Int
    maximum_length::Int
    components::Tuple{Vararg{Tuple{Int,Int}}}

    function NPALevel(
        token::_ValidatedConstructorToken,
        label::AbstractString,
        base_level::Int,
        maximum_length::Int,
        components::Tuple{Vararg{Tuple{Int,Int}}},
    )
        _require_validated_constructor_token(token)
        base_level >= 0 || throw(ArgumentError("base NPA level must be nonnegative"))
        maximum_length >= base_level ||
            throw(ArgumentError("maximum NPA word length is inconsistent"))
        all(
            component ->
                component[1] >= 0 && component[2] >= 0 && sum(component) <= maximum_length,
            components,
        ) || throw(ArgumentError("NPA level components are inconsistent"))
        return new(String(label), base_level, maximum_length, components)
    end
end

function NPALevel(level::Integer)
    checked = _nonlocal_nonnegative_integer(level, "level")
    return NPALevel(_VALIDATED_CONSTRUCTOR_TOKEN, string(checked), checked, checked, ())
end

function NPALevel(level::AbstractString)
    tokens = split(lowercase(strip(level)), '+')
    isempty(tokens) && throw(ArgumentError("NPA level string must not be empty"))
    base = tryparse(Int, strip(first(tokens)))
    base === nothing &&
        throw(ArgumentError("the first NPA level component must be a nonnegative integer"))
    base = _nonlocal_nonnegative_integer(base, "base NPA level")
    components = Tuple{Int,Int}[]
    maximum_length = base
    for token in Iterators.drop(tokens, 1)
        cleaned = strip(token)
        isempty(cleaned) && throw(ArgumentError("NPA level components must not be empty"))
        all(character -> character in ('a', 'b'), cleaned) ||
            throw(ArgumentError("intermediate NPA components may contain only 'a' and 'b'"))
        count_a = count(==('a'), cleaned)
        count_b = count(==('b'), cleaned)
        push!(components, (count_a, count_b))
        maximum_length = max(maximum_length, length(cleaned))
    end
    return NPALevel(
        _VALIDATED_CONSTRUCTOR_TOKEN, strip(level), base, maximum_length, Tuple(components)
    )
end

NPALevel(level::NPALevel) = level

function _npa_word_allowed(word::NPAWord, level::NPALevel)
    length(word) <= level.base_level && return true
    count_a = count(letter -> letter.party == 0x01, word.letters)
    count_b = length(word) - count_a
    return any(
        component ->
            length(word) <= sum(component) &&
            count_a <= component[1] &&
            count_b <= component[2],
        level.components,
    )
end

function _npa_alphabet(scenario::BellScenario)
    alphabet = NPALetter[]
    for x in 1:scenario.alice_settings, a in 1:(scenario.alice_outputs - 1)
        push!(alphabet, NPALetter(:alice, x, a))
    end
    for y in 1:scenario.bob_settings, b in 1:(scenario.bob_outputs - 1)
        push!(alphabet, NPALetter(:bob, y, b))
    end
    return alphabet
end

function _npa_catalog(
    scenario::BellScenario,
    level::NPALevel;
    max_words=512,
    max_word_generation_work=2_000_000,
)
    word_limit = _nonlocal_positive_integer(max_words, "max_words")
    work_limit = _nonlocal_positive_integer(
        max_word_generation_work, "max_word_generation_work"
    )
    alphabet = _npa_alphabet(scenario)
    alphabet_count = length(alphabet)
    estimated_work = sum(big(alphabet_count)^degree for degree in 0:level.maximum_length)
    estimated_work <= work_limit || throw(
        ArgumentError(
            "NPA word generation needs up to $estimated_work raw products, " *
            "exceeding max_word_generation_work=$work_limit",
        ),
    )
    catalog = NPAWord[NPAWord()]
    seen = Set(catalog)
    frontier = NPAWord[NPAWord()]
    work = 0
    for _ in 1:level.maximum_length
        next_frontier = NPAWord[]
        for prefix in frontier, letter in alphabet
            work += 1
            word = npa_word((prefix.letters..., letter))
            word === nothing && continue
            push!(next_frontier, word)
            _npa_word_allowed(word, level) || continue
            if !(word in seen)
                length(catalog) < word_limit || throw(
                    ArgumentError(
                        "NPA catalog exceeds max_words=$word_limit before allocation"
                    ),
                )
                push!(seen, word)
                push!(catalog, word)
            end
        end
        frontier = unique(next_frontier)
    end
    return Tuple(catalog), work, estimated_work
end

"""
    NPAStatus

Certificate-aware NPA outcome. Basic probability/no-signalling violations can
be exact certificates. Floating SDP feasibility or infeasibility remains
numerical evidence unless independently certified.
"""
@enum NPAStatus::UInt8 begin
    NPABasicConditionsSatisfied
    NPABasicConditionsViolated
    NPANumericallyFeasible
    NPANumericallyInfeasible
    NPANumericalBoundary
    NPAResourceLimit
    NPABackendUnavailable
    NPABackendFailure
    NPAInvalidCertificate
end

"""
    NPAProblem

Package-owned NPA feasibility or Bell-objective SDP. `words` records the
canonical moment basis and `program` contains only numeric affine data.
"""
struct NPAProblem{L,W,P,B,O,M}
    scenario::BellScenario
    level::L
    words::W
    program::P
    behavior::B
    objective_kind::Symbol
    objective_metadata::O
    moment_view::M
    word_generation_work::Int

    function NPAProblem(
        token::_ValidatedConstructorToken,
        scenario::BellScenario,
        level::L,
        words::W,
        program::P,
        behavior::B,
        objective_kind::Symbol,
        objective_metadata::O,
        moment_view::M,
        word_generation_work::Int,
    ) where {L,W,P,B,O,M}
        _require_validated_constructor_token(token)
        level isa NPALevel || throw(ArgumentError("NPA problem level is invalid"))
        words isa Tuple || throw(ArgumentError("NPA problem words must be immutable"))
        !isempty(words) && first(words) == NPAWord() ||
            throw(ArgumentError("an NPA word catalog must start with the identity"))
        all(word -> word isa NPAWord, words) ||
            throw(ArgumentError("an NPA catalog may contain only NPA words"))
        length(Set(words)) == length(words) ||
            throw(ArgumentError("an NPA word catalog must not contain duplicates"))
        objective_kind in (:feasibility, :bell_upper_bound) ||
            throw(ArgumentError("invalid NPA objective kind"))
        word_generation_work >= 0 ||
            throw(ArgumentError("NPA word-generation work must be nonnegative"))
        behavior === nothing ||
            behavior.scenario == scenario ||
            throw(DimensionMismatch("NPA behavior and scenario do not match"))
        return new{L,W,P,B,O,M}(
            scenario,
            level,
            words,
            program,
            behavior,
            objective_kind,
            objective_metadata,
            moment_view,
            word_generation_work,
        )
    end
end

"""
    NPAResult

Status-rich hierarchy result. `verdict` is set only for exact basic-condition
branches. A numerical feasible point supplies `moment_matrix` and residuals
but is not promoted to an exact quantum-membership statement.
"""
struct NPAResult{P,M,O,R,T}
    status::NPAStatus
    verdict::Union{Nothing,Bool}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    level::NPALevel
    problem::P
    moment_matrix::M
    optimization_result::O
    residuals::R
    tolerance::T
    message::String

    function NPAResult(
        token::_ValidatedConstructorToken,
        status::NPAStatus,
        verdict::Union{Nothing,Bool},
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        level::NPALevel,
        problem,
        moment_matrix,
        optimization_result,
        residuals,
        tolerance,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        if certified
            status === NPABasicConditionsSatisfied ||
                throw(ArgumentError("only exact basic conditions can be certified"))
            verdict === true ||
                throw(ArgumentError("a certified basic-condition result must be true"))
            certificate_kind === :exact_probability_and_no_signalling_conditions ||
                throw(ArgumentError("invalid NPA certificate kind"))
        else
            verdict === nothing ||
                throw(ArgumentError("an uncertified NPA result cannot carry a verdict"))
            certificate_kind === nothing ||
                throw(ArgumentError("an uncertified NPA result cannot name a certificate"))
        end
        status === NPABasicConditionsSatisfied &&
            !certified &&
            throw(ArgumentError("exact basic-condition success must be certified"))
        tolerance isa Real &&
        !(tolerance isa Bool) &&
        isfinite(tolerance) &&
        tolerance >= zero(tolerance) ||
            throw(ArgumentError("NPA tolerance must be finite and nonnegative"))
        owned_moment = _nonlocal_owned_evidence(moment_matrix)
        return new{
            typeof(problem),
            typeof(owned_moment),
            typeof(optimization_result),
            typeof(residuals),
            typeof(tolerance),
        }(
            status,
            verdict,
            certified,
            certificate_kind,
            level,
            problem,
            owned_moment,
            optimization_result,
            residuals,
            tolerance,
            String(message),
        )
    end
end

function Base.show(io::IO, result::NPAResult)
    return print(
        io,
        "NPAResult(status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", level=",
        repr(result.level.label),
        ")",
    )
end

function _npa_entry_affine(
    matrix::HermitianAffineMatrix{T}, row::Int, column::Int, component::Symbol
) where {T<:Real}
    component in (:real, :imag) || throw(ArgumentError("component must be :real or :imag"))
    entry = matrix.constant[row, column]
    constant = component === :real ? real(entry) : imag(entry)
    variables = Int[]
    coefficients = T[]
    for term in matrix.terms
        value = term.coefficient[row, column]
        coefficient = component === :real ? real(value) : imag(value)
        iszero(coefficient) && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    return AffineScalar(constant, variables, coefficients, matrix.variable_count)
end

function _npa_affine_difference(
    left::AffineScalar{T}, right::AffineScalar{T}; right_sign=one(T)
) where {T<:Real}
    coefficients = copy(left.coefficients)
    coefficients -= right_sign * right.coefficients
    return AffineScalar(
        left.constant - right_sign * right.constant,
        coefficients;
        variable_count=length(coefficients),
    )
end

function _npa_affine_shift(function_data::AffineScalar{T}, target) where {T<:Real}
    return AffineScalar(
        function_data.constant - convert(T, target),
        function_data.coefficients;
        variable_count=length(function_data.coefficients),
    )
end

function _npa_known_value(
    word::NPAWord, behavior::Union{Nothing,CollinsGisinBehavior}, scenario::BellScenario
)
    isempty(word) &&
        return one(behavior === nothing ? Float64 : eltype(behavior.coefficients))
    behavior === nothing && return nothing
    if length(word) == 1
        letter = first(word.letters)
        return if letter.party == 0x01
            behavior.coefficients[
                _nonlocal_aindex(letter.outcome, letter.setting, scenario), 1
            ]
        else
            behavior.coefficients[
                1, _nonlocal_bindex(letter.outcome, letter.setting, scenario)
            ]
        end
    elseif length(word) == 2
        first_letter, second_letter = word.letters
        if first_letter.party == 0x01 && second_letter.party == 0x02
            return behavior.coefficients[
                _nonlocal_aindex(first_letter.outcome, first_letter.setting, scenario),
                _nonlocal_bindex(second_letter.outcome, second_letter.setting, scenario),
            ]
        end
    end
    return nothing
end

function _npa_word_orientation(word::NPAWord)
    reversed = _npa_adjoint(word)
    reversed === nothing && return word, false
    if _npa_word_isless(reversed, word)
        return reversed, true
    end
    return word, false
end

function _npa_moment_constraints(
    moment::HermitianAffineMatrix{T},
    words,
    scenario::BellScenario,
    behavior::Union{Nothing,CollinsGisinBehavior},
) where {T<:Real}
    equalities = AffineEquality{T}[]
    representatives = Dict{NPAWord,Tuple{Int,Int,Bool}}()
    for row in eachindex(words), column in row:length(words)
        word = _npa_moment_word(words[row], words[column])
        real_entry = _npa_entry_affine(moment, row, column, :real)
        imag_entry = _npa_entry_affine(moment, row, column, :imag)
        if word === nothing
            push!(
                equalities,
                AffineEquality(real_entry, Symbol(:npa_zero_real_, row, :_, column)),
            )
            push!(
                equalities,
                AffineEquality(imag_entry, Symbol(:npa_zero_imag_, row, :_, column)),
            )
            continue
        end
        known = _npa_known_value(word, behavior, scenario)
        if known !== nothing
            push!(
                equalities,
                AffineEquality(
                    _npa_affine_shift(real_entry, known),
                    Symbol(:npa_known_real_, row, :_, column),
                ),
            )
            push!(
                equalities,
                AffineEquality(imag_entry, Symbol(:npa_known_imag_, row, :_, column)),
            )
            continue
        end
        representative, conjugated = _npa_word_orientation(word)
        if representative == _npa_adjoint(representative)
            push!(
                equalities,
                AffineEquality(
                    imag_entry, Symbol(:npa_self_adjoint_imag_, row, :_, column)
                ),
            )
        end
        if !haskey(representatives, representative)
            representatives[representative] = (row, column, conjugated)
            continue
        end
        ref_row, ref_column, ref_conjugated = representatives[representative]
        ref_real = _npa_entry_affine(moment, ref_row, ref_column, :real)
        ref_imag = _npa_entry_affine(moment, ref_row, ref_column, :imag)
        push!(
            equalities,
            AffineEquality(
                _npa_affine_difference(real_entry, ref_real),
                Symbol(:npa_word_real_, row, :_, column),
            ),
        )
        sign = conjugated == ref_conjugated ? one(T) : -one(T)
        push!(
            equalities,
            AffineEquality(
                _npa_affine_difference(imag_entry, ref_imag; right_sign=sign),
                Symbol(:npa_word_imag_, row, :_, column),
            ),
        )
    end
    return equalities
end

function _npa_objective_from_cg(
    moment::HermitianAffineMatrix{T}, words, scenario::BellScenario, cg_coefficients
) where {T<:Real}
    lookup = Dict{NPAWord,Int}(word => index for (index, word) in enumerate(words))
    objective = AffineScalar(
        convert(T, cg_coefficients[1, 1]),
        zeros(T, moment.variable_count);
        variable_count=moment.variable_count,
    )
    coefficients = copy(objective.coefficients)
    constant = objective.constant
    identity_index = lookup[NPAWord()]
    function add_entry!(word::NPAWord, coefficient)
        iszero(coefficient) && return nothing
        if haskey(lookup, word)
            entry = _npa_entry_affine(moment, identity_index, lookup[word], :real)
        elseif length(word) == 2
            first_word = NPAWord((word.letters[1],))
            second_word = NPAWord((word.letters[2],))
            haskey(lookup, first_word) && haskey(lookup, second_word) ||
                throw(ArgumentError("NPA level does not contain the required Bell moment"))
            entry = _npa_entry_affine(
                moment, lookup[first_word], lookup[second_word], :real
            )
        else
            throw(ArgumentError("NPA level does not contain the required Bell moment"))
        end
        constant += convert(T, coefficient) * entry.constant
        coefficients .+= convert(T, coefficient) .* entry.coefficients
        return nothing
    end
    for x in 1:scenario.alice_settings, a in 1:(scenario.alice_outputs - 1)
        add_entry!(
            NPAWord((NPALetter(:alice, x, a),)),
            cg_coefficients[_nonlocal_aindex(a, x, scenario), 1],
        )
    end
    for y in 1:scenario.bob_settings, b in 1:(scenario.bob_outputs - 1)
        add_entry!(
            NPAWord((NPALetter(:bob, y, b),)),
            cg_coefficients[1, _nonlocal_bindex(b, y, scenario)],
        )
    end
    for x in 1:scenario.alice_settings, y in 1:scenario.bob_settings
        for a in 1:(scenario.alice_outputs - 1), b in 1:(scenario.bob_outputs - 1)
            add_entry!(
                NPAWord((NPALetter(:alice, x, a), NPALetter(:bob, y, b))),
                cg_coefficients[
                    _nonlocal_aindex(a, x, scenario), _nonlocal_bindex(b, y, scenario)
                ],
            )
        end
    end
    return AffineScalar(constant, coefficients; variable_count=moment.variable_count)
end

function _npa_build_problem(
    scenario::BellScenario,
    level::NPALevel;
    behavior::Union{Nothing,CollinsGisinBehavior}=nothing,
    cg_objective=nothing,
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    level.maximum_length > 0 ||
        throw(ArgumentError("an SDP NPA problem requires level at least one"))
    words, work, _ = _npa_catalog(
        scenario,
        level;
        max_words=max_words,
        max_word_generation_work=max_word_generation_work,
    )
    dimension = length(words)
    variable_count_big = big(dimension)^2
    variable_count_big <= typemax(Int) ||
        throw(ArgumentError("the NPA moment variable count is too large for Int"))
    variable_count_big <= limits.max_variables || throw(
        ArgumentError(
            "the NPA moment model needs $variable_count_big variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    dimension <= limits.max_psd_dimension || throw(
        ArgumentError(
            "the NPA moment matrix has dimension $dimension, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    variable_count = Int(variable_count_big)
    coefficient_type = if behavior !== nothing
        eltype(behavior.coefficients)
    elseif cg_objective !== nothing
        eltype(cg_objective)
    else
        Float64
    end
    coefficient_type <: Real ||
        throw(ArgumentError("NPA coefficients must have a real element type"))
    moment = hermitian_variable(
        :npa_moment,
        dimension;
        variable_count=variable_count,
        coefficient_type=coefficient_type,
    )
    equalities = _npa_moment_constraints(moment, words, scenario, behavior)
    objective = if cg_objective === nothing
        AffineScalar(
            zero(coefficient_type),
            zeros(coefficient_type, variable_count);
            variable_count=variable_count,
        )
    else
        _npa_objective_from_cg(moment, words, scenario, cg_objective)
    end
    feasible = zeros(coefficient_type, variable_count)
    # Identity moment. This is only a warm start, not a known feasible point
    # for arbitrary supplied behavior.
    feasible[1] = one(coefficient_type)
    program = SemidefiniteProgram(
        cg_objective === nothing ? :npa_membership : :npa_bell_upper_bound,
        cg_objective === nothing ? :feasibility : :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=[moment],
        primal_views=[moment],
        initial_point=feasible,
        known_feasible_point=nothing,
        limits=limits,
        metadata=(
            formulation=:navascues_pironio_acin_moment_relaxation,
            level=level.label,
            scenario=Tuple(scenario),
            word_count=dimension,
            word_generation_work=work,
            real_moment_variables=variable_count,
        ),
    )
    return NPAProblem(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        scenario,
        level,
        words,
        program,
        behavior,
        cg_objective === nothing ? :feasibility : :bell_upper_bound,
        if cg_objective === nothing
            NamedTuple()
        else
            (cg_coefficients=_nonlocal_read_only(cg_objective),)
        end,
        moment,
        work,
    )
end

"""
    npa_problem(behavior; level=1, ...)

Build a solver-neutral NPA membership SDP from a validated full-probability or
Collins--Gisin behavior.
"""
function npa_problem(
    behavior::CollinsGisinBehavior;
    level=1,
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    parsed = NPALevel(level)
    return _npa_build_problem(
        behavior.scenario,
        parsed;
        behavior=behavior,
        max_words=max_words,
        max_word_generation_work=max_word_generation_work,
        limits=limits,
    )
end

function npa_problem(behavior::FullProbabilityBehavior; kwargs...)
    return npa_problem(collins_gisin_behavior(behavior); kwargs...)
end

function npa_problem(
    correlations::AbstractMatrix{<:Real},
    scenario::BellScenario;
    atol=nothing,
    rtol=nothing,
    kwargs...,
)
    return npa_problem(
        CollinsGisinBehavior(correlations, scenario; atol=atol, rtol=rtol); kwargs...
    )
end

function npa_problem(correlations::AbstractMatrix{<:Real}, desc::AbstractVector; kwargs...)
    return npa_problem(correlations, BellScenario(desc); kwargs...)
end

function _npa_status_from_optimization(optimization)
    if optimization.status in (OptimizationOptimal, OptimizationFeasible)
        return NPANumericallyFeasible
    elseif optimization.status === OptimizationInfeasible
        return NPANumericallyInfeasible
    elseif optimization.status === OptimizationBackendUnavailable
        return NPABackendUnavailable
    elseif optimization.status === OptimizationLimit
        return NPAResourceLimit
    else
        return NPABackendFailure
    end
end

function _npa_result_from_problem(problem::NPAProblem, backend)
    optimization = solve_optimization(problem.program, backend)
    status = _npa_status_from_optimization(optimization)
    moment = if optimization.primal === nothing
        nothing
    elseif haskey(optimization.primal.views, :npa_moment)
        Matrix(optimization.primal.views.npa_moment)
    else
        nothing
    end
    tolerance = if problem.behavior === nothing
        backend isa JuMPBackend ? backend.atol + backend.rtol : 0.0
    else
        problem.behavior.tolerance
    end
    residuals = (
        primal=optimization.primal_residual,
        dual=optimization.dual_residual,
        absolute_gap=optimization.absolute_gap,
        relative_gap=optimization.relative_gap,
    )
    message = if status === NPANumericallyFeasible
        "a numerical moment matrix satisfies the requested finite NPA relaxation; " *
        "this is not an exact quantum-membership certificate"
    elseif status === NPANumericallyInfeasible
        "the backend reported finite-level infeasibility; inspect its dual status " *
        "and certificate before using it as exclusion evidence"
    elseif status === NPABackendUnavailable
        "the package-owned NPA model was built, but no optimization backend was available"
    elseif status === NPAResourceLimit
        "the optimization backend reached a configured resource limit"
    else
        "the optimization backend did not return independently validated NPA evidence"
    end
    return NPAResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        nothing,
        false,
        nothing,
        problem.level,
        problem,
        moment,
        optimization,
        residuals,
        tolerance,
        message,
    )
end

"""
    npa_membership(behavior; level=1, backend=NoOptimizationBackend(), ...)

Check basic level-zero probability/no-signalling conditions exactly or build
and solve a finite NPA relaxation. Numerical solver output is never rounded to
a Boolean.
"""
function npa_membership(
    behavior::Union{FullProbabilityBehavior,CollinsGisinBehavior};
    level=1,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    parsed = NPALevel(level)
    full = if behavior isa FullProbabilityBehavior
        behavior
    else
        full_probability_behavior(behavior)
    end
    if parsed.maximum_length == 0
        if full.boundary
            return NPAResult(
                _VALIDATED_CONSTRUCTOR_TOKEN,
                NPANumericalBoundary,
                nothing,
                false,
                nothing,
                parsed,
                nothing,
                nothing,
                nothing,
                (
                    normalization=full.normalization_residual,
                    no_signalling=full.no_signalling_residual,
                    minimum_probability=full.minimum_probability,
                ),
                full.tolerance,
                "basic conditions hold only within a nonzero tolerance band",
            )
        end
        return NPAResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            NPABasicConditionsSatisfied,
            true,
            true,
            :exact_probability_and_no_signalling_conditions,
            parsed,
            nothing,
            nothing,
            nothing,
            (
                normalization=full.normalization_residual,
                no_signalling=full.no_signalling_residual,
                minimum_probability=full.minimum_probability,
            ),
            full.tolerance,
            "the full behavior exactly satisfies positivity, normalization, and " *
            "no-signalling",
        )
    end
    problem = npa_problem(
        behavior;
        level=parsed,
        max_words=max_words,
        max_word_generation_work=max_word_generation_work,
        limits=limits,
    )
    return _npa_result_from_problem(problem, backend)
end

function npa_membership(
    correlations::AbstractMatrix{<:Real},
    scenario::BellScenario;
    atol=nothing,
    rtol=nothing,
    kwargs...,
)
    behavior = CollinsGisinBehavior(correlations, scenario; atol=atol, rtol=rtol)
    return npa_membership(behavior; kwargs...)
end

function npa_membership(
    correlations::AbstractMatrix{<:Real}, desc::AbstractVector; kwargs...
)
    return npa_membership(correlations, BellScenario(desc); kwargs...)
end
