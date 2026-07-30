# Status-aware exact, no-signalling, XOR, and NPA Bell-value routes.
#
# This file composes the validated scenario layer and the solver-neutral NPA
# layer. No optimizer is selected implicitly and no numerical optimum is
# promoted to theorem-level exactness.

@enum NonlocalValueStatus::UInt8 begin
    NonlocalValueExact
    NonlocalValueNumericalLowerBound
    NonlocalValueNumericalUpperBound
    NonlocalValueNumericalInterval
    NonlocalValueResourceLimit
    NonlocalValueBackendUnavailable
    NonlocalValueBackendFailure
    NonlocalValueInvalidCertificate
end

"""
    DeterministicStrategy

Owned one-based deterministic outputs for every Alice and Bob setting.
"""
struct DeterministicStrategy{A<:AbstractVector{Int},B<:AbstractVector{Int}}
    alice_outputs::A
    bob_outputs::B

    function DeterministicStrategy(
        token::_ValidatedConstructorToken, alice_outputs::A, bob_outputs::B
    ) where {A<:AbstractVector{Int},B<:AbstractVector{Int}}
        _require_validated_constructor_token(token)
        alice_outputs isa _NonlocalReadOnlyArray &&
        bob_outputs isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("deterministic strategies require read-only storage"))
        !isempty(alice_outputs) && !isempty(bob_outputs) ||
            throw(ArgumentError("deterministic strategies must answer every setting"))
        all(>(0), alice_outputs) && all(>(0), bob_outputs) ||
            throw(ArgumentError("deterministic outputs must be positive"))
        return new{A,B}(alice_outputs, bob_outputs)
    end
end

function DeterministicStrategy(
    alice_outputs::AbstractVector{<:Integer}, bob_outputs::AbstractVector{<:Integer}
)
    Base.require_one_based_indexing(alice_outputs, bob_outputs)
    any(value -> value isa Bool, alice_outputs) &&
        throw(ArgumentError("Alice outputs must be non-Boolean integers"))
    any(value -> value isa Bool, bob_outputs) &&
        throw(ArgumentError("Bob outputs must be non-Boolean integers"))
    alice = Int.(alice_outputs)
    bob = Int.(bob_outputs)
    return DeterministicStrategy(
        _VALIDATED_CONSTRUCTOR_TOKEN, _nonlocal_read_only(alice), _nonlocal_read_only(bob)
    )
end

"""
    NonlocalValueResult

Status-rich value or bound. `value` is populated only for exact exhaustive
branches or when independently validated lower and upper bounds coincide.
Numerical relaxation values remain in `lower_bound` or `upper_bound` with
their certification flags set explicitly.
"""
struct NonlocalValueResult{T,L,U,S,P,O,W,D}
    status::NonlocalValueStatus
    regime::Symbol
    value::T
    lower_bound::L
    upper_bound::U
    exact::Bool
    certified_lower::Bool
    certified_upper::Bool
    certificate_kind::Union{Nothing,Symbol}
    strategy::S
    problem::P
    optimization_result::O
    witness::W
    diagnostics::D
    message::String

    function NonlocalValueResult(
        token::_ValidatedConstructorToken,
        status::NonlocalValueStatus,
        regime::Symbol,
        value,
        lower_bound,
        upper_bound,
        exact::Bool,
        certified_lower::Bool,
        certified_upper::Bool,
        certificate_kind::Union{Nothing,Symbol},
        strategy,
        problem,
        optimization_result,
        witness,
        diagnostics,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        if exact
            status === NonlocalValueExact ||
                throw(ArgumentError("an exact nonlocal value needs exact status"))
            value !== nothing && lower_bound == value && upper_bound == value ||
                throw(ArgumentError("exact value and bounds must agree"))
            certified_lower && certified_upper ||
                throw(ArgumentError("an exact value certifies both bounds"))
            certificate_kind !== nothing ||
                throw(ArgumentError("an exact value must name its certificate"))
            strategy isa DeterministicStrategy ||
                throw(ArgumentError("an exact value needs a deterministic strategy"))
        else
            value === nothing ||
                throw(ArgumentError("a nonexact result cannot populate value"))
            !certified_lower && !certified_upper || throw(
                ArgumentError("numerical bounds are not exact mathematical certificates"),
            )
            certificate_kind === nothing ||
                throw(ArgumentError("a numerical result cannot name an exact certificate"))
            strategy === nothing || throw(
                ArgumentError("a numerical value result cannot carry an exact strategy")
            )
        end
        status === NonlocalValueNumericalUpperBound &&
            upper_bound === nothing &&
            throw(ArgumentError("a numerical-upper status requires an upper bound"))
        status in (
                NonlocalValueResourceLimit,
                NonlocalValueBackendUnavailable,
                NonlocalValueBackendFailure,
                NonlocalValueInvalidCertificate,
            ) &&
            (lower_bound !== nothing || upper_bound !== nothing) &&
            throw(ArgumentError("failed value results cannot carry bounds"))
        owned_witness = _nonlocal_owned_evidence(witness)
        return new{
            typeof(value),
            typeof(lower_bound),
            typeof(upper_bound),
            typeof(strategy),
            typeof(problem),
            typeof(optimization_result),
            typeof(owned_witness),
            typeof(diagnostics),
        }(
            status,
            regime,
            value,
            lower_bound,
            upper_bound,
            exact,
            certified_lower,
            certified_upper,
            certificate_kind,
            strategy,
            problem,
            optimization_result,
            owned_witness,
            diagnostics,
            String(message),
        )
    end
end

function Base.show(io::IO, result::NonlocalValueResult)
    return print(
        io,
        "NonlocalValueResult(status=",
        result.status,
        ", regime=",
        result.regime,
        ", value=",
        result.value,
        ", lower_bound=",
        result.lower_bound,
        ", upper_bound=",
        result.upper_bound,
        ")",
    )
end

function _nonlocal_status_from_optimization(optimization, bound_kind::Symbol)
    if optimization.status in (OptimizationOptimal, OptimizationFeasible)
        return if bound_kind === :lower
            NonlocalValueNumericalLowerBound
        else
            NonlocalValueNumericalUpperBound
        end
    elseif optimization.status === OptimizationLimit
        return NonlocalValueResourceLimit
    elseif optimization.status === OptimizationBackendUnavailable
        return NonlocalValueBackendUnavailable
    else
        return NonlocalValueBackendFailure
    end
end

function _nonlocal_solver_upper(optimization)
    for candidate in (optimization.objective_bound, optimization.dual_objective_value)
        candidate === nothing && continue
        isfinite(candidate) || continue
        return candidate
    end
    return nothing
end

function _nonlocal_exact_result(
    value,
    regime::Symbol,
    strategy,
    diagnostics;
    certificate_kind=:exhaustive_deterministic_strategies,
    message="exhaustive deterministic-strategy optimization completed",
)
    return NonlocalValueResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        NonlocalValueExact,
        regime,
        value,
        value,
        value,
        true,
        true,
        true,
        certificate_kind,
        strategy,
        nothing,
        nothing,
        strategy,
        diagnostics,
        message,
    )
end

function _nonlocal_probability_index(a::Int, b::Int, x::Int, y::Int, scenario::BellScenario)
    oa, ob, ma, _ = Tuple(scenario)
    return a + oa * ((b - 1) + ob * ((x - 1) + ma * (y - 1)))
end

function _nonlocal_probability_affine(scenario::BellScenario, variable_count::Int, terms)
    variables = Int[]
    coefficients = Float64[]
    for (index, coefficient) in terms
        iszero(coefficient) && continue
        push!(variables, index)
        push!(coefficients, Float64(coefficient))
    end
    return AffineScalar(0.0, variables, coefficients, variable_count)
end

function _nonlocal_lp_coefficient_type(functional::BellFunctional)
    T = eltype(functional.coefficients)
    T in (Float32, Float64) || throw(
        ArgumentError(
            "the optional optimization route currently requires Float32 or Float64 " *
            "Bell coefficients; exact inputs remain supported by the classical route",
        ),
    )
    return T
end

function _nonlocal_no_signalling_problem(
    functional::BellFunctional; limits::OptimizationLimits=OptimizationLimits()
)
    scenario = functional.scenario
    oa, ob, ma, mb = Tuple(scenario)
    variable_count_big = prod(BigInt, Tuple(scenario))
    variable_count_big <= typemax(Int) ||
        throw(ArgumentError("the no-signalling variable count is too large for Int"))
    variable_count_big <= limits.max_variables || throw(
        ArgumentError(
            "the no-signalling model needs $variable_count_big variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    equality_count =
        BigInt(ma) * mb +
        BigInt(ma) * oa * max(mb - 1, 0) +
        BigInt(mb) * ob * max(ma - 1, 0)
    equality_count <= limits.max_equalities || throw(
        ArgumentError(
            "the no-signalling model needs $equality_count equalities, exceeding " *
            "max_equalities=$(limits.max_equalities)",
        ),
    )
    variable_count_big <= limits.max_intervals || throw(
        ArgumentError(
            "the no-signalling model needs $variable_count_big intervals, exceeding " *
            "max_intervals=$(limits.max_intervals)",
        ),
    )
    variable_count = Int(variable_count_big)
    T = _nonlocal_lp_coefficient_type(functional)
    objective_variables = collect(1:variable_count)
    objective_coefficients = Vector{T}(undef, variable_count)
    for y in 1:mb, x in 1:ma, b in 1:ob, a in 1:oa
        index = _nonlocal_probability_index(a, b, x, y, scenario)
        objective_coefficients[index] = functional.coefficients[a, b, x, y]
    end
    objective = AffineScalar(
        zero(T), objective_variables, objective_coefficients, variable_count
    )
    equalities = AffineEquality{T}[]
    for x in 1:ma, y in 1:mb
        variables = Int[
            _nonlocal_probability_index(a, b, x, y, scenario) for a in 1:oa for b in 1:ob
        ]
        coefficients = fill(one(T), length(variables))
        push!(
            equalities,
            AffineEquality(
                AffineScalar(-one(T), variables, coefficients, variable_count),
                Symbol(:normalization_, x, :_, y),
            ),
        )
    end
    for x in 1:ma, a in 1:oa, y in 2:mb
        variables = Int[]
        coefficients = T[]
        for b in 1:ob
            push!(variables, _nonlocal_probability_index(a, b, x, y, scenario))
            push!(coefficients, one(T))
            push!(variables, _nonlocal_probability_index(a, b, x, 1, scenario))
            push!(coefficients, -one(T))
        end
        push!(
            equalities,
            AffineEquality(
                AffineScalar(zero(T), variables, coefficients, variable_count),
                Symbol(:alice_no_signalling_, a, :_, x, :_, y),
            ),
        )
    end
    for y in 1:mb, b in 1:ob, x in 2:ma
        variables = Int[]
        coefficients = T[]
        for a in 1:oa
            push!(variables, _nonlocal_probability_index(a, b, x, y, scenario))
            push!(coefficients, one(T))
            push!(variables, _nonlocal_probability_index(a, b, 1, y, scenario))
            push!(coefficients, -one(T))
        end
        push!(
            equalities,
            AffineEquality(
                AffineScalar(zero(T), variables, coefficients, variable_count),
                Symbol(:bob_no_signalling_, b, :_, x, :_, y),
            ),
        )
    end
    intervals = AffineInterval{T}[]
    for variable in 1:variable_count
        push!(
            intervals,
            AffineInterval(
                AffineScalar(zero(T), [variable], [one(T)], variable_count),
                zero(T),
                nothing,
                Symbol(:probability_nonnegative_, variable),
            ),
        )
    end
    feasible = fill(one(T) / convert(T, oa * ob), variable_count)
    program = SemidefiniteProgram(
        :no_signalling_bell_bound,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        intervals=intervals,
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:full_probability_no_signalling_linear_program,
            scenario=Tuple(scenario),
            probability_axis_order=(
                :alice_output, :bob_output, :alice_setting, :bob_setting
            ),
        ),
    )
    return (
        program=program,
        scenario=scenario,
        functional=functional,
        probability_shape=Tuple(scenario),
    )
end

function _nonlocal_behavior_from_coordinates(problem, coordinates)
    scenario = problem.scenario
    T = eltype(coordinates)
    probabilities = Array{T}(undef, Tuple(scenario))
    oa, ob, ma, mb = Tuple(scenario)
    for y in 1:mb, x in 1:ma, b in 1:ob, a in 1:oa
        probabilities[a, b, x, y] = coordinates[_nonlocal_probability_index(
            a, b, x, y, scenario
        )]
    end
    return probabilities
end

function _nonlocal_npa_bound_problem(
    functional::BellFunctional;
    level=1,
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    parsed = NPALevel(level)
    parsed.maximum_length > 0 ||
        throw(ArgumentError("a Bell quantum upper bound requires NPA level at least one"))
    cg = _nonlocal_full_functional_to_cg(functional)
    return _npa_build_problem(
        functional.scenario,
        parsed;
        cg_objective=cg,
        max_words=max_words,
        max_word_generation_work=max_word_generation_work,
        limits=limits,
    )
end

function _nonlocal_solve_upper_problem(problem, backend, regime::Symbol)
    program = problem isa NPAProblem ? problem.program : problem.program
    optimization = solve_optimization(program, backend)
    status = _nonlocal_status_from_optimization(optimization, :upper)
    upper = if optimization.status in (OptimizationOptimal, OptimizationFeasible)
        _nonlocal_solver_upper(optimization)
    else
        nothing
    end
    status === NonlocalValueNumericalUpperBound &&
        upper === nothing &&
        (status = NonlocalValueInvalidCertificate)
    witness = if optimization.primal === nothing
        nothing
    elseif problem isa NPAProblem
        if haskey(optimization.primal.views, :npa_moment)
            Matrix(optimization.primal.views.npa_moment)
        else
            nothing
        end
    elseif hasproperty(problem, :scenario)
        _nonlocal_behavior_from_coordinates(problem, optimization.primal.coordinates)
    elseif haskey(optimization.primal.views, :xor_gram)
        Matrix(optimization.primal.views.xor_gram)
    else
        nothing
    end
    message = if status === NonlocalValueNumericalUpperBound
        "the finite relaxation produced a numerical upper bound; it is not " *
        "promoted to an exact value or a theorem-level certificate"
    elseif status === NonlocalValueInvalidCertificate
        "the solver returned a primal relaxation value without an objective " *
        "bound or dual objective; a feasible maximization point is not an upper bound"
    elseif status === NonlocalValueBackendUnavailable
        "the package-owned optimization problem was built, but no backend was available"
    elseif status === NonlocalValueResourceLimit
        "the optimizer reached a configured resource limit"
    else
        "the optimizer did not return usable upper-bound evidence"
    end
    return NonlocalValueResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        regime,
        nothing,
        nothing,
        upper,
        false,
        false,
        false,
        nothing,
        nothing,
        problem,
        optimization,
        witness,
        (
            termination_status=optimization.termination_status,
            primal_status=optimization.primal_status,
            dual_status=optimization.dual_status,
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
            absolute_gap=optimization.absolute_gap,
            relative_gap=optimization.relative_gap,
            primal_objective=optimization.objective_value,
            solver_objective_bound=optimization.objective_bound,
            dual_objective=optimization.dual_objective_value,
        ),
        message,
    )
end

function _nonlocal_functional_game(functional::BellFunctional)
    scenario = functional.scenario
    ma, mb = scenario.alice_settings, scenario.bob_settings
    probability = Rational{BigInt}(1, ma * mb)
    probabilities = fill(probability, ma, mb)
    payoff_type = promote_type(eltype(functional.coefficients), Rational{BigInt})
    payoff = Array{payoff_type}(undef, Tuple(scenario))
    scale = ma * mb
    for index in CartesianIndices(payoff)
        payoff[index] = scale * functional.coefficients[index]
    end
    return NonlocalGame(probabilities, payoff; atol=0, rtol=0)
end

"""
    bell_inequality_bound(coefficients, scenario; notation=:cg,
                          regime=:classical, ...)

Compute an exact guarded classical maximum, a no-signalling LP bound, or a
finite-level NPA quantum upper bound.
"""
function bell_inequality_bound(
    coefficients,
    scenario::BellScenario;
    notation::Symbol=:collins_gisin,
    regime::Symbol=:classical,
    level=1,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_strategies=1_000_000,
    max_work=100_000_000,
    max_entries=1_000_000,
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    functional = BellFunctional(
        coefficients, scenario; notation=notation, max_entries=max_entries
    )
    if regime === :classical
        exact = _nonlocal_classical_optimum(
            _nonlocal_functional_game(functional);
            max_strategies=max_strategies,
            max_work=max_work,
        )
        strategy = DeterministicStrategy(exact.alice_strategy, exact.bob_strategy)
        return _nonlocal_exact_result(
            exact.value,
            :classical,
            strategy,
            (
                enumerated_party=exact.enumerated_party,
                strategies_evaluated=exact.strategies_evaluated,
                joint_strategy_count=exact.joint_strategy_count,
                scalar_work=exact.scalar_work,
                source_notation=notation,
            ),
        )
    elseif regime in (:no_signalling, :nosignal)
        problem = _nonlocal_no_signalling_problem(functional; limits=limits)
        return _nonlocal_solve_upper_problem(problem, backend, :no_signalling)
    elseif regime === :quantum
        problem = _nonlocal_npa_bound_problem(
            functional;
            level=level,
            max_words=max_words,
            max_word_generation_work=max_word_generation_work,
            limits=limits,
        )
        return _nonlocal_solve_upper_problem(problem, backend, :quantum_npa)
    end
    return throw(
        ArgumentError("regime must be :classical, :quantum, :no_signalling, or :nosignal")
    )
end

function bell_inequality_bound(coefficients, desc::AbstractVector; kwargs...)
    return bell_inequality_bound(coefficients, BellScenario(desc); kwargs...)
end

"""
    nonlocal_game_value(game; regime=:classical, ...)

Independently specified replacement for the `NonlocalGameValue` dependency
missing from the pinned QETLAB tree. Classical, no-signalling, and finite-NPA
branches share the reviewed Bell-functional implementations.
"""
function nonlocal_game_value(
    game::NonlocalGame;
    regime::Symbol=:classical,
    level=1,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_strategies=1_000_000,
    max_work=100_000_000,
    max_words=512,
    max_word_generation_work=2_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    if regime === :classical
        exact = _nonlocal_classical_optimum(
            game; max_strategies=max_strategies, max_work=max_work
        )
        return _nonlocal_exact_result(
            exact.value,
            :classical,
            DeterministicStrategy(exact.alice_strategy, exact.bob_strategy),
            (
                enumerated_party=exact.enumerated_party,
                strategies_evaluated=exact.strategies_evaluated,
                joint_strategy_count=exact.joint_strategy_count,
                scalar_work=exact.scalar_work,
            ),
        )
    end
    functional = BellFunctional(game)
    if regime in (:no_signalling, :nosignal)
        return _nonlocal_solve_upper_problem(
            _nonlocal_no_signalling_problem(functional; limits=limits),
            backend,
            :no_signalling,
        )
    elseif regime === :quantum
        return _nonlocal_solve_upper_problem(
            _nonlocal_npa_bound_problem(
                functional;
                level=level,
                max_words=max_words,
                max_word_generation_work=max_word_generation_work,
                limits=limits,
            ),
            backend,
            :quantum_npa,
        )
    end
    return throw(
        ArgumentError("regime must be :classical, :quantum, :no_signalling, or :nosignal")
    )
end

function _xor_validate(probabilities, winning_parity; atol, rtol, max_entries)
    probabilities isa AbstractMatrix{<:Real} ||
        throw(ArgumentError("XOR probabilities must be a real matrix"))
    winning_parity isa AbstractMatrix ||
        throw(ArgumentError("winning_parity must be a matrix"))
    _nonlocal_check_real_array(probabilities, "XOR probabilities")
    Base.require_one_based_indexing(winning_parity)
    size(probabilities) == size(winning_parity) || throw(
        DimensionMismatch("XOR probabilities and winning_parity must have equal size")
    )
    _nonlocal_check_budget(length(probabilities), max_entries, "max_entries")
    all(value -> value in (false, true, 0, 1), winning_parity) ||
        throw(ArgumentError("winning_parity must contain only zero/one entries"))
    isempty(probabilities) && throw(ArgumentError("XOR games need at least one question"))
    tolerance = _nonlocal_tolerance(
        eltype(probabilities), atol, rtol, maximum(abs, probabilities)
    )
    minimum(probabilities) < -tolerance &&
        throw(DomainError(minimum(probabilities), "XOR probabilities must be nonnegative"))
    total = sum(probabilities)
    abs(total - one(total)) <= tolerance || throw(
        ArgumentError(
            "XOR probabilities must sum to one; they are never normalized implicitly"
        ),
    )
    return copy(probabilities), Int.(winning_parity), tolerance
end

function _xor_classical(
    probabilities, parity; max_strategies=1_000_000, max_work=100_000_000
)
    ma, mb = size(probabilities)
    C = probabilities .* (one(eltype(probabilities)) .- 2 .* parity)
    enumerate_alice = ma <= mb
    count = big(2)^(enumerate_alice ? ma : mb)
    _nonlocal_check_budget(count, max_strategies, "max_strategies")
    work = count * ma * mb
    _nonlocal_check_budget(work, max_work, "max_work")
    T = eltype(C)
    best_bias = nothing
    best_alice = Int[]
    best_bob = Int[]
    if enumerate_alice
        for encoded in BigInt(0):(count - 1)
            alice_bits = _nonlocal_decode_strategy(encoded, 2, ma) .- 1
            alice_signs = one(T) .- 2 .* alice_bits
            bob_signs = Vector{T}(undef, mb)
            bias = zero(T)
            for y in 1:mb
                coefficient = sum(C[x, y] * alice_signs[x] for x in 1:ma)
                bob_signs[y] = coefficient >= zero(coefficient) ? one(T) : -one(T)
                bias += abs(coefficient)
            end
            if best_bias === nothing || bias > best_bias
                best_bias = bias
                best_alice = alice_bits .+ 1
                best_bob = Int.((one(T) .- bob_signs) ./ 2) .+ 1
            end
        end
    else
        for encoded in BigInt(0):(count - 1)
            bob_bits = _nonlocal_decode_strategy(encoded, 2, mb) .- 1
            bob_signs = one(T) .- 2 .* bob_bits
            alice_signs = Vector{T}(undef, ma)
            bias = zero(T)
            for x in 1:ma
                coefficient = sum(C[x, y] * bob_signs[y] for y in 1:mb)
                alice_signs[x] = coefficient >= zero(coefficient) ? one(T) : -one(T)
                bias += abs(coefficient)
            end
            if best_bias === nothing || bias > best_bias
                best_bias = bias
                best_alice = Int.((one(T) .- alice_signs) ./ 2) .+ 1
                best_bob = bob_bits .+ 1
            end
        end
    end
    value = (one(best_bias) + best_bias) / 2
    return (
        value=value,
        bias=best_bias,
        strategy=DeterministicStrategy(best_alice, best_bob),
        strategies_evaluated=count,
        enumerated_party=enumerate_alice ? :alice : :bob,
        joint_strategy_count=big(2)^(ma + mb),
        scalar_work=work,
    )
end

function _xor_quantum_problem(
    probabilities, parity; limits::OptimizationLimits=OptimizationLimits()
)
    T = eltype(probabilities)
    T in (Float32, Float64) || throw(
        ArgumentError("the XOR SDP currently requires Float32 or Float64 probabilities")
    )
    ma, mb = size(probabilities)
    dimension_big = BigInt(ma) + mb
    dimension_big <= typemax(Int) ||
        throw(ArgumentError("the XOR Gram dimension is too large for Int"))
    dimension_big <= limits.max_psd_dimension || throw(
        ArgumentError(
            "the XOR Gram matrix has dimension $dimension_big, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    variable_count_big = dimension_big^2
    variable_count_big <= limits.max_variables || throw(
        ArgumentError(
            "the XOR SDP needs $variable_count_big variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    dimension = Int(dimension_big)
    variable_count = Int(variable_count_big)
    gram = hermitian_variable(
        :xor_gram, dimension; variable_count=variable_count, coefficient_type=T
    )
    equalities = AffineEquality{T}[]
    for index in 1:dimension
        diagonal = _npa_entry_affine(gram, index, index, :real)
        push!(
            equalities,
            AffineEquality(
                _npa_affine_shift(diagonal, one(T)), Symbol(:xor_unit_diagonal_, index)
            ),
        )
    end
    coefficients = spzeros(T, variable_count)
    constant = one(T) / 2
    for x in 1:ma, y in 1:mb
        signed = probabilities[x, y] * (one(T) - 2 * parity[x, y]) / 2
        entry = _npa_entry_affine(gram, x, ma + y, :real)
        constant += signed * entry.constant
        coefficients .+= signed .* entry.coefficients
    end
    objective = AffineScalar(constant, coefficients; variable_count=variable_count)
    feasible = zeros(T, variable_count)
    feasible[1:dimension] .= one(T)
    program = SemidefiniteProgram(
        :xor_quantum_value,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=[gram],
        primal_views=[gram],
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(formulation=:tsirelson_xor_gram_sdp, alice_settings=ma, bob_settings=mb),
    )
    return (program=program, gram=gram, probabilities=probabilities, parity=parity)
end

"""
    xor_game_value(probabilities, winning_parity; regime=:classical, ...)

Compute the exact guarded classical value or the Tsirelson SDP quantum value
of a binary XOR game. The quantum result is retained as a numerical upper
bound with its Gram matrix and solver diagnostics.
"""
function xor_game_value(
    probabilities,
    winning_parity;
    regime::Symbol=:classical,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_strategies=1_000_000,
    max_work=100_000_000,
    max_entries=1_000_000,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_probabilities, parity, tolerance = _xor_validate(
        probabilities, winning_parity; atol=atol, rtol=rtol, max_entries=max_entries
    )
    if regime === :classical
        exact = _xor_classical(
            checked_probabilities, parity; max_strategies=max_strategies, max_work=max_work
        )
        return _nonlocal_exact_result(
            exact.value,
            :classical,
            exact.strategy,
            (
                bias=exact.bias,
                strategies_evaluated=exact.strategies_evaluated,
                enumerated_party=exact.enumerated_party,
                joint_strategy_count=exact.joint_strategy_count,
                scalar_work=exact.scalar_work,
                tolerance=tolerance,
            );
            certificate_kind=:exhaustive_xor_strategy_reduction,
        )
    elseif regime === :quantum
        problem = _xor_quantum_problem(checked_probabilities, parity; limits=limits)
        result = _nonlocal_solve_upper_problem(problem, backend, :quantum_xor)
        gram =
            if result.optimization_result === nothing ||
                result.optimization_result.primal === nothing
                nothing
            elseif haskey(result.optimization_result.primal.views, :xor_gram)
                Matrix(result.optimization_result.primal.views.xor_gram)
            else
                nothing
            end
        return NonlocalValueResult(
            _VALIDATED_CONSTRUCTOR_TOKEN,
            result.status,
            result.regime,
            result.value,
            result.lower_bound,
            result.upper_bound,
            result.exact,
            result.certified_lower,
            result.certified_upper,
            result.certificate_kind,
            result.strategy,
            problem,
            result.optimization_result,
            gram,
            merge(result.diagnostics, (tolerance=tolerance,)),
            result.message,
        )
    end
    return throw(ArgumentError("XOR regime must be :classical or :quantum"))
end
