# Bounded see-saw lower bounds, BCS composition, and the corrected rectangular
# qubit Bell relaxation.
#
# Source contracts:
#   NonlocalGameLB.m           afc7c3d4b6bdfb0d946269575e70f179fa253bb7c2c1c4361be69c8fa0d387c0
#   BCSGameLB.m                f0c3e510f545cd902cfc212835bd17d294c978dc142828f4cbd10d9548f4a298
#   BCSGameValue.m             42f9188009aeb5aa81609234b0fa5a26c07f7ba29d4a4fccb7de93f170833e63
#   BellInequalityMaxQubits.m  4dfe24e80ec4c2767a1f42a9bfbc3d6b1fa08682998bbdd911c12be362f2f90e
# at QETLAB revision d8589610f00cff106537268dee2e2a1153f3a601.
#
# QETLAB is BSD-2-Clause; see licenses/QETLAB-LICENSE.txt.

using Random: AbstractRNG, randn

"""
    NonlocalLowerBoundResult

Attained see-saw candidate evidence. `lower_bound` is populated only from
explicit returned state/measurement data whose numerical feasibility and
objective were checked. Because tolerance-feasible matrices are not an exact
feasibility certificate, `certified_lower` remains false. The candidate is
never called an exact quantum value.
"""
struct NonlocalLowerBoundResult{L,G,A,B,H,O,D,T}
    status::NonlocalValueStatus
    lower_bound::L
    certified_lower::Bool
    local_dimension::Int
    game::G
    alice_operators::A
    bob_measurements::B
    objective_history::H
    optimization_history::O
    iterations::Int
    converged::Bool
    diagnostics::D
    tolerance::T
    message::String

    function NonlocalLowerBoundResult(
        token::_ValidatedConstructorToken,
        status::NonlocalValueStatus,
        lower_bound,
        certified_lower::Bool,
        local_dimension::Int,
        game,
        alice_operators,
        bob_measurements,
        objective_history,
        optimization_history,
        iterations::Int,
        converged::Bool,
        diagnostics,
        tolerance,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        local_dimension > 0 || throw(ArgumentError("the local dimension must be positive"))
        iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
        certified_lower && throw(
            ArgumentError(
                "a tolerance-feasible see-saw candidate is not an exact lower-bound certificate",
            ),
        )
        if lower_bound === nothing
            alice_operators === nothing && bob_measurements === nothing || throw(
                ArgumentError(
                    "a result without a candidate cannot carry strategy operators"
                ),
            )
        else
            lower_bound isa Real && !(lower_bound isa Bool) && isfinite(lower_bound) ||
                throw(ArgumentError("a numerical lower-bound candidate must be finite"))
            alice_operators !== nothing && bob_measurements !== nothing || throw(
                ArgumentError("a lower-bound candidate must retain both operator families"),
            )
        end
        status === NonlocalValueNumericalLowerBound &&
            lower_bound === nothing &&
            throw(ArgumentError("a numerical-lower status requires a candidate"))
        converged &&
            status !== NonlocalValueNumericalLowerBound &&
            throw(ArgumentError("only a numerical-lower result may report convergence"))
        owned_alice = if alice_operators === nothing
            nothing
        else
            _nonlocal_read_only_nested(alice_operators)
        end
        owned_bob = if bob_measurements === nothing
            nothing
        else
            _nonlocal_read_only_nested(bob_measurements)
        end
        owned_objectives = _nonlocal_read_only(collect(objective_history))
        owned_optimizations = _nonlocal_read_only(collect(optimization_history))
        iterations <= length(owned_objectives) ||
            throw(ArgumentError("iteration count exceeds objective history"))
        tolerance isa Real &&
        !(tolerance isa Bool) &&
        isfinite(tolerance) &&
        tolerance >= zero(tolerance) ||
            throw(ArgumentError("see-saw tolerance must be finite and nonnegative"))
        return new{
            typeof(lower_bound),
            typeof(game),
            typeof(owned_alice),
            typeof(owned_bob),
            typeof(owned_objectives),
            typeof(owned_optimizations),
            typeof(diagnostics),
            typeof(tolerance),
        }(
            status,
            lower_bound,
            false,
            local_dimension,
            game,
            owned_alice,
            owned_bob,
            owned_objectives,
            owned_optimizations,
            iterations,
            converged,
            diagnostics,
            tolerance,
            String(message),
        )
    end
end

function Base.show(io::IO, result::NonlocalLowerBoundResult)
    return print(
        io,
        "NonlocalLowerBoundResult(status=",
        result.status,
        ", lower_bound=",
        result.lower_bound,
        ", iterations=",
        result.iterations,
        ", converged=",
        result.converged,
        ")",
    )
end

function _nonlocal_affine_matrix_combination(
    name::Symbol, matrices, coefficients; constant=nothing
)
    isempty(matrices) && throw(ArgumentError("at least one affine matrix is required"))
    length(matrices) == length(coefficients) || throw(
        DimensionMismatch("matrix and coefficient collections must have equal length")
    )
    first_matrix = first(matrices)
    T = _optimization_complex_real_type(eltype(first_matrix.constant))
    dimension = first_matrix.dimension
    variable_count = first_matrix.variable_count
    all(
        matrix -> matrix.dimension == dimension && matrix.variable_count == variable_count,
        matrices,
    ) || throw(DimensionMismatch("affine matrices are incompatible"))
    initial_constant = if constant === nothing
        spzeros(Complex{T}, dimension, dimension)
    else
        _optimization_check_hermitian(constant, "affine combination constant")
        size(constant) == (dimension, dimension) ||
            throw(DimensionMismatch("affine combination constant has the wrong size"))
        sparse(Complex{T}.(constant))
    end
    variables = Int[]
    coefficient_matrices = SparseMatrixCSC{Complex{T},Int}[]
    for (matrix, scalar) in zip(matrices, coefficients)
        scalar isa Real && !(scalar isa Bool) && isfinite(scalar) ||
            throw(ArgumentError("affine combination coefficients must be finite real"))
        for term in matrix.terms
            push!(variables, term.variable)
            push!(coefficient_matrices, convert(T, scalar) .* term.coefficient)
        end
        initial_constant += convert(T, scalar) .* matrix.constant
    end
    return HermitianAffineMatrix(
        name, initial_constant, variables, coefficient_matrices, variable_count
    )
end

function _nonlocal_write_hermitian_coordinates!(
    coordinates, first_variable::Int, matrix::AbstractMatrix
)
    dimension = size(matrix, 1)
    size(matrix, 2) == dimension || throw(DimensionMismatch("matrix must be square"))
    ishermitian(matrix) || throw(ArgumentError("matrix must be exactly Hermitian"))
    cursor = first_variable
    for index in 1:dimension
        coordinates[cursor] = real(matrix[index, index])
        cursor += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        coordinates[cursor] = real(matrix[row, column])
        cursor += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        coordinates[cursor] = imag(matrix[row, column])
        cursor += 1
    end
    return cursor
end

function _nonlocal_trace_objective(
    views, fixed_matrices, weights, variable_count::Int, ::Type{T}
) where {T<:Real}
    length(views) == length(fixed_matrices) == length(weights) ||
        throw(DimensionMismatch("objective collections must have equal length"))
    coefficients = spzeros(T, variable_count)
    constant = zero(T)
    for (view, fixed, weight) in zip(views, fixed_matrices, weights)
        for term in view.terms
            value = convert(T, weight * real(dot(fixed, term.coefficient)))
            iszero(value) && continue
            coefficients[term.variable] += value
        end
        constant += convert(T, weight * real(dot(fixed, view.constant)))
    end
    return AffineScalar(constant, coefficients; variable_count=variable_count)
end

function _nonlocal_float_game(game::NonlocalGame)
    T = promote_type(eltype(game.probabilities), eltype(game.payoff))
    T in (Float32, Float64) || throw(
        ArgumentError(
            "the see-saw SDP currently requires Float32 or Float64 game data; " *
            "the exact classical route supports exact inputs",
        ),
    )
    return T
end

function _nonlocal_random_povms(
    rng::AbstractRNG,
    ::Type{T},
    dimension::Int,
    outcomes::Int,
    settings::Int;
    max_initialization_attempts::Int,
) where {T<:AbstractFloat}
    for attempt in 1:max_initialization_attempts
        measurements = Array{Matrix{Complex{T}}}(undef, outcomes, settings)
        success = true
        for y in 1:settings
            raw = Matrix{Complex{T}}[]
            total = zeros(Complex{T}, dimension, dimension)
            for _ in 1:outcomes
                ginibre =
                    randn(rng, T, dimension, dimension) +
                    complex(zero(T), one(T)) .* randn(rng, T, dimension, dimension)
                effect = ginibre * adjoint(ginibre)
                push!(raw, effect)
                total += effect
            end
            decomposition = eigen(Hermitian(total))
            minimum(decomposition.values) > sqrt(eps(T)) || begin
                success = false
                break
            end
            inverse_root =
                decomposition.vectors *
                Diagonal(inv.(sqrt.(decomposition.values))) *
                adjoint(decomposition.vectors)
            for b in 1:outcomes
                measurements[b, y] = Matrix(Hermitian(inverse_root * raw[b] * inverse_root))
            end
        end
        success && return measurements, attempt
    end
    return nothing, max_initialization_attempts
end

function _nonlocal_seesaw_preflight(
    game::NonlocalGame,
    dimension::Int,
    initialization_attempts::Int;
    max_dense_entries,
    max_initialization_work,
    limits::OptimizationLimits,
)
    oa, ob, ma, mb = Tuple(game.scenario)
    dimension_big = BigInt(dimension)
    square = dimension_big^2
    cube = dimension_big^3
    alice_variables = (BigInt(oa) * ma + 1) * square
    bob_variables = BigInt(ob) * mb * square
    maximum_variables = max(alice_variables, bob_variables)
    maximum_variables <= limits.max_variables || throw(
        ArgumentError(
            "the see-saw model needs up to $maximum_variables variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    dimension <= limits.max_psd_dimension || throw(
        ArgumentError(
            "the see-saw PSD blocks have dimension $dimension, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    alice_psd_blocks = BigInt(oa) * ma + 1
    bob_psd_blocks = BigInt(ob) * mb
    maximum_psd_blocks = max(alice_psd_blocks, bob_psd_blocks)
    maximum_psd_blocks <= limits.max_psd_blocks || throw(
        ArgumentError(
            "the see-saw model needs up to $maximum_psd_blocks PSD blocks, exceeding " *
            "max_psd_blocks=$(limits.max_psd_blocks)",
        ),
    )
    alice_equalities = BigInt(ma) * square + 1
    bob_equalities = BigInt(mb) * square
    maximum_equalities = max(alice_equalities, bob_equalities)
    maximum_equalities <= limits.max_equalities || throw(
        ArgumentError(
            "the see-saw model needs up to $maximum_equalities equalities, exceeding " *
            "max_equalities=$(limits.max_equalities)",
        ),
    )
    hermitian_term_entries = 2 * square - dimension_big
    alice_stored_entries =
        alice_variables +
        BigInt(ma) * (oa + 1) * square +
        dimension_big +
        2 * alice_psd_blocks * hermitian_term_entries
    bob_stored_entries =
        bob_variables +
        BigInt(mb) * ob * square +
        2 * bob_psd_blocks * hermitian_term_entries
    real_block_entries_per_matrix = (2 * dimension_big) * (2 * dimension_big + 1) ÷ 2
    maximum_model_entries = max(
        alice_stored_entries,
        bob_stored_entries,
        alice_psd_blocks * real_block_entries_per_matrix,
        bob_psd_blocks * real_block_entries_per_matrix,
    )
    maximum_model_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "the see-saw model needs up to $maximum_model_entries stored or " *
            "real-block entries, exceeding max_model_entries=" *
            "$(limits.max_model_entries)",
        ),
    )
    dense_entries = square * (BigInt(ob) * mb + ob + 4)
    _nonlocal_check_budget(dense_entries, max_dense_entries, "max_dense_entries")
    initialization_work = BigInt(initialization_attempts) * mb * max(ob, 1) * cube
    _nonlocal_check_budget(
        initialization_work, max_initialization_work, "max_initialization_work"
    )
    return (
        dense_entries=dense_entries,
        initialization_work=initialization_work,
        maximum_variables=maximum_variables,
        maximum_model_entries=maximum_model_entries,
    )
end

function _nonlocal_alice_problem(
    game::NonlocalGame,
    bob_measurements,
    dimension::Int,
    ::Type{T};
    limits::OptimizationLimits,
) where {T<:AbstractFloat}
    oa, ob, ma, mb = Tuple(game.scenario)
    block = dimension^2
    variable_count = (oa * ma + 1) * block
    views = Matrix{HermitianAffineMatrix{T}}(undef, oa, ma)
    for x in 1:ma, a in 1:oa
        block_index = (x - 1) * oa + a
        views[a, x] = hermitian_variable(
            Symbol(:alice_, a, :_, x),
            dimension;
            first_variable=(block_index - 1) * block + 1,
            variable_count=variable_count,
            coefficient_type=T,
        )
    end
    rho = hermitian_variable(
        :shared_state,
        dimension;
        first_variable=oa * ma * block + 1,
        variable_count=variable_count,
        coefficient_type=T,
    )
    equalities = AffineEquality{T}[]
    for x in 1:ma
        combination = _nonlocal_affine_matrix_combination(
            Symbol(:alice_sum_minus_state_, x),
            vcat([views[a, x] for a in 1:oa], [rho]),
            vcat([one(T) for _ in 1:oa], [-one(T)]),
        )
        append!(
            equalities,
            hermitian_equalities(combination; name_prefix=Symbol(:alice_complete_, x)),
        )
    end
    push!(
        equalities,
        AffineEquality(_npa_affine_shift(trace_affine(rho), one(T)), :shared_state_trace),
    )
    objective_views = HermitianAffineMatrix{T}[]
    fixed = Matrix{Complex{T}}[]
    weights = T[]
    for x in 1:ma, y in 1:mb, a in 1:oa, b in 1:ob
        weight = convert(T, game.probabilities[x, y] * game.payoff[a, b, x, y])
        iszero(weight) && continue
        push!(objective_views, views[a, x])
        push!(fixed, bob_measurements[b, y])
        push!(weights, weight)
    end
    objective = _nonlocal_trace_objective(
        objective_views, fixed, weights, variable_count, T
    )
    feasible = zeros(T, variable_count)
    rho_value = Matrix{Complex{T}}(I, dimension, dimension) / convert(T, dimension)
    for x in 1:ma, a in 1:oa
        block_index = (x - 1) * oa + a
        _nonlocal_write_hermitian_coordinates!(
            feasible, (block_index - 1) * block + 1, rho_value / convert(T, oa)
        )
    end
    _nonlocal_write_hermitian_coordinates!(feasible, oa * ma * block + 1, rho_value)
    all_views = vcat(vec(views), [rho])
    program = SemidefiniteProgram(
        :nonlocal_seesaw_alice_state_step,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=all_views,
        primal_views=all_views,
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:fixed_bob_subnormalized_alice_state_effects,
            local_dimension=dimension,
            scenario=Tuple(game.scenario),
        ),
    )
    return (program=program, views=views, state_view=rho)
end

function _nonlocal_bob_problem(
    game::NonlocalGame,
    alice_operators,
    dimension::Int,
    ::Type{T};
    limits::OptimizationLimits,
) where {T<:AbstractFloat}
    oa, ob, ma, mb = Tuple(game.scenario)
    block = dimension^2
    variable_count = ob * mb * block
    views = Matrix{HermitianAffineMatrix{T}}(undef, ob, mb)
    for y in 1:mb, b in 1:ob
        block_index = (y - 1) * ob + b
        views[b, y] = hermitian_variable(
            Symbol(:bob_, b, :_, y),
            dimension;
            first_variable=(block_index - 1) * block + 1,
            variable_count=variable_count,
            coefficient_type=T,
        )
    end
    equalities = AffineEquality{T}[]
    identity_matrix = Matrix{Complex{T}}(I, dimension, dimension)
    for y in 1:mb
        sum_view = _nonlocal_affine_matrix_combination(
            Symbol(:bob_sum_, y), [views[b, y] for b in 1:ob], fill(one(T), ob)
        )
        append!(
            equalities,
            hermitian_equalities(
                sum_view; target=identity_matrix, name_prefix=Symbol(:bob_complete_, y)
            ),
        )
    end
    objective_views = HermitianAffineMatrix{T}[]
    fixed = Matrix{Complex{T}}[]
    weights = T[]
    for x in 1:ma, y in 1:mb, a in 1:oa, b in 1:ob
        weight = convert(T, game.probabilities[x, y] * game.payoff[a, b, x, y])
        iszero(weight) && continue
        push!(objective_views, views[b, y])
        push!(fixed, alice_operators[a, x])
        push!(weights, weight)
    end
    objective = _nonlocal_trace_objective(
        objective_views, fixed, weights, variable_count, T
    )
    feasible = zeros(T, variable_count)
    effect = identity_matrix / convert(T, ob)
    for y in 1:mb, b in 1:ob
        block_index = (y - 1) * ob + b
        _nonlocal_write_hermitian_coordinates!(
            feasible, (block_index - 1) * block + 1, effect
        )
    end
    all_views = vec(views)
    program = SemidefiniteProgram(
        :nonlocal_seesaw_bob_step,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=all_views,
        primal_views=all_views,
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:fixed_alice_bob_povm,
            local_dimension=dimension,
            scenario=Tuple(game.scenario),
        ),
    )
    return (program=program, views=views)
end

function _nonlocal_extract_views(optimization, views)
    optimization.primal === nothing && return nothing
    result = similar(views, Matrix{Complex{eltype(optimization.primal.coordinates)}})
    for index in eachindex(views)
        name = views[index].name
        haskey(optimization.primal.views, name) || return nothing
        result[index] = Matrix(getproperty(optimization.primal.views, name))
    end
    return result
end

function _nonlocal_attained_value(game::NonlocalGame, alice, bob)
    oa, ob, ma, mb = Tuple(game.scenario)
    T = promote_type(
        eltype(game.probabilities),
        eltype(game.payoff),
        eltype(first(alice)),
        eltype(first(bob)),
    )
    value = zero(typeof(real(zero(T))))
    for x in 1:ma, y in 1:mb, a in 1:oa, b in 1:ob
        value +=
            game.probabilities[x, y] *
            game.payoff[a, b, x, y] *
            real(dot(bob[b, y], alice[a, x]))
    end
    return value
end

function _nonlocal_lb_failure(
    status,
    game,
    dimension,
    tolerance,
    message;
    objective_history=Float64[],
    optimization_history=Any[],
    iterations=0,
    diagnostics=NamedTuple(),
)
    return NonlocalLowerBoundResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        nothing,
        false,
        dimension,
        game,
        nothing,
        nothing,
        objective_history,
        optimization_history,
        iterations,
        false,
        diagnostics,
        tolerance,
        message,
    )
end

"""
    nonlocal_game_lower_bound(rng, local_dimension, game; backend, ...)

Bounded alternating SDP heuristic. The explicit RNG is used only to initialize
Bob's POVMs. Each returned candidate is re-evaluated from the owned
tolerance-feasible operators; solver objectives alone are never accepted.
Numerical feasibility is not promoted to an exact lower-bound certificate.
"""
function nonlocal_game_lower_bound(
    rng::AbstractRNG,
    local_dimension,
    game::NonlocalGame;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_iterations=100,
    max_initialization_attempts=4,
    max_dense_entries=1_000_000,
    max_initialization_work=100_000_000,
    atol=1.0e-7,
    rtol=1.0e-7,
    limits::OptimizationLimits=OptimizationLimits(),
)
    dimension = _nonlocal_positive_integer(local_dimension, "local_dimension")
    iterations_limit = _nonlocal_nonnegative_integer(max_iterations, "max_iterations")
    initialization_limit = _nonlocal_positive_integer(
        max_initialization_attempts, "max_initialization_attempts"
    )
    if iterations_limit == 0
        T = promote_type(eltype(game.probabilities), eltype(game.payoff))
        tolerance = _nonlocal_tolerance(T <: Real ? T : Float64, atol, rtol)
        return _nonlocal_lb_failure(
            NonlocalValueResourceLimit,
            game,
            dimension,
            tolerance,
            "max_iterations=0 requests no heuristic work; the RNG was not consumed",
        )
    elseif backend isa NoOptimizationBackend
        T = promote_type(eltype(game.probabilities), eltype(game.payoff))
        tolerance = _nonlocal_tolerance(T <: Real ? T : Float64, atol, rtol)
        return _nonlocal_lb_failure(
            NonlocalValueBackendUnavailable,
            game,
            dimension,
            tolerance,
            "an explicit optimization backend is required; the RNG was not consumed",
        )
    end
    T = _nonlocal_float_game(game)
    tolerance = _nonlocal_tolerance(T, atol, rtol)
    preflight = _nonlocal_seesaw_preflight(
        game,
        dimension,
        initialization_limit;
        max_dense_entries=max_dense_entries,
        max_initialization_work=max_initialization_work,
        limits=limits,
    )
    bob, initialization_attempts = _nonlocal_random_povms(
        rng,
        T,
        dimension,
        game.scenario.bob_outputs,
        game.scenario.bob_settings;
        max_initialization_attempts=initialization_limit,
    )
    bob === nothing && return _nonlocal_lb_failure(
        NonlocalValueBackendFailure,
        game,
        dimension,
        tolerance,
        "bounded random POVM initialization failed";
        diagnostics=(initialization_attempts=initialization_attempts,),
    )
    objective_history = T[]
    optimization_history = Any[]
    best_value = nothing
    best_alice = nothing
    best_bob = nothing
    previous = nothing
    converged = false
    for iteration in 1:iterations_limit
        alice_problem = _nonlocal_alice_problem(game, bob, dimension, T; limits=limits)
        alice_optimization = solve_optimization(alice_problem.program, backend)
        push!(optimization_history, alice_optimization)
        if !(alice_optimization.status in (OptimizationOptimal, OptimizationFeasible))
            status = if alice_optimization.status === OptimizationLimit
                NonlocalValueResourceLimit
            else
                NonlocalValueBackendFailure
            end
            return NonlocalLowerBoundResult(
                _VALIDATED_CONSTRUCTOR_TOKEN,
                status,
                best_value,
                false,
                dimension,
                game,
                best_alice,
                best_bob,
                objective_history,
                optimization_history,
                iteration - 1,
                false,
                (
                    initialization_attempts=initialization_attempts,
                    failed_step=:alice_state,
                    termination_status=alice_optimization.termination_status,
                    primal_status=alice_optimization.primal_status,
                    dual_status=alice_optimization.dual_status,
                ),
                tolerance,
                "the Alice/state optimization step did not return a validated feasible point",
            )
        end
        alice = _nonlocal_extract_views(alice_optimization, alice_problem.views)
        alice === nothing && return _nonlocal_lb_failure(
            NonlocalValueInvalidCertificate,
            game,
            dimension,
            tolerance,
            "the backend omitted an Alice primal view";
            objective_history=objective_history,
            optimization_history=optimization_history,
            iterations=iteration - 1,
        )
        bob_problem = _nonlocal_bob_problem(game, alice, dimension, T; limits=limits)
        bob_optimization = solve_optimization(bob_problem.program, backend)
        push!(optimization_history, bob_optimization)
        if !(bob_optimization.status in (OptimizationOptimal, OptimizationFeasible))
            status = if bob_optimization.status === OptimizationLimit
                NonlocalValueResourceLimit
            else
                NonlocalValueBackendFailure
            end
            return NonlocalLowerBoundResult(
                _VALIDATED_CONSTRUCTOR_TOKEN,
                status,
                best_value,
                false,
                dimension,
                game,
                best_alice,
                best_bob,
                objective_history,
                optimization_history,
                iteration - 1,
                false,
                (
                    initialization_attempts=initialization_attempts,
                    failed_step=:bob,
                    termination_status=bob_optimization.termination_status,
                    primal_status=bob_optimization.primal_status,
                    dual_status=bob_optimization.dual_status,
                ),
                tolerance,
                "the Bob optimization step did not return a validated feasible point",
            )
        end
        candidate_bob = _nonlocal_extract_views(bob_optimization, bob_problem.views)
        candidate_bob === nothing && return _nonlocal_lb_failure(
            NonlocalValueInvalidCertificate,
            game,
            dimension,
            tolerance,
            "the backend omitted a Bob primal view";
            objective_history=objective_history,
            optimization_history=optimization_history,
            iterations=iteration - 1,
        )
        attained = convert(T, _nonlocal_attained_value(game, alice, candidate_bob))
        isfinite(attained) || return _nonlocal_lb_failure(
            NonlocalValueInvalidCertificate,
            game,
            dimension,
            tolerance,
            "the reconstructed strategy has a nonfinite objective";
            objective_history=objective_history,
            optimization_history=optimization_history,
            iterations=iteration - 1,
        )
        push!(objective_history, attained)
        if best_value === nothing || attained > best_value
            best_value = attained
            best_alice = map(copy, alice)
            best_bob = map(copy, candidate_bob)
        end
        bob = candidate_bob
        if previous !== nothing &&
            abs(attained - previous) <= tolerance * max(one(T), abs(attained))
            converged = true
            break
        end
        previous = attained
    end
    status = converged ? NonlocalValueNumericalLowerBound : NonlocalValueResourceLimit
    message = if converged
        "the bounded see-saw converged; the returned attained value is a " *
        "heuristic lower bound, not the exact quantum value"
    else
        "the bounded see-saw reached max_iterations; the best attained " *
        "strategy remains a valid numerical lower bound"
    end
    return NonlocalLowerBoundResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        best_value,
        false,
        dimension,
        game,
        best_alice,
        best_bob,
        objective_history,
        optimization_history,
        length(objective_history),
        converged,
        (
            initialization_attempts=initialization_attempts,
            final_change=if length(objective_history) >= 2
                abs(objective_history[end] - objective_history[end - 1])
            else
                nothing
            end,
            max_iterations=iterations_limit,
            preflight=preflight,
        ),
        tolerance,
        message,
    )
end

function nonlocal_game_lower_bound(
    rng::AbstractRNG,
    local_dimension,
    probabilities::AbstractMatrix,
    payoff::AbstractArray{<:Real,4};
    kwargs...,
)
    return nonlocal_game_lower_bound(
        rng, local_dimension, NonlocalGame(probabilities, payoff); kwargs...
    )
end

"""
    bcs_game_lower_bound(rng, local_dimension, constraints; ...)

Convert a validated BCS game and run [`nonlocal_game_lower_bound`](@ref).
"""
function bcs_game_lower_bound(
    rng::AbstractRNG,
    local_dimension,
    constraints;
    max_entries=1_000_000,
    coefficient_type::Type{T}=Float64,
    kwargs...,
) where {T<:AbstractFloat}
    T in (Float32, Float64) ||
        throw(ArgumentError("coefficient_type must be Float32 or Float64"))
    bcs = if constraints isa BCSGame
        constraints
    else
        BCSGame(constraints; max_entries=max_entries)
    end
    exact_game = nonlocal_game(bcs; max_entries=max_entries)
    # The BCS conversion is exact rational arithmetic. Entering a numerical
    # conic solver necessarily chooses a floating coefficient type, so that
    # choice is an explicit public keyword rather than a hidden promotion.
    game = NonlocalGame(
        T.(exact_game.probabilities),
        T.(exact_game.payoff);
        atol=zero(T),
        rtol=sqrt(eps(T)),
        max_entries=max_entries,
    )
    return nonlocal_game_lower_bound(rng, local_dimension, game; kwargs...)
end

"""
    bcs_game_value(constraints; regime=:classical, ...)

Classical/no-signalling/NPA dispatcher for BCS games. This independently
supplies the `NonlocalGameValue` capability missing from the pinned QETLAB
source tree.
"""
function bcs_game_value(
    constraints; regime::Symbol=:classical, max_entries=1_000_000, kwargs...
)
    bcs = if constraints isa BCSGame
        constraints
    else
        BCSGame(constraints; max_entries=max_entries)
    end
    game = nonlocal_game(bcs; max_entries=max_entries)
    return nonlocal_game_value(game; regime=regime, kwargs...)
end

"""
    BellQubitProblem

Corrected rectangular-setting qubit Bell PPT relaxation. The distinguished
Alice/Bob clone pair is a four-dimensional subsystem, so their shared
entanglement is retained; all remaining two-dimensional clone systems are
subject to one representative from every complementary PPT bipartition.
"""
struct BellQubitProblem{F,A,B,P,V,D}
    functional::F
    alice_values::A
    bob_values::B
    program::P
    state_view::V
    grouped_dimensions::D
    ppt_partitions::Tuple

    function BellQubitProblem(
        token::_ValidatedConstructorToken,
        functional::F,
        alice_values::A,
        bob_values::B,
        program::P,
        state_view::V,
        grouped_dimensions::D,
        ppt_partitions::Tuple,
    ) where {F,A,B,P,V,D}
        _require_validated_constructor_token(token)
        alice_values isa _NonlocalReadOnlyArray && bob_values isa _NonlocalReadOnlyArray ||
            throw(ArgumentError("qubit outcome values require read-only storage"))
        length(alice_values) == 2 && length(bob_values) == 2 ||
            throw(DimensionMismatch("qubit outcome values must each have length two"))
        !isempty(grouped_dimensions) && first(grouped_dimensions) == 4 ||
            throw(ArgumentError("the distinguished qubit pair must have dimension four"))
        all(>(0), grouped_dimensions) ||
            throw(ArgumentError("grouped dimensions must be positive"))
        length(Set(ppt_partitions)) == length(ppt_partitions) ||
            throw(ArgumentError("PPT partitions must not contain duplicates"))
        return new{F,A,B,P,V,D}(
            functional,
            alice_values,
            bob_values,
            program,
            state_view,
            grouped_dimensions,
            ppt_partitions,
        )
    end
end

"""
    BellQubitResult

Numerical upper bound and relaxation state. The returned state is the PPT
relaxation variable, not a physical two-qubit state or a strategy.
"""
struct BellQubitResult{U,P,S,O,D}
    status::NonlocalValueStatus
    upper_bound::U
    certified_upper::Bool
    problem::P
    relaxation_state::S
    optimization_result::O
    diagnostics::D
    message::String

    function BellQubitResult(
        token::_ValidatedConstructorToken,
        status::NonlocalValueStatus,
        upper_bound,
        certified_upper::Bool,
        problem,
        relaxation_state,
        optimization_result,
        diagnostics,
        message::AbstractString,
    )
        _require_validated_constructor_token(token)
        certified_upper && throw(
            ArgumentError("a numerical qubit relaxation is not an exact upper certificate"),
        )
        status === NonlocalValueNumericalUpperBound &&
            upper_bound === nothing &&
            throw(ArgumentError("a numerical-upper status requires an upper bound"))
        status !== NonlocalValueNumericalUpperBound &&
            upper_bound !== nothing &&
            throw(ArgumentError("an inconclusive qubit result cannot carry an upper bound"))
        owned_state = _nonlocal_owned_evidence(relaxation_state)
        return new{
            typeof(upper_bound),
            typeof(problem),
            typeof(owned_state),
            typeof(optimization_result),
            typeof(diagnostics),
        }(
            status,
            upper_bound,
            false,
            problem,
            owned_state,
            optimization_result,
            diagnostics,
            String(message),
        )
    end
end

function _bell_qubit_local_projector(
    grouped_dimensions::Tuple,
    subsystem::Int,
    outcome::Int,
    party_in_pair::Symbol,
    ::Type{T},
) where {T<:Real}
    total = prod(grouped_dimensions)
    diagonal = zeros(T, total)
    layout = SubsystemLayout(grouped_dimensions)
    for linear in 1:total
        basis = linear_to_basis(linear, layout)
        selected = if subsystem == 1
            pair_value = basis[1] - 1
            party_in_pair === :alice ? div(pair_value, 2) + 1 : mod(pair_value, 2) + 1
        else
            basis[subsystem]
        end
        selected == outcome && (diagonal[linear] = one(T))
    end
    return spdiagm(0 => Complex{T}.(diagonal))
end

function _bell_qubit_partitions(group_count::Int; max_ppt_constraints)
    group_count <= 1 && return ()
    total = big(2)^(group_count - 1) - 1
    _nonlocal_check_budget(total, max_ppt_constraints, "max_ppt_constraints")
    partitions = Vector{Tuple{Vararg{Int}}}()
    full_mask = big(2)^group_count - 1
    for mask in BigInt(1):(full_mask - 1)
        complement = xor(mask, full_mask)
        mask < complement || continue
        subset = Tuple(
            index for index in 1:group_count if !iszero(mask & (big(1) << (index - 1)))
        )
        push!(partitions, subset)
    end
    return Tuple(partitions)
end

function _bell_qubit_problem(
    joint_coefficients::AbstractMatrix{<:Real},
    alice_coefficients::AbstractVector{<:Real},
    bob_coefficients::AbstractVector{<:Real},
    alice_values::AbstractVector{<:Real},
    bob_values::AbstractVector{<:Real};
    max_dimension=256,
    max_ppt_constraints=127,
    limits::OptimizationLimits=OptimizationLimits(),
)
    for (value, name) in (
        (joint_coefficients, "joint_coefficients"),
        (alice_coefficients, "alice_coefficients"),
        (bob_coefficients, "bob_coefficients"),
        (alice_values, "alice_values"),
        (bob_values, "bob_values"),
    )
        _nonlocal_check_real_array(value, name)
    end
    ma, mb = size(joint_coefficients)
    length(alice_coefficients) == ma ||
        throw(DimensionMismatch("alice_coefficients must have length $ma"))
    length(bob_coefficients) == mb ||
        throw(DimensionMismatch("bob_coefficients must have length $mb"))
    length(alice_values) == 2 && length(bob_values) == 2 || throw(
        DimensionMismatch("qubit Bell relaxations require exactly two outcome values")
    )
    ma > 0 && mb > 0 || throw(ArgumentError("each party needs at least one setting"))
    T = promote_type(
        eltype(joint_coefficients),
        eltype(alice_coefficients),
        eltype(bob_coefficients),
        eltype(alice_values),
        eltype(bob_values),
    )
    T in (Float32, Float64) || throw(
        ArgumentError("the qubit Bell SDP currently requires Float32 or Float64 data")
    )
    grouped_dimensions = Tuple(vcat([4], fill(2, ma + mb - 2)))
    total_dimension_big = prod(BigInt, grouped_dimensions)
    _nonlocal_check_budget(total_dimension_big, max_dimension, "max_dimension")
    total_dimension_big <= limits.max_psd_dimension || throw(
        ArgumentError(
            "the qubit relaxation PSD dimension $total_dimension_big exceeds " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    variable_count_big = total_dimension_big^2
    variable_count_big <= limits.max_variables || throw(
        ArgumentError(
            "the qubit relaxation needs $variable_count_big variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    total_dimension = Int(total_dimension_big)
    variable_count = Int(variable_count_big)
    partitions = _bell_qubit_partitions(
        length(grouped_dimensions); max_ppt_constraints=max_ppt_constraints
    )
    BigInt(length(partitions)) + 1 <= limits.max_psd_blocks || throw(
        ArgumentError(
            "the qubit relaxation needs $(length(partitions) + 1) PSD blocks, " *
            "exceeding max_psd_blocks=$(limits.max_psd_blocks)",
        ),
    )
    state = hermitian_variable(
        :bell_qubit_relaxation_state,
        total_dimension;
        variable_count=variable_count,
        coefficient_type=T,
    )
    alice_effects = Matrix{SparseMatrixCSC{Complex{T},Int}}(undef, 2, ma)
    bob_effects = Matrix{SparseMatrixCSC{Complex{T},Int}}(undef, 2, mb)
    for x in 1:ma, outcome in 1:2
        subsystem = x == 1 ? 1 : x
        alice_effects[outcome, x] = _bell_qubit_local_projector(
            grouped_dimensions, subsystem, outcome, :alice, T
        )
    end
    for y in 1:mb, outcome in 1:2
        subsystem = y == 1 ? 1 : ma + y - 1
        bob_effects[outcome, y] = _bell_qubit_local_projector(
            grouped_dimensions, subsystem, outcome, :bob, T
        )
    end
    objective_matrix = spzeros(Complex{T}, total_dimension, total_dimension)
    for x in 1:ma, y in 1:mb, a in 1:2, b in 1:2
        coefficient = joint_coefficients[x, y] * alice_values[a] * bob_values[b]
        y == 1 && (coefficient += alice_coefficients[x] * alice_values[a])
        x == 1 && (coefficient += bob_coefficients[y] * bob_values[b])
        iszero(coefficient) && continue
        objective_matrix +=
            convert(T, coefficient) .* (alice_effects[a, x] * bob_effects[b, y])
    end
    objective = _nonlocal_trace_objective(
        [state], [objective_matrix], [one(T)], variable_count, T
    )
    equalities = [
        AffineEquality(
            _npa_affine_shift(trace_affine(state), one(T)), :bell_qubit_state_trace
        ),
    ]
    psd_constraints = HermitianAffineMatrix{T}[state]
    for (index, partition) in enumerate(partitions)
        push!(
            psd_constraints,
            partial_transpose_affine(
                state,
                grouped_dimensions;
                systems=partition,
                name=Symbol(:bell_qubit_ppt_, index),
            ),
        )
    end
    feasible = zeros(T, variable_count)
    _nonlocal_write_hermitian_coordinates!(
        feasible,
        1,
        Matrix{Complex{T}}(I, total_dimension, total_dimension) /
        convert(T, total_dimension),
    )
    program = SemidefiniteProgram(
        :bell_inequality_qubit_ppt_bound,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=psd_constraints,
        primal_views=[state],
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:rectangular_qubit_clone_ppt_relaxation,
            alice_settings=ma,
            bob_settings=mb,
            grouped_dimensions=grouped_dimensions,
            ppt_constraint_count=length(partitions),
            upstream_rectangular_loop_corrected=true,
        ),
    )
    scenario = BellScenario(2, 2, ma, mb)
    full_coefficients = zeros(T, 2, 2, ma, mb)
    for x in 1:ma, y in 1:mb, a in 1:2, b in 1:2
        full_coefficients[a, b, x, y] =
            joint_coefficients[x, y] * alice_values[a] * bob_values[b] +
            alice_coefficients[x] * alice_values[a] / mb +
            bob_coefficients[y] * bob_values[b] / ma
    end
    functional = BellFunctional(full_coefficients, scenario; notation=:full_probability)
    return BellQubitProblem(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        functional,
        _nonlocal_read_only(collect(alice_values)),
        _nonlocal_read_only(collect(bob_values)),
        program,
        state,
        grouped_dimensions,
        partitions,
    )
end

"""
    bell_inequality_qubit_bound(joint, alice, bob, alice_values, bob_values;
                                backend, ...)

Build and solve the corrected rectangular-setting qubit PPT relaxation. Unlike
the pinned source, Bob's loops and dimensions use Bob's setting count.
"""
function bell_inequality_qubit_bound(
    joint_coefficients,
    alice_coefficients,
    bob_coefficients,
    alice_values,
    bob_values;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    max_dimension=256,
    max_ppt_constraints=127,
    limits::OptimizationLimits=OptimizationLimits(),
)
    problem = _bell_qubit_problem(
        joint_coefficients,
        alice_coefficients,
        bob_coefficients,
        alice_values,
        bob_values;
        max_dimension=max_dimension,
        max_ppt_constraints=max_ppt_constraints,
        limits=limits,
    )
    optimization = solve_optimization(problem.program, backend)
    status = _nonlocal_status_from_optimization(optimization, :upper)
    upper = if optimization.status in (OptimizationOptimal, OptimizationFeasible)
        _nonlocal_solver_upper(optimization)
    else
        nothing
    end
    status === NonlocalValueNumericalUpperBound &&
        upper === nothing &&
        (status = NonlocalValueInvalidCertificate)
    state =
        if optimization.primal !== nothing &&
            haskey(optimization.primal.views, :bell_qubit_relaxation_state)
            Matrix(optimization.primal.views.bell_qubit_relaxation_state)
        else
            nothing
        end
    message = if status === NonlocalValueInvalidCertificate
        "the solver returned a primal relaxation value without an objective " *
        "bound or dual objective; it was not mislabeled as an upper bound"
    elseif upper === nothing
        "the qubit PPT relaxation did not return usable upper-bound evidence"
    else
        "the corrected rectangular qubit PPT relaxation returned a numerical " *
        "upper bound; it is not an exact qubit value"
    end
    return BellQubitResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        upper,
        false,
        problem,
        state,
        optimization,
        (
            alice_settings=size(joint_coefficients, 1),
            bob_settings=size(joint_coefficients, 2),
            grouped_dimensions=problem.grouped_dimensions,
            ppt_constraint_count=length(problem.ppt_partitions),
            upstream_rectangular_loop_corrected=true,
            termination_status=optimization.termination_status,
            primal_status=optimization.primal_status,
            dual_status=optimization.dual_status,
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
            primal_objective=optimization.objective_value,
            solver_objective_bound=optimization.objective_bound,
            dual_objective=optimization.dual_objective_value,
        ),
        message,
    )
end
