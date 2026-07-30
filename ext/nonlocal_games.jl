# Optional JuMP convenience dispatch for package-owned nonlocal-game models.
# All model construction, limits, status semantics, and evidence interpretation
# remain in the core source; this file only makes an explicit JuMPBackend
# available positionally after the package extension has loaded.

const _NONLOCAL_GAMES_EXTENSION_LOADED = true

function QET.npa_membership(behavior, backend::QET.JuMPBackend; kwargs...)
    return QET.npa_membership(behavior; backend=backend, kwargs...)
end

function QET.xor_game_value(
    probabilities, winning_parity, backend::QET.JuMPBackend; kwargs...
)
    return QET.xor_game_value(probabilities, winning_parity; backend=backend, kwargs...)
end

function QET.bell_inequality_bound(
    coefficients, scenario, backend::QET.JuMPBackend; kwargs...
)
    return QET.bell_inequality_bound(coefficients, scenario; backend=backend, kwargs...)
end

function QET.bell_inequality_qubit_bound(
    joint_coefficients,
    alice_coefficients,
    bob_coefficients,
    alice_values,
    bob_values,
    backend::QET.JuMPBackend;
    kwargs...,
)
    return QET.bell_inequality_qubit_bound(
        joint_coefficients,
        alice_coefficients,
        bob_coefficients,
        alice_values,
        bob_values;
        backend=backend,
        kwargs...,
    )
end

function QET.nonlocal_game_value(game, backend::QET.JuMPBackend; kwargs...)
    return QET.nonlocal_game_value(game; backend=backend, kwargs...)
end

function QET.nonlocal_game_lower_bound(
    rng, local_dimension, game, backend::QET.JuMPBackend; kwargs...
)
    return QET.nonlocal_game_lower_bound(
        rng, local_dimension, game; backend=backend, kwargs...
    )
end

function QET.bcs_game_lower_bound(
    rng, local_dimension, constraints, backend::QET.JuMPBackend; kwargs...
)
    return QET.bcs_game_lower_bound(
        rng, local_dimension, constraints; backend=backend, kwargs...
    )
end

function QET.bcs_game_value(constraints, backend::QET.JuMPBackend; kwargs...)
    return QET.bcs_game_value(constraints; backend=backend, kwargs...)
end
