# Optional positional-backend convenience for the package-owned S(k)-norm
# hierarchy. The actual solver translation is the generic, audited
# SemidefiniteProgram-to-JuMP bridge.

function QET.sk_operator_norm(
    rng::QET.Random.AbstractRNG, operator, backend::QET.JuMPBackend; kwargs...
)
    return QET.sk_operator_norm(rng, operator; backend, kwargs...)
end
