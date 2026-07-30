# Optional positional-backend convenience for the package-owned block-
# positivity analysis.

function QET.is_block_positive(
    rng::QET.Random.AbstractRNG, operator, backend::QET.JuMPBackend; kwargs...
)
    return QET.is_block_positive(rng, operator; backend, kwargs...)
end
