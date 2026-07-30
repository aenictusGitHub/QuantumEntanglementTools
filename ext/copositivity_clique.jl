# Row-specific optional-JuMP entry points. Model translation remains the
# generic package-owned SemidefiniteProgram path in this extension.

function QET.copositivity_criterion(
    rng, matrix::AbstractMatrix, backend::QET.JuMPBackend; kwargs...
)
    return QET.copositivity_criterion(rng, matrix; backend=backend, kwargs...)
end

function QET.clique_number_bounds(
    rng, adjacency::AbstractMatrix, backend::QET.JuMPBackend; kwargs...
)
    return QET.clique_number_bounds(rng, adjacency; backend=backend, kwargs...)
end

const _COPOSITIVITY_CLIQUE_EXTENSION_LOADED = true
