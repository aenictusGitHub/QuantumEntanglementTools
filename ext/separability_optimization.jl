# Row-specific optional-backend entry points. This file is included from the
# package's JuMP extension only; the core formulations remain solver-neutral.

using Random: AbstractRNG

const _SEPARABILITY_OPTIMIZATION_EXTENSION_LOADED = true

function QET.local_distinguishability(states, dims, backend::QET.JuMPBackend; kwargs...)
    return QET.local_distinguishability(states, dims; backend=backend, kwargs...)
end

function QET.is_separable(
    rho::AbstractMatrix{<:Number}, dims, backend::QET.JuMPBackend; kwargs...
)
    return QET.is_separable(rho, dims; backend=backend, kwargs...)
end

function QET.is_separable(
    rng::AbstractRNG,
    rho::AbstractMatrix{<:Number},
    dims,
    backend::QET.JuMPBackend;
    kwargs...,
)
    return QET.is_separable(rng, rho, dims; backend=backend, kwargs...)
end

function QET.upb_sep_distinguishable(local_factors, backend::QET.JuMPBackend; kwargs...)
    return QET.upb_sep_distinguishable(local_factors; backend=backend, kwargs...)
end
