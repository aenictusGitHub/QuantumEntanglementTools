# Row-specific optional-extension entry point. The actual model translation is
# the generic package-owned SemidefiniteProgram path in this extension.

function QET.polynomial_sos_bounds(
    rng, polynomial::QET.HomogeneousPolynomial, backend::QET.JuMPBackend; kwargs...
)
    return QET.polynomial_sos_bounds(rng, polynomial; backend=backend, kwargs...)
end

const _POLYNOMIAL_SOS_EXTENSION_LOADED = true

function QET.polynomial_sos_bounds(
    rng,
    coefficients::AbstractVector,
    variables,
    half_degree,
    level,
    backend::QET.JuMPBackend;
    kwargs...,
)
    return QET.polynomial_sos_bounds(
        rng, coefficients, variables, half_degree, level; backend=backend, kwargs...
    )
end
