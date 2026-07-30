```@meta
CurrentModule = QuantumEntanglementTools
DocTestSetup = quote
    using QuantumEntanglementTools
    using LinearAlgebra
end
```

# Minimum-error state discrimination

`state_distinguishability` computes or bounds the largest success probability
for identifying one member of a known quantum-state ensemble. It accepts
either a tuple/vector of density matrices or a matrix whose columns are pure
states.

For normalized states $\rho_j$ with priors $p_j$, the optimization is

```math
\begin{aligned}
\text{maximize}\quad & \sum_j p_j\,\mathrm{tr}(M_j\rho_j),\\
\text{subject to}\quad & M_j \succeq 0,\\
& \sum_j M_j = I.
\end{aligned}
```

The effects $M_j$ form a positive-operator-valued measure (POVM). The
solver-neutral builder `state_discrimination_problem` returns this affine SDP
as a package-owned `SemidefiniteProgram`; it never stores JuMP variables or
optimizer objects.

## Solver-free certificates

Two states use the Helstrom formula

```math
p_{\mathrm{success}}
= \frac{1+\lVert p_1\rho_1-p_2\rho_2\rVert_1}{2}.
```

The returned POVM is the positive spectral projector of
$p_1\rho_1-p_2\rho_2$ and its complement. A deterministic prior and an
exactly orthogonal support decomposition also produce solver-free
certificates.

```jldoctest
ket0 = ComplexF64[1, 0]
ketplus = ComplexF64[1, 1] / sqrt(2)
result = state_distinguishability(hcat(ket0, ketplus))

(
    result.status ==
        QuantumEntanglementTools.StateDiscriminationHelstromOptimal,
    isapprox(result.success_probability, (1 + inv(sqrt(2))) / 2),
)

# output
(true, true)
```

`result.measurement` is returned only after independent positivity,
Hermiticity, completeness, objective, and probability checks.

## Three or more nonorthogonal states

The general branch needs an explicit optional backend:

```julia
using Hypatia, JuMP, QuantumEntanglementTools

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)

omega = cis(2pi / 3)
trine = hcat(
    ComplexF64[1, 1] / sqrt(2),
    ComplexF64[1, omega] / sqrt(2),
    ComplexF64[1, omega^2] / sqrt(2),
)
result = state_distinguishability(trine; backend)

result.lower_bound
result.upper_bound
result.measurement
result.optimization_result
```

For the trine ensemble the optimum is $2/3$. The dedicated optional test
environment checks the POVM, primal and dual residuals, and this analytic
value with Hypatia, then cross-checks it with SCS.

Without a backend, the model is still inspectable and the outcome is explicit:

```jldoctest
trine = hcat(
    ComplexF64[1, 1] / sqrt(2),
    ComplexF64[1, cis(2pi / 3)] / sqrt(2),
    ComplexF64[1, cis(4pi / 3)] / sqrt(2),
)
unavailable = state_distinguishability(trine)

(
    unavailable.status,
    unavailable.success_probability,
    unavailable.problem.program.variable_count,
)

# output
(QuantumEntanglementTools.StateDiscriminationBackendUnavailable, nothing, 12)
```

## Reading the result

`StateDiscriminationResult` deliberately separates:

- `lower_bound`: the success probability of a residual-checked POVM;
- `upper_bound`: an analytic value or a validated dual/optimizer bound;
- `success_probability`: present only when the available bounds agree within
  the recorded tolerance;
- `measurement`: a residual-checked POVM;
- `dual_operator`: reconstructed dual data and slack diagnostics when
  available;
- `optimization_result`: termination, primal/dual statuses, bounds, gaps,
  residuals, iterations, timing, optimizer versions, and options.

A solver-reported optimum is not relabeled as an exact mathematical
certificate. Limits, malformed optimizers, missing backends, and numerical
boundaries retain structured inconclusive statuses.

## Validation, storage, and limits

Density matrices must be finite, normalized, exactly Hermitian, positive
semidefinite, square, and equal in dimension. Pure columns must be finite and
normalized. Priors must be finite, nonnegative, have the same length as the
ensemble, and sum to one within the stated tolerance.

No input is normalized, symmetrized, clipped, or projected. Sparse inputs
require `allow_densify=true` because state validation and POVM spectral checks
use dense eigendecompositions. The optional real-block complex-PSD translation
also requires a backend with `allow_densify=true`.

`max_dimension`, `max_states`, `max_variables`, and `OptimizationLimits`
preflight model growth before allocation. For Hilbert-space dimension $d$
and $n$ states, the primal uses $n d^2$ real coordinates, $n$ PSD blocks
of complex dimension $d$, and $d^2$ scalar completeness equalities.

## QETLAB migration and provenance

This API covers the executable contract of pinned QETLAB
`Distinguishability.m`. `MATLABCompat.Distinguishability` delegates to the same
reviewed kernel.

The native and compatibility APIs intentionally do not reproduce QETLAB's
silent normalization of density matrices and pure columns. They also retain
optimizer statuses and residuals instead of returning an unchecked CVX scalar.

The pinned two-pure-state scalar formula has a separate unequal-prior defect:
it uses

```math
\frac{1}{2}
+\frac{1}{2}\sqrt{
2(p_1^2+p_2^2)-4p_1p_2
\lvert\langle\psi_1,\psi_2\rangle\rvert^2
},
```

which agrees with Helstrom only for equal priors and can exceed one. Its own
returned POVM attains the correct smaller value. The committed Octave fixture
records both numbers; the Julia APIs return the residual-checked Helstrom
value.

The solver-independent formulation and validation were checked against
[John Watrous, *The Theory of Quantum Information*, Section
3.1](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).

The complete docstrings are collected in the
[Julia-native optimization API](api/native_optimization.md).
