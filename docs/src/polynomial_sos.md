# Polynomial SOS hierarchy

[`polynomial_sos_bounds`](@ref) implements the SDP hierarchy used by QETLAB
`PolynomialSOS` for real, even-degree homogeneous polynomials on the unit
sphere. It keeps the conic outer relaxation separate from sampled feasible
inner values.

## A complete small example

Consider

```math
p(x_1,x_2)=x_1^2+2x_1x_2+3x_2^2.
```

Its maximum and minimum on $x_1^2+x_2^2=1$ are $2+\sqrt{2}$ and
$2-\sqrt{2}$.

```julia
using Hypatia
using JuMP
using QuantumEntanglementTools
using Random

p = HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)

result = polynomial_sos_bounds(
    MersenneTwister(7),
    p,
    backend;
    level=0,
    sense=:max,
    inner_samples=64,
)
```

Read the two sides independently:

```julia
result.outer_bound
result.inner_bound
result.best_point
result.moment_matrix
result.optimization_result
```

For maximization, the SOS value is an outer upper bound and every returned
sample is an attained lower bound. For minimization, those directions reverse.
The sampled point is feasible evidence, but sampling quality is heuristic.

## Solver-neutral model

[`polynomial_sos_problem`](@ref) returns a
[`SemidefiniteProgram`](@ref) without selecting a solver. If $p$ has degree
$2d$ and the hierarchy level is $k$, let $V$ embed the normalized
degree-$(d+k)$ symmetric monomial basis into
$(\mathbb{R}^n)^{\otimes(d+k)}$. The primal moment matrix $\rho$ obeys

```math
\rho\succeq0,\qquad
\mathop{\mathrm{tr}}\rho=1,\qquad
\left(V\rho V^\mathsf{T}\right)^{T_1}
=V\rho V^\mathsf{T}.
```

The objective is the trace pairing between $\rho$ and the compact symmetric
matrix returned by [`polynomial_as_matrix`](@ref). `sense=:max` maximizes that
pairing; `sense=:min` minimizes it.

The package reconstructs the moment matrix from primal coordinates and
retains the generic optimizer result: termination, primal and dual statuses,
objective/bound/gap values, residuals, solver metadata, limits, and warnings.

## Status and certificate semantics

No optimizer is implicit. With no backend, `result.status` is
`OptimizationBackendUnavailable`; requested samples can still produce an
inner value.

`OptimizationOptimal` means the solver reported a primal-dual optimum
consistent with the package's available residual checks. Floating conic
evidence is not automatically a theorem-level SOS certificate, so
`result.certified_outer` remains `false`. Limited, inaccurate, inconsistent,
unsupported, malformed, and failed outcomes remain distinguishable, and do
not fabricate an outer value.

A degree-zero homogeneous polynomial is handled exactly without a solver.

## Explicit randomness and targets

The first argument is an explicit `rng::AbstractRNG`. `inner_samples` is an
exact deterministic count; it is never derived from elapsed solver time. The
global random stream is not touched.

If a computed outer bound already proves the requested `target`, sampling is
skipped without consuming the supplied RNG. The result records
`target_status`, `samples_requested`, and `samples_evaluated`.

## Resource policy

Before allocating a symmetric basis or affine SDP, the builder checks:

- $n^{d+k}$ against `max_full_dimension`;
- symmetric-basis stored entries;
- moment variables and PSD dimension;
- partial-transpose equality count;
- affine construction work;
- the standard `OptimizationLimits`;
- polynomial degree, monomial, exponent-table, dense-entry, nonzero, sample,
  and work limits.

The optional extension uses the documented complex-Hermitian real-block
translation and therefore requires `allow_densify=true`. Although the
polynomial data are real, the package's central Hermitian coordinate model is
used; partial-transpose invariance removes irrelevant imaginary directions.

## Differences from QETLAB

This implementation follows the pinned primal hierarchy and coefficient
ordering, but deliberately replaces two unsafe behaviors:

- randomized inner-bound work has a fixed sample count instead of a
  wall-clock-derived loop;
- solver failures and numerical boundaries return structured statuses instead
  of an unchecked scalar.

Compatibility mode preserves the useful argument order while requiring an
explicit RNG and backend for the corresponding branches.
