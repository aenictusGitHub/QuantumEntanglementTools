# Matsumoto fidelity: numeric value and SDP model

The Matsumoto fidelity of density matrices $\rho$ and $\sigma$ is

```math
F_{\mathrm{M}}(\rho,\sigma)
=
\mathop{\mathrm{tr}}(\rho\mathbin{\#}\sigma),
```

where $\#$ is the Kubo--Ando matrix geometric mean. It is not the Uhlmann
fidelity returned by [`fidelity`](@ref).

## Direct numeric evaluation

Use [`matsumoto_fidelity`](@ref) when both states are known numerically:

```julia
using LinearAlgebra
using QuantumEntanglementTools

rho = Diagonal([0.7, 0.3])
sigma = Diagonal([0.4, 0.6])
value = matsumoto_fidelity(rho, sigma)
```

The numeric algorithm uses support-aware factorizations and does not form an
explicit inverse or add a regularizing identity. It validates both density
matrices without normalizing, symmetrizing, clipping, or repairing them.
General sparse inputs require `allow_densify=true`.

## Explicit semidefinite lift

QETLAB also accepts CVX expressions and silently opens a nested model.
QuantumEntanglementTools exposes that operation as package-owned data:

```math
F_{\mathrm{M}}(\rho,\sigma)
=
\max_X\left\{
\mathop{\mathrm{tr}}X:
\begin{bmatrix}
\rho & X\\
X & \sigma
\end{bmatrix}\succeq0,\quad X=X^\dagger
\right\}.
```

For fixed numeric states, construct and solve the SDP explicitly:

```julia
using Hypatia
using JuMP

problem = matsumoto_fidelity_problem(rho, sigma)
backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)
result = solve_optimization(problem, backend)
```

No optimizer is selected by default. Without `backend`, the result has status
`OptimizationBackendUnavailable`; it does not fabricate a value.

Inspect solver evidence before using the result:

```julia
result.status
result.objective_value
result.objective_bound
result.primal_residual
result.dual_residual
result.absolute_gap
result.primal.views.matsumoto_fidelity_coupling
```

An `OptimizationOptimal` status means the reported solver outcome is
consistent with the package's available residual checks. It is not
automatically a theorem-level certificate.

## Composing with affine state variables

[`matsumoto_fidelity_model`](@ref) also accepts two
[`HermitianAffineMatrix`](@ref) objects with the same coordinate count. It
appends the Hermitian coupling coordinates and returns:

- the affine objective `tr(X)`;
- the block PSD constraint;
- the named coupling view;
- the old and lifted coordinate counts.

[`matsumoto_fidelity_problem`](@ref) can lift caller-supplied affine
equalities, intervals, PSD constraints, and primal views automatically. The
caller must impose positivity, trace, and any other density-matrix conditions
on affine inputs; symbolic model data cannot be validated as a particular
state before solving.

## Resource and upstream behavior

The model checks scalar-variable, PSD-dimension, stored-entry, and real-block
allocation limits before construction. The optional extension represents a
complex block of dimension $2d$ as a real block of dimension $4d$, so
`allow_densify=true` is explicit at the backend.

This model is the no-loss Julia-native replacement for the CVX-expression
branch of QETLAB `MatsumotoFidelity.m` at the pinned revision. The numeric
Julia implementation deliberately corrects QETLAB's additive
$10^{-8}I$ regularization and explicit-inverse behavior, which changes
singular-state mathematics. The correction and the solver formulation are
tested separately.
