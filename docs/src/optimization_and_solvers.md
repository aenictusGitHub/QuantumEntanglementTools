# Optimization and solvers

The core now has a solver-neutral affine SDP representation and stable
status-rich results. JuMP is an optional translation layer; concrete solvers
remain user choices and test-environment dependencies. Completed rows built on
this layer still supply their own formulations, limits, certificate checks,
and evidence; the absolute-PPT ordering-realization model is one such use.

## Mathematical representation

A [`SemidefiniteProgram`](@ref) uses real coordinates $x\in\mathbb{R}^n$,
a scalar affine objective, scalar equalities and intervals, and Hermitian
affine PSD blocks:

```math
\begin{aligned}
\mathop{\mathrm{minimize\ or\ maximize}}\quad&
    c_0+c^\mathsf{T}x,\\
\text{subject to}\quad&
    a_{j0}+a_j^\mathsf{T}x=0,\\
&   \ell_k\leq b_{k0}+b_k^\mathsf{T}x\leq u_k,\\
&   H_{r0}+\sum_i x_iH_{ri}\succeq0.
\end{aligned}
```

The numeric matrices and sparse coefficient vectors are package-owned. No
JuMP expression or optimizer object appears in this model.

Complex Hermitian matrices use one reviewed portable convention. For
$H=A+iB$,

```math
H\succeq0
\quad\Longleftrightarrow\quad
\mathcal{R}(H)=
\begin{bmatrix}A&-B\\B&A\end{bmatrix}\succeq0.
```

The eigenvalues of $\mathcal{R}(H)$ are those of $H$, each repeated twice. If
$Y$ is a real-block dual, the core applies the adjoint embedding

```math
\mathcal{R}^*(Y)
=(Y_{11}+Y_{22})+i(Y_{21}-Y_{12})
```

before exposing a Hermitian dual block. The extension must allocate a dense
matrix of JuMP affine expressions for each real block, so an SDP backend
requires `allow_densify=true`. Model limits are checked first.

## Core building blocks

The solver-neutral layer provides:

- `AffineScalar`, `AffineEquality`, and `AffineInterval`;
- `HermitianAffineMatrix` and `hermitian_variable`;
- `trace_affine`, `partial_trace_affine`, and
  `partial_transpose_affine`;
- `tensor_affine` and `choi_trace_preserving_affine`;
- `hermitian_equalities`;
- `real_block_embedding` and `hermitian_dual_from_real_block`;
- `primal_residual` and `dual_stationarity_residual`;
- `OptimizationLimits`, `SemidefiniteProgram`, `JuMPBackend`, and
  `OptimizationResult`.

Hermitian coordinates are diagonal entries, real upper-triangle entries, then
imaginary upper-triangle entries. `hermitian_variable(:X, d)` therefore uses
$d^2$ real coordinates. Trace and subsystem maps act coefficientwise and
preserve sparse coefficient storage.

## Explicit backend

There is no default optimizer:

```julia
using QuantumEntanglementTools

result = solve_optimization(problem)
@assert result.status ==
        QuantumEntanglementTools.OptimizationBackendUnavailable
```

After adding JuMP and an SDP solver to the active environment, pass their
factory explicitly. Hypatia is illustrative; it is not a core dependency.
From the Julia package prompt opened with `]`, run

```text
add JuMP Hypatia
```

Then construct and inspect the backend explicitly:

```julia
using QuantumEntanglementTools
using JuMP
using Hypatia

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    optimizer_options=(iter_limit=1_000,),
    time_limit_seconds=60,
    allow_densify=true,
    atol=1e-7,
    rtol=1e-7,
)
readiness = backend_status(backend)
@assert readiness.configured

result = solve_optimization(problem, backend)
```

Optimizer options pass through JuMP's public attribute API. A coordinate
`initial_point` in the problem becomes a JuMP warm start. Unsupported options,
malformed factories, missing cones, and backend exceptions become structured
results rather than mathematical negatives.

If this setup does not solve, inspect `backend_status()` and
`backend_status(backend)` before changing the mathematical model. The first
reports whether the JuMP extension is ready for configuration; the second
reports whether the explicit backend is configured. Solver availability is
confirmed only by attempting a solve.

## Result contract

[`OptimizationResult`](@ref) retains:

- package status and the raw termination, primal, and dual statuses;
- primal objective, dual objective, and objective bound when available;
- absolute and relative gaps when available;
- independent primal feasibility and dual stationarity residuals;
- iterations and solve time when exposed by the optimizer;
- primal coordinates, named Hermitian views, scalar duals, real PSD duals, and
  pulled-back Hermitian PSD duals;
- modeling-layer, MOI, configured/reported solver, version, option, tolerance,
  model-limit, embedding, and raw-status metadata;
- warnings and a status-specific explanation.

`OptimizationOptimal` means that the optimizer reported an optimum and the
available package residual checks are consistent at
`atol + rtol * max(1, abs(objective))`. It is still not automatically a
theorem-level certificate. `OptimizationInconsistent` means a nominally
successful solver status conflicts with an independent residual. Inaccurate,
limited, numerical-failure, unsupported, malformed, and missing-backend
outcomes stay distinct.

An `INFEASIBLE` status becomes `OptimizationInfeasible` only with the solver's
infeasibility-certificate result status. A `DUAL_INFEASIBLE` status becomes
`OptimizationUnbounded` only if the primal ray status is present and a supplied
`known_feasible_point` independently passes the package residual check.
Otherwise the result is `OptimizationUnknown`.

## Resource policy

[`OptimizationLimits`](@ref) bounds:

- real scalar variables;
- scalar equalities and intervals;
- PSD block count and original complex dimension;
- stored sparse coefficient entries;
- triangle entries in the doubled real PSD blocks.

Counts use checked `BigInt` arithmetic before JuMP construction. PSD residuals
require a dense Hermitian eigensolve and therefore the same explicit
densification permission. No input is normalized, Hermitian-symmetrized,
projected, clipped, or repaired.

## Validation evidence

The isolated `test/extensions/jump_optimization` environment contains JuMP,
Hypatia, SCS, and no core solver dependency. Its 92 assertions pass on Julia
1.10.11 and Julia 1.12.6. They cover:

- a complex two-by-two PSD optimum with primal and dual residuals;
- a cross-check of that optimum with SCS;
- infeasible and unbounded solver outcomes;
- iteration limit and numerical failure;
- unsupported cone/model operations;
- missing and malformed backends;
- inaccurate feasible and deliberately inconsistent primal evidence;
- sparse/type-preserving affine maps, complex real-block equivalence, Choi
  trace preservation, model validation, and pre-allocation limits.

The tested versions were JuMP 1.31.1, MOI 1.51.2, Hypatia 0.10.0, and SCS
2.6.4. Solver versions and tolerances are runtime evidence and must be recorded
again for each optimization-backed public row.

The absolute-PPT slice adds 23 focused assertions in the same isolated
environment. Hypatia and SCS independently propose a feasible two-qubit
log-gap point; the core then converts that point to exact rationals and
rechecks every strict product-order inequality before returning a negative
absolute-PPT certificate. Missing and malformed backends remain structured
inconclusive outcomes.

The symmetric-extension slice builds outer full-space and bosonic-compressed
models, together with the Navascués--Owari--Plenio inner approximation.
Hypatia tests validate reconstructed primal extensions and normalized dual
separators; SCS independently exercises a feasible outer model. A negative
outer dual is accepted as an entanglement witness only after the complete dual
residual and normalization checks pass. A negative inner dual is deliberately
labeled only as a separator from the selected inner cone. Missing, malformed,
limited, and tolerance-boundary outcomes remain inconclusive. See
[Symmetric-extension hierarchies](symmetric_extensions.md).

The completed model-expression foundations add:

- explicit affine PSD constraint construction for `IsPSD`;
- fixed and composable affine Matsumoto-fidelity block hypographs;
- exact top-`k`, `p` conic epigraph atoms for rectangular complex affine data;
- minimum-error state-discrimination POVM models and residual-checked bounds;
- the symmetric-moment `PolynomialSOS` hierarchy with explicit RNG and sample
  counts for attained inner evidence.

Each has dedicated Julia 1.10/1.12 core tests. Hypatia and SCS exercise the
applicable optimizer paths, including limits and malformed-backend behavior.
See [Positive-semidefinite model constraints](positive_semidefinite_constraints.md),
[Matsumoto-fidelity SDP models](matsumoto_fidelity_model.md),
[Top-k p-norms in optimization models](top_k_p_norm_models.md),
[Minimum-error state discrimination](state_discrimination.md), and
[Polynomial SOS hierarchy](polynomial_sos.md).

## Formulation gate for later rows

Before marking an upstream CVX routine complete:

1. document its primal and, where available, dual;
2. state cones, subsystem order, and complex embedding;
3. return package-owned model and result types;
4. set model-size and work limits before allocation;
5. test small analytic or independent primal/dual instances;
6. cover missing, infeasible, unbounded, inaccurate, limited, unsupported, and
   failed statuses that are meaningful for that formulation;
7. explain whether each reported number is exact, a certified bound, a
   relaxation, a heuristic, or unknown.
