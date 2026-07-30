# Optimization architecture

Optimization-backed functionality is split into a dependency-free numeric
model, an optional translator, and a user-selected solver. The split keeps
solver objects and modeling expressions out of the package's public problem
and result contracts.

```text
package-owned affine SDP
          |
          | explicit JuMPBackend
          v
optional JuMP/MOI extension
          |
          | explicit optimizer factory
          v
user-installed conic solver
```

The accepted design is recorded in
[ADR 0005](https://github.com/aenictusGitHub/QuantumEntanglementTools/blob/main/docs/adr/0005-optimization-layer.md).
The user-facing model, backend, result, and resource contracts are described in
[Optimization and solvers](optimization_and_solvers.md).

## Layer responsibilities

The core owns:

- scalar-affine and Hermitian-affine numeric data;
- real coordinate conventions and subsystem-aware linear maps;
- SDP validation and pre-allocation limits;
- backend-independent statuses, metadata, primal and dual evidence;
- independent primal and dual residual calculations.

The optional `QuantumEntanglementToolsJuMPExt` extension owns:

- translation to documented JuMP and MathOptInterface APIs;
- the explicit real-block representation of complex Hermitian PSD cones;
- optimizer option, warm-start, and time-limit forwarding;
- extraction and normalization of solver statuses and evidence.

The caller owns the concrete optimizer choice. `JuMPBackend` requires an
optimizer factory; the package has no global or default solver.

## Complex Hermitian cones

For a Hermitian matrix $H=A+iB$, the extension constructs the real symmetric
block

```math
\mathcal{R}(H)=
\begin{bmatrix}
A&-B\\
B&A
\end{bmatrix}.
```

The constraints $H\succeq0$ and $\mathcal{R}(H)\succeq0$ are equivalent. This
doubles each PSD block dimension and produces dense JuMP affine-expression
matrices. The caller must therefore set `allow_densify=true`, and the core
checks all configured model limits before the extension allocates the block.

## Status and evidence boundary

Raw solver outcomes never become Boolean mathematical conclusions. An optimal
solver status must agree with the available package residual checks. An
infeasible result requires an infeasibility-certificate status. An unbounded
result requires a primal-ray status and an independently checked feasible
anchor. Limited, inaccurate, inconsistent, unsupported, malformed, and failed
outcomes remain distinct.

`OptimizationResult` records the raw and package statuses, objectives and
bounds, residuals, iterations, solve time, reconstructed primal and dual
objects, solver versions and options, tolerances, model limits, warnings, and
an explanation. This evidence is not automatically a theorem-level
certificate.

## Completion boundary

This architecture completes no solver-backed QETLAB row by itself. Each
upstream routine still needs a reviewed formulation, row-specific provenance,
analytic or independent validation, solver-status tests, documentation, and a
precise certificate interpretation before its ledger status can change.
