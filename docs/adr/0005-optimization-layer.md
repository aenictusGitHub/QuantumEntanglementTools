# ADR 0005: package-owned conic models with an optional JuMP/MOI translator

- Status: Accepted
- Date: 2026-07-30

## Context

Several pinned QETLAB capabilities use CVX. The Julia package needs portable
SDP modeling without making any optimizer, modeling layer, or global solver
state part of the dependency-free core. Solver objects also cannot define the
long-lived public result contract: their status enums, raw objects, and
available attributes differ.

The decision gate in the proposed version of this ADR was exercised with JuMP
1.31.1, MathOptInterface (MOI) 1.51.2, Hypatia 0.10.0, and SCS 2.6.4 on Julia
1.10.11 and Julia 1.12.6. Hypatia was the primary SDP test backend and SCS
independently cross-checked the representative complex SDP. This is
implementation evidence, not a promise that these exact patch versions are the
only compatible versions.

## Alternatives

| Alternative | Maintenance and portability | Complex PSD and duals | Cost |
|---|---|---|---|
| Direct MOI construction | Precise but repeats low-level function/set/index plumbing in every formulation | Full access, but manual bridges, result indices, and matrix-cone scalar products are easy to get wrong | Lowest modeling dependency; highest package maintenance burden |
| JuMP plus MOI | Stable public model/status/dual APIs and broad optimizer portability | JuMP has Hermitian cones and public dual extraction; explicit real blocks work with ordinary real PSD solvers | Optional compile/load cost; macros and model objects must remain inside the extension |
| Only a project-owned conic layer | Best stable core boundary and easiest deterministic preflight | Cannot solve by itself; still needs a maintained translator and status adapter | More package data types, but no solver leakage |

Direct MOI gave no demonstrated benefit large enough to justify duplicating
JuMP's affine-expression, bridge, warm-start, and result-query machinery.
Depending only on JuMP would instead leak backend expressions into native
problem specifications.

## Decision

Use a hybrid of the latter two alternatives:

1. The core owns immutable scalar-affine, Hermitian-affine, SDP, limits,
   backend-configuration, metadata, primal/dual, status, and result types.
   Native formulations contain numeric data, never JuMP or MOI expressions.
2. JuMP is a weak dependency loaded through
   `QuantumEntanglementToolsJuMPExt`. The extension translates the
   package-owned model using documented JuMP and MOI APIs.
3. A `JuMPBackend` always carries an explicit optimizer factory and options.
   There is no global or package default optimizer. Concrete solvers remain in
   dedicated test or example environments.
4. Complex Hermitian PSD constraints use the centralized real block
   embedding

   ```math
   H=A+iB \succeq 0
   \quad\Longleftrightarrow\quad
   \begin{bmatrix}A&-B\\B&A\end{bmatrix}\succeq0.
   ```

   The extension does not rely on solver-native complex cones. This doubles
   each PSD dimension, so the caller must opt into densification and the core
   checks variable, constraint, PSD-dimension, coefficient-entry, and
   real-block scalarization limits before model construction.
5. The core centralizes Hermitian coordinates, trace, partial trace, partial
   transpose, Choi trace-preservation, tensor maps, real-block dual pullback,
   and primal/dual residual calculations.
6. Solver results retain termination, primal and dual status; primal, dual, and
   bound objectives when available; absolute and relative gaps; residuals;
   iterations; solve time; raw status; layer/interface/solver versions; options;
   and reconstructed witnesses. A solver optimum is not automatically a
   mathematical certificate.
7. An infeasibility or unbounded classification requires the corresponding
   solver certificate status. A dual-infeasible result is classified as
   unbounded only when a package-owned feasible point is supplied and passes
   the independent residual check. Inaccurate, limited, inconsistent,
   unsupported, malformed, failed, and missing-backend paths remain explicit.

## Consequences

- Core loading and solver-free APIs remain independent of JuMP and all solvers.
- New optimization-backed QETLAB rows must first define package-owned numeric
  model data, then add a thin extension translation or compose the accepted SDP
  representation. They may not return raw JuMP expressions.
- Model authors can use coordinate warm starts without backend coupling.
- Portable real blocks favor auditability over the smaller native Hermitian
  cone representation. A future accepted ADR may add a proven equivalent
  native-complex translation, but it must preserve residual and dual semantics.
- The infrastructure alone completes no QETLAB solver-backed row. Each row
  still needs its formulation, provenance, tests, documentation, and
  certificate interpretation.
