# Entanglement backends

No optional entanglement backend is currently declared production-ready.

## Result semantics

High-level analysis uses package-owned `EntanglementReport` values, with one
`EntanglementAttempt` per method run. A report distinguishes:

- `entangled`: backed by the stated witness or sufficient criterion;
- `separable`: backed by a valid sufficient criterion or decomposition;
- `unknown`: no requested certificate was obtained.

A failed necessary test is generally not a separability certificate. An
exception, time limit, numerical failure, or undocumented backend boolean must
not be translated into `separable`.

Every report records whether the conclusion is certified, its named
`certificate_kind`, method/backend, package-owned evidence, ordered attempt
history, and a message. Necessary-criterion details use `CriterionResult`,
whose `CriterionStatus` is `CriterionEntanglementDetected`,
`CriterionSatisfied`, or `CriterionUnknown`, with the measured value, exact
boundary, numerical tolerance, and optional witness. A satisfied necessary
condition is not silently relabeled `separable`.

## Native pipeline

The dependency-free backend is represented by `NativeEntanglementBackend`;
`backend_capabilities` reports its package-owned metadata, and
`available_entanglement_backends()` discovers it without loading optional
packages. `AbstractEntanglementBackend` and `AbstractEntanglementMethod` are
dispatch interfaces for future extensions. `NativePPT` is the validated
configuration for the currently supported explicit
`detect_entanglement(state, dims, method)` call.

For a density matrix, `analyze_entanglement(...;
strategy=:certificates_first)` runs:

1. `ppt_criterion`;
2. `realignment_criterion`;
3. `reduction_criterion` for bipartite layouts.

Every attempted method is retained. A robust violation stops with a named
entanglement certificate. A PPT pass certifies separability only in the exact
bipartite `2×2` or `2×3` domain, where PPT is sufficient. In higher dimensions
a PPT pass is `unknown`; realignment and reduction passes are also `unknown`.
Tolerance-boundary values remain `unknown`. For bipartite pure vectors, a
trailing Schmidt coefficient above the configured threshold certifies
entanglement. Separability is certified only when the computed trailing
coefficients are exactly zero; a tolerance-defined rank-one result with
nonzero trailing coefficients remains `unknown`.

These APIs are intentionally narrower than a general separability solver. They
do not export or claim a complete replacement for QETLAB `IsSeparable`.
`MATLABCompat.IsPPT` is also intentionally tri-state: it accepts a finite
Hermitian operator without requiring unit trace but returns `CriterionResult`
rather than collapsing a numerical-boundary case to a Boolean.

The native pipeline passes 68 focused assertions. The underlying scalar
measures and criteria, including all 11 Tier D compatibility wrappers, pass
162 focused assertions. The Tier D supplemental Octave/QETLAB artifact has 13
fixtures and 34 passing comparisons; its SHA-256 is
`ad0cdc45077390fc1eb736fc7c7ff1ec41696c796a508b536774cb6e0020160a`.
This evidence does not establish complete QETLAB parity or supported-platform
coverage.

## EntanglementDetection.jl 0.2.2

The pinned audit target is version 0.2.2 at
`5f60da1ceef6442acb669e10acc2fa47670bab06`. It must remain a weak dependency
loaded through a Julia extension.

The upstream audit found integration hazards:

- `separable_distance` calls `Random.seed!(0)`;
- a supplied logfile can redirect global stdout, and the routine logs/flushes;
- `AlternatingSeparableLMO(..., parallelism=true)` globally sets BLAS threads to
  one.

Until these are isolated without monkey-patching private APIs, and regression
tests prove the caller's RNG stream, stdout, logging, and BLAS configuration are
unchanged, the extension must not be described as safe or complete.

Required fresh-process tests cover: core package alone, both load orders,
precompilation/repeated loading, absence of method ambiguities or type piracy,
backend exceptions, and global-state preservation.
