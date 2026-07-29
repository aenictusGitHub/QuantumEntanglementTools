# Entanglement backends

The native certificate pipeline is dependency-free. The optional
EntanglementDetection.jl adapter is available for exact version 0.2.2 as a
locally validated heuristic candidate generator; it is not a source of
package-certified conclusions.

## Result semantics

High-level analysis uses package-owned `EntanglementReport` values, with one
`EntanglementAttempt` per method run. A report distinguishes:

- `entangled`: backed by the stated witness or sufficient criterion;
- `separable`: backed by a valid sufficient criterion or decomposition;
- `unknown`: no requested certificate was obtained.

Passing a necessary condition is generally not a separability certificate. An
exception, time limit, numerical failure, or undocumented backend boolean must
not be translated into `separable`.

Every report records whether the conclusion is certified, its named
`certificate_kind`, method/backend, package-owned evidence, ordered attempt
history, and a message. Necessary-criterion details use `CriterionResult`,
whose `CriterionStatus` is `CriterionEntanglementDetected`,
`CriterionSatisfied`, or `CriterionUnknown`, with the measured value, exact
boundary, numerical tolerance, and optional witness. A satisfied necessary
condition is not silently relabeled `separable`.

For practical constructions and complete output examples, start with
[Separability by example](separability_examples.md). There is intentionally no
general Boolean `is_separable`; use the structured status and certification
fields instead.

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

The pipeline does not call `in_separable_ball`. That function is a separate
sufficient test with its own `SeparableBallResult` statuses; invoke it
explicitly when a density matrix may be close enough to the maximally mixed
state.

These APIs are intentionally narrower than a general separability solver. They
do not export or claim a complete replacement for QETLAB `IsSeparable`.
`MATLABCompat.IsPPT` is also intentionally tri-state: it accepts a finite
Hermitian operator without requiring unit trace but returns `CriterionResult`
rather than collapsing a numerical-boundary case to a Boolean.

The native pipeline passes 68 focused assertions. The underlying scalar
measures and criteria, including all 11 Tier D compatibility wrappers, pass
168 focused assertions. The Tier D supplemental Octave/QETLAB artifact has 13
fixtures and 34 passing comparisons; its SHA-256 is
`ad0cdc45077390fc1eb736fc7c7ff1ec41696c796a508b536774cb6e0020160a`.
This evidence does not establish complete QETLAB parity or supported-platform
coverage.

## EntanglementDetection.jl 0.2.2

The pinned audit target is version 0.2.2 at
`5f60da1ceef6442acb669e10acc2fa47670bab06`. It is an exact-version weak
dependency loaded through a Julia extension; core loading and `Pkg.test()` do
not install it. The runtime checks the package version but not the source-tree
hash, so controlled test and release environments must resolve the registered
release or the recorded checkout.

The upstream audit found integration hazards:

- `separable_distance` calls `Random.seed!(0)`;
- a supplied logfile can redirect global stdout, and the routine logs/flushes;
- `AlternatingSeparableLMO(..., parallelism=true)` globally sets BLAS threads to
  one.

`EntanglementDetectionSearch` contains those effects by launching the documented
public backend call in a bounded child Julia process. It exposes neither logfile
nor parallelism controls, preserves the caller's RNG/stdout/logger/BLAS state,
and converts every backend `true`, `false`, or `nothing` value into package-owned
candidate evidence inside an uncertified `unknown` report. Timeouts, backend
exceptions, process failures, and invalid structured responses are also
`unknown`, never mathematical negatives.

The dedicated suite passes 125/125 assertions locally on Julia 1.12.6,
including both load orders, method ambiguities, a live search, explicit
real-to-complex representation conversion, caller-state preservation, bounded
TERM-to-KILL escalation, interrupt cleanup, response validation, and
output/read limits. The core pipeline suite separately checks dependency
absence. EntanglementDetection 0.2.2 currently resolves only on Julia 1.11 or
later because of Ket 0.9 registry compatibility. The six-job Julia 1.11/1.12
Linux/macOS/Windows workflow passed at predecessor commit `6bf8d61`; an exact
release-candidate rerun remains required.

See [EntanglementDetection.jl extension](entanglement_detection_extension.md)
for installation, execution, failure semantics, IPC trust boundaries, and
remaining limitations.
