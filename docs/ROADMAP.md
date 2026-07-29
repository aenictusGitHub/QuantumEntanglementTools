# Roadmap

The roadmap is gated by evidence rather than dates. Detailed current status is in
`PORTING_STATUS.md`.

## M0 — audit and scaffold

Pin sources and license hashes; review the function inventory; establish
provenance, conventions, package/docs/test environments, basic CI, and
maintainer handoff records.

## M1 — subsystem kernel

Verify tensor products, permutations, partial trace, partial transpose,
realignment, plans, sparse behavior, generic types, properties, independent
formulations, and representative benchmarks.

## M2 — states, operators, and random constructors

Complete the reviewed inventory slice with explicit RNG APIs, named-family
normalizations, analytic fixtures, docs, and benchmarks.

## M3 — channels and maps

Freeze Choi/Kraus/superoperator conventions through round trips and physicality
tests, then implement the mapped channel and positive-map inventory.

## M4 — measures and non-optimization criteria

Validate norms, entropies, fidelity conventions, entanglement/coherence
measures, and certificate versus necessary-condition semantics.

## M5 — optimization

Accept an architecture ADR only after representative modeling prototypes.
Implement documented primal/dual formulations through optional, status-aware
solver extensions.

## M6 — EntanglementDetection integration

Maintain the exact-version, load-order-safe adapter through a fresh child
process. The local 0.2.2 gate now covers caller RNG, stdout/logger, BLAS state,
timeouts, forced termination, interrupts, malformed responses, and conservative
candidate translation. Completion still requires the configured Julia
1.11/1.12 Linux/macOS/Windows workflow to pass remotely; every future backend
version must be re-audited before widening compatibility.

## M7 — QETLAB completeness sweep

Account for every reviewed public upstream function and example with
implementation, compatibility mapping, documented supersession, or explicit
blocker. Complete provenance and validation matrices.

## M8 — release candidate

Review API stability, cross-platform/minimum-version CI, documentation, legal
notices, registry readiness, allocation/performance evidence, and SemVer policy.
Do not tag 1.0 solely because the function count is high.
