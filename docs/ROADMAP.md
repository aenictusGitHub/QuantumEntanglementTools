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

Build an optional load-order-safe adapter. Resolve or isolate the pinned
version's global RNG, stdout/logging, and BLAS-thread side effects before calling
the integration safe.

## M7 — QETLAB completeness sweep

Account for every reviewed public upstream function and example with
implementation, compatibility mapping, documented supersession, or explicit
blocker. Complete provenance and validation matrices.

## M8 — release candidate

Review API stability, cross-platform/minimum-version CI, documentation, legal
notices, registry readiness, allocation/performance evidence, and SemVer policy.
Do not tag 1.0 solely because the function count is high.

## Separate QUBIT4MATLAB stage

This stage is blocked before source inspection by a material license conflict.
If authoritative clarification resolves it, work starts from a tested main commit
on an isolated branch, performs a per-file audit, reuses overlapping core
operations, and ends in human technical/legal review rather than an automatic
merge.
