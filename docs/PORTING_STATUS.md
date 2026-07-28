# Porting status

Last updated: 2026-07-28.

This is the human-readable status summary. The generated upstream inventory is
the authoritative function ledger once reviewed. A function is not complete
merely because a similarly named Julia method exists.

## Snapshot

| Area | Evidence-based status |
|---|---|
| Repository/package shell | Present; pre-alpha |
| Neutral name audit | Provisional checks found no exact GitHub or local General registry collision; JuliaHub result unverified |
| QETLAB pin | Inspected at `d8589610f00cff106537268dee2e2a1153f3a601` |
| Raw QETLAB source count | 127 root `.m` files and 36 helper `.m` files |
| Authoritative inventory | Generated (163 files, 503 dependency edges, no detected cycles); 76 rows manually reviewed and 87 pending |
| Tier A native kernel | Implemented; 304/304 focused local tests pass on Julia 1.12.6 |
| Tier A QETLAB entry-point mappings | 12 implemented and manually reviewed mappings |
| Tier B operators/states/random | 22 complete mappings plus the scalar bipartite Werner path; 941/941 focused local tests pass |
| Tier B deferred scope | Multipartite Werner is partial; `RandomSuperoperator` and `RandomPPTState` are explicitly deferred and not exported |
| Tier C channels/maps | 24 native public bindings and 11 `MATLABCompat` wrappers have 35 provenance rows; 154/154 focused local assertions pass |
| Tier C QETLAB entry-point mappings | Five channel/map constructors are implemented; six general-map entry points are partial because two-sided Kraus and/or rectangular row/column operator-space forms remain absent |
| Tier D measures/criteria | 22 native measures/criteria bindings and 11 `MATLABCompat` wrappers are recorded; 162/162 focused local assertions pass |
| Tier D native pipeline | 10 project-native orchestration exports are recorded; 68/68 certificate, tri-state, low-dimensional PPT, conservative pure-Schmidt, and failure-semantics assertions pass |
| Tier D QETLAB entry-point mappings | Ten mappings are implemented; `Entropy` is partial because only the verified von Neumann `ALPHA=1` branch is exposed |
| Tier E coherence | Three reviewed mappings are implemented; 52/52 focused local assertions and 25/25 supplemental fixture assertions pass |
| Tier E product analysis | 10 native bindings and six `MATLABCompat` entry points are recorded; 175/175 native and 52/52 compatibility assertions pass |
| Tier E product mappings | Four mappings are implemented and two are partial; the supplemental artifact passes 68/68 assertions over 14 fixtures |
| Tier E matrix analysis | Four native bindings and four `MATLABCompat` entry points are recorded; 127/127 native and 32/32 compatibility assertions pass |
| Tier E matrix mappings | Four mappings are implemented with native/compatibility distinctions for strong versus weak majorization and compound boundary shapes/errors; the supplemental artifact passes 59/59 assertions over 22 fixtures |
| Tier E matrix predicates | Nine native bindings and four structured-result compatibility entry points pass 166/166 native and 37/37 compatibility assertions; `IsPSD` is partial because the pinned CVX symbolic branch is omitted; no MATLAB-family predicate oracle has been run |
| Recorded full local package run | The integrated 2,270-assertion package corpus passes on Julia 1.12.6 and the minimum supported Julia 1.10.11 |
| Reviewed status totals | 63 implemented, 11 partial, 2 deferred, and 87 pending inventory rows |
| Upstream differential validation | Tier A: 13 Octave/QETLAB fixtures and 52 assertions pass; Tier B: 18 deterministic fixtures and 72 assertions pass; Tier C: 7 fixtures and 28 assertions pass; Tier D: 13 fixtures and 34 assertions pass; Tier E coherence: 6 fixtures and 25 assertions pass; Tier E product: 14 fixtures and 68 assertions pass; Tier E matrix: 22 fixtures and 59 assertions pass. MATLAB not run; Octave evidence is function-specific and supplemental |
| Benchmark smoke | 42 quick cases ran locally, including three product-analysis and three matrix-analysis cases; not a regression baseline or comparative performance claim |
| Native entanglement orchestration | Dependency-free certificate-first pipeline implemented; a satisfied necessary criterion remains `unknown` except for the exact `2×2`/`2×3` PPT theorem, and no full `IsSeparable` mapping is claimed |
| EntanglementDetection extension | Not implemented/safe; global RNG, stdout, logging, and BLAS-thread hazards require isolation tests |
| QUBIT4MATLAB v6.5 | **Blocked** by material archive-license conflict; no implementation `.m` files inspected |

## Milestones

| Milestone | Status | Exit evidence still required |
|---|---|---|
| M0 — audit and scaffold | In progress | Review the remaining 87 inventory rows and obtain passing remote basic CI |
| M1 — subsystem kernel | In progress; local implementation/tests and reviewed mappings pass | Cross-version/platform CI, MATLAB differential review, and reviewed benchmark baseline |
| M2 — states/operators/random | In progress; reviewed core slice locally validated | Add reviewed benchmarks, supported-platform CI, MATLAB evidence, and decide/develop deferred multipartite Werner/random capabilities in their proper milestones |
| M3 — channels/maps | In progress; representation, physicality, constructor, wrapper, and focused oracle tests pass locally | Add two-sided Kraus and rectangular row/column operator-space support or retain explicit partial statuses; add supported-platform CI, MATLAB evidence, and reviewed benchmarks |
| M4 — measures/criteria | In progress; reviewed scalar measures, tri-state necessary criteria, native orchestration, focused tests, supplemental oracle fixtures, and benchmark smoke pass locally | Implement Rényi entropy or retain the explicit partial status; add supported-platform CI and MATLAB evidence; do not claim a general `IsSeparable` implementation |
| M5 — optimization | Not started | Documented formulations, optional solvers, status-aware tests |
| M6 — EntanglementDetection | Native backend interface implemented; optional adapter blocked on safe design | Optional load-order tests and containment of observed global side effects |
| M7 — QETLAB parity sweep | Not started | Every reviewed inventory row implemented, mapped, superseded, or explicitly blocked |
| M8 — release candidate | Not started | Full CI, docs, legal review, benchmarks, API and registry review |
| QUBIT4MATLAB branch | Blocked | Human license clarification before source inspection, then isolated audited branch |

## Completion vocabulary

- `implemented`: code exists; this alone is not a completion claim.
- `verified`: specification, provenance, implementation, applicable tests, and
  documentation pass.
- `compatibility_alias`: tested wrapper delegates to a verified native API.
- `superseded_with_documented_mapping`: no capability loss and migration is
  documented.
- `blocked_with_explicit_reason`: evidence and required human/external action are
  recorded.

`forgotten`, blank, or an unreviewed generated classification is never an
acceptable final status.

## Immediate gates

1. Review the remaining 87 generated root/helper rows and dependency edges.
2. Keep `UpstreamManifest.toml`, `PROVENANCE.toml`, and reviewed status overlays
   synchronized as later milestones export new bindings.
3. Run the Tier A--Tier E reviewed slices and native pipeline on Julia 1.10 and
   the supported-platform CI matrix.
4. Add MATLAB differential evidence where MATLAB becomes available; keep
   analytic/property/independent formulations as primary checks.
5. Turn the current reviewed benchmark smoke cases into a regression baseline
   without claiming comparative speedups.
6. Keep optional integrations optional and preserve tri-state certification
   semantics.
7. Resolve the QUBIT4MATLAB license conflict before any implementation source
   inspection.
