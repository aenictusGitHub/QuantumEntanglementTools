# Porting status

Last updated: 2026-07-29.

This is the human-readable status summary. The generated upstream inventory is
the authoritative function ledger once reviewed. A function is not complete
merely because a similarly named Julia method exists.

## Snapshot

| Area | Evidence-based status |
|---|---|
| Repository/package shell | Present; experimental `v0.1.0` release candidate |
| Neutral name audit | The 2026-07-28 local General snapshot (12,263 non-JLL names) has no exact case-insensitive name or UUID collision, and AutoMerge 1.0.0's similarity check passes; GitHub and JuliaHub checks remain provisional and must be rerun before public registration |
| QETLAB pin | Inspected at `d8589610f00cff106537268dee2e2a1153f3a601` |
| Raw QETLAB source count | 127 root `.m` files and 36 helper `.m` files |
| Authoritative inventory | Generated (163 files, 503 dependency edges, no detected cycles); 76 rows manually reviewed and 87 pending |
| Tier A native kernel | Implemented; 345/345 focused local tests pass on Julia 1.12.6 and Julia 1.10.11, including constructor-integrity and non-one-based-array rejection |
| Tier A QETLAB entry-point mappings | 12 implemented and manually reviewed mappings |
| Tier B operators/states/random | 22 complete mappings plus the scalar bipartite Werner path; 954/954 focused local tests pass, including capped pre-allocation guards and compatibility-keyword forwarding for Brauer-state combinatorics |
| Tier B deferred scope | Multipartite Werner is partial; `RandomSuperoperator` and `RandomPPTState` are explicitly deferred and not exported |
| Tier C channels/maps | 24 native public bindings and 11 `MATLABCompat` wrappers have 35 provenance rows; 177/177 focused local assertions pass, including representation-constructor and post-construction storage-invariant safety |
| Tier C QETLAB entry-point mappings | Five channel/map constructors are implemented; six general-map entry points are partial because two-sided Kraus and/or rectangular row/column operator-space forms remain absent |
| Tier D measures/criteria | 22 native measures/criteria bindings and 11 `MATLABCompat` wrappers are recorded; 168/168 focused local assertions pass |
| Tier D native pipeline | 10 project-native orchestration exports are recorded; 68/68 certificate, tri-state, low-dimensional PPT, conservative pure-Schmidt, and failure-semantics assertions pass |
| Tier D QETLAB entry-point mappings | Ten mappings are implemented; `Entropy` is partial because only the verified von Neumann `ALPHA=1` branch is exposed |
| Tier E coherence | Three reviewed mappings are implemented; 54/54 focused local assertions and 25/25 supplemental fixture assertions pass |
| Tier E product analysis | 10 native bindings and six `MATLABCompat` entry points are recorded; 179/179 native and 52/52 compatibility assertions pass |
| Tier E product mappings | Four mappings are implemented and two are partial; the supplemental artifact passes 68/68 assertions over 14 fixtures |
| Tier E matrix analysis | Four native bindings and four `MATLABCompat` entry points are recorded; 131/131 native and 34/34 compatibility assertions pass |
| Tier E matrix mappings | Four mappings are implemented with native/compatibility distinctions for strong versus weak majorization and compound boundary shapes/errors; the supplemental artifact passes 59/59 assertions over 22 fixtures |
| Tier E matrix predicates | Nine native bindings and four structured-result compatibility entry points pass 170/170 native and 37/37 compatibility assertions; `IsPSD` is partial because the pinned CVX symbolic branch is omitted; no MATLAB-family predicate oracle has been run |
| Recorded full local package run | The integrated 2,417-assertion package corpus (2,369 core plus 48 executable tutorials) passes on Julia 1.12.6 and the minimum supported Julia 1.10.11 |
| Browser code generator | Original BSD-licensed documentation tool covers nine curated state families, six analysis/measure routes, guarded decomposable PPT witnesses, and the five-qubit published symmetric witness; 71 deterministic JavaScript assertions and all nine generated Julia branches pass locally on Julia 1.12.6 and 1.10.11 |
| Reviewed status totals | 63 implemented, 11 partial, 2 deferred, and 87 pending inventory rows |
| Upstream differential validation | Tier A: 13 Octave/QETLAB fixtures and 52 assertions pass; Tier B: 18 deterministic fixtures and 72 assertions pass; Tier C: 7 fixtures and 28 assertions pass; Tier D: 13 fixtures and 34 assertions pass; Tier E coherence: 6 fixtures and 25 assertions pass; Tier E product: 14 fixtures and 68 assertions pass; Tier E matrix: 22 fixtures and 59 assertions pass. MATLAB not run; Octave evidence is function-specific and supplemental |
| Benchmark smoke | 42 quick cases ran locally, including three product-analysis and three matrix-analysis cases; not a regression baseline or comparative performance claim |
| Native entanglement orchestration | Dependency-free certificate-first pipeline implemented; a satisfied necessary criterion remains `unknown` except for the exact `2×2`/`2×3` PPT theorem, and no full `IsSeparable` mapping is claimed |
| EntanglementDetection extension | Exact 0.2.2 weak-dependency adapter implemented with bounded child-process isolation and conservative `unknown` reports; 125/125 focused assertions pass locally on Julia 1.12.6; the six-job Julia 1.11/1.12 Linux/macOS/Windows workflow passed at `6bf8d61`, while an exact release-candidate rerun remains required |
| Prior remote evidence | Core and documentation workflows passed at `485b6a3`; the optional-extension matrix passed at `6bf8d61`; those predecessor runs are not evidence for the unpushed release-candidate tree |
| Exact release-candidate remote evidence | Pending: Quality previously failed because its ignored QETLAB checkout was absent, and the previously green Coverage job hid a rejected upload; both workflows are corrected locally but require successful runs and accepted Codecov ingestion on the exact candidate commit |

## Milestones

| Milestone | Status | Exit evidence still required |
|---|---|---|
| M0 — audit and scaffold | In progress | Review the remaining 87 inventory rows and obtain passing remote CI on the exact candidate commit |
| M1 — subsystem kernel | In progress; local implementation/tests and reviewed mappings pass | Cross-version/platform CI, MATLAB differential review, and reviewed benchmark baseline |
| M2 — states/operators/random | In progress; reviewed core slice locally validated | Add reviewed benchmarks, supported-platform CI, MATLAB evidence, and decide/develop deferred multipartite Werner/random capabilities in their proper milestones |
| M3 — channels/maps | In progress; representation, physicality, constructor, wrapper, and focused oracle tests pass locally | Add two-sided Kraus and rectangular row/column operator-space support or retain explicit partial statuses; add supported-platform CI, MATLAB evidence, and reviewed benchmarks |
| M4 — measures/criteria | In progress; reviewed scalar measures, tri-state necessary criteria, native orchestration, focused tests, supplemental oracle fixtures, and benchmark smoke pass locally | Implement Rényi entropy or retain the explicit partial status; add supported-platform CI and MATLAB evidence; do not claim a general `IsSeparable` implementation |
| M5 — optimization | Not started | Documented formulations, optional solvers, status-aware tests |
| M6 — EntanglementDetection | Local adapter, load-order checks, lifecycle hardening, caller-state isolation, and tri-state translation pass for exact 0.2.2; the six-job remote matrix passed at `6bf8d61` | Rerun the configured Julia 1.11/1.12 Linux/macOS/Windows matrix on the exact candidate; reassess every widened backend version and retain the trusted-worker/non-certificate limitations |
| M7 — QETLAB parity sweep | Not started | Every reviewed inventory row implemented, mapped, superseded, or explicitly blocked |
| M8 — release candidate | In progress; scoped `v0.1.0` release candidate preparation | Complete remote CI, documentation, legal/API/archive review, non-delegable human review of Codex-assisted work, and the explicit visibility/rename decision required for General registration |

The `v0.1.0` candidate is scoped to the package-owned API whose provenance,
documentation, and tests are recorded. It does not complete M7 or change the
fact that 87 QETLAB inventory rows remain pending. A private GitHub release and
General registration have separate gates in `docs/RELEASE_CHECKLIST.md`.

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
6. Rerun the EntanglementDetection Julia 1.11/1.12 platform matrix on the exact
   release candidate and keep its exact-version, child-process, trusted-worker,
   and uncertified-result boundaries explicit.
