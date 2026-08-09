# Porting status

Last updated: 2026-08-09.

This is the human-readable status summary. The generated upstream inventory is
the authoritative function ledger once reviewed. A function is not complete
merely because a similarly named Julia method exists.

## Snapshot

<!-- qetlab-current-claims: begin -->

| Area | Evidence-based status |
|---|---|
| Repository state | `main` and `origin/main` are at `c69185f46c7027906e07d4965ed58dbfeab6457f`; five GitHub Pages implementation/evidence files, two status records, and the reconciled generated snapshot are uncommitted locally. No commit, push, tag, release, visibility change, or branch-setting change was made by the coding agent; the repository owner enabled GitHub Actions as the Pages source |
| QETLAB source | Pinned and clean at `d8589610f00cff106537268dee2e2a1153f3a601`; 163 MATLAB files, 127 public functions, 36 private helpers, 503 dependency edges, and no detected cycle |
| Strict completion ledger | 127/127 public rows are verified with final status; completion queue contains 0 public rows; 36/36 internal helpers have terminal dispositions; 0 required helpers remain; 0 static completion failures |
| Public API/provenance | 477 runtime exports (347 native/module and 130 `MATLABCompat`) match 477 provenance records |
| Full package corpus | 9,484/9,484 assertions—9,400 core plus 84 executable-tutorial assertions—pass on Julia 1.12.6 and the installed Julia 1.10.11 |
| Optional optimization | The complete JuMP/Hypatia/SCS environment passes 847/847 assertions on both installed Julia lines; solver output remains status-rich and is not automatically a certificate |
| EntanglementDetection.jl | Exact 0.2.2 adapter passes 141/141 assertions on Julia 1.12.6; its effective resolver floor remains Julia 1.11 because of Ket 0.9 |
| User-facing diagnostics | Additive result interpretation, density validation, strategy/backend discovery, symmetric multiqubit/multiqudit coordinates, and task-first documentation preserve existing function names and certificate boundaries; the Pages source is enabled, but the first successful post-fix deployment remains pending |
| Interactive generator | Schema 2 provides 11 bounded families, nine curated presets, five additional analysis routes, four fixed core separability profiles, resource planning, and strict versioned JSON portability; 178 JavaScriptCore checks pass, and the 13-module generated Julia bundle completes on Julia 1.12.6 and 1.10.11 |
| Independent predicates | Randomized spectrum/minor validation passes 130/130 assertions on both installed Julia lines |
| Executable tutorials | Seven standalone scripts pass 84/84 assertions on both installed Julia lines; 36 assertions cover the new seeded Schmidt and exact Tiles-UPB workflows |
| Quick benchmarks | All 114 declared quick benchmark cases completed without failure on the current uncommitted worktree; this is local smoke evidence only, and targeted paired minima remain diagnostics rather than a stable comparative baseline |
| Remote/release evidence | Dirty-worktree preflight and isolated archive smoke pass on both installed Julia lines; exact-commit supported-platform CI, Codecov ingestion, maintainer review, and publication decisions remain open |

The strict result is a bounded local repository-evidence claim. It does not
establish complete MATLAB/QETLAB parity, supported-platform CI, comparative
performance, API stability, release approval, or the maintainer's
non-delegable review.

<!-- qetlab-current-claims: end -->

The usability layer is additive: existing status fields and computational entry
points remain available. The new helpers explain rather than reclassify
results, density validation never repairs input, and readiness diagnostics do
not load optional packages or turn solver availability into a mathematical
conclusion. The browser generator uses those same conservative semantics:
necessary tests remain one-sided, `unknown` remains explicit, optional
backends are only inspected, and generated separability searches use fixed
dependency-free strategy tuples and resource limits.

## Milestones

| Milestone | Status | Exit evidence still required |
|---|---|---|
| M0--M6 — architecture and work packages | Local implementation and focused evidence complete | Exact-tree remote/platform evidence and human review remain release gates |
| M7 — QETLAB completeness sweep | Local static objective reached: 127/127 public rows and 36/36 helpers are terminal | MATLAB equivalence is not claimed; keep source-free fixtures supplemental to analytic/property evidence |
| M8 — release convergence | In progress; result ergonomics, validation/discovery helpers, schema-2 generator, task-first docs, the Pages workflow, and local exact-tree gates are integrated; the Pages source is enabled | A successful exact-commit Pages deployment, exact committed-tree CI, legal/API review, and non-delegable maintainer review |

## Generator capability evidence (2026-07-31)

The browser-local generator now supports 11 curated families, including an
exact `BigInt`/`Rational{BigInt}` Tiles-UPB complement and locally seeded
Hilbert--Schmidt or Bures random density matrices. Five new routes expose
non-mutating validation, marginals, pure-state Schmidt data, bounded
separability reports, and backend readiness. Four separability profiles expand
only to explicit dependency-free strategies; no free-form Julia identifier,
optional solver, or implicit random source can enter generated code.

Schema-1 configurations migrate with every new analysis and output option
disabled. Schema 2 rejects unknown keys, non-finite or out-of-range values,
future versions, and JSON over 16 KiB measured as UTF-8. The UI adds a live
resource estimate, accessible field errors, stale-output protection, and local
copy/download/load controls without browser persistence or network transfer.

The JavaScriptCore/JXA harness passes 178 checks. Its 13 generated Julia modules
contain 130 structural and runtime assertions and complete on Julia 1.12.6 and
Julia 1.10.11. Strict Documenter builds also pass on both lines; only the
previously recorded HTML/search-index size warnings remain.

## Completion vocabulary

- `verified`: final public-row status after the static ledger finds synchronized
  native and compatibility mappings, provenance, implementation, tests, and
  documentation.
- `compatibility_alias`: tested wrapper delegates to a verified native API.
- `superseded_with_documented_mapping`: no capability loss and migration is
  documented.
- Private-helper dispositions record whether a helper was replaced internally,
  intentionally excluded, or subsumed by a completed parent. They are not
  public API promises.

The strict checker is deliberately static. Passing it does not prove semantic
correctness, MATLAB parity, or solver-certificate validity.

## Symmetric multiqubit and multiqudit tools (2026-08-09)

The project-native API now includes ten occupation-coordinate operations for
symmetric states. Exact dimension/rank arithmetic, multiqudit Dicke
construction, compressed product coordinates, collective operators, split
isometries, reductions, and coordinate/ambient maximally mixed states are
integrated with provenance, conceptual/API documentation, and resource guards.
The focused suite passes 1,014/1,014 on Julia 1.12.6 and 1.10.11. This addition
does not change the pinned QETLAB inventory or claim a QETLAB compatibility
mapping for these project-native functions.

## GitHub Pages deployment repair (2026-08-09)

The project and code-generator URLs returned 404 because Pages was not yet
provisioned and the latest Documentation workflow failed before upload and
deployment. The owner has now selected GitHub Actions as the Pages source. The
local repair removes the platform-sensitive numerical dependency from the live
Tiles example by retaining the exact boundary-state proof while applying the
independent realignment calculation to a fixed, exactly defined full-rank
depolarized neighbor. It also corrects the documented workflow artifact name.

The focused example and all 84 tutorial assertions pass on Julia 1.12.6 and
1.10.11. The 178-check JavaScriptCore/JXA generator harness, its generated
13-module Julia smoke on both Julia lines, strict CI-shaped Documenter builds
on both lines, the formatter gate, and `git diff --check` pass. The remote
branch remains at `c69185f46c7027906e07d4965ed58dbfeab6457f`; the five
implementation/evidence files, two status records, and reconciled generated
snapshot are uncommitted, so the live site is expected to remain unavailable
until an authorized push triggers a successful Documentation deployment.

## Remaining gates

The formatter, local documentation, benchmark smoke, claim reconciliation,
two-version runtime checks, dirty-worktree distribution preflight, and isolated
archive smoke have been completed for this worktree.

1. Have the maintainer perform the non-delegable mathematical, API, provenance,
   licensing, and generated-change review.
2. Only with explicit authorization, commit and push the Pages repair, require
   the exact-SHA Documentation deployment and supported-platform CI, verify the
   live root and `/code_generator/` routes, and confirm Codecov ingestion.
3. Treat tagging, a private GitHub release, public visibility, and General
   registration as separate, explicitly authorized decisions.

## Twirl completion evidence (2026-07-30)

The pinned `Twirl` row is implemented for Werner, isotropic, real, and Pauli
branches with strict copy/dimension validation, sparse spanning operators,
exact Gram-dependence elimination, and explicit combinatorial, dense-entry,
nonzero, and work guards. The focused suite passes 149/149 on Julia 1.12.6 and
1.10.11. Six source-free pinned Octave/QETLAB fixtures pass 18/18 native,
compatibility, and idempotence comparisons; MATLAB was not run.

## Absolute-PPT completion evidence (2026-07-30)

The pinned `AbsPPTConstraints` and `IsAbsPPT` rows are implemented through
package-owned numeric and affine-LMI representations, explicit dimensions,
typed sparse opt-in, bounded monotone/criss-cross ordering enumeration, and
certificate-aware tri-state results. The native API distinguishes analytic,
sufficient, exhaustive, negative-certificate, capped, boundary, and backend
outcomes. A negative result is issued only after exact rational realization of
the candidate ordering; optional solver work remains isolated in the JuMP
extension.

On both Julia 1.12.6 and Julia 1.10.11, the focused suites pass 197/197 native,
24/24 compatibility, 28/28 committed Octave-fixture, and 23/23 dedicated JuMP
extension assertions. The complete JuMP-extension runner passes 115/115 on
both Julia lines. The fixture digest is
`17b825d6e15fae1f2391162137255603653de245fea39661c18ba1db65e74f6c`;
MATLAB was not run.

After regeneration, the static ledgers pass on both Julia lines with 316 public
bindings, 316 provenance entries, 163 reviewed upstream rows, 503 dependency
edges, zero cycles, and zero completion-checker failures. The completion view
contains 127 public rows: 95 implementation-complete, five partial, one
deferred, and 26 blocked, with 32 queued rows. The documentation-math scan
passes over 45 Markdown files, 210 inline spans, and 99 display blocks, and the
79-case quick benchmark smoke includes the two absolute-PPT cases. These are
local implementation and validation claims, not a general QETLAB parity or
performance claim.

## Symmetric-extension and random-PPT completion evidence (2026-07-30)

The pinned `SymmetricExtension` and `SymmetricInnerExtension` rows now use
package-owned solver-neutral SDP models with an explicit optional backend.
Outer models support full permutation-invariant and bosonic-compressed
extensions, representative PPT constraints, solver-free theorem branches, and
independently validated primal or dual certificates. Inner models implement
the Navascués--Owari--Plenio transform and preserve the essential warning that
a negative inner-cone dual is not automatically an entanglement witness. The
private exact Jacobi recurrence supersedes the pinned `jacobi_poly` helper
without exporting it.

`RandomPPTState` is now a bounded, mandatory-explicit-RNG construction. A
shifted induced state covers full-rank requests and a convex mixture of random
product projectors covers low-rank requests. A candidate becomes `state` only
after independent Hermiticity, PSD, trace-one, PPT, and rank-bound checks;
resource and zero-iteration exits occur before consuming the RNG.

On Julia 1.12.6 and Julia 1.10.11, the focused suites pass 162/162 native,
31/31 compatibility, 85/85 dedicated Hypatia/SCS extension, and 31/31
committed Octave-fixture assertions. The source-free fixture covers five
solver-free extension decisions, three Jacobi coefficient vectors, and
random-PPT properties, with SHA-256
`572dc3d846d1adb0b2ee64f9ee5a3b7a6903eb86d56927ff76e5e8fa520551f0`;
MATLAB and CVX were not run.

Regenerated ledgers pass with 331 public bindings and provenance entries, 163
reviewed upstream rows, 503 dependency edges, zero cycles, 98
implementation-complete public rows, 29 queued public rows, 25 terminal helper
dispositions, and zero static-evidence failures. The 82-case quick benchmark
smoke passes, including analytic outer-extension, bosonic PPT model-building,
and random-PPT cases. The documentation-math scan passes 48 Markdown files,
244 inline spans, and 109 display blocks. The repository-wide Documenter run
was blocked only by the concurrently pending state-discrimination module
integration; no symmetric-extension or random-PPT documentation failure was
reported.
