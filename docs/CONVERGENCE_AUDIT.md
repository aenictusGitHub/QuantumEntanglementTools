# Convergence audit

Evidence date: 2026-07-29.

This audit establishes a bounded starting point for the QETLAB convergence
work. It is an evidence record, not a parity declaration, release approval, or
substitute for the authoritative inventory, provenance, and completion-policy
ledgers.

<!-- qetlab-current-claims: begin -->
## Superseding local completion evidence (2026-07-30)

The later strict static repository audit at pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601` records all 127/127 public rows as
verified with final status and all 36/36 private helpers with terminal
dispositions, with zero queued tasks and zero static-evidence failures. The
strict completion checker passes locally on Julia 1.12.6 and the installed
Julia 1.10.0.

This supersedes only the historical completion-queue conclusions below; their
dated audit counts remain as the baseline evidence actually observed on
2026-07-29. The new result is local static repository evidence, not
MATLAB/QETLAB parity, remote-CI evidence, a performance claim, an API-stability
guarantee, release approval, or the required human review.
<!-- qetlab-current-claims: end -->

## Evidence boundaries

Three states must remain distinct:

| Evidence layer | Observed state | What it supports |
|---|---|---|
| Initial exact baseline | At the start of this audit, local `main`, `origin/main`, and `HEAD` were clean and equal to `2d965bfbbcb6f3af350e293f1034b1fdcd4f937e` | Exact local and remote evidence explicitly tied to that commit |
| Current prospective worktree | Phase 0 corrections and convergence work are uncommitted | Local development only; an exact committed-tree and remote rerun is still required |
| Publication state | The GitHub repository is private, with no tag, GitHub release, or GitHub Pages deployment | No published `v0.1.0` release exists |

Results from the initial baseline do not automatically validate later
uncommitted changes. Likewise, a successful local run does not replace a
successful remote run at the exact commit proposed for release.

## Upstream and ledger baseline

- The exact QETLAB source pin is
  `d8589610f00cff106537268dee2e2a1153f3a601`.
- The generated inventory contains 163 source rows and 503 dependency edges,
  with zero automatically detected cycles.
- The inventory distinguishes 127 public entry points from 36 private helpers.
  Classification is not implementation, and implementation is not verified
  behavioral equivalence.
- The public-API/provenance consistency check covers 214 exported bindings and
  214 corresponding provenance records.
- The differential fixture corpus contains 338 assertions. It is
  function-specific evidence and does not establish general MATLAB
  equivalence.

Per-row status and completion decisions remain owned by the inventory and
completion-policy ledgers. This document deliberately does not promote any row
from partial, deferred, or blocked status.

## Recorded local validation

The following results form the recorded local comparison baseline:

| Check | Recorded result |
|---|---|
| Package tests | 2,417 assertions on Julia 1.12.6 and Julia 1.10.11 |
| Executable tutorials | 48 assertions on Julia 1.12.6 and Julia 1.10.11 |
| EntanglementDetection.jl extension | 125 assertions with the exact 0.2.2 integration |
| Inventory generator | 163 rows, 503 dependency edges, zero detected cycles |
| Public API and provenance | 214 exports and 214 provenance records |
| Independent matrix-predicate validation | 130 assertions |
| Aqua | 11 assertions |
| JET | 25 representative probes |
| Benchmark smoke | 42 cases; not a performance or regression claim |
| Documentation | Strict builds passed on Julia 1.12.6 and Julia 1.10.11 |
| QETLAB differential fixtures | 338 assertions against the exact upstream pin |

The optional extension has an effective Julia 1.11 resolver floor because of
its dependency graph. Its 125-assertion result was recorded on Julia 1.12.6;
it is not a Julia 1.10 extension claim.

Because the current Phase 0 worktree is uncommitted, these counts must be rerun
and attached to the exact resulting commit before they can serve as current
release evidence.

## Exact remote baseline

Remote workflow evidence was inspected at the initial exact commit
`2d965bfbbcb6f3af350e293f1034b1fdcd4f937e`:

- Core passed.
- Documentation passed.
- Quality failed only because the tracked editor-workspace file caused the
  release-integrity check to reject the tree.
- Coverage's package-test steps passed, but the Codecov upload failed with
  `Repository not found`. A passing test step is not evidence of successful
  coverage ingestion.

No exact remote rerun has yet validated the current uncommitted Phase 0
worktree.

## Repository controls and external blockers

- The effective `main` controls observed protect against deletion and
  non-fast-forward updates.
- Required status checks and required pull-request reviews are absent.
- Codecov repository access or configuration must be corrected and a later
  upload must be observed as accepted.
- No release tag, GitHub release, or Pages deployment exists.
- A non-delegable maintainer review remains necessary before publication.

## Source-scope audit

A scoped scan found no occurrence of the excluded-source identifier in the
current tracked tree or current commit archive. A history-only occurrence
remains reachable in prior commits. No history rewrite was performed.

This is a source-scope observation, not a license or secret-scanning guarantee.
Any destructive history rewrite requires separate, explicit authorization and
coordination.

## Evidence required after Phase 0

Before any release claim:

1. Commit the Phase 0 changes as a coherent, reviewable tree.
2. Rerun the applicable local gates on that exact commit and record the commit
   and dirty-worktree state.
3. Push only with explicit authorization, then require successful Core,
   Documentation, Quality, Coverage, and optional-extension evidence at that
   exact SHA.
4. Confirm that Codecov accepted the upload, rather than relying on the test
   step alone.
5. Add required status checks and review requirements if branch policy is to
   enforce those gates.
6. Complete the per-row implementation, verification, and human-review work
   recorded in the authoritative ledgers.
