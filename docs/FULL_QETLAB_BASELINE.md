# Full-QETLAB baseline

Evidence date: 2026-07-29.

This is the reproducible starting point for the program to cover all 127
public rows at the pinned QETLAB revision. It is not a completeness or parity
claim. The generated per-row plan is
[`QETLAB_COMPLETION_PLAN.md`](QETLAB_COMPLETION_PLAN.md), and the broader
governance evidence is in [`CONVERGENCE_AUDIT.md`](CONVERGENCE_AUDIT.md).

<!-- qetlab-current-claims: begin -->
## Superseding local completion evidence (2026-07-30)

The later strict static repository audit at pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601` records all 127/127 public rows as
verified with final status and all 36/36 private helpers with terminal
dispositions, with zero queued tasks and zero static-evidence failures. The
strict completion checker passes locally on Julia 1.12.6 and the installed
Julia 1.10.0.

This supersedes only the historical completion-queue conclusions below; their
dated counts remain as the baseline evidence actually observed on 2026-07-29.
The new result is local static repository evidence, not MATLAB/QETLAB parity,
remote-CI evidence, a performance claim, an API-stability guarantee, release
approval, or the required human review.
<!-- qetlab-current-claims: end -->

## Repository states

The audit began from clean local and remote `main` at
`2d965bfbbcb6f3af350e293f1034b1fdcd4f937e`. During the audit, the maintainer
committed and pushed the bounded governance corrections as
`ec9094dad43a7531b16b1f1d282a490ccba0c543` (`Baseline corrections`). The
current convergence implementation is an uncommitted worktree on top of that
maintainer commit. No agent commit, tag, release, history rewrite, Pages
deployment, or repository-visibility change was made.

Results tied to the initial exact commit and results from the prospective
worktree are deliberately kept separate. Before release, every applicable
gate must be rerun on the exact committed tree that is proposed for
publication.

## Inventory and implementation snapshot

| Measure | Count |
|---|---:|
| Pinned QETLAB source rows | 163 |
| Public rows | 127 |
| Internal helpers | 36 |
| Static source-call edges | 503 |
| Automatically detected source-call cycles | 0 |
| Public rows marked implemented | 63 |
| Public rows marked partial | 15 |
| Public rows marked deferred | 19 |
| Public rows blocked with an explicit reason | 30 |
| Public rows still requiring completion | 64 |
| Exported package and compatibility bindings | 214 |
| Public provenance records | 214 |

The 163-row source-disposition classification is complete. Per-function
implementation, semantic validation, MATLAB/QETLAB differential evidence, and
maintainer release review are not complete. A generated non-strict checker
reports the current queue; its strict mode is expected to fail until all 127
public rows and their required private helpers satisfy the completion
contract.

## Reproduced local gates

These commands were run from the repository root during baseline reproduction.
Counts are observations from the run, not copied expectations.

| Command or gate | Julia | Exit/result |
|---|---|---|
| `Pkg.instantiate(); Pkg.precompile(); Pkg.test()` | 1.12.6 | Passed, 2,417 assertions |
| `Pkg.instantiate(); Pkg.test()` | 1.10.11 | Passed, 2,417 assertions |
| `tutorials/runtests.jl` | 1.12.6 and 1.10.11 | Passed, 48 assertions on each |
| `scripts/build_upstream_inventory.jl --check` | 1.12.6 and 1.10.11 | Passed, 163 rows, 503 edges, zero cycles |
| `scripts/check_public_api.jl` | 1.12.6 and 1.10.11 | Passed, 214 exports and 214 provenance rows |
| `scripts/validate_matrix_predicates.jl` | 1.12.6 | Passed, 130 assertions |
| `quality/run_quality.jl` | 1.12.6 | Aqua 11/11 and JET 25/25 passed |
| `benchmark/benchmarks.jl --quick --no-save` | 1.12.6 | All 42 smoke cases ran |
| `docs/make.jl` | 1.12.6 and 1.10.11 | Strict builds and doctests passed |
| Optional EntanglementDetection suite | 1.12.6 | Passed, 125 assertions with exact version 0.2.2 |
| Seven source-free QETLAB fixture comparisons | 1.12.6 | Passed, 338 assertions in total |

The fixture total comprises Tier A 52, Tier B 72, Tier C 28, Tier D 34,
coherence 25, product analysis 68, and matrix analysis 59 assertions. Octave
was available; MATLAB was not. These fixtures are supplemental,
function-specific evidence and cannot be generalized to untested branches.

The quick benchmark is a smoke run, not a performance comparison or regression
baseline. The optional extension resolves only on Julia 1.11 or later because
of its Ket 0.9 dependency and remains outside the core test target.

## Baseline failures and environment boundaries

- The release-integrity gate at the initial exact commit rejected an
  editor-workspace path. The maintainer's later `ec9094d` commit removes that
  path and adjusts unreleased metadata. An isolated clean-clone check of
  `ec9094d` still rejects the line-wrapped README parity warning; the
  prospective worktree fixes that checker pattern. The final prospective tree
  still needs its own exact archive check.
- One quality run initially encountered sandbox/cache-write restrictions.
  The permitted rerun passed; this was an execution-environment failure, not a
  package assertion failure.
- At the initial exact commit, remote Core and Documentation passed. Remote
  Quality reflected the editor-workspace release-integrity failure. Coverage
  executed the tests successfully but Codecov rejected its upload with
  `Repository not found`.
- During this audit, the local GitHub CLI credential became invalid and the
  connected GitHub app could not read the private repository. Exact remote
  workflow evidence for `ec9094d` is therefore unverified until authentication
  is restored.
- Required status checks and pull-request reviews are not enforced by the
  observed branch rules. Only deletion and non-fast-forward protection were
  observed.

## Current completion gate

Run:

```sh
julia --startup-file=no --project=. scripts/check_qetlab_completion.jl --check
julia --startup-file=no --project=. scripts/check_qetlab_completion.jl --strict
```

The first command must validate the generated 127-row plan and print the honest
current counts. The second command must remain nonzero while any public row is
partial, deferred, blocked, pending, stubbed, undocumented, untested, or
missing required provenance, or while a required private helper lacks a tested
implementation or no-loss supersession.

Strict completion is necessary but not sufficient for release: supported
platform CI, accepted coverage ingestion, applicable authoritative oracle
evidence, legal review, performance review, and the maintainer's
non-delegable review remain separate gates.
