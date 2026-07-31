# AGENTS.md

This file is the operational entry point for maintainers and coding agents working
in this repository. Read it together with
[`docs/SESSION_HANDOFF.md`](docs/SESSION_HANDOFF.md) and
[`docs/PORTING_STATUS.md`](docs/PORTING_STATUS.md) before changing code.

## Project state

`QuantumEntanglementTools` is unreleased pre-1.0 software whose development goal
is complete behavioral coverage of the documented public API at the pinned
QETLAB revision. Classification, source disposition, and the strict static
implementation ledger are complete. Semantic validation, oracle evidence,
supported-platform remote CI, and QETLAB/MATLAB behavioral parity remain
separate milestones. Never infer completeness from a file existing or a symbol
being exported, and do not widen a milestone's verified scope beyond functions
backed by the provenance, test, and documentation ledgers.

The package name is provisional. It is intentionally neutral and must not be
described as an official QETLAB project.

<!-- qetlab-current-claims: begin -->

- The inventory is pinned to QETLAB revision
  `d8589610f00cff106537268dee2e2a1153f3a601`. Its strict static ledger reports
  127/127 public rows are verified with the required final status, 36/36 internal
  helpers have terminal dispositions, the completion queue contains 0 public
  rows, 0 required internal helpers remain, and there are 0 static completion
  failures.
- The public API contains 467 public bindings (337 native/module and 130
  `MATLABCompat`) with matching provenance entries.
- The 8,360-assertion full package suite passed 8,360/8,360: 8,276 core
  assertions plus 84 executable-tutorial assertions, including 36
  for the two QETLAB-introduction workflows. It passes on Julia 1.12.6 and the
  installed Julia 1.10.0. The full optional JuMP suite passed 846/846 on both
  Julia lines. The exact
  EntanglementDetection.jl 0.2.2 extension passed 130/130 focused assertions on
  the current compatible Julia.
- All 114 declared quick benchmark cases completed for the clean `f32dd233`
  baseline and the candidate now committed as `a50f516`; targeted paired
  observations remain local
  diagnostics, not a stable performance baseline.
- These local and static results do not establish QETLAB/MATLAB parity,
  supported-platform remote CI, comparative performance, API stability,
  release approval, or non-delegable human review. No version has been tagged
  or published.

<!-- qetlab-current-claims: end -->

## Sources of truth

- `porting/qetlab_inventory.toml` is the completeness ledger once generated.
- `UpstreamManifest.toml` records upstream revisions and license evidence.
- `PROVENANCE.toml` records the origin, implementation kind, tests, and docs for
  each public function.
- `docs/src/conventions.md` records user-visible mathematical conventions.
- `docs/PORTING_STATUS.md` summarizes verified milestone status.
- `docs/SESSION_HANDOFF.md` records the current working state and blockers.
- `porting/qetlab_completion_plan.toml` records the generated all-public-row
  completion program.
- `artifacts/convergence/current_snapshot.toml` records reconciled current
  counts and static-audit candidates.

If these disagree, stop making broad claims, resolve the discrepancy from primary
evidence, and update all affected records.

## Non-negotiable engineering rules

1. Inspect `git status` before editing and preserve unrelated or pre-existing
   work. Do not reset or clean the worktree.
2. Do not export placeholders, fabricated results, or unchecked numerical
   conclusions.
3. A public function is not complete until its provenance, specification,
   implementation, tests, and documentation are all present.
4. Preserve meaningful numeric element types and sparse structure. Any
   densification must be explicit, documented, guarded, and tested.
5. Randomized APIs accept an explicit `rng::AbstractRNG`; never mutate the
   caller's global random stream.
6. Distinguish certificates from necessary tests and heuristics. An
   inconclusive computation is `unknown`, not a negative mathematical result.
7. Optional solvers and integrations belong in package extensions. Do not use
   private dependency APIs, type piracy, runtime method injection, or
   `Requires.jl`.
8. Do not claim parity or performance without recorded validation or benchmark
   evidence.
9. Keep source-derived licensing and attribution at file/function granularity.

## Julia conventions

- Minimum supported Julia version: 1.10. The optional
  EntanglementDetection.jl environment has an effective Julia 1.11 resolver
  floor because of its Ket 0.9 dependency.
- Prefer one top-level module and standard `AbstractVector`/`AbstractMatrix`
  inputs.
- Use lowercase `snake_case` for the Julia-native API. MATLAB-compatible names
  belong in a separate compatibility namespace.
- Use `adjoint` where conjugation is required and `A \ b` or factorizations
  instead of `inv(A) * b`.
- Validate subsystem dimensions and indices centrally. Never silently normalize,
  symmetrize, clip, or repair user input.
- Document subsystem order, vectorization, normalization, shapes, tolerances,
  failure behavior, and complexity for every important operation.

## Expected checks

Run the smallest relevant test while iterating, then the full applicable checks:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --startup-file=no --project=. tutorials/runtests.jl
julia --startup-file=no --project=. scripts/check_release.jl
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
julia --project=benchmark benchmark/benchmarks.jl --quick --no-save
julia --startup-file=no --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
```

Instantiate dedicated optional-extension environments as documented in
`CONTRIBUTING.md`; do not add an optional backend to the core test target.

Before handing work off, update the status and handoff documents with commands
actually run, failures, uncommitted files, and external blockers. Do not leave
critical information only in chat.
