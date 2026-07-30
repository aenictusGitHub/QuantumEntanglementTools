# Session handoff

Snapshot date: 2026-07-30. This file describes the local convergence worktree,
not a committed or published release.

## Repository state

<!-- qetlab-current-claims: begin -->

- Branch: `main`.
- `HEAD` and `origin/main`:
  `ec9094dad43a7531b16b1f1d282a490ccba0c543`.
- Remote: private `origin` at
  `https://github.com/aenictusGitHub/QuantumEntanglementTools.git`.
- Worktree: intentionally dirty with 66 tracked paths modified and 224
  untracked paths at this snapshot. The changes span the QETLAB completion
  implementation, optional extensions, tests, fixtures, generated ledgers,
  documentation, CI definitions, and convergence evidence. Use
  `git status --short` for the exact live list; do not reset or clean it.
- No commit, push, tag, release, repository-visibility change, branch-setting
  change, history rewrite, or registry submission was made.
- Pinned QETLAB source:
  `d8589610f00cff106537268dee2e2a1153f3a601`, with 163 MATLAB files,
  127 public functions, 36 internal helpers, 503 dependency edges, and no
  detected cycle.
- The strict terminal-coverage gate passes: 127/127 public rows are verified,
  the completion queue contains 0 public rows, 36/36 internal helpers have
  terminal dispositions, 0 required helpers remain, and there are 0 static
  completion failures.
- The runtime API contains 458 exports (328 native/module and 130
  `MATLABCompat`) matched by 458 provenance records.
- The exact current-tree package corpus passes 8,117/8,117 assertions,
  including 48 executable-tutorial assertions, on Julia 1.12.6 and the
  installed Julia 1.10.0.
- The complete JuMP/Hypatia/SCS extension passes 836/836 assertions on both
  installed Julia lines. The exact EntanglementDetection.jl 0.2.2 extension
  passes 125/125 assertions on Julia 1.12.6; its dependency resolver floor is
  Julia 1.11.
- All 27 source-free QETLAB comparator groups pass 1,347/1,347 assertions on
  both installed Julia lines, with fixture SHA-256 files verified. MATLAB was
  not run.
- Strict documentation builds pass on both installed Julia lines. The math
  compatibility scan covers 55 Markdown files, 301 inline spans, and 128
  display blocks. The browser generator passes 71/71 JavaScriptCore checks,
  and its nine generated examples pass 35/35 Julia assertions on both lines.
- Aqua passes 11/11, the representative JET set passes 25/25, independent
  randomized matrix-predicate validation passes 130/130 on both lines, and the
  formatter gate passes.
- All 114 declared quick benchmark cases completed with `--no-save`. This is
  execution smoke evidence, not a comparative-performance result or regression
  baseline.

These are bounded local implementation and validation claims. They do not
establish general MATLAB/QETLAB parity, supported-platform remote CI,
comparative performance, API stability, release approval, or the maintainer's
non-delegable review.

<!-- qetlab-current-claims: end -->

## Completion result

All 127 public rows in `porting/qetlab_status.toml` now carry the final
`verified` coverage status required by the execution prompt. Strict mode
deliberately rejects the former nonterminal `implemented` label. All 36 private
helper rows are either internally superseded or intentionally excluded with a
terminal disposition; none required by a public row remains deferred.

The generated sources of truth are:

- `porting/qetlab_inventory.toml` and `.csv`;
- `porting/qetlab_completion_plan.toml`;
- `porting/qetlab_completion_queue.toml`;
- `docs/QETLAB_COMPLETION_PLAN.md`;
- `PROVENANCE.toml`;
- `artifacts/convergence/current_snapshot.toml`.

The public surface uses Julia-native structured results, explicit RNGs and
resource limits, guarded densification, optional solver extensions, and
status-rich numerical outcomes. Necessary criteria, relaxations, heuristics,
and numerical feasibility evidence are not promoted to unchecked Boolean
certificates.

In particular, the UPB separable-discrimination path does not claim a
reachable exact positive certificate for floating solver input. Residual-checked
floating feasible points and computed-cone Farkas evidence remain structured
inconclusive results. The Bell-qubit numerical relaxation likewise cannot
enter a legacy certified tuple path.

## Final local validation commands

The following commands passed on the current convergence tree unless a
qualification is stated:

```sh
julia --compiled-modules=no --startup-file=no --project=. test/runtests.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. test/runtests.jl

julia --compiled-modules=no --startup-file=no \
  --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no \
  --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl

julia --compiled-modules=no --startup-file=no --project=docs docs/make.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=docs docs/make.jl

julia --compiled-modules=no --startup-file=no --project=quality \
  quality/format.jl
julia --compiled-modules=no --startup-file=no --project=quality \
  quality/run_quality.jl

julia --compiled-modules=no --startup-file=no --project=benchmark \
  benchmark/benchmarks.jl --quick --no-save

julia --compiled-modules=no --startup-file=no --project=. \
  scripts/check_qetlab_completion.jl --strict
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/check_public_api.jl
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/build_upstream_inventory.jl --check
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/build_qetlab_completion_plan.jl --check
```

The full oracle invocation sets
`JULIA_LOAD_PATH="$PWD/test/oracle:@:@stdlib"` and executes all 27
`test/oracle/compare_*_oracle.jl` files against their explicitly named frozen
fixtures on both Julia lines. Each fixture digest is checked before and after
the run.

The final claim reconciliation, dirty-worktree release preflight, and isolated
archive smoke passed on the current tree. They must be rerun after any
subsequent edit:

```sh
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/reconcile_project_claims.jl
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/reconcile_project_claims.jl --check --skip-tests
julia --compiled-modules=no --startup-file=no --project=. \
  scripts/check_release.jl --allow-dirty --archive-smoke
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. \
  scripts/check_release.jl --allow-dirty --archive-smoke
```

## Environment notes

- Host: macOS on Apple arm64, Apple M4 Pro.
- Installed Julia lines used for final local evidence: 1.12.6 and 1.10.0.
  The installed minimum-line binary is 1.10.0, despite older documentation
  referring to 1.10.11.
- The EntanglementDetection.jl 0.2.2 isolated environment resolves only on
  Julia 1.11 or newer because of Ket 0.9; this does not raise the core package
  minimum.
- Sandboxed Julia initially failed while attempting to write
  `~/.julia/logs/manifest_usage.toml.pid`; permitted reruns passed.
- The first fresh-depot archive smoke was blocked only by sandbox DNS. The
  permitted reruns downloaded registry metadata into disposable depots and
  passed on both Julia 1.12.6 and 1.10.0.
- The ignored docs manifest required a version-specific refresh for each Julia
  line. It currently reflects Julia 1.12.6 and is not a tracked source change.
- Node.js was unavailable. The documented JavaScriptCore/JXA fallback passed
  the generator checks.
- `cffconvert` was not installed, so schema-tool validation of `CITATION.cff`
  remains a release-checklist item. The package release preflight's metadata
  consistency checks passed.
- Documentation emitted only size warnings: five API/migration pages exceeded
  100 KiB but remained below the 200 KiB hard limit; `search_index.js`
  exceeded the 500 KiB warning threshold.
- Claim reconciliation emitted Julia 1.12 world-age warnings for anonymous
  test bindings created by repeated `include` calls. The measured core and
  extension testsets still passed; no package runtime failure was observed.
- MATLAB and proprietary/commercial optimization backends were not used.

## Remaining external gates

There is no remaining public QETLAB row or required private helper in the local
completion queue. Remaining work requires authority or infrastructure outside
this implementation task:

1. maintainer mathematical, API, provenance, licensing, generated-ledger, and
   release review;
2. an explicitly authorized commit and push;
3. exact-commit supported-platform CI, Codecov ingestion, and archive
   validation;
4. separate explicit decisions for tagging, a private GitHub release, public
   visibility, and Julia General registration;
5. authoritative MATLAB comparison if broader MATLAB/QETLAB parity is ever to
   be claimed.

Do not describe this package as an official QETLAB project, and do not turn
the local terminal ledger into a universal mathematical, parity, performance,
or production-readiness claim.
