# Session handoff

Snapshot date: 2026-08-09. This file describes the local convergence worktree,
not a committed or published release.

## Repository state

<!-- qetlab-current-claims: begin -->

- Branch: `main`.
- `HEAD` and `origin/main`:
  `8b2fcbafc646c528aade8bbf695b1112e45511a9`.
- Remote: private `origin` at
  `https://github.com/aenictusGitHub/QuantumEntanglementTools.git`.
- Worktree: intentionally dirty with release-hardening changes layered on the
  previously committed tutorials, user-experience work, and schema-2 browser
  generator. The current pass adds stricter release/oracle gates, provenance
  enforcement, bounded combinatorial resources, non-projecting Choi-to-Kraus
  conversion, bounded optional-backend output, and GitHub maintenance/release
  scaffolding, plus ten project-native symmetric-state APIs with exact
  occupation indexing, multiqudit Dicke states, compressed collective
  operators, split isometries, and direct reductions. Use `git status --short`
  for the exact live list; do not reset or clean it.
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
- The runtime API contains 477 exports (347 native/module and 130
  `MATLABCompat`) matched by 477 provenance records.
- The exact current-tree package corpus passes 9,484/9,484 assertions: 9,400
  core assertions plus 84 executable-tutorial assertions, on
  Julia 1.12.6 and the installed Julia 1.10.11. Seven tutorial scripts are
  covered; 36 assertions exercise the new seeded Schmidt-decomposition and
  exact Tiles-UPB bound-entanglement workflows.
- The complete JuMP/Hypatia/SCS extension passes 847/847 assertions on both
  installed Julia lines. The exact EntanglementDetection.jl 0.2.2 extension
  passes 141/141 assertions on Julia 1.12.6; its dependency resolver floor is
  Julia 1.11.
- All 27 source-free QETLAB comparator groups pass 1,347/1,347 assertions on
  both installed Julia lines, with fixture SHA-256 files verified. MATLAB was
  not run.
- Strict documentation builds pass on both installed Julia lines. The math
  compatibility scan covers 55 Markdown files, 302 inline spans, and 128
  display blocks. The schema-2 browser generator passes 178/178 JavaScriptCore
  checks; its 13-module Julia bundle contains 130 structural and runtime
  assertions and completes on both lines.
- Aqua passes 11/11, the representative JET set passes 25/25, independent
  randomized matrix-predicate validation passes 130/130 on both lines, and the
  formatter gate passes.
- All 114 declared quick benchmark cases completed without failure on the
  current uncommitted tree with one Julia and one BLAS thread. This is local
  smoke evidence only; targeted minima remain diagnostics rather than a stable
  comparative-performance result or regression baseline.

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

The performance pass preserves every public name and signature. Homogeneous
dense BLAS-float Kraus and paired operator-sum factors now form Choi matrices
from compact factor-column products; mixed, exact, arbitrary-precision, and
sparse inputs retain the previous termwise path. Compound matrices use
allocation-light, partially pivoted kernels for standard floating-point
two-by-two and three-by-three minors, exact formulas for widened exact types,
and the standard-library determinant for other numeric types and larger
minors. See `docs/BENCHMARK_REPORT.md` for the bounded measurements and
numerical cross-checks.

The tutorial pass adds two independently written, core-only workflows inspired
by QETLAB's introductory examples. The seeded Schmidt workflow reconstructs a
random two-qutrit pure state manually and through `tensor_sum`. The Tiles
workflow verifies the UPB and complementary projector in exact arithmetic,
proves PPT by an exact partial-transpose identity, derives entanglement from the
range criterion, and separately records the floating PPT boundary as `unknown`
before realignment certifies entanglement. Both are live Documenter examples
and retain the pinned QETLAB source and license relationship in
`PROVENANCE.toml`.

The user-experience pass keeps all existing entry points and detailed status
fields while adding a shared interpretation layer:
`conclusion`, `is_conclusive`, `is_certified`, and `explain`.
`validate_density_matrix` reports shape, dimension, finiteness, trace,
Hermiticity, positivity, sparse-densification, and numerical-boundary details
without repairing input. `available_separability_strategies`,
`describe_strategy`, and `backend_status` make method selection and optional
integration readiness discoverable without loading a package or solving a
problem. Rich terminal displays expose evidence and warnings, while compact
single-line displays remain unchanged.

The documentation now starts from user tasks, orders tutorials progressively,
and checks that every copyable `JuMPBackend` example imports JuMP explicitly.
The `main` documentation workflow packages and deploys the rendered
Documenter site, including equations and the browser-local code generator.
Deployment still requires the repository owner to select **GitHub Actions** as
the Pages source after this work is committed and pushed. On a private
repository, Pages availability and site visibility depend on the GitHub plan
and repository settings.

The generator now has 11 bounded state families and nine curated workflows.
The new families are an exact-arithmetic Tiles-UPB complement and seeded
Hilbert--Schmidt/Bures density matrices; the random diagnostics preset uses a
full-rank state so its purity and entropy routes remain demonstrative without
clipping a rank-boundary eigenvalue. New analysis choices expose density
validation, marginals, pure-state Schmidt data, bounded core-only separability
profiles, and side-effect-free backend readiness. Result-helper output keeps
the raw status visible, while optional generated assertions are structural and
do not hard-code parameter-dependent entanglement conclusions.

Portable schema-2 JSON is canonical, capped at 16 KiB as UTF-8, and rejects
unknown keys, malformed or future schemas, free-form identifiers, and
non-finite or out-of-range values. Schema-1 input migrates with every new
analysis and output option disabled. The UI adds accessible errors, resource
estimates, stale-output protection, and local copy/download/load controls; it
does not execute Julia, upload data, use browser storage, or regenerate code
until the user requests it.

In particular, the UPB separable-discrimination path does not claim a
reachable exact positive certificate for floating solver input. Residual-checked
floating feasible points and computed-cone Farkas evidence remain structured
inconclusive results. The Bell-qubit numerical relaxation likewise cannot
enter a legacy certified tuple path.

## Legal and provenance hardening (2026-07-31)

The uncommitted release-hardening tree now gives each previously identified
QETLAB-informed nonlocal-game and optimization source file, including the JuMP
PSD materialization, an exact upstream-filename/revision BSD-2-Clause header.
`PROVENANCE.toml` maps the PSD extension as an additional implementation file,
and `scripts/check_public_api.jl` now rejects every QETLAB-attributed Julia file
whose preamble omits QETLAB, the pinned revision, BSD-2-Clause, or the retained
license path. The checker passes for all 477 public bindings on Julia 1.12.6 and
1.10.11.

`UpstreamManifest.toml` now pins the inspected
PermutationalInvariantDynamics.jl generator reference at commit
`49c64b1c0fc5b301531582d470144c5b6b3d4030`, including tree, license, REUSE,
and four inspected-path hashes. No reference source or template was copied.
The manifest also records the DOI, publisher-PDF hash, copyright notice, and
reference-only relationship for the Louvet--Serrano-Ensástiga--Bastin--Martin
SAPPT paper. The executable tutorial carries the matching citation and
non-redistribution notice. Generated Julia examples now carry BSD-3-Clause
SPDX headers, complete-license location, user-material warning, and the full
paper DOI where the symmetric witness coefficients are emitted.

Focused checks passed: the 178-check JavaScriptCore generator harness, its
current-Julia generated bundle, the symmetric-SAPPT tutorial, the strict
Documenter build, and both Julia-line public-API/provenance checks. The full
repository formatter was then applied and its repository-wide gate passed.

## Reachable-history visibility blocker (2026-07-31)

The candidate tree and its source archive contain no `QUBIT4MATLAB` path or
text reference. The name is nevertheless reachable in the existing Git
history through commits `9b0d0d3` and `6bf8d61`, including former documentation
and license paths. This does not contaminate the candidate archive, but it means
that changing the current private repository to public visibility would expose
the removed history.

Do not rewrite or replace history automatically. Before public visibility, the
maintainer must explicitly choose and record either a coordinated history
rewrite/clean repository migration or a reviewed legal decision accepting
tree-only removal. A private GitHub prerelease and a public repository remain
separate decisions.

## Final local validation commands

The following commands passed on the current convergence tree unless a
qualification is stated:

```sh
julia --compiled-modules=no --startup-file=no --project=. test/runtests.jl
julia +1.10 --compiled-modules=no --startup-file=no --project=. test/runtests.jl

julia --compiled-modules=no --startup-file=no --project=. \
  -e 'using Test, QuantumEntanglementTools; include("test/symmetric_states.jl")'
julia +1.10 --compiled-modules=no --startup-file=no --project=. \
  -e 'using Test, QuantumEntanglementTools; include("test/symmetric_states.jl")'

julia --compiled-modules=no --startup-file=no --project=. tutorials/runtests.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. tutorials/runtests.jl

julia --compiled-modules=no --startup-file=no \
  --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no \
  --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl

julia --compiled-modules=no --startup-file=no \
  --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl

julia --compiled-modules=no --startup-file=no --project=docs docs/make.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=docs docs/make.jl

julia --compiled-modules=no --startup-file=no --project=. \
  scripts/check_optional_backend_docs.jl

osascript -l JavaScript docs/test/code_generator_jxa.js <temporary-directory>
julia --compiled-modules=no --startup-file=no --project=docs \
  <temporary-directory>/generated_smoke.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=docs \
  <temporary-directory>/generated_smoke.jl

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
- Installed Julia lines used for final local evidence: 1.12.6 and 1.10.11.
  The installed minimum-line binary is 1.10.11, despite older documentation
  referring to 1.10.0.
- The EntanglementDetection.jl 0.2.2 isolated environment resolves only on
  Julia 1.11 or newer because of Ket 0.9; this does not raise the core package
  minimum.
- Sandboxed Julia initially failed while attempting to write
  `~/.julia/logs/manifest_usage.toml.pid`; permitted reruns passed.
- The sandboxed Aqua persistent-task probe also could not create Julia registry
  and manifest-usage pidfiles. Its permitted rerun passed Aqua 11/11 and the
  representative JET checks 25/25.
- The first fresh-depot archive smoke was blocked only by sandbox DNS. The
  permitted reruns downloaded registry metadata into disposable depots and
  passed on both Julia 1.12.6 and 1.10.11.
- The ignored docs manifest required a version-specific refresh for each Julia
  line. It currently reflects Julia 1.12.6 and is not a tracked source change.
- The first Julia 1.10 documentation attempt found the Julia 1.12-flavored
  ignored manifest without its version-specific OpenSSL artifact. Refreshing
  that environment with Julia 1.10 resolved the artifact; the strict build
  passed, and the ignored manifest was then restored with Julia 1.12.6.
- Node.js was unavailable. The documented JavaScriptCore/JXA fallback passed
  all 178 generator checks and emitted the two-version smoke bundle.
- `cffconvert` was not installed, so schema-tool validation of `CITATION.cff`
  remains a release-checklist item. The package release preflight's metadata
  consistency checks passed.
- No Pages deployment was attempted because this worktree was not committed or
  pushed. After pushing, select **Settings → Pages → GitHub Actions**. Treat
  site visibility as a separate decision for this private repository.
- Documentation emitted only size warnings: five API/migration pages exceeded
  100 KiB but remained below the 200 KiB hard limit; `search_index.js`
  exceeded the 500 KiB warning threshold.
- The exact current dirty tree completed all 114 declared quick benchmark
  cases with `--quick --no-save`; those timings remain local smoke evidence,
  not a comparative-performance claim.
- Claim reconciliation emitted Julia 1.12 world-age warnings for anonymous
  test bindings created by repeated `include` calls. The measured core and
  extension testsets still passed; no package runtime failure was observed.
- MATLAB and proprietary/commercial optimization backends were not used.

## Symmetric-state toolkit addition (2026-08-09)

Ten project-native public functions now cover exact symmetric-subspace sizes,
occupation enumeration and rank/unrank, multiqudit Dicke states, product-state
coordinates, collective one-body operators, bipartition isometries, direct
compressed reductions, and maximally mixed symmetric states. Their 1,014
focused assertions pass on Julia 1.12.6 and 1.10.11; both complete package runs
pass 9,484/9,484. The public-API/provenance gate, generated completion-artifact
freshness check, strict completion checker, and strict Documenter build pass.

Occupation-tuple materialization has a hard 1,000-level safety cap for Julia
1.10 compiler stability. Large sparse wrappers may still require a temporary
canonical sparse copy before resource rejection, and nontrivial split/reduction
construction materializes guarded occupation workspaces. These are documented
resource limitations, not silent densification or correctness fallbacks.

## Remaining external gates

There is no remaining public QETLAB row or required private helper in the local
completion queue. Remaining work requires authority or infrastructure outside
this implementation task:

1. maintainer mathematical, API, provenance, licensing, generated-ledger, and
   release review;
2. an explicitly authorized commit and push;
3. selection of **GitHub Actions** as the repository's Pages source, followed
   by verification of the deployed site's intended access level;
4. exact-commit supported-platform CI, Codecov ingestion, and archive
   validation;
5. separate explicit decisions for tagging, a private GitHub release, public
   visibility, and Julia General registration;
6. authoritative MATLAB comparison if broader MATLAB/QETLAB parity is ever to
   be claimed.

Do not describe this package as an official QETLAB project, and do not turn
the local terminal ledger into a universal mathematical, parity, performance,
or production-readiness claim.
