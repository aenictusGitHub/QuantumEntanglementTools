# Session handoff

Snapshot date: 2026-09-08. This file describes the local convergence worktree,
not a committed or published release.

## Repository state

<!-- qetlab-current-claims: begin -->

- Branch: `main`.
- `HEAD` and `origin/main`:
  `1d611e4f61f2d740602018dce05afd47f2ecf620`.
- Remote: private `origin` at
  `https://github.com/aenictusGitHub/QuantumEntanglementTools.git`.
- Worktree: intentionally dirty with a bounded performance, numerical-
  stability, regression-test, CI, documentation, and evidence update. Use
  `git status --short` for the exact live list; do not reset or clean it.
- No commit, push, tag, release, repository-visibility change, branch-setting
  change, history rewrite, or registry submission was made in this pass. The
  preceding Pages repair is committed at `1d611e4f61f2d740602018dce05afd47f2ecf620`
  and was successfully deployed.
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
- The exact current-tree package corpus passes 9,759/9,759 assertions: 9,675
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
  compatibility scan covers 56 Markdown files, 324 inline spans, and 132
  display blocks. The schema-2 browser generator passes 178/178 JavaScriptCore
  checks; its 13-module Julia bundle contains 130 structural and runtime
  assertions and completes on both lines.
- Aqua passes 11/11, the representative JET set passes 25/25, independent
  randomized matrix-predicate validation passes 130/130 on both lines, and the
  formatter gate passes.
- All 120 declared quick benchmark cases completed without failure on the
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

The current performance and stability pass preserves every public name and
signature. Diagonal density measures and coherence now avoid dense spectral
factorizations; coefficient-only Schmidt operations use `svdvals`; and dense
order-two additive compounds use direct lexicographic pair indexing. Logical
sparse wrappers remain sparse through subsystem transformations. Fixed-width
integer tensor products, partial traces, and parallel repetitions now reject
overflow instead of wrapping; Boolean partial traces consistently return the
additive `Int` codomain; and large `Float16` equal-superposition constructors
can no longer silently normalize every amplitude to zero. Generic, sparse,
higher-order, exact, and full-vector-producing paths retain their documented
roles. See `docs/BENCHMARK_REPORT.md` for bounded smoke results and explicitly
non-baseline local timing observations.

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
The repository owner selected **GitHub Actions** as the Pages source on
2026-08-09. The Pages repair was committed, its Documentation deployment
succeeded, and the root and `/code_generator/` routes returned HTTP 200 in the
preceding session. On a private repository, Pages availability and site
visibility continue to depend on the GitHub plan and repository settings.

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
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. test/runtests.jl

julia --compiled-modules=no --startup-file=no --project=. \
  -e 'using Test, QuantumEntanglementTools; include("test/symmetric_states.jl")'
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. \
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

julia --compiled-modules=no --startup-file=no scripts/build_docs.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no scripts/build_docs.jl

julia --compiled-modules=no --startup-file=no --project=. \
  scripts/check_optional_backend_docs.jl

osascript -l JavaScript docs/test/code_generator_jxa.js <temporary-directory>
julia --compiled-modules=no --startup-file=no --project=. \
  <temporary-directory>/generated_smoke.jl
/Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia \
  --compiled-modules=no --startup-file=no --project=. \
  <temporary-directory>/generated_smoke.jl

julia --compiled-modules=no --startup-file=no --project=quality \
  quality/format.jl
julia --compiled-modules=no --startup-file=no --project=quality \
  quality/run_quality.jl

julia --compiled-modules=no --startup-file=no --project=benchmark \
  benchmark/benchmarks.jl --quick \
  --output=/private/tmp/qet-performance-final-20260908

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
- A direct Julia 1.10 documentation attempt exposed an ignored Julia
  1.12-generated manifest whose bundled OpenSSL entry had no Julia 1.10 source.
  A fresh Julia 1.10 environment built the strict docs successfully. The new
  `scripts/build_docs.jl` helper always resolves a disposable environment for
  the running Julia version and therefore does not reuse or rewrite the
  ignored `docs/Manifest.toml`.
- Node.js was unavailable. The documented JavaScriptCore/JXA fallback passed
  all 178 generator checks and emitted the two-version smoke bundle.
- `cffconvert` was not installed, so schema-tool validation of `CITATION.cff`
  remains a release-checklist item. The package release preflight's metadata
  consistency checks passed.
- The repository owner selected **Settings → Pages → GitHub Actions**. The
  repaired Documentation workflow subsequently deployed commit `1d611e4`, and
  the root and code-generator routes were verified with HTTP 200. Treat site
  visibility as a separate decision for this private repository.
- Documentation emitted only size warnings: five API/migration pages exceeded
  100 KiB but remained below the 200 KiB hard limit; `search_index.js`
  exceeded the 500 KiB warning threshold.
- The exact current dirty tree completed all 120 declared quick benchmark
  cases with results saved under
  `/private/tmp/qet-performance-final-20260908`; those timings remain local
  smoke evidence, not a comparative-performance claim.
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
pass 9,759/9,759. The public-API/provenance gate, generated completion-artifact
freshness check, strict completion checker, and strict Documenter build pass.

Occupation-tuple materialization has a hard 1,000-level safety cap for Julia
1.10 compiler stability. Large sparse wrappers may still require a temporary
canonical sparse copy before resource rejection, and nontrivial split/reduction
construction materializes guarded occupation workspaces. These are documented
resource limitations, not silent densification or correctness fallbacks.

## GitHub Pages deployment repair (2026-08-09)

The public project root and code-generator route returned HTTP 404 because the
Pages site had not been provisioned and the latest Documentation workflow
failed before artifact upload and deployment. Both supported Julia jobs failed
in the live Tiles tutorial when a platform-dependent boundary classification
left `entanglement_report.evidence` equal to `nothing`.

The local repair keeps the exact rank-four PPT/range-criterion proof unchanged
and keeps its floating PPT result explicitly `unknown`. It runs the independent
realignment certificate on a fixed full-rank state obtained by mixing exactly
`1//1024` of the separable maximally mixed state into the Tiles state. The
neighbor remains PPT by exact partial-transpose equality, has a safely positive
eigenvalue floor, and retains a robust realignment violation. Assertions now
validate report status, method, and evidence type before dereferencing the
evidence. The generator documentation also names the workflow's actual
run-specific `documentation-<run-id>-<run-attempt>` artifact.

Focused Tiles runs and the complete 84-assertion tutorial suite pass on Julia
1.12.6 and 1.10.11. The 178-check JavaScriptCore/JXA generator harness and its
13-module generated Julia smoke pass on both Julia lines. Strict CI-shaped
Documenter builds pass on both lines through disposable per-version docs
environments; they emit only the previously recorded HTML/search-index size
warnings. The formatter gate and `git diff --check` pass. Node.js was
unavailable, so the documented JavaScriptCore/JXA fallback was used.

The repair is committed at
`1d611e4f61f2d740602018dce05afd47f2ecf620`. Its Documentation workflow
deployed successfully, and the live root and `/code_generator/` routes both
returned HTTP 200 in the preceding session.

## Performance and stability pass (2026-09-08)

This uncommitted pass adds structure-aware fast paths for diagonal scalar and
coherence measures, coefficient-only Schmidt analysis, and dense order-two
additive compounds. It also closes silent-result hazards around sparse array
wrappers, fixed-width integer overflow, Boolean partial-trace accumulation,
and large `Float16` equal-superposition normalization. Tests exercise dense,
sparse, wrapped, exact, generic-float, overflow, boundary, nonmutation, and
resource-policy behavior on both supported local Julia lines.

The complete core/tutorial corpus passes 9,759/9,759, both optional JuMP
suites pass 847/847, the compatible EntanglementDetection.jl suite passes
141/141, and all 27 source-free comparator groups pass 1,347/1,347 on both
Julia lines. All 120 quick benchmark cases completed on Julia 1.12.6 with one
Julia and one BLAS thread. Targeted minima remain local diagnostics; they are
not speed guarantees or a comparative regression baseline. Checked arithmetic
currently covers real fixed-width integer element types, not `Complex{Int}`;
ordinary floating-point underflow remains governed by the requested type.

## Remaining external gates

There is no remaining public QETLAB row or required private helper in the local
completion queue. Remaining work requires authority or infrastructure outside
this implementation task:

1. maintainer mathematical, API, provenance, licensing, generated-ledger, and
   release review;
2. an explicitly authorized commit and push of the current performance and
   stability candidate;
3. exact-commit supported-platform CI, Documentation deployment, Codecov
   ingestion, and archive
   validation;
4. verification of the deployed site's intended access level;
5. separate explicit decisions for tagging, a private GitHub release, public
   visibility, and Julia General registration;
6. authoritative MATLAB comparison if broader MATLAB/QETLAB parity is ever to
   be claimed.

Do not describe this package as an official QETLAB project, and do not turn
the local terminal ledger into a universal mathematical, parity, performance,
or production-readiness claim.
