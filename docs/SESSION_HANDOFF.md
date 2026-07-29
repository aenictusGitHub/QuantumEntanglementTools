# Session handoff

Snapshot date: 2026-07-29. Refresh this file at the end of each substantive
session; do not assume the repository state below remains current.

## Repository state

- Branch: `main`.
- Candidate code commit:
  `0b63359159e1c0c1527c8753f78b61940701eb25`
  (`release: prepare experimental v0.1.0 candidate`). A documentation-only
  evidence commit follows it; use `git rev-parse HEAD` for the final local
  handoff revision.
- Candidate base: `485b6a3` (`docs: add entanglement example code generator`),
  which is also the current `origin/main` revision before release preparation.
- Remote: private `origin` at
  `https://github.com/aenictusGitHub/QuantumEntanglementTools.git`.
- Package: `QuantumEntanglementTools`, UUID
  `45675e5b-5c8b-4983-b92d-4c3725d56c4e`, experimental version `0.1.0`.
- Package author metadata: `John MARTIN <jmartin@uliege.be>`.
- Expected post-commit worktree: the user-owned staged
  `docs/src/QuantumEntanglementTools.code-workspace` remains outside the
  inventory-review commit, together with ignored local manifests, built
  documentation, benchmark output, oracle output, and development checkouts.
  Never include the workspace file in an agent-authored commit.
- No push, tag, GitHub release, visibility change, branch-rule change, or
  General-registry submission was made during release preparation.

## Candidate scope and evidence

This is an experimental `v0.1.0` candidate for the package-owned, documented
API. It is not a full QETLAB port, a parity claim, or a production-stability
claim. The generated inventory contains 163 upstream files and 503 dependency
edges with no automatically detected cycle. All 163 rows now have source-level
dispositions and none remains pending. Among the 127 public rows, 63 are
implemented, 15 partial, 19 deferred, and 30 blocked with explicit reasons.
The 36 private helpers are recorded separately as 13 internally superseded, two
unreferenced and intentionally excluded, four partially covered, and 17
deferred with their parent scope. The public API and provenance gate covers 214
bindings.

The 2026-07-29 inventory sweep inspected all 87 formerly pending MATLAB files
against the pinned source, current API, provenance, tests, and documentation.
It added a per-row reviewed overlay entry and the human-readable
`docs/INVENTORY_REVIEW.md`. The generator now reports
`review_status = "source_review_complete"`, and the public-API checker fails if
any future upstream row lacks an overlay entry or remains automatically
pending. Source review is a classification milestone, not implementation or
human-release-review evidence.

The integrated local package corpus contains 2,417 assertions: 2,369 core
assertions plus 48 executable-tutorial assertions. It passes on Julia 1.12.6
and the minimum supported Julia 1.10.11. Focused counts are:

- Tier A subsystem kernel: 345/345.
- Tier B operators, states, and explicit-RNG random objects: 954/954.
- Tier C channels and maps: 177/177.
- Tier D measures and necessary criteria: 168/168.
- Native certificate-first entanglement pipeline: 68/68.
- Tier E coherence: 54/54.
- Tier E product analysis: 179/179 native and 52/52 compatibility.
- Tier E matrix analysis: 131/131 native and 34/34 compatibility.
- Tier E matrix predicates: 170/170 native and 37/37 compatibility.
- Standalone executable tutorials: 48/48.

The exact EntanglementDetection.jl 0.2.2 integration remains a weak-dependency
extension. Its local Julia 1.12.6 suite passes 125/125 assertions. Every search
runs in a bounded child process, and upstream results remain uncertified
candidate evidence reported as `unknown`. The compatible dependency graph has
an effective Julia 1.11 resolver floor because of Ket 0.9; the core package
continues to support Julia 1.10.

Source-free Octave/QETLAB comparisons pass 338 assertions in total: Tier A 52,
Tier B 72, Tier C 28, Tier D 34, Tier E coherence 25, Tier E product 68, and
Tier E matrix 59. Matrix predicates additionally pass 130/130 independent
randomized eigenspectrum/minor checks. MATLAB is absent, and Octave evidence is
function-specific rather than a general MATLAB-equivalence claim.

All 42 quick benchmark cases completed. Aqua passes 11/11, all 25
representative JET probes pass, the formatter gate passes, the inventory
check remains at 163 files/503 edges/zero detected cycles, and the
public-API/provenance gate remains at 214 bindings. CFFConvert 2.0.0 validates
`CITATION.cff` against schema 1.2.0. The offline upstream audit passes all
10 pin/license checks. A scoped history scan found no common token,
private-key, or suspicious credential-filename signatures; this is a heuristic
check, not proof that history contains no secret.

Strict Documenter builds, doctests, exported-doc checks, and live tutorials
pass on Julia 1.12.6 and Julia 1.10.11. The API reference is split into native
and compatibility pages; the native page is about 157 KiB, below the 200 KiB
hard limit. A Julia 1.12.6 CI-mode pretty-URL build passed. The browser-local
code generator passes 71 deterministic JavaScriptCore assertions, and the
generated smoke program covers all nine state families on Julia 1.12.6 and
Julia 1.10.11.

The exact Git archive of candidate code commit
`0b63359159e1c0c1527c8753f78b61940701eb25` passes the release check and
fresh-depot smoke on Julia 1.12.6 and Julia 1.10.11. Both runs produced tar
SHA-256 `226302eb31191306111fad3aae293548beceec9c5c45640ddda1a110b22c2396`.
After this evidence was recorded in a documentation-only commit, the same
two archive gates were rerun successfully against the final local `HEAD`.

## Release hardening completed

- Public plan objects now own read-only lookup arrays, and Kraus
  representations own read-only operator collections. Validated inner
  constructors prevent bypassing those invariants.
- Public array APIs reject non-one-based axes before using positional indexing;
  copy-only compatibility paths may preserve caller axes without indexing.
  Channel consumers validate cached Kraus and matrix representation axes,
  shapes, and finite entries before use.
- `brauer_states` rejects excessive matching and nonzero counts during capped
  counting, before enumeration or allocation. Its `max_matchings=100_000` and
  `max_nonzeros=1_000_000` guards can be disabled only explicitly with
  `nothing`.
- Standalone docstrings were added for Choi, superoperator, and Kraus
  representation conversion.
- Release metadata, citation files, security policy, changelog, contribution
  guidance, status pages, and the release checklist now describe the scoped
  experimental candidate and avoid unsupported parity claims.
- The release checker validates worktree state, metadata consistency,
  extension entry points, forbidden paths and byte signatures, symlinks,
  submodules, Git LFS pointers, tag shape/target, and the exact Git archive. Its
  optional archive smoke uses a fresh depot.
- GitHub Actions are pinned to full immutable commit SHAs, checkout credentials
  are not persisted, tag triggers are present, core coverage includes Julia
  1.10 macOS and Windows, optional-extension path triggers include all source
  changes, Quality obtains its pinned comparison checkout, and coverage uses
  OIDC with upload failure made fatal.

## Remote evidence and external state

- Core and documentation workflows passed at predecessor commit `485b6a3`.
- The six-job Julia 1.11/1.12 Linux/macOS/Windows optional-extension matrix
  passed at predecessor commit `6bf8d61`.
- The predecessor Quality workflow failed only because its ignored upstream
  comparison checkout was absent. The predecessor Coverage job was green even
  though Codecov rejected the tokenless upload. Both configurations are
  corrected locally, but neither correction has remote evidence on the exact
  candidate commit.
- The GitHub repository currently has deletion and non-fast-forward protection
  on `main`, but no required status checks. Repository Actions policy permits
  unpinned actions even though this candidate pins its own workflows.
  Vulnerability alerts, automated security fixes, secret scanning, and code
  scanning were observed disabled.
- No tag or GitHub release exists.

## Failures encountered and resolved

- The formatter initially reported drift in
  `scripts/build_upstream_inventory.jl` and the new release files; formatting
  was applied and the full-tree gate then passed.
- A documentation environment manifest resolved by one Julia version caused a
  cross-version local docs run to fail. The Julia 1.12 copy was moved to
  `/private/tmp/qet-docs-manifest-julia112-20260729.toml`; developing the local
  package and resolving a fresh ignored Julia 1.10 docs manifest fixed the run.
  No manifest is committed.
- A sandboxed documentation run could not write Julia's package-usage log.
  Rerunning the same build with the required filesystem permission completed
  successfully; this was an environment restriction, not a documentation or
  package failure.
- Initial sandboxed package, focused-test, and quality invocations could not
  write Julia compiled-cache pidfiles. The same commands passed with permission
  to use Julia's local cache; this was an environment restriction.
- Node.js was unavailable for the browser-generator test. The documented macOS
  JavaScriptCore fallback passed 71 assertions and emitted a smoke program that
  passed on Julia 1.12.6 and Julia 1.10.11.
- `cffconvert` was not on the interactive shell's `PATH`; the isolated
  CFFConvert 2.0.0 validation environment at
  `/private/tmp/qet-cff-validator-20260729` validated the citation against
  schema 1.2.0.
- The prior Quality and Codecov remote shortcomings above remain unverified
  until an exact-candidate push and remote run.
- During the completed-inventory sweep, sandboxed Quality and docs invocations
  again failed only because Julia could not create cache or usage-log pidfiles.
  The permitted reruns passed. The first permitted docs build then exposed an
  ignored stale cross-version `docs/Manifest.toml` entry for `MbedTLS_jll`;
  developing the local package, resolving, and instantiating the ignored docs
  environment removed that stale entry, and the exact docs build passed.

## Inventory-sweep validation

The following commands passed after the 87-row source review:

```sh
julia --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
julia --startup-file=no --project=. scripts/check_public_api.jl
julia +1.10 --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
julia +1.10 --startup-file=no --project=. scripts/check_public_api.jl
julia --startup-file=no --project=quality quality/format.jl
julia --startup-file=no --project=quality quality/run_quality.jl
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=docs docs/make.jl
```

Results were 163 inventoried and reviewed rows, 503 dependency edges, zero
detected cycles, 214 API/provenance bindings, Aqua 11/11, JET 25/25, the full
2,417-assertion package corpus, and a strict Documenter build. Documenter
retained the known non-fatal warning that the approximately 158 KiB native API
page is above the 100 KiB warning threshold but below the 200 KiB hard limit.

## Remaining release gates

1. Have a maintainer perform the non-delegable human review required for
   Codex-assisted registry submissions, including mathematical conclusions,
   public API, licenses, and generated changes.
2. With explicit authorization, push the candidate and require successful
   Core, Docs, Quality, Coverage, and optional-extension runs at that exact SHA;
   confirm Codecov accepted the upload rather than relying only on a green job.
3. Decide separately whether this remains a private experimental release or
   becomes a public General-registry candidate. General registration requires a
   public repository and the conventional `.jl` repository URL; the current
   private URL intentionally does not satisfy that preflight.
4. Before any public visibility change, decide whether to preserve or rewrite
   historical commits that contain references to material removed from the
   current tree. History rewriting is destructive and was not authorized.
5. Only after the exact remote gates and human review pass, update the release
   date if necessary, create an annotated `v0.1.0` tag, rerun tagged-archive
   checks, and create the GitHub release. Do not infer authorization to publish
   from this handoff.

## Commands to rerun

```sh
git status --short
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia +1.10 --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=. tutorials/runtests.jl
julia --startup-file=no --project=docs docs/make.jl
julia +1.10 --startup-file=no --project=docs docs/make.jl
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
julia --startup-file=no --project=quality quality/run_quality.jl
julia --startup-file=no --project=quality quality/format.jl
julia --startup-file=no --project=benchmark benchmark/benchmarks.jl --quick --no-save
julia --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
julia --startup-file=no --project=. scripts/check_public_api.jl
julia --startup-file=no --project=. scripts/validate_matrix_predicates.jl
julia --startup-file=no --project=. scripts/check_release.jl --allow-dirty
julia +1.10 --startup-file=no --project=. scripts/check_release.jl --allow-dirty
julia --startup-file=no --project=. scripts/check_release.jl --archive-smoke
julia +1.10 --startup-file=no --project=. scripts/check_release.jl --archive-smoke
```

MATLAB is absent locally. Octave 11.3.0 is available but must not be treated as
an equivalent oracle without function-specific compatibility evidence.

## GitHub documentation-preview repair

The commit containing this section repairs the two repository-preview failures
reported from GitHub:

- 110 Documenter-only inline double-backtick equations across four pages now
  use `$...$`, which is rendered by GitHub and by Documenter. The 49 fenced
  display equations are unchanged, and table-cell bipartitions use `\vert` so
  that GitHub does not parse them as column separators.
- The full code-generator form moved from an `@raw html` Markdown block into
  the existing deferred UI asset. The source page now renders cleanly on
  GitHub, explains the platform limitation, and points private-repository users
  to the seven-day Documentation workflow artifact, local build commands, and
  executable tutorials. The built Documenter page still inserts the complete
  interactive form after the `Interactive generator` heading.
- `scripts/check_docs_math_compat.jl` is invoked by `docs/make.jl`, and changes
  to that script trigger the Documentation workflow. It rejects the incompatible
  inline form, unbalanced dollar delimiters, `\operatorname`, literal math
  pipes in table cells, and unclosed fences.

Focused validation completed:

```sh
julia --startup-file=no --project=. scripts/check_docs_math_compat.jl
julia --startup-file=no --project=. -e 'include("scripts/check_docs_math_compat.jl"); mktemp() do path, io; write(io, "bad ``x``\n"); close(io); result = DocsMathCompatibility.check_file(path); @assert !isempty(result.errors); end; println("Checker rejection self-test passed.")'
/usr/bin/osascript -l JavaScript docs/test/code_generator_jxa.js \
  /private/tmp/qet-generator-smoke-20260729-final
julia --startup-file=no --project=docs \
  /private/tmp/qet-generator-smoke-20260729-final/generated_smoke.jl
julia --startup-file=no --project=docs docs/make.jl
git diff --check
```

The compatibility scan passed 25 pages, 110 inline spans, and 49 display
blocks. Its synthetic double-backtick rejection test passed. JavaScriptCore
passed 71 generator assertions, all nine emitted Julia branches passed, and the
strict Documenter build/doctests passed. A live GitHub Markdown API render also
returned math-renderer elements for representative bra-ket, floor, `\mathrm`,
and table expressions. Node.js remains unavailable locally, so the documented
macOS JavaScriptCore fallback was used.

The first sandboxed Documenter rerun failed only when Julia attempted to write
`~/.julia/logs/manifest_usage.toml.pid`; the permitted rerun passed. The sole
build warning remains the known approximately 158 KiB native API page.

No GitHub setting, Pages site, branch, or deployment was changed. The existing
workflow produces an authenticated artifact, not a hosted site. GitHub Pages is
not configured; publishing Pages from this private personal repository would
make the documentation site public and therefore still requires an explicit
maintainer visibility decision. No push was performed as part of this repair.
