# Session handoff

Snapshot date: 2026-07-29. Refresh this file at the end of each substantive
session; do not assume the worktree details below remain current after parallel
work continues.

## Repository state

- Branch: `main`.
- Milestone: `docs: add entanglement example code generator`. This handoff is
  part of that milestone commit; use `git rev-parse HEAD` for its exact
  revision.
- Remote: private `origin` at
  `https://github.com/aenictusGitHub/QuantumEntanglementTools.git`; `main`
  tracks `origin/main`. The separability documentation commits are present on
  the remote. After this milestone, local `main` is two commits ahead of
  `origin/main`; neither the GitHub-math compatibility commit nor the browser
  code-generator commit has been pushed.
- Package: `QuantumEntanglementTools`, UUID
  `45675e5b-5c8b-4983-b92d-4c3725d56c4e`, version `0.1.0`.
- Package author metadata: `John MARTIN <jmartin@uliege.be>`.
- Worktree: clean after the milestone commit, apart from ignored local
  manifests, generated documentation, benchmark output, oracle output, and
  upstream development checkouts.

## Evidence secured

- QETLAB inspected at
  `d8589610f00cff106537268dee2e2a1153f3a601` on 2026-07-28.
- QETLAB root BSD-2 license and Bruno Luong helper BSD-2 notices have recorded
  hashes in `docs/LEGAL.md`.
- EntanglementDetection.jl 0.2.2/tag SHA and MIT license hash are recorded.
  Its observed global RNG, stdout/logging, and BLAS-thread side effects are
  contained by the completed fresh-child-process adapter and covered by the
  dedicated extension suite.
- The state family and rounded witness matrices in Phys. Rev. A 111, 042418
  (2025) were independently audited against the supplied journal PDF. The
  executable tutorial records the GHZ phase mismatch explicitly, gives an
  exact Bernstein-basis block-positivity proof for the printed five-qubit
  witness, and does not claim to rerun the source paper's SDP.
- The documentation includes an original browser-local code generator for nine
  curated state families, six analysis/measure routes, decomposable PPT
  witnesses, and the five-qubit published symmetric witness. Its DOM-free core
  passes 71 deterministic validation assertions, and all nine generated Julia
  branches pass on Julia 1.12.6 and 1.10.11. The strict documentation build
  copies all three local assets. The GPL-3.0-only reference generator was used
  only for interaction-pattern inspiration; no source, styling, prose, or
  templates were reused, as recorded in `docs/LEGAL.md`.
- Markdown and matching API-docstring equations avoid GitHub's unsupported
  operator-name macro; the replacement roman-text forms pass the strict local
  Documenter build.
- Local Julia/OS/BLAS/tool evidence is in `docs/BUILD_ENVIRONMENT.md`.

## Current status

Documentation, governance, legal, citation, ADR, and CI scaffolding are present.
The generated QETLAB inventory contains 127 root-level and 36 helper MATLAB
files, 503 dependency edges, and no automatically detected cycle. The reviewed
overlay covers 76 rows: 63 implemented mappings, 11 partial mappings, and two
explicit deferrals; 87 rows remain pending.

Tier A is implemented with native and `MATLABCompat` APIs for layouts/basis
indices, tensor products/sums, permutations/swaps, partial trace/transpose,
realignment, and symmetric/antisymmetric projectors/bases. Its focused local
Julia 1.12.6 suite passes 304/304.

Tier B implements 22 complete QETLAB mappings for operator bases, named states,
and explicit-RNG random objects, plus the verified scalar bipartite Werner
path. Its focused suite passes 941/941. All six native and six compatibility
random constructors require a leading explicit RNG and pass a global-stream
isolation regression. Multipartite Werner is partial; `RandomSuperoperator` and
`RandomPPTState` are deferred and not exported.

Tier C passes 154/154 focused assertions. Tier D passes 162/162
measures/criteria and 68/68 native-pipeline assertions. Tier E coherence passes
52/52. Tier E product analysis passes 175/175 native plus 52/52 compatibility
assertions; four of its six reviewed mappings are implemented and two remain
partial. Tier E matrix analysis passes 127/127 native plus 32/32 compatibility
assertions across four implemented mappings. Tier E matrix predicates pass
166/166 native plus 37/37 compatibility assertions. Their APIs retain
three-valued boundary outcomes, witnesses, explicit sparse-densification gates,
and combinatorial limits. `IsPSD` remains partial because the pinned CVX
symbolic branch is omitted.

Five deterministic executable tutorials cover subsystem reductions and
ordering, local channel action and representation conversion, separability
certificates, symmetric SAPPT states and witnesses, and certificate-aware
entanglement analysis. Each script runs independently from `tutorials/`;
`tutorials/runtests.jl` supplies a 48/48
automated gate included by `Pkg.test()`. `docs/src/tutorials.md` executes those
same scripts as live Documenter examples instead of publishing copied output.
`docs/src/separability_examples.md` also provides direct, copyable examples for
pure products, mixed product-state decompositions, low-dimensional PPT
certification, higher-dimensional inconclusive pipelines, and the separate
separable-ball result vocabulary. The paper-specific page and tutorial add a
finite five-qubit separable decomposition, same-spectrum SAPPT representatives,
and published plus decomposable witness constructions while keeping the
paper's numerical boundary distinct from package-owned certificates.

The browser-local entanglement example generator complements those fixed
tutorials. Its bounded form emits dense, executable Julia scripts for product,
Bell, diagonal-mixture, GHZ, Dicke, isotropic, Werner, Horodecki, and symmetric
SAPPT states. It rejects arbitrary code and non-finite inputs, preserves
certificate versus necessary-test semantics, never executes or persists user
parameters, and offers copy/download controls. Documentation CI syntax-checks
the JavaScript, runs the deterministic core cases, executes the generated smoke
bundle on both documentation Julia versions, and checks the rendered assets.

The optional EntanglementDetection.jl 0.2.2 integration is implemented as a
weak-dependency Julia extension. Every heuristic search runs in a fresh child
process so upstream RNG, stdout/logging, and BLAS-thread mutations do not alter
the caller. Backend conclusions are retained only as candidate evidence:
adapter reports remain `:unknown` and uncertified, including when the backend
suggests `:entangled` or `:separable`. The dedicated Julia 1.12.6 suite passes
125/125 focused assertions. The adapter checks the loaded package version, not
its source-tree hash; the recorded local checkout and controlled CI resolution
provide source provenance. Its isolated environment has a Julia 1.11 resolver
floor because the compatible Ket 0.9 releases require Julia 1.11; the core
package continues to support Julia 1.10. The configured Julia 1.11/1.12
Linux/macOS/Windows extension workflow has not yet run remotely.

Source-free Octave/QETLAB fixtures pass Tier A 52/52, Tier B 72/72, Tier C
28/28, Tier D 34/34, Tier E coherence 25/25, Tier E product 68/68, and Tier E
matrix 59/59. The matrix artifact has 22 fixtures, including five reviewed
discrepancies for weak majorization and compound boundary behavior; its SHA-256
is `e37685c262ce5982d10dd705cef8c172d49d9c55c89a0a67d4de729a5068f540`.
Octave is not treated as generally MATLAB-equivalent. Matrix predicates have
analytic/property evidence, including 130/130 independent randomized
eigenspectrum/minor cross-checks, but no MATLAB-family oracle.

All 42 quick benchmark cases completed, including exactly three product and
three matrix cases. Aqua passes 11/11 and all 25 representative JET probes
pass; those probes do not cover matrix predicates. The integrated package
corpus contains 2,318 assertions—2,270 core plus 48 executable-tutorial
assertions—and passes on Julia 1.12.6 and Julia 1.10.11. Strict Documenter,
formatter, inventory (163 files/503 edges/zero cycles), and
public-API/provenance (214 bindings) gates pass. Every committed oracle
comparator passes, including the 59-assertion matrix comparator. Remote core,
documentation, quality, and optional-extension CI evidence remains pending.

## Next safe work

1. Review the remaining 87 inventory rows and manually correct classifications.
2. Keep `UpstreamManifest.toml`, `PROVENANCE.toml`, status overlays, and
   generated inventory synchronized as later slices land.
3. Run the remote core/docs/quality matrices and the optional-extension
   Julia 1.11/1.12 platform matrix.
4. Add authoritative MATLAB differential checks when MATLAB is available.
5. Preserve the EntanglementDetection adapter's child-process boundary and
   uncertified-candidate semantics while reviewing remote platform results and
   future pinned-upstream changes.
6. Keep this handoff current with subsequent commits, failures, and
   uncommitted files.
7. Push the two local documentation commits only when explicitly authorized,
   then inspect the remote documentation matrix and generated-site artifact.

## Commands to rerun

```sh
git status --short
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia +1.10 --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --startup-file=no --project=. tutorials/subsystem_reductions.jl
julia --startup-file=no --project=. tutorials/local_channel_noise.jl
julia --startup-file=no --project=. tutorials/entanglement_certificates.jl
julia --startup-file=no --project=. tutorials/separability_examples.jl
julia --startup-file=no --project=. tutorials/symmetric_sappt_witnesses.jl
julia --startup-file=no --project=. tutorials/runtests.jl
/usr/bin/osascript -l JavaScript docs/test/code_generator_jxa.js \
  /tmp/qet-generator-smoke
julia --startup-file=no --project=docs /tmp/qet-generator-smoke/generated_smoke.jl
julia +1.10 --startup-file=no --project=docs \
  /tmp/qet-generator-smoke/generated_smoke.jl
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
julia --project=test/extensions/entanglement_detection -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
julia --startup-file=no --project=benchmark benchmark/benchmarks.jl --quick --no-save
julia --startup-file=no --project=quality quality/run_quality.jl
julia --startup-file=no --project=quality quality/format.jl
julia --startup-file=no --project=test/oracle test/oracle/compare_tier_e_product_oracle.jl \
  test/oracle/fixtures/tier_e_product_octave_11_3_qetlab_d858961.json
julia --startup-file=no --project=test/oracle test/oracle/compare_tier_e_matrix_analysis_oracle.jl \
  test/oracle/fixtures/tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json
julia --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
julia --startup-file=no --project=. scripts/check_public_api.jl
julia --startup-file=no --project=. scripts/validate_matrix_predicates.jl
```

MATLAB is absent locally. Octave 11.3.0 is available but must not be treated as
an equivalent oracle without function-specific compatibility evidence.
