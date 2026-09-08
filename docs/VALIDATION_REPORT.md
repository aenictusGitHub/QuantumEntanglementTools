# Validation report

Evidence date: 2026-09-08.

Status: local completion and release-gate evidence for the uncommitted
convergence worktree; no MATLAB parity, remote supported-platform, production
backend, performance, or release-approval claim.

## Current summary

<!-- qetlab-current-claims: begin -->

- The strict completion checker passes with 127/127 public rows verified with
  final status, completion queue containing 0 public rows, 36/36 internal
  helpers with terminal dispositions, 0 required helpers remaining, and 0
  static completion failures.
- The API checker matches 477 runtime exports (347 native/module and 130
  `MATLABCompat`) to 477 provenance records.
- The full package corpus passes 9,759/9,759 assertions: 9,675
  core assertions plus 84 executable-tutorial assertions, on Julia 1.12.6 and
  the installed Julia 1.10.11.
- Seven standalone tutorials pass 84/84 and independent matrix-predicate
  validation passes 130/130 on both installed Julia lines.
- The complete JuMP/Hypatia/SCS environment passes 847/847 assertions on both
  installed Julia lines. Numerical solver outcomes retain status, residual,
  bound, and certificate metadata and are not automatically mathematical
  certificates.
- The EntanglementDetection.jl 0.2.2 adapter passes 141/141 assertions on Julia
  1.12.6. Its candidate evidence remains conservatively `unknown`.
- All 27 source-free oracle comparators pass 1,347/1,347 assertions on both
  installed Julia lines, and every fixture SHA-256 is verified. These
  function-specific Octave/QETLAB fixtures are supplemental evidence; MATLAB
  was not run.
- Aqua passes 11/11 and the representative JET set passes 25/25. The formatter
  gate passes after applying the repository formatter.
- All 120 declared quick benchmark cases completed without failure on the
  current uncommitted worktree, with one Julia and one BLAS thread. The
  targeted paired observations remain local diagnostics, not a stable
  comparative-performance or regression-baseline claim.
- The dirty-worktree distribution preflight and an isolated fresh-depot archive
  load smoke pass on both installed Julia lines. Dirty mode uses a temporary
  Git index/object store and does not change the repository index; it is not
  exact committed-tree release evidence.

This evidence applies to a dirty local worktree based on
`1d611e4f61f2d740602018dce05afd47f2ecf620`. It is not exact committed-tree
remote CI, API-stability evidence, release approval, or the maintainer's
non-delegable review.

<!-- qetlab-current-claims: end -->

## Historical baseline summary (2026-07-29)

The following section is retained as dated pre-convergence evidence. Its
counts are superseded by the current summary above.

The Tier A subsystem kernel passes 345/345 tests and the Tier B
operators/states/random slice passes 954/954 tests locally on Julia 1.12.6.
The Tier C channels/maps slice passes 177/177 focused package assertions, the
Tier D measures/criteria slice passes 168/168, and the project-native
entanglement pipeline passes 68/68. Tier E coherence passes 54/54, product
analysis passes 179/179 native plus 52/52 compatibility assertions, matrix
analysis passes 131/131 native plus 34/34 compatibility assertions, and matrix
predicates pass 170/170 native plus 37/37 compatibility assertions. The
integrated 2,417-assertion package corpus, comprising 2,369 core assertions
plus 48 executable-tutorial assertions, passes on Julia 1.12.6 and the minimum
supported Julia 1.10.11. The standalone tutorial runner also passes 48/48 on
both Julia versions, and the public-API/provenance gate passes over 214 public
bindings.

The optional EntanglementDetection.jl 0.2.2 focused suite passes 125/125 on
Julia 1.12.6. Its isolated dependency environment has an effective Julia 1.11
resolver floor because compatible Ket 0.9 releases require Julia 1.11; this
does not raise the core package's Julia 1.10 minimum. Every upstream candidate
conclusion remains `unknown` and uncertified at the package boundary. This
establishes an experimental release-candidate baseline, not complete QETLAB
parity or a production-support claim. The predecessor core and documentation
workflows passed at `485b6a3`, and the six-job optional-extension matrix passed
at `6bf8d61`. Those runs do not validate the unpushed release-candidate tree;
exact-candidate remote evidence and accepted Codecov ingestion remain pending.
MATLAB is absent.

| Validation class | Current evidence |
|---|---|
| Package load | Passed as part of `Pkg.test()` |
| Core unit tests | Tier A 345/345, Tier B 954/954, Tier C 177/177, Tier D measures/criteria 168/168, native pipeline 68/68, Tier E coherence 54/54, product analysis 179/179 native plus 52/52 compatibility, matrix analysis 131/131 native plus 34/34 compatibility, and matrix predicates 170/170 native plus 37/37 compatibility pass locally; these 2,369 core assertions plus 48 tutorial assertions give an integrated total of 2,417/2,417 on Julia 1.12.6 and 1.10.11 |
| Executable tutorials | Five repository-native scripts run standalone, through a 48-assertion tutorial gate, from `Pkg.test()`, and as live Documenter examples; the standalone gate passes 48/48 on Julia 1.12.6 and 1.10.11 |
| Tier A analytic tests | Bell reduction/PT spectrum, exact bases/projectors, tensor identities |
| Tier B analytic tests | Operator-basis identities, named-state support/normalization, mixed-state PSD/trace/PPT properties |
| Tier C analytic tests | Kraus/Choi/superoperator round trips, channel application, CP/TP/unital diagnostics, Hilbert--Schmidt duality, complementary/partial maps, and analytic channel/positive-map formulas |
| Tier D analytic tests | Schatten/Ky Fan norms, purity/entropy/fidelity/trace distance, negativity, Schmidt reconstruction/rank, two-qubit concurrence, and PPT/realignment/reduction criterion boundaries and witnesses |
| Tier E coherence analytic tests | Pure and mixed-state coherence formulas, basis transforms, generic precision, strict validation, and the reviewed upstream rank-counting discrepancy |
| Tier E product analytic tests | Operator-Schmidt reconstruction and rank, multipartite product factors and residuals, pure/mixed entanglement of formation, and structured separable-ball certificate boundaries |
| Tier E matrix analytic tests | Strong and weak majorization contracts, exact elementary symmetric polynomials, Cauchy--Binet and compound identities, additive finite differences/eigenvalue sums, sparse paths, and reviewed boundary-shape discrepancies |
| Tier E matrix-predicate analytic tests | Exact and floating PSD, Hermiticity boundaries, principal-submatrix witnesses, strict total positivity, total nonsingularity, structured tri-state results, sparse gates, and combinatorial limits |
| Native entanglement pipeline | PPT certificates, exact `2×2`/`2×3` PPT separability theorem, higher-dimensional `unknown`, pure-state Schmidt certificates, certificate-first attempt order, backend metadata, and invalid/failure semantics |
| Property tests | Index and representation round trips, permutation inverse, trace preservation, PT involution, realignment inverse, projector identities, map duality, physicality properties, fidelity symmetry, Schmidt reconstruction, certificate witness checks, and matrix-predicate invariance/scaling/minor properties |
| Independent formulations | Explicit permutation/projector reconstructions, direct operator/state formulas and perfect matchings, direct Kraus application, Choi block reconstruction, analytic channel/map actions, direct singular-value/eigenvalue identities, phase-independent Schmidt reconstruction, and 130 randomized PSD/minor cross-checks against eigenspectra and compound matrices |
| QETLAB differential tests | Tier A: 13 source-free Octave/QETLAB fixtures pass 52 assertions. Tier B: 18 fixtures pass 72 assertions. Tier C: 7 fixtures pass 28. Tier D: 13 fixtures pass 34. Tier E coherence: 6 fixtures pass 25. Tier E product: 14 fixtures pass 68. Tier E matrix: 22 fixtures pass 59. Matrix predicates have analytic/property evidence but no MATLAB-family oracle. MATLAB unavailable |
| Sparse/generic-number tests | Sparse vectors/matrices and representative `Float32`, rational, `BigFloat`, complex, and abstract `Number` paths pass; later slices add sparse representation reshuffles, structure-aware spectral paths, exact symmetric polynomials/minors, and explicit dense-SVD gates |
| Invalid/adversarial inputs | Dimension/index/repetition/shape, physical-parameter, RNG-option, invalid-label, nonfinite map/state, non-CP recovery, implicit sparse densification, tolerance boundary, unsupported compatibility forms, exact-arithmetic overflow, compound boundary-order, and combinatorial predicate-guard cases covered |
| Explicit RNG safety | Six Tier B native/wrapper random constructors and the compatibility random `PauliChannel(rng, Q)` form are seeded/property tested; regressions verify the global stream is unchanged |
| MATLAB compatibility wrappers | Tier A, Tier B, supported Tier C, all 11 Tier D wrappers, six Tier E product entry points, four Tier E matrix entry points, and four structured matrix-predicate entry points pass locally; explicit partial statuses and structured-result differences remain documented |
| Benchmark smoke | 42 quick cases ran locally, including three product-analysis and three matrix-analysis cases; no regression threshold or comparative performance claim |
| Optional extension/load order | EntanglementDetection.jl 0.2.2 is integrated through a weak-dependency extension and an isolated child process; 125/125 focused assertions pass on Julia 1.12.6, including load order, lifecycle, failure, caller-state, timeout, IPC, and conservative-result checks. The compatible dependency graph resolves on Julia 1.11+; the six-job Julia 1.11/1.12 Linux/macOS/Windows matrix passed at predecessor commit `6bf8d61`, with an exact-candidate rerun pending |
| Optimization statuses | Pending |
| Quality and API consistency | Aqua passes 11/11, all 25 representative JET probes pass, and the provenance consistency gate passes over 214 public bindings; the JET set does not cover matrix predicates |
| Doctests/docs build | Strict Documenter builds and live tutorial examples passed locally on Julia 1.12.6 and 1.10.11; the split native reference is about 157 KiB, below the 200 KiB hard limit, and a Julia 1.12.6 CI-mode pretty-URL build passed |
| Julia 1.10/stable/nightly CI | Local Julia 1.10.11 and 1.12.6 core suites pass. The predecessor core workflow passed at `485b6a3`; the hardened exact-candidate workflow and nightly schedule still require remote runs. The optional extension's effective floor is Julia 1.11 |
| Linux/macOS/Windows CI | The optional Julia 1.11/1.12 six-job platform matrix passed at `6bf8d61`. The expanded core Julia 1.10 macOS/Windows jobs and all exact-candidate workflows remain pending |

## Historical commands and results (2026-07-29)

The historical 2,270-assertion core baseline was run on the exact staged source
tree committed as implementation milestone
`9b0d0d3b8177eded5b8743d250044d5428625b1b`. The release audit is based on
`485b6a3` and adds 99 regression assertions: 81 for constructor integrity,
post-construction plan/Kraus storage safety, and non-one-based-array rejection;
12 for the Brauer-state capped pre-allocation guards and compatibility
keywords; and six for post-construction channel-representation invariants.
It does not change the committed oracle artifacts. The environment is recorded
in `BUILD_ENVIRONMENT.md`:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

Results on Julia 1.12.6 and Julia 1.10.11: every package testset passed. The
core corpus comprises Tier A `345`, Tier B `577 + 214 + 133 + 30 = 954`,
Tier C `177`, Tier D measures/criteria
`6 + 25 + 30 + 11 + 15 + 11 + 38 + 32 = 168`, the native pipeline `68`,
coherence `54`, product `179 + 52`, matrix analysis `131 + 34`, and matrix
predicates `170 + 37`, for `2,369 / 2,369`. The executable tutorial gate adds
`48 / 48`, giving the integrated package result `2,417 / 2,417`.
The Julia 1.10 invocation warned that the ignored development
`Manifest.toml` had been resolved by Julia 1.12 and that project compatibility
had changed; `Pkg.test()` still resolved its temporary test environment and
passed. No root manifest is committed, so CI resolves against each configured
Julia version.

The exact tutorial scripts were also run through their standalone gate:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

Result: `48 passed / 48 total` on Julia 1.12.6 and Julia 1.10.11. The runner
executes `tutorials/subsystem_reductions.jl`,
`tutorials/local_channel_noise.jl`, and
`tutorials/entanglement_certificates.jl`, plus
`tutorials/separability_examples.jl` and
`tutorials/symmetric_sappt_witnesses.jl`; `Pkg.test()` includes the same
runner, and the strict documentation build evaluates the same calculations as
live examples.

The optional EntanglementDetection.jl environment and focused suite were run
with:

```sh
julia --project=test/extensions/entanglement_detection -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
```

Result on Julia 1.12.6: `125 passed / 125 total` against
EntanglementDetection.jl 0.2.2. The 125-assertion extension suite covers both
load orders, repeated loading, method-ambiguity checks, configuration and
density validation, caller RNG/stdout/logging/BLAS preservation, backend
failures, malformed IPC, bounded reads, timeouts, forced termination, and child
cleanup, including a post-launch output-stream close failure. The 68-assertion
core pipeline suite separately verifies dependency absence and its actionable
error.

Independent smoke calls covered `Float32`, `Float64`, `ComplexF32`,
`ComplexF64`, and multipartite input; those smoke calls are supplementary and
are not counted as separate assertions in the committed focused suite.

Searches run in a disposable child Julia process and exchange data through
trusted local Julia `Serialization`; this is an isolation boundary for audited
upstream side effects, not a security sandbox. Response and captured-output
reads are bounded. Tested timeout and injected-wait-error paths terminate and
reap their children; if bounded forced reaping ever fails, the adapter instead
returns `:termination_failed` and leaves the temporary directory for
process-exit cleanup.

The recorded ignored test manifest used the audited local checkout at
`dev/upstream/EntanglementDetection.jl`, pinned to
`5f60da1ceef6442acb669e10acc2fa47670bab06`. That checkout is not part of a
fresh clone; the setup command above and CI resolve the registered release
under the exact `0.2.2` compatibility bound. The runtime checks the loaded
package version but does not verify a source-tree hash, so a modified path
dependency retaining version 0.2.2 is trusted rather than detected.

The adapter exposes backend output only as candidate evidence. Whether
EntanglementDetection.jl suggests entangled, separable, or inconclusive, the
package-owned report remains `status = :unknown` and `certified = false`.
The isolated environment's compatible Ket 0.9 dependency makes Julia 1.11 the
effective resolver floor. This report cites the Julia 1.12.6 focused local run;
the Julia 1.11/1.12 Linux/macOS/Windows matrix passed at predecessor commit
`6bf8d61`, while an exact release-candidate rerun remains required.

The Tier C slice was also run directly through the package test environment:

```sh
julia --startup-file=no --project=. -e \
  'using Test, QuantumEntanglementTools; include("test/tier_c_channels_maps.jl")'
```

Result: `Tier C channel/map representations | 177 passed / 177 total`.
Coverage includes the 24 native public map bindings/types and all 11 Tier C
`MATLABCompat` wrappers. The supported native representation model permits
unequal input/output Hilbert-space dimensions, but the compatibility layer
does not yet implement QETLAB's two-sided left/right Kraus-cell form or
independent rectangular row/column operator spaces.

The development-only oracle generated 13 JSON fixtures from the pinned QETLAB
checkout with Octave 11.3.0. The committed source-free artifact has a sibling
SHA-256 record; its comparison run passed 52/52 native/wrapper assertions for
tensor operations, permutations/swaps, trace, partial transpose (including a
complex non-Hermitian input), realignment, and symmetric/antisymmetric
projectors. This is function-specific supplementary evidence, not a blanket
claim that Octave reproduces MATLAB/QETLAB behavior.

```sh
scripts/matlab_oracle/run_tier_a_oracle.sh --engine octave
julia --project=test/oracle -e \
  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=test/oracle test/oracle/compare_tier_a_oracle.jl \
  test/oracle/fixtures/tier_a_octave_11_3_qetlab_d858961.json
```

Result: `QETLAB Tier A differential fixture | 52 passed / 52 total`. MATLAB,
CVX, and a CVX solver were not detected; the fixture metadata records that
absence.

The analogous deterministic Tier B artifact was generated from the same pinned
QETLAB checkout with Octave 11.3.0:

```sh
scripts/matlab_oracle/run_tier_b_oracle.sh --engine octave
julia --project=test/oracle test/oracle/compare_tier_b_oracle.jl \
  test/oracle/fixtures/tier_b_octave_11_3_qetlab_d858961.json
```

Result: `QETLAB Tier B differential fixture | 72 passed / 72 total` across 18
fixtures. Both the native and `MATLABCompat` path are checked for each fixture.
The committed artifact is
`test/oracle/fixtures/tier_b_octave_11_3_qetlab_d858961.json`, with SHA-256
`61f9e585b6a608f593272a4c62819cc5865da0c4a7e8d79113a9123087c840e5`.
Randomized functions are intentionally absent because equal numeric seeds do
not imply equal streams across Julia and MATLAB-family engines.

The Tier C artifact was generated from the same pinned checkout with Octave
11.3.0:

```sh
scripts/matlab_oracle/run_tier_c_oracle.sh --engine octave
julia --project=test/oracle test/oracle/compare_tier_c_oracle.jl \
  test/oracle/fixtures/tier_c_octave_11_3_qetlab_d858961.json
```

Result: `QETLAB Tier C differential fixture | 28 passed / 28 total` across 7
fixtures. Both native and `MATLABCompat` results are checked for depolarizing,
dephasing, and Pauli channels, the Choi and reduction maps, whole-map
application, and partial-map application. The committed artifact is
`test/oracle/fixtures/tier_c_octave_11_3_qetlab_d858961.json`, with SHA-256
`46b31802d97ee8366da163d7c723408ac708c94325aa96c83a34307355a15c03`.
Representation round trips and canonical Kraus recovery use analytic/property
tests because Kraus bases are not unique. The fixture is supplemental
function-specific evidence; MATLAB was not run.

The Tier D measures/criteria slice and the project-native pipeline were also
run directly through the package test environment:

```sh
julia --startup-file=no --project=. -e \
  'using Test, QuantumEntanglementTools; include("test/tier_d_measures_criteria.jl")'
julia --startup-file=no --project=. -e \
  'using Test, QuantumEntanglementTools; include("test/tier_d_entanglement_pipeline.jl")'
```

Results: `Tier D scalar measures and criteria | 168 passed / 168 total` and
`Tier D entanglement pipeline | 68 passed / 68 total`. The first count includes
all 11 Tier D `MATLABCompat` entry points. The pipeline treats PPT, realignment,
and reduction passes as inconclusive unless a separately stated theorem is
sufficient. It certifies PPT separability only in bipartite `2×2` and `2×3`.
For pure vectors, a trailing Schmidt coefficient above tolerance certifies
entanglement; separability requires the computed trailing coefficients to be
exactly zero, while a tolerance-defined rank-one result with nonzero trailing
coefficients remains `unknown`. The pipeline does not claim full QETLAB
`IsSeparable`.

The Tier D artifact was generated from the same pinned checkout with Octave
11.3.0:

```sh
scripts/matlab_oracle/run_tier_d_oracle.sh --engine octave
julia --project=test/oracle test/oracle/compare_tier_d_oracle.jl \
  test/oracle/fixtures/tier_d_octave_11_3_qetlab_d858961.json
```

Result: `QETLAB Tier D differential fixture | 34 passed / 34 total` across 13
fixtures. The committed artifact is
`test/oracle/fixtures/tier_d_octave_11_3_qetlab_d858961.json`, with SHA-256
`ad0cdc45077390fc1eb736fc7c7ff1ec41696c796a508b536774cb6e0020160a`.
It covers trace/Schatten/Ky Fan norms, purity, the von Neumann `Entropy`
`ALPHA=1` branch, unsquared root fidelity, negativity, Schmidt coefficients
and rank, concurrence, two PPT outcomes, and the realignment trace norm.
PPT fixture booleans are compared to the corresponding tri-state meaning
without coercing `CriterionResult` or asserting separability. Schmidt vectors
are validated by phase-independent reconstruction rather than entrywise
comparison. This is supplemental function-specific evidence; MATLAB was not
run.

The Tier E focused suites and supplemental artifacts were run with:

```sh
julia --project=. -e \
  'using QuantumEntanglementTools, Test; include("test/tier_e_product_analysis.jl")'
julia --project=. -e \
  'using QuantumEntanglementTools, Test; include("test/tier_e_product_compat.jl")'
julia --project=test/oracle test/oracle/compare_tier_e_product_oracle.jl \
  test/oracle/fixtures/tier_e_product_octave_11_3_qetlab_d858961.json
julia --project=test/oracle test/oracle/compare_tier_e_coherence_oracle.jl \
  test/oracle/fixtures/tier_e_coherence_octave_11_3_qetlab_d858961.json
julia --project=. -e \
  'using QuantumEntanglementTools, Test; include("test/tier_e_matrix_analysis.jl"); include("test/tier_e_matrix_analysis_compat.jl")'
julia --project=test/oracle test/oracle/compare_tier_e_matrix_analysis_oracle.jl \
  test/oracle/fixtures/tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json
```

Product analysis passes 179/179 native and 52/52 compatibility assertions.
Its 14-fixture artifact passes 68/68 with SHA-256
`ab6414c1a684141db74782616d4c18e79c8e6039aad53a695d8c723eed598d85`.
Matrix analysis passes 131/131 native and 34/34 compatibility assertions. Its
22-fixture artifact passes 59/59 with SHA-256
`e37685c262ce5982d10dd705cef8c172d49d9c55c89a0a67d4de729a5068f540`;
17 fixtures are agreements and five preserve reviewed QETLAB discrepancies
without turning them into native expected values. The earlier Tier E coherence
slice passes 54/54 locally and 25/25 against six fixtures (SHA-256
`11bcaaee88fac8a595e9a4eff164432dbaa4e141cdd26554da2522e22981811a`).
MATLAB was not run.

The matrix-predicate slice was run both through the integrated package suite
and its focused files:

```sh
julia --startup-file=no --project=. -e \
  'using QuantumEntanglementTools, Test; include("test/tier_e_matrix_predicates.jl"); include("test/tier_e_matrix_predicates_compat.jl")'
julia --startup-file=no --project=. scripts/validate_matrix_predicates.jl
```

Results: 170/170 native and 37/37 compatibility assertions. An independent
randomized check compared PSD outcomes with eigenspectra and all-minor outcomes
with compound-matrix determinants for 130/130 assertions. Exact
integer/rational, `Float32`, `Float64`, `BigFloat`, complex nonsingularity,
structured diagonal, explicit sparse-densification, tolerance-boundary, witness,
and combinatorial-guard paths are covered. No MATLAB-family predicate fixture
has been generated, and the `IsPSD` inventory mapping remains partial because
the pinned CVX symbolic branch is omitted.

```sh
julia --project=quality quality/run_quality.jl
```

Result: Aqua passes 11/11 and all 25 representative JET inference probes pass,
including four product-analysis calls, four matrix-analysis calls, and the
package-owned `EntanglementDetectionSearch` configuration constructor. The JET
smoke does not load the optional dependency and does not currently include the
matrix-predicate slice.

```sh
julia --startup-file=no --project=. scripts/check_public_api.jl
```

Result: the public-API/provenance consistency gate passes over 214 public
bindings. This validates ledger consistency, not behavioral parity for every
binding.

```sh
julia --project=benchmark benchmark/benchmarks.jl --quick --no-save
```

Result: all 42 quick benchmark cases completed. The product and matrix slices
each add exactly three representative cases. The matrix cases cover two dense
`32×32` matrices compared through singular values, an `8×8` order-three
compound, and a `16×16`
order-two additive compound. This is a local allocation/runtime smoke run, not
a saved regression baseline or a comparative performance claim.

```sh
julia --project=docs docs/make.jl
```

Results on Julia 1.12.6 and Julia 1.10.11: Documenter completed doctests,
cross-references, strict exported-doc checks, live tutorial examples, and HTML
rendering without errors. The API reference is split into native and
compatibility pages; the native page is about 157 KiB, below the 200 KiB hard
limit, and emits only a non-failing size warning. A Julia 1.12.6 CI-mode build
also produced the expected pretty URLs and local code-generator assets.

Release-specific local gates were also run:

```sh
julia --startup-file=no --project=. scripts/check_release.jl --allow-dirty
julia +1.10 --startup-file=no --project=. scripts/check_release.jl --allow-dirty
julia --startup-file=no --project=quality quality/format.jl
julia --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
julia --startup-file=no --project=. scripts/check_public_api.jl
```

The dirty-worktree release preflight passed on Julia 1.12.6 and Julia 1.10.11.
CFFConvert 2.0.0 validated `CITATION.cff` against schema 1.2.0; the offline
upstream audit passed 10/10 pin and license checks; the inventory remained at
163 files, 503 edges, and zero detected cycles; and the public API remained at
214 provenance-covered bindings. The exact Git archive of candidate code
commit `0b63359159e1c0c1527c8753f78b61940701eb25` passed the release check and
fresh-depot smoke on Julia 1.12.6 and Julia 1.10.11; both runs produced tar
SHA-256 `226302eb31191306111fad3aae293548beceec9c5c45640ddda1a110b22c2396`.
The same two archive gates were rerun successfully on the final
documentation-only evidence `HEAD`.

## Acceptance rule

Update this report only with commands and artifacts actually run. For each
validated function or group, record:

- repository commit and dirty state;
- environment/manifest;
- test and oracle command;
- input families, sizes, types, and sparsity;
- comparison definition and tolerance rationale;
- pass/fail counts and unresolved discrepancies;
- links to fixtures without copying restricted source or assets.

An upstream match alone is not sufficient when independent mathematical
evidence is practical. A discrepancy must be investigated, not hidden by a
broader tolerance.

## Known gaps

The local static implementation ledger is complete, but broader release
evidence is not:

- MATLAB/CVX was not run, and Octave fixtures do not establish general MATLAB
  equivalence.
- Remote Linux/macOS/Windows CI and accepted Codecov ingestion have not run on
  an exact commit containing this convergence work.
- The 120-case quick benchmark is a smoke suite, not a reviewed regression
  baseline or performance comparison.
- The optional EntanglementDetection environment has an effective Julia 1.11
  resolver floor and treats the trusted child worker as local IPC, not a
  security boundary.
- The worktree is uncommitted, so exact archive and fresh-depot release gates
  cannot yet identify the final source state.
- The maintainer's non-delegable mathematical, API, licensing, provenance, and
  generated-change review remains open.
