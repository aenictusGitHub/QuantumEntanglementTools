# Session handoff

Snapshot date: 2026-07-28. Refresh this file at the end of each substantive
session; do not assume the worktree details below remain current after parallel
work continues.

## Repository state

- Branch: `main`.
- Base commit: none; the branch was unborn when this scaffold was created.
- Remote: none configured at that snapshot.
- Package: `QuantumEntanglementTools`, UUID
  `45675e5b-5c8b-4983-b92d-4c3725d56c4e`, version `0.1.0`.
- Package author metadata: `John MARTIN <jmartin@uliege.be>`.
- Worktree: intentionally uncommitted initial project work was present; inspect
  `git status` before proceeding and do not discard it.

## Evidence secured

- QETLAB inspected at
  `d8589610f00cff106537268dee2e2a1153f3a601` on 2026-07-28.
- QETLAB root BSD-2 license and Bruno Luong helper BSD-2 notices have recorded
  hashes in `docs/LEGAL.md`.
- EntanglementDetection.jl 0.2.2/tag SHA and MIT license hash are recorded.
  Its observed global RNG, stdout/logging, and BLAS-thread side effects are
  explicit M6 gates.
- The separately supplied QUBIT4MATLAB license is preserved verbatim.
- The author-hosted QUBIT4MATLAB v6.5 archive and bundled license are hashed.
  Their terms materially conflict, and no `.m` implementation files were
  inspected.
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

Source-free Octave/QETLAB fixtures pass Tier A 52/52, Tier B 72/72, Tier C
28/28, Tier D 34/34, Tier E coherence 25/25, Tier E product 68/68, and Tier E
matrix 59/59. The matrix artifact has 22 fixtures, including five reviewed
discrepancies for weak majorization and compound boundary behavior; its SHA-256
is `e37685c262ce5982d10dd705cef8c172d49d9c55c89a0a67d4de729a5068f540`.
Octave is not treated as generally MATLAB-equivalent. Matrix predicates have
analytic/property evidence, including 130/130 independent randomized
eigenspectrum/minor cross-checks, but no MATLAB-family oracle.

All 42 quick benchmark cases completed, including exactly three product and
three matrix cases. Aqua passes 11/11 and all 24 representative JET probes
pass; those probes do not cover matrix predicates. The integrated 2,270-assertion
corpus passes on Julia 1.12.6 and Julia 1.10.11. Strict Documenter, formatter,
inventory (163 files/503 edges/zero cycles), and public-API/provenance
(212 bindings) gates pass. Every committed oracle comparator passes, including
the 59-assertion matrix comparator.

## Blocking issue

QUBIT4MATLAB source-informed work is blocked. The bundled v6.5 license differs
materially from the separately supplied license. Required action: obtain
authoritative human/rightsholder clarification and preserve it as evidence.
Until then, do not inspect archive implementation `.m` files, do not claim the
supplied license governs the archive, and do not create a derived port.

## Next safe work

1. Review the remaining 87 inventory rows and manually correct classifications.
2. Keep `UpstreamManifest.toml`, `PROVENANCE.toml`, status overlays, and
   generated inventory synchronized as Tier C and later slices land.
3. Run the remote core/docs matrix on Julia 1.10/stable and supported platforms.
4. Add authoritative MATLAB differential checks when MATLAB is available.
5. Design the optional EntanglementDetection adapter around observed global
   side effects; do not expose it as safe before fresh-process regression tests.
6. Refresh this handoff with resulting commits, failures, and uncommitted files.

## Commands to rerun

```sh
git status --short
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia +1.10 --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
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
