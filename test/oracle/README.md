# Development-only QETLAB oracle

This harness generates small numerical fixtures from the exact QETLAB revision
recorded in `UpstreamManifest.toml`. MATLAB is optional and is never a package
or public-CI dependency. Octave 11.3.0 generated the committed Tier A--Tier E
artifacts described below, but their results are function-specific supplemental
evidence only—not a claim that Octave is generally equivalent to MATLAB.

The JSON schema stores every array as its row count, column count, and separate
column-major real and imaginary vectors. That representation preserves complex
arrays without relying on engine-specific JSON handling. Metadata records the
engine, engine version, platform, QETLAB commit, and any detected CVX version
and solver. A sibling `.sha256` file protects each generated document.

`fixtures/tier_a_octave_11_3_qetlab_d858961.json` is the small, source-free
supplemental fixture generated during the 2026-07-28 audit. The comparison
script uses a freshly generated local fixture when one exists, otherwise it
checks this committed artifact.

The analogous Tier B files cover 18 deterministic operator/state families
(with separate `3 x 3` and `2 x 4` Horodecki fixtures). Randomized functions
are validated by seeded analytic/property tests instead of comparing unrelated
random streams across languages.

The Tier C fixture covers five channel/map Choi constructors plus whole-map and
partial-map application. Representation round trips and canonical Kraus
recovery remain independently tested in Julia because Kraus bases are not
unique.

The Tier D fixture covers 13 deterministic norm, state-measure, Schmidt,
concurrence, PPT, and realignment cases. It checks 34 native assertions.
Phase-ambiguous Schmidt vectors are validated by analytic reconstruction in
the Julia test suite rather than compared entry-by-entry. QETLAB `IsPPT`
booleans are compared only with the corresponding mathematical meaning of the
native tri-state result; a passed necessary condition is never relabeled as a
separability certificate.

The Tier E coherence fixture covers six deterministic `l1`, relative-entropy,
and coherence-rank cases and checks 25 assertions. Two rank fixtures record the
pinned QETLAB implementation's zero-counting bug; the Julia assertions
deliberately follow the documented nonzero-coefficient definition instead.
Its SHA-256 is
`11bcaaee88fac8a595e9a4eff164432dbaa4e141cdd26554da2522e22981811a`.

The Tier E product-analysis fixture covers 14 deterministic operator-Schmidt,
product-vector/operator, entanglement-of-formation, and separable-ball cases.
It checks 68 native, compatibility, reconstruction, and discrepancy assertions.
Product Booleans are compared only away from tolerance boundaries, and failure
of the sufficient separable-ball condition is never treated as entanglement.
The fixture records pinned rectangular/Hermitian operator-Schmidt failures and
the zero-concurrence `EntFormation` `NaN` without reproducing them in Julia. Its
SHA-256 is
`ab6414c1a684141db74782616d4c18e79c8e6039aad53a695d8c723eed598d85`.

The Tier E matrix-analysis fixture covers 22 deterministic majorization,
elementary-symmetric-polynomial, compound, and additive-compound cases. Its 17
agreement fixtures and five reviewed discrepancy fixtures pass 59 assertions,
including checks that the compatibility layer preserves reviewed weak
majorization and compound boundary behavior without making those behaviors
native expected values. Its SHA-256 is
`e37685c262ce5982d10dd705cef8c172d49d9c55c89a0a67d4de729a5068f540`.

The matrix-predicate slice passes 170 native and 37 compatibility assertions in
the package test suite. It has no MATLAB-family oracle artifact; do not infer
upstream parity from its analytic, exact, property, or randomized checks.

The parallel-repetition fixture records the one-copy and two-copy branches of
the pinned `ParallelRepetition` entry point for a rectangular four-axis game.
The full tensor shape and column-major flattened coefficients verify the
packed-copy order. This is supplemental Octave evidence; the allocation,
overflow, exact-type, and invalid-input contracts are independently tested in
Julia. Its SHA-256 is
`fc2d9df8cfa32827a9461b543f5395cb768a2e7391d73545986a639b1eb81774`.

The entangled-subspace fixture covers equal and unequal local dimensions,
nonmaximal prefix selection, and the `r=2` branch of the pinned
`EntangledSubspace` construction. Exact coefficient and support comparison
also verifies the documented tensor-factor order. This is supplemental Octave
evidence; exact arithmetic, Schmidt-rank guarantees, allocation guards, and
invalid-input behavior are independently tested in Julia. Its SHA-256 is
`d71e8001d9df7844b468e4fa257bf7fe1513d7d974d46c7914048df11ab82d69`.

The entangling-gate fixture records pinned flags for identity, swap, CNOT, and
controlled-Z gates plus the two returned upstream witness candidates. The CNOT
candidate is a valid normalized product input with entangled output. The
controlled-Z candidate exposes a pinned defect: it is left unnormalized and
its output remains product despite the true flag. The Julia implementation
retains the flag but supplies a validated four-support phase-grid witness
instead. Its SHA-256 is
`35f8d5a20171c5da8a986311e21db95cd0d85cb53a0d39218d142e3c2472d617`.

The IsUPB fixture covers the two-qutrit Tiles and three-qubit Shifts UPBs,
an extendible family, a complete product basis, a nonorthogonal unextendible
family, and a complex extension witness. Its 33 assertions check conclusive
native and compatibility behavior, exact witness orthogonality, and three
reviewed pinned discrepancies: missing full-definition checks and a
nonconjugating complex nullspace. Its SHA-256 is
`56f2f9bd0a3fa7d18ed31f45f896926f5e2c6df3b2400b59b0295f4fae0119e1`.

The minimum-UPB-size fixture covers every branch of the pinned theorem and
exception table through 11 exact dimension families, plus the reviewed
`(2,3,4)` unknown/error route. The native API returns that last case as an
explicit structured `:unknown` instead of fabricating a minimum.
Its SHA-256 is
`8ed38bab7ad6c1222de6174d964a9df53cd3af471c3192637af752d318584b26`.

The UPB-catalog fixture compares 17 deterministic named and dimension-driven
families with the pinned source up to independent phase on every local vector.
Its 137 parity/invariant assertions are supplemented by eight assertions that
retain the pinned nonorthogonal `John2^4k(8)` reshape-order result and validate
the corrected native graph construction. Randomized families use seeded
theorem/property tests in Julia instead of comparing unrelated random streams.
The fixture metadata records both reviewed source-file hashes. Its SHA-256 is
`cac29476f30f3a901e7c7ae36877f54097ffe1c6f91d049180b2bb1113585a93`.

The multipartite-Werner fixture records the pinned three-party loop-overwrite
defect, the valid one-entry-vector bipartite route, and the invalid parameter
length. The Julia native and compatibility APIs implement the documented sum
over every permutation, enforce Hermiticity and positivity, and deliberately
do not reproduce the last-parameter-only output. Its SHA-256 is
`9e19082444c97468756b4a7ba60be29885f6eb6928c5c7530bccf3063047b427`.

The pure-state `k`-coherence robustness fixture covers eight normalized,
sorted, nonnegative theorem-domain calls and records both pinned outputs.
Three correction fixtures retain the unsafe pinned behavior for an unsorted
vector, a complex phased vector, and a nonnormalized vector. The Julia native
and compatibility APIs sort magnitudes, validate normalization, and return a
real phase-invariant theorem value; the discrepancy values are never used as
native expectations. Its SHA-256 is
`551f529a86960a43e8765e2934007a1a640b86f8327a2f19f79860d9d5ce51a6`.

The Twirl fixture covers Werner, isotropic, real, and Pauli two-copy
projections plus three-copy Werner and real projections of general complex
operators. Native, compatibility, and idempotence checks pass 18 assertions.
Strict invalid-copy/dimension behavior is independently tested because the
pinned routine under-validates values below two. Its SHA-256 is
`acefd7ad638fa8e0eb0e30812e975f53c864de85e8abc990c88c832d00dc9e77`.

The random-superoperator fixture records five seeded QETLAB Choi matrices
spanning unconstrained, trace-preserving, unital, equal-dimensional
bistochastic, real, complex, rank-controlled, and unequal proportional-output
branches. The comparator uses them as property evidence rather than comparing
unrelated random streams entrywise, and records the permissive-flag,
oversized-rank, and unequal-unitality upstream discrepancies. Octave 11.3
cannot execute the pinned complex constrained branch because of an upstream
`PartialTrace`/`cellfun` incompatibility, so complex coverage is unconstrained
and constrained fixtures are real. The 51 pinned, 45 native, and 10
compatibility assertions pass on Julia 1.12.6 and 1.10.11. Its SHA-256 is
`2bf98f8cbbc188be160755d8e866634658de109a0a92b402a353287c1ad9ca7a`.

The absolute-PPT fixture records the pinned $p=2$ and $p=3$ Hildebrand
matrices, the first three `LIM`-capped $p=4$ matrices, the
`ESC_IF_NPOS` early exit, and five legacy `IsAbsPPT` outputs. Its 28 integrated
native and compatibility assertions also retain the two pinned scalar-DIM
conventions and the floating boundary that the native API intentionally maps
to `unknown`. Octave 11.3.0 generated the committed source-free artifact from
the pinned revision; MATLAB was not run. Its SHA-256 is
`17b825d6e15fae1f2391162137255603653de245fea39661c18ba1db65e74f6c`.

The symmetric-extension fixture records five solver-free outer-extension
decisions, three exact `jacobi_poly` coefficient vectors used by the inner
hierarchy, and one seeded upstream random-PPT property record. The Julia
comparison checks extension decisions and recurrence coefficients directly,
but checks the random state only through PSD, unit-trace, PPT, and rank
properties because the languages do not share an RNG. The private Julia
Jacobi recurrence supersedes the upstream helper without exporting it.
Octave 11.3.0 generated the artifact; MATLAB and CVX were not run. Its
SHA-256 is
`572dc3d846d1adb0b2ee64f9ee5a3b7a6903eb86d56927ff76e5e8fa520551f0`.

The state-discrimination fixture covers one-state, two-state pure and mixed,
and mutually orthogonal multi-state solver-free branches. It records both the
pinned scalar and the objective attained by the pinned POVM for unequal pure
priors, where the upstream closed formula can exceed one. The Julia comparator
uses the correct Helstrom value, validates every returned measurement, and
rejects pinned silent state normalization and negative-prior acceptance. Its
119 assertions pass against the source-free Octave 11.3.0 artifact with
SHA-256
`88f343739fc3e6618760cee454ad738dd63fb4bfa154aadb80e37393121691fd`.

The coherence-optimization fixture covers deterministic conclusive,
inconclusive, and invalid-input routes for the two coherence criteria and the
three optimization problems. The 51 comparator assertions preserve the pinned
`IsAbskIncoh` bandwidth-helper defect as discrepancy evidence while validating
the corrected bounded graph-bandwidth search used by the Julia API. Solver
results are checked separately in the optional JuMP extension tests, including
termination, primal, and dual statuses. Octave 11.3.0 generated the source-free
artifact; its SHA-256 is
`2b44fc70b33183c01c8c8dad4c71e122ed182fa776b6e95cc15ac693eb499dad`.

The channel-optimization fixture covers 11 deterministic solver-free
`DiamondNorm`, `CBNorm`, `ChannelDistinguishability`, and
`MaximumOutputFidelity` cases. Its 122 assertions validate native analytic
certificates and retain three pinned discrepancies: the missing affine
channel Holevo--Helstrom conversion, accepted negative priors, and
unequal-Kraus-rank truncation in maximum output fidelity. General SDP evidence
comes from separate Hypatia/SCS extension tests because CVX was not available
to Octave. The source-free fixture SHA-256 is
`4855dae732c48b4e8e31e357e6ae1fcf392d62572bf6e6092d278ed30fe327c7`.

The S(`k`)-norm fixture covers four numeric `kpNormDual` plateau/endpoint
cases, two solver-free exact `SkOperatorNorm` cases, and four
`IsBlockPositive` certificate or boundary cases. Its 26 assertions use only
pinned branches that do not invoke CVX or randomized iteration; the affine
dual-norm atom and SDP relaxations are validated independently with Hypatia
and SCS. The source-free Octave 11.3.0 artifact has SHA-256
`c21811376362955cb96b179f3a0eb1f1c37aa0db0e5445cff4eabbc02c1deabe`.

The copositivity/clique fixture covers exact direct certificates and
solver-free polynomial branches for `IsCopositive` and `CliqueNumber`.
Its 80 assertions compare pinned outputs while separately enforcing the
native tri-state and certified-graph-bound contracts. Hypatia/SCS evidence is
kept in the optional extension suite. The source-free Octave 11.3.0 artifact
has SHA-256
`31806fca115a245a4669e72a7f5bcafcf04b4532102f51c2acd2dada33d81a55`.

From the repository root:

```sh
scripts/matlab_oracle/run_tier_a_oracle.sh --engine auto
julia --project=test/oracle -e \
  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=test/oracle test/oracle/compare_tier_a_oracle.jl \
  test/oracle/fixtures/tier_a_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_b_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_b_oracle.jl \
  test/oracle/fixtures/tier_b_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_c_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_c_oracle.jl \
  test/oracle/fixtures/tier_c_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_d_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_d_oracle.jl \
  test/oracle/fixtures/tier_d_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_e_coherence_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_coherence_oracle.jl \
  test/oracle/fixtures/tier_e_coherence_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_e_product_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_product_oracle.jl \
  test/oracle/fixtures/tier_e_product_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_tier_e_matrix_analysis_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_matrix_analysis_oracle.jl \
  test/oracle/fixtures/tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_parallel_repetition_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_parallel_repetition_oracle.jl \
  test/oracle/fixtures/parallel_repetition_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_entangled_subspace_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_entangled_subspace_oracle.jl \
  test/oracle/fixtures/entangled_subspace_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_entangling_gate_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_entangling_gate_oracle.jl \
  test/oracle/fixtures/entangling_gate_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_is_upb_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_is_upb_oracle.jl \
  test/oracle/fixtures/is_upb_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_minimum_upb_size_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_minimum_upb_size_oracle.jl \
  test/oracle/fixtures/minimum_upb_size_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_upb_catalog_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_upb_catalog_oracle.jl \
  test/oracle/fixtures/upb_catalog_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_multipartite_werner_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_multipartite_werner_oracle.jl \
  test/oracle/fixtures/multipartite_werner_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_robk_coherence_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_robk_coherence_oracle.jl \
  test/oracle/fixtures/robk_coherence_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_twirl_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_twirl_oracle.jl \
  test/oracle/fixtures/twirl_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_random_superoperator_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_random_superoperator_oracle.jl \
  test/oracle/fixtures/random_superoperator_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_absolute_ppt_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_absolute_ppt_oracle.jl \
  test/oracle/fixtures/absolute_ppt_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_symmetric_extensions_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_symmetric_extensions_oracle.jl \
  test/oracle/fixtures/symmetric_extensions_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_state_discrimination_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_state_discrimination_oracle.jl \
  test/oracle/fixtures/state_discrimination_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_coherence_optimization_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_coherence_optimization_oracle.jl \
  test/oracle/fixtures/coherence_optimization_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_channel_optimization_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_channel_optimization_oracle.jl \
  test/oracle/fixtures/channel_optimization_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_sk_norms_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_sk_norms_oracle.jl \
  test/oracle/fixtures/sk_norms_octave_11_3_qetlab_d858961.json
scripts/matlab_oracle/run_copositivity_clique_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_copositivity_clique_oracle.jl \
  test/oracle/fixtures/copositivity_clique_octave_11_3_qetlab_d858961.json
```

The explicit fixture arguments above verify the committed artifacts even when
an ignored `test/oracle/generated/` file exists. Omit the argument only when
intentionally comparing a freshly generated fixture.

Use `--engine matlab` for authoritative MATLAB/QETLAB differential evidence or
`--engine octave` for the explicitly limited supplemental check. `--qetlab`
selects another checkout, but generation fails unless its Git revision is the
pinned commit. `--output` can select a fixture path.

Only the generated numerical values and metadata may be committed. Do not copy
QETLAB source into this directory, and never put MATLAB, CVX, or solver binaries
in the repository.
