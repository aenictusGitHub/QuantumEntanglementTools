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
