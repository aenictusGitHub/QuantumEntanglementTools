# Contributing

Thank you for helping build `QuantumEntanglementTools`. The `0.1.x` API is
experimental, so correctness, traceable provenance, and explicit conventions
take priority over API breadth.

## Before starting

Read `AGENTS.md`, `docs/PORTING_STATUS.md`, `docs/SESSION_HANDOFF.md`, and
`docs/src/conventions.md`. Check the upstream inventory and provenance records
before selecting a function. Open an issue before making a broad convention,
dependency, or public-API change when an issue tracker is available.

## Development setup

Install Julia 1.10 or later, then instantiate and test from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The current full corpus contains 2,417 assertions—2,369 core plus 48
executable-tutorial assertions—and has passed locally on Julia 1.12.6 and Julia
1.10.11. Those runs are development evidence, not a substitute for
supported-platform CI or MATLAB validation.

Run the tutorials independently with:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

Build the documentation with:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Run the 42-case non-recording quick benchmark smoke suite with:

```sh
julia --project=benchmark benchmark/benchmarks.jl --quick --no-save
```

Omit `--no-save` only when you intend to keep a local raw artifact under the
ignored `benchmark/results/local/` directory. Do not present planned commands
as passing checks or the smoke suite as a performance comparison.

Set up and run the optional EntanglementDetection.jl 0.2.2 extension tests from
the repository root with:

```sh
julia --project=test/extensions/entanglement_detection -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
```

This fresh-clone path resolves the registered release under the exact `0.2.2`
compatibility bound. Maintainers may instead develop the ignored audited
checkout at `dev/upstream/EntanglementDetection.jl` after verifying its commit
against `UpstreamManifest.toml`. The runtime adapter enforces the package
version, not a source-tree hash, so controlled validation environments are
responsible for source integrity. Run this optional environment on Julia 1.11
or later; its Ket 0.9 dependency does not currently resolve on the core
package's Julia 1.10 minimum. The dedicated extension suite passes 125/125
locally on Julia 1.12.6. The six-job Julia 1.11/1.12 Linux/macOS/Windows
workflow passed at predecessor commit `6bf8d61`; an exact release-candidate
rerun remains required.

Run the package-quality and ledger checks with:

```sh
julia --project=quality quality/run_quality.jl
julia --project=quality quality/format.jl
julia --project=. scripts/build_upstream_inventory.jl --check
julia --project=. scripts/check_public_api.jl
julia --project=. scripts/validate_matrix_predicates.jl
```

The current quality run reports Aqua 11/11 and 25 representative JET probes.
Those JET probes do not cover the matrix-predicate functions. The inventory
currently has 76 manually reviewed rows: 63 implemented, 11 partial, two
deferred, and 87 pending.

Run the optional development oracle against a pinned QETLAB checkout with:

```sh
scripts/matlab_oracle/run_tier_a_oracle.sh --engine auto
julia --project=test/oracle -e \
  'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=test/oracle test/oracle/compare_tier_a_oracle.jl
scripts/matlab_oracle/run_tier_b_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_b_oracle.jl
scripts/matlab_oracle/run_tier_c_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_c_oracle.jl
scripts/matlab_oracle/run_tier_d_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_d_oracle.jl
scripts/matlab_oracle/run_tier_e_coherence_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_coherence_oracle.jl
scripts/matlab_oracle/run_tier_e_product_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_product_oracle.jl
scripts/matlab_oracle/run_tier_e_matrix_analysis_oracle.sh --engine auto
julia --project=test/oracle test/oracle/compare_tier_e_matrix_analysis_oracle.jl
```

MATLAB is preferred when present. Octave results are labeled supplemental and
function-specific. The committed matrix-analysis artifact passes 59 assertions
over 22 fixtures. Matrix predicates have 170 native and 37 compatibility
assertions but no MATLAB-family oracle.

## Completing a public function

A contribution adding a public operation should include:

- an upstream or mathematical specification and provenance classification;
- a reviewed Julia-native signature and any compatibility mapping;
- documented dimensions, ordering, normalization, tolerances, errors, and
  complexity;
- type-generic and sparse behavior, or an explicit justified limitation;
- analytic, property, invalid-input, and independent/differential tests as
  applicable;
- runnable documentation and an inventory status update;
- benchmarks for hot operations, without unsubstantiated speed claims.

Never convert solver failure, a heuristic result, or an unmet necessary
condition into a definitive mathematical conclusion.

## Style and scope

- Follow Julia base style and use lowercase `snake_case` for native APIs.
- Prefer ordinary array abstractions and avoid type piracy.
- Keep optional integrations in Julia package extensions.
- Accept an explicit RNG in randomized APIs.
- Keep pull requests focused; preserve unrelated work.
- Add dependencies only with an architecture and license rationale.
- Do not paste upstream prose or code without recording permission and
  attribution.

## Changes and review

Describe the mathematical behavior, provenance, test commands, numerical
tolerances, allocation/densification implications, and known limitations in the
change summary. Update `CHANGELOG.md` under **Unreleased** for user-visible
changes.

By contributing, you agree that your contribution may be distributed under the
project's BSD 3-Clause License and that you have the right to submit it. Report
security-sensitive issues according to `SECURITY.md`, not in a public ticket.
