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

<!-- qetlab-current-claims: begin -->

The pinned inventory revision is
`d8589610f00cff106537268dee2e2a1153f3a601`. Its strict static ledger reports
127/127 public rows are verified with the required final status, 36/36 internal
helpers have terminal dispositions, the completion queue contains 0 public
rows, 0 required internal helpers remain, and 0 static completion failures.
The public API has 477 public bindings (347 native/module and 130
`MATLABCompat`) with matching provenance entries.

The 9,484-assertion full package suite passed 9,484/9,484: 9,400 core
assertions plus 84 executable-tutorial assertions, including 36 for the
two QETLAB-introduction workflows. It passes on Julia 1.12.6 and the installed
Julia 1.10.11. The full optional JuMP suite passed 847/847 on both Julia lines. The
EntanglementDetection.jl extension passed 141/141 focused assertions on the
current compatible Julia. All 114 declared quick benchmark cases completed
without failure on the current uncommitted release-hardening worktree based on
`8b2fcbaf`. This is local smoke evidence only; targeted paired observations
remain diagnostics, not a stable comparative-performance baseline.

These runs and static checks are development evidence, not QETLAB/MATLAB
parity, supported-platform remote CI, comparative performance, API stability,
release approval, or non-delegable human review. No version has been tagged or
published.

<!-- qetlab-current-claims: end -->

Run the tutorials independently with:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

Build the documentation with:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The benchmark source declares 114 quick cases. Run its non-recording smoke
command with:

```sh
julia --project=benchmark benchmark/benchmarks.jl --quick --no-save
```

Omit `--no-save` only when you intend to keep a local raw artifact under the
ignored `benchmark/results/local/` directory. Do not present planned commands
as passing checks, the declared count as execution evidence, or the smoke suite
as a performance comparison.

Set up and run the optional JuMP extension suite from the repository root with:

```sh
julia --project=test/extensions/jump_optimization -e '
    using Pkg
    Pkg.develop(PackageSpec(path=pwd()))
    Pkg.instantiate()
'
julia --startup-file=no --project=test/extensions/jump_optimization \
  test/extensions/jump_optimization/runtests.jl
```

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
package's Julia 1.10 minimum. The dedicated extension suite passes 141/141
locally on Julia 1.12.6. The six-job Julia 1.11/1.12 Linux/macOS/Windows
workflow passed at predecessor commit `6bf8d61`; a rerun on the exact current
commit remains required.

Run the package-quality and ledger checks with:

```sh
julia --project=quality quality/run_quality.jl
julia --project=quality quality/format.jl
julia --project=. scripts/build_upstream_inventory.jl --check
julia --project=. scripts/check_public_api.jl
julia --project=. scripts/validate_matrix_predicates.jl
```

The strict static completion ledger is a separate gate from quality, oracle,
documentation, release, and supported-platform CI checks. Record each command
as passing only after it has run against the exact tree under review.

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
