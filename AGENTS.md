# AGENTS.md

This file is the operational entry point for maintainers and coding agents working
in this repository. Read it together with
[`docs/SESSION_HANDOFF.md`](docs/SESSION_HANDOFF.md) and
[`docs/PORTING_STATUS.md`](docs/PORTING_STATUS.md) before changing code.

## Project state

`QuantumEntanglementTools` is a pre-alpha Julia package intended to become an
idiomatic, independently maintained successor to selected QETLAB functionality.
The public API, upstream inventory, and validation baseline are not complete.
Never infer completeness from a file existing or a symbol being exported.

The package name is provisional. It is intentionally neutral and must not be
described as an official QETLAB project.

## Sources of truth

- `porting/qetlab_inventory.toml` is the completeness ledger once generated.
- `UpstreamManifest.toml` records upstream revisions and license evidence.
- `PROVENANCE.toml` records the origin, implementation kind, tests, and docs for
  each public function.
- `docs/src/conventions.md` records user-visible mathematical conventions.
- `docs/PORTING_STATUS.md` summarizes verified milestone status.
- `docs/SESSION_HANDOFF.md` records the current working state and blockers.

If these disagree, stop making broad claims, resolve the discrepancy from primary
evidence, and update all affected records.

## Non-negotiable engineering rules

1. Inspect `git status` before editing and preserve unrelated or pre-existing
   work. Do not reset or clean the worktree.
2. Do not export placeholders, fabricated results, or unchecked numerical
   conclusions.
3. A public function is not complete until its provenance, specification,
   implementation, tests, and documentation are all present.
4. Preserve meaningful numeric element types and sparse structure. Any
   densification must be explicit, documented, guarded, and tested.
5. Randomized APIs accept an explicit `rng::AbstractRNG`; never mutate the
   caller's global random stream.
6. Distinguish certificates from necessary tests and heuristics. An
   inconclusive computation is `unknown`, not a negative mathematical result.
7. Optional solvers and integrations belong in package extensions. Do not use
   private dependency APIs, type piracy, runtime method injection, or
   `Requires.jl`.
8. Do not claim parity or performance without recorded validation or benchmark
   evidence.
9. Keep source-derived licensing and attribution at file/function granularity.

## Julia conventions

- Minimum supported Julia version: 1.10. The optional
  EntanglementDetection.jl environment has an effective Julia 1.11 resolver
  floor because of its Ket 0.9 dependency.
- Prefer one top-level module and standard `AbstractVector`/`AbstractMatrix`
  inputs.
- Use lowercase `snake_case` for the Julia-native API. MATLAB-compatible names
  belong in a separate compatibility namespace.
- Use `adjoint` where conjugation is required and `A \ b` or factorizations
  instead of `inv(A) * b`.
- Validate subsystem dimensions and indices centrally. Never silently normalize,
  symmetrize, clip, or repair user input.
- Document subsystem order, vectorization, normalization, shapes, tolerances,
  failure behavior, and complexity for every important operation.

## Expected checks

Run the smallest relevant test while iterating, then the full applicable checks:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --startup-file=no --project=. tutorials/runtests.jl
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
julia --project=benchmark benchmark/benchmarks.jl --quick --no-save
julia --startup-file=no --project=test/extensions/entanglement_detection \
  test/extensions/entanglement_detection/runtests.jl
```

Instantiate dedicated optional-extension environments as documented in
`CONTRIBUTING.md`; do not add an optional backend to the core test target.

Before handing work off, update the status and handoff documents with commands
actually run, failures, uncommitted files, and external blockers. Do not leave
critical information only in chat.
