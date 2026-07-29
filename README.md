# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is an independent Julia package for
quantum-information and entanglement calculations. Development currently
targets an unreleased, experimental `0.1.0` milestone with a scoped,
type-generic, sparse-aware API, explicit QETLAB migration helpers, and optional
backend integrations.

> [!WARNING]
> The unreleased `0.1.0` development milestone is experimental. It is not a
> claim of complete QETLAB parity, API stability, comparative performance, or
> registry availability. Check the [porting status](docs/PORTING_STATUS.md) and
> [validation report](docs/VALIDATION_REPORT.md) before relying on a migration
> mapping or numerical certificate.

## Current status

- The local package corpus contains 2,417 passing assertions on Julia 1.10.11
  and 1.12.6: 2,369 core assertions plus 48 assertions that execute the
  published tutorials.
- All 163 inventoried QETLAB files have source-reviewed dispositions. Among the
  127 public rows, 63 mappings are implemented, 15 are partial, 19 are
  deferred, and 30 are blocked with explicit reasons; none remains pending.
  The 36 private helpers are tracked separately.
- The implemented scope covers subsystem operations, states and random objects,
  channels, scalar measures and criteria, certificate-aware entanglement
  analysis, coherence, product analysis, and selected matrix analysis.
- The exact EntanglementDetection.jl 0.2.2 integration is optional and isolated
  in a child process. Its heuristic output remains uncertified candidate
  evidence and never becomes a package-owned separability certificate.
- Full test, oracle, platform, benchmark, and limitation details live in the
  [porting status](docs/PORTING_STATUS.md),
  [inventory source review](docs/INVENTORY_REVIEW.md),
  [validation report](docs/VALIDATION_REPORT.md), and
  [benchmark report](docs/BENCHMARK_REPORT.md).

## Installation

No version has been released or registered. Clone the repository and use
Julia's package manager from the checkout:

```julia
using Pkg
Pkg.develop(path="/path/to/QuantumEntanglementTools")
using QuantumEntanglementTools
```

## Small tested example

The Tier A suite checks this Bell-state reduction using the Tier B constructor:

```julia
using QuantumEntanglementTools
using LinearAlgebra

ψ = bell_state()
ρA = partial_trace(ψ, (2, 2); trace_out = (2,))

@assert isapprox(ρA, Matrix{ComplexF64}(I, 2, 2) / 2)
```

Pure-vector reduction returns an operator matrix. Subsystem labels are one-based,
the first tensor factor is most significant, and tracing every subsystem returns
a `1 × 1` matrix. See [`docs/src/conventions.md`](docs/src/conventions.md).

### Separability quick start

Construct product states explicitly so the subsystem order remains visible:

```julia
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
ψ01 = tensor_product(ket0, ket1)

report = analyze_entanglement(ψ01, (2, 2))

@assert report.status === :separable
@assert report.certified
@assert report.certificate_kind === :pure_product_decomposition
```

There is intentionally no general Boolean `is_separable`: mixed-state
separability is hard, and `:unknown` must not be confused with entanglement.
The [separability examples](docs/src/separability_examples.md) build explicit
mixed states, compare the native pipeline with `in_separable_ball`, and explain
every status and certificate. The
[symmetric SAPPT example](docs/src/paper_symmetric_separability.md) adds an
explicit five-qubit separable decomposition and constructive entanglement
witnesses for the state family in Phys. Rev. A 111, 042418 (2025).
The browser-local
[entanglement example code generator](docs/src/code_generator.md) turns
curated state families, criteria, measures, and witness choices into complete
downloadable Julia scripts without executing or uploading the selected
parameters.

Randomized APIs never choose an implicit process-global stream:

```julia
using Random

rng = Xoshiro(0x514554)
ρ = random_density_matrix(rng, 4; rank = 2)
```

The `MATLABCompat` randomized spellings also require a leading RNG. See
[states, operators, and random objects](docs/src/states_operators_random.md)
for examples and current limitations.

Five deterministic repository tutorials run as ordinary scripts and as part of
`Pkg.test()`:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

See [Executable tutorials](docs/src/tutorials.md) for the standalone subsystem,
channel, separability, symmetric-witness, and entanglement-certificate
workflows.

## Design commitments

The core API accepts ordinary Julia vectors and matrices, preserves useful
numeric types and sparsity, makes subsystem and normalization conventions
explicit, and reports inconclusive numerical or solver outcomes honestly.
MATLAB-style compatibility names remain isolated from the Julia-native API.
Heavy solvers and third-party detection packages will remain optional.

## Documentation and development

Start with:

- [Getting started](docs/src/getting_started.md)
- [Entanglement example code generator](docs/src/code_generator.md)
- [Separability by example](docs/src/separability_examples.md)
- [Symmetric SAPPT states and witnesses](docs/src/paper_symmetric_separability.md)
- [Mathematical conventions](docs/src/conventions.md)
- [API reference](docs/src/api/index.md)
- [Executable tutorials](docs/src/tutorials.md)
- [EntanglementDetection extension](docs/src/entanglement_detection_extension.md)
- [Migration ledger](docs/src/migration_from_qetlab.md)
- [Contributing](CONTRIBUTING.md)
- [Release checklist](docs/RELEASE_CHECKLIST.md)
- [Legal and provenance status](docs/LEGAL.md)

Build the local documentation with:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Citation

Development citation metadata is available in [`CITATION.cff`](CITATION.cff)
and [`CITATION.bib`](CITATION.bib). No version or archival DOI has been
published; include the exact commit used.

## Development and review disclosure

This repository has received substantial assistance from OpenAI Codex in its
implementation, tests, documentation, and release preparation. The maintainer
is responsible for every distributed result. This disclosure does not claim
that the non-delegable human review required for General-registry submission has
been completed; registration must wait until that review is recorded in the
[release checklist](docs/RELEASE_CHECKLIST.md).

## License and non-affiliation

New project code and independently written documentation are licensed under the
BSD 3-Clause License; see [`LICENSE`](LICENSE). Source-derived portions may carry
additional compatible terms recorded in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md),
[`NOTICE`](NOTICE), and `PROVENANCE.toml`.

This project is independently maintained. It is not affiliated with, endorsed by,
or officially supported by QETLAB, its maintainers, or their contributors.
