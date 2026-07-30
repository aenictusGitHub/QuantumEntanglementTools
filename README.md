# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is an independent Julia package for
quantum-information and entanglement calculations. Development currently
targets an unreleased, experimental `0.1.0` milestone toward complete
behavioral coverage of the public API at a pinned QETLAB revision. It provides
a type-generic, sparse-aware API, explicit QETLAB migration helpers, and
optional backend integrations.

> [!WARNING]
> The unreleased `0.1.0` development milestone is experimental.
> Static implementation completion is not QETLAB/MATLAB behavioral parity. It
> does not establish supported-platform remote CI, API stability, comparative
> performance, release approval, human review, or registry availability. Check
> the [porting status](docs/PORTING_STATUS.md) and
> [validation report](docs/VALIDATION_REPORT.md) before relying on a migration
> mapping or numerical certificate.

## Current status

<!-- qetlab-current-claims: begin -->

- The inventory is pinned to QETLAB revision
  `d8589610f00cff106537268dee2e2a1153f3a601`. All 163 inventoried QETLAB files
  are source-reviewed. The strict static ledger reports 127/127 public rows
  are verified with the required final status, 36/36 internal helpers have
  terminal dispositions, the completion queue contains 0 public rows, 0
  required internal helpers remain, and there are 0 static completion failures.
- The package exports 458 public bindings, each with a matching provenance
  entry.
- The 8,133-assertion full package suite passed 8,133/8,133, including 48
  executable-tutorial assertions, on Julia 1.12.6 and the installed Julia
  1.10.0. The full optional JuMP suite passed 836/836 on both Julia lines.
- The exact EntanglementDetection.jl 0.2.2 integration remains optional and
  child-process isolated. The EntanglementDetection.jl extension passed 125/125
  focused assertions on the current compatible Julia; heuristic output remains
  uncertified candidate evidence.
- All 114 declared quick benchmark cases completed for a clean `f32dd233`
  baseline and the candidate; targeted paired observations remain local
  diagnostics, not a stable performance baseline.
- These are static-ledger and local-test results, not a claim of complete
  QETLAB parity or MATLAB parity, supported-platform remote CI, comparative
  performance, API stability, release approval, or non-delegable human review.
  No version has been tagged or published.

<!-- qetlab-current-claims: end -->

Test, oracle, platform, benchmark, and limitation details live in the
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

`is_separable(rho, dims; strategies=...)` is deliberately a structured,
certificate-first API rather than a bare Boolean predicate. It returns an
`EntanglementReport` with ordered evidence and retains inconclusive outcomes as
`:unknown`, which must not be confused with entanglement.
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
