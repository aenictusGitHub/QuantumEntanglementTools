# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is a provisional, independent Julia package for
quantum-information and entanglement calculations. Its long-term goal is an
idiomatic, type-generic, sparse-aware Julia successor to the documented public
functionality of QETLAB, with explicit migration helpers and genuinely optional
backend integrations.

> [!WARNING]
> This repository is pre-alpha. Upstream inventory, implementation, validation,
> and API design are in progress. The local evidence below is not a claim of
> complete QETLAB parity, supported-release validation, a performance advantage,
> registry availability, or API stability.

## Current status

- Full local package suite: 2,295/2,295 assertions pass on Julia 1.12.6 and
  Julia 1.10.11: 2,270 core assertions plus 25 assertions that execute the
  published tutorials. This is a development baseline, not supported-platform
  or MATLAB validation.
- Inventory: 76 of 163 QETLAB rows have been manually reviewed: 63 are marked
  implemented, 11 partial, and two deferred; 87 remain pending.
- Tier A subsystem kernel: 304/304 focused assertions pass, covering
  dimensions/basis indices, tensor products and sums, subsystem
  permutation/swap, partial trace/transpose, realignment,
  symmetric/antisymmetric projectors, reusable plans, sparse/generic paths, and
  MATLAB-compatible wrappers.
- Tier B operators/states/random slice: 941/941 focused local tests pass. It
  includes operator bases, deterministic named states, and six randomized
  constructors whose native and compatibility APIs require an explicit RNG.
  The scalar bipartite Werner family is implemented, its multipartite form is
  partial, and `RandomSuperoperator` plus `RandomPPTState` are deferred.
- Tier C channels/maps: 154/154 focused assertions pass for the representation,
  conversion, application, physicality, constructor, and reviewed
  compatibility slice.
- Tier D measures/criteria and native entanglement pipeline: 162/162 and 68/68
  focused assertions pass, respectively. Necessary tests retain structured
  inconclusive outcomes rather than being reported as separability.
- Tier E: coherence passes 52/52; product analysis passes 175/175 native plus
  52/52 compatibility assertions; matrix analysis passes 127/127 native plus
  32/32 compatibility assertions; and matrix predicates pass 166/166 native
  plus 37/37 compatibility assertions. The `IsPSD` mapping remains partial
  because the pinned CVX symbolic branch is not implemented.
- Differential evidence: source-free fixtures generated with Octave 11.3.0
  from the pinned QETLAB revision provide function-specific supplemental
  checks. The matrix-analysis artifact passes 59/59 assertions across 22
  fixtures. No MATLAB-family oracle has been run for the matrix predicates, and
  Octave evidence is not general MATLAB equivalence.
- Quality and performance smoke: Aqua passes 11/11 and 25 representative JET
  probes pass; those JET probes do not cover the matrix predicates. All 42
  non-recording quick benchmark cases complete, without a regression threshold
  or comparative performance claim.
- QETLAB parity: not claimed; inventory/provenance review and MATLAB
  differential validation remain incomplete. Consult
  [`docs/PORTING_STATUS.md`](docs/PORTING_STATUS.md).
- Optional backends: the exact EntanglementDetection.jl 0.2.2 weak-dependency
  adapter passes 125/125 focused assertions on Julia 1.12.6. Searches run in a
  bounded child process so the audited backend's RNG/stdout/BLAS side effects
  do not escape into the caller. Backend conclusions remain uncertified
  candidate evidence and the package report stays `unknown`. The optional
  environment currently resolves on Julia 1.11 or later; its configured remote
  platform matrix has not run.

## Installation

The package is not registered. For local development, clone the repository and
use Julia's package manager from the checkout:

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

Randomized APIs never choose an implicit process-global stream:

```julia
using Random

rng = Xoshiro(0x514554)
ρ = random_density_matrix(rng, 4; rank = 2)
```

The `MATLABCompat` randomized spellings also require a leading RNG. See
[states, operators, and random objects](docs/src/states_operators_random.md)
for examples and current limitations.

Three deterministic tutorials run as ordinary scripts and as part of
`Pkg.test()`:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

See [Executable tutorials](docs/src/tutorials.md) for the standalone subsystem,
channel, and entanglement-certificate workflows.

## Design commitments

The core API will accept ordinary Julia vectors and matrices, preserve useful
numeric types and sparsity, make subsystem and normalization conventions
explicit, and report inconclusive numerical or solver outcomes honestly.
MATLAB-style compatibility names will be isolated from the Julia-native API.
Heavy solvers and third-party detection packages will remain optional.

## Documentation and development

Start with:

- [Getting started](docs/src/getting_started.md)
- [Mathematical conventions](docs/src/conventions.md)
- [States, operators, and random objects](docs/src/states_operators_random.md)
- [Product structure and separable-ball certificates](docs/src/product_analysis.md)
- [Matrix analysis](docs/src/matrix_analysis.md)
- [Matrix predicates](docs/src/matrix_predicates.md)
- [Executable tutorials](docs/src/tutorials.md)
- [EntanglementDetection extension](docs/src/entanglement_detection_extension.md)
- [Migration ledger](docs/src/migration_from_qetlab.md)
- [Contributing](CONTRIBUTING.md)
- [Legal and provenance status](docs/LEGAL.md)

Build the local documentation with:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Citation

The software is not yet associated with an archival release or DOI. Provisional
citation metadata is available in [`CITATION.cff`](CITATION.cff) and
[`CITATION.bib`](CITATION.bib); cite the exact commit used until a release is
published.

## License and non-affiliation

New project code and independently written documentation are licensed under the
BSD 3-Clause License; see [`LICENSE`](LICENSE). Source-derived portions may carry
additional compatible terms recorded in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md),
[`NOTICE`](NOTICE), and `PROVENANCE.toml`.

This project is independently maintained. It is not affiliated with, endorsed by,
or officially supported by QETLAB, its maintainers, or their contributors.
