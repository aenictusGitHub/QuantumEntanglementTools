# QuantumEntanglementTools.jl

[![Live documentation](https://img.shields.io/badge/docs-live-blue.svg)](https://aenictusgithub.github.io/QuantumEntanglementTools/)
[![Documentation workflow](https://github.com/aenictusGitHub/QuantumEntanglementTools/actions/workflows/docs.yml/badge.svg?branch=main)](https://github.com/aenictusGitHub/QuantumEntanglementTools/actions/workflows/docs.yml)

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

## Five-minute start

QuantumEntanglementTools requires Julia 1.10 or later. No version has been
released or registered, so install the current development branch directly
from GitHub. Access to the private repository and working Git credentials are
required:

```julia
using Pkg
Pkg.add(
    url="https://github.com/aenictusGitHub/QuantumEntanglementTools.git",
    rev="main",
)
using QuantumEntanglementTools
```

For package development, clone the repository and replace `Pkg.add(...)` with
`Pkg.develop(path="/path/to/QuantumEntanglementTools")`.

Here is a complete separability calculation:

```julia
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
psi = tensor_product(ket0, ket1)

report = analyze_entanglement(psi, (2, 2))
@assert conclusion(report) === :separable
@assert is_certified(report)
println(report)
# EntanglementReport(status=separable, certified=true, method=pure_schmidt, attempts=1)
```

The conclusion is certified because the exact trailing Schmidt coefficients
vanish. A result with `status === :unknown` means that the requested methods
found neither a separability certificate nor an entanglement certificate; it
is not a negative answer.

## Choose an API

Start from the question you want to answer:

| Task | Recommended entry point | What to inspect |
|---|---|---|
| Reduce or rearrange subsystems | `partial_trace`, `partial_transpose`, `permute_subsystems` | Returned array and documented subsystem order |
| Work with symmetric multiqubit or multiqudit states | `generalized_dicke_state`, `symmetric_product_coordinates`, `symmetric_reduced_state` | Occupation ordering, representation, and resource guards |
| Classify a bipartite pure state | `analyze_entanglement(psi, dims)` | `EntanglementReport.status` and `certificate_kind` |
| Apply one necessary entanglement criterion | `ppt_criterion`, `realignment_criterion`, `reduction_criterion` | `CriterionResult.status`, witness, and tolerance |
| Run the dependency-free criterion pipeline | `analyze_entanglement(rho, dims)` | Ordered `attempts`; a pass can still lead to `:unknown` |
| Seek a separability or entanglement certificate | `is_separable(rho, dims; strategies=...)` | `certified`, `certificate_kind`, and `attempts` |
| Apply the sufficient separable-ball test | `in_separable_ball(rho, dims)` | `:separable_certified`, `:outside_ball`, or boundary status |
| Use an optional detection package | `detect_entanglement(rho, dims, method)` | Backend metadata and whether evidence is certified |
| Interpret any common structured result | `conclusion`, `is_conclusive`, `is_certified`, `explain` | Conservative conclusion and human-readable reason |

For one criterion, `CriterionEntanglementDetected` is a positive entanglement
detection; `CriterionSatisfied` only means that the state passed that
criterion; and `CriterionUnknown` records a tolerance or numerical boundary.

The [separability examples](docs/src/separability_examples.md) build explicit
mixed states and explain each conclusion. The
[interactive code generator](https://aenictusgithub.github.io/QuantumEntanglementTools/code_generator/)
creates complete, downloadable Julia scripts in the browser without executing
Julia or uploading parameters.

## A subsystem example

```julia
using QuantumEntanglementTools
using LinearAlgebra

psi = bell_state()
rho_A = partial_trace(psi, (2, 2); trace_out=(2,))

@assert isapprox(rho_A, Matrix{ComplexF64}(I, 2, 2) / 2)
```

Pure-vector reduction returns an operator matrix. Subsystem labels are
one-based, the first tensor factor is most significant, and tracing every
subsystem returns a `1 × 1` matrix. See
[Mathematical conventions](docs/src/conventions.md).

Randomized APIs never choose an implicit process-global stream:

```julia
using Random

rng = Xoshiro(0x514554)
ρ = random_density_matrix(rng, 4; rank = 2)
```

The `MATLABCompat` randomized spellings also require a leading RNG. See
[states, operators, and random objects](docs/src/states_operators_random.md)
for examples and current limitations.

Seven deterministic repository tutorials run as ordinary scripts and as part of
`Pkg.test()`:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

See [Executable tutorials](docs/src/tutorials.md) for the standalone subsystem,
channel, separability, symmetric-witness, entanglement-certificate,
Schmidt-decomposition, and Tiles-UPB bound-entanglement workflows.

## Current status

<!-- qetlab-current-claims: begin -->

- The inventory is pinned to QETLAB revision
  `d8589610f00cff106537268dee2e2a1153f3a601`. All 163 inventoried QETLAB files
  are source-reviewed. The strict static ledger reports 127/127 public rows
  are verified with the required final status, 36/36 internal helpers have
  terminal dispositions, the completion queue contains 0 public rows, 0
  required internal helpers remain, and there are 0 static completion failures.
- The package exports 477 public bindings (347 native/module and 130
  `MATLABCompat`), each with a matching provenance entry.
- The 9,759-assertion full package suite passed 9,759/9,759: the core accounts
  for 9,675 assertions, and the seven executable tutorials account for 84/84
  assertions (48 existing plus 36 for the two new workflows), on Julia
  1.12.6 and the installed Julia 1.10.11. The full optional JuMP suite passed
  847/847 on both Julia lines.
- The exact EntanglementDetection.jl 0.2.2 integration remains optional and
  child-process isolated. The EntanglementDetection.jl extension passed 141/141
  focused assertions on the current compatible Julia; heuristic output remains
  uncertified candidate evidence.
- All 120 declared quick benchmark cases completed without failure on the
  current uncommitted performance-and-stability worktree based on `1d611e4`.
  This is
  local smoke evidence only; targeted paired observations remain diagnostics,
  not a stable comparative-performance baseline.
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

## Design commitments

The core API accepts ordinary Julia vectors and matrices, preserves useful
numeric types and sparsity, makes subsystem and normalization conventions
explicit, and reports inconclusive numerical or solver outcomes honestly.
MATLAB-style compatibility names remain isolated from the Julia-native API.
Heavy solvers and third-party detection packages will remain optional.

## Documentation and development

Start with:

- [Live rendered documentation](https://aenictusgithub.github.io/QuantumEntanglementTools/)
- [Five-minute quick start](docs/src/getting_started.md)
- [Interactive entanglement code generator](https://aenictusgithub.github.io/QuantumEntanglementTools/code_generator/)
- [Separability by example](docs/src/separability_examples.md)
- [Symmetric SAPPT states and witnesses](docs/src/paper_symmetric_separability.md)
- [Symmetric multiqubit and multiqudit states](docs/src/symmetric_states.md)
- [Mathematical conventions](docs/src/conventions.md)
- [API reference](docs/src/api/index.md)
- [Executable tutorials](docs/src/tutorials.md)
- [EntanglementDetection extension](docs/src/entanglement_detection_extension.md)
- [Migration ledger](docs/src/migration_from_qetlab.md)
- [Contributing](CONTRIBUTING.md)
- [Support and correctness reports](SUPPORT.md)
- [Security policy](SECURITY.md)
- [Release checklist](docs/RELEASE_CHECKLIST.md)
- [Legal and provenance status](docs/LEGAL.md)

Build the local documentation with:

```sh
julia --startup-file=no scripts/build_docs.jl
```

The helper resolves the docs dependencies in a temporary environment for the
running Julia version and leaves the rendered site in `docs/build`. It exits
with an error if dependency resolution or the strict Documenter build fails.

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
