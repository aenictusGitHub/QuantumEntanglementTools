# Five-minute quick start

QuantumEntanglementTools works with ordinary Julia vectors and matrices while
keeping subsystem dimensions and certificate meaning explicit. It requires
Julia 1.10 or later.

## 1. Install the development version

No version has been released or registered. If your GitHub account has access
to the private repository, install its current development branch with:

```julia
using Pkg
Pkg.add(
    url="https://github.com/aenictusGitHub/QuantumEntanglementTools.git",
    rev="main",
)
using QuantumEntanglementTools
```

For package development, clone the repository and use
`Pkg.develop(path="/path/to/QuantumEntanglementTools")` instead.

## 2. Obtain a separability certificate

Constructing the tensor factors explicitly makes the subsystem order visible:

```@example quick-start-product
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
psi = tensor_product(ket0, ket1)

report = analyze_entanglement(psi, (2, 2))
@assert conclusion(report) === :separable
@assert is_certified(report)
report
```

The report is `:separable` and certified by an exact pure-product
decomposition. Its compact display also identifies the method and number of
attempts.

## 3. Interpret the result

For `EntanglementReport`, always read `status` together with `certified`:

| Status | Meaning |
|---|---|
| `:separable` with `certified == true` | The report contains a valid separability certificate for the stated inputs and tolerances. |
| `:entangled` with `certified == true` | The report contains a valid entanglement certificate or witness. |
| `:unknown` | The selected methods established neither conclusion. It is not evidence for either side. |

Use `conclusion(report)`, `is_conclusive(report)`, `is_certified(report)`, and
`explain(report)` for a consistent interface across the most common
entanglement, criterion, separable-ball, and optimization results. The original
status and all detailed fields remain available. Inspect
`report.certificate_kind` and `report.attempts` when the route to the
conclusion matters.

Single criteria use a related vocabulary:
`CriterionEntanglementDetected` is a positive detection,
`CriterionSatisfied` means only that the state passed that necessary test, and
`CriterionUnknown` records a tolerance or numerical boundary.

## 4. Choose the right entry point

| Question | Entry point | Result to inspect |
|---|---|---|
| How do I reduce, transpose, or reorder subsystems? | `partial_trace`, `partial_transpose`, `permute_subsystems` | Returned array and the subsystem convention |
| Is this bipartite pure state a product state? | `analyze_entanglement(psi, dims)` | `EntanglementReport` |
| Does one criterion detect entanglement? | `ppt_criterion`, `realignment_criterion`, `reduction_criterion` | `CriterionResult` |
| What do the dependency-free criteria conclude? | `analyze_entanglement(rho, dims)` | Ordered `attempts` in `EntanglementReport` |
| Can a configured strategy certify either conclusion? | `is_separable(rho, dims; strategies=...)` | Certificate, evidence, and ordered attempts |
| Is the state inside the sufficient separable ball? | `in_separable_ball(rho, dims)` | `SeparableBallResult`; `:outside_ball` is not entanglement |
| Should an optional detector analyze it? | `detect_entanglement(rho, dims, method)` | Backend metadata and certification flag |

## Diagnose an input without changing it

When a density matrix is rejected, inspect its dimensions and numerical
residuals without normalizing or repairing it:

```@example quick-start-validation
using QuantumEntanglementTools

rho = [0.5 0.1; 0.1 0.5]
diagnostics = validate_density_matrix(rho, (2,))

@assert diagnostics.valid
(
    diagnostics.status,
    diagnostics.trace_residual,
    diagnostics.hermiticity_residual,
    diagnostics.minimum_eigenvalue,
)
```

For a sparse non-diagonal matrix, the spectral check remains incomplete unless
`allow_densify=true` is supplied explicitly. The report explains this rather
than converting storage silently.

!!! note "Common surprises"
    - `prod(dims)` must match the vector length or matrix dimension.
    - Subsystem indices are one-based and tensor factors are ordered from left
      to right, with the first factor most significant.
    - `unknown`, `CriterionSatisfied`, and `:outside_ball` are inconclusive;
      none means “not entangled.”
    - Randomized methods require an explicit `rng::AbstractRNG`.
    - Sparse spectral work and SDP translation require explicit densification
      permission.
    - JuMP, concrete optimizers, and EntanglementDetection.jl are optional and
      are never loaded or selected implicitly.

For explicit mixed separable states, higher-dimensional inconclusive outcomes,
and the separable-ball test, continue with
[Separability by example](separability_examples.md). For a five-qubit family
with a finite decomposition and constructive witnesses, see
[Symmetric SAPPT states and witnesses](paper_symmetric_separability.md).

## 5. Try a subsystem calculation

```@example quick-start-subsystems
using QuantumEntanglementTools
using LinearAlgebra

psi = bell_state()
rho_A = partial_trace(psi, (2, 2); trace_out=(2,))
@assert isapprox(rho_A, Matrix{ComplexF64}(I, 2, 2) / 2)
rho_A
```

Pure-vector reduction returns an operator matrix. Subsystem labels are
one-based and the first tensor factor is most significant. The
[mathematical conventions](conventions.md) page illustrates all ordering
choices.

## Where to go next

- [Executable tutorials](tutorials.md) progress from subsystem operations and
  Schmidt decomposition to separability certificates and bound entanglement.
- The [interactive code generator](code_generator.md) assembles a complete
  local Julia script from a curated state family and analysis route.
- [Migration from QETLAB](migration_from_qetlab.md) maps familiar QETLAB names
  to the Julia-native and `MATLABCompat` APIs.
- [API overview](api/index.md) provides a task-based map of the full reference.

## Randomized operations require an RNG

No randomized package method draws from an implicit global stream:

```julia
using QuantumEntanglementTools
using Random
using LinearAlgebra

rng = Xoshiro(2026)
p = random_probabilities(rng, 6)
U = random_unitary(rng, 4)
effects = random_povm(rng, 3, 4)

@assert sum(p) ≈ 1
@assert U' * U ≈ I
@assert sum(effects) ≈ I
```

The compatibility namespace follows the same safety rule:

```julia
rng = Xoshiro(2026)
ρ = MATLABCompat.RandomDensityMatrix(rng, 4, 0, 2, "hs")
```

Read [States, operators, and random objects](states_operators_random.md) before
using the Tier B slice; it records physical parameter ranges, sparse behavior,
oracle evidence, resource limits, and compatibility differences.

The core supports Julia 1.10. The optional exact-version
[EntanglementDetection.jl extension](entanglement_detection_extension.md)
currently resolves on Julia 1.11 or later because of its backend dependency
graph; it is not installed by the core workflow above.

## Development status and validation

The experimental `0.1.0` milestone is unreleased and unregistered. At pinned
QETLAB revision `d8589610f00cff106537268dee2e2a1153f3a601`, the strict static
completion checker passes with 127/127 public rows carrying final `verified`
status, 36/36 internal helpers assigned terminal dispositions, no queued rows,
and 477 exported bindings with matching provenance entries. The API may still
change before publication or in later `0.x` releases.

The direct local full corpus passes 9,759/9,759 assertions: the core accounts
for 9,675 assertions, and the seven standalone executable tutorials
pass 84/84 assertions (48 existing plus 36 for the two new workflows), on Julia
1.12.6 and the installed Julia 1.10.11. This is local implementation and
validation evidence, not MATLAB/QETLAB parity, remote supported-platform CI,
comparative performance, API stability, release approval, or human review.

Maintainers can reproduce the core and tutorial gates from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --startup-file=no --project=. tutorials/runtests.jl
```

## Build these docs

```sh
julia --startup-file=no scripts/build_docs.jl
```

The helper copies `docs/Project.toml` into a temporary environment, resolves it
for the running Julia version, and leaves the rendered site in `docs/build`.
It does not reuse or write the ignored `docs/Manifest.toml`, and exits with an
error when dependency resolution or the strict build fails.

## Before relying on a function

Check all of the following:

1. its inventory row records the completed mapping and its evidence;
2. `PROVENANCE.toml` records its specification and implementation origin;
3. the API reference describes shapes, conventions, tolerances, errors, and
   certification meaning;
4. applicable analytic, property, invalid-input, sparse/generic, and independent
   or differential tests pass;
5. any optional backend reports its version and solver status.

For separability and nonlocal optimization in particular, inspect the
structured status rather than treating a numerical value or necessary-test
pass as a theorem; see
[Separability and local discrimination](separability_optimization.md) and
[Bell inequalities and nonlocal games](nonlocal_optimization.md).
