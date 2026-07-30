# Symmetric-extension hierarchies

`symmetric_extension` and `symmetric_inner_extension` expose the outer and
inner semidefinite hierarchies as status-rich Julia APIs. They do not turn a
solver failure, iteration limit, or tolerance boundary into a Boolean answer.

The optional optimization layer is explicit:

```julia
using Hypatia
using QuantumEntanglementTools

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
    atol=5e-7,
    rtol=5e-7,
)
```

Hypatia is installed only in the dedicated optimization environment. It is not
a core dependency and there is no global default solver.

## Outer extensions

For an operator $X_{AB}$ and hierarchy order $k$, an outer extension is a
positive operator on

```math
\mathcal{H}_A \otimes \mathcal{H}_{B_1} \otimes \cdots
\otimes \mathcal{H}_{B_k}
```

whose $AB_1$ marginal is $X_{AB}$ and which is invariant under permutations of
the $B$ copies. The implementation orders tensor factors as
`A, B₁, ..., Bₖ`, with subsystem one the slowest-varying index. Adjacent-swap
equalities impose permutation invariance directly. A bosonic model instead
uses an orthonormal occupation-number basis and returns the lifted ambient
extension.

```julia
rho = Matrix{Float64}(I, 4, 4) / 4

result = symmetric_extension(
    rho;
    dims=(2, 2),
    order=2,
    ppt=true,
    bosonic=true,
    backend=backend,
    prefer_analytic=false,
    allow_densify=true,
)

@assert result.verdict === true
@assert result.optimization_result.status == OptimizationOptimal
@assert result.residuals.valid
@assert result.residuals.marginal_residual <= result.tolerance
@assert result.residuals.permutation_residual <= result.tolerance
@assert result.residuals.bosonic_support_residual <= result.tolerance
@assert result.residuals.ppt_violation <= result.tolerance
```

When `ppt=true`, the model includes the representative partial transposes over
$B_1$, then $B_1B_2$, and so on. Symmetry makes these representatives
sufficient for the corresponding copy subsets. Every returned primal
candidate is independently checked for:

- Hermiticity and positive semidefiniteness;
- trace and $AB_1$ marginal;
- all adjacent-copy symmetries;
- bosonic support when requested;
- every modeled partial-transpose constraint.

The complete package-owned `OptimizationResult` is retained, including raw
termination/primal/dual status symbols, objective data, residuals, optimizer
identity and version, options, work limits, and reconstructed primal/dual
values.

### Solver-free theorem branches

With `NoOptimizationBackend()` and `prefer_analytic=true`, the following
branches avoid an SDP:

- `order=1`: a PSD input is its own extension; with `ppt=true`, the input
  itself is returned only after the requested partial transpose is also
  checked, in any bipartite dimensions;
- two qubits, `order=2`, without PPT: the exact analytic criterion

```math
\mathrm{Tr}(\rho_B^2)
\geq
\mathrm{Tr}(\rho_{AB}^2)-4\sqrt{\det(\rho_{AB})}
```

  is used;
- total bipartite dimension at most six with `ppt=true`: the
  Peres--Horodecki theorem decides existence at every order;
- in every bipartite dimension, a robust negative partial transpose excludes
  any requested PPT extension and yields a normalized decomposable witness.

Floating values inside the declared tolerance band return
`SymmetricExtensionNumericalBoundary`, with `verdict === nothing`. A theorem
decision need not contain an explicit extension matrix. Supply a backend and
set `prefer_analytic=false` when a reconstructed primal extension is required.

For a robust NPT state, the solver-free branch constructs a decomposable
witness. For a solver-proved outer infeasibility, the marginal-equality dual is
converted to a Hermitian separator, normalized to

```math
\mathrm{Tr}(W X)=-1,
```

and accepted only if the complete dual cone/stationarity residual and
normalization residual pass independently.

## Inner extensions

`symmetric_inner_extension` implements the
Navascués--Owari--Plenio inner approximation to the separable cone. Its
positive variable is compressed to the symmetric subspace. For the non-PPT
hierarchy it imposes

```math
X_{AB}
=
\frac{k}{k+d_B}\sigma_{AB}
+
\frac{1}{k+d_B}\sigma_A \otimes I_B.
```

The PPT hierarchy uses the published Jacobi-root mixing coefficient. A private
descending-coefficient Jacobi recurrence, with independent polynomial and root
tests, supersedes QETLAB's internal `jacobi_poly` helper. No helper name is
exported.

```julia
inner = symmetric_inner_extension(
    rho;
    dims=(2, 2),
    order=2,
    ppt=true,
    backend=backend,
    allow_densify=true,
)

@assert inner.verdict === true
@assert inner.residuals.valid
@assert inner.problem.mixing_parameter ≈ 1 - inv(sqrt(3))
```

!!! warning "A negative inner result is not an entanglement certificate"
    A checked negative dual object separates the input from the selected inner
    approximation only. The inner cone is a subset of the separable cone, so
    this object need not be nonnegative on every separable state. Accordingly,
    `result.witness.entanglement_witness` is always `false` for the inner
    hierarchy. This preserves the warning in the pinned QETLAB documentation.

## Problem builders and limits

`symmetric_extension_problem` and `symmetric_inner_extension_problem` return
package-owned models containing `SemidefiniteProgram` values and affine
Hermitian views. They never contain JuMP expressions or solver-owned objects.
This permits model inspection and backend-independent testing.

The high-level APIs check `max_order`, `max_dense_entries`, and
`OptimizationLimits` before model allocation. Limit exhaustion is a structured
`SymmetricExtensionResourceLimit` with no verdict. Sparse input requires
`allow_densify=true` because input spectral validation and the portable complex
PSD real-block embedding are dense operations.

Malformed dimensions, nonfinite entries, non-Hermitian matrices, and robustly
non-PSD inputs are rejected. A negative input eigenvalue that remains inside
the tolerance boundary produces an inconclusive high-level result; low-level
problem builders reject it. The implementation never normalizes, symmetrizes,
clips, or projects the supplied operator.

## References and upstream scope

The formulations are independently checked against:

- [Doherty, Parrilo, and Spedalieri, *A complete family of separability criteria*](https://arxiv.org/abs/quant-ph/0308032);
- [Navascués, Owari, and Plenio, *A complete criterion for separability detection*](https://arxiv.org/abs/0906.2735);
- [Chen et al., *Symmetric Extension of Two-Qubit States*](https://arxiv.org/abs/1310.3530).

The compatibility surface maps the pinned `K`, `DIM`, `PPT`, `BOS`, and `TOL`
arguments, while requiring an explicit backend whenever QETLAB would have
selected CVX. Solver statuses are not collapsed to `0` or `1`.
