# Product structure and separable-ball certificates

This Tier E slice covers operator Schmidt decompositions, tolerance-aware
product analysis, exact entanglement-of-formation formulas in their known
closed-form domains, and the Gurvits--Barnum separable ball. The APIs are
dependency-free and deliberately distinguish a numerical classification from
an exact algebraic statement or entanglement certificate.

## Operator Schmidt decomposition

For a bipartite operator

```math
X : \mathbb{C}^{c_A}\otimes\mathbb{C}^{c_B}
    \longrightarrow
    \mathbb{C}^{r_A}\otimes\mathbb{C}^{r_B},
```

`operator_schmidt_decomposition(X, (rA, rB), (cA, cB))` returns coefficients
and Frobenius-orthonormal local factors satisfying

```math
X = \sum_k s_k A_k \otimes B_k
```

up to the numerical SVD reconstruction error. The one-layout form
`operator_schmidt_decomposition(X, dims)` is the square-operator convenience
method. All thin-decomposition coefficients are returned, including numerical
zeros. `operator_schmidt_coefficients` avoids constructing the local factors,
and `operator_schmidt_rank` returns a tolerance-defined numerical rank rather
than an exact symbolic rank. The typed result records `factor_convention` as
`:general` or `:hermitian`. For the Hermitian convention it also reports the
imaginary Hermitian-coordinate residual and the derived tolerance that bounded
its removal.

```@example product-analysis
using QuantumEntanglementTools

A = [1.0 2.0; 0.0 1.0]
B = [0.0 1.0; 1.0 0.0]
X = tensor_product(A, B)
decomposition = operator_schmidt_decomposition(X, (2, 2))
reconstruction = tensor_sum(
    decomposition.left_factors,
    decomposition.right_factors;
    weights = decomposition.coefficients,
)
isapprox(reconstruction, X)
```

The default `factor_convention=:general` path performs a full dense complex SVD
after realignment. Singular-vector phases and bases inside degenerate singular
subspaces are not unique. Compare coefficients, reconstruction, factor Gram
matrices, or subspace projectors rather than comparing factors entrywise.

For an exactly Hermitian operator with locally square layouts,
`hermitian_factors=true` changes to real orthonormal Hermitian operator bases:

```@example product-analysis
using LinearAlgebra

H = ComplexF64[
    1 0 0 0 0 0
    0 2 im 0 0 0
    0 -im 3 0 0 0
    0 0 0 4 0 0
    0 0 0 0 5 2im
    0 0 0 0 -2im 6
]
hermitian_decomposition = operator_schmidt_decomposition(
    H,
    (2, 3);
    hermitian_factors = true,
)
(
    hermitian_decomposition.factor_convention,
    all(ishermitian, hermitian_decomposition.left_factors),
    all(ishermitian, hermitian_decomposition.right_factors),
    isapprox(
        tensor_sum(
            hermitian_decomposition.left_factors,
            hermitian_decomposition.right_factors;
            weights = hermitian_decomposition.coefficients,
        ),
        H,
    ),
)
```

This construction works for unequal local dimensions such as `(2, 3)`. It
corrects the pinned routine's accidental linear indexing of its `DIM` matrix,
which used the first local dimension twice. The real coordinate extraction is
allowed only after its imaginary residual is below
`8 * max(dA²,dB²) * eps(T) * max(1, coordinate_scale)`; both values are stored
in the result. A near-Hermitian input is not projected or symmetrized and
raises when the Hermitian convention is explicitly requested. Locally
rectangular factor spaces cannot contain Hermitian factors and therefore
reject that explicit request.

Both conventions use dense SVDs. Sparse inputs therefore require
`allow_densify=true`, and only `Float32`, `Float64`, `ComplexF32`, and
`ComplexF64` are accepted. The operator is never normalized, symmetrized,
clipped, or mutated.

## Product-vector and product-operator analysis

`is_product_vector` and `is_product_operator` return
[`ProductAnalysisResult`](@ref), not a `Bool`. They recursively compare the
full nonleading Schmidt-tail norm at every cut and then check the complete
product reconstruction. Residuals and thresholds are reported in the original
input's units.

The three statuses are:

- `:within_tolerance`: every cut and the final reconstruction are below their
  requested thresholds;
- `:outside_tolerance`: at least one residual is above threshold;
- `:boundary`: no residual is above threshold and at least one positive
  residual equals its threshold.

```@example product-analysis
product_state = tensor_product(
    ComplexF64[1, 2im],
    ComplexF64[3, -1],
    ComplexF64[2, im],
)
analysis = is_product_vector(product_state, (2, 2, 2))
(analysis.status, length(analysis.factors), analysis.reconstruction_residual)
```

The result is a tolerance-aware numerical classification, not a symbolic
factorization proof. Factors are returned even for an outside-tolerance input
so that the leading product approximation is inspectable. Zero vectors and
zero operators are rejected because their product factors are nonunique.
Multipartite operators may use separate row and column layouts, so rectangular
local factors are supported.

## Entangling-gate certificates

[`is_entangling_gate`](@ref) returns an
[`EntanglingGateResult`](@ref), never an unchecked Boolean. It accepts a shared
subsystem layout, a scalar first bipartite dimension, or a two-row matrix of
output and input local dimensions. The three statuses are:

- `:not_entangling`, with local factors and a subsystem permutation that
  reconstruct the gate;
- `:entangling`, with a normalized product input whose output has an
  outside-tolerance product analysis;
- `:unknown`, when unitarity or product structure lies on a numerical
  boundary or the bounded witness grid is inconclusive.

```@example product-analysis
using LinearAlgebra

controlled_z = Matrix(Diagonal(ComplexF64[1, 1, 1, -1]))
gate_result = is_entangling_gate(controlled_z, (2, 2))
(
    status = gate_result.status,
    witness_support = count(!iszero, gate_result.product_witness),
    output_status = gate_result.certificate_analysis.status,
)
```

Before witness search, every output-subsystem permutation is tested for a
local tensor factorization. Factorial permutation work is guarded. The search
then checks the pinned one- and two-support product inputs and, if necessary,
a corrected finite phase grid of local basis vectors and normalized
`eᵢ + z*eⱼ` vectors for `z in (1, -1, im, -im)`. Candidate count and dense
matrix-vector work are preflighted separately. Sparse gates require explicit
densification, and no gate is normalized or projected onto a unitary.

The pinned routine's witness loop is not reliable for all entangling gates.
In particular, its controlled-Z call returns the true entangling flag but
falls through with an unnormalized two-support candidate whose output is
still product. The Julia implementation records that discrepancy and returns
a validated four-support product witness instead. A six-fixture source-free
artifact generated by Octave 11.3.0 from the pinned revision passes 13
comparison and defect assertions; authoritative MATLAB execution was not run.

## Minimum UPB size

[`minimum_upb_size`](@ref) is an exact-or-unknown lookup for the orthogonal
unextendible-product-basis theorem table reviewed at the pinned QETLAB
revision. Its [`MinimumUPBSizeResult`](@ref) contains the exact minimum when
known, the universal counting lower bound, canonical sorted dimensions, and a
primary-source citation. An uncovered multipartite family returns
`:unknown` with `size === nothing`; it is not assigned a guessed value.

```@example product-analysis
known_upb_minimum = minimum_upb_size((4, 6))
unknown_upb_minimum = minimum_upb_size((2, 3, 4))
(
    known = (
        known_upb_minimum.status,
        known_upb_minimum.size,
        known_upb_minimum.reference_key,
    ),
    unresolved = (
        unknown_upb_minimum.status,
        unknown_upb_minimum.size,
        unknown_upb_minimum.lower_bound,
    ),
)
```

The table is grounded in the primary results of
[Alon and Lovász](https://doi.org/10.1006/jcta.2000.3122),
[DiVincenzo et al.](https://arxiv.org/abs/quant-ph/9908070),
[Feng](https://doi.org/10.1016/j.dam.2005.10.011),
[Chen and Johnston](https://arxiv.org/abs/1301.1406), and
[Johnston](https://arxiv.org/abs/1302.1604). The API deliberately describes
unknown as “not determined by the reviewed table,” rather than asserting that
no later or more specialized result exists.

`MATLABCompat.MinUPBSize` preserves the known-case integer return and optional
reference output. It raises a `DomainError` carrying the native structured
result for an unresolved case, matching the pinned error route without losing
the lower bound. Eleven exact table branches and one unknown branch pass a
source-free Octave 11.3.0 fixture generated from the pinned QETLAB revision;
authoritative MATLAB execution was not run.

To analyze a supplied family rather than a dimension table, use
[`is_upb`](@ref). It checks product structure, mutual orthogonality,
incompleteness, and unextendibility, and returns either an exhaustion
certificate, an explicit product extension, or a structured `:unknown` when a
numerical or resource boundary prevents a conclusion. The exact, numerical,
complex-witness, sparse, and combinatorial contracts are documented in
[UPB certificates and extension witnesses](is_upb.md).

## Entanglement of formation

`entanglement_of_formation` implements only domains with a verified exact
closed form:

- normalized bipartite pure vectors in arbitrary local dimensions, using the
  entropy of squared Schmidt coefficients;
- normalized rank-one density matrices in arbitrary local dimensions, using
  the same pure-state formula after recovering the projected vector;
- normalized two-qubit density matrices, using Wootters concurrence and binary
  entropy.

The logarithm base defaults to two and may be changed explicitly. An exactly
idempotent rank-one projector is recognized without a numerical projection.
For a floating matrix that lies only within the numerical rank-one boundary,
set `rank_boundary_policy=:project` to explicitly discard the bounded
spectral tail and use the normalized leading eigenvector. Sparse spectral work
requires explicit densification.

All matrix boundary policies default to `:reject`.
`psd_boundary_policy=:project`, `rank_boundary_policy=:project`, and
`range_boundary_policy=:project` permit only the corresponding bounded
correction already proved to lie inside the requested numerical tolerance.
Mixed density matrices outside the two-qubit domain remain unsupported.

```@example product-analysis
bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
rectangular_projector = zeros(6, 6)
rectangular_projector[[1, 5], [1, 5]] .= 0.5
(
    entanglement_of_formation(bell, (2, 2)),
    entanglement_of_formation(bell * bell', (2, 2)),
    entanglement_of_formation(rectangular_projector, (2, 3)),
)
```

## Gurvits--Barnum separable ball

For a normalized bipartite state of total dimension `D > 1`,
`in_separable_ball` tests the sufficient condition

```math
\mathrm{tr}(\rho^2) \leq \frac{1}{D-1}.
```

The result is a [`SeparableBallResult`](@ref):

- `:separable_certified` is a sufficient separability certificate;
- `:outside_ball` means only that this test did not certify the state;
- `:unknown` records a numerical boundary in input validation or the ball
  inequality.

An outside-ball result is never evidence of entanglement. The native API
requires an already normalized density matrix or normalized eigenvalue vector
and never divides by the trace. `Diagonal` matrices and supplied eigenvalues
use structure-aware linear-time paths; a general matrix requires dense
spectral validation.

```@example product-analysis
using LinearAlgebra

maximally_mixed = Diagonal(fill(0.25, 4))
certificate = in_separable_ball(maximally_mixed, (2, 2))
(certificate.status, certificate.purity, certificate.boundary)
```

## QETLAB compatibility and evidence

The `MATLABCompat` namespace exposes
`OperatorSchmidtDecomposition`, `OperatorSchmidtRank`, `IsProductVector`,
`IsProductOperator`, `IsEntanglingGate`, `IsUPB`, `MinUPBSize`,
`EntFormation`, and `InSeparableBall`. The wrappers retain reviewed dimension and selection
defaults where doing so is safe, but return structured Julia results when a
Boolean would erase a numerical boundary or a failed sufficient certificate.

Notable intentional differences are:

- `OperatorSchmidtDecomposition` returns a three-field named tuple that also
  supports positional destructuring as `s, U, V = ...`. `K=0`, `K=-1`, and
  positive `K` select numerically nonzero, full-thin, and leading terms,
  respectively.
- For exactly Hermitian input and locally square `DIM`, the wrapper
  automatically requests Hermitian factors. It corrects the pinned
  unequal-local-dimension indexing failure and reapplies `K` after repair; the
  pinned branch could accidentally ignore `K`. Automatic mode never treats a
  merely near-Hermitian input as Hermitian. `hermitian_factors=true` can require
  the convention, while `false` disables it.
- General rectangular row/column layouts use the valid thin decomposition
  instead of reproducing the pinned reshape failures recorded by the committed
  Octave discrepancy fixtures. If a globally Hermitian matrix is described by
  locally rectangular factor spaces, automatic mode returns general factors
  because rectangular matrices cannot themselves be Hermitian.
- `IsProductVector` and `IsProductOperator` return
  `ProductAnalysisResult`; boundary cases are not collapsed to `true` or
  `false`.
- `IsEntanglingGate` returns the same certificate-bearing structured result as
  the native API. It corrects the pinned controlled-Z witness fall-through and
  retains numerical boundary cases as `:unknown`.
- `IsUPB` returns the reviewed Boolean for conclusive calls, optionally returns
  compact witness factors, and exposes the native result with
  `structured=true`. It raises on `:unknown` rather than turning an
  inconclusive bounded search into `false`.
- `MinUPBSize` returns the pinned integer for known theorem-table cases and
  raises with the native structured result for unresolved dimensions; the
  native API returns that case directly as `:unknown`.
- `EntFormation` returns the mathematical zero for zero-concurrence two-qubit
  states instead of reproducing the pinned implementation's `NaN`. It covers
  the pinned pure-vector, arbitrary-dimensional rank-one-density, and mixed
  two-qubit branches, while requiring explicit opt-in for numerical boundary
  projections.
- `InSeparableBall` retains QETLAB's positive-trace normalization convenience
  but validates finite positive semidefinite input and returns a structured
  sufficient-certificate result.

The committed Octave 11.3.0/QETLAB artifact contains 14 deterministic fixtures
and passes 68 native, compatibility, reconstruction, and discrepancy
assertions. Its SHA-256 is
`ab6414c1a684141db74782616d4c18e79c8e6039aad53a695d8c723eed598d85`.
Product Booleans are compared only away from tolerance boundaries. This is
function-specific supplemental evidence; MATLAB was not run, and factor
vectors are validated by reconstruction rather than entrywise phase
comparison.
