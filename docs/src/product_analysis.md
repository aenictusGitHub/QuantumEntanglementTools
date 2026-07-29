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
than an exact symbolic rank.

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

The implementation performs a full dense SVD after realignment. Sparse inputs
therefore require `allow_densify=true`, and only BLAS floating element types
are accepted. It does not normalize the operator or repair singular-vector
phases. In particular, Hermitian input does not imply Hermitian local factors.

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

## Entanglement of formation

`entanglement_of_formation` implements only domains with a verified exact
closed form:

- normalized bipartite pure vectors in arbitrary local dimensions, using the
  entropy of squared Schmidt coefficients;
- normalized two-qubit density matrices, using Wootters concurrence and binary
  entropy.

The logarithm base defaults to two and may be changed explicitly. A rank-one
density matrix in dimensions other than `(2, 2)` is not silently converted to
a pure vector: pass the state vector when using the arbitrary-dimensional pure
formula. Sparse spectral work requires explicit densification.

For two-qubit density matrices, `psd_boundary_policy=:project` and
`range_boundary_policy=:project` permit only bounded corrections already
proved to lie inside the requested numerical tolerances. Use `:reject` to
refuse those boundary corrections. Inputs are never normalized.

```@example product-analysis
bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
(
    entanglement_of_formation(bell, (2, 2)),
    entanglement_of_formation(bell * bell', (2, 2)),
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
`IsProductOperator`, `EntFormation`, and `InSeparableBall`. The wrappers retain
reviewed dimension and selection defaults where doing so is safe, but return
structured Julia results when a Boolean would erase a numerical boundary or a
failed sufficient certificate.

Notable intentional differences are:

- `OperatorSchmidtDecomposition` returns a named result and does not reproduce
  QETLAB's Hermitian-factor repair; the pinned implementation's failures for
  rectangular and unequal-local-dimension Hermitian inputs are recorded, not
  emulated.
- `IsProductVector` and `IsProductOperator` return
  `ProductAnalysisResult`; boundary cases are not collapsed to `true` or
  `false`.
- `EntFormation` returns the mathematical zero for zero-concurrence two-qubit
  states instead of reproducing the pinned implementation's `NaN`, and it does
  not auto-convert higher-dimensional rank-one density matrices.
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
