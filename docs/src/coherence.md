# Coherence

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

The dependency-free coherence slice implements four basis-dependent
quantities for normalized pure vectors and validated density matrices:

- `l1_coherence`;
- `relative_entropy_coherence`;
- `coherence_rank` for pure vectors; and
- `pure_k_coherence_robustness` for pure vectors.

Inputs are never normalized, repaired, or clipped. Pure-vector paths preserve
generic real or complex floating-point precision, including `BigFloat`, and do
not densify sparse vectors. Matrix entropy and positivity validation require a
BLAS floating element type because the core deliberately has no
arbitrary-precision eigensolver dependency.

## Definitions and example

For a density matrix $\rho$ in the selected computational basis,

```math
C_{l_1}(\rho) = \sum_{i \ne j} |\rho_{ij}|,
\qquad
C_{\mathrm{rel}}(\rho) = S(\mathrm{diag}(\rho)) - S(\rho).
```

The logarithm base for `relative_entropy_coherence` is required explicitly.
For a pure vector $\psi$, coherence rank is the number of basis
coefficients above the requested numerical threshold.

```jldoctest coherence
julia> ψ = ComplexF64[1, 1] / sqrt(2);

julia> isapprox(l1_coherence(ψ), 1)
true

julia> isapprox(relative_entropy_coherence(ψ; base = 2), 1)
true

julia> coherence_rank(ψ)
2

julia> pure_k_coherence_robustness(ψ, 2).value
0.0
```

A supplied `basis` for `coherence_rank` is a finite unitary matrix whose
columns are basis vectors. After the unitarity check, coordinates use the
adjoint basis action directly, never an explicit inverse or generic
factorization. A sparse basis transform may become dense and therefore
requires `allow_densify=true`.

The pure-state robustness operation first sorts coefficient magnitudes, so
coordinate permutations and complex basis phases do not change its value.
Its structured result records the largest branch selected by Theorem 1,
including the theorem's non-strict equality rule, together with adjacent gaps
that distinguish stable and tolerance-near boundaries. The complete formula,
branch contract, and precision behavior are documented in
[Pure-state robustness of k-coherence](pure_k_coherence_robustness.md).

## Validation and numerical policy

- Pure vectors must be finite and normalized within `atol` and `rtol`.
- Density matrices must be finite, square, normalized, Hermitian, and positive
  semidefinite within the documented tolerances.
- Sparse density matrices require `allow_densify=true` because validation uses
  a dense eigendecomposition.
- No entropy probability, eigenvalue, coherence value, or rank coefficient is
  clipped.
- Rank is a tolerance-defined numerical count, not an exact symbolic rank.

The pure `l1` path uses
$(\sum_i |\psi_i|)^2-\lVert\psi\rVert_2^2$ and is linear in the vector
length. Pure relative entropy is the Shannon entropy of `abs2.(ψ)`. Matrix
relative entropy requires a dense Hermitian eigendecomposition and therefore
has cubic time and quadratic workspace in the matrix dimension.

## QETLAB migration and reviewed discrepancy

`MATLABCompat.L1NormCoherence`, `MATLABCompat.RelEntCoherence`,
`MATLABCompat.CoherenceRank`, and `MATLABCompat.RobkCohValue` preserve the
reviewed call shape where it agrees with the strict native contracts. The
relative-entropy wrapper keeps QETLAB's base-two default; the native API
requires `base`. `RobkCohValue` preserves the pinned two outputs while
intentionally correcting unsorted, complex-phase, and nonnormalized inputs.

The pinned `CoherenceRank.m` documentation defines the number of nonzero basis
coefficients, but its implementation increments the result for coefficients
at or below the tolerance. The Julia-native function and compatibility entry
point implement the documented mathematical definition and intentionally do
not reproduce that zero-counting bug.

The committed source-free Octave 11.3.0/QETLAB artifact records six
deterministic fixtures and passes 25 assertions. Two fixtures preserve the
actual buggy upstream `CoherenceRank` outputs as reviewed discrepancy evidence;
they are not used as mathematical expected values. The fixture SHA-256 is
`11bcaaee88fac8a595e9a4eff164432dbaa4e141cdd26554da2522e22981811a`.
Octave evidence is function-specific and supplemental; MATLAB was not run.

The solver-backed and theorem-based `k`-incoherence, absolute incoherence,
robustness, trace-distance, and generalized robustness APIs are documented
separately in
[Coherence criteria and optimization](coherence_optimization.md).

The separate pure-state robustness artifact contains eight theorem-domain
agreement fixtures and three reviewed correction cases. The focused suite
passes 82 assertions and the artifact comparator passes 58 on Julia 1.12.6 and
Julia 1.10.11. It passes the
normalized real inputs through the pinned QETLAB routine, while retaining the
unsafe unsorted, complex, and nonnormalized results only as discrepancy
evidence. Its SHA-256 is
`551f529a86960a43e8765e2934007a1a640b86f8327a2f19f79860d9d5ce51a6`.
