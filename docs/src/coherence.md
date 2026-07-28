# Coherence

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

The dependency-free coherence slice implements three basis-dependent
quantities for normalized pure vectors and validated density matrices:

- `l1_coherence`;
- `relative_entropy_coherence`;
- `coherence_rank` for pure vectors.

Inputs are never normalized, repaired, or clipped. Pure-vector paths preserve
generic real or complex floating-point precision, including `BigFloat`, and do
not densify sparse vectors. Matrix entropy and positivity validation require a
BLAS floating element type because the core deliberately has no
arbitrary-precision eigensolver dependency.

## Definitions and example

For a density matrix ``\rho`` in the selected computational basis,

```math
C_{l_1}(\rho) = \sum_{i \ne j} |\rho_{ij}|,
\qquad
C_{\mathrm{rel}}(\rho) = S(\operatorname{diag}(\rho)) - S(\rho).
```

The logarithm base for `relative_entropy_coherence` is required explicitly.
For a pure vector ``\psi``, coherence rank is the number of basis
coefficients above the requested numerical threshold.

```jldoctest coherence
julia> ψ = ComplexF64[1, 1] / sqrt(2);

julia> isapprox(l1_coherence(ψ), 1)
true

julia> isapprox(relative_entropy_coherence(ψ; base = 2), 1)
true

julia> coherence_rank(ψ)
2
```

A supplied `basis` for `coherence_rank` is a finite unitary matrix whose
columns are basis vectors. Coordinates use a linear solve, never an explicit
inverse. A sparse basis transform may become dense and therefore requires
`allow_densify=true`.

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
``(\sum_i |\psi_i|)^2-\lVert\psi\rVert_2^2`` and is linear in the vector
length. Pure relative entropy is the Shannon entropy of `abs2.(ψ)`. Matrix
relative entropy requires a dense Hermitian eigendecomposition and therefore
has cubic time and quadratic workspace in the matrix dimension.

## QETLAB migration and reviewed discrepancy

`MATLABCompat.L1NormCoherence`, `MATLABCompat.RelEntCoherence`, and
`MATLABCompat.CoherenceRank` preserve the reviewed call shape where it agrees
with the strict native contracts. The relative-entropy wrapper keeps QETLAB's
base-two default; the native API requires `base`.

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
