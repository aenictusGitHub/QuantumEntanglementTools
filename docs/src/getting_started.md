# Getting started

## Requirements

The minimum supported Julia version is 1.10. Development currently targets the
repository checkout; the package is not registered and has no stable release.

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The package test entry point executes the three published tutorials. They can
also be run directly:

```sh
julia --startup-file=no --project=. tutorials/runtests.jl
```

To use the checkout from another Julia environment:

```julia
using Pkg
Pkg.develop(path="/path/to/QuantumEntanglementTools")
using QuantumEntanglementTools
```

## Locally tested deterministic operations

```julia
using QuantumEntanglementTools
using LinearAlgebra

ψ = bell_state()
ρA = partial_trace(ψ, (2, 2); trace_out = (2,))
@assert isapprox(ρA, Matrix{ComplexF64}(I, 2, 2) / 2)
```

This behavior passes the local Tier A/Tier B test suite. It is not a claim of
full QETLAB parity or cross-platform validation.

Other implemented Tier A families include tensor products/sums, subsystem
permutations and swaps, partial transpose, realignment/inverse realignment,
symmetric and antisymmetric projectors/bases, reusable plans, and a
`MATLABCompat` namespace. Consult the generated API reference for signatures and
the status ledger before relying on a migration mapping.

Tier B also supplies operator bases and named-state constructors:

```julia
X = pauli(:X)
F = fourier_matrix(4)
ghz = ghz_state(2, 3)
ρ = horodecki_state(0.3; dims = (3, 3))
```

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
oracle evidence, and deliberately deferred functions.

The core supports Julia 1.10. The optional exact-version
[EntanglementDetection.jl extension](entanglement_detection_extension.md)
currently resolves on Julia 1.11 or later because of its backend dependency
graph; it is not installed by the core workflow above.

## Build these docs

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The first command records the local package path in the disposable docs
environment. `docs/Manifest.toml` is intentionally not versioned.

## Before relying on a function

Check all of the following:

1. its inventory row is reviewed and does not say merely `implemented`;
2. `PROVENANCE.toml` records its specification and implementation origin;
3. the API reference describes shapes, conventions, tolerances, errors, and
   certification meaning;
4. applicable analytic, property, invalid-input, sparse/generic, and independent
   or differential tests pass;
5. any optional backend reports its version and solver status.

Until those gates exist, treat results as development output.
