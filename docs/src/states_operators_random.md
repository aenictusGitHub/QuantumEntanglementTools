# States, operators, and random objects

Tier B provides a reviewed implementation slice for common operator bases,
named states, and randomized quantum objects. The `0.1.x` API remains
experimental: the evidence below applies to the exact functions and inputs
tested, not to all of QETLAB.

## Operator bases

The native constructors use zero-based basis labels where QETLAB does:

```julia
using QuantumEntanglementTools
using LinearAlgebra

X = pauli(:X)
Z = pauli(:Z)
XZ = pauli((:X, :Z); sparse_output = true)
W12 = generalized_pauli(1, 2, 3)
λ8 = gell_mann(8)
F4 = fourier_matrix(4)

@assert X * X == I
@assert size(XZ) == (4, 4)
@assert isapprox(W12' * W12, I)
@assert ishermitian(λ8)
@assert isapprox(F4' * F4, I)
```

`pauli` rejects unknown labels instead of treating them as identity.
Generalized Pauli, generalized Gell-Mann, and Fourier constructors accept a
floating-point type keyword for higher-precision output. Sparse output is
selected explicitly and multi-qubit sparse Pauli products are constructed
without first materializing the dense tensor product.

## Named states

Pure-state constructors make normalization and sparse output explicit:

```julia
using QuantumEntanglementTools
using LinearAlgebra

bell = bell_state()
ghz = ghz_state(2, 3)
w = w_state(4)
dicke = dicke_state(4, 1)

@assert norm(bell) ≈ 1
@assert norm(ghz) ≈ 1
@assert w == dicke
```

Supplied GHZ and W coefficients are used exactly; the package never silently
normalizes user data. `dicke_state` counts its combinations with `BigInt`
before enumeration. By default it rejects more than 1,000,000 nonzeros,
10,000,000 dense entries, or 100,000,000 estimated scalar operations; each
guard may be set to `nothing` only after the caller reviews the requested
cost. `MATLABCompat.DickeState` keeps sparse compatibility output and forwards
the same limits.

For multiqudits and scalable occupation-coordinate calculations, use the
project-native [symmetric-state toolkit](symmetric_states.md). It adds exact
occupation rank/unrank, generalized Dicke states, product coordinates,
collective one-body operators, symmetric bipartition isometries, direct
reduced states, and symmetric maximally mixed states. These additions do not
create new `MATLABCompat` mappings or widen QETLAB parity claims.

Mixed-state families validate their documented physical parameter ranges
instead of clipping inputs:

```julia
ρiso = isotropic_state(3, 1 // 4)
ρwer = werner_state(3, 0.2)
ρwer3 = werner_state(2, [0.05, 0.04, 0.03, 0.03, 0.02])
ρhor = horodecki_state(0.3; dims = (3, 3))

@assert tr(ρiso) == 1
@assert tr(ρwer) ≈ 1
@assert tr(ρwer3) ≈ 1
@assert tr(ρhor) ≈ 1
```

Also implemented are maximally entangled, Bell, Gisin, Breuer, Brauer, and
chessboard constructors. A Werner parameter vector must contain `p! - 1`
entries and constructs the normalized documented sum over the nonidentity
permutations of `p` equal-dimensional subsystems. Coefficients paired by
inverse permutations must be exact conjugates. Positivity is certified either
by `sum(abs, alpha) <= 1` or by an explicitly bounded dense eigendecomposition;
inputs are never symmetrized, clipped, or repaired. Permutation, sparse-entry,
dense-validation, and work guards are configurable.

The pinned multipartite implementation accidentally reinitializes its
accumulator on every permutation, so its output retains only the final
parameter. The Julia native and compatibility APIs implement the documented
sum instead. A source-free Octave 11.3.0 fixture records that upstream output
and verifies the correction; MATLAB was not run.
Because the number of Brauer columns is `(2*pairs-1)!!`,
`brauer_states` rejects more than 100,000 matchings or 1,000,000 stored entries
by default, before enumerating them. The `max_matchings` and `max_nonzeros`
keywords may be raised—or explicitly set to `nothing`—only after reviewing the
memory cost.

## Entangled subspaces

[`entangled_subspace`](@ref) constructs the pinned sparse
diagonal-Vandermonde basis for a bipartite subspace in which every nonzero
vector has Schmidt rank at least `r + 1`. For local dimensions `(d₁, d₂)`, the
sharp maximal subspace dimension is `(d₁-r) * (d₂-r)`.

```@example entangled-subspace
using QuantumEntanglementTools
using LinearAlgebra
using SparseArrays

basis = entangled_subspace(2, (2, 3); r=1, coefficient_type=BigInt)
first_coefficient_matrix = reshape(basis[:, 1], 3, 2)

(size=size(basis), sparse=issparse(basis), schmidt_rank=rank(Float64.(first_coefficient_matrix)))
```

Subsystem one is the slowest-varying tensor factor, so a basis column reshapes
to a `d₂ × d₁` coefficient matrix. The returned columns are linearly
independent but deliberately not normalized or orthogonal. Exact grid tests on
small maximal subspaces verify the promised lower Schmidt rank for every
nonzero tested linear combination, rather than checking only the generator
columns.

The default coefficient type is `Float64`, matching QETLAB. The native keyword
accepts other concrete numeric types, including `BigInt`,
`Rational{BigInt}`, `Float32`, and `BigFloat`, without an implicit conversion.
`max_nonzeros` and `max_work` preflight the sparse construction and may be
disabled only explicitly with `nothing`. `MATLABCompat.EntangledSubspace`
preserves the pinned positional argument order and sparse `Float64` output.
Four source-free fixtures generated by Octave 11.3.0 from the pinned QETLAB
revision exactly match the equal-dimension, rectangular, nonmaximal-prefix,
and `r=2` branches. This is supplemental evidence; authoritative MATLAB
execution has not been run.

## Mandatory explicit RNG

Every randomized API requires an `rng::AbstractRNG` as its first argument:

```julia
using QuantumEntanglementTools
using Random
using LinearAlgebra

rng = Xoshiro(0x514554)
p = random_probabilities(rng, 5)
ψ = random_state_vector(rng, (3, 4); schmidt_rank = 2)
ρ = random_density_matrix(rng, 4; rank = 2)
U = random_unitary(rng, 4)
A = random_graph(rng, 8; edge_probability = 0.25)
effects = random_povm(rng, 3, 4)

@assert sum(p) ≈ 1
@assert norm(ψ) ≈ 1
@assert tr(ρ) ≈ 1
@assert U' * U ≈ I
@assert A == A'
@assert sum(effects) ≈ I
```

There are intentionally no convenience methods that draw from
`Random.default_rng()`. The `MATLABCompat` randomized names also require a
leading RNG, which is a documented safety deviation from QETLAB:

```julia
rng = Xoshiro(7)
ρ = MATLABCompat.RandomDensityMatrix(rng, 4, 0, 2, "hs")
```

The test suite seeds and exercises every native and compatibility random
constructor between draws from Julia's global stream and verifies that the
global sequence is unchanged.

`random_povm` uses an independent documented Haar-isometry construction.
[`random_superoperator`](@ref) is the separate bounded random-map API: it
returns a structured result with a Kraus complete-positivity certificate,
marginal residuals, and explicit termination evidence. See
[Bounded random completely positive maps](random_superoperators.md).

[`random_ppt_state`](@ref) constructs a bipartite PPT density matrix with a
mandatory RNG and independently checked PSD, trace, partial-transpose, and rank
evidence. Full-rank requests use a bounded shifted-induced construction;
rank-limited requests use a bounded random separable mixture. The structured
result distinguishes verified construction, resource/iteration limits,
numerical failure, and verification failure, and never exposes an unverified
candidate as `state`. See [Random PPT states](random_ppt_states.md).

## Differential evidence and limits

The committed Tier B oracle artifact contains 18 deterministic operator/state
fixtures generated by Octave 11.3.0 from QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`. The comparison suite checks both
the native and `MATLABCompat` paths, yielding 72 assertions. This is
function-specific supplemental evidence only. Octave is not assumed to be
generally MATLAB-compatible, MATLAB was unavailable, and randomized streams
are validated with seeded analytic/property tests rather than unrelated
cross-language draws.

See [Reproducibility](reproducibility.md) for the metadata contract and
[Migration from QETLAB](migration_from_qetlab.md) for argument-level mappings.
