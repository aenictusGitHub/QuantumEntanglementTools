# Guarded group twirls

`twirl` applies the four deterministic averaging channels exposed by the
pinned QETLAB `Twirl` entry point. It accepts a general finite square operator:
the input need not be Hermitian, positive, or trace one, and the implementation
does not silently impose any of those properties.

```julia
using LinearAlgebra
using QuantumEntanglementTools
using SparseArrays

rho = sparse([
    0.4 0.0 0.0 0.1
    0.0 0.2 0.0 0.0
    0.0 0.0 0.1 0.0
    0.1 0.0 0.0 0.3
])

rho_iso = twirl(rho; kind=:isotropic)
@assert tr(rho_iso) ≈ tr(rho)
@assert twirl(rho_iso; kind=:isotropic) ≈ rho_iso
@assert issparse(rho_iso)
```

## Four supported projections

Let `X` act on `p` equal local systems of dimension `d`, so its matrix size is
`d^p` by `d^p`. The local dimension must be an exact integer; it is never
obtained by rounding a floating-point root.

| `kind` | Copies | Fixed-point spanning family | Basis size |
|---|---:|---|---:|
| `:werner` | any positive `p` | subsystem-permutation operators | `p!` |
| `:isotropic` | exactly `2` | maximally entangled projector and its complement | `2` |
| `:real` | any positive `p` | Brauer operators indexed by perfect matchings | `(2p-1)!!` |
| `:pauli` | exactly `2`, with `d=2^q` | Bell projectors from vectorized Pauli strings | `4^q=d^2` |

The Werner map is

```math
\mathcal{T}_{\mathrm{U}}(X)
=
\int
U^{\otimes p} X \left(U^\dagger\right)^{\otimes p}
\mathrm{d}U .
```

Its fixed-point space is represented without forming tensor powers of sampled
unitaries. Instead, `twirl` projects onto the sparse permutation-operator
span. For the real branch, the same construction uses the Brauer invariant
span of the orthogonal group. The compact Gram matrices have exact integer
entries. If a spanning family is linearly dependent, exact row elimination
selects an independent subfamily before the coefficient solve.

For the isotropic branch, with normalized
$|\Phi_d\rangle=d^{-1/2}\sum_j |j,j\rangle$ and
$P_\Phi=|\Phi_d\rangle\langle\Phi_d|$,

```math
\mathcal{T}_{\mathrm{iso}}(X)
=
\langle\Phi_d|X|\Phi_d\rangle P_\Phi
+
\frac{\mathrm{tr}(X)-\langle\Phi_d|X|\Phi_d\rangle}{d^2-1}
(I-P_\Phi).
```

The one-dimensional case is handled directly, without the removable
`d^2-1` singularity. The Pauli branch is the finite bilateral average

```math
\mathcal{T}_{\mathrm{P}}(X)
=
\frac{1}{4^q}
\sum_{Q\in\mathcal{P}_q}
(Q\otimes Q)X(Q\otimes Q)^\dagger ,
```

implemented as a projection onto mutually orthogonal Bell projectors. Global
Pauli phases do not affect the result.

These maps are trace preserving, completely positive, and idempotent in exact
arithmetic. Floating-point results can contain ordinary roundoff; no
symmetrization, eigenvalue clipping, trace normalization, or idempotence repair
is performed.

## Sparse storage and explicit budgets

```julia
projected = twirl(
    rho;
    kind=:real,
    copies=2,
    max_basis_size=256,
    max_nonzeros=5_000_000,
    max_dense_entries=1_000_000,
    max_work=1_000_000_000,
)
```

Sparse input produces sparse output by default. The input itself is never
densified: permutation, Brauer, maximally-entangled, and Bell projectors are
assembled sparsely. A compact dense Gram matrix is used only in coefficient
space. To request a dense result from sparse input, set both
`sparse_output=false` and `allow_densify=true`.

The guards are checked before the corresponding enumeration or allocation:

- `max_basis_size` caps `p!`, `(2p-1)!!`, or `4^q`;
- `max_nonzeros` caps a conservative spanning-family storage estimate;
- `max_dense_entries` caps both the compact Gram/elimination matrices and an
  explicitly requested dense output;
- `max_work` caps a conservative construction and compact-solve estimate.

Each guard accepts `nothing` as an explicit opt-out. Disabling a guard does not
make an allocation feasible and does not suppress integer-overflow checks.

For a spanning-family size `b` and operator dimension `n=d^p`, the Werner and
real branches use `O(b^2)` compact Gram storage, `O(b^3)` elimination/solve
work, and at most `O(bn)` sparse construction work. The Pauli branch has
`b=n=d^2` and an `O(n^2)` stored-entry/work estimate. The isotropic branch is
linear in `n` apart from scanning a dense input.

## Exact arithmetic and numerical boundaries

Integer and rational inputs are widened to `Rational{BigInt}` so division by
Gram coefficients remains exact:

```julia
X = Rational{Int}[
    1 2 0 1
    3 5 1 0
    0 2 4 1
    1 0 3 6
]

Y = twirl(X; kind=:werner)
@assert eltype(Y) == Rational{BigInt}
@assert twirl(Y; kind=:werner) == Y
@assert tr(Y) == tr(X)
```

`Float32`, `Float64`, their complex counterparts, and generic precision
supported by Julia's dense linear solve retain their precision. There is no
pseudoinverse threshold: dependence is decided from the exact integer Gram
matrix. A singularity or unsupported scalar operation is reported rather than
hidden by down-conversion.

## Migration from QETLAB

The compatibility spelling keeps QETLAB's positional defaults,
case-insensitive type strings, and sparse result storage:

```julia
using QuantumEntanglementTools.MATLABCompat

Y = Twirl(Matrix(rho), "PaUlI", 2)
@assert issparse(Y)
```

The native spelling deliberately uses strict lowercase symbols:

```julia
Y = twirl(Matrix(rho); kind=:pauli, copies=2)
```

There are two intentional validation corrections:

- `:isotropic` and `:pauli` require exactly two copies. The pinned routine
  rejects only values greater than two and otherwise reaches dimensionally
  inconsistent formulas.
- Matrix size, integer local dimension, Pauli power-of-two dimension, finite
  entries, and all resource limits are validated before construction.

This implementation is source-informed by QETLAB `Twirl.m` at revision
`d8589610f00cff106537268dee2e2a1153f3a601`, whose recorded SHA-256 is
`a2ed4de3c937cb3e451d0b8e84508cdbe4b3bc23013558015d0e2b972853cd33`.
QETLAB is BSD-2-Clause licensed; see the repository license ledger.

The committed source-free Octave 11.3.0 fixture checks all four two-copy
branches plus the three-copy Werner and real branches. Its 18 native,
compatibility, and idempotence comparisons pass with SHA-256
`acefd7ad638fa8e0eb0e30812e975f53c864de85e8abc990c88c832d00dc9e77`.
This is function-specific supplemental evidence against the pinned QETLAB
checkout. MATLAB was not run, and no general Octave/MATLAB equivalence is
claimed.
