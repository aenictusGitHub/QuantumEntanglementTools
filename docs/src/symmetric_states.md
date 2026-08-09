# Symmetric multiqubit and multiqudit states

This page is the practical guide to the package's occupation-number tools for
permutation-symmetric states. The API works for qubits and arbitrary local
dimension, stays in the compressed symmetric basis whenever possible, and
checks combinatorial costs before allocation.

These functions are project-native additions. They complement the full-space
[`symmetric_subspace_basis`](@ref) and [`symmetric_projector`](@ref); they are
not MATLAB/QETLAB compatibility claims.

## Choose a representation

For $N$ parties of local dimension $d$, the ambient tensor space has
dimension $d^N$, while the symmetric subspace has dimension

```math
D_{d,N} = \binom{N+d-1}{N}.
```

Occupation coordinates label a basis vector by
$(n_1,\ldots,n_d)$, where $n_j$ counts local level $|j-1\rangle$ and
$\sum_j n_j=N$. The ordering matches the columns of
[`symmetric_subspace_basis`](@ref): earlier occupations descend first.

```@example symmetric-occupations
using QuantumEntanglementTools

dimension = symmetric_subspace_dimension(3, 4)
occupations = symmetric_occupations(3, 2)
round_trip = [
    symmetric_basis_occupation(symmetric_basis_index(n), 3, 2)
    for n in occupations
]

(dimension=dimension, occupations=occupations, round_trip=round_trip)
```

The exact dimension and basis indices are `BigInt`, so merely asking for a
combinatorial size cannot overflow an `Int`. Functions that allocate arrays
perform a separate representability check.

`parties=0` is a one-dimensional vacuum convention for the coordinate APIs on
this page. The older ambient `symmetric_subspace_basis` and
`symmetric_projector` constructors require at least one tensor factor.

## Generalized Dicke states

[`generalized_dicke_state`](@ref) constructs the normalized multiqudit state

```math
|D_{\boldsymbol n}\rangle =
\frac{1}{\sqrt{M_{\boldsymbol n}}}
\sum_{\mathrm{type}(x)=\boldsymbol n}|x_1\cdots x_N\rangle,
\qquad
M_{\boldsymbol n}=\frac{N!}{\prod_j n_j!}.
```

Subsystem 1 is the slowest-varying tensor factor. Trailing zeros in an
occupation are significant because they fix the local dimension.

```@example symmetric-dicke
using QuantumEntanglementTools
using LinearAlgebra
using SparseArrays

psi_qutrit = generalized_dicke_state((1, 1, 1))
psi_qubit = generalized_dicke_state((3, 2))

(
    qutrit_ambient_length=length(psi_qutrit),
    qutrit_support=nnz(psi_qutrit),
    qutrit_norm=norm(psi_qutrit),
    agrees_with_qubit_dicke=psi_qubit == dicke_state(5, 2),
)
```

`normalized=false` places `one(T)` on every word. No supplied coefficient or
state is silently normalized.

## Symmetric product states

For a local vector $z$, [`symmetric_product_coordinates`](@ref) returns the
coordinates of $z^{\otimes N}$ directly:

```math
c_{\boldsymbol n} =
\sqrt{M_{\boldsymbol n}}\prod_j z_j^{n_j}.
```

```@example symmetric-products
using QuantumEntanglementTools
using LinearAlgebra

z = ComplexF64[1, im, 2]
z /= norm(z)
c = symmetric_product_coordinates(z, 5)

(coordinate_count=length(c), norm=norm(c))
```

The local vector is never normalized by the function. Widened intermediates
avoid forming enormous multinomial coefficients as `Float64` before the final
conversion; use `BigFloat` or `Complex{BigFloat}` input when the requested
output itself is outside machine range. The `N=0` and `N=1` branches preserve
exact input types. For higher `N`, a nonzero coefficient that underflows the
requested output type raises `OverflowError` instead of silently becoming
zero.

## Collective one-body observables

[`symmetric_collective_operator`](@ref) represents
$\sum_{r=1}^N A^{(r)}$ without constructing the ambient operator. For
$a\ne b$, the local entry $A_{ab}$ moves one particle from level $b$ to
level $a$ with coefficient $\sqrt{n_b(n_a+1)}$.

```@example symmetric-collective
using QuantumEntanglementTools
using LinearAlgebra

sx = ComplexF64[0 1; 1 0] / 2
sy = ComplexF64[0 -im; im 0] / 2
sz = ComplexF64[1 0; 0 -1] / 2

Jx = symmetric_collective_operator(sx, 6)
Jy = symmetric_collective_operator(sy, 6)
Jz = symmetric_collective_operator(sz, 6)

commutator_residual = norm(Jx * Jy - Jy * Jx - im * Jz)
```

The input operator is not Hermitized. A non-Hermitian local operator produces
the corresponding non-Hermitian collective operator. The `N=0`, `N=1`, and
diagonal-only branches avoid unnecessary radical promotion; fixed-width
integer overflow is reported instead of wrapping.

## Bipartitions and reduced states

[`symmetric_split_isometry`](@ref) exposes the exact relationship between a
global occupation basis and a $k|(N-k)$ tensor product of two symmetric
bases. If $a+c=n$, its coefficient is

```math
\sqrt{\frac{M_a M_c}{M_n}}.
```

The kept factor is the slow, left tensor factor, consistent with the package's
subsystem convention.

```@example symmetric-split
using QuantumEntanglementTools
using LinearAlgebra

V = symmetric_split_isometry(3, 4, 2)
isometry_residual = norm(adjoint(V) * V - I)
```

The identity relation is exact at the `keep=0` and `keep=N` endpoints. For a
nontrivial split it holds up to the precision of `T`; construct `T=BigFloat`
inside a suitable `setprecision` block when more than the current `BigFloat`
precision is required.

[`symmetric_reduced_state`](@ref) uses the same coefficients directly and
avoids the ambient $d^N\times d^N$ density matrix. A vector input is treated
as a pure state; a matrix input is treated linearly as an operator.

```@example symmetric-reduction
using QuantumEntanglementTools
using LinearAlgebra

# The occupation (3,1) is column 2 for four qubits.
c = zeros(5)
c[2] = 1
rho1 = symmetric_reduced_state(c, 2, 4; keep=1)

(rho1=Matrix(rho1), trace=tr(rho1))
```

`keep=0` returns an explicit `1 × 1` norm-squared or trace. `keep=N` returns
the pure projector or the original operator in symmetric coordinates. Inputs
are never normalized, Hermitized, projected, or repaired. Existing high-
precision `BigFloat` input precision is retained through nontrivial reductions,
and nonfinite output caused by a narrow floating type is reported.

## Maximally mixed symmetric state

[`symmetric_maximally_mixed_state`](@ref) returns $I_D/D$ in compressed
coordinates by default. Select `representation=:ambient` explicitly to obtain
$P_{\mathrm{sym}}/D$ in the full tensor space.

```@example symmetric-mixed
using QuantumEntanglementTools
using LinearAlgebra

rho_coordinates = symmetric_maximally_mixed_state(3, 3)
rho_ambient = symmetric_maximally_mixed_state(
    2,
    3;
    representation=:ambient,
    T=Rational{BigInt},
)

(
    coordinate_size=size(rho_coordinates),
    coordinate_trace=tr(rho_coordinates),
    ambient_size=size(rho_ambient),
    ambient_trace=tr(rho_ambient),
)
```

This state is maximally mixed only inside the symmetric subspace. The
constructor makes no separability claim. Rational output has exact unit trace;
floating output stores the rounded value `1/D`, so its trace is approximate.
Accumulate very large low-precision diagonals in a wider type when checking
normalization (for example, `sum(Float64, nonzeros(rho))` for `Float16`).

## A terminology warning

Two notions that are often both called "symmetric" are different for mixed
states:

- permutation invariance means $U_\pi\rho U_\pi^\dagger=\rho$;
- symmetric-subspace support means
  $P_{\mathrm{sym}}\rho P_{\mathrm{sym}}=\rho$.

The identity operator is permutation invariant but is not supported only on
the symmetric subspace. The coordinate APIs on this page represent the second
notion. They do not silently project a general permutation-invariant operator
into that subspace.

## Resource and numeric behavior

Occupation counts, multinomial support sizes, split dimensions, candidate
nonzeros, dense entries, and conservative work estimates are computed with
`BigInt` before allocation. The main guards are:

- `max_occupations` for compressed sectors;
- `max_entries` for either the returned tuple from
  [`symmetric_basis_occupation`](@ref) or the total stored tuple fields from
  [`symmetric_occupations`](@ref);
- `max_nonzeros` for sparse construction or possible sparse output;
- `max_dense_entries` when dense output is requested;
- `max_work` for conservative scalar-operation estimates.

Set a guard to `nothing` only after reviewing the requested cost. Array-length
representability checks remain mandatory. Sparse input to
[`symmetric_reduced_state`](@ref) produces sparse output by default; ambient
construction and dense output are always explicit choices.

Operations that must materialize occupation tuples use at most 1,000 local
levels. This hard safety limit avoids compiler instability from enormous Julia
tuple types; it does not limit the particle count, and the exact dimension
helper plus the `N=0` and `N=1` endpoint paths do not need such tuples.

The exact dimension helper returns an unguarded `BigInt`; at astronomical
balanced dimensions, that scalar integer can itself be expensive. Nontrivial
split and reduction paths currently materialize guarded occupation and weight
workspaces. Sparse matrix wrappers may also be canonicalized into a temporary
sparse matrix so that their logical entries—not hidden parent storage—are
validated correctly. Review `max_work` as well as output limits for very large
jobs.

The focused tests compare occupation columns with the existing full symmetric
basis, multiqudit reductions with ambient [`partial_trace`](@ref), collective
operators with explicit tensor sums, and qubit generalized Dicke states with
[`dicke_state`](@ref). These are local analytic and property checks, not a
claim of external QETLAB/MATLAB parity.
