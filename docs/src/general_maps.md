# General maps and rectangular operator spaces

`QuantumEntanglementTools` distinguishes a general linear map from a quantum
channel. A map may act between rectangular matrix spaces,

```math
\Phi : \mathbb{C}^{m\times n} \longrightarrow
       \mathbb{C}^{p\times q},
```

whereas complete positivity is meaningful only when the domain and codomain
are square matrix algebras.

## Describe the operator space

`OperatorSpace(m, n, p, q)` records all four dimensions. The tuple constructor
is equivalent:

```@example general_maps
using QuantumEntanglementTools
using LinearAlgebra

space = OperatorSpace((2, 3), (3, 4))
@assert input_size(space) == (2, 3)
@assert output_size(space) == (3, 4)
nothing
```

`input_dimension` and `output_dimension` remain available for square matrix
spaces. They throw `ArgumentError` for a genuinely rectangular space rather
than discarding a distinct column dimension.

All dimensions must be positive integers. Constructors also check the products
needed by Choi and superoperator storage for `Int` overflow.

## General two-sided sums

`OperatorSumRepresentation(left, right)` represents

```math
\Phi(X)=\sum_{\ell=1}^{r} A_\ell X B_\ell^\dagger.
```

For the operator space above, each `A` has size `3 × 2` and each `B` has size
`4 × 3`. The constructor validates matching nonempty collections, one-based
axes, shapes, finite numeric entries, and concrete element types. It copies
the factors and preserves sparse factor storage. Use `operator_sum_factors`
to obtain independent mutable copies.

```@example general_maps
A = [1.0 0.0; 0.0 1.0; 1.0 -1.0]
B = [1.0 0.0 0.0; 0.0 1.0 1.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
map = OperatorSumRepresentation([A], [B])
X = reshape(1.0:6.0, 2, 3)

@assert apply_channel(X, map) == A * X * B'
@assert operator_space(map) == space
nothing
```

This type makes no positivity or channel claim. A
`KrausRepresentation([K₁, ...])` is the completely positive specialization
with identical left and right factors.

Direct `apply_channel(X, map)` dispatches on the typed representation.
Homogeneous strided operator-sum factors use a preallocated `mul!` kernel with
one reusable intermediate and output accumulator. Sparse or mixed-storage
factors use ordinary Julia matrix products so sparse structure is not
discarded. Neither path constructs a global superoperator.

For migration, `MATLABCompat.ApplyMap(X, PHI)` accepts the pinned raw forms:

- a generalized Choi matrix;
- a tuple/vector or one-column matrix of CP Kraus operators;
- a one-row matrix containing more than two CP Kraus operators; or
- an $r \times 2$ factor-cell analogue whose rows contain
  $(A_\ell,B_\ell)$.

Ambiguous or wider factor matrices are rejected explicitly. All valid forms
delegate to `apply_channel`.

## Choi and superoperator forms

For matrix units `Eᵢⱼ` in the input space, the generalized Choi matrix uses
input-first ordering:

```math
J(\Phi)_{(i,a),(j,b)} = \Phi(E_{ij})_{a,b}.
```

It has size `(m*p) × (n*q)`. The superoperator has size
`(p*q) × (m*n)` and is defined by

```math
\mathrm{vec}(\Phi(X))=S(\Phi)\mathrm{vec}(X).
```

A paired term contributes
`vec(A) * vec(B)'` to the Choi matrix and `kron(conj(B), A)` to the
superoperator. Conversions are index reshuffles and preserve sparse storage.

```@example general_maps
choi = choi_representation(map)
superop = superoperator_representation(map)

@assert size(choi_matrix(choi)) == (6, 12)
@assert size(superoperator_matrix(superop)) == (12, 6)
@assert vec(apply_channel(X, map)) == superoperator_matrix(superop) * vec(X)
@assert choi_matrix(choi_representation(superop)) == choi_matrix(choi)
nothing
```

Use `ChoiRepresentation(matrix, space)` or
`SuperoperatorRepresentation(matrix, space)` when a raw matrix represents a
rectangular operator space. The one-argument constructors infer only square
matrix-algebra dimensions.

`MATLABCompat.ChoiMatrix(PHI, SYS)` preserves the two QETLAB orderings for
factor input. `SYS=2` is the package's input-first convention. `SYS=1`
independently swaps the row and column tensor factors, which remains valid
when those factors have unequal dimensions. As in the pinned source, a numeric
matrix is treated as an already formed Choi matrix and copied before `SYS` is
examined.

## Hilbert--Schmidt dual

The package uses

```math
\langle Y,\Phi(X)\rangle_{\mathrm{HS}}
=\mathrm{tr}\!\left(Y^\dagger\Phi(X)\right).
```

For a paired term, the dual acts as
`A' * Y * B`. It swaps the input and output matrix spaces and is available for
operator-sum, Choi, and superoperator representations.

```@example general_maps
dual = dual_channel(map)
Y = reshape(ComplexF64.(1:12), 3, 4)
Xc = ComplexF64.(X)

@assert dot(vec(Y), vec(apply_channel(Xc, map))) ≈
        dot(vec(apply_channel(Y, dual)), vec(Xc))
nothing
```

`MATLABCompat.DualMap` preserves the representation style supplied by
migration code. It adjoints every raw factor-cell entry without changing that
cell array's shape. For a raw Choi matrix with unequal or rectangular spaces,
pass QETLAB's `DIM`: a two-entry vector for square algebras or a $2 \times 2$
matrix whose first row contains input/output row dimensions and whose second
row contains input/output column dimensions.

## Apply a map to one subsystem

The five-argument native form

```julia
partial_map(X, map, subsystem, row_dims, column_dims)
```

allows the selected map to change both the row and column dimension of one
local matrix space. It permutes the chosen local row and column axes into a
block position, applies the map blockwise, and reverses the permutation. It
does not form $I \otimes S(\Phi)$. The four-argument overload uses the same
layout for rows and columns.

`MATLABCompat.PartialMap(X, PHI, SYS, DIM)` additionally accepts QETLAB's
scalar, dimension-vector, and two-row rectangular `DIM` forms plus the raw map
forms listed above. A scalar `DIM` is valid only for a square input whose
dimension it divides exactly. Fully sparse factor and input data remain
sparse.

```@example general_maps
first = reshape(1.0:6.0, 2, 3)
second = reshape(2.0:7.0, 2, 3)
local_map = OperatorSumRepresentation(
    [[1.0 0.0; 0.0 1.0; 1.0 -1.0]],
    [[1.0 0.0 0.0; 0.0 1.0 1.0; 1.0 0.0 -1.0; 0.0 1.0 0.0]],
)

mapped = partial_map(kron(first, second), local_map, 2, (2, 2), (3, 3))
@assert mapped == kron(first, apply_channel(second, local_map))
nothing
```

## Hermiticity preservation

For maps between square matrix algebras, [`is_hermiticity_preserving`](@ref)
uses the theorem that $\Phi$ preserves Hermiticity exactly when its Choi
matrix is Hermitian. It returns a [`MatrixPredicateResult`](@ref): satisfied,
violated with a decisive-entry witness, or unknown when a nonzero defect lies
inside the requested numerical tolerance. Exact integer/rational data are
decided exactly. A genuinely rectangular matrix space is rejected as
non-applicable rather than mislabeled `false`.

`MATLABCompat.IsHermPreserving` accepts raw Choi and factor data and uses
QETLAB's default absolute tolerance $\epsilon^{3/4}$. Unlike QETLAB's
Boolean comparison, a nonzero boundary residual remains `unknown`; no input is
silently symmetrized. A nonsquare raw Choi matrix retains the pinned negative
conclusion because it cannot be Hermitian.

## Complete positivity

[`is_completely_positive`](@ref) returns a
[`MatrixPredicateResult`](@ref), not a Boolean. A supplied
[`KrausRepresentation`](@ref) is satisfied by construction. Other
representations must first have an exactly Hermitian Choi matrix outside any
numerical boundary and must then pass the structured positive-semidefiniteness
test:

```@example general_maps
identity_map = KrausRepresentation([Matrix{Float64}(I, 2, 2)])
cp_result = is_completely_positive(identity_map)

@assert cp_result.status === MatrixPredicateSatisfied
nothing
```

A robust Choi-Hermiticity or negative-eigenvalue witness is
`MatrixPredicateViolated`. A nonzero defect or negative eigenvalue inside the
requested tolerance is `MatrixPredicateUnknown`: it is not silently promoted
to a completely positive map. Exact integer and rational Choi matrices are
decided exactly and require zero tolerances. Sparse spectral work requires
`allow_densify=true`.

`MATLABCompat.IsCP(PHI, TOL)` accepts raw Choi matrices and one- or two-column
factor data using QETLAB's scalar absolute tolerance. It retains the same
structured three-way result, so migration code must inspect `.status` instead
of using the result as a Boolean.

## Recover canonical factors

`operator_sum_decomposition` computes a thin SVD of the Choi matrix and
returns an `OperatorSumDecompositionResult`. The result retains the singular
values, cutoff, retained rank, and Frobenius norm of discarded singular
values; truncation is therefore visible rather than a silent repair.

The standard-library SVD path accepts `Float32`, `Float64`, and their complex
counterparts. It rejects exact and arbitrary-precision matrices instead of
downcasting them. Sparse Choi data requires `allow_densify=true`.

```@example general_maps
decomposition = operator_sum_decomposition(choi; atol=0, rtol=0)
reconstructed = choi_matrix(decomposition.representation)

@assert reconstructed ≈ choi_matrix(choi)
@assert decomposition.discarded_frobenius_norm == 0
nothing
```

[`canonical_map_decomposition`](@ref) reproduces the three mathematical
branches behind QETLAB's `KrausOperators` while keeping their meanings
explicit:

1. established complete positivity gives equal left/right factors;
2. established Hermiticity preservation gives positive pairs followed by
   signed negative pairs; and
3. a general map gives paired factors from the unmodified Choi SVD.

Its [`CanonicalMapDecompositionResult`](@ref) records the branch,
structured CP and Hermiticity diagnostics, cutoff, retained rank, discarded
Frobenius norm, and actual Choi reconstruction residual. Only the first branch
is a Kraus representation. Accordingly, [`kraus_operators`](@ref) returns
matrices only when complete positivity is established; otherwise it throws
`DomainError` and directs callers to the paired result.

`MATLABCompat.KrausOperators` preserves QETLAB's raw output shapes: a vector
for the CP branch and an $r \times 2$ matrix for the signed-Hermitian or
general branch. The wrapper supports the full rectangular
`DIM=[input_rows output_rows; input_columns output_columns]` case; this also
corrects the pinned `pad_array` concatenation failure for unequal Choi row and
column sizes. Set `diagnostics=true` to obtain the native structured result.

The decomposition is nonunique: singular subspaces may rotate and individual
factor pairs may change phase. Validate reconstruction and subspaces, not
entrywise factor equality.

## Complementary maps

[`complementary_channel`](@ref) preserves a supplied dilation. If
$K_r : \mathbb{C}^{d_\mathrm{in}}\to\mathbb{C}^{d_\mathrm{out}}$ are the
supplied Kraus matrices, its factors satisfy
$L_a[r,i]=K_r[a,i]$; redundant Kraus matrices therefore produce a larger,
but intentional, environment.

For paired factors, $L_a[r,i]=A_r[a,i]$ and
$R_a[r,j]=B_r[a,j]$. The input operator space may be rectangular, while the
output must be square so the shared output index $a$ is defined. Unequal
output row and column dimensions are rejected explicitly. This is an
intentional correction of the pinned routine's failing
`cell2mat`/`mat2cell` branch, not an omitted valid complement.

```@example general_maps
redundant_dilation = KrausRepresentation([
    Matrix{Float64}(I, 2, 2) / sqrt(2),
    Matrix{Float64}(I, 2, 2) / sqrt(2),
])
complement = complementary_channel(redundant_dilation)

@assert output_dimension(complement) == 2
@assert choi_matrix(complementary_channel(redundant_dilation)) ==
        choi_matrix(complement)
nothing
```

Choi and superoperator inputs first use the canonical factorization and retain
their representation kind. `MATLABCompat.ComplementaryMap` returns raw factor
data for raw factor input and a raw Choi matrix for numeric Choi input.

## Applicability and storage

`is_completely_positive`, `is_trace_preserving`, `is_unital`,
`kraus_representation`, and quantum-channel operations require square input
and output algebras. They reject a genuinely rectangular operator space with a
precise `ArgumentError`; returning `false` would incorrectly describe a
non-applicable property as disproved. Paired complementary maps are the
exception described above: they allow rectangular input spaces but require a
square output space.

Direct application of an operator sum does not materialize its superoperator.
Explicit conversion may allocate an object with `(p*q) * (m*n)` entries.
Sparse Choi/superoperator reshuffles preserve sparse storage, while
factorization-based decompositions have an explicit densification opt-in.

See [Mathematical conventions](conventions.md),
[Sparse and large-scale computations](sparse_and_large_scale.md), and
[ADR 0006](https://github.com/aenictusGitHub/QuantumEntanglementTools/blob/main/docs/adr/0006-general-operator-space-maps.md)
for the governing
contracts.
