# ADR 0006: rectangular operator spaces and two-sided map sums

- Status: Accepted and implemented
- Date: 2026-07-29

## Context

The current map types model linear maps between square matrix algebras,
including unequal input and output Hilbert-space dimensions. That is enough for
maps from `m × m` matrices to `p × p` matrices, but it is not enough for the
full pinned QETLAB contracts. In particular, the reviewed `ApplyMap`,
`ChoiMatrix`, `DualMap`, and `PartialMap` sources also accept:

- maps from `m × n` matrices to `p × q` matrices;
- paired left/right operator data; and
- independent row and column subsystem dimensions.

Calling a pair `(Aᵢ, Bᵢ)` a set of Kraus operators would incorrectly imply
complete positivity. Adding independent exceptions to the existing channel
functions would also leave Choi, superoperator, dual, partial-map, and
diagnostic conventions inconsistent.

The pinned sources were checked at revision
`d8589610f00cff106537268dee2e2a1153f3a601`. Their paired form acts as

```math
\Phi(X) = \sum_i A_i X B_i^\dagger.
```

For a completely positive map, the right factors equal the left factors.
Canonical recovery uses an eigendecomposition for Hermiticity-preserving maps
and an SVD for general maps.

## Decision

### Operator-space dimensions

Introduce one immutable, validated descriptor for

```math
\Phi : \mathbb{C}^{m\times n} \longrightarrow
       \mathbb{C}^{p\times q}.
```

It stores four positive checked dimensions:

- input rows `m`;
- input columns `n`;
- output rows `p`; and
- output columns `q`.

Public accessors will expose input and output matrix sizes. Existing scalar
`input_dimension` and `output_dimension` behavior remains available for maps
between square matrix algebras; it will not silently discard a distinct column
dimension.

### Representation types

Keep `KrausRepresentation` as the completely positive specialization

```math
\Phi(X) = \sum_i K_i X K_i^\dagger
```

between square matrix algebras. Add a separately named general operator-sum
representation for

```math
\Phi(X) = \sum_i A_i X B_i^\dagger.
```

For each term, `Aᵢ` has size `p × m` and `Bᵢ` has size `q × n`. Constructors
validate all axes, shapes, finite entries, nonempty paired collections, and
checked dimensions. They copy and read-only-wrap the collections, while
preserving dense or sparse matrix storage. Internal validated constructors use
the package's existing constructor token.

The general representation is not called Kraus and carries no positivity,
Hermiticity-preservation, or channel claim. A verified equal-pair operator sum
may be converted to `KrausRepresentation`; malformed or merely numerically
similar pairs are not silently promoted.

`ChoiRepresentation` and `SuperoperatorRepresentation` will carry the same
four-dimensional descriptor. The existing square-algebra constructors remain
source compatible.

### Choi and vectorization convention

Column-major vectorization remains the only native convention. For matrix
units `Eᵢⱼ` in the `m × n` input space, define the generalized Choi matrix by

```math
J(\Phi)_{(i,a),(j,b)} = \Phi(E_{ij})_{a,b}.
```

It has size `(m*p) × (n*q)` and need not be square. A two-sided term contributes

```math
\mathrm{vec}(A_i)\mathrm{vec}(B_i)^\dagger.
```

The superoperator is the `(p*q) × (m*n)` matrix satisfying

```math
\mathrm{vec}(\Phi(X)) = S(\Phi)\mathrm{vec}(X),
```

so one paired term contributes

```math
\overline{B_i}\otimes A_i.
```

Choi/superoperator conversion is an index reshuffle and must preserve sparse
storage. Constructors require explicit dimensions whenever matrix shape alone
does not determine the operator spaces.

### Application, duals, and partial maps

Direct operator-sum application is primary. It accumulates `Aᵢ * X * Bᵢ'`
with reusable workspaces and `mul!` where supported; it does not construct a
superoperator. Explicit Choi input may use a reshuffle or a guarded
superoperator path, chosen by measured storage and allocation behavior.

The Hilbert--Schmidt dual maps `p × q` matrices to `m × n` matrices:

```math
\Phi^\dagger(Y) = \sum_i A_i^\dagger Y B_i.
```

It therefore swaps the input/output matrix spaces and adjoints both stored
factor collections. Tests must verify the Hilbert--Schmidt identity for complex
rectangular spaces, not only representation round trips.

`partial_map` will accept independent row and column subsystem layouts. The
selected row and column dimensions must match `m` and `n`; the output layouts
replace them by `p` and `q`. The implementation will use axis
reshape/permutation and block contractions, with a sparse path or an explicit
guarded limitation. It will not construct a full tensor-product
superoperator.

### Diagnostics and decompositions

Complete positivity and Hermiticity preservation are properties of maps
between square *-algebras. Therefore:

- for `m = n` and `p = q`, Hermiticity preservation is equivalent to a
  Hermitian Choi matrix and complete positivity to a positive-semidefinite Choi
  matrix;
- for genuinely rectangular operator spaces, these diagnostics return a typed
  `not_applicable` result or reject the request explicitly;
- tolerance-aware diagnostics retain residuals and raw extremal values instead
  of silently symmetrizing or clipping; and
- Boolean convenience methods are offered only away from an explicitly
  reported numerical boundary.

Canonical decomposition of a square Hermiticity-preserving map uses a Hermitian
eigendecomposition. Positive eigenvalues give equal factor pairs and negative
eigenvalues give sign-opposed pairs. A general Choi matrix uses an SVD,
returning paired left/right factors. A completely positive result may return
the existing one-list `KrausRepresentation`. Reconstruction, subspace freedom,
zero maps, tolerance thresholds, and unavoidable dense factorizations are
documented and tested.

### Migration and implementation order

Implementation proceeds as one representation change, in this order:

1. the dimension descriptor and general operator-sum type;
2. direct application and Choi/superoperator conversions;
3. Hilbert--Schmidt dual and rectangular partial-map action;
4. Hermiticity-preserving and complete-positivity diagnostics;
5. canonical CP/HP/general decompositions;
6. complementary-map semantics for supplied CP dilations; and
7. explicit-RNG random map construction with bounded normalization.

The affected QETLAB rows remain partial or deferred until their native API,
compatibility mapping, provenance, analytic/property/invalid/sparse/generic
tests, documentation, and applicable oracle evidence all pass.

## Implementation

Implemented on 2026-07-30. The package now provides the checked
`OperatorSpace` descriptor and separately named `OperatorSumRepresentation`,
rectangular Choi and superoperator conversions, direct two-sided application,
Hilbert--Schmidt duals, rectangular partial-map action, typed physicality
diagnostics, and canonical CP, Hermiticity-preserving, and general
decompositions. The implementation is covered by the corresponding
provenance, documentation, analytic, property, invalid-input, sparse, and
compatibility evidence.

The strict static completion ledger at the pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601` now records the affected public rows
as implementation-complete. This is local static repository evidence; it does
not by itself establish MATLAB/QETLAB parity, remote-CI success, performance,
API stability, release approval, or human review.

## Alternatives rejected

- **Treat paired factors as Kraus operators.** This makes a false complete
  positivity claim.
- **Represent every map only by a dense superoperator.** This loses sparse and
  matrix-free application and scales poorly for local maps.
- **Keep square Choi types and add raw-matrix overloads for each function.**
  This duplicates dimension inference and permits convention drift.
- **Use untyped tuples or two-column arrays.** This does not make ownership,
  dimensions, factor roles, or mathematical status explicit.
- **Return `false` for rectangular CP/HP queries.** The property is not
  applicable, not disproved.

## Consequences

- Existing CP channel construction remains available and becomes an explicit
  specialization of the broader map model.
- Some public structs need an internal field migration, but existing
  constructors for square matrix algebras can remain compatible.
- General Choi matrices can be rectangular, so code may no longer infer that a
  value called a Choi representation is Hermitian or physical.
- `KrausOperators` compatibility requires a typed paired-factor result for
  non-CP maps.
- Sparse preservation and densification budgets become representation-level
  obligations rather than per-wrapper conventions.
- The decision unblocks the nine first-priority map rows; it does not mark any
  of them complete by itself.
