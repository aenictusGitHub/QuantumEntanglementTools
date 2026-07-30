# Mathematical conventions

Status: Tier A decisions below have executable local coverage in the 345-test
subsystem suite. Channel conventions have Kraus/Choi/superoperator round-trip
and physicality coverage in the Tier C suite, and the general-map work package
adds rectangular basis-action, reshuffle, and Hilbert--Schmidt-dual properties.
Exact-candidate cross-platform and authoritative MATLAB validation remain
pending. Named-basis normalizations not listed below remain unresolved.

## Subsystems and indices

Dimensions are written `(d₁, …, dₙ)` in left-to-right tensor-factor order.
Subsystem labels are one-based. Subsystem 1 is most significant and changes
slowest; subsystem `n` is least significant and changes fastest.

For one-based component indices `(i₁, …, iₙ)`, the one-based linear basis index
is

```math
1 + \sum_{k=1}^{n} (i_k - 1) \prod_{j=k+1}^{n} d_j.
```

For example, with dimensions `(2, 3)`, the basis order is
`(1,1), (1,2), (1,3), (2,1), (2,2), (2,3)`.

Dimension entries must be positive integers whose product matches the relevant
vector length or matrix dimension. Repeated, zero, negative, or out-of-range
subsystem labels are errors. Whether an API canonicalizes unsorted unique labels
or preserves their order must be documented per operation.

Native positional kernels require axes that begin at one. Offset or otherwise
non-one-based arrays are rejected explicitly before subsystem, channel, state,
random-object, or matrix-analysis indexing begins. Compatibility methods that
only copy a caller's array may preserve its axes without indexing it.

## Tensor products and vectorization

Tensor factors use the same left-to-right convention as `kron(A, B)`.
Vectorization is Julia/MATLAB column-major `vec`: rows vary fastest within each
matrix column. Consequently, the intended identity is

```math
\mathrm{vec}(A X B) = (B^\mathsf{T} \otimes A)
\mathrm{vec}(X).
```

The transpose on `B` is ordinary transpose, not adjoint.

## Partial trace

Reducing a pure vector means forming the mathematical projector `ψψ†` and
returning an operator matrix for the retained systems. The implementation should
avoid materializing the full projector when a direct contraction is available.

Tracing all subsystems returns an explicit `1 × 1` matrix; local tests cover both
pure-vector and matrix inputs.

Partial trace does not silently normalize the input.

## Partial transpose

Partial transpose swaps row and column indices for each selected subsystem. It
does not conjugate values. Applying it twice to the same systems must return the
original mathematical operator, subject only to documented numerical behavior.

## Adjoint and transpose

Use `adjoint` (`'`) for bra/ket and Hermitian-conjugation operations. Use
`transpose` only where the mathematical definition is non-conjugating, such as
the selected-index operation above. Complex inputs make this distinction
observable and therefore mandatory in tests.

## Choi matrices

For a map from `m × n` matrices to `p × q` matrices, let `Eᵢⱼ` denote the
`m × n` matrix unit. The generalized, unnormalized Choi representation is

```math
J(\Phi)_{(i,a),(j,b)}=\Phi(E_{ij})_{a,b}.
```

It has size `(m*p) × (n*q)`, with input index first and output index
second. For maps between square matrix algebras this is equivalently

```math
J(\Phi) = \sum_{i,j} |i\rangle\langle j| \otimes
          \Phi(|i\rangle\langle j|).
```

Under this convention:

- a two-sided term `A * X * B'` contributes `vec(A) * vec(B)'`;
- the matching superoperator term is `kron(conj(B), A)`;
- trace preservation means `tr_output(J) = I_input` for square algebras;
- `tr(J) = d_input` for a trace-preserving channel; and
- conversions must state all four operator-space dimensions when shape does
  not determine them.

Kraus, Choi, and superoperator round-trip and physicality tests exercise this
public contract. The complete runnable examples are in
[General maps and rectangular operator spaces](general_maps.md).

## General operator sums and duals

A general paired representation means

```math
\Phi(X)=\sum_\ell A_\ell X B_\ell^\dagger.
```

It is not called a Kraus representation and does not imply complete
positivity. The Hilbert--Schmidt inner product is

```math
\langle Y,X\rangle_{\mathrm{HS}}=\mathrm{tr}(Y^\dagger X),
```

so the adjoint map is

```math
\Phi^\dagger(Y)=\sum_\ell A_\ell^\dagger Y B_\ell.
```

Complete positivity, Hermiticity preservation, trace preservation, and
unitality are channel-algebra diagnostics. They are not silently extended to
genuinely rectangular matrix spaces.

## Normalization and tolerance

The default maximally entangled state is normalized. Normalizations for Pauli,
generalized Pauli, Gell-Mann, spin, Bell, graph, and other named families remain
to be fixed with their inventory slices.

No operation may silently normalize, symmetrize, replace non-finite values, or
clip eigenvalues. Any explicit roundoff repair must be bounded by a tolerance
derived from numeric type, problem scale, dimension, and conditioning, and must
be reported in result metadata.
