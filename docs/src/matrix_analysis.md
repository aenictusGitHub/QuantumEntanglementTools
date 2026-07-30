# Matrix analysis

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

The matrix-analysis slice provides standard strong majorization, elementary
symmetric polynomials, multiplicative and additive compound matrices, and
bounded numerical commutants. Its contracts prefer mathematical shape
conventions, exact arithmetic where possible, and explicit sparse behavior over
reproducing incidental behavior of the pinned MATLAB implementation.

## Strong majorization

`majorizes(a, b)` compares finite real vectors after padding both inputs with
zeros to a common length and sorting them in descending order. It checks every
proper prefix and also requires equal total sums within one symmetric
tolerance:

```math
\mathrm{atol}
+ \mathrm{rtol}
  \max\!\left(\lVert a\rVert_1,\lVert b\rVert_1,1\right).
```

Exact vector inputs default to zero relative tolerance. Floating inputs default
to `sqrt(eps(T))`; callers can select `atol` and `rtol` explicitly. Matrix
inputs are compared through their singular values, including a `1×n`
`AbstractMatrix`. Sparse inputs are rejected because singular-value work and
sorting could otherwise densify them implicitly. Prefix arithmetic widens
integer and rational inputs to `BigInt` and `Rational{BigInt}` respectively,
while homogeneous `Float32` and `Float64` inputs retain their floating type.

```jldoctest matrix-analysis
julia> majorizes([4, 1, 1], [3, 2, 1]; rtol = 0)
true

julia> majorizes([2, 0], [1, 0]; rtol = 0)
false

julia> majorizes([1, -1], [1, 0, -1]; rtol = 0)
true
```

`MATLABCompat.Majorizes` preserves the pinned `Majorizes.m` convention instead:
it implements weak majorization without equal-total validation, treats row and
column matrices as vectors, sorts before appending padding zeros, and uses the
one-sided allowance
`atol + rtol * norm(A_values)`. Its default `rtol=eps(Float64)^(3/4)` matches
the pinned routine; both tolerances are exposed explicitly. Sparse inputs
require `allow_densify=true`. Its prefix accumulation uses the same
exact-widening and floating-type-preservation policy as the native path.

Consequently, `MATLABCompat.Majorizes([2,0], [1,0]; rtol=0)` returns `true`
while native `majorizes` returns `false`. The compatibility wrapper also
reproduces the reviewed negative-padding and one-row-vector results. Callers
should choose the native API for the standard strong definition.

## Elementary symmetric polynomials

`elementary_symmetric_polynomial(values, order)` uses a descending dynamic
program with `O(length(values) * order)` arithmetic and `O(order)` storage.
It accepts sparse vectors without densifying them, uses `e₀ = 1`, and retains
integer, rational, and compatible complex arithmetic. Intermediate exact
arithmetic is widened; narrowing that cannot represent the result raises an
`OverflowError`.

```jldoctest matrix-analysis
julia> elementary_symmetric_polynomial([1, 2, 3, 4], 2)
35

julia> elementary_symmetric_polynomial(
           Rational{Int}[1//2, 2//3, 3//4],
           2,
       )
29//24
```

The reviewed upstream file and function are named `ElemSymPoly.m` and
`ElemSymPoly`. `MATLABCompat.ElemSymPoly` preserves that spelling, accepts
Julia vectors and MATLAB-shaped row or column matrices, and delegates to the
native exact implementation. Exact result types can therefore differ from the
pinned floating-point projection path.

## Multiplicative compounds

For an `m×n` matrix, `compound_matrix(A, k)` returns a
`binomial(m,k) × binomial(n,k)` matrix of minors in lexicographic subset order.
Order zero is the `1×1` multiplicative identity. If `k` exceeds only one input
dimension, the zero dimension is retained rather than collapsing the result:

```jldoctest matrix-analysis
julia> compound_matrix([1 2; 3 5; 7 11], 2)
3×1 Matrix{Int64}:
 -1
 -3
 -2

julia> size(compound_matrix(zeros(Int, 2, 3), 3))
(0, 1)
```

Pinned QETLAB returns `0×0` whenever `k > min(m,n)`, and
`MATLABCompat.CompoundMatrix` preserves that compatibility shape. The native
function retains the more informative mathematical shape. For orders where
values exist, the wrapper delegates to the native implementation: dense inputs
return dense matrices and sparse inputs return sparse matrices by default,
with `sparse_output` overriding either choice. Exact minors use widened integer
or rational arithmetic with checked narrowing rather than the floating
projection path used for ordinary numeric QETLAB inputs.

## Additive compounds

`additive_compound_matrix(A, k)` requires a square matrix and represents the
induced action on the `k`th exterior power. Order zero is the `1×1` additive
zero, order one is a copy of `A`, and an order above the dimension returns a
`0×0` matrix. The same explicit sparse-output and checked exact-arithmetic
policy applies.

```jldoctest matrix-analysis
julia> additive_compound_matrix([1 2; 3 4], 0)
1×1 Matrix{Int64}:
 0

julia> additive_compound_matrix([1 2; 3 4], 2)
1×1 Matrix{Int64}:
 5
```

The pinned `AdditiveCompoundMatrix.m` dependency path errors at order zero;
`MATLABCompat.AdditiveCompoundMatrix` reproduces that reviewed error. Use the
native function for the mathematically defined `1×1` result.

## Matrix commutants

For one square matrix or a nonempty collection of equally sized square
matrices, `commutant(generators)` returns matrices `X₁, …, Xᵣ` spanning every
matrix that commutes with all generators. The basis is orthonormal in the
Hilbert--Schmidt inner product. With column-major vectorization, the
implementation stacks the equations

```math
\left(I_n \otimes A - A^\mathsf{T} \otimes I_n\right)
\mathrm{vec}(X) = 0
```

and obtains their numerical null space by SVD. The ordinary transpose in this
equation is intentional, including for complex generators.

```jldoctest matrix-analysis
julia> A = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 4.0];

julia> basis = commutant(A);

julia> length(basis)
3

julia> maximum(maximum(abs, A * X - X * A) for X in basis) < 1e-12
true
```

The basis itself is not unique: signs, complex phases, and unitary mixing
inside a degenerate null space may change between factorization libraries.
Tests and downstream code should compare its dimension, commutator residuals,
or the projector formed from `hcat(vec.(basis)...)`, never basis entries.
`atol` and `rtol` control the numerical rank threshold explicitly.

Sparse generators produce a sparse Kronecker commutator first, but the
dependency-free core has no sparse null-space factorization. They therefore
require `allow_densify=true`, after which `max_entries` and `max_work` still
guard the dense allocation and SVD. Dense inputs are subject to the same
budgets. The native result is dense because a commutant basis is generally
dense; `MATLABCompat.Commutant` converts the basis back to sparse storage when
all generators were sparse, preserving the pinned storage convention at the
cost of potentially inefficient sparse matrices.

The supported factorization types are `Float32`, `Float64`, `ComplexF32`, and
`ComplexF64`; ordinary integer inputs promote to `Float64`. `BigFloat` and
`Complex{BigFloat}` are rejected instead of being silently narrowed because
the dependency-free environment has no generic-precision SVD. With `g`
generators of size `n×n`, dense storage is `O(g n⁴)` and the conservative SVD
work estimate is `O(g n⁶)`.

Full signatures are collected in the [API reference](api/index.md).
