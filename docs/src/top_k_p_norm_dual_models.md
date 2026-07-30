# Dual top-`k`, `p`-norm model expressions

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

The numeric function
[`top_k_p_norm_dual`](@ref) evaluates the dual of the top-`k`, `p`-norm
directly from vector magnitudes or matrix singular values. The separate
[`top_k_p_norm_dual_epigraph`](@ref) API provides the missing
model-expression capability: it creates a solver-neutral convex epigraph for
a [`ComplexAffineMatrix`](@ref).

## A small affine example

This example describes a complex $2 \times 2$ matrix with eight real model
coordinates and constructs the dual top-two Euclidean-norm atom. It does not
create a solver or mutate a model.

```jldoctest
julia> X = complex_affine_variable(:X, 2, 2);

julia> atom = top_k_p_norm_dual_epigraph(X, 2, 2);

julia> (atom.k, atom.p, atom.conjugate_order, atom.singular_value_count)
(2, 2, 2.0, 2)

julia> (atom.psd_block_count, atom.auxiliary_variable_count)
(4, 41)
```

After loading JuMP and a conic optimizer, call
`add_top_k_p_norm_dual_epigraph!(model, atom, coordinates;
allow_densify=true)`. The returned handle has an `epigraph` field that can
appear in an objective or another constraint. No optimizer is selected
implicitly.

## What the lifted model means

Let $s_1 \ge \cdots \ge s_n \ge 0$ majorize the singular values of the
matrix expression, and let $q=p/(p-1)$. For $1<p<\infty$, the dual norm is the
matrix version of the vector $k$-support norm. Its homogeneous perspective
epigraph is

```math
\begin{aligned}
0 &\le z_i \le t, \\
\sum_i z_i &\le kt, \\
u_i^{1/q}z_i^{1-1/q} &\ge s_i, \\
\sum_i u_i &\le t.
\end{aligned}
```

The third line is a standard three-dimensional power cone. The singular-value
majorants use exact Ky Fan semidefinite constraints. Consequently the
formulation is composable: the input may contain model variables, rather than
being evaluated before optimization.

The endpoints simplify:

| Requested order | Dual spectral gauge |
|---|---|
| $p=1$ | $\max\{\max_i s_i,\ \sum_i s_i/k\}$ |
| $1<p<\infty$ | perspective power-cone formulation above |
| $p=\infty$ | $\sum_i s_i$ |

A one-row or one-column affine input is diagonalized internally, so the
singular values are its entry magnitudes. This matches the numeric vector
contract and avoids the pinned signed-vector ordering defect.

## Resource and failure behavior

Construction checks variable, scalar-constraint, PSD-block, PSD-dimension, and
real-block-entry limits before the modeling extension allocates a dense
expression matrix. Materialization requires
`allow_densify=true`; this is never inferred from a dense-looking input.
Unsupported optimizer cones, optimizer limits, and numerical failures remain
backend outcomes. They do not change the mathematical meaning of the atom or
produce a fabricated norm value.

The source-free QETLAB fixture covers numeric plateau and endpoint branches.
The symbolic branch is validated independently with Hypatia and SCS because
the pinned implementation depends on CVX.

```@meta
DocTestSetup = nothing
```
