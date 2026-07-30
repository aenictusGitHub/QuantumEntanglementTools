# Top-k p-norms in optimization models

For a matrix $X$ with singular values
$s_1(X)\geq\cdots\geq s_r(X)\geq0$, the top-$k$, $p$-norm is

```math
\lVert X\rVert_{(k,p)}
=
\left(\sum_{j=1}^{\min(k,r)}s_j(X)^p\right)^{1/p},
\qquad 1\leq p<\infty.
```

For $p=\infty$, it is $s_1(X)$. For vectors, the same definition uses the
largest entry magnitudes rather than the single singular value of a row or
column matrix.

## Numeric arrays

Use [`top_k_p_norm`](@ref) when the input is known:

```julia
using QuantumEntanglementTools

top_k_p_norm([3.0, -4.0, 2.0], 2, 2)
```

Vector entries are ordered by magnitude. This intentionally corrects the
pinned QETLAB numeric branch, which sorts signed vector entries and can return
a value that is not a norm.

## Solver-neutral affine variables

Use [`ComplexAffineMatrix`](@ref) and
[`top_k_p_norm_epigraph`](@ref) when the input is part of an optimization
model:

```julia
X = complex_affine_variable(:X, 2, 3)
atom = top_k_p_norm_epigraph(X, 2, 3)
```

`ComplexAffineMatrix` stores arbitrary rectangular complex affine data over
real coordinates. It is distinct from [`HermitianAffineMatrix`](@ref), which
is square and Hermitian by construction. A complex variable uses real entries
first in column-major order, then imaginary entries.

The atom is package-owned data. It records the Hermitian dilation, clipped
$k$, norm order, exact auxiliary-model dimensions, and allocation limits. It
does not select a solver or mutate an ambient model.

## JuMP materialization

After loading JuMP and a conic solver, materialize the atom explicitly:

```julia
using Hypatia
using JuMP

model = JuMP.Model(Hypatia.Optimizer)
x = @variable(model, [1:X.variable_count])
handle = add_top_k_p_norm_epigraph!(
    model,
    atom,
    x;
    allow_densify=true,
)

# Add constraints that define x, then use:
@objective(model, Min, handle.epigraph)
```

The handle also exposes the singular-value majorant, auxiliary thresholds,
PSD constraint references, and formulation identifier. Solver status remains
owned by the surrounding JuMP model; inspect `termination_status`,
`primal_status`, `dual_status`, objective values/bounds when supported, and
the model's residual evidence. A missing, failed, or limited solve is not a
norm value.

## Exact conic formulation

The extension forms the Hermitian dilation

```math
\mathcal{D}(X)=
\begin{bmatrix}
0&X\\
X^\dagger&0
\end{bmatrix}.
```

Its largest $j$ eigenvalues are the largest $j$ singular values of $X$ for
$j\leq r$. For each $j$, the Ky Fan epigraph uses

```math
\sum_{i=1}^{j}s_i(X)
=
\min_{\alpha,Z}
\left\{
j\alpha+\mathop{\mathrm{tr}}Z:
Z\succeq0,\quad
Z-\mathcal{D}(X)+\alpha I\succeq0
\right\}.
```

Ordered nonnegative majorants are then constrained by an MOI `NormCone`.
The cases $p=1$ and $p=\infty$ use direct linear epigraphs. This is a convex
epigraph, not a randomized estimate.

One-row and one-column affine inputs are first embedded as diagonal matrices,
so their singular values are exactly the entry magnitudes.

## Limits and type behavior

The core checks the total scalar variables, PSD block count and dimension,
stored sparse entries, and dense real-block triangle entries before the
extension allocates anything. `allow_densify=true` is required at
materialization because each complex Hermitian PSD constraint is translated
to a dense real symmetric block.

Core affine constants and coefficients preserve their real component type and
sparse storage. MathOptInterface's general `NormCone` records $p$ as a
`Float64`; this affects the optional cone parameter, not numeric
[`top_k_p_norm`](@ref).

This atom is the no-loss Julia-native replacement for the CVX-expression
branch of QETLAB `kpNorm.m` at the pinned revision. It uses public JuMP/MOI
APIs and does not open a hidden nested optimization model.
