# Solver-free polynomial foundations

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

This page covers the numeric, dependency-free foundations corresponding to
QETLAB's `CopositivePolynomial`, `PolynomialAsMatrix`, and
`PolynomialOptimize` routines. It does not add a symbolic polynomial parser,
CVX expressions, or a semidefinite-programming backend.

## Coefficients and their order

`HomogeneousPolynomial` stores the coefficients of a homogeneous polynomial
without changing their numeric type. Sparse input stays sparse. The fixed
ordering is `:qetlab_lexicographic`.

For three variables and degree two:

```jldoctest polynomial-foundations
julia> monomial_exponents(3, 2)
6×3 Matrix{Int64}:
 2  0  0
 1  1  0
 1  0  1
 0  2  0
 0  1  1
 0  0  2
```

Thus `[a,b,c,d,e,f]` represents

```math
a x_1^2 + b x_1x_2 + c x_1x_3
+ d x_2^2 + e x_2x_3 + f x_3^2.
```

There are $\binom{n+m-1}{m}$ coefficients for an $n$-variable polynomial of
degree $m$. The constructor checks this count before copying the input.

```jldoctest polynomial-foundations
julia> using LinearAlgebra, SparseArrays

julia> p = HomogeneousPolynomial(sparsevec([1, 6], [2, 3], 6), 3, 2)
HomogeneousPolynomial(3 variables, degree 2, sparse Int64 coefficients)

julia> evaluate_polynomial(p, [1//2, 1//3, 1//4])
11//16
```

The coefficient vector is owned by the polynomial. Julia arrays are mutable,
so every consuming routine also revalidates the stored length and finiteness.

## A copositivity quartic

For a real symmetric matrix $C$, `copositive_polynomial(C)` constructs

```math
p(x) = y^\mathsf{T} C y,\qquad y_i=x_i^2.
```

The matrix is copositive exactly when this quartic is nonnegative for every
real $x$.

```jldoctest polynomial-foundations
julia> C = [1 2; 2 3];

julia> p = copositive_polynomial(C)
HomogeneousPolynomial(2 variables, degree 4, sparse Int64 coefficients)

julia> collect(p.coefficients)
5-element Vector{Int64}:
 1
 0
 4
 0
 3

julia> x = [2//3, -3//5];

julia> evaluate_polynomial(p, x) == dot(x .^ 2, C * (x .^ 2))
true
```

The pinned MATLAB routine silently applies `(C+C')/2`. The Julia API does not:
it rejects a nonsymmetric matrix. If taking a symmetric part is intended, make
that modeling decision explicitly before calling the function.

## Compact symmetric matrix representation

Let $p$ have degree $2d$, let $K$ be the hierarchy level, and put $k=d+K$.
For every weak composition $\alpha$ of $k$, define the normalized symmetric
monomial

```math
v_\alpha(x)=
\sqrt{\frac{k!}{\alpha_1!\cdots\alpha_n!}}\,
x_1^{\alpha_1}\cdots x_n^{\alpha_n}.
```

`polynomial_as_matrix(p; level=K)` returns the fully symmetric matrix $M$ in
the same ordering as `monomial_exponents(n,k)`, satisfying

```math
v(x)^\mathsf{T} M v(x)
=
\left(x_1^2+\cdots+x_n^2\right)^K p(x).
```

For a quadratic, this recovers the familiar symmetric coefficient matrix:

```jldoctest polynomial-foundations
julia> p = HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2);

julia> Matrix(polynomial_as_matrix(p))
2×2 Matrix{Float64}:
 1.0  1.0
 1.0  3.0
```

Sparse output is the default. Pass `sparse_output=false` to request
densification explicitly; `max_dense_entries` is checked before that dense
matrix is allocated. Integer and rational coefficients remain exact in the
polynomial object, but the square-root normalization of the compact matrix
generally requires floating entries. `Float32` and `Float64` inputs retain
their precision.

## Bounds on the unit sphere

`polynomial_bounds` applies the solver-free generalized-eigenvalue hierarchy
to a real, even-degree polynomial on
$x_1^2+\cdots+x_n^2=1$. It requires an explicit random-number generator even
when no inner samples are requested.

```julia
using Random

p = HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
result = polynomial_bounds(
    MersenneTwister(2026),
    p;
    level=1,
    sense=:max,
    inner_samples=500,
    allow_densify=true,
)

@assert result.outer_bound >= 2 + sqrt(2)
@assert result.inner_bound <= 2 + sqrt(2)
@assert result.samples_evaluated == 500
@assert isapprox(norm(result.best_point), 1)
```

The two values have deliberately different meanings:

| Field | Meaning for `sense=:max` | Meaning for `sense=:min` |
|---|---|---|
| `outer_bound` | hierarchy upper bound | hierarchy lower bound |
| `inner_bound` | best sampled feasible lower value | best sampled feasible upper value |
| `outer_kind` | `:hierarchy_outer_bound` | `:hierarchy_outer_bound` |
| `inner_kind` | `:sampled_feasible` when requested | `:sampled_feasible` when requested |

`outer_uncertainty` records the numerical reconstruction allowance used to pad
the computed generalized eigenvalue. A sampled inner value is attained at
`best_point`, but sampling does not certify global optimality.

`inner_samples` is an exact count. It is never derived from a wall clock, and
only the supplied `rng` is consumed. If the outer bound already establishes a
requested `target`, `inner_kind` is `:skipped_by_outer_target` and no random
draw occurs.

Maximization is evaluated directly. This corrects a pinned recursion defect:
the MATLAB routine negates the polynomial to call minimization but forwards
the target without changing its sign. The Julia result reports the target
conclusion explicitly as `:outer_proves_at_most`,
`:outer_proves_at_least`, or `:not_proven`.

## Resource and solver boundaries

The following guards are checked before the corresponding expensive work:

- `max_terms` bounds coefficient and exponent tables;
- `max_exponent_entries` bounds each dense monomial-exponent table, including
  the variable axis that a term count alone does not control;
- `max_degree` bounds factorial scaling work;
- `max_dimension` bounds the compact matrix dimension;
- `max_dense_entries` bounds each explicitly requested or internally required
  dense matrix before conversion;
- `max_nonzeros` bounds sparse matrix construction;
- `max_work` bounds deterministic construction, eigensolver, and sampling
  estimates;
- `max_samples` bounds the requested exact sample count.

`polynomial_bounds` uses the standard-library dense symmetric eigensolver and
therefore requires `allow_densify=true`. `BigFloat`, symbolic expressions, CVX
objects, and solver-backed SOS models are outside this solver-free method. They
must be supplied by a separately tested optional integration rather than
silently converted.

## API

```@docs
HomogeneousPolynomial
PolynomialOptimizationResult
monomial_exponents
evaluate_polynomial
copositive_polynomial
polynomial_as_matrix
polynomial_bounds
```
