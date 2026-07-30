# Pure-state robustness of k-coherence

`pure_k_coherence_robustness(state, k)` evaluates a solver-free closed formula
for a normalized pure state. It returns the common standard and generalized
robustness of $k$-coherence proved in
[Johnston *et al.*, Phys. Rev. A 98, 022328 (2018)](https://doi.org/10.1103/PhysRevA.98.022328).

This API is deliberately stricter than the pinned `RobkCohValue` routine:

- it accepts real or complex floating-point coefficients and sorts their
  magnitudes in descending order;
- it checks that the vector is finite, nonempty, one-based, and normalized
  without changing it;
- it requires $2 \leq k \leq n$ as in Theorem 1; and
- it returns structured branch diagnostics rather than hiding a
  tolerance-sensitive comparison.

## Formula

Write the sorted magnitudes as
$a_1 \geq a_2 \geq \cdots \geq a_n \geq 0$ and define
$s_j = \sum_{i=j}^{n} a_i$. The branch index $\ell$ is the largest integer in
$\{2,\ldots,k\}$ satisfying

```math
a_{\ell-1} \geq \frac{s_\ell}{k-\ell+1}.
```

If no such index exists, $\ell=1$. The returned value is

```math
R_k =
\frac{s_\ell^2}{k-\ell+1}
- \sum_{i=\ell}^{n} a_i^2.
```

The comparison is non-strict. An exact equality therefore selects that branch,
and the descending search ensures that the largest admissible index is used.

## Reading the result

```julia
psi = ComplexF64[1, im, 1, -im] / 2
result = pure_k_coherence_robustness(psi, 2)

result.value
result.branch_index
result.branch_status
```

`branch_status` is `:stable`, `:exact_equality`, or `:near_boundary`.
`selected_gap` records how far the selected comparison succeeds, while
`next_gap` records how far the next larger branch fails. A `nothing` gap means
that the corresponding adjacent branch does not exist.

The `atol` and `rtol` keywords serve two purposes: they validate the input
normalization and set the reported branch-stability threshold. They never
normalize the vector and never change which branch the exact input-arithmetic
comparison selects.

## Corrected compatibility behavior

The pinned MATLAB loop assumes that its input already consists of sorted,
nonnegative, normalized coefficients. Applying it directly to an unsorted or
complex state can select the wrong branch or even return a complex quantity.
The Julia compatibility wrapper delegates to the validated native operation,
so it intentionally corrects those unsafe cases while retaining the two pinned
outputs `(robustness, branch_index)`.

Sorting dominates the cost: $O(n\log n)$ time and $O(n)$ workspace. The
calculation retains `Float32`, `Float64`, or `BigFloat` arithmetic from the
state. Integer and rational coefficient vectors are rejected rather than
silently converted.
