# S(`k`) operator norms and block positivity

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
    using LinearAlgebra
    using Random
end
```

[`sk_operator_norm`](@ref) returns lower and upper bounds on a bipartite
S(`k`) operator norm. [`is_block_positive`](@ref) uses those bounds without
collapsing an inconclusive computation into a Boolean answer.

Both native functions require an explicit random-number generator. Exact
branches and target decisions that finish before a search do not consume it.

## Exact and bounded results

If $k$ reaches the smaller local dimension, the S(`k`) norm is the ordinary
operator norm:

```jldoctest
julia> X = Diagonal([4.0, 3.0, 2.0, 1.0]);

julia> result = sk_operator_norm(
           MersenneTwister(1), X; k=2, dims=(2, 2)
       );

julia> (result.status, result.lower_bound, result.upper_bound, result.exact)
(QuantumEntanglementTools.SKOperatorNormExact, 4.0, 4.0, true)

julia> result.lower_witness.validated
true
```

For a general input, inspect the interval and its provenance:

```jldoctest
julia> X = Diagonal([4.0, 3.0, 2.0, 1.0]);

julia> bounded = sk_operator_norm(
           MersenneTwister(2), X;
           k=1, dims=(2, 2), strength=0
       );

julia> (bounded.lower_bound, bounded.upper_bound)
(4.0, 4.0)

julia> (bounded.exact, bounded.lower_bound_kind, bounded.upper_bound_kind)
(false, :explicit_witness, :operator_norm)
```

The numbers happen to coincide here because a computational product vector is
optimal. The result is still not labeled analytically exact: equality of
floating-point bounds is evidence, not a proof that an exact-case theorem was
used.

`lower_witness.left` and `lower_witness.right` have Schmidt rank at most `k`
and directly verify

```math
\left|\langle u, Xv\rangle\right|.
```

For a positive-semidefinite operator the witness is quadratic, with $u=v$.
Random projected iterations can improve such lower bounds, but they remain
witness evidence rather than global-optimality certificates.

## Optional upper relaxations

[`sk_operator_norm_problem`](@ref) builds a package-owned
[`SemidefiniteProgram`](@ref):

- at level one and $k=1$, it uses the PPT relaxation;
- at level one and $k>1$, it uses
  $k\,\mathrm{Tr}_B(\rho)\otimes I_B-\rho \succeq 0$;
- at higher levels for $k=1$, it uses a bosonic-compressed PPT symmetric
  extension.

An explicit [`JuMPBackend`](@ref) activates the optional modeling extension.
The result keeps the backend termination status, objective bound, primal and
dual residuals, solver name and version, options, and model limits. A
numerically checked dual relaxation is stored in
`relaxation_upper_bound` and `upper_witness`; it is not relabeled as an exact
symbolic value.

Higher `strength` requests more work, subject to `max_hierarchy_level`,
[`OptimizationLimits`](@ref), `max_dense_entries`, and `max_work`.
The unbounded upstream `strength=-1` behavior is intentionally rejected.

## Block positivity: three outcomes

A negative product-vector expectation is a direct certificate:

```jldoctest
julia> W = Diagonal([-1.0, 2.0, 2.0, 2.0]);

julia> answer = is_block_positive(
           MersenneTwister(3), W; dims=(2, 2)
       );

julia> (answer.status, answer.verdict, answer.certified)
(QuantumEntanglementTools.BlockPositivityViolated, false, true)

julia> (answer.witness.expectation, answer.witness.schmidt_rank)
(-1.0, 1)
```

Positive definiteness is an immediate positive certificate:

```jldoctest
julia> answer = is_block_positive(
           MersenneTwister(4), Matrix{Float64}(I, 4, 4);
           dims=(2, 2)
       );

julia> (answer.status, answer.verdict)
(QuantumEntanglementTools.BlockPositivityCertified, true)
```

For an indefinite $X$, the implementation uses

```math
X\text{ is }k\text{-block positive}
\quad\Longleftrightarrow\quad
\left\|\,\|X\|I-X\,\right\|_{S(k)}\le \|X\|.
```

A quadratic lower witness above the threshold certifies a negative
expectation. A validated upper relaxation strictly below the threshold
certifies block positivity. If the bounds overlap, the optimizer stops at a
limit, or the result lies in the tolerance band, `verdict` is `nothing` and
the status is unknown, boundary, resource-limited, or backend-failure as
appropriate.

In particular, zero-eigenvalue PSD decisions are not silently rounded to
positive. Inputs must be exactly Hermitian, and no path symmetrizes,
normalizes, clips, projects, or repairs the operator. Sparse inputs require
explicit `allow_densify=true` and remain subject to the same allocation
limits.

## Reproducibility checklist

- Pass `rng` explicitly, even when an exact case is expected.
- Record `dims`, `k`, `strength`, tolerances, and all resource limits.
- Treat `lower_witness.validated` as lower-bound evidence only.
- Treat a solver relaxation as numerical certificate evidence at its recorded
  residual tolerance, not as exact arithmetic.
- Branch on `verdict === nothing`; never coerce it to `false`.

```@meta
DocTestSetup = nothing
```
