# Copositivity and clique-number bounds

This page covers the native replacements for QETLAB's `IsCopositive` and
`CliqueNumber` entry points. Both APIs expose the evidence they actually have:
an exact witness or graph certificate can decide a question, while a floating
SOS value remains numerical evidence.

The implementation follows the pinned QETLAB polynomial construction, but it
does not reproduce three unsafe behaviors:

- no fixed decimal offset turns a numerical boundary into `true` or changes an
  integer bound;
- no wall-clock loop decides how many random points to draw; and
- an unavailable or failed optimizer is not reported as a negative result.

The mathematical background is the Motzkin--Straus theorem and copositive-cone
hierarchies. See the
[original Motzkin--Straus paper](https://doi.org/10.4153/CJM-1965-053-6)
and the
[de Klerk--Pasechnik copositive formulation](https://doi.org/10.1137/S1052623401383248).

## Copositivity: three possible conclusions

A real symmetric matrix $C$ is copositive when

```math
y^\mathsf{T} C y \geq 0
\qquad\text{for every }y\geq 0.
```

By homogeneity, it is enough to consider the nonnegative simplex. With
$y_i=x_i^2$, the package constructs the quartic

```math
p_C(x)=\sum_{i,j} C_{ij}x_i^2x_j^2
```

and bounds its minimum on the real unit sphere.

`copositivity_criterion` can return:

- `verdict == true` with an entrywise-nonnegative or exact-PSD sufficient
  certificate;
- `verdict == false` with an exact rational nonnegative vector whose quadratic
  value is negative outside the tolerance band; or
- `verdict == nothing` when only hierarchy evidence, a boundary value, a
  resource limit, or a backend failure is available.

Here is an exact positive example with negative off-diagonal entries:

```@example copositive_clique
using QuantumEntanglementTools
using LinearAlgebra
using Random

C = Rational{Int}[1 -1; -1 1]
positive = copositivity_criterion(MersenneTwister(1), C)

@assert positive.verdict === true
@assert positive.certificate_kind === :exact_positive_semidefinite
(positive.status, positive.lower_bound, positive.upper_bound)
```

The next matrix has a direct two-coordinate counterexample:

```@example copositive_clique
Cbad = [1 -2; -2 1]
negative = copositivity_criterion(MersenneTwister(2), Cbad)

@assert negative.verdict === false
@assert negative.witness.simplex_vector == [1 // 2, 1 // 2]
@assert negative.witness.value == -1 // 2
(negative.status, negative.certificate_kind)
```

The witness uses exact `Rational{BigInt}` arithmetic, even when the original
matrix used a binary floating-point type. Its simplex coordinates are an
owned, read-only array: copy them if you need mutable working storage.
Exact coordinate, pair, entrywise-nonnegative, and exact-PSD branches run
before polynomial construction, so `negative.polynomial === nothing` here.

### A boundary is not a proof

The pinned routine declares success when its numerical lower bound is at least
`-1e-9`. That can label a small negative diagonal as copositive. The native
contract instead retains the boundary:

```@example copositive_clique
Cboundary = [-5.0e-10 0.0; 0.0 1.0]
boundary = copositivity_criterion(MersenneTwister(3), Cboundary)

@assert boundary.verdict === nothing
@assert boundary.status ===
    QuantumEntanglementTools.CopositivityNumericalBoundary
(boundary.tolerance, boundary.message)
```

Set `atol` and `rtol` explicitly when the data have a known error model.
Integer and rational matrices are exact inputs and therefore require both
tolerances to be zero.

## Selecting a hierarchy

The default `hierarchy=:sos` constructs the package-owned
`SemidefiniteProgram`. No solver is selected implicitly:

```@example copositive_clique
H = [
     1 -1  1  1 -1
    -1  1 -1  1  1
     1 -1  1 -1  1
     1  1 -1  1 -1
    -1  1  1 -1  1
]

without_backend = copositivity_criterion(MersenneTwister(4), H)
@assert without_backend.verdict === nothing
@assert without_backend.status ===
    QuantumEntanglementTools.CopositivityBackendUnavailable
without_backend.hierarchy_result.status
```

After loading JuMP and a supported optimizer, pass a backend explicitly:

```julia
using Hypatia
using JuMP

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)
result = copositivity_criterion(MersenneTwister(5), H; backend)
```

An optimal floating solve supplies a numerical SOS lower bound, not an exact
copositivity certificate. Inspect `result.hierarchy_result`, its primal and
dual residuals, and the returned moment matrix.

The solver-free alternative is selected with `hierarchy=:nosdp`. It performs a
guarded dense generalized-eigenvalue calculation, so densification must be
acknowledged:

```@example copositive_clique
solver_free = copositivity_criterion(
    MersenneTwister(6),
    H;
    hierarchy=:nosdp,
    allow_densify=true,
    inner_samples=16,
)

@assert solver_free.samples_evaluated == 16
(solver_free.lower_bound, solver_free.upper_bound, solver_free.verdict)
```

For minimization, the hierarchy outer value is a lower bound and the best
sample is an attained upper bound. `inner_samples` is the exact number of
draws. Passing zero requests no samples, and a branch that is already decided
does not consume the caller's RNG.

## Clique-number bounds

`clique_number_bounds` accepts only a nonempty simple undirected graph:

- the adjacency matrix is exactly symmetric;
- every diagonal entry is zero; and
- every entry is exactly zero or one.

Weighted, directed, looped, nonfinite, and silently rounded graphs are
rejected.

For an adjacency matrix $A$, the Motzkin--Straus theorem states

```math
\max_{\substack{y\geq 0\\\sum_i y_i=1}} y^\mathsf{T}Ay
=1-\frac{1}{\omega(G)},
```

where $\omega(G)$ is the clique number. The implementation preserves the
continuous hierarchy values, but its returned integer interval is certified
independently:

- a deterministic greedy clique proves the lower bound;
- edge count and maximum degree prove elementary upper bounds;
- a deterministic proper coloring proves another upper bound; and
- an attained hierarchy point may improve the lower bound only after exact
  rational normalization and evaluation.

The result retains all three component upper bounds in `upper_certificate`,
including the full proper coloring. `upper_certificate_kinds` lists the
components that attain the final minimum, not merely all computations that
were attempted. The returned clique, coloring, and exact simplex witness use
owned, read-only storage.

The five-cycle is a useful example:

```@example copositive_clique
function cycle_adjacency(n)
    A = zeros(Int, n, n)
    for vertex in 1:n
        neighbor = mod1(vertex + 1, n)
        A[vertex, neighbor] = 1
        A[neighbor, vertex] = 1
    end
    return A
end

C5 = cycle_adjacency(5)
bounds = clique_number_bounds(MersenneTwister(7), C5)

@assert bounds.lower_bound == 2
@assert bounds.upper_bound == 3
@assert bounds.bounds_certified
@assert !bounds.exact
@assert bounds.upper_certificate_kinds ==
    (:edge_count, :maximum_degree, :greedy_coloring)
(
    bounds.status,
    bounds.best_clique,
    bounds.upper_certificate.greedy_coloring_bound,
)
```

The missing SOS backend does not invalidate those graph-theoretic bounds. A
solver-free hierarchy run retains additional evidence:

```@example copositive_clique
sampled_bounds = clique_number_bounds(
    MersenneTwister(8),
    C5;
    hierarchy=:nosdp,
    allow_densify=true,
    inner_samples=16,
)

@assert sampled_bounds.samples_evaluated == 16
@assert sum(sampled_bounds.motzkin_straus_witness.simplex_vector) == 1
(
    sampled_bounds.continuous_lower_bound,
    sampled_bounds.continuous_upper_bound,
    sampled_bounds.uncertified_upper_candidate,
)
```

The field `uncertified_upper_candidate` shows what the floating continuous
outer value suggests after outward tolerance padding. It never replaces the
certified integer `upper_bound`.

Edgeless and complete graphs are handled exactly without invoking a hierarchy
or consuming the RNG. Their `polynomial` field is `nothing`, because no
polynomial model is needed:

```@example copositive_clique
empty_graph = clique_number_bounds(MersenneTwister(9), zeros(Int, 4, 4))
complete_graph = clique_number_bounds(
    MersenneTwister(10),
    ones(Int, 4, 4) - Matrix{Int}(I, 4, 4),
)

@assert (empty_graph.lower_bound, empty_graph.upper_bound) == (1, 1)
@assert (complete_graph.lower_bound, complete_graph.upper_bound) == (4, 4)
@assert empty_graph.polynomial === nothing
@assert complete_graph.polynomial === nothing
nothing
```

## Resource and storage policy

Both entry points validate matrix dimensions and logical entry counts before
cheap exact certificates and construct a polynomial only when hierarchy work
is still needed. The graph routine also guards deterministic clique and
coloring work. Polynomial term counts, tensor dimensions, dense entries,
nonzeros, model variables, PSD blocks, samples, and exact-certificate bit
sizes have separate limits.

Sparse input is retained by the quartic polynomial construction. Operations
that genuinely need dense linear algebra require `allow_densify=true`; no
matrix is silently symmetrized or converted from weighted to unweighted form.
