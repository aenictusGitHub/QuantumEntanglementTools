```@meta
CurrentModule = QuantumEntanglementTools
DocTestSetup = quote
    using QuantumEntanglementTools
    using LinearAlgebra
    using Random
end
```

# Separability and local discrimination

The APIs on this page preserve three possible conclusions:

- `:separable` means that a sufficient certificate was checked;
- `:entangled` means that a witness or necessary-condition violation was
  checked;
- `:unknown` means that the selected methods did not decide the question.

In particular, passing PPT or an outer symmetric-extension test does not prove
separability in general dimensions. Conversely, failing an inner
approximation does not prove entanglement. `is_separable` returns an
`EntanglementReport`, not a `Bool`, so these distinctions cannot be silently
discarded.

## Explicit separable and entangled examples

A density matrix diagonal in the computational product basis has its
decomposition written directly in its diagonal entries. For example,

```math
\rho_{\mathrm{sep}}
= 0.4\lvert 00\rangle\langle 00\rvert
+ 0.1\lvert 01\rangle\langle 01\rvert
+ 0.2\lvert 10\rangle\langle 10\rvert
+ 0.3\lvert 11\rangle\langle 11\rvert .
```

```@example is-separable-diagonal
using LinearAlgebra
using QuantumEntanglementTools

rho_sep = Diagonal([0.4, 0.1, 0.2, 0.3])
separable = is_separable(rho_sep, (2, 2))

@assert separable.status === :separable
@assert separable.certified
@assert separable.certificate_kind ===
    :computational_basis_product_decomposition

(
    status=separable.status,
    certificate=separable.certificate_kind,
    weights=[component.weight for component in separable.evidence.components],
)
```

The Bell state

```math
\lvert\Phi^+\rangle
= \frac{\lvert 00\rangle+\lvert 11\rangle}{\sqrt{2}}
```

has a negative partial transpose. That violation supplies a decomposable
entanglement witness:

```@example is-separable-bell
using QuantumEntanglementTools

# Write the exactly represented projector directly. Constructing it through
# `1 / sqrt(2)` can leave a nonzero floating trace residual.
rho_bell = ComplexF64[
    0.5 0 0 0.5
    0   0 0 0
    0   0 0 0
    0.5 0 0 0.5
]
entangled = is_separable(rho_bell, (2, 2); strategies=(:ppt,))

@assert entangled.status === :entangled
@assert entangled.certified
@assert entangled.certificate_kind === :negative_partial_transpose_witness

(
    status=entangled.status,
    method=entangled.method,
    certificate=entangled.certificate_kind,
)
```

Now consider the two-qutrit isotropic state

```math
\rho(q)
= (1-q)\frac{I_9}{9}
+q\lvert\Phi_3\rangle\langle\Phi_3\rvert,
\qquad
\lvert\Phi_3\rangle
= \frac{1}{\sqrt{3}}\sum_{j=0}^{2}\lvert jj\rangle .
```

At `q = 0.1`, PPT alone is inconclusive:

```@example is-separable-unknown
using LinearAlgebra
using QuantumEntanglementTools

phi3 = vec(Matrix{Float64}(I, 3, 3)) / sqrt(3)
rho3 = 0.9 * Matrix{Float64}(I, 9, 9) / 9 + 0.1 * (phi3 * phi3')
unknown = is_separable(rho3, (3, 3); strategies=(:ppt,))

@assert unknown.status === :unknown
@assert !unknown.certified

(
    status=unknown.status,
    attempted_method=only(unknown.attempts).method,
    attempted_status=only(unknown.attempts).status,
)
```

`unknown` is the correct conclusion here: this call established a necessary
condition and no sufficient certificate.

## Selecting methods

The default `strategies=:qetlab_deterministic` runs a bounded,
certificate-first sequence inspired by the pinned QETLAB `IsSeparable`
routine. Individual methods can be selected explicitly:

```julia
report = is_separable(
    rho,
    dims;
    strategies=(
        :ppt,
        :realignment,
        :centered_realignment,
        :reduction,
        :positive_maps,
        :separable_ball,
    ),
)
```

The supported deterministic strategies include PPT and low-rank PPT theorems,
realignment and centered realignment, reduction, low-dimensional
qubit--qudit criteria, the $3\times3$ rank-four Chow test, the
Gurvits--Barnum ball, positive-map tests, filter covariance, and outer and
inner symmetric extensions. Each `EntanglementAttempt` retains the method,
status, certification flag, raw evidence, backend kind, and explanatory
message. Returned attempt histories use immutable tuple storage. Numeric arrays
inside package-owned decisive evidence are copied into read-only views, so
mutating a returned object cannot silently invalidate its certificate.

The product-subtraction strategy is heuristic, so it requires an explicit
random-number generator and explicit finite budgets:

```julia
using Random

rng = Xoshiro(2026)
report = is_separable(
    rng,
    rho,
    dims;
    strategies=(:randomized_subtraction,),
    max_subtractions=12,
    max_product_restarts=8,
    max_product_iterations=32,
    max_work=10_000_000,
)
```

The method never uses Julia's global random stream. Exhausting a search or
work budget produces `:unknown`, not a negative mathematical result.

## Symmetric-extension certificates

Outer and inner hierarchies need an explicit optional optimization backend:

```julia
using Hypatia
using QuantumEntanglementTools

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
    atol=5e-7,
    rtol=5e-7,
)

outer = is_separable(
    rho,
    dims;
    strategies=(:symmetric_extension,),
    extension_orders=(2, 3),
    extension_ppt=true,
    backend=backend,
)

inner = is_separable(
    rho,
    dims;
    strategies=(:symmetric_inner_extension,),
    extension_orders=(2,),
    extension_ppt=true,
    backend=backend,
)
```

Interpret these results in their proper direction:

- a validated outer-hierarchy separator certifies entanglement;
- a feasible outer extension is only a necessary-condition pass and remains
  `:unknown`;
- validated membership in the inner cone certifies separability;
- nonmembership in the inner cone remains `:unknown`.

The raw hierarchy result, primal or dual data, residuals, optimizer identity,
and options are retained in the corresponding attempt.

## Local state discrimination

`local_distinguishability` builds a separable-measurement **outer**
relaxation. For states $\rho_j$, priors $p_j$, and effects $M_j$, the objective
is

```math
\sum_j p_j\,\mathrm{tr}(M_j\rho_j).
```

Every effect is represented as the first-copy marginal of a positive,
permutation-invariant extension. Optional PPT and bosonic constraints tighten
or compress the model. The package also retains the always-valid separable
strategy that guesses the most probable state.

The state traces (or pure-column squared norms) and the represented prior sum
must be exactly one. A value accepted only through `atol`/`rtol` is rejected
before the analytic single-guess branch or the universal probability-one bound.
In particular, a prior vector such as `[1.0, 1e-10]` is not a deterministic
ensemble, and a deterministic prior requires every other represented prior to
be exactly zero.

This model can be inspected without installing a solver:

```@example local-discrimination-model
using QuantumEntanglementTools

states = hcat(
    ComplexF64[1, 0, 0, 0],
    ComplexF64[0, 0, 0, 1],
)
problem = local_distinguishability_problem(
    states,
    (2, 2);
    order=2,
    ppt=true,
    bosonic=true,
)

(
    variables=problem.program.variable_count,
    state_count=length(problem.states),
    dimensions=problem.dimensions,
)
```

With an optional backend:

```julia
result = local_distinguishability(
    states,
    (2, 2);
    order=2,
    ppt=true,
    bosonic=true,
    backend=backend,
)

result.separable_lower_bound
result.separable_upper_bound
result.relaxation_value
result.measurement
result.residuals
```

Read the fields separately:

- `separable_lower_bound` is achieved by an actual separable POVM;
- `separable_upper_bound` is the universal probability-one bound unless a
  validated dual bound for the configured outer relaxation supplies a tighter
  value;
- `relaxation_value` is populated only when checked primal and dual bounds
  agree;
- `measurement` is feasible for the outer hierarchy, but is not thereby
  certified separable or LOCC.
- `dual_solution` retains the analytic or solver completeness dual and raw
  dual data; it is not by itself an entanglement witness.

Thus the relaxation can bound separable discrimination without pretending
that its returned primal effects form a separable measurement.

## Perfect separable discrimination of a UPB

For a minimal unextendible product basis, `upb_sep_distinguishable` enumerates
the finitely many isolated replacement product vectors $q_r$ and asks whether

```math
I = \sum_r \lambda_r q_r q_r^\dagger,
\qquad \lambda_r \geq 0 .
```

The two-qutrit Tiles UPB has 30 replacement candidates:

```@example upb-separable-discrimination-model
using QuantumEntanglementTools

tiles = upb(:tiles)
result = upb_sep_distinguishable(tiles.local_factors)

@assert result.status === :backend_unavailable
@assert result.separably_distinguishable === nothing
@assert result.candidates_generated == 30

(
    status=result.status,
    candidates=result.candidates_generated,
    partitions=result.partitions_examined,
)
```

Supplying `backend=backend` solves the package-owned nonnegative conic
feasibility problem. A small floating-point reconstruction residual is
reported as `feasibility == :numerically_feasible` while
`separably_distinguishable` remains `nothing`. The current replacement
generator uses floating SVDs, so even a bitwise reconstruction of its computed
projectors is not by itself a proof about the exact input UPB. Likewise, a
Farkas functional can validate separation from the computed floating cone but
is not promoted to `false` without a rigorous margin covering replacement
generation error. Because both the current public input domain and replacement
generator are floating-point only, there is no reachable exact linked
positive-certificate path either. The Boolean field therefore always remains
`nothing`; compatibility mode with `structured=false` always raises
`DomainError`.

The input must first pass the full UPB analysis: product structure, mutual
orthogonality, incompleteness, and unextendibility. Complex overlaps use the
conjugating Hilbert inner product.

## Input validation and resource limits

Density matrices must already be finite, exactly represented with unit trace,
Hermitian, positive semidefinite, and dimensionally consistent. The
separability and local-discrimination certificate APIs deliberately reject a
nonzero normalization residual even when it lies inside `atol`/`rtol`. State
vectors and UPB factors are also checked and are never normalized implicitly.
No API symmetrizes, clips, projects, or repairs its input.

Dense spectral work and hierarchy growth are guarded by
`max_dense_entries`, `max_work`, `max_filter_work`, `max_extension_order`, and
`OptimizationLimits`. Randomized subtraction also has separate restart,
iteration, and subtraction limits. UPB enumeration additionally checks
`max_partitions`, `max_candidates`, `max_states`, and its own scalar-work
bound. Its work counter charges each recursive partition-search node before
visiting it, then preflights the estimated SVD, candidate-construction, and
duplicate-comparison work before evaluating a complete partition. Exhausting
any of these valid resource caps returns `status == :resource_limit` with no
Boolean conclusion; malformed cap values still raise `ArgumentError`.

## References and pinned-source differences

The formulations and certificate directions were checked against:

- [Horodecki, Horodecki, and Horodecki, *Separability of mixed states:
  necessary and sufficient conditions*](https://arxiv.org/abs/quant-ph/9605038);
- [Doherty, Parrilo, and Spedalieri, *A complete family of separability
  criteria*](https://arxiv.org/abs/quant-ph/0308032);
- [Cosentino, *Positive-partial-transpose-indistinguishable
  states via semidefinite programming*](https://arxiv.org/abs/1205.1031);
- [Bandyopadhyay et al., *Limitations on separable measurements by
  convex optimization*](https://arxiv.org/abs/1408.6981);
- [Horodecki, Lewenstein, Vidal, and Cirac, *Operational criterion and
  constructive checks for the separability of low-rank density
  matrices*](https://arxiv.org/abs/quant-ph/0002089).

These native APIs preserve the broad executable intent of pinned QETLAB
`LocalDistinguishability.m`, `IsSeparable.m`, and
`UPBSepDistinguishable.m`, with deliberate corrections:

- state inputs are rejected rather than silently normalized;
- randomized subtraction uses only a caller-supplied RNG and fixes the pinned
  routine's stale-update step;
- floating eigensolver equality is never promoted to the input-level exact
  rank-one identity-perturbation hypothesis;
- complex UPB replacement vectors use conjugating transposes;
- CVX or linear-program termination, floating-cone Farkas separation, and
  small residuals remain structured numerical evidence instead of unchecked
  Boolean answers.
