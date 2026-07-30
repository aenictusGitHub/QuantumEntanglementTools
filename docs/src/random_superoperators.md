# Bounded random completely positive maps

[`random_superoperator`](@ref) generates a dense, rank-controlled completely
positive map between square matrix algebras. It covers the executable
`RandomSuperoperator` branches at the pinned QETLAB revision while making the
RNG, representation, marginal constraints, numerical diagnostics, and resource
limits explicit.

## A trace-preserving channel

The native function returns a [`RandomSuperoperatorResult`](@ref), not an
unchecked matrix:

```@example random-superoperators
using QuantumEntanglementTools
using LinearAlgebra
using Random

channel = random_superoperator(
    Xoshiro(0x5153),
    (2, 3);
    trace_preserving=true,
    kraus_rank=2,
    representation=:choi,
)

@assert channel.status == :success
@assert channel.complete_positivity_guaranteed
@assert channel.trace_preservation_guaranteed
@assert is_trace_preserving(channel.representation)

(
    size=size(choi_matrix(channel.representation)),
    numerical_rank=channel.numerical_kraus_rank,
    tp_residual=channel.trace_preservation_residual,
)
```

`dim=(d_in,d_out)` means a map from `d_in × d_in` matrices to
`d_out × d_out` matrices. A scalar uses the same input and output dimension.
The Choi matrix follows the package's unnormalized, input-first convention, so
its trace is `d_in` for a trace-preserving map.

The generator starts from `kraus_rank` Ginibre factors. It therefore preserves
complete positivity and the rank upper bound algebraically; it does not form a
matrix and clip negative eigenvalues. A generic draw has exactly the requested
rank with probability one. Because a concrete floating-point draw can be
ill-conditioned, the result separately records:

- `factor_singular_values`, `rank_threshold`, and
  `numerical_kraus_rank`;
- the TP, unital, and proportional-unital residuals and their tolerances;
- the balancing iteration count, residual history, minimum encountered
  marginal eigenvalue, and largest cumulative filter condition number;
- attempts, conservative work, failure location, status, and message.

The construction guarantee and the numerical diagnostic are deliberately not
the same claim.

## Unital and bistochastic maps

Set `unital=true` to require $\Phi(I_{\mathrm{in}})=I_{\mathrm{out}}$.
Trace preservation and unitality can both hold only when the two dimensions
agree, because their Choi marginal traces must agree.

```@example random-superoperators
bistochastic = random_superoperator(
    Xoshiro(7),
    3;
    trace_preserving=true,
    unital=true,
    kraus_rank=1,
)

@assert bistochastic.status == :success
@assert bistochastic.trace_preservation_guaranteed
@assert bistochastic.unitality_guaranteed
@assert is_trace_preserving(bistochastic.representation)
@assert is_unital(bistochastic.representation)

(bistochastic.status, bistochastic.balancing_iterations)
```

Simultaneous constraints use a factor-space Sinkhorn iteration. Each sweep
alternately applies positive-definite inverse square roots to the input and
output Kraus marginals. This is algebraically equivalent to local filtering of
the Choi matrix, but avoids constructing a dense Choi matrix when
`representation=:kraus`.

`max_iterations` bounds every attempt and `max_attempts` bounds redrawing.
Singular or ill-conditioned marginals are not projected or regularized.

## The unequal-dimensional correction

The pinned routine accepts `TP=1, UN=1` with unequal dimensions, warns that
unitality is impossible, and then produces the proportional condition

```math
\Phi(I_{\mathrm{in}})
=
\frac{d_{\mathrm{in}}}{d_{\mathrm{out}}}I_{\mathrm{out}}.
```

Calling that map unital would be mathematically false. The native API therefore
rejects `trace_preserving=true, unital=true` for unequal dimensions. The
executable construction remains available under an explicit corrected name:

```@example random-superoperators
proportional = random_superoperator(
    Xoshiro(8),
    (2, 3);
    trace_preserving=true,
    unital=false,
    proportional_unital=true,
    kraus_rank=2,
)

@assert proportional.status == :success
@assert proportional.trace_preservation_guaranteed
@assert proportional.proportional_unitality_guaranteed
@assert !proportional.unitality_guaranteed
@assert proportional.proportional_unital_factor == 2 / 3

identity_output = apply_channel(
    Matrix{Float64}(I, 2, 2),
    proportional.kraus_representation,
)
@assert isapprox(identity_output, (2 / 3) * Matrix{Float64}(I, 3, 3))

proportional.proportional_unitality_residual
```

This is a documented upstream correction, not an omitted branch.

## Controlled failure

Only `status == :success` establishes the requested TP or unital guarantee.
Exhausted attempts and work limits are structured outcomes:

```@example random-superoperators
bounded = random_superoperator(
    Xoshiro(9),
    3;
    unital=true,
    kraus_rank=2,
    max_attempts=1,
    max_iterations=0,
)

@assert bounded.status == :max_attempts
@assert !bounded.succeeded
@assert bounded.complete_positivity_guaranteed
@assert !bounded.trace_preservation_guaranteed
@assert !bounded.unitality_guaranteed

(bounded.status, bounded.last_attempt_status)
```

The returned map is still completely positive because it retains the final
Kraus draw. It must not be used as though the unsuccessful marginal
normalization had been established.

`max_dimension` caps each local dimension. `max_entries` guards a conservative
peak-storage estimate, including a requested dense Choi or superoperator.
`max_work` bounds a conservative arithmetic estimate across draws and
balancing sweeps. Each limit may be set to `nothing` only after reviewing the
cost.

Random maps are dense almost surely, so there is no misleading sparse-output
option. The dependency-free spectral path supports `T=Float32` and
`T=Float64`; `real=false` uses the corresponding complex type. Unsupported
arbitrary precision is rejected rather than silently down-converted.

## QETLAB compatibility

The compatibility signature keeps the pinned positional order but still
requires an explicit leading RNG:

```julia
MATLABCompat.RandomSuperoperator(
    rng, DIM, TP=1, UN=0, RE=0, KR=prod(DIM);
    diagnostics=false,
    allow_proportional_unital=false,
    T=Float64,
    atol=0,
    rtol=nothing,
    max_attempts=8,
    max_iterations=1000,
    max_condition_number=nothing,
    max_dimension=4096,
    max_entries=10_000_000,
    max_work=1_000_000_000,
)
```

On success, the default return is the raw Choi matrix expected by QETLAB code.
Set `diagnostics=true` to retain the native structured result. A failed native
status raises `DomainError` on the raw-output path.

The compatibility layer accepts `Bool`, `0`, or `1` for the three flags and
rejects other values. This corrects the pinned behavior in which, for example,
`TP=2` silently selects an unconstrained branch. It also rejects
`KR > prod(DIM)`, where the documentation's almost-sure exact-rank statement is
impossible. For the unequal `TP=UN=1` branch, callers must explicitly set
`allow_proportional_unital=true`; the diagnostic result continues to report
that strict unitality was not established.

Both API layers leave Julia's global random stream untouched.

## Pinned-source evidence

The source-free fixture
`random_superoperator_octave_11_3_qetlab_d858961.json` records five seeded
QETLAB outputs spanning unconstrained, TP-only, unital-only, equal-dimensional
bistochastic, real, complex, rank-controlled, and unequal proportional
branches. Its SHA-256 is
`2bf98f8cbbc188be160755d8e866634658de109a0a92b402a353287c1ad9ca7a`.

The comparator validates positivity, marginal residuals, traces, rank bounds,
real output, the unequal-dimensional defect, permissive flag handling, and
the oversized-rank discrepancy. It then checks the corresponding native
and compatibility properties with independent Julia RNG draws. The focused
suite passes 139 assertions and the comparator passes 106 assertions on Julia
1.12.6 and 1.10.11. Random matrices are not compared entrywise across
unrelated RNG algorithms. The committed artifact was
generated with Octave 11.3.0 and is supplemental, function-specific evidence;
MATLAB was not available.
