# Bounded operator Sinkhorn scaling

`operator_sinkhorn` computes local invertible filters that make every
one-party marginal proportional to the identity. It is a bounded numerical
normalization routine, not a separability test and not a guarantee that every
positive-semidefinite input admits a nonsingular scaling.

## Result and reconstruction contract

For subsystem dimensions $(d_1,\ldots,d_m)$, the result stores dense local
filters $F_1,\ldots,F_m$ and a scaled operator satisfying

```math
\sigma \approx
\left(\bigotimes_{j=1}^{m} F_j\right)
\rho
\left(\bigotimes_{j=1}^{m} F_j\right)^\dagger.
```

The original positive trace is preserved. The iteration explicitly uses the
trace-one work matrix $\rho/\mathrm{tr}(\rho)$, then restores the original
trace algebraically, so this gauge choice does not silently replace the
user's input or change the returned normalization.

`OperatorSinkhornResult` contains:

- `scaled_operator` and the authoritative multipartite `local_filters`;
- `left_filter` and `right_filter`, aliases for the first two local filters;
- `iterations`, `status`, `converged`, and `message`;
- the aggregate `residual_history` and individual
  `final_marginal_residuals`;
- final local-filter condition numbers, their checked maximum, minimum
  encountered marginal eigenvalues, and the largest checked marginal
  Hermiticity residual;
- the convergence threshold, condition limit, original trace, failed
  subsystem, and conservative work count.

Only `status == :converged` is success. Other statuses are
`:max_iterations`, `:singular_marginal`, `:ill_conditioned`, `:work_limit`,
and `:numerical_failure`.

## A complete bipartite example

The product of two full-rank diagonal states is scaled to the maximally mixed
state in one complete sweep:

```@example operator-sinkhorn
using LinearAlgebra
using QuantumEntanglementTools

rho = tensor_product(
    Diagonal([0.8, 0.2]),
    Diagonal([0.6, 0.3, 0.1]),
)
result = operator_sinkhorn(rho, (2, 3))

@assert result.status == :converged
@assert result.iterations == 1

global_filter = tensor_product(result.local_filters...)
@assert isapprox(
    result.scaled_operator,
    global_filter * rho * global_filter';
    atol=2e-14,
)
@assert isapprox(
    partial_trace(result.scaled_operator, (2, 3); trace_out=(2,)),
    Matrix{Float64}(I, 2, 2) / 2;
    atol=2e-14,
)
@assert isapprox(
    partial_trace(result.scaled_operator, (2, 3); trace_out=(1,)),
    Matrix{Float64}(I, 3, 3) / 3;
    atol=2e-14,
)

(result.status, result.iterations)
```

The same API supports more than two parties. `local_filters` then contains
one matrix per entry of `dims`; `left_filter` and `right_filter` remain
convenient aliases for the first two.

## Controlled failure is data

A product projector has singular one-party marginals. The native API returns
the checked state of the iteration instead of throwing away the evidence:

```@example operator-sinkhorn
projector = zeros(4, 4)
projector[1, 1] = 1
failure = operator_sinkhorn(projector, (2, 2))

@assert failure.status == :singular_marginal
@assert failure.failed_subsystem == 1
@assert !failure.converged

(failure.status, failure.failed_subsystem)
```

Low rank alone is not automatically failure. For example, a Bell-state
projector already has identity one-party marginals and converges with zero
sweeps. Failure is tied to the encountered marginal factorization and the
explicit conditioning limit.

## Numerical and resource policy

Every complete sweep uses the Hermitian positive-definite inverse square root
of the current marginal,

```math
T_j = \frac{M_j^{-1/2}}{\sqrt{d_j}},
```

constructed from a Hermitian eigendecomposition and scalar reciprocal square
roots. No matrix inverse or unrestricted matrix square root is formed.
Marginal averaging is confined to an eigensolver work matrix only after its
Hermiticity defect passes a dimension- and scale-aware roundoff check; the
input and returned operator are never symmetrized.

The aggregate stopping residual is

```math
r = \sum_j \left\|M_j-\frac{I_{d_j}}{d_j}\right\|_{\mathrm F},
```

with threshold
$\mathtt{atol}+\mathtt{rtol}\sum_j 1/\sqrt{d_j}$.
The defaults are `atol=0` and `rtol=sqrt(eps(R))`.
`max_iterations` bounds complete sweeps. `max_condition_number` defaults to
`1/sqrt(eps(R))`. `max_entries` guards the dense full-operator work matrix,
and `max_work` bounds a conservative spectral, partial-trace, and dense
congruence estimate.

The dependency-free implementation supports `Float32`, `Float64`,
`ComplexF32`, and `ComplexF64`. It does not silently down-convert
`BigFloat`. Sparse input is rejected unless `allow_densify=true`, and the
result is dense.

The input must be finite, square, exactly Hermitian, positive semidefinite,
and have positive trace. Negative eigenvalues are rejected rather than
clipped. This strict boundary is intentional: an inconclusive spectral
boundary is not relabeled as a valid density operator.

## QETLAB compatibility

`MATLABCompat.OperatorSinkhorn` keeps the reviewed default, scalar, and
multipartite-vector `DIM` forms. `TOL` is an absolute aggregate residual
threshold and also sets the pinned `1/TOL` conditioning heuristic. On success
it returns exactly two named fields that destructure positionally:

```@example operator-sinkhorn
sigma, filters = QuantumEntanglementTools.MATLABCompat.OperatorSinkhorn(
    Matrix{Float64}(I, 4, 4) / 4,
)

@assert length(filters) == 2
@assert sigma == Matrix{Float64}(I, 4, 4) / 4
```

For every nonconverged native status, the wrapper throws a `DomainError`
containing the structured result. It never changes process-wide warning
settings. Unlike the pinned loop, it has mandatory iteration and work
bounds, does not call `inv` or unrestricted `sqrtm`, and does not apply final
symmetrization or trace repair. These are deliberate safety corrections to
the routine at pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`.
