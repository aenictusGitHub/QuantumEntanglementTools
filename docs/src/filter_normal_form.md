# Bounded filter normal form

`filter_normal_form` computes a bipartite filter normal form without hiding
iteration failure, rank deficiency, conditioning, or resource exhaustion
behind a Boolean or an exception from a matrix inverse. It combines the
bounded [`operator_sinkhorn`](@ref) implementation with the reviewed
Hermitian-factor [`operator_schmidt_decomposition`](@ref).

## What a successful result certifies

For

```julia
result = filter_normal_form(rho, (dA, dB))
```

with `result.status == :converged`, define

```math
n = d_A d_B,\qquad
\tau = \mathtt{result.normal\_form\_trace},\qquad
K = F_A\otimes F_B.
```

The result has checked both identities

```math
\widetilde{\rho}
\approx K\rho K^\dagger
\approx
\frac{1}{n}\left(
  \tau I +
  \sum_k \xi_k G_{A,k}\otimes G_{B,k}
\right).
```

The corresponding fields are:

- `filtered_operator` for $\widetilde{\rho}$;
- `coefficients` for the complete thin vector $\xi$;
- `left_operators` and `right_operators` for the Hermitian,
  Frobenius-orthonormal $G_{A,k}$ and $G_{B,k}$;
- `left_filter`, `right_filter`, and `local_filters` for $F_A$ and
  $F_B$;
- `normal_form_reconstruction_residual` and `filter_identity_residual` for
  the two independently checked equations.

The formula includes $\tau$ so it remains correct for a positive operator
whose trace is not one. The implementation never silently normalizes the
input or repairs the returned trace.

Here is an executable two-qubit example.

```jldoctest filter-normal-form
julia> using QuantumEntanglementTools, LinearAlgebra

julia> bell = [1.0, 0, 0, 1] / sqrt(2);

julia> rho = 0.8 * (bell * bell') + 0.2 * Matrix{Float64}(I, 4, 4) / 4;

julia> result = filter_normal_form(rho, (2, 2));

julia> (result.status, result.input_numerical_rank,
        result.coefficient_numerical_rank)
(:converged, 4, 3)

julia> result.normal_form_reconstruction_residual <=
       result.normal_form_reconstruction_tolerance
true

julia> result.filter_identity_residual <= result.filter_identity_tolerance
true
```

The normal-form reconstruction can also be checked directly.

```jldoctest filter-normal-form
julia> reconstructed =
           (result.normal_form_trace * Matrix{Float64}(I, 4, 4) +
            tensor_sum(result.left_operators, result.right_operators;
                       weights=result.coefficients)) / 4;

julia> isapprox(reconstructed, result.filtered_operator; atol=2e-15)
true
```

## Status is not a bare success flag

The native routine always returns a `FilterNormalFormResult` after a valid
input has entered the bounded algorithm:

| `status` | Meaning |
|:--|:--|
| `:converged` | Sinkhorn balancing and both reconstruction checks passed |
| `:max_iterations` | The requested number of complete Sinkhorn sweeps was exhausted |
| `:singular_marginal` | A local inverse square root does not exist |
| `:ill_conditioned` | A cumulative local filter reached the condition-number limit |
| `:work_limit` | The next validation, sweep, or decomposition stage exceeds `max_work` |
| `:numerical_failure` | A checked spectral, Hermiticity, or reconstruction boundary failed |

Only `:converged` sets `converged=true`. When Sinkhorn started,
`sinkhorn_result` retains its residual history, local marginal eigenvalues,
condition numbers, failed subsystem, and work accounting.

For example, the pure product state $|00\rangle$ has singular one-party
marginals.

```jldoctest
julia> using QuantumEntanglementTools

julia> product_state = zeros(4, 4); product_state[1, 1] = 1;

julia> failure = filter_normal_form(product_state, (2, 2));

julia> (failure.status, failure.failed_subsystem, failure.converged)
(:singular_marginal, 1, false)
```

This is an algorithmic failure to construct the requested form. It is not a
separability or entanglement conclusion.

## Rank and coefficient diagnostics

`input_numerical_rank` counts input eigenvalues strictly above

```math
\mathtt{rank\_atol}
+ \mathtt{rank\_rtol}\lambda_{\max}(\rho).
```

The default relative tolerance is `size(rho, 1) * eps(R)`.
`input_minimum_eigenvalue`, `input_maximum_eigenvalue`, and
`input_numerically_full_rank` expose the underlying decision. A negative
input eigenvalue is rejected; it is never clipped.

Low rank does not automatically mean failure. A maximally entangled pure
state has rank one but full-rank balanced marginals:

```jldoctest
julia> using QuantumEntanglementTools

julia> bell = ComplexF64[1, 0, 0, 1] / sqrt(2);

julia> result = filter_normal_form(bell * bell', (2, 2));

julia> (result.status, result.input_numerical_rank,
        result.input_numerically_full_rank)
(:converged, 1, false)
```

The native `coefficients` field retains every thin operator-Schmidt term,
including numerical zeros. `coefficient_numerical_rank` instead counts terms
above

```math
\mathtt{coefficient\_atol}
+ \mathtt{coefficient\_rtol}\max_k \xi_k.
```

Its default relative tolerance is
`max(dA^2, dB^2) * eps(R)`, the reviewed scale of QETLAB's default
operator-Schmidt truncation. Keeping the full vector makes truncation
visible rather than silently changing the decomposition.

## Numerical and resource policy

The dependency-free path accepts `Float32`, `Float64`, `ComplexF32`, and
`ComplexF64` and preserves that precision. Integer and `BigFloat` matrices
are rejected rather than converted. Sparse matrices require the explicit
opt-in `allow_densify=true`; the result is dense.

The routine requires an exactly Hermitian, finite, positive-semidefinite
input with positive trace. It does not:

- normalize or symmetrize the input;
- clip a negative eigenvalue;
- use `inv`, `sqrtm`, or global warning state;
- repair the returned filtered operator or its trace;
- discard operator-Schmidt terms in the native result.

Floating-point congruences can produce a tiny skew-Hermitian residual. The
returned `centered_operator` remains untouched. Only a derived
decomposition work matrix is averaged with its adjoint, and only after
`decomposition_hermiticity_residual` is below the reported
`decomposition_hermiticity_tolerance`.

The important bounds are:

- `max_iterations` for complete Sinkhorn sweeps;
- `max_condition_number` for cumulative local filters;
- `max_entries` for a conservative count of explicitly materialized dense
  operator and Hermitian-basis entries;
- `max_work` for input rank analysis, Sinkhorn work, decomposition, and both
  reconstruction checks.

`work_used`, `estimated_postprocessing_work`, and
`required_next_stage_work` explain where a `:work_limit` result stopped. Pass
`nothing` for an allocation or work bound only when unbounded resource use
is intentional.

## QETLAB compatibility

```julia
xi, GA, GB, FA, FB =
    MATLABCompat.FilterNormalForm(rho, (dA, dB), tol)
```

The compatibility wrapper:

- accepts the pinned default, scalar, singleton, and two-element `DIM`
  forms;
- returns exactly five positionally destructurable fields
  `(xi, GA, GB, FA, FB)`;
- returns only terms above the reviewed QETLAB operator-Schmidt threshold;
- raises a `DomainError` containing the native result unless the complete
  checked pipeline converged.

The pinned `FilterNormalForm.m` parses `TOL` but does not pass it to
`OperatorSinkhorn`, so the supplied value has no effect. The wrapper
intentionally corrects that omission: `TOL` is the absolute aggregate
Sinkhorn residual threshold and sets the reviewed conditioning limit
`1/TOL`. A nonpositive tolerance is rejected because the replacement
iteration is always bounded.

The pinned routine also turns every `OperatorSinkhorn:LowRank` error into
`FilterNormalForm:NoFNF`. The Julia result preserves the distinct singular,
ill-conditioned, iteration, work, and numerical statuses instead. The
native trace-aware identity above also corrects the pinned prose formula,
which implicitly assumes a trace-one density matrix.
