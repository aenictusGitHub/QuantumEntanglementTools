# Certifying an unextendible product basis

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

An unextendible product basis (UPB) is an **incomplete, mutually orthogonal
family of product states** whose orthogonal complement contains no product
state. `is_upb` checks all parts of that definition. It does not label a
nonorthogonal family or a complete product basis as a UPB merely because no
additional product vector is orthogonal to it.

The native API returns a [`UPBAnalysisResult`](@ref), not a bare Boolean:

```julia
result = is_upb(local_factors; normalization=:require)

result.status   # :upb, :not_upb, or :unknown
result.is_upb   # true, false, or nothing
result.reason   # the certificate, counterexample, boundary, or limit reason
```

This follows the original definition in Bennett *et al.*, *Phys. Rev. Lett.*
82, 5385 (1999), [arXiv:quant-ph/9808030](https://arxiv.org/abs/quant-ph/9808030).
The finite partition search follows the documented behavior of the pinned
QETLAB `IsUPB` routine, with the validation and complex-inner-product
corrections described below.

## Example 1: certify the two-qutrit Tiles UPB

The columns below are unnormalized local factors of the five Tiles states.
Using `normalization=:allow` is explicit: the input is not changed, and only
zero sets, spans, and orthogonality are used.

```jldoctest is_upb_witness
julia> left = [
           1  1  0  0  1
           0 -1  0  1  1
           0  0  1 -1  1
       ];

julia> right = [
           1  0  0  1  1
          -1  0  1  0  1
           0  1 -1  0  1
       ];

julia> result = is_upb(left, right; normalization=:allow);

julia> (result.status, result.reason, result.certificate_kind)
(:upb, :unextendible, :exact)

julia> result.partitions_examined
20
```

Because the entries are integers, the implementation promotes them to exact
`Rational{BigInt}` work arithmetic. Exhausting all 20 admissible ordered
partitions is therefore an exact certificate for these supplied states.

With normalized floating-point columns, the same call returns
`certificate_kind == :tolerance_robust`. Such a result is a declared numerical
classification, not an exact symbolic proof.

## Example 2: obtain and verify an extension witness

These two orthogonal two-qubit states are extendible:

```julia
julia> using LinearAlgebra

julia> left = [1 0; 0 1];

julia> right = [1 1; 0 0];

julia> result = is_upb(left, right; materialize_witness=true);

julia> (result.status, result.reason, result.is_upb)
(:not_upb, :extension_witness, false)

julia> result.witness_factors
(Rational{BigInt}[0, 1], Rational{BigInt}[0, 1])

julia> result.witness_vector
4-element Vector{Rational{BigInt}}:
 0
 0
 0
 1

julia> [
           prod(dot(result.witness_factors[j], factors[:, k])
                for (j, factors) in enumerate((left, right)))
           for k in axes(left, 2)
       ]
2-element Vector{Rational{BigInt}}:
 0
 0
```

`witness_factors` is the primary certificate. It stays compact in large tensor
products. `witness_vector` is allocated only when
`materialize_witness=true`, and its length is checked against
`max_dense_entries`.

`witness_partition[j]` lists the supplied states whose vanishing overlap is
certified by local witness factor `j`. Together the parts contain every state
index exactly once.

## Example 3: analyze global product-vector columns

The second form accepts a vector or matrix in the package tensor convention:

```julia
global_states = hcat(
    [tensor_product(left[:, k], right[:, k]) for k in axes(left, 2)]...,
)
result = is_upb(global_states, (3, 3); normalization=:allow)
```

Each global column is first factorized as a rank-one tensor. A column that is
definitely not a product vector returns

```text
status = :not_upb
reason = :not_product
```

and records its one-based column index in `offending_state`. A factorization
residual inside the tolerance boundary band returns `:unknown`.
`product_residual` records the largest factorization residual separately from
the orthogonality residual of any returned extension witness.

Subsystem 1 is the slowest-varying factor, consistently with
`tensor_product(a, b) == kron(a, b)`.

## What is certified

For supplied product states

```math
\lvert u_k\rangle =
\lvert u_{1,k}\rangle \mathbin{\otimes} \cdots
\mathbin{\otimes} \lvert u_{p,k}\rangle,
```

a product extension

```math
\lvert x\rangle =
\lvert x_1\rangle \mathbin{\otimes} \cdots
\mathbin{\otimes} \lvert x_p\rangle
```

is orthogonal to state `k` exactly when at least one local overlap vanishes:

```math
\langle x \mid u_k\rangle =
\prod_{j=1}^{p} \langle x_j \mid u_{j,k}\rangle = 0.
```

Assign each state index to one subsystem where its local overlap vanishes.
An extension exists exactly when there is an ordered partition
`(S₁, …, Sₚ)` for which every local family
`{u[j,k] : k in Sⱼ}` has rank smaller than `dⱼ`. A nonzero vector in each
local orthogonal complement gives the returned witness.

When the state count is at least `sum(dⱼ - 1)`, it is sufficient to enumerate
parts of size at least `dⱼ - 1`: any smaller deficient part can be enlarged to
that size while it still contains fewer than `dⱼ` vectors. When fewer states
were supplied, a deficient partition and witness exist immediately.

The implementation enumerates these partitions lazily and never stores the
full combinatorial list. This tested private algorithm supersedes the pinned
`vec_partitions` helper.

## Result meanings

| `status` | Typical `reason` | Meaning |
|---|---|---|
| `:upb` | `:unextendible` | The input is an incomplete orthogonal product family and every admissible partition was exhausted. |
| `:not_upb` | `:extension_witness` | A product vector in the orthogonal complement was constructed. |
| `:not_upb` | `:not_product` | A global input column is not a product vector. |
| `:not_upb` | `:not_orthogonal` | Two supplied product states are not orthogonal. |
| `:not_upb` | `:complete_basis` | The orthogonal product family fills the ambient space and is not incomplete. |
| `:unknown` | `:numerical_boundary` | A required floating-point zero or rank decision lies in the boundary band. |
| `:unknown` | `:state_limit` | `max_states` stopped input-dependent pair checks and recursive search depth. |
| `:unknown` | `:partition_limit` | `max_partitions` stopped the exhaustive search. |
| `:unknown` | `:work_limit` | `max_work` stopped local rank/complement work. |

An `:unknown` result is never converted into `false`.

## Normalization and tolerances

The default `normalization=:require` checks that each global product column has
unit norm. It never modifies input. Clearly non-unit input raises an
`ArgumentError`; a floating value inside the normalization boundary band
returns `:unknown`.

Use `normalization=:allow` to accept arbitrary nonzero scaling explicitly.
Floating analysis then uses private normalized local columns so rank decisions
are scale independent. This is recorded as `analysis_rescaled=true`; the
caller's arrays are untouched. Exact integer and rational paths do not
normalize and use exact elimination.

For floating inputs, the lower zero threshold at scale `n` is

```math
\mathit{atol} + n\,\mathit{rtol}.
```

Values strictly below it are treated as zero, values through
`boundary_factor` times that threshold are `:unknown`, and larger values are
treated as nonzero. The default is `atol=0`, `rtol=8eps(T)`, and
`boundary_factor=8`. Exact inputs require zero tolerances.

## Resource and storage controls

- `max_partitions` bounds complete partition candidates.
- `max_work` bounds a deterministic scalar-work estimate for local rank and
  complement calculations.
- `max_states=256` bounds pairwise validation work and recursive search depth.
- `max_dense_entries` guards private dense workspaces and optional full-witness
  materialization.
- `allow_densify=false` rejects sparse local or global input. Set it to `true`
  explicitly after choosing an appropriate dense-entry budget.

`max_states`, `max_partitions`, and `max_work` may be set to `nothing` to
disable those guards. A search limit produces a structured `:unknown`;
validation and pre-allocation violations throw before an unsafe operation.

The partition search is exponential in the worst case. For `s` states and `p`
parties it visits at most `p^s` assignments, with infeasible minimum-size
branches pruned. Each leaf performs bounded local rank work. The result reports
`partitions_examined` and `work_used`.

## Deliberate differences from the pinned routine

The pinned QETLAB implementation checks only extendibility. It can therefore
return true for a nonorthogonal unextendible product set or for a complete
product basis. The native API checks the full UPB definition.

The pinned witness calculation uses a nonconjugating transpose. For complex
local factors, the returned vector need not be orthogonal under the Hilbert
space inner product. This implementation solves the Hermitian-orthogonality
equations and verifies every returned witness.

The compatibility entry point is a thin
`MATLABCompat.IsUPB(U, V, ...; return_witness=false)` wrapper. It will opt into
`normalization=:allow` by default, returns a Boolean only for conclusive native
results, returns `(boolean, witness_factors)` when `return_witness=true`, and
raises a status-rich `DomainError` for `:unknown`. Set `structured=true` to
receive the full native result.

## Evidence

- Pinned source revision:
  `d8589610f00cff106537268dee2e2a1153f3a601`.
- Pinned files:
  `IsUPB.m` and `helpers/vec_partitions.m`.
- Primary definition:
  [Bennett *et al.* (1999)](https://arxiv.org/abs/quant-ph/9808030).
- Focused tests cover exact Tiles certification, floating and `BigFloat`
  behavior, global factorization, exact complex witnesses, nonproduct,
  nonorthogonal and complete-basis rejection, tolerance boundaries, sparse
  opt-in, invalid inputs, and deterministic limits.
- The function-specific Octave fixture compares conclusive pinned branches and
  records the pinned complex-witness and missing-definition-check deviations.
  Octave evidence is supplemental; MATLAB has not been run.
