# Matrix predicates

```@meta
DocTestSetup = quote
    using QuantumEntanglementTools
end
```

The matrix-predicate slice maps QETLAB `IsPSD`, `IsLocallyPSD`,
`IsTotallyPositive`, and `IsTotallyNonsingular` to structured Julia results.
Every call returns a `MatrixPredicateResult` whose status is
`MatrixPredicateSatisfied`, `MatrixPredicateViolated`, or
`MatrixPredicateUnknown`. Floating-point values inside an explicit tolerance
band are never collapsed to a Boolean.

## Positive semidefiniteness

`is_positive_semidefinite(A)` requires a nonempty square matrix. Integer and
rational inputs are checked exactly. Floating BLAS types use a Hermitian
eigendecomposition; `BigFloat` and other supported types use a generic LDL
congruence calculation. The input is not symmetrized: a Hermiticity defect
inside tolerance is `unknown`, and one outside tolerance violates the
predicate.

```jldoctest matrix-predicates
julia> is_positive_semidefinite([2.0 1.0; 1.0 2.0]).status === MatrixPredicateSatisfied
true

julia> is_positive_semidefinite(zeros(2, 2)).status === MatrixPredicateUnknown
true

julia> is_positive_semidefinite(Rational{Int}[1 1; 1 1]).status === MatrixPredicateSatisfied
true
```

`is_locally_positive_semidefinite(A, k)` applies the same contract to all
`binomial(size(A,1), k)` principal submatrices. Its witness records the
principal indices and nested PSD result.

## All-minor predicates

`is_totally_positive(A; orders=nothing)` requires strict positivity of every
requested square minor of a real matrix. `is_totally_nonsingular` accepts real
or complex matrices and requires nonzero determinants. Exact integer/rational
arithmetic is widened. A represented singular floating minor is a violation;
a nonzero determinant inside tolerance is `unknown`.

```jldoctest matrix-predicates
julia> is_totally_positive(Rational{Int}[1 1; 1 2]).status === MatrixPredicateSatisfied
true

julia> is_totally_nonsingular([1.0 1.0; 1.0 1.0 + 1e-12];
           orders=2, atol=1e-8, rtol=0).status === MatrixPredicateUnknown
true
```

All-minor work is combinatorial:
`sum(binomial(m,k) * binomial(n,k) for k in orders)` determinants, with
`O(k^3)` arithmetic per determinant. `max_minors=100_000` and
`max_submatrices=100_000` guard that work; `nothing` explicitly disables a
guard. Sparse input requires `allow_densify=true`; diagonal PSD input has a
structure-aware linear path.

## Compatibility entry points

`MATLABCompat.IsPSD`, `IsLocallyPSD`, `IsTotallyPositive`, and
`IsTotallyNonsingular` preserve QETLAB spellings and positional arguments but
return the same structured result. They do not reproduce boundary-collapsing
Booleans. `MATLABCompat.IsPSD` also omits the CVX branch and never substitutes
the Hermitian part. The all-minor wrappers preserve QETLAB's single default
determinant tolerance while treating singular minors according to the strict
mathematical definitions.

Full signatures are in the [API reference](api/index.md).
