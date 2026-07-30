# Absolute PPT from a spectrum

An operator is **absolutely PPT** when its partial transpose remains positive
semidefinite after every global unitary conjugation. Unlike ordinary PPT, this
is a property of the spectrum and the bipartite dimensions: changing the
eigenvectors cannot change the answer.

`abs_ppt_constraints` constructs the finite spectral LMI family used by
QETLAB's `AbsPPTConstraints`, while `is_abs_ppt` evaluates that family and
returns a certificate-aware [`IsAbsPPTResult`](@ref). The native API never
reduces an inconclusive computation to `false`.

Absolute PPT and separability are different notions. Absolute PPT rules out
the creation of NPT entanglement by a global unitary; it does not, in general,
certify separability or absolute separability.

## The finite LMI criterion

Let the real spectrum be ordered as

```math
\lambda_1\geq\lambda_2\geq\cdots\geq\lambda_D\geq0,
\qquad D=d_A d_B,
```

and put $p=\min(d_A,d_B)$. Hildebrand's criterion considers the possible
orderings of products $x_i x_j$ with $1\leq i\leq j\leq p$. Each admitted
ordering determines a matrix $\Lambda$. The associated constraint is

```math
\Lambda+\Lambda^{\mathsf T}\succeq0.
```

The diagonal and upper-triangle entries use the smallest eigenvalues, in
reverse product-order rank. The strictly lower-triangle contribution
subtracts the largest eigenvalues in the induced relative ordering of the
off-diagonal products. An [`AbsPPTOrdering`](@ref) records both orderings
explicitly:

- `positive_pairs[r] == (i, j)` means that $x_i x_j$ has rank `r`, with the
  largest product first;
- `negative_pairs` contains the pairs with $i<j$, in the same relative order.

The enumeration uses QETLAB's monotonicity and criss-cross exclusion rules. At
local dimensions $p=2,3,4,5,6$, it produces respectively
`1, 2, 10, 114, 2612` matrices. The $p=6$ QETLAB family contains four
[documented redundant orderings](https://njohnston.ca/2014/02/counting-the-possible-orderings-of-pairwise-multiplication/)
beyond the 2608 realizable orderings. This is why a negative result is accepted
only after the offending ordering has also received an exact realization
certificate.

The mathematical source is R. Hildebrand,
[“Positive partial transpose from spectra”](https://doi.org/10.1103/PhysRevA.76.052325),
Phys. Rev. A **76**, 052325 (2007); the
[author manuscript](https://arxiv.org/abs/quant-ph/0502170) is openly
available. The port follows QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`.

## Two-qubit example

For $p=2$, there is one matrix:

```math
M(\lambda)=
\begin{bmatrix}
2\lambda_4 & \lambda_3-\lambda_1\\
\lambda_3-\lambda_1 & 2\lambda_2
\end{bmatrix}.
```

Thus the two-qubit absolute-PPT condition is

```math
\lambda_1\leq\lambda_3+2\sqrt{\lambda_2\lambda_4}.
```

The state with spectrum $(0.45,0.35,0.10,0.10)$ satisfies the inequality
strictly:

```@example absolute-ppt-two-qubit
using QuantumEntanglementTools

spectrum = [0.45, 0.35, 0.10, 0.10]
family = abs_ppt_constraints(spectrum; dims=(2, 2))
constraint = family.constraints[1]

@assert family.exhaustive
@assert length(family.constraints) == 1
@assert constraint.matrix ≈ [
    0.20 -0.35
   -0.35  0.70
]

result = is_abs_ppt(spectrum; dims=(2, 2))
@assert result.status ===
        QuantumEntanglementTools.AbsolutePPTExhaustiveCertified
@assert result.verdict === true

(
    matrix=constraint.matrix,
    status=result.status,
    verdict=result.verdict,
    minimum_margin=result.margin,
)
```

For the pure spectrum $(1,0,0,0)$, the same matrix has a negative eigenvalue.
The result retains both that PSD diagnostic and an exact rational realization
of the product ordering:

```@example absolute-ppt-two-qubit-negative
using QuantumEntanglementTools

result = is_abs_ppt([1.0, 0.0, 0.0, 0.0]; dims=(2, 2))

@assert result.status === QuantumEntanglementTools.AbsolutePPTCertifiedNot
@assert result.verdict === false
@assert result.violating_constraint.diagnostic.value < 0
@assert result.ordering_certificate.minimum_gap > 0
@assert result.ordering_certificate.minimum_product_margin > 0

(
    status=result.status,
    violating_eigenvalue=result.violating_constraint.diagnostic.value,
    certificate_source=result.ordering_certificate.source,
)
```

The certificate uses log gaps
$g_i=\log x_i-\log x_{i+1}$. Every strict product-order inequality is
rechecked in rational arithmetic before the negative verdict is returned.

## A three-by-three construction

For $p=3$, exactly two criss-cross orderings survive. Using the descending
spectrum $(9,8,\ldots,1)$ makes the indexing visible:

```@example absolute-ppt-three-by-three
using QuantumEntanglementTools

spectrum9 = collect(9:-1:1)
family3 = abs_ppt_constraints(
    spectrum9;
    dims=(3, 3),
    max_constraints=nothing,
    max_work=nothing,
    max_entries=nothing,
)

@assert family3.exhaustive
@assert length(family3.constraints) == 2
@assert family3.known_criss_cross_count == 2

matrices = Set(constraint.matrix for constraint in family3.constraints)
@assert matrices == Set([
    [2 -7 -4; -7 6 -2; -4 -2 12],
    [2 -7 -5; -7 8 -2; -5 -2 12],
])

[
    (
        positive_pairs=collect(constraint.ordering.positive_pairs),
        negative_pairs=collect(constraint.ordering.negative_pairs),
        matrix=constraint.matrix,
    ) for constraint in family3.constraints
]
```

These integer eigenvalues are only a transparent construction example; they
are not normalized. Neither function requires trace one because the
absolute-PPT inequalities are homogeneous.

## Reading `is_abs_ppt`

[`IsAbsPPTResult`](@ref) separates theorem-level conclusions from resource and
numerical outcomes:

| Status | `verdict` | Meaning |
|---|---:|---|
| `AbsolutePPTAnalyticCertified` | `true` | Zero operator or one-dimensional local factor |
| `AbsolutePPTSufficientTestPassed` | `true` | A scale-invariant sufficient test passed |
| `AbsolutePPTExhaustiveCertified` | `true` | Every member of the exhaustively enumerated family passed |
| `AbsolutePPTCertifiedNot` | `false` | A robust negative LMI and an exactly verified realizable ordering were found |
| `AbsolutePPTCappedUnknown` | `nothing` | A work, storage, enumeration, or deterministic-realization budget ended first |
| `AbsolutePPTNumericalBoundary` | `nothing` | Input validation or an LMI lies inside the floating tolerance boundary |
| `AbsolutePPTBackendUnavailable` | `nothing` | A realization solve was needed but no optional backend was configured |
| `AbsolutePPTBackendFailure` | `nothing` | The backend failed or its point did not pass exact re-verification |

Inspect `status` and `verdict` together. In particular, neither a capped family
with no observed violation nor a floating boundary is a positive result.

The default sufficient branches are the scale-invariant
Gurvits--Barnum separable-ball inequality and a Gershgorin bound for all
Hildebrand matrices. Set `use_sufficient_tests=false` when the purpose is to
inspect the complete LMI path:

```@example absolute-ppt-statuses
using QuantumEntanglementTools

fast = is_abs_ppt(fill(0.25, 4); dims=(2, 2))
@assert fast.status ===
        QuantumEntanglementTools.AbsolutePPTSufficientTestPassed
@assert fast.certificate_kind === :gurvits_barnum_separable_ball

capped = is_abs_ppt(
    fill(1 / 9, 9);
    dims=(3, 3),
    use_sufficient_tests=false,
    max_constraints=1,
)
@assert capped.status === QuantumEntanglementTools.AbsolutePPTCappedUnknown
@assert capped.verdict === nothing
@assert capped.checked_constraints == 1
@assert capped.planned_orderings == 2

(fast=fast, capped=capped)
```

Exact rational data can certify an equality exactly. The corresponding
floating computation is conservatively reported as a numerical boundary:

```@example absolute-ppt-boundary
using QuantumEntanglementTools

exact = is_abs_ppt(
    [1 // 2, 1 // 6, 1 // 6, 1 // 6];
    dims=(2, 2),
)
floating = is_abs_ppt(
    [0.5, 1 / 6, 1 / 6, 1 / 6];
    dims=(2, 2),
)

@assert exact.verdict === true
@assert floating.status ===
        QuantumEntanglementTools.AbsolutePPTNumericalBoundary
@assert floating.verdict === nothing

(exact=exact.status, floating=floating.status)
```

## Affine spectra and optional solvers

An affine eigenvalue model returns package-owned
[`HermitianAffineMatrix`](@ref) values, never JuMP expressions:

```@example absolute-ppt-affine
using QuantumEntanglementTools

variables = 4
affine_spectrum = [
    AffineScalar(
        0.0,
        [index == variable ? 1.0 : 0.0 for index in 1:variables],
    ) for variable in 1:variables
]

affine_family = abs_ppt_constraints(
    affine_spectrum;
    dims=(2, 2),
    assume_ordered=true,
)

point = [0.4, 0.3, 0.2, 0.1]
evaluated = Matrix(
    QuantumEntanglementTools.evaluate_affine(
        affine_family.constraints[1].matrix,
        point,
    ),
)

@assert evaluated ≈ [0.2 -0.2; -0.2 0.6]
(symbolic=affine_family.symbolic, matrix=evaluated)
```

The surrounding model must impose
$\lambda_1\geq\cdots\geq\lambda_D\geq0$; symbolic expressions cannot be
silently sorted.

Most negative-ordering certificates are found by a deterministic bounded
search. If that search ends first, an optional JuMP backend can solve the
package-owned feasibility model:

```julia
using QuantumEntanglementTools
using Hypatia

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)

result = is_abs_ppt(
    [1.0, 0.0, 0.0, 0.0];
    dims=(2, 2),
    max_realization_iterations=0,
    backend,
)
```

The optional solver only proposes log-gap coordinates. The core converts them
to exact rationals and checks every ordering inequality independently. A
nominally feasible solver status by itself cannot produce
`AbsolutePPTCertifiedNot`.

## Limits, dimensions, types, and sparse inputs

Enumeration is combinatorial. [`AbsPPTConstraintFamily`](@ref) records:

- the monotone-ordering upper count and, when known, the QETLAB criss-cross
  count;
- candidates checked, work charged, and matrix entries stored;
- the active `max_constraints`, `max_work`, and `max_entries` limits;
- whether enumeration was exhaustive;
- the first robust violation or numerical boundary when
  `stop_on_violation=true`.

Defaults are deliberately finite. Passing `nothing` disables a particular
limit and should be an explicit decision:

```@example absolute-ppt-limits
using QuantumEntanglementTools

partial = abs_ppt_constraints(
    collect(9.0:-1:1);
    dims=(3, 3),
    max_constraints=1,
)

@assert partial.status ===
        QuantumEntanglementTools.AbsPPTEnumerationConstraintLimit
@assert !partial.exhaustive
@assert length(partial.constraints) == 1

(
    status=partial.status,
    checked=partial.candidate_orderings_checked,
    work=partial.work_used,
    stored=partial.entries_stored,
)
```

Dimension tuples use the package's left-to-right tensor-factor convention.
`dims=nothing` is accepted only when the total dimension is a perfect square.
A scalar `dims=d` means `(d, D ÷ d)` in both native functions. This consistent
native rule deliberately corrects the incompatible scalar-dimension behavior
of the two pinned MATLAB routines.

Numeric spectra are sorted in descending order. They must be real, finite, and
one-based; no value is normalized or clipped. Matrix input must be square,
finite, and Hermitian. A tiny Hermiticity defect within floating tolerance
returns `AbsolutePPTNumericalBoundary` from `is_abs_ppt`, while a larger defect
is an error. Wrap a matrix in `Hermitian` only when that is the caller's
intended mathematical assertion.

`Rational`, `Float32`, and `BigFloat` spectrum vectors preserve their meaningful
element type. Sparse spectra and general sparse Hermitian matrices require
`allow_densify=true`; a sparse diagonal matrix can be read without that
permission. Use `sparse_output=true` to request sparse numeric LMI matrices.

## Differences from the pinned MATLAB behavior

The native API intentionally makes several ambiguities observable:

- scalar dimensions have one consistent meaning, described above;
- the $p=1$ family is empty but exhaustive, and `is_abs_ppt` returns an analytic
  certificate instead of failing on an empty minimum;
- the $p=6$ QETLAB count of 2612 is retained, while negative verdicts require
  exact ordering realizability;
- complex “eigenvalue vectors” are rejected instead of silently discarding
  imaginary parts;
- general sparse matrix eigendecomposition requires explicit densification;
- floating equality is a boundary, not an unconditional Boolean pass; and
- all limits, backend failures, and incomplete searches have structured
  statuses.

For migrating QETLAB code, `MATLABCompat.AbsPPTConstraints` keeps the
positional `DIM`, `ESC_IF_NPOS`, and `LIM` arguments and returns a vector of
matrices. Its scalar `DIM=d` deliberately retains the pinned `(d,d)` meaning.
`MATLABCompat.IsAbsPPT` retains the other pinned scalar rule and maps certified
yes, certified no, and every inconclusive status to `1`, `0`, and `-1`
respectively. Both wrappers accept `structured=true` to expose the native
evidence. The plain constraint wrapper raises if an implicit resource guard
produces an unexpected partial family.

The source-free oracle fixture records pinned QETLAB matrices and legacy
Boolean outputs separately from these documented native corrections.
