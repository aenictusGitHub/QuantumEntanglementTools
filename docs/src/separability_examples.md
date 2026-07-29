# Separability by example

A bipartite state is separable when it can be written as a convex mixture of
product states. For a pure state, this reduces to being a single product
vector. Deciding separability for an arbitrary mixed state is difficult, so
this package returns structured conclusions instead of a potentially
misleading `Bool`.

There is intentionally no general `is_separable` function. A result of
`unknown` means that the methods run did not produce a certificate; it does
not mean that the state is entangled or separable.

## Choosing an API

| Input and question | Start with | A conclusive result | An inconclusive result |
|---|---|---|---|
| Bipartite pure vector | `analyze_entanglement(psi, dims)` | `status == :separable` with an exact product-decomposition certificate | `status == :unknown` |
| Density matrix in `2×2` or `2×3` | `analyze_entanglement(rho, dims)` | `status == :separable` when the low-dimensional PPT theorem applies | `status == :unknown` |
| Bipartite density matrix in any supported dimensions | `in_separable_ball(rho, dims)` | `status == :separable_certified` inside the Gurvits--Barnum ball | `:outside_ball` or `:unknown`; neither implies entanglement |
| Higher-dimensional density matrix, looking for entanglement | `analyze_entanglement(rho, dims)` | A criterion violation can return certified `:entangled` | Passing the requested necessary criteria returns `:unknown` |

`analyze_entanglement` returns an
[`EntanglementReport`](@ref), whose status is `:separable`, `:entangled`, or
`:unknown`. Inspect `certified`, `certificate_kind`, and `message` before using
the conclusion. `in_separable_ball` returns a [`SeparableBallResult`](@ref)
with the separate status vocabulary `:separable_certified`, `:outside_ball`,
or `:unknown`.

The density-matrix pipeline currently tries PPT, realignment, and reduction
criteria. It does **not** call `in_separable_ball`; run the ball test explicitly
when it is relevant.

## 1. A pure product state

The vector ``|0\rangle\otimes|1\rangle`` is visibly a product state. Subsystem
dimensions are supplied from left to right, so `(2, 2)` assigns two levels to
each factor.

```@example separability-pure
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
psi01 = tensor_product(ket0, ket1)

report = analyze_entanglement(psi01, (2, 2))

@assert report.status === :separable
@assert report.certified
@assert report.certificate_kind === :pure_product_decomposition

(
    status=report.status,
    certified=report.certified,
    certificate=report.certificate_kind,
    message=report.message,
)
```

Here the package computes a Schmidt decomposition. Exactly zero trailing
Schmidt coefficients provide the named pure-product certificate.

## 2. An explicit mixed separable state

Consider the full-rank two-qubit mixture

```math
\rho =
\frac{3}{8}|00\rangle\langle 00|
+ \frac{1}{8}|01\rangle\langle 01|
+ \frac{1}{8}|10\rangle\langle 10|
+ \frac{3}{8}|11\rangle\langle 11|.
```

The displayed convex decomposition already makes the state separable. The code
constructs that decomposition directly and then obtains two package-owned
sufficient certificates.

```@example separability-mixed-two-qubit
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
product_basis = [
    tensor_product(ket0, ket0),
    tensor_product(ket0, ket1),
    tensor_product(ket1, ket0),
    tensor_product(ket1, ket1),
]
weights = [3 / 8, 1 / 8, 1 / 8, 3 / 8]
rho = sum(
    weights[index] * (product_basis[index] * product_basis[index]')
    for index in eachindex(weights)
)

pipeline = analyze_entanglement(rho, (2, 2); atol=1e-12, rtol=0)
ball = in_separable_ball(rho, (2, 2); atol=1e-12, rtol=0)

@assert pipeline.status === :separable
@assert pipeline.certified
@assert pipeline.certificate_kind === :ppt_low_dimension_theorem
@assert ball.status === :separable_certified

(
    pipeline_status=pipeline.status,
    pipeline_certificate=pipeline.certificate_kind,
    ball_status=ball.status,
    purity=ball.purity,
    ball_boundary=ball.boundary,
)
```

The two status spellings belong to different structured result types. The
pipeline certifies the state through the exact PPT theorem for two qubits. The
ball test independently certifies it because its purity is strictly below the
reported sufficient bound.

## 3. A higher-dimensional state needs the right certificate

In the ordered product basis
``|00\rangle,|01\rangle,\ldots,|22\rangle``, the following diagonal
``3\times3`` bipartite state is an explicit mixture of nine product-basis
projectors:

```math
\rho_3 = \operatorname{diag}
\left(\frac18,\frac18,\frac18,\frac18,\frac18,\frac18,\frac18,
\frac1{16},\frac1{16}\right).
```

The native pipeline cannot infer separability merely because its necessary
criteria pass in these dimensions. The separable-ball test supplies the
missing sufficient certificate.

```@example separability-three-by-three
using LinearAlgebra
using QuantumEntanglementTools

weights3 = [fill(1 / 8, 7); 1 / 16; 1 / 16]
rho3 = Matrix(Diagonal(weights3))

pipeline3 = analyze_entanglement(rho3, (3, 3); atol=1e-12, rtol=0)
ball3 = in_separable_ball(rho3, (3, 3); atol=1e-12, rtol=0)

@assert pipeline3.status === :unknown
@assert !pipeline3.certified
@assert ball3.status === :separable_certified

attempts = [attempt.method => attempt.status for attempt in pipeline3.attempts]
(
    pipeline_status=pipeline3.status,
    pipeline_attempts=attempts,
    ball_status=ball3.status,
    purity=ball3.purity,
    ball_boundary=ball3.boundary,
)
```

These results are not contradictory. They apply different criteria: the
pipeline found no certificate, while the ball test satisfied a sufficient
separability condition.

## 4. Outside the separable ball does not mean entangled

The separable ball is centered on the maximally mixed state and is only a
sufficient test. A pure product density matrix lies outside that ball even
though its underlying vector has an exact product certificate.

```@example separability-outside-ball
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
psi00 = tensor_product(ket0, ket0)
rho00 = psi00 * psi00'

pure_report = analyze_entanglement(psi00, (2, 2))
ball_result = in_separable_ball(rho00, (2, 2); atol=1e-12, rtol=0)

@assert pure_report.status === :separable
@assert pure_report.certified
@assert ball_result.status === :outside_ball

(
    known_state=pure_report.status,
    known_certificate=pure_report.certificate_kind,
    ball_status=ball_result.status,
    ball_message=ball_result.message,
)
```

Never negate a one-sided certificate: `:outside_ball` says only that this
particular test did not certify the state. The same rule applies to `unknown`
pipeline and optional-backend results.

For the complete criterion ordering, tolerance behavior, and witness semantics,
see [Entanglement backends](entanglement_backends.md). For the separable-ball
formula and structured input paths, see
[Product structure and separable-ball certificates](product_analysis.md).
