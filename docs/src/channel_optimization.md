```@meta
CurrentModule = QuantumEntanglementTools
DocTestSetup = quote
    using QuantumEntanglementTools
    using LinearAlgebra
end
```

# Channel norms, discrimination, and output fidelity

This page covers four optimization-backed QETLAB entry points with one
solver-neutral native design:

- `diamond_norm` for the completely bounded trace norm;
- `cb_norm` for the completely bounded operator norm;
- `channel_distinguishability` for the optimal one-use success probability;
- `maximum_output_fidelity` for the largest root fidelity between channel
  outputs.

Every native call returns a package-owned result. A scalar is present only
after an analytic theorem applies or residual-checked lower and upper bounds
agree. Missing solvers, iteration limits, malformed optimizer factories, and
inconsistent numerical evidence remain explicit statuses.

## Diamond norm

For a map $\Phi:\mathrm{L}(\mathcal{X})\to\mathrm{L}(\mathcal{Y})$, the
diamond norm is the induced trace norm of
$\Phi\otimes\mathrm{id}_{\mathcal{X}}$. The solver-neutral builder uses the
Watrous primal SDP

```math
\begin{aligned}
\text{maximize}\quad
& \mathrm{Re}\langle J(\Phi),X\rangle,\\
\text{subject to}\quad
&
\begin{bmatrix}
\rho_0\otimes I_{\mathcal{Y}} & X\\
X^\dagger & \rho_1\otimes I_{\mathcal{Y}}
\end{bmatrix}\succeq0,\\
& \rho_0,\rho_1\succeq0,\\
& \mathrm{tr}(\rho_0)=\mathrm{tr}(\rho_1)=1.
\end{aligned}
```

The package Choi convention is input-first:
$J_{(i,a),(j,b)}=\Phi(E_{ij})_{a,b}$. This is why the diagonal blocks contain
$\rho\otimes I_{\mathcal{Y}}$. Input and output dimensions may differ, and
general `OperatorSumRepresentation` data are accepted, but both domain and
codomain must be square matrix algebras.

Completely positive maps use

```math
\lVert\Phi\rVert_\diamond
=\lVert\Phi^\ast(I_{\mathcal{Y}})\rVert_\infty
```

without a solver when complete positivity is established. Exact zero maps and
one-dimensional input or output spaces also have solver-free branches.

```jldoctest
identity_channel = KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
result = diamond_norm(identity_channel)

(
    result.status,
    result.value,
    result.certificate_kind,
)

# output
(QuantumEntanglementTools.ChannelOptimizationAnalyticOptimal, 1.0, :completely_positive_adjoint_identity)
```

`diamond_norm_problem(map)` always builds the inspectable affine SDP, even
when `diamond_norm(map)` could take an analytic shortcut.

## Completely bounded norm

Finite-dimensional duality gives

```math
\lVert\Phi\rVert_{\mathrm{cb}}
=\lVert\Phi^\ast\rVert_\diamond.
```

`cb_norm` applies `dual_channel`, delegates to the same reviewed diamond-norm
kernel, and records the adjoint transformation in the returned certificate.

```jldoctest
reset = KrausRepresentation([
    ComplexF64[1 0; 0 0],
    ComplexF64[0 1; 0 0],
])
result = cb_norm(reset)

(result.quantity, result.value)

# output
(:completely_bounded, 2.0)
```

The value is not generally one for a channel. In this example the reset
channel sends the identity to $2\lvert0\rangle\langle0\rvert$, so its
completely bounded norm is two.

## Channel discrimination

For normalized priors $p_1,p_2$, the channel Holevo--Helstrom theorem gives

```math
p_{\mathrm{success}}
=\frac{
p_1+p_2+
\lVert p_1\Phi_1-p_2\Phi_2\rVert_\diamond
}{2}.
```

Both inputs must be completely positive, trace preserving, and equal in input
and output dimensions. Priors must be finite, nonnegative, have length two,
and sum to one within the recorded tolerance. They are never normalized.

Identical channels and a zero prior have solver-free certificates:

```jldoctest
identity_channel = KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
result = channel_distinguishability(
    identity_channel,
    identity_channel;
    priors=[0.8, 0.2],
)

(
    result.status,
    result.success_probability,
    result.certificate_kind,
)

# output
(QuantumEntanglementTools.ChannelOptimizationAnalyticOptimal, 0.8, :identical_channels)
```

Without a backend, a nontrivial call still retains the always-guess lower
bound $\max(p_1,p_2)$ and the probability upper bound one.

## Maximum output fidelity

The native convention is root fidelity:

```math
F_{\max}(\Phi,\Psi)
=\max_{\rho,\sigma}
\left\lVert
\sqrt{\Phi(\rho)}\sqrt{\Psi(\sigma)}
\right\rVert_1.
```

The two input states are independent. The direct SDP is

```math
\begin{aligned}
\text{maximize}\quad
& \mathrm{Re}\,\mathrm{tr}(X),\\
\text{subject to}\quad
&
\begin{bmatrix}
\Phi(\rho) & X\\
X^\dagger & \Psi(\sigma)
\end{bmatrix}\succeq0,\\
& \rho,\sigma\succeq0,\\
& \mathrm{tr}(\rho)=\mathrm{tr}(\sigma)=1.
\end{aligned}
```

Identical channels, exact replacer channels, and channels with a common
computational-basis output have solver-free certificates. For example, the
identity and the phase-flip channel with Kraus operators $0.8I$ and $0.6Z$
both output $\lvert0\rangle\langle0\rvert$ on input
$\lvert0\rangle\langle0\rvert$:

```jldoctest
identity_channel = KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
phase_flip = KrausRepresentation([
    0.8Matrix{ComplexF64}(I, 2, 2),
    0.6ComplexF64[1 0; 0 -1],
])
result = maximum_output_fidelity(identity_channel, phase_flip)

(
    result.value,
    result.certificate_kind,
    result.output_states[1] == result.output_states[2],
)

# output
(1.0, :common_basis_output, true)
```

The default contract requires trace-preserving channels. The mathematically
broader completely-positive-map definition is available only through
`require_trace_preserving=false`; its value need not lie in $[0,1]$.

## Using an optional solver

General branches require JuMP plus an explicit conic optimizer. No optimizer
is selected globally or added to the core dependency graph.

```julia
using Hypatia, JuMP, QuantumEntanglementTools

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)

swap = ComplexF64[
    1 0 0 0
    0 0 1 0
    0 1 0 0
    0 0 0 1
]
transpose_map = ChoiRepresentation(swap, 2, 2)
result = diamond_norm(transpose_map; backend)

result.value
result.lower_bound
result.upper_bound
result.density_operators
result.witness_operator
result.optimization_result
```

The qubit transpose has diamond norm two. The dedicated extension tests check
this value with Hypatia and SCS, along with a rectangular two-sided elementary
map, channel discrimination against the completely depolarizing channel, and
the direct maximum-output-fidelity SDP.

## Reading bounds and statuses

`ChannelNormResult`, `ChannelDistinguishabilityResult`, and
`MaximumOutputFidelityResult` separate:

- a residual-checked primal lower bound;
- an available dual or optimizer upper bound;
- a scalar value only when bounds agree;
- reconstructed density operators, output states, block couplings, or norm
  witnesses;
- generic termination, primal, and dual statuses;
- primal and dual residuals, gaps, iterations, timing, optimizer versions, and
  options;
- theorem-level analytic certificates, separately from numerical solver agreement.

The last distinction matters: a successful numerical solve is not marked
`certified=true`. It remains a residual-checked numerical optimum.

## Validation, storage, and limits

No map, prior, density operator, or Choi matrix is normalized, symmetrized,
clipped, padded, or projected.

- `max_input_dimension`, `max_output_dimension`, `max_choi_dimension`, and
  `max_choi_entries` guard representation growth before model construction.
- `max_variables` guards affine coordinates.
- `OptimizationLimits` guards equalities, PSD blocks, block dimensions,
  stored coefficients, and the real-block cone translation.
- `max_dense_entries` guards spectral validation.
- Sparse spectral shortcuts require `allow_densify=true`. The JuMP backend
  also requires its own explicit `allow_densify=true` because a complex
  Hermitian PSD cone becomes a dense real block.

Numerical complete-positivity or trace-preservation boundaries are not
silently repaired into channels.

## Pinned QETLAB deviations

The committed solver-free Octave fixture records two upstream defects.

First, pinned `ChannelDistinguishability.m` returns only
$\lVert p_1\Phi_1-p_2\Phi_2\rVert_\diamond$ in its non-deterministic branch.
For identical channels with priors $(0.8,0.2)$ it therefore returns $0.6$
while labeling the result a probability; the theorem and this package return
$0.8$.

Second, pinned `MaximumOutputFidelity.m` retains only the first
`min(Kraus ranks)` canonical Kraus operators when constructing its cross map.
For the identity and the phase-flip channel above it returns $0.8$, even
though their displayed common output certifies the exact value one. The native
implementation uses the direct fidelity SDP and is invariant under Kraus
representation.

The pinned channel-discrimination validation also accepts a negative prior
such as $p=(1.1,-0.1)$. The native API rejects it. General CVX branches are not
claimed as Octave evidence; they are exercised independently with Hypatia and
SCS.

The formulations and the corrected channel-discrimination probability follow
[John Watrous, *The Theory of Quantum Information*, Sections
3.3.3--3.3.4](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).

The complete docstrings are collected in the
[Julia-native optimization API](api/native_optimization.md).
