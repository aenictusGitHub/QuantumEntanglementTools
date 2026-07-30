# Coherence optimization and absolute incoherence

This page covers five status-rich native operations:

- `is_k_incoherent` tests whether a state has coherence number at most `k`;
- `is_absolutely_k_incoherent` asks whether every unitary conjugate is
  `k`-incoherent;
- `robustness_coherence` computes robustness against incoherent states;
- `trace_distance_coherence` minimizes trace distance to an incoherent state;
- `generalized_robustness_k_coherence` computes robustness against states with
  coherence number at most `k`.

The first two return `CoherenceCriterionResult`. The three measures return
`CoherenceOptimizationResult`. These objects keep theorem provenance, solver
status, residuals, free states, noise, and decompositions separate. A missing
backend, iteration limit, or numerical boundary is never converted into a
mathematical `false` or a fabricated scalar.

## Input contract

`state` may be a normalized floating-point vector or a floating-point density
matrix. Matrices must be square, finite, exactly Hermitian, positive
semidefinite within the requested tolerance, and trace one within that
tolerance. Vectors must have squared norm one. Validation never normalizes,
symmetrizes, clips eigenvalues, or projects an input onto the state space.

`k` is an integer from `1` through the Hilbert-space dimension. The free sets
are nested:

```math
\mathcal{C}_1 \subseteq \mathcal{C}_2 \subseteq \cdots
\subseteq \mathcal{C}_d.
```

Here $\mathcal{C}_1$ is the set of diagonal density matrices and
$\mathcal{C}_d$ is the full state space.

Sparse inputs require `allow_densify=true`, because density validation and the
portable complex-Hermitian SDP representation perform dense spectral work.
This explicit opt-in prevents an accidental large allocation.

## Analytic example

For the maximally coherent four-level state
$|\psi\rangle=(1,1,1,1)^\mathsf{T}/2$, the exact branches need no solver:

```julia
using QuantumEntanglementTools

psi = fill(0.5 + 0im, 4)

robustness = robustness_coherence(psi)
robustness.value                 # 3.0
robustness.method                # :pure_state_k_support_norm
robustness.free_state            # optimal diagonal density matrix
robustness.noise_state           # normalized optimal noise

distance = trace_distance_coherence(psi)
distance.value                   # 1.5
distance.free_state              # I / 4

level_two = generalized_robustness_k_coherence(psi, 2)
level_two.value                  # 1.0
level_two.decomposition.supports # supports of the optimal 2-sparse blocks
```

For all three robustness calls, `unnormalized_noise` is retained and satisfies
the cone equality with the returned free state. When the robustness is exactly
zero, `noise_state` is `nothing` and `unnormalized_noise` is the zero matrix.
The implementation never evaluates `zero_matrix / 0`.

## The reconstructed free cone

A positive matrix is `k`-incoherent exactly when it has factor width at most
`k`. Equivalently, it can be decomposed as

```math
Y=\sum_{\substack{S\subseteq\{1,\ldots,d\}\\|S|=k}}
E_S A_S E_S^*,\qquad A_S\succeq 0,
```

where $E_S$ embeds the coordinates indexed by $S$. This definition supplies
the capability that the pinned `GenRobustnesskCoherence.m` source expected
from an absent `IskCoherent.m` file.

The exact feasibility model for `is_k_incoherent` sets $Y=\rho$. Generalized
robustness solves

```math
\begin{aligned}
\mathop{\mathrm{minimize}}\quad & \mathrm{tr}(X)\\
\mathrm{subject\ to}\quad
& \rho+X=\sum_S E_S A_S E_S^*,\\
& X\succeq0,\qquad A_S\succeq0.
\end{aligned}
```

The result retains the unnormalized and normalized free blocks, the optimal
noise, the solver result, and the complete solver-neutral
`SemidefiniteProgram`.

For pure states, the implementation uses the exact `k`-support-norm theorem.
It also constructs an optimizer: a deterministic fixed-cardinality marginal
decomposition produces `k`-sparse atoms whose covariance is the positive
noise. Thus the analytic branch returns more than the scalar theorem value.

## Testing `k`-incoherence

With `strategy=:auto`, `is_k_incoherent` tries the following theorem-backed
routes before building the general SDP:

1. diagonal and full-dimension endpoint characterizations;
2. the comparison-matrix criterion, which is necessary and sufficient for
   `k == 2` and sufficient for larger `k`;
3. the purity-ball sufficient condition;
4. an exact but bounded graph-bandwidth recognition search;
5. the dephasing positive-semidefinite sufficient condition;
6. the factor-width block feasibility SDP.

A sufficient theorem can certify `true`, but failure of a one-sided theorem
does not certify `false`. If the SDP cannot be run, `verdict` remains
`nothing`.

### Corrected band-ordering helper

The pinned private `has_band_k_ordering.m` defines its reversal helper with
arguments `(unselected, candidate, num_placed)` but calls it as
`(unselected, num_placed, candidate)`. This is not just cosmetic: for

```math
A=\begin{bmatrix}
0&1&1\\
1&0&0\\
1&0&0
\end{bmatrix},
```

the pinned helper reports no one-based bandwidth-2 ordering, although
`[2, 1, 3]` is one.

The native implementation does not translate the faulty pruning rule. It
performs an independently specified exact layout search, bounded by
`max_band_search_nodes`, and verifies the returned ordering. If the budget is
exhausted, the band test is inconclusive and the function continues to another
certificate or the SDP.

The one-based convention is deliberate: a zero matrix has bandwidth `0`, a
diagonal matrix has bandwidth `1`, and an entry at positions `i,j` requires
$|i-j|<k$ for bandwidth at most $k$.

## Absolute `k`-incoherence

Absolute `k`-incoherence depends only on the spectrum. The automatic
classifier implements the primary-source conclusions separately:

- only the maximally mixed state is absolutely 1-incoherent;
- every state is absolutely `d`-incoherent;
- rank at most $d-k$ rules the property out;
- the tight-rank state with equal nonzero eigenvalues certifies it;
- a maximal eigenvalue at most $1/(d-k+1)$ is sufficient;
- the purity inequality is necessary and sufficient for `k == 2` in
  dimensions at most three and sufficient in larger dimensions;
- for `k == d-1`, a theorem-backed real positive-semidefinite feasibility
  model is necessary and sufficient.

For other parameter ranges, the cited spectral results are one-sided. The
correct result may therefore be `status == :unknown` and `verdict === nothing`.
This is an intentional mathematical distinction.

## Trace-distance formulation

The general trace-distance problem introduces positive semidefinite matrices
$P,N$ and a diagonal probability vector $p$:

```math
\begin{aligned}
\mathop{\mathrm{minimize}}\quad & \mathrm{tr}(P)+\mathrm{tr}(N)\\
\mathrm{subject\ to}\quad
& \rho-\mathrm{diag}(p)=P-N,\\
& P\succeq0,\quad N\succeq0,\quad p\geq0,\quad
\sum_i p_i=1.
\end{aligned}
```

`positive_part`, `negative_part`, and `free_state` retain that primal
decomposition. Pure states use the exact closed form, including the closest
diagonal state and branch diagnostics. Qubits use
$2|\rho_{12}|$. The native qubit result consistently returns a density matrix;
the pinned MATLAB branch returns a diagonal matrix even though its help text
describes a vector.

## Selecting a backend

The core package never chooses or loads a solver implicitly. After loading
JuMP and a supported conic optimizer, construct an explicit backend:

```julia
using QuantumEntanglementTools
using JuMP
using Hypatia

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
    atol=1e-7,
    rtol=1e-7,
)

rho = ComplexF64[
    0.40  0.05   0.02im
    0.05  0.35   0.03
   -0.02im 0.03  0.25
]

result = robustness_coherence(rho; strategy=:sdp, backend)
result.status
result.value
result.optimization.primal_residual
result.optimization.optimizer
```

`strategy=:sdp` is useful for cross-checking an analytic branch. Hypatia and
SCS are tested independently, but neither is a core dependency.

Model construction is guarded by `OptimizationLimits` and, for factor-width
models, `max_subsets`. Exceeding a preflight limit returns
`status == :resource_limit` before allocating a solver model.

## Result interpretation

Common criterion statuses include:

- `:certified_true` and `:certified_false` for theorem conclusions;
- `:solver_feasible` and `:solver_infeasible` for numerical SDP outcomes;
- `:boundary` and `:unknown` for unresolved mathematical cases;
- `:backend_unavailable`, `:solver_limit`, `:malformed_backend`,
  `:inconsistent`, and `:resource_limit` for operational failures.

Common measure statuses include `:analytic_exact`, `:optimal`, `:feasible`,
and the same operational failure statuses. A consistent numerical optimum is
not marked `exact`; raw backend termination, primal/dual status, objective
bounds, residuals, options, and version metadata remain in `optimization`.

## References

The formulations and theorem branches follow:

- Ringbauer *et al.*, *Certification and Quantification of Multilevel Quantum
  Coherence*, Physical Review X **8**, 041007 (2018);
- Johnston, Moein, Pereira, and Plosker, *Absolutely k-Incoherent Quantum
  States and Spectral Inequalities for Factor Width of a Matrix*, Physical
  Review A **106**, 052417 (2022);
- Napoli *et al.*, *Robustness of Coherence*, Physical Review Letters **116**,
  150502 (2016);
- Piani *et al.*, *Robustness of Asymmetry and Coherence of Quantum States*,
  Physical Review A **93**, 042107 (2016);
- Chen *et al.*, *Quantifying the Coherence of Pure Quantum States*, Physical
  Review A **94**, 042313 (2016);
- Johnston *et al.*, *Some Notes on the Robustness of k-Coherence and
  k-Entanglement*, Physical Review A **98**, 022328 (2018).

The implementation is source-informed but independent. QETLAB source
dispositions, exact pinned hashes, solver-free fixtures, and the corrected
band-ordering counterexample are recorded in the repository evidence files.

The complete docstrings are collected in the
[Julia-native optimization API](api/native_optimization.md).
