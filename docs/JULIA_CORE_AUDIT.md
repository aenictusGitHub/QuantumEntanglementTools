# Julia core audit

Audit date: 2026-07-29.

This report reviews the Julia-native core, compatibility layer, optional
EntanglementDetection extension, focused tests, benchmark evidence, and the
architecture/conventions documents. It is a release-hardening audit, not a
claim of complete QETLAB parity. A public symbol being present does not by
itself establish that the corresponding inventory row is complete.

<!-- qetlab-current-claims: begin -->
## Superseding local completion evidence (2026-07-30)

The later strict static repository audit at pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601` records all 127/127 public rows as
verified with final status and all 36/36 private helpers with terminal
dispositions, with zero queued tasks and zero static-evidence failures. The
strict completion checker passes locally on Julia 1.12.6 and the installed
Julia 1.10.0.

This supersedes only the historical completion and implementation-priority
conclusions below; the dated audit findings remain as the evidence actually
observed on 2026-07-29. The new result is local static repository evidence,
not MATLAB/QETLAB parity, remote-CI evidence, a performance claim, an
API-stability guarantee, release approval, or the required human review.
<!-- qetlab-current-claims: end -->

## Executive ranking

| Priority | Category | State | Finding |
|---|---|---|---|
| P0 | Architecture | Accepted design; implementation pending | General rectangular operator spaces and two-sided map sums require the representation change specified by ADR 0006. |
| P0 | Correctness and soundness | Resolved in the current worktree | `MATLABCompat.IsPPT` previously symmetrized a slightly non-Hermitian input and could report a false `CriterionSatisfied`; it now returns `CriterionUnknown` with residual evidence. |
| P0/P1 | Correctness policy | Open | Choi-to-Kraus conversion can project an approximately Hermitian, approximately positive Choi matrix to a different CP map without returning repair metadata. |
| P1 | Resource safety | Open | Symmetric/antisymmetric projector construction and Dicke-state enumeration have integer-overflow checks but no practical combinatorial work/nonzero guards. |
| P1 | Resource safety | Open | Multiplicative and additive compound matrices have no practical output/minor/work cap. |
| P1 | Classification completeness | Partially resolved | The generic floating LDL path now scans for robustly negative diagonal witnesses and continues across exactly decoupled boundary pivots; a general pivoted method for coupled boundary pivots remains open. |
| P1 | Performance and resource planning | Open | Tensor and channel kernels create full intermediate terms without a common output/work budget. |
| P2 | Generic arithmetic | Resolved in the current worktree | Native and compatibility majorization now retain floating accumulator types while preserving widened exact arithmetic. |
| P2 | Performance | Partially resolved | `coherence_rank` now uses the validated unitary-adjoint path; pure-state negativity still takes the full density-matrix route. |
| P2 | Extension hardening | Open and documented | EntanglementDetection configuration has no seed field, child output files are not size-bounded while running, and high-level evidence fields use `Any`. |

## Architecture

### P0: rectangular operator spaces and two-sided maps

The current representation hierarchy models maps between square matrix
algebras. [`AbstractMapRepresentation`](../src/channels/channels.jl#L32-L45)
stores one input and one output Hilbert-space dimension, and
[`KrausRepresentation`](../src/channels/channels.jl#L47-L79) specifically means

```math
\Phi(X)=\sum_i K_i X K_i^\dagger.
```

Unequal input and output Hilbert dimensions are supported, but the input and
output operators themselves remain square. The current model therefore does
not describe a general map

```math
\Phi:\mathbb{C}^{m\times n}\longrightarrow\mathbb{C}^{p\times q},
\qquad
\Phi(X)=\sum_i A_i X B_i^\dagger.
```

This limitation is explicit in the compatibility layer:

- [`_compat_map`](../src/compat/MATLABCompat.jl#L857-L879) rejects paired
  left/right operator collections;
- [`ApplyMap`](../src/compat/MATLABCompat.jl#L881-L898) requires square input;
  and
- [`PartialMap`](../src/compat/MATLABCompat.jl#L970-L1024) rejects independent
  row and column subsystem dimensions.

[ADR 0006](adr/0006-general-operator-space-maps.md) is accepted and gives the
right package-wide design: a checked four-dimension operator-space descriptor,
a separately named two-sided operator-sum representation, rectangular Choi and
superoperator shapes, direct application, rectangular Hilbert--Schmidt duals,
and typed `not_applicable` physicality diagnostics. Implementation is still
pending. The affected map rows must remain partial or deferred until the native
API, compatibility mapping, provenance, analytic/property/invalid/sparse tests,
documentation, and applicable oracle evidence all agree.

This should be implemented as one representation change. Adding isolated raw
matrix overloads would duplicate dimension inference and create inconsistent
Choi, dual, partial-map, and sparsity conventions.

## Correctness and soundness

### Resolved P0: compatibility PPT Hermiticity boundary

The previous compatibility implementation accepted a nonzero Hermiticity
residual inside `TOL`, replaced the input by its Hermitian part, and could then
return `CriterionSatisfied`. PPT is not defined for that unrepaired input, so a
successful necessary-condition result was unsound.

[`MATLABCompat.IsPPT`](../src/compat/MATLABCompat.jl#L2034-L2121) now has the
following boundary:

- a Hermiticity residual above `TOL` raises;
- a nonzero residual inside `TOL` returns `CriterionUnknown`;
- the witness records `kind=:hermiticity_boundary`, the affected indices,
  defect value, and maximum residual;
- the input is neither mutated nor replaced by its Hermitian part; and
- only an exactly Hermitian input proceeds to the partial-transpose
  eigendecomposition.

The focused regression uses a one-sided `1e-10` perturbation with `TOL=1e-8`
and verifies the unknown status, evidence, and input immutability at
[`test/tier_d_measures_criteria.jl`](../test/tier_d_measures_criteria.jl#L338-L347).
The complete focused Tier D file passed 174/174 assertions, including 38/38
compatibility assertions.

### Open P0/P1: Choi-to-Kraus projection policy

[`kraus_representation`](../src/channels/channels.jl#L478-L538) converts a Choi
or superoperator representation by:

1. accepting approximate Hermiticity;
2. diagonalizing the explicitly formed Hermitian part;
3. allowing negative eigenvalues down to `-threshold`; and
4. retaining only eigenvalues strictly greater than `threshold`.

Thus small negative eigenvalues are clipped and small positive eigenvalues are
dropped. The returned `KrausRepresentation` can encode a different map, but the
return type contains no projection magnitude, discarded spectrum, or
reconstruction residual. A focused diagnostic with Choi spectrum
`[1, -1e-10, 0, 0]` and `atol=1e-8` returned a map with spectrum
`[1, 0, 0, 0]`.

The tolerance behavior is documented, so this is not a hidden implementation
accident. It remains inconsistent with the stricter package rule that a repair
must be explicit and report what changed. Recommended resolution:

- make strict rejection the default;
- require an explicit projection policy for bounded repair; and
- return a typed conversion result containing the original residuals,
  threshold, discarded eigenvalues, reconstruction error, and resulting
  representation.

Existing channel tests cover exact round trips and grossly non-CP or
non-Hermitian inputs at
[`test/tier_c_channels_maps.jl`](../test/tier_c_channels_maps.jl#L156-L205) and
[`test/tier_c_channels_maps.jl`](../test/tier_c_channels_maps.jl#L245-L272).
They do not yet lock the within-tolerance conversion boundary.

### Soundness controls already working

The native criterion path is conservative:

- [`ppt_criterion`](../src/entanglement/criteria.jl#L152-L198),
  [`realignment_criterion`](../src/entanglement/criteria.jl#L213-L253), and
  [`reduction_criterion`](../src/entanglement/criteria.jl#L287-L348) force
  structural input uncertainty to `CriterionUnknown`;
- a passed necessary condition is not described as a separability certificate;
- grossly invalid density matrices raise rather than being normalized or
  clipped; and
- optional-backend candidates remain uncertified evidence rather than being
  promoted to mathematical conclusions.

No solver-backed core routine currently collapses time limits, numerical
failure, or unavailable dependencies to a negative mathematical result. The
solver-independent optimization layer is still a proposed architecture rather
than a misleading partial implementation.

## Classification completeness

### Partially resolved P1: generic floating LDL boundary pivots

BLAS floating matrices use a Hermitian eigendecomposition, exact
integer/rational matrices use exact congruence arithmetic, and `Diagonal`
matrices have a structure-aware path. Other floating types, including
`BigFloat`, use
[`_matrix_predicate_psd_ldl`](../src/linear_algebra/matrix_predicates.jl#L281-L361).

The original method returned `MatrixPredicateUnknown` immediately when a pivot
lay inside tolerance. For example, a dense

```julia
BigFloat[0 0; 0 -1]
```

previously returned `Unknown` at the first pivot and never observed the
decisive second negative diagonal entry.

The current worktree first checks every diagonal entry for a robustly negative
basis-vector witness. It also records, rather than immediately returning at, a
boundary pivot whose remaining row and column are exactly decoupled. The LDL
scan then continues through the trailing block, and a violation takes
precedence over the recorded uncertainty. Regressions cover both a later
negative diagonal and a later indefinite `2 × 2` block; the focused native
matrix-predicate suite passes 178/178 assertions.

A boundary pivot coupled to the trailing block remains conservatively
`MatrixPredicateUnknown`; a pivoted congruence method would be needed to
improve that case without dividing by an uncertain pivot. Matrix predicates
also currently have no benchmark cases, as recorded in
[`docs/BENCHMARK_REPORT.md`](BENCHMARK_REPORT.md#L7-L12).

## Resource safety

### P1: projectors and Dicke states

The projector implementation recursively enumerates occupation sequences and
materializes every unique orbit at
[`src/subsystem/projectors.jl`](../src/subsystem/projectors.jl#L7-L64).
[`_projection_arguments`](../src/subsystem/projectors.jl#L75-L80) checks only
that `local_dimension^copies` fits `Int`. The basis constructors then collect
all orbit entries, while the projector constructors collect every pair within
each orbit at
[`src/subsystem/projectors.jl`](../src/subsystem/projectors.jl#L96-L225).

Similarly, [`dicke_state`](../src/states/states.jl#L153-L200) materializes all
`binomial(parties, excitations)` indices. The ambient power of two can fit
`Int` even when that binomial count is far beyond available memory.

Add preflight counts using `BigInt` and practical `max_columns`,
`max_nonzeros`, and `max_work` keywords. `nothing` may explicitly disable a
guard after caller review. The guarded
[`brauer_states`](../src/states/states.jl#L383-L490) API is the local model for
this policy.

Current projector tests use local dimensions and copy counts only through
three. The quick benchmark has one symmetric `d=4`, copies-four case, already
recording 335,840 bytes and 1,627 allocations.

### P1: compound matrices

[`_matrix_analysis_binomial`](../src/linear_algebra/matrix_analysis.jl#L309-L319)
only checks whether a combination count fits `Int`.
[`compound_matrix`](../src/linear_algebra/matrix_analysis.jl#L378-L427) checks
the product of its output dimensions only for integer overflow, then
materializes row and column combination tables and evaluates every minor. Its
sparse-output mode still enumerates every minor.

[`additive_compound_matrix`](../src/linear_algebra/matrix_analysis.jl#L454-L540)
likewise has no practical count cap. It materializes all combinations and a
content-keyed dictionary before enumerating exterior-action updates.

Add `max_entries`, `max_minors`, and/or `max_work` preflight guards computed
with `BigInt`. Arithmetic-overflow tests are already present, but practical
resource-limit tests are not. The quick benchmark shows why the distinction
matters:

- dense `8×8`, order-three compound: 1,738,768 bytes and 53,367 allocations;
- dense `16×16`, order-two additive compound: 427,872 bytes and 7,242
  allocations.

## Performance and output planning

### P1: tensor and channel intermediates

The following paths build full-size intermediates without a common output or
workspace budget:

- [`tensor_product`](../src/tensor_products.jl#L19-L43) and
  [`tensor_power`](../src/tensor_products.jl#L52-L60) repeatedly rebuild a
  growing `kron` result;
- [`tensor_sum`](../src/tensor_products.jl#L212-L255) materializes each full
  tensor term and stores expected shapes in an abstract `Vector{Tuple}`;
- Kraus-to-Choi and Kraus-to-superoperator conversion first materialize every
  full term at [`src/channels/channels.jl`](../src/channels/channels.jl#L406-L417);
- Kraus application first materializes every `K*X*K'` term at
  [`src/channels/channels.jl`](../src/channels/channels.jl#L594-L600); and
- [`partial_map`](../src/channels/representations.jl#L105-L156) permutes the
  whole operator and maps every block into another full output.

The benchmark report records about 16.78 MB for each dimension-16 Kraus matrix
conversion, while tensor-product coverage contains only one dense two-factor
case. First add checked result/work estimates and representative dense, sparse,
many-factor, and local-map benchmarks. Then compare direct accumulation,
preallocation, `mul!`, contractions, and sparse accumulation. Optimization
should follow measurements rather than changing algorithms speculatively.

### Resolved P2: majorization accumulator promotion

Native and compatibility majorization previously initialized scale or prefix
accumulators with `BigInt(0)`. In Julia, adding either `Float32` or `Float64`
to that value promotes the computation to `BigFloat`.

Both paths now choose an accumulator from the actual value types.
Integer-only inputs still widen to `BigInt`, rational-only inputs widen to
`Rational{BigInt}`, and homogeneous `Float32` or `Float64` inputs retain that
floating type. Focused type regressions accompany the existing exact,
tolerance, and weak-versus-strong semantic tests. The native suite passes
134/134 assertions and the compatibility suite passes 37/37.

### Partially resolved P2: structure-aware fast paths

After verifying that a supplied basis is unitary,
[`coherence_rank`](../src/coherence/coherence.jl) now computes coordinates
with `basis' * state`, avoiding the previous generic linear solve. The focused
coherence suite passes 54/54 assertions. The existing dimension-64 benchmark
should be rerun on the final tree to record the allocation effect.

The pure-state negativity path forms a full outer-product density matrix at
[`src/measures/scalar_measures.jl`](../src/measures/scalar_measures.jl#L500-L505)
and then follows dense validation, partial transpose, and eigendecomposition.
For a selected bipartition, reshaping the state and using Schmidt singular
values avoids the `O(n^2)` density allocation and the dense `O(n^3)` spectral
route. These are performance gaps, not current numerical-correctness defects.

## EntanglementDetection extension hardening

The extension boundary already isolates a known stateful backend in a fresh
child process, enforces a wall timeout, validates response schemas, caps
serialized response size, truncates returned output excerpts, and maps every
backend failure or candidate conclusion to uncertified package-owned evidence.
Tests verify caller RNG, logger, stdout, and BLAS-thread preservation.

Three limitations remain:

1. [`EntanglementDetectionSearch`](../src/entanglement/backend_interface.jl#L40-L126)
   has no seed field, and the serialized request at
   [`ext/QuantumEntanglementToolsEntanglementDetectionExt.jl`](../ext/QuantumEntanglementToolsEntanglementDetectionExt.jl#L293-L326)
   carries no seed. Child isolation protects caller state, but the experimental
   configuration does not encode or report an explicit random seed. If the
   audited backend exposes no public seed parameter, that limitation should be
   represented explicitly rather than bypassed through private APIs.
2. [`_output_excerpt`](../ext/QuantumEntanglementToolsEntanglementDetectionExt.jl#L72-L79)
   limits only how much of the completed stdout/stderr file is read. The child
   can grow those files without bound while running. Use bounded pipes/ring
   buffers or monitor file size and terminate the child. OS memory/CPU
   sandboxing remains absent and is honestly reported as
   `resource_sandboxed=false` by
   [`backend_capabilities`](../src/entanglement/backend_interface.jl#L304-L322).
3. [`EntanglementAttempt.raw_result`](../src/entanglement/backend_interface.jl#L167-L201)
   and [`EntanglementReport.evidence`](../src/entanglement/backend_interface.jl#L212-L252)
   are `Any`. They are not in a numerical hot loop, but parametric evidence
   fields would preserve extensibility while improving inference and making
   package-owned result structure more explicit.

These are hardening items, not evidence that the current adapter falsely
certifies backend output.

## Positive core findings

The audit did not find the following commonly risky patterns:

- Core randomized APIs require an explicit `rng::AbstractRNG`; no core
  `Random.seed!` or global-stream mutation was found.
- Sparse spectral operations consistently require explicit
  `allow_densify=true`, and sparse-aware subsystem kernels preserve sparse
  structure.
- No `inv(A) * b` pattern or unrestricted general matrix square root was found.
  PSD square roots use an explicit Hermitian eigendecomposition and a stated
  boundary policy.
- Subsystem permutation, partial trace, partial transpose, and realignment use
  reusable index plans instead of constructing full permutation matrices in
  their kernels. Explicit permutation-operator constructors are intentional
  public matrix constructors.
- Important kernels accept `AbstractVector`/`AbstractMatrix`; focused tests
  exercise structured matrices, sparse arrays, views, exact arithmetic, and
  `BigFloat` where supported.
- No hidden conversion of generic inputs to `ComplexF64` was found. Some
  spectral entry points deliberately restrict accepted element types instead of
  silently changing precision.
- Native entanglement criteria and the high-level pipeline use structured
  three-valued outcomes. An inconclusive calculation remains `unknown`, and a
  passed necessary test is not called a separability certificate.
- Most compatibility functions delegate to a native implementation. The few
  independent numerical bodies preserve explicitly documented compatibility
  semantics rather than silently becoming a second native API.

## Recommended closure order

1. Implement ADR 0006 as a coherent map-representation migration.
2. Decide and test the strict/projecting Choi-to-Kraus conversion contract.
3. Add practical combinatorial guards to projectors, Dicke states, and compound
   matrices.
4. Improve generic PSD completeness around boundary pivots.
5. Add output/work planning and representative allocation benchmarks before
   rewriting tensor/channel kernels.
6. Remove unintended `BigFloat` promotion and implement the two clear
   structure-aware fast paths.
7. Harden optional-backend reproducibility and child output/resource handling,
   then parameterize high-level evidence containers.

Release claims should be updated only after the corresponding implementation,
tests, documentation, provenance, and benchmark evidence are recorded together.
