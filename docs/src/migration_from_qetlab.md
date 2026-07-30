# Migration from QETLAB

<!-- qetlab-current-claims: begin -->
At pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`, the strict static completion
checker passes with 127/127 public rows verified with final status,
36/36 internal helpers assigned terminal dispositions, no queued rows, and
458 exported bindings with matching provenance entries. The direct local full
corpus passes 8,117/8,117 assertions, including 48 executable-tutorial
assertions, on Julia 1.12.6 and the installed Julia 1.10.0.

This is package-local implementation and validation evidence. It does not
establish MATLAB/QETLAB parity, remote supported-platform CI, comparative
performance, API stability, release approval, or human review.
<!-- qetlab-current-claims: end -->

The mappings below preserve the reviewed argument, convention, correction, and
row-level evidence for each completed slice. Focused and oracle counts recorded
alongside those slices remain useful evidence, but the concise block above is
the current package-wide aggregate.

## Tier A subsystem kernel

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `Tensor` | `tensor_product`, `tensor_power` | `MATLABCompat.Tensor` | Native repeated product uses `copies` or `tensor_power` | Implemented; locally tested; parity review pending |
| `TensorSum` | `tensor_sum` | `MATLABCompat.TensorSum` | Native weights are a `weights` keyword | Implemented; locally tested; parity review pending |
| `KroneckerSum` | `kronecker_sum` | `MATLABCompat.KroneckerSum` | Native repeated form uses `copies` | Implemented; locally tested; parity review pending |
| `PermuteSystems` | `permute_subsystems` | `MATLABCompat.PermuteSystems` | Native API uses dimension/permutation keywords or reusable plans; booleans replace MATLAB flags | Implemented; locally tested; parity review pending |
| `Swap` | `swap_subsystems` | `MATLABCompat.Swap` | Native API takes explicit dimensions and subsystem labels | Implemented; locally tested; parity review pending |
| `PermutationOperator` | `permutation_operator` | `MATLABCompat.PermutationOperator` | Native output is sparse by default and uses keywords for inverse/output form | Implemented; locally tested; parity review pending |
| `SwapOperator` | `swap_operator` | `MATLABCompat.SwapOperator` | Native subsystem pair and output form are keywords | Implemented; locally tested; parity review pending |
| `PartialTrace` | `partial_trace` | `MATLABCompat.PartialTrace` | Native `trace_out` keyword; implementation selected by dispatch rather than `MODE` | Implemented; locally tested; parity review pending |
| `PartialTranspose` | `partial_transpose` | `MATLABCompat.PartialTranspose` | Native `systems` keyword; rectangular form uses explicit row/column layouts or a plan | Implemented; locally tested; parity review pending |
| `Realignment` | `realign` | `MATLABCompat.Realignment` | Native generalized party selection uses `systems`; reusable plan available | Implemented; locally tested; parity review pending |
| `SymmetricProjection` | `symmetric_projector`, `symmetric_subspace_basis` | `MATLABCompat.SymmetricProjection` | Projector and orthonormal-basis outputs are separate native functions | Implemented; locally tested; parity review pending |
| `AntisymmetricProjection` | `antisymmetric_projector`, `antisymmetric_subspace_basis` | `MATLABCompat.AntisymmetricProjection` | Projector and orthonormal-basis outputs are separate native functions | Implemented; locally tested; parity review pending |

The compatibility namespace also exposes tested `BasisToLinear`,
`LinearToBasis`, and inverse-realignment spellings used by the Tier A wrappers.
Their upstream public/internal classification must follow the reviewed inventory
rather than this explanatory note.

## Tier B operators and states

The deterministic Tier B oracle contains 18 fixtures (the two Horodecki local
dimension choices are separate) and checks both API layers for 72 assertions.
It was generated with Octave 11.3.0 and is supplemental function-specific
evidence, not general MATLAB equivalence.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `Pauli` | `pauli` | `MATLABCompat.Pauli` | Native labels are validated strictly; native dense output is the default | Implemented; local and Octave fixture tests pass |
| `GenPauli` | `generalized_pauli` | `MATLABCompat.GenPauli` | Zero-based indices retained; native booleans/keywords control sparse output and precision | Implemented; local and Octave fixture tests pass |
| `GellMann` | `gell_mann` | `MATLABCompat.GellMann` | Native keyword selects sparse output and floating-point precision | Implemented; local and Octave fixture tests pass |
| `GenGellMann` | `generalized_gell_mann` | `MATLABCompat.GenGellMann` | Zero-based indices retained; native keyword selects output form/precision | Implemented; local and Octave fixture tests pass |
| `FourierMatrix` | `fourier_matrix` | `MATLABCompat.FourierMatrix` | Native `T` keyword selects real floating-point precision | Implemented; local and Octave fixture tests pass |
| `MaxEntangled` | `maximally_entangled` | `MATLABCompat.MaxEntangled` | Native normalization, sparse output, and type are keywords | Implemented; local and Octave fixture tests pass |
| `Bell` | `bell_state` | `MATLABCompat.Bell` | Native index is restricted to `0:3`; compatibility indices retain modulo-four behavior | Implemented; local and Octave fixture tests pass |
| `GHZState` | `ghz_state` | `MATLABCompat.GHZState` | Native coefficients and output form are keywords; supplied coefficients are not normalized | Implemented; local and Octave fixture tests pass |
| `WState` | `w_state` | `MATLABCompat.WState` | Native coefficients and output form are keywords; supplied coefficients are not normalized | Implemented; local and Octave fixture tests pass |
| `DickeState` | `dicke_state` | `MATLABCompat.DickeState` | Native normalization and output form are keywords | Implemented; local and Octave fixture tests pass |
| `IsotropicState` | `isotropic_state` | `MATLABCompat.IsotropicState` | Native constructor enforces the physical positivity interval and supports exact rationals | Implemented; local and Octave fixture tests pass |
| `WernerState` | `werner_state` | `MATLABCompat.WernerState` | Scalar form enforces the bipartite physical range; vector form implements the normalized lexicographic permutation sum with exact inverse-coefficient Hermiticity, PSD validation, and resource guards | Implemented; intentionally corrects the pinned multipartite loop-overwrite defect, which is retained as an Octave discrepancy fixture |
| `HorodeckiState` | `horodecki_state` | `MATLABCompat.HorodeckiState` | Local dimensions use a `dims` keyword in the native API | Implemented; both `3×3` and `2×4` local/Octave fixtures pass |
| `GisinState` | `gisin_state` | `MATLABCompat.GisinState` | Mixing probability is validated; angle is in radians | Implemented; local and Octave fixture tests pass |
| `BreuerState` | `breuer_state` | `MATLABCompat.BreuerState` | Even dimension and convex weight are validated; sparse output is a keyword | Implemented; local and Octave fixture tests pass |
| `BrauerStates` | `brauer_states` | `MATLABCompat.BrauerStates` | Returns a sparse matrix; numeric output type and pre-allocation combinatorial guards are keywords | Implemented; local guard/property and exact Octave fixture tests pass |
| `ChessboardState` | `chessboard_state` | `MATLABCompat.ChessboardState` | Optional `s`/`t` are native keywords; construction does not issue an implicit PPT verdict | Implemented; local and Octave fixture tests pass |
| `EntangledSubspace` | `entangled_subspace` | `MATLABCompat.EntangledSubspace` | Native output is the sparse diagonal-Vandermonde basis with explicit coefficient type, exact nonzero/work guards, and `r=0` support; columns are not silently normalized | Implemented; sharp dimension bounds, exact small-subspace Schmidt-rank properties, rectangular dimensions, types, guards, and Julia 1.10 are tested |

## Tier B randomized constructors

Every native and compatibility signature below adds a mandatory leading
`rng::AbstractRNG`. This deliberately changes QETLAB call syntax to guarantee
that package calls do not read or mutate Julia's process-global stream.
Randomized streams are tested with seeded properties rather than compared
across languages.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `RandomProbabilities` | `random_probabilities` | `MATLABCompat.RandomProbabilities` | Mandatory leading RNG; returns a Dirichlet-one simplex draw | Implemented; seeded/simplex/global-stream tests pass |
| `RandomStateVector` | `random_state_vector` | `MATLABCompat.RandomStateVector` | Mandatory RNG; `real` and `schmidt_rank` are native keywords | Implemented; norm/rank/reproducibility/global-stream tests pass |
| `RandomDensityMatrix` | `random_density_matrix` | `MATLABCompat.RandomDensityMatrix` | Mandatory RNG; rank and distribution are validated keywords | Implemented; PSD/trace/rank/Bures/global-stream tests pass |
| `RandomUnitary` | `random_unitary` | `MATLABCompat.RandomUnitary` | Mandatory RNG; `real=true` selects Haar orthogonal output | Implemented; unitary/orthogonal/reproducibility/global-stream tests pass |
| `RandomGraph` | `random_graph` | `MATLABCompat.RandomGraph` | Mandatory RNG; native output is a `BitMatrix` and probability is a keyword | Implemented; graph-property/reproducibility/global-stream tests pass |
| `RandomPOVM` | `random_povm` | `MATLABCompat.RandomPOVM` | Mandatory RNG; native implementation uses a Haar isometry rather than an unfinished random superoperator | Implemented; PSD/completeness/reproducibility/global-stream tests pass |
| `RandomSuperoperator` | `random_superoperator` | `MATLABCompat.RandomSuperoperator` | Mandatory RNG; native output is a structured result with Kraus certificate and bounded failure; compatibility output is a raw Choi matrix on success; unequal TP-plus-unital requests require an explicit corrected proportional-output opt-in | Implemented; 139/139 focused and 106/106 supplemental property/oracle assertions pass on Julia 1.12.6 and 1.10.11 |
| `RandomPPTState` | `random_ppt_state` | `MATLABCompat.RandomPPTState` | Mandatory RNG and dimensions; native output is a structured, independently verified result. Full-rank requests use a bounded shifted-induced construction and low-rank requests a bounded separable mixture; the wrapper returns a matrix only on verified success | Implemented with seeded PSD, trace, PPT, rank, type, resource, RNG-isolation, compatibility, and source-free property-oracle tests |

## Tier C channels and maps

Tier C uses a project-native representation layer below
`AbstractMapRepresentation`. `KrausRepresentation` records a completely
positive dilation; `OperatorSumRepresentation` records general paired
left/right factors; `ChoiRepresentation` and `SuperoperatorRepresentation`
record matrix forms; and `OperatorSpace` retains all four dimensions for maps
between rectangular matrix spaces. See
[General maps and rectangular operator spaces](general_maps.md) for the full
ownership, ordering, applicability, and densification contracts.

The representation convention is an unnormalized, input-first Choi matrix

```math
J(\Phi) = \sum_{i,j} E_{ij} \otimes \Phi(E_{ij}),
```

and Julia's column-major vectorization, so
`vec(Φ(X)) == superoperator_matrix(Φ) * vec(X)`. Choi/superoperator
conversion is an index reshuffle and preserves sparse storage. Canonical
Choi-to-Kraus recovery is defined only when complete positivity is
established. `canonical_map_decomposition` separately exposes CP,
Hermiticity-preserving, and general paired-factor branches, including the
cutoff and reconstruction residual. The physicality functions are structured
numerical diagnostics at the requested tolerances, not optimization
certificates, and do not repair or symmetrize the stored map.

The committed Tier C Octave/QETLAB artifact has 7 deterministic fixtures and
checks both native and compatibility results for 28 assertions. Its SHA-256 is
`46b31802d97ee8366da163d7c723408ac708c94325aa96c83a34307355a15c03`.
This is supplemental function-specific evidence, not general MATLAB
equivalence.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `ApplyMap` | `apply_channel` | `MATLABCompat.ApplyMap` | Typed operator sums and raw one-column, one-row CP, paired-factor, and generalized-Choi forms support independent input/output row and column dimensions | Implemented; analytic, property, sparse, exact, generic, wrapper, and prior CP-slice Octave tests pass |
| `ChoiMatrix` | `choi_matrix` | `MATLABCompat.ChoiMatrix` | Native output uses input-first ordering; the wrapper implements `SYS=1/2`, including independent row/column swaps, and retains the pinned numeric early return | Implemented; rectangular, paired, ordering, sparse, exact, and wrapper tests pass |
| `KrausOperators` | `kraus_operators`, `canonical_map_decomposition` | `MATLABCompat.KrausOperators` | Native `kraus_operators` is CP-only; the structured canonical result distinguishes CP, signed Hermitian, and general paired factors. The wrapper retains QETLAB's vector/two-column shapes and fixes its full-rectangular `pad_array` failure | Implemented with corrected boundaries; reconstruction, ordering, zero, rectangular, sparse-policy, and numerical-boundary tests pass |
| `ComplementaryMap` | `complementary_channel` | `MATLABCompat.ComplementaryMap` | Supplied dilation is preserved. Paired maps allow rectangular input and require square output; unequal output row/column dimensions are rejected as an intentional correction of the pinned failing reshape branch | Implemented with corrected boundaries; dilation, paired, representation-kind, exact, generic, sparse, and invalid-shape tests pass |
| `DualMap` | `dual_channel` | `MATLABCompat.DualMap` | Hilbert--Schmidt duality covers all typed forms and raw factor cells without changing the raw cell shape | Implemented; complex-duality, rectangular, raw-shape, and representation-kind tests pass |
| `PartialMap` | `partial_map` | `MATLABCompat.PartialMap` | Native and wrapper forms accept independent row/column subsystem layouts and apply the selected map blockwise without constructing a global tensor-product superoperator | Implemented; subsystem, dimension-changing, paired, rectangular, sparse, and invalid-layout tests pass |
| `DepolarizingChannel` | `depolarizing_channel` | `MATLABCompat.DepolarizingChannel` | Native constructor enforces the CPTP interval; the wrapper preserves QETLAB's raw Choi formula for every finite real `P` | Implemented; local and Octave fixture tests pass |
| `DephasingChannel` | `dephasing_channel` | `MATLABCompat.DephasingChannel` | Native constructor enforces the CP interval; the wrapper preserves QETLAB's raw Choi formula for every finite real `P` | Implemented; local and Octave fixture tests pass |
| `PauliChannel` | `pauli_channel` | `MATLABCompat.PauliChannel` | Native probabilities are validated without clipping/normalization; random compatibility form is `PauliChannel(rng::AbstractRNG, Q)` | Implemented; deterministic fixture plus explicit-RNG/global-stream tests pass |
| `ChoiMap` | `choi_map` | `MATLABCompat.ChoiMap` | Native output is a typed general-map representation; wrapper output is the raw Choi matrix | Implemented; local and exact Octave fixture tests pass |
| `ReductionMap` | `reduction_map` | `MATLABCompat.ReductionMap` | Native output is a typed general-map representation; wrapper output is the raw Choi matrix | Implemented; local and exact Octave fixture tests pass |
| `IsHermPreserving` | `is_hermiticity_preserving` | `MATLABCompat.IsHermPreserving` | Returns `MatrixPredicateResult`; an exactly zero Choi-Hermiticity defect is satisfied, a robust defect is violated, and a nonzero tolerance-boundary defect is unknown. Exact data are decided exactly | Implemented with corrected boundaries; CP, signed Hermitian, non-Hermitian, exact, sparse, BigFloat, raw-wrapper, and invalid-tolerance tests pass |
| `IsCP` | `is_completely_positive` | `MATLABCompat.IsCP` | Returns `MatrixPredicateResult`: a nonzero Choi-Hermiticity or PSD tolerance boundary is `unknown`, not QETLAB's Boolean positive. The wrapper accepts QETLAB's scalar tolerance and raw paired factors | Implemented with corrected boundaries; CP, violation, unknown, exact, BigFloat, sparse-policy, rectangular-applicability, raw-wrapper, and invalid-tolerance tests pass |
| `Twirl` | `twirl` | `MATLABCompat.Twirl` | Native `kind` is a strict lowercase symbol and output storage follows the input by default. The wrapper keeps case-insensitive positional `TYPE, P` and sparse output. Isotropic/Pauli calls require exactly two copies; exact local-dimension roots and explicit factorial, double-factorial, exponential, sparse, compact-dense, and work guards replace pinned rounding and uncapped pseudoinverse construction | Implemented for Werner, isotropic, real, and Pauli twirls; 149 local assertions and 18 six-fixture pinned Octave comparisons pass on Julia 1.12.6 and 1.10.11; MATLAB not run |

“Rectangular channel” and “rectangular operator space” remain different. A
channel between unequal Hilbert dimensions maps, for example, `2×2` matrices
to `3×3` matrices and may be completely positive. A general rectangular map
may instead map `M_{m,n}` to `M_{p,q}` with four independent dimensions;
complete positivity is then non-applicable, while application, Choi and
superoperator conversion, duality, decomposition, and partial mapping remain
defined.

The twirl implementation projects onto sparse permutation, maximally
entangled, Brauer, or Bell-projector spans without densifying its input.
Dependent Werner and real spanning families are reduced using exact
integer-Gram elimination before the coefficient solve, so integer and rational
inputs retain exact `Rational{BigInt}` arithmetic. No trace normalization,
symmetrization, positivity clipping, or pseudoinverse tolerance is applied.
See [Guarded group twirls](twirls.md) for definitions, complexity, and resource
contracts.

## Tier D measures and entanglement criteria

Tier D includes the scalar-measure and criterion bindings below, their reviewed
QETLAB entry-point wrappers, and a separate project-native orchestration layer.
The committed Octave 11.3.0/QETLAB artifact contains 13 fixtures, passes 34
comparisons, and has SHA-256
`ad0cdc45077390fc1eb736fc7c7ff1ec41696c796a508b536774cb6e0020160a`.
This is supplemental function-specific evidence, not general MATLAB
equivalence.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `TraceNorm` | `trace_norm` | `MATLABCompat.TraceNorm` | Sparse input requires explicit `allow_densify=true`; dependency-free implementation computes a full SVD | Implemented; local and Octave fixture tests pass |
| `SchattenNorm` | `schatten_norm` | `MATLABCompat.SchattenNorm` | Native order is validated as real `p ≥ 1`, including `Inf`; sparse densification is opt-in | Implemented; local and Octave fixture tests pass |
| `KyFanNorm` | `ky_fan_norm` | `MATLABCompat.KyFanNorm` | Native `k` is validated against `min(size(X)...)`; implementation computes the full spectrum | Implemented; local and Octave fixture tests pass |
| `Purity` | `purity` | `MATLABCompat.Purity` | Native function requires a validated density matrix; wrapper deliberately preserves QETLAB's unchecked `real(tr(RHO^2))` scalar operation | Implemented; distinct native/wrapper semantics are locally and fixture tested |
| `Entropy` | `renyi_entropy`; `von_neumann_entropy` for order one | `MATLABCompat.Entropy` | Native logarithm base is required; wrapper keeps `BASE=2`; orders `0`, `1`, finite `α ≥ 0`, and `Inf` use explicit endpoint formulas and strict density validation | Implemented for numeric density matrices; analytic endpoint, continuity, precision, sparse, invalid-input, and wrapper tests pass |
| `kpNorm` | `top_k_p_norm`; `top_k_p_norm_epigraph` for affine model data | `MATLABCompat.kpNorm` | Julia vectors use entry magnitudes and matrices use singular values; `k` is clipped to the available spectrum. The model path returns an explicit solver-neutral epigraph and never starts a nested solve | Implemented for numeric arrays and rectangular complex affine expressions; the pinned signed-vector sorting defect is deliberately corrected |
| `kpNormDual` | `top_k_p_norm_dual`; `top_k_p_norm_dual_epigraph` for affine model data | `MATLABCompat.kpNormDual` | Numeric arrays use the verified plateau formula and exact `p=1`/`Inf` endpoints. The affine path returns an exact solver-neutral Ky Fan/k-support epigraph with optional power-cone materialization instead of opening a hidden nested model | Implemented for numeric arrays and rectangular complex affine expressions; brute-force duality plus Hypatia/SCS model checks pass |
| `SkOperatorNorm` | `sk_operator_norm`; `sk_operator_norm_problem` | `MATLABCompat.SkOperatorNorm` | Requires a leading explicit RNG and returns witnessed lower and theorem/relaxation upper bounds with exactness, target, resource, hierarchy, and solver statuses kept separate. The wrapper preserves positional `K,DIM,STR,TARGET,TOL` but defaults to the structured result | Implemented with exact cases, bounded projected search, level-one PPT/reduction relaxations, and a bosonic PPT hierarchy for `k=1` |
| `IsBlockPositive` | `is_block_positive` | `MATLABCompat.IsBlockPositive` | Requires a leading explicit RNG and returns a tri-state certificate result. Legacy-shaped output maps certified true/false/unknown to `1/0/-1` and includes a validated negative witness when available | Implemented without collapsing overlapping bounds, tolerance boundaries, backend failures, or resource limits to false |
| `InducedMatrixNorm` | `induced_matrix_norm` | `MATLABCompat.InducedMatrixNorm` | Requires a leading explicit RNG and returns `InducedMatrixNormResult`; only proved closed-form branches have `bound_kind=:exact`, while every alternating branch remains a witnessed lower bound even after numerical convergence | Implemented for the complete pinned numeric contract with deterministic iteration/work limits and no global RNG mutation |
| `InducedSchattenNorm` | `induced_schatten_lower_bound` | `MATLABCompat.InducedSchattenNorm` | Requires a leading explicit RNG and returns a matrix witness in `InducedSchattenNormResult`; the `2 -> 2` transfer-matrix branch is exact and every other alternating branch is a lower bound | Implemented for the complete pinned numeric contract with deterministic iteration/work/dense-entry limits and explicit zero-gradient status |
| `Fidelity` | `fidelity` | `MATLABCompat.Fidelity` | Default is QETLAB's unsquared Uhlmann root fidelity; native `squared=true` is explicit; density inputs are strictly validated | Implemented; local and Octave fixture tests pass |
| `MatsumotoFidelity` | `matsumoto_fidelity`; `matsumoto_fidelity_model` and `matsumoto_fidelity_problem` for affine/model data | `MATLABCompat.MatsumotoFidelity` | Uses a support-aware matrix geometric mean with solves instead of explicit inverses or additive identity regularization. The symbolic path returns an explicit block-PSD hypograph rather than opening a hidden nested model | Implemented for numeric matrices and affine expressions; independent numeric, model, Hypatia, and SCS checks pass |
| `Negativity` | `negativity` | `MATLABCompat.Negativity` | Native dimensions and transposed subsystem set are explicit; wrapper retains equal-dimension inference and subsystem 2 | Implemented; local and Octave fixture tests pass |
| `SchmidtDecomposition` | `schmidt_decomposition` | `MATLABCompat.SchmidtDecomposition` | Native output is `SchmidtDecompositionResult`; wrapper returns a named tuple and preserves reviewed `K` selection rather than MATLAB output arity | Implemented; reconstruction and coefficient fixtures pass |
| `SkVectorNorm` | `schmidt_k_norm` | `MATLABCompat.SkVectorNorm` | Native subsystem dimensions are explicit; the wrapper preserves QETLAB's nearest-integer square-root default and scalar `DIM`, including a rectangular default when it divides the vector length; `k` is clipped to the smaller subsystem | Implemented; rectangular, sparse, generic full-spectrum, invalid-input, and wrapper tests pass |
| `SchmidtRank` | `schmidt_rank` | `MATLABCompat.SchmidtRank` | Native absolute/relative tolerances are explicit; result is a tolerance-defined numerical rank | Implemented; local and Octave fixture tests pass |
| `Concurrence` | `concurrence` | `MATLABCompat.Concurrence` | Domain is explicitly normalized two-qubit pure vectors or `4×4` density matrices | Implemented; local and Octave fixture tests pass |
| `IsPPT` | `ppt_criterion` | `MATLABCompat.IsPPT` | Native input is a density matrix; wrapper also accepts finite Hermitian operators. Both return structured tri-state `CriterionResult`, intentionally not QETLAB's boundary-collapsing Boolean. A nonzero Hermiticity residual inside the wrapper tolerance is `unknown` with residual evidence; the input is never symmetrized | Implemented with documented tri-state result; local and Octave fixture tests pass |
| `Distinguishability` | `state_distinguishability`; `state_discrimination_problem` | `MATLABCompat.Distinguishability` | States and priors are never normalized implicitly. Exact theorem branches and solver-derived bounds remain distinct; the general POVM SDP requires an explicit backend and retains measurement, dual, residual, and status evidence | Implemented with Helstrom/orthogonal certificates, Hypatia/SCS tests, and a source-free fixture recording the pinned unequal-prior formula defect |

The ordinary norm entry points above accept numeric arrays. The `kpNorm`
model branch uses `ComplexAffineMatrix` and `TopKPNormEpigraph`; JuMP
materialization is explicit and optional. General sparse matrix norms require
explicit `allow_densify=true`; sparse vectors remain sparse. `Diagonal`
matrices use their diagonal values directly, which preserves `BigFloat`
arithmetic without a dense factorization. General dense SVD/eigenvalue paths
retain the standard-library BLAS element-type boundary and never convert
precision implicitly.

### Induced matrix norms: exact values versus lower bounds

`induced_matrix_norm(rng, X, p; q=p, ...)` returns an
`InducedMatrixNormResult`, not a bare scalar. The maximum absolute column sum
(`1 -> 1`), largest singular value (`2 -> 2`), and maximum absolute row sum
(`Inf -> Inf`) branches are exact. The zero operator is also recognized
exactly. Every other order follows the pinned alternating Hölder-equality
updates and has `bound_kind=:lower_bound`, including
`status=:converged_lower_bound`. Numerical stationarity is not promoted to a
global optimum.

The result carries the normalized right-multiplication `witness`, its input
and output norms, normalization and value residuals, the last objective-change
residual, iteration count, deterministic work count, and termination status.
`max_iterations` counts completed alternating updates. `max_work` uses a
documented size-based accounting model for initialization, matrix-vector
products, and vector updates; neither budget consults elapsed time. A work or
iteration limit returns the best witnessed lower bound already available.
A budget too small either to complete a selected exact branch or to construct
the initial lower-bound evidence raises before calculation.

An explicit `rng::AbstractRNG` is mandatory. Exact branches and calls with a
supplied `initial_vector` do not consume it. Calls without a start draw only
from that RNG, so Julia's global random stream is untouched. The compatibility
form keeps QETLAB's positional `(P,Q,TOL,V0)` order, accepts
case-insensitive `"fro"` as order two, treats scalar `V0` as the random-start
sentinel, and vectorizes row or column starts. A malformed non-scalar start
raises rather than silently warning and switching to randomness.

Finite `Float32`, `Float64`, `ComplexF32`, and `ComplexF64` matrices are
supported without precision conversion. General alternating updates and the
exact column/row-sum branches preserve sparse matrix storage. A sparse
`2 -> 2` call requires `allow_densify=true` because it uses a full SVD.
The exact infinity-norm witness phase-aligns the maximizing row; this repairs
the pinned routine's all-ones optional witness for signed or complex rows
without changing its exact norm value.

`induced_schatten_lower_bound(rng, map, p; q=p, ...)` applies the analogous
contract to linear maps between square matrix algebras. Its `2 -> 2` value is
the exact largest singular value of the project transfer matrix. Other orders
use alternating matrix supporting functionals and remain
`bound_kind=:lower_bound`, including `status=:converged_lower_bound`.
`:stationary_zero` records that the chosen start reached a zero output or
dual gradient; it is deliberately not an exact-zero conclusion.

The result retains the matrix witness, its Schatten input/output norms,
residuals, objective history, and deterministic work accounting. General
branches use dense SVDs on input- and output-sized matrices; sparse work
matrices require `allow_densify=true`, and `max_dense_entries` checks each
conversion. The compatibility wrapper keeps QETLAB's positional
`(P,Q,DIM,TOL,X0)` order and `"fro"` alias while adding the mandatory RNG and
rejecting malformed starts instead of silently replacing them. See
[Induced Schatten norm lower bounds](induced_schatten_norm.md) for a runnable
example and complete status semantics.

The focused suite passes 135 formula, witness, residual, budget, RNG,
Float32, complex, sparse, compatibility, and invalid-input assertions on
Julia 1.12 and 1.10. A live Octave 11.3.0 check at the pinned revision with
`X = [1 2; -3 0.5; 0.25 -1]`, `p=3`, `q=2`, `TOL=1e-10`, and
`V0=[1;-2]` returned lower bound `3.207347868940717` and witness
`[0.9666227862946062; -0.4591963115719015]`; the Julia compatibility path
reproduces both to displayed precision. This is branch-specific evidence, not
a claim that a randomized lower bound equals the induced norm.

### Matsumoto fidelity and singular support

`matsumoto_fidelity(rho, sigma)` computes
`tr(rho # sigma)`, where `#` is the Kubo--Ando geometric mean. It is not an
alias for Uhlmann [`fidelity`](@ref). For positive-definite inputs, the
implementation uses Cholesky solves and a Hermitian eigendecomposition. For
singular inputs it computes the intersection of the proved spectral supports,
forms each state's shorted operator on that intersection, and takes the
positive-definite geometric mean there. No full pseudoinverse, explicit matrix
inverse, or `rho + epsilon * I` regularizer is formed.

This support definition has an important boundary consequence: distinct pure
states have zero Matsumoto fidelity even when their vector overlap is nonzero.
A small positive eigenvalue or near-coincident singular supports can make the
exact support numerically ambiguous. The default
`support_boundary_policy=:reject` reports that ambiguity.
`:project` is an explicit opt-in to discard only a bounded positive spectral
tail and canonically identify near-coincident principal support directions.
Small negative eigenvalues are always rejected, and input matrices must be
exactly Hermitian; the operation never normalizes or symmetrizes them.
General sparse input requires `allow_densify=true`. Two `Diagonal` states use
the exact commuting formula directly and retain `BigFloat` arithmetic.

The focused suite passes 62 analytic, singular-support, precision, sparse,
boundary, invalid-input, independent-limit, and wrapper assertions on Julia
1.12 and 1.10. A separate real symmetric JuMP/Hypatia check in the optional
extension environment reproduced the full-rank `2 x 2` value within
`1.0e-9`; its singular boundary solve agreed within `2.4e-5`, with the looser
agreement recorded rather than promoted to core oracle evidence. The
independent test also verifies convergence of
`(rho + epsilon I) # (sigma + epsilon I)` to the support-aware singular
answer.

Live Octave checks of the pinned numeric branch agree away from singular
boundaries. They also expose the reason that its `1e-8 I` regularizer is not
reproduced: two distinct pure states return about `1.4142e-4` upstream instead
of the exact zero, an identical pure state returns `1.000000005`, and a
reviewed singular-intersection case differs by about `7.2e-5`. The
[matrix-geometric-mean definition and SDP characterization](https://arxiv.org/abs/2006.06918)
support the exact limiting semantics used here. The pinned CVX-variable branch
is represented by package-owned affine SDP data, so the dependency-free core
can construct and inspect it without selecting a solver.

The separability and local-discrimination rows use structured native results
and safe compatibility surfaces:

| QETLAB function | Julia-native function | Compatibility name | Changed arguments and result semantics |
|---|---|---|---|
| `IsSeparable` | `is_separable` | `MATLABCompat.IsSeparable` | Native output is an `EntanglementReport` with ordered exact, sufficient, necessary, heuristic, and hierarchy attempts. The compatibility form requires a leading explicit RNG and returns `1` or `0` only for a validated certificate; `unknown` throws when `structured=false` |
| `LocalDistinguishability` | `local_distinguishability_problem`; `local_distinguishability` | `MATLABCompat.LocalDistinguishability` | Native output separates the always-achievable separable lower quantity from the numerical outer-hierarchy upper quantity and hierarchy POVM. The latter is relaxation evidence, not a certified separable or LOCC measurement. Legacy `(dist, meas, dual_sol)` output is available only after an accepted solve and checked residuals |
| `UPBSepDistinguishable` | `upb_sep_distinguishability_problem`; `upb_sep_distinguishable` | `MATLABCompat.UPBSepDistinguishable` | Complex overlaps use the conjugating Hilbert inner product. Only an exact reconstruction rigorously linked to the represented input can become `true`. Floating-SVD replacement vectors, numerical feasibility, and a separator for their floating cone remain structured `unknown`; no floating result becomes `false` |

The focused WP6 suites pass 211/211 native, 42/42 compatibility, and 72/72
optional Hypatia/SCS extension assertions on Julia 1.12.6 and the installed
Julia 1.10.0. The source-free QETLAB comparators pass 26/26 separability and
18/18 local-discrimination assertions on both Julia lines. The oracle covers
deterministic solver-free branches only; MATLAB and CVX were not run. See
[Separability and local discrimination](separability_optimization.md) for
examples and complete certificate semantics.

The Bell, nonlocal-game, and NPA rows use the same structured-result boundary:

| QETLAB function | Julia-native function | Compatibility name | Changed arguments and result semantics |
|---|---|---|---|
| `NPAHierarchy` | `npa_problem`; `npa_membership` | `MATLABCompat.NPAHierarchy` | Level zero can return an exact basic-condition verdict. Positive hierarchy levels retain the finite moment model, solver status, moment matrix, and residuals as numerical evidence and are never rounded to a Boolean membership certificate |
| `NonlocalGameLB` | `nonlocal_game_lower_bound` | `MATLABCompat.NonlocalGameLB` | Requires a leading explicit RNG plus deterministic iteration, dense-entry, model-entry, and scalar-work limits. The see-saw candidate has `certified_lower=false`, so the compatibility scalar form rejects it |
| `XORGameValue` | `xor_game_value` | `MATLABCompat.XORGameValue` | Guarded classical enumeration is exact. The quantum Gram SDP retains only an uncertified numerical upper quantity obtained from objective-bound or dual evidence |
| `BellInequalityMax` | `bell_inequality_bound` | `MATLABCompat.BellInequalityMax` | Collins--Gisin, full-correlator, and full-probability inputs are converted through validated typed representations. The guarded classical branch is exact; no-signalling and finite-NPA branches remain numerical upper quantities |
| `BellInequalityMaxQubits` | `bell_inequality_qubit_bound` | `MATLABCompat.BellInequalityMaxQubits` | Corrects the pinned rectangular-setting loop defect and exposes a solver-neutral PPT-relaxation problem. Its returned state is a relaxation variable, not a physical strategy, and its numerical upper quantity is uncertified |
| `BCSGameLB` | `bcs_game_lower_bound` | `MATLABCompat.BCSGameLB` | Uses an exact typed BCS-to-nonlocal conversion followed by the same explicit-RNG bounded see-saw. Numerical candidates are never collapsed to scalars |
| `BCSGameValue` | `bcs_game_value` | `MATLABCompat.BCSGameValue` | Supplies the `NonlocalGameValue` capability missing from the pinned source tree. Its classical branch is exact; no-signalling and finite-NPA branches retain numerical upper quantities and solver diagnostics |

See [Bell inequalities and nonlocal games](nonlocal_optimization.md) for
executable examples, model limits, and the distinction between exact values,
numerical upper relaxations, and attained numerical candidates.

The project-native scalar additions `trace_distance`,
`logarithmic_negativity`, and `schmidt_coefficients` have no fabricated
one-to-one QETLAB mapping. Likewise, `CriterionStatus`, `CriterionResult`,
`realignment_criterion`, and `reduction_criterion` are independently designed
structured APIs. A criterion violation can certify entanglement; a satisfied
necessary condition is not generally a separability certificate, and a
tolerance-boundary result is `unknown`.

The project-native orchestration exports
`AbstractEntanglementMethod`, `AbstractEntanglementBackend`,
`NativeEntanglementBackend`, `NativePPT`, `EntanglementAttempt`,
`EntanglementReport`, `detect_entanglement`, `analyze_entanglement`,
`is_separable`,
`backend_capabilities`, and `available_entanglement_backends`. The
certificate-first matrix path tries PPT, realignment/CCNR, and reduction
criteria. It certifies separability from PPT only in the exact bipartite
`2×2`/`2×3` domain. For pure vectors, a trailing Schmidt coefficient above the
configured threshold certifies entanglement, but separability is certified only
when the computed trailing coefficients are exactly zero; a tolerance-defined
rank-one result with nonzero trailing coefficients is `unknown`. The composite
`is_separable` API adds the pinned behavioral strategy family while preserving
these certificate boundaries and explicit resource limits.

## Channel norms and optimization

The status-aware channel-optimization slice passes 124 focused native
assertions, 21 compatibility assertions, 74 optional Hypatia/SCS assertions,
and 122 source-free oracle assertions on Julia 1.12 and Julia 1.10. Its
committed Octave 11.3.0/QETLAB fixture has SHA-256
`4855dae732c48b4e8e31e357e6ae1fcf392d62572bf6e6092d278ed30fe327c7`.
MATLAB and CVX were not run; general SDP evidence comes independently from
the Julia optional-solver environment.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `DiamondNorm` | `diamond_norm`; `diamond_norm_problem` | `MATLABCompat.DiamondNorm` | Native maps carry their operator-space dimensions; raw compatibility inputs retain positional `DIM`. Analytic certificates, numerical lower/upper bounds, missing backends, limits, and failures are distinct, and sparse densification is explicit | Implemented for square input/output algebras, unequal dimensions, and one- or two-sided map representations with a solver-neutral Watrous primal |
| `CBNorm` | `cb_norm` | `MATLABCompat.CBNorm` | The result records the Hilbert--Schmidt-adjoint reduction and retains the full diamond-norm evidence; a compatibility scalar is returned only after a conclusive result | Implemented with analytic and Hypatia/SCS branches |
| `ChannelDistinguishability` | `channel_distinguishability` | `MATLABCompat.ChannelDistinguishability` | Priors are finite, nonnegative, sum to one, and are never normalized. The native and compatibility routes apply the full channel Holevo--Helstrom conversion `(p₁+p₂+norm)/2`, correcting the pinned routine's non-deterministic probability defect | Implemented with exact identical/deterministic certificates and status-rich weighted diamond-norm bounds |
| `MaximumOutputFidelity` | `maximum_output_fidelity`; `maximum_output_fidelity_problem` | `MATLABCompat.MaximumOutputFidelity` | The convention is root fidelity. A direct SDP is invariant under Kraus representation and does not reproduce the pinned unequal-Kraus-rank truncation; identical, common-output, and replacer cases have analytic certificates | Implemented with direct solver-neutral modeling and optional Hypatia/SCS evidence |

The compatibility wrappers return structured results by default. Set
`structured=false` only when a QETLAB-shaped scalar is required; an
inconclusive result then raises instead of discarding its bounds or status.
See [Channel norms, discrimination, and output fidelity](channel_optimization.md).

## Tier E coherence slice

The focused coherence suite passes 54 local assertions. Its committed
Octave 11.3.0/QETLAB artifact contains six deterministic fixtures and passes
25 native/compatibility/discrepancy assertions, with SHA-256
`11bcaaee88fac8a595e9a4eff164432dbaa4e141cdd26554da2522e22981811a`.
MATLAB was not run.

The separate pure-state robustness suite passes 82 focused assertions on Julia
1.12.6 and Julia 1.10.11. Its source-free Octave artifact adds 58 assertions
over eight theorem-domain agreement fixtures and three reviewed correction
cases, with SHA-256
`551f529a86960a43e8765e2934007a1a640b86f8327a2f19f79860d9d5ce51a6`.

The status-aware coherence-optimization slice passes 386 shared core
assertions, 88 optional Hypatia/SCS assertions, and 51 source-free oracle
assertions on Julia 1.12 and 1.10. Its fixture SHA-256 is
`2b44fc70b33183c01c8c8dad4c71e122ed182fa776b6e95cc15ac693eb499dad`.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `L1NormCoherence` | `l1_coherence` | `MATLABCompat.L1NormCoherence` | Native pure vectors and density matrices are strictly validated; sparse density validation requires explicit densification permission | Implemented; local and Octave fixture tests pass |
| `RelEntCoherence` | `relative_entropy_coherence` | `MATLABCompat.RelEntCoherence` | Native logarithm base is required; wrapper keeps QETLAB's base-two default | Implemented; local and Octave fixture tests pass |
| `CoherenceRank` | `coherence_rank` | `MATLABCompat.CoherenceRank` | Native tolerances and optional basis are keywords. Both Julia paths count nonzero coefficients as documented rather than reproducing the pinned implementation's zero-counting bug | Implemented with documented upstream bug fix and committed discrepancy fixtures |
| `RobkCohValue` | `pure_k_coherence_robustness` | `MATLABCompat.RobkCohValue` | Native output is a structured result with the theorem branch and adjacent stability gaps; compatibility returns `(robustness, branch_index)`. Both sort coefficient magnitudes, support complex phases, require normalization, and enforce $2 \leq k \leq n$ | Implemented from Theorem 1; eight pinned-domain fixtures agree and three unsafe upstream behaviors are retained as correction evidence |
| `IskIncoherent` | `is_k_incoherent` | `MATLABCompat.IskIncoherent` | The native result distinguishes exact/sufficient theorems, bounded band search, factor-width SDP evidence, limits, boundaries, and unknown. The wrapper returns `1`/`0` only when conclusive | Implemented with a corrected bounded graph-layout search and Hypatia/SCS factor-width tests |
| `IsAbskIncoh` | `is_absolutely_k_incoherent` | `MATLABCompat.IsAbskIncoh` | Spectral theorem branches retain one-sided versus necessary-and-sufficient semantics; the $k=d-1$ model is explicit and status-aware | Implemented with strict state validation, boundary/unknown results, and optional solver evidence |
| `RobustnessCoherence` | `robustness_coherence` | `MATLABCompat.RobustnessCoherence` | Exact pure/qubit/diagonal branches and the general SDP return free/noise state evidence; scalar compatibility output is available only when a value exists | Implemented with explicit backend, model limits, residuals, and structured failure |
| `TraceDistanceCoherence` | `trace_distance_coherence` | `MATLABCompat.TraceDistanceCoherence` | The closest state is consistently a density matrix natively and a diagonal vector in legacy-shaped output, correcting the pinned qubit inconsistency | Implemented for exact pure/qubit and general status-aware SDP branches |
| `GenRobustnesskCoherence` | `generalized_robustness_k_coherence` | `MATLABCompat.GenRobustnesskCoherence` | The absent pinned `IskCoherent` dependency is replaced by the factor-width cone. Zero robustness has `noise_state=nothing` rather than division by zero | Implemented with exact pure-state and general SDP branches |

Pure-vector paths preserve generic floating-point precision and sparse storage.
Matrix entropy/positivity work has the same explicit BLAS-type and
`allow_densify` policy as the Tier D spectral routines. See
[Coherence](coherence.md) and
[Pure-state robustness of k-coherence](pure_k_coherence_robustness.md), and
[Coherence criteria and optimization](coherence_optimization.md).

The `RobkCohValue` theorem is permutation- and phase-invariant, but the pinned
loop applies its comparisons directly to the supplied entries. Its unsorted
fixture selects a different branch and value, its complex fixture returns a
complex quantity, and its nonnormalized fixture returns a scale-dependent
number. The Julia-native and compatibility paths intentionally correct all
three cases rather than treating them as parity targets. Equality uses the
theorem's non-strict comparison and selects the largest admissible branch;
tolerances report branch sensitivity but never move the input between
branches.

## Tier E product analysis and separable-ball slice

The core product-analysis slice adds 10 native bindings and six reviewed
QETLAB entry-point wrappers. Its native focused suite passes 197 assertions
and its compatibility suite passes 61. The committed Octave 11.3.0/QETLAB artifact contains 14
deterministic fixtures and passes 68 native, wrapper, reconstruction, and
discrepancy assertions, with SHA-256
`ab6414c1a684141db74782616d4c18e79c8e6039aad53a695d8c723eed598d85`.
MATLAB was not run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `OperatorSchmidtDecomposition` | `operator_schmidt_decomposition`, `operator_schmidt_coefficients` | `MATLABCompat.OperatorSchmidtDecomposition` | Native result is typed, supports separate rectangular row/column layouts and a requested real Hermitian-basis convention, and returns the full thin decomposition. The wrapper returns a positionally destructurable three-field named tuple; it preserves `K=-1`, `K=0`, and positive-`K` selection and automatically uses Hermitian factors for exactly Hermitian input with locally square dimensions | Implemented for finite BLAS floating matrices with documented fixes for pinned rectangular reshapes, unequal-local-dimension Hermitian basis indexing, and repaired-branch `K` handling; local invariant tests and supported Octave coefficient fixtures pass, MATLAB not run |
| `OperatorSchmidtRank` | `operator_schmidt_rank` | `MATLABCompat.OperatorSchmidtRank` | Native rank uses explicit absolute/relative tolerances; wrapper preserves QETLAB's reviewed default absolute threshold | Implemented; local and Octave fixture tests pass |
| `IsProductVector` | `is_product_vector` | `MATLABCompat.IsProductVector` | Returns `ProductAnalysisResult` with factors, scaled cut residuals, reconstruction residual, thresholds, and a three-way status instead of collapsing a tolerance boundary to a Boolean | Implemented with intentional structured result; local and Octave fixture tests pass away from boundaries |
| `IsProductOperator` | `is_product_operator` | `MATLABCompat.IsProductOperator` | Supports multipartite rectangular local operators and returns the same structured numerical analysis instead of a Boolean | Implemented with intentional structured result; local and Octave fixture tests pass away from boundaries |
| `IsEntanglingGate` | `is_entangling_gate` | `MATLABCompat.IsEntanglingGate` | Returns `EntanglingGateResult` with either a local-factor/permutation certificate, a validated product-state witness, or `:unknown`; factorial permutations, witness candidates, and work are guarded | Implemented with a corrected finite phase-grid witness search; pinned identity/swap/CNOT/controlled-Z flags agree, while the invalid controlled-Z upstream witness fall-through is recorded rather than reproduced |
| `IsUPB` | `is_upb` | `MATLABCompat.IsUPB` | Native output is a `UPBAnalysisResult` that checks the complete definition, retains exact or tolerance-aware exhaustion evidence, and returns compact product-extension factors. Boolean compatibility output is available only for conclusive results | Implemented with lazy guarded partition search and Hermitian complex witnesses; 146 focused assertions and 33 committed-oracle assertions pass, while pinned nonorthogonal/complete-basis false positives and a complex-witness defect are recorded rather than reproduced |
| `MinUPBSize` | `minimum_upb_size` | `MATLABCompat.MinUPBSize` | Native output is an exact-or-unknown `MinimumUPBSizeResult` with the counting lower bound and primary reference; compatibility returns the known integer and carries the structured result in the unresolved-case error | Implemented for every pinned theorem-table branch; 11 exact families and one unknown route match the pinned Octave fixture |
| `UPB` | `upb` | `MATLABCompat.UPB` | Native output is a stable `UPBConstruction` containing owned local factors, global vectors, family/reference metadata, validation residuals, and resource/RNG diagnostics. Numeric dispatch preserves caller subsystem order. The compatibility `output` keyword replaces MATLAB `nargout`; randomized routes require a leading explicit RNG | Every executable named and dimension-driven catalog family is implemented with guarded allocations and searches. The invalid pinned `GenTiles1(2)` singleton and nonorthogonal `John2^4k` reshape branch are corrected and recorded as upstream discrepancies |
| `EntFormation` | `entanglement_of_formation` | `MATLABCompat.EntFormation` | Covers pure vectors and rank-one density matrices in arbitrary bipartite dimensions plus mixed two-qubit states; wrapper preserves the rounded-square-root dimension default, while all bounded PSD/rank/range projections require explicit opt-in | Implemented with intentional safety and zero-limit corrections; pinned zero-concurrence `NaN` is returned as mathematical zero |
| `InSeparableBall` | `in_separable_ball` | `MATLABCompat.InSeparableBall` | Native input must already be normalized and the result distinguishes a certificate, outside-ball failure, and numerical unknown; wrapper retains safe positive-trace normalization | Implemented with structured sufficient-certificate semantics; local and Octave fixture tests pass |

An additional live Octave 11.3.0 check against the pinned `EntFormation.m`
confirmed its arbitrary-dimensional rank-one matrix conversion: a `2 x 3`
maximally entangled projector returned `1`, and a nonuniform `2 x 3` pure
state returned the same value (`0.8812908992306926`) in vector and projector
forms. The same run reconfirmed the reviewed upstream `NaN` at zero
concurrence.

`ProductAnalysisResult` is a numerical classification, not an exact symbolic
proof. `SeparableBallResult(:outside_ball, ...)` is a failed sufficient test,
not an entanglement conclusion. General matrix spectral work requires explicit
sparse densification; diagonal matrices and supplied eigenvalues have
structure-aware paths. See
[Product structure and separable-ball certificates](product_analysis.md).

The separate six-fixture `IsEntanglingGate` artifact passes 13 assertions
against Octave 11.3.0 and the pinned QETLAB revision. It records that the
pinned controlled-Z call returns an unnormalized product-preserving final
candidate despite its true entangling flag; the Julia result instead carries
a validated four-support witness. MATLAB was not run.

The `IsUPB` artifact covers the Tiles and Shifts UPBs, an extendible set, a
complete product basis, a nonorthogonal unextendible set, and a complex
extension witness. It passes 33 native, compatibility, agreement, and
discrepancy assertions. The native API enforces the full UPB definition and
Hermitian orthogonality, deliberately correcting the pinned false positives
and nonconjugating witness calculation. The committed artifact SHA-256 is
`56f2f9bd0a3fa7d18ed31f45f896926f5e2c6df3b2400b59b0295f4fae0119e1`;
MATLAB was not run. See
[UPB certificates and extension witnesses](is_upb.md).

The `MinUPBSize` artifact covers 11 exact theorem/exception families and the
pinned `(2,3,4)` unknown error route. The native result retains the counting
lower bound when the exact size is unresolved. This evidence was generated by
Octave 11.3.0 from the pinned revision; MATLAB was not run.

The separate `UPB` catalog artifact compares local factors up to independent
phase for 17 deterministic named/dimension routes and passes 137 parity and
invariant assertions. Eight additional assertions retain the pinned
`John2^4k(8)` reshape-order failure and validate the corrected native
construction. Randomized Alon–Lovász and Chen–Johnston routes are covered by
seeded theorem, full-spark, orthogonality, unextendibility, and RNG-isolation
tests rather than unrelated entrywise random draws. The committed artifact
SHA-256 is
`cac29476f30f3a901e7c7ae36877f54097ffe1c6f91d049180b2bb1113585a93`;
MATLAB was not run. See
[UPB construction catalog](upb_catalog.md).

## Absolute-PPT spectral criteria

The absolute-PPT slice ports the finite Hildebrand spectral LMI construction
without returning CVX or JuMP-owned expressions. Numeric and affine spectra
produce package-owned matrices and ordering metadata; optional JuMP backends
are used only to propose a product-order realization after the deterministic
bounded search ends. Every solver point is rechecked in exact rational
arithmetic before it can support a negative certificate.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `AbsPPTConstraints` | `abs_ppt_constraints` | `MATLABCompat.AbsPPTConstraints` | Native output is an `AbsPPTConstraintFamily` with explicit criss-cross orderings, exhaustion state, counts, work/storage limits, and numeric or affine Hermitian matrices. Native scalar `dims=d` consistently infers the second factor; the wrapper retains the pinned `(d,d)` interpretation. An unexpected resource-capped plain compatibility result raises instead of masquerading as the full family | Implemented; local numeric, exact, affine, sparse, ownership, p=1…6 count, limit, compatibility, and source-free pinned-Octave fixture tests pass |
| `IsAbsPPT` | `is_abs_ppt` | `MATLABCompat.IsAbsPPT` | Native output distinguishes analytic, sufficient, exhaustive, certified-negative, capped, numerical-boundary, unavailable-backend, and failed-backend outcomes. The wrapper maps certified yes/no/inconclusive to `1/0/-1` and preserves the pinned scalar/default rectangular dimension inference | Implemented; local unitary-orbit, certificate, boundary, resource, malformed-backend, Hypatia, SCS, compatibility, and pinned-Octave tests pass |

The $p=6$ QETLAB criss-cross family contains 2612 matrices, including four
documented redundant orderings beyond the 2608 realizable orderings. A robust
negative matrix therefore becomes `false` only with an exact ordering
realization. Conversely, passing a capped family is always inconclusive.
Complex eigenvalue vectors are rejected, general sparse spectral work requires
explicit densification, and no input is normalized, clipped, symmetrized, or
trace-repaired. See [Absolute PPT from a spectrum](absolute_ppt.md) for the
criterion, examples, statuses, and solver-independent affine construction.

## Symmetric-extension hierarchies

The outer and inner hierarchies use package-owned affine SDP models and require
an explicit optional backend for nonanalytic instances. Solver termination is
never itself a mathematical Boolean: reconstructed primal extensions and dual
separators must pass independent residual checks.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `SymmetricExtension` | `symmetric_extension` | `MATLABCompat.SymmetricExtension` | Native output is a structured result with explicit theorem, backend, limit, boundary, primal, and dual evidence. `symmetric_extension_problem` exposes the solver-neutral SDP. The compatibility wrapper returns `1` or `0` only for conclusive results and otherwise raises unless `structured=true` | Implemented with analytic order-one, two-qubit, low-dimensional PPT, universal NPT, full and bosonic SDP, validated primal/dual, compatibility, Hypatia/SCS, and source-free oracle tests |
| `SymmetricInnerExtension` | `symmetric_inner_extension` | `MATLABCompat.SymmetricInnerExtension` | Native output preserves the NOP mixing parameter and explicitly labels a negative dual as a separator from the inner cone, not automatically an entanglement witness. `symmetric_inner_extension_problem` exposes the solver-neutral SDP | Implemented with non-PPT and PPT formulations, exact private Jacobi recurrence, validated primal/dual, compatibility, Hypatia, and source-free oracle tests |

Tensor factors are ordered `A, B₁, ..., Bₖ`. Inputs are not normalized,
symmetrized, clipped, or repaired. Sparse spectral work and the portable
complex real-block SDP embedding require explicit densification. See
[Symmetric-extension hierarchies](symmetric_extensions.md) for the precise
certificate and inconclusive-status contract.

## Bounded iterative operator scaling

The focused `OperatorSinkhorn` suite passes 124 assertions on Julia 1.12.6
and Julia 1.10.11. These are independent reconstruction, marginal, status,
conditioning, sparse-policy, type, and resource-bound tests; no
function-specific MATLAB-family oracle has been recorded.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `OperatorSinkhorn` | `operator_sinkhorn` | `MATLABCompat.OperatorSinkhorn` | Native output is an `OperatorSinkhornResult` with multipartite filters, first-two left/right aliases, residual history, convergence status, conditioning diagnostics, failure subsystem, and work accounting. The compatibility result has exactly two positionally destructurable fields `(sigma, filters)` and is returned only after checked convergence | Implemented for finite BLAS floating density operators with explicit sparse densification, iteration, conditioning, allocation, and work limits. The implementation intentionally replaces pinned `inv`/`sqrtm`, unbounded iteration, global warning mutation, catch-all failure, final symmetrization, and trace repair with factorizations and auditable statuses; MATLAB not run |

The native algorithm uses a documented trace-one work gauge while preserving
the original trace and direct local-filter reconstruction of the returned
operator. Inputs are never normalized in place, symmetrized, or spectrally
clipped. Singular and ill-conditioned marginals are controlled nonconverged
results rather than Boolean or negative mathematical conclusions. See
[Bounded operator Sinkhorn scaling](operator_sinkhorn.md).

## Bounded filter normal form

The focused `FilterNormalForm` suite passes 148 assertions on Julia 1.12.6
and Julia 1.10.11. It covers both reconstruction identities, equal and
rectangular dimensions, a full-rank complex input, a rank-one state whose
balanced marginals still admit a form, every native failure status, all
reviewed compatibility dimension forms, the corrected `TOL` route,
precision and sparse policy, and entry/work guards. No function-specific
MATLAB-family oracle has been recorded.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `FilterNormalForm` | `filter_normal_form` | `MATLABCompat.FilterNormalForm` | Native output is a `FilterNormalFormResult` containing the complete thin coefficient vector, Hermitian orthonormal local operators, filters, numerical input/coefficient ranks, nested Sinkhorn diagnostics, two reconstruction residuals, and total work accounting. Compatibility returns exactly five positionally destructurable fields `(xi, GA, GB, FA, FB)` and truncates only at the reviewed QETLAB operator-Schmidt threshold | Implemented for finite BLAS floating bipartite positive operators with explicit sparse densification, iteration, conditioning, allocation, and work limits. `TOL` is intentionally forwarded to Sinkhorn, correcting the pinned routine's parsed-but-unused argument. Distinct singular, ill-conditioned, iteration, work, and numerical statuses replace the pinned catch-and-rethrow collapse; MATLAB not run |

The native identity includes the filtered operator's actual trace, so it is
valid beyond trace-one density matrices. Rank deficiency is reported but is
not treated as an automatic failure: a low-rank state with nonsingular
balanced marginals may have a valid form. No input or output is silently
normalized, symmetrized, clipped, or trace-repaired. A derived
operator-Schmidt work matrix is averaged with its adjoint only inside a
reported roundoff boundary. See
[Bounded filter normal form](filter_normal_form.md).

## Tier E matrix analysis

The native matrix-analysis suite passes 131 local assertions, and the separate
compatibility suite passes 34. The committed Octave 11.3.0/QETLAB artifact
contains 22 deterministic fixtures: 17 agreement fixtures and five reviewed
semantic discrepancies. Its 59 assertions pass with SHA-256
`e37685c262ce5982d10dd705cef8c172d49d9c55c89a0a67d4de729a5068f540`.
MATLAB was not run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `Majorizes` | `majorizes` | `MATLABCompat.Majorizes` | Native uses standard strong majorization, equal totals, symmetric tolerances, padding before sorting, and singular values for every matrix shape. The compatibility entry point preserves pinned weak totals, row/column vectors, sort-before-padding, and its one-sided tolerance while exposing safe tolerance and sparse-densification keywords | Implemented; native and compatibility contracts plus all reviewed discrepancies are locally tested |
| `ElemSymPoly` | `elementary_symmetric_polynomial` | `MATLABCompat.ElemSymPoly` | Native name is descriptive Julia `snake_case`; row/column matrices are accepted by the wrapper. The dynamic program preserves sparse vectors and checked exact arithmetic rather than enumerating combinations in floating arithmetic | Implemented; local and Octave fixture tests pass |
| `CompoundMatrix` | `compound_matrix` | `MATLABCompat.CompoundMatrix` | Native retains `binomial(m,k) × binomial(n,k)` shapes when one dimension is zero; the wrapper preserves pinned `0×0` for `k > min(m,n)`. Exact minor arithmetic and sparse output selection are retained where values agree | Implemented; local and reviewed shape-discrepancy tests pass |
| `AdditiveCompoundMatrix` | `additive_compound_matrix` | `MATLABCompat.AdditiveCompoundMatrix` | Native defines order zero as the `1×1` additive zero; the wrapper reproduces the reviewed pinned dependency-path error. Positive orders use checked exact arithmetic and explicit sparse output | Implemented; local and reviewed order-zero-discrepancy tests pass |
| `Commutant` | `commutant` | `MATLABCompat.Commutant` | A matrix or nonempty Julia tuple/vector replaces a MATLAB matrix or cell. The native result is a dense Hilbert--Schmidt-orthonormal basis with explicit rank tolerances, densification permission, and allocation/work budgets. The wrapper preserves sparse output storage when all inputs were sparse | Implemented for finite numeric matrices supported by the standard-library dense SVD; subspace, residual, type, sparse-policy, and guard tests pass; MATLAB not run |

The five discrepancy fixtures are evidence, not native expected values: weak
totals, negative-vector padding order, row-matrix vector semantics,
rectangular high-order compound shape, and the additive order-zero error. See
[Matrix analysis](matrix_analysis.md) for definitions and complexity.

The `Commutant` basis is nonunique, so its focused tests compare null-space
dimensions, Hilbert--Schmidt orthonormality, commutator residuals, and subspace
projectors under a unitary basis change. Sparse input still requires explicit
permission for the dense SVD. The compatibility wrapper's sparse basis storage
matches the pinned contract, but its individual sparse matrices may contain
many nonzeros. The source-informed null-space policy retains Bruno Luong's
separate bundled-`spnull` attribution; see the repository-level
`docs/LEGAL.md`.

## Solver-free polynomial foundations

The polynomial layer fixes one public coefficient-order contract and keeps
numeric construction separate from symbolic modeling or solver integration.
Its focused suite passes 181 assertions on Julia 1.12.6 and Julia 1.10.11.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `CopositivePolynomial` | `copositive_polynomial` | `MATLABCompat.CopositivePolynomial` | Native returns an owned `HomogeneousPolynomial`; the wrapper returns QETLAB-ordered coefficients and preserves sparse storage unless `dense_output=true`. Both reject nonsymmetric input instead of silently applying `(C+C')/2` | Implemented for finite real symmetric numeric matrices with exact-type and combinatorial-overflow tests |
| `PolynomialAsMatrix` | `polynomial_as_matrix` | `MATLABCompat.PolynomialAsMatrix` | The fixed `:qetlab_lexicographic` coefficient order is explicit. Sparse output is the default, densification is opt-in and guarded, and the compact normalized monomial basis is documented | Implemented for finite numeric coefficients; CVX and other symbolic expressions are intentionally outside this solver-free API |
| `PolynomialOptimize` | `polynomial_bounds` | `MATLABCompat.PolynomialOptimize` | An explicit `rng::AbstractRNG` and exact `inner_samples` replace global RNG and wall-clock sampling. The structured result distinguishes hierarchy outer bounds from sampled feasible values. Maximization handles `target` directly, correcting the pinned target-sign recursion defect | Implemented for finite real even-degree numeric polynomials using the standard-library generalized eigensolver with explicit densification and work limits; symbolic and solver-backed branches are not claimed |
| `PolynomialSOS` | `polynomial_sos_problem`; `polynomial_sos_bounds` | `MATLABCompat.PolynomialSOS` | The primal symmetric-moment SDP is solver-neutral; an explicit backend supplies outer relaxation evidence. A mandatory RNG and exact sample count replace global, wall-clock-derived sampling, and the result keeps outer and attained inner bounds separate | Implemented with missing-backend, target, RNG-isolation, resource, Hypatia, SCS, limit, and failure tests on Julia 1.10 and 1.12 |
| `IsCopositive` | `copositivity_criterion` | `MATLABCompat.IsCopositive` | Requires a leading explicit RNG. Exact sufficient certificates and exact rational negative witnesses produce true/false; floating hierarchy values, boundaries, backend failures, and limits remain inconclusive. Invalid modes are rejected | Implemented with lazy exact branches, owned read-only witnesses, bounded SOS/solver-free routes, and checked `structured=false` compatibility output |
| `CliqueNumber` | `clique_number_bounds` | `MATLABCompat.CliqueNumber` | Requires a leading explicit RNG and returns graph-certified integer bounds separately from floating hierarchy candidates. Legacy `(ub,lb)` output uses only certified edge, degree, coloring, clique, or exact Motzkin--Straus evidence | Implemented with exhaustive graphs through five vertices, guarded exact large-graph branches, and inspectable read-only certificate payloads |

The compatibility optimizer therefore also requires a leading RNG and returns
a `PolynomialOptimizationResult`, rather than collapsing bounds, provenance,
sampling evidence, and target status into two positional scalars. See
[Solver-free polynomial foundations](polynomial_foundations.md) for the
ordering identity, runnable examples, and all resource guards.

## Nonlocal-game tensor construction

| QETLAB function | Julia-native API | Compatibility API | Status and migration notes |
|---|---|---|---|
| `ParallelRepetition` | `parallel_repetition(game, repetitions; max_entries, max_work)` | `MATLABCompat.ParallelRepetition(V, REPT; max_entries, max_work)` | Implemented for finite four-dimensional numeric coefficient tensors. Packed copy order agrees with QETLAB's tensor convention, exact coefficient types are preserved, and output allocation/work are checked before execution. Nonpositive repetitions are rejected instead of reproducing the pinned routine's accidental unchanged-input behavior. |

This row constructs a repeated coefficient tensor only. It does not compute a
classical, no-signalling, or quantum game value. See
[Nonlocal games](nonlocal_games.md) for the axis convention, a runnable
two-copy example, and exponential resource limits.

## Tier E matrix predicates

The native predicate suite passes 170 assertions and its compatibility suite
passes 37. No MATLAB predicate oracle has been run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `IsPSD` | `is_positive_semidefinite`; `positive_semidefinite_constraint` for affine data | `MATLABCompat.IsPSD` | Numeric input returns a three-valued `MatrixPredicateResult` without symmetrization. Affine input returns owned solver-neutral cone data and can be materialized explicitly through the JuMP extension | Implemented for numeric and affine inputs with copy/nonmutation, unavailable-backend, Hypatia-status, limit, and failure tests |
| `IsLocallyPSD` | `is_locally_positive_semidefinite` | `MATLABCompat.IsLocallyPSD` | Preserves principal-index witnesses and guards `binomial(n,k)` work | Implemented; analytic, exact, boundary, sparse, and guard tests pass |
| `IsTotallyPositive` | `is_totally_positive` | `MATLABCompat.IsTotallyPositive` | Strictly singular minors violate the mathematical predicate; nonzero boundary determinants are `unknown`, unlike QETLAB's permissive Boolean comparison | Implemented with documented strict/boundary divergence |
| `IsTotallyNonsingular` | `is_totally_nonsingular` | `MATLABCompat.IsTotallyNonsingular` | Exact inputs are exact; represented singularity is distinguished from a nonzero boundary determinant; all-minor work is guarded | Implemented with intentional structured result |

The compatibility all-minor entry points keep QETLAB's positional
`SUB_SIZES, TOL` order and single default determinant tolerance. They never
coerce `unknown` to `true` or `false`. See
[Matrix predicates](matrix_predicates.md).

## Migration policy

- Julia-native functions use lowercase `snake_case`, ordinary Julia arrays,
  multiple dispatch, and keywords.
- QETLAB spellings belong in a separate `MATLABCompat` namespace and delegate to
  tested native methods where the contracts agree. Direct compatibility entry
  points such as unchecked `Purity` and operator-valued tri-state `IsPPT`
  document why strict native delegation would change the reviewed contract.
- A wrapper documents any ordering, shape, tolerance, output-arity, or default
  difference. It does not silently guess an ambiguous MATLAB call.
- Random compatibility wrappers require a leading explicit RNG even though the
  original MATLAB entry points draw from MATLAB's implicit global stream. This
  includes the scalar random `PauliChannel(rng, Q)` compatibility form.
- Deprecated upstream functions receive an alias or an explicit documented
  replacement where useful.
- An operation becomes `compatibility_alias` only after both native and wrapper
  paths are tested.

MATLAB and Julia are both column-major, but that fact alone does not establish
equivalent reshape, permutation, adjoint, or subsystem-label semantics. See
[Mathematical conventions](conventions.md).
