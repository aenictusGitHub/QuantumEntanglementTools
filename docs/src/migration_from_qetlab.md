# Migration from QETLAB

The mappings below are covered by focused local suites: 304 Tier A assertions,
941 Tier B assertions, 154 Tier C assertions, 162 Tier D measures/criteria
assertions, and 68 project-native entanglement-pipeline assertions. The
integrated corpus passes 2,295 assertions (2,270 core plus 25 executable
tutorials) on Julia 1.12.6 and 1.10.11. Reviewed
inventory/provenance rows exist for these slices, including
explicit partial statuses where the Tier C compatibility surface does not yet
cover the full pinned QETLAB entry point and where Tier D `Entropy` implements
only `ALPHA=1`. MATLAB differential validation and supported-platform CI remain
incomplete, so “locally tested” is not a full parity claim. The table grows
from the reviewed inventory rather than an optimistic list of names.

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
| `WernerState` | `werner_state` | `MATLABCompat.WernerState` | Only the verified scalar bipartite form is exposed; physical range is enforced | **Partial:** scalar bipartite local/Octave tests pass; multipartite vector form deferred |
| `HorodeckiState` | `horodecki_state` | `MATLABCompat.HorodeckiState` | Local dimensions use a `dims` keyword in the native API | Implemented; both `3×3` and `2×4` local/Octave fixtures pass |
| `GisinState` | `gisin_state` | `MATLABCompat.GisinState` | Mixing probability is validated; angle is in radians | Implemented; local and Octave fixture tests pass |
| `BreuerState` | `breuer_state` | `MATLABCompat.BreuerState` | Even dimension and convex weight are validated; sparse output is a keyword | Implemented; local and Octave fixture tests pass |
| `BrauerStates` | `brauer_states` | `MATLABCompat.BrauerStates` | Returns a sparse matrix; numeric output type is a keyword | Implemented; local and exact Octave fixture tests pass |
| `ChessboardState` | `chessboard_state` | `MATLABCompat.ChessboardState` | Optional `s`/`t` are native keywords; construction does not issue an implicit PPT verdict | Implemented; local and Octave fixture tests pass |

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

`RandomSuperoperator` and `RandomPPTState` are not exported. They remain
explicitly deferred until channel representations and optimization/backend
semantics are complete.

## Tier C channels and maps

Tier C introduces a project-native representation layer around three concrete
types: `KrausRepresentation`, `ChoiRepresentation`, and
`SuperoperatorRepresentation`, all below `AbstractMapRepresentation`. The
project-native accessors and conversions are `input_dimension`,
`output_dimension`, `kraus_representation`, `choi_representation`,
`superoperator_representation`, and `superoperator_matrix`; the
tolerance-controlled numerical diagnostics are `is_completely_positive`,
`is_trace_preserving`, and `is_unital`. These 13 bindings are recorded as
project-native rather than being assigned fabricated one-to-one QETLAB
origins.

The representation convention is an unnormalized, input-first Choi matrix

```math
J(\Phi) = \sum_{i,j} E_{ij} \otimes \Phi(E_{ij}),
```

and Julia's column-major vectorization, so
`vec(Φ(X)) == superoperator_matrix(Φ) * vec(X)`. Choi/superoperator
conversion is an index reshuffle and preserves sparse storage. Canonical
Choi-to-Kraus recovery is defined only for completely positive maps and
explicitly densifies for a tolerance-controlled Hermitian eigendecomposition.
The three physicality functions are numerical tests at the requested
tolerances, not optimization certificates, and do not repair or symmetrize the
stored map.

The committed Tier C Octave/QETLAB artifact has 7 deterministic fixtures and
checks both native and compatibility results for 28 assertions. Its SHA-256 is
`46b31802d97ee8366da163d7c723408ac708c94325aa96c83a34307355a15c03`.
This is supplemental function-specific evidence, not general MATLAB
equivalence.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `ApplyMap` | `apply_channel` | `MATLABCompat.ApplyMap` | Native input is a typed representation or one-sided CP Kraus collection; two-column left/right Kraus cells and rectangular row/column operator spaces are rejected | **Partial:** supported square-operator-space slice passes local and Octave fixture tests |
| `ChoiMatrix` | `choi_matrix` | `MATLABCompat.ChoiMatrix` | Native output uses the package's input-first convention; the wrapper implements `SYS=1/2` for supported CP Kraus input and copies raw Choi input | **Partial:** two-sided Kraus and rectangular row/column operator-space forms are absent |
| `KrausOperators` | `kraus_operators` | `MATLABCompat.KrausOperators` | Only canonical CP recovery is exposed; general non-CP QETLAB left/right outputs are not fabricated | **Partial:** CP slice and explicit non-CP rejection are locally tested |
| `ComplementaryMap` | `complementary_channel` | `MATLABCompat.ComplementaryMap` | Kraus input fixes the dilation; matrix representations use canonical CP recovery and retain their representation kind | **Partial:** CP square-operator-space slice is locally tested |
| `DualMap` | `dual_channel` | `MATLABCompat.DualMap` | Native dispatch preserves representation kind for square matrix algebras | **Partial:** two-sided Kraus cells and rectangular row/column operator spaces are absent |
| `PartialMap` | `partial_map` | `MATLABCompat.PartialMap` | Native API takes one subsystem and a dimension tuple; unequal input/output Hilbert dimensions are supported, but independent row/column subsystem layouts are not | **Partial:** supported multipartite square-operator-space slice passes local and Octave fixture tests |
| `DepolarizingChannel` | `depolarizing_channel` | `MATLABCompat.DepolarizingChannel` | Native constructor enforces the CPTP interval; the wrapper preserves QETLAB's raw Choi formula for every finite real `P` | Implemented; local and Octave fixture tests pass |
| `DephasingChannel` | `dephasing_channel` | `MATLABCompat.DephasingChannel` | Native constructor enforces the CP interval; the wrapper preserves QETLAB's raw Choi formula for every finite real `P` | Implemented; local and Octave fixture tests pass |
| `PauliChannel` | `pauli_channel` | `MATLABCompat.PauliChannel` | Native probabilities are validated without clipping/normalization; random compatibility form is `PauliChannel(rng::AbstractRNG, Q)` | Implemented; deterministic fixture plus explicit-RNG/global-stream tests pass |
| `ChoiMap` | `choi_map` | `MATLABCompat.ChoiMap` | Native output is a typed general-map representation; wrapper output is the raw Choi matrix | Implemented; local and exact Octave fixture tests pass |
| `ReductionMap` | `reduction_map` | `MATLABCompat.ReductionMap` | Native output is a typed general-map representation; wrapper output is the raw Choi matrix | Implemented; local and exact Octave fixture tests pass |

“Rectangular channel” and “rectangular operator space” are different here.
The native representation layer supports channels between unequal Hilbert
dimensions, for example maps from `2×2` matrices to `3×3` matrices. The
remaining compatibility limitation is a general linear map from independently
rectangular matrix spaces, such as `M_{m,n}` to `M_{p,q}`, encoded by separate
row and column dimension arrays in QETLAB.

## Tier D measures and entanglement criteria

Tier D adds 22 native scalar-measure/criterion bindings, 11 reviewed QETLAB
entry-point wrappers, and a separate 10-binding project-native orchestration
layer. The focused suites pass 162 measures/criteria assertions and 68
pipeline assertions. The committed Octave 11.3.0/QETLAB artifact contains 13
fixtures, passes 34 comparisons, and has SHA-256
`ad0cdc45077390fc1eb736fc7c7ff1ec41696c796a508b536774cb6e0020160a`.
This is supplemental function-specific evidence, not general MATLAB
equivalence.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `TraceNorm` | `trace_norm` | `MATLABCompat.TraceNorm` | Sparse input requires explicit `allow_densify=true`; dependency-free implementation computes a full SVD | Implemented; local and Octave fixture tests pass |
| `SchattenNorm` | `schatten_norm` | `MATLABCompat.SchattenNorm` | Native order is validated as real `p ≥ 1`, including `Inf`; sparse densification is opt-in | Implemented; local and Octave fixture tests pass |
| `KyFanNorm` | `ky_fan_norm` | `MATLABCompat.KyFanNorm` | Native `k` is validated against `min(size(X)...)`; implementation computes the full spectrum | Implemented; local and Octave fixture tests pass |
| `Purity` | `purity` | `MATLABCompat.Purity` | Native function requires a validated density matrix; wrapper deliberately preserves QETLAB's unchecked `real(tr(RHO^2))` scalar operation | Implemented; distinct native/wrapper semantics are locally and fixture tested |
| `Entropy` | `von_neumann_entropy` | `MATLABCompat.Entropy` | Native logarithm base is required; wrapper keeps `BASE=2`; every Rényi `ALPHA != 1` request is rejected explicitly | **Partial:** verified von Neumann `ALPHA=1` branch only |
| `Fidelity` | `fidelity` | `MATLABCompat.Fidelity` | Default is QETLAB's unsquared Uhlmann root fidelity; native `squared=true` is explicit; density inputs are strictly validated | Implemented; local and Octave fixture tests pass |
| `Negativity` | `negativity` | `MATLABCompat.Negativity` | Native dimensions and transposed subsystem set are explicit; wrapper retains equal-dimension inference and subsystem 2 | Implemented; local and Octave fixture tests pass |
| `SchmidtDecomposition` | `schmidt_decomposition` | `MATLABCompat.SchmidtDecomposition` | Native output is `SchmidtDecompositionResult`; wrapper returns a named tuple and preserves reviewed `K` selection rather than MATLAB output arity | Implemented; reconstruction and coefficient fixtures pass |
| `SchmidtRank` | `schmidt_rank` | `MATLABCompat.SchmidtRank` | Native absolute/relative tolerances are explicit; result is a tolerance-defined numerical rank | Implemented; local and Octave fixture tests pass |
| `Concurrence` | `concurrence` | `MATLABCompat.Concurrence` | Domain is explicitly normalized two-qubit pure vectors or `4×4` density matrices | Implemented; local and Octave fixture tests pass |
| `IsPPT` | `ppt_criterion` | `MATLABCompat.IsPPT` | Native input is a density matrix; wrapper also accepts finite Hermitian operators. Both return structured tri-state `CriterionResult`, intentionally not QETLAB's boundary-collapsing Boolean | Implemented with documented tri-state result; local and Octave fixture tests pass |

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
`backend_capabilities`, and `available_entanglement_backends`. The
certificate-first matrix path tries PPT, realignment/CCNR, and reduction
criteria. It certifies separability from PPT only in the exact bipartite
`2×2`/`2×3` domain. For pure vectors, a trailing Schmidt coefficient above the
configured threshold certifies entanglement, but separability is certified only
when the computed trailing coefficients are exactly zero; a tolerance-defined
rank-one result with nonzero trailing coefficients is `unknown`. These bindings
do not implement or claim parity with QETLAB `IsSeparable`.

## Tier E coherence slice

The focused coherence suite passes 52 local assertions. Its committed
Octave 11.3.0/QETLAB artifact contains six deterministic fixtures and passes
25 native/compatibility/discrepancy assertions, with SHA-256
`11bcaaee88fac8a595e9a4eff164432dbaa4e141cdd26554da2522e22981811a`.
MATLAB was not run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `L1NormCoherence` | `l1_coherence` | `MATLABCompat.L1NormCoherence` | Native pure vectors and density matrices are strictly validated; sparse density validation requires explicit densification permission | Implemented; local and Octave fixture tests pass |
| `RelEntCoherence` | `relative_entropy_coherence` | `MATLABCompat.RelEntCoherence` | Native logarithm base is required; wrapper keeps QETLAB's base-two default | Implemented; local and Octave fixture tests pass |
| `CoherenceRank` | `coherence_rank` | `MATLABCompat.CoherenceRank` | Native tolerances and optional basis are keywords. Both Julia paths count nonzero coefficients as documented rather than reproducing the pinned implementation's zero-counting bug | Implemented with documented upstream bug fix and committed discrepancy fixtures |

Pure-vector paths preserve generic floating-point precision and sparse storage.
Matrix entropy/positivity work has the same explicit BLAS-type and
`allow_densify` policy as the Tier D spectral routines. See
[Coherence](coherence.md).

## Tier E product analysis and separable-ball slice

This slice adds 10 native bindings and six reviewed QETLAB entry-point
wrappers. The native focused suite passes 175 assertions and the compatibility
suite passes 52. The committed Octave 11.3.0/QETLAB artifact contains 14
deterministic fixtures and passes 68 native, wrapper, reconstruction, and
discrepancy assertions, with SHA-256
`ab6414c1a684141db74782616d4c18e79c8e6039aad53a695d8c723eed598d85`.
MATLAB was not run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `OperatorSchmidtDecomposition` | `operator_schmidt_decomposition`, `operator_schmidt_coefficients` | `MATLABCompat.OperatorSchmidtDecomposition` | Native result is typed, supports separate rectangular row/column layouts, returns the full thin decomposition, and performs no Hermitian-factor repair; wrapper preserves reviewed `K` selection | **Partial compatibility:** supported coefficients/reconstruction pass; pinned rectangular and unequal-dimension Hermitian failures are recorded rather than reproduced |
| `OperatorSchmidtRank` | `operator_schmidt_rank` | `MATLABCompat.OperatorSchmidtRank` | Native rank uses explicit absolute/relative tolerances; wrapper preserves QETLAB's reviewed default absolute threshold | Implemented; local and Octave fixture tests pass |
| `IsProductVector` | `is_product_vector` | `MATLABCompat.IsProductVector` | Returns `ProductAnalysisResult` with factors, scaled cut residuals, reconstruction residual, thresholds, and a three-way status instead of collapsing a tolerance boundary to a Boolean | Implemented with intentional structured result; local and Octave fixture tests pass away from boundaries |
| `IsProductOperator` | `is_product_operator` | `MATLABCompat.IsProductOperator` | Supports multipartite rectangular local operators and returns the same structured numerical analysis instead of a Boolean | Implemented with intentional structured result; local and Octave fixture tests pass away from boundaries |
| `EntFormation` | `entanglement_of_formation` | `MATLABCompat.EntFormation` | Exact pure bipartite and mixed two-qubit domains only; native base is explicit with default two and bounded projection policies are named; higher-dimensional rank-one matrices are not auto-converted | **Partial compatibility:** verified exact domains pass; pinned zero-concurrence `NaN` is corrected to mathematical zero |
| `InSeparableBall` | `in_separable_ball` | `MATLABCompat.InSeparableBall` | Native input must already be normalized and the result distinguishes a certificate, outside-ball failure, and numerical unknown; wrapper retains safe positive-trace normalization | Implemented with structured sufficient-certificate semantics; local and Octave fixture tests pass |

`ProductAnalysisResult` is a numerical classification, not an exact symbolic
proof. `SeparableBallResult(:outside_ball, ...)` is a failed sufficient test,
not an entanglement conclusion. General matrix spectral work requires explicit
sparse densification; diagonal matrices and supplied eigenvalues have
structure-aware paths. See
[Product structure and separable-ball certificates](product_analysis.md).

## Tier E matrix analysis

The native matrix-analysis suite passes 127 local assertions, and the separate
compatibility suite passes 32. The committed Octave 11.3.0/QETLAB artifact
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

The five discrepancy fixtures are evidence, not native expected values: weak
totals, negative-vector padding order, row-matrix vector semantics,
rectangular high-order compound shape, and the additive order-zero error. See
[Matrix analysis](matrix_analysis.md) for definitions and complexity.

## Tier E matrix predicates

The native predicate suite passes 166 assertions and its compatibility suite
passes 37. No MATLAB predicate oracle has been run.

| QETLAB function | Julia-native function | Compatibility name | Changed arguments/conventions | Status |
|---|---|---|---|---|
| `IsPSD` | `is_positive_semidefinite` | `MATLABCompat.IsPSD` | Returns a three-valued `MatrixPredicateResult`; neither Julia path silently symmetrizes input, and the QETLAB CVX branch is omitted | Partial: numeric matrices are implemented and tested; the CVX symbolic branch is not |
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
