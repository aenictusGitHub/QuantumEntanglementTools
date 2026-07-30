# API reference

At pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`, the strict static completion
checker passes with 127/127 public rows carrying final `verified` status, 36/36
internal helpers assigned terminal dispositions, no queued rows, and 458
exported bindings with matching provenance entries. The direct local full
corpus passes 8,133/8,133 assertions, including 48 executable-tutorial
assertions, on Julia 1.12.6 and the installed Julia 1.10.0.

The reference below covers all exported native operations and compatibility
wrappers. These results establish local implementation and validation evidence,
not MATLAB/QETLAB parity, remote supported-platform CI, comparative
performance, API stability, release approval, or human review.

Randomized constructors and randomized lower-bound routines such as
`random_superoperator`, `induced_matrix_norm`, and
`induced_schatten_lower_bound` require an explicit
`rng::AbstractRNG`, including their
`MATLABCompat` spellings. See
[States, operators, and random objects](../states_operators_random.md) for
examples, physical parameter validation, resource limits, and compatibility
differences.
Channel representation conventions, general-map limitations, raw compatibility
constructors, and the explicit-RNG `PauliChannel` form are documented in
[Migration from QETLAB](../migration_from_qetlab.md#tier-c-channels-and-maps).
The four-dimensional operator-space descriptor, two-sided sums, rectangular
Choi/superoperator forms, dual convention, structured complete-positivity
diagnostics, canonical CP/Hermitian/general factor branches, and complementary
map boundary are documented in
[General maps and rectangular operator spaces](../general_maps.md).
Exact-versus-lower-bound semantics, matrix witnesses, zero-gradient status,
and deterministic SVD/work guards for induced Schatten norms are documented
in [Induced Schatten norm lower bounds](../induced_schatten_norm.md).
Explicit marginal constraints, Kraus complete-positivity certificates,
bounded factor-space balancing, corrected unequal-dimensional semantics, and
resource guards for random maps are documented in
[Bounded random completely positive maps](../random_superoperators.md).
Sparse permutation/Brauer/projector structure, exact dependence elimination,
strict dimension/copy validation, and factorial/double-factorial/exponential
guards for the four supported twirls are documented in
[Guarded group twirls](../twirls.md).
Type-preserving homogeneous coefficients, the fixed monomial order, compact
symmetric matrix representations, deterministic sampling, corrected target
semantics, and explicit solver/resource boundaries are documented in
[Solver-free polynomial foundations](../polynomial_foundations.md).
The status-aware symmetric-moment hierarchy, explicit RNG/sample accounting,
and outer-versus-attained-inner semantics are documented in
[Polynomial SOS hierarchy](../polynomial_sos.md).
Exact copositivity witnesses, numerical hierarchy boundaries, and independently
inspectable graph certificates are documented in
[Copositivity and clique-number bounds](../copositivity_and_cliques.md).
Minimum-error POVM construction, Helstrom and orthogonal certificates, and
solver-derived primal/dual bounds are documented in
[Minimum-error state discrimination](../state_discrimination.md).
Completely bounded channel norms, the channel Holevo--Helstrom reduction,
direct maximum-output-root-fidelity models, bounds, statuses, and reviewed
upstream defect corrections are documented in
[Channel norms, discrimination, and output fidelity](../channel_optimization.md).
Exact and witnessed S(`k`) bounds, explicit-RNG projected searches,
solver-neutral PPT/reduction/symmetric-extension relaxations, and tri-state
block-positivity conclusions are documented in
[S(k) operator norms and block positivity](../sk_operator_norm_and_block_positivity.md).
The composable affine-expression branch of `kpNormDual` is documented in
[Dual top-k p-norm model expressions](../top_k_p_norm_dual_models.md).
The solver-free pure-state `k`-coherence formula, exact equality convention,
branch-stability diagnostics, and corrected compatibility behavior are
documented in
[Pure-state robustness of k-coherence](../pure_k_coherence_robustness.md).
Full-definition UPB validation, exact and tolerance-aware partition
certificates, compact extension witnesses, and explicit combinatorial limits
are documented in
[UPB certificates and extension witnesses](../is_upb.md).
The complete executable UPB construction catalog, stable local/global result,
dimension dispatch, explicit-RNG randomized branches, family provenance, and
allocation/work limits are documented in
[UPB construction catalog](../upb_catalog.md).
Tier D tri-state criteria, certificate semantics, numeric Rényi entropy,
structured exact-versus-lower-bound induced norms, support-aware Matsumoto
fidelity, and the native pipeline are documented in
[Entanglement backends](../entanglement_backends.md) and
[Migration from QETLAB](../migration_from_qetlab.md#tier-d-measures-and-entanglement-criteria).
The composite certificate-first `is_separable` report, local-discrimination
models, and exact-versus-inconclusive decision boundaries are documented in
[Separability and local discrimination](../separability_optimization.md).
The exact-version optional adapter, its Julia 1.11 resolver floor, and its
uncertified child-process result semantics are documented in
[EntanglementDetection.jl extension](../entanglement_detection_extension.md).
Tier E operator Schmidt decompositions, structured product and
entangling-gate analysis, exact-or-unknown UPB-size lookup, closed-form
entanglement of formation, and sufficient separable-ball certificates are documented in
[Product structure and separable-ball certificates](../product_analysis.md).
The bounded multipartite operator-scaling result, residual and conditioning
diagnostics, successful-only compatibility outputs, and explicit resource
limits are documented in
[Bounded operator Sinkhorn scaling](../operator_sinkhorn.md). The
rank-aware bipartite decomposition, reconstruction checks, and five-output
compatibility wrapper are documented in
[Bounded filter normal form](../filter_normal_form.md).
Strong and weak majorization conventions, exact symmetric polynomials,
compound-matrix boundary shapes, and bounded numerical commutants are
documented in
[Matrix analysis](../matrix_analysis.md).
The checked, type-preserving full-probability tensor construction for
parallel nonlocal-game repetition is documented in
[Nonlocal games](../nonlocal_games.md).
Exact classical enumeration, NPA and no-signalling relaxations, and
explicit-RNG attained candidates are documented in
[Bell inequalities and nonlocal games](../nonlocal_optimization.md).
Structured predicate outcomes, witnesses, tolerance boundaries, and
combinatorial guards are documented in
[Matrix predicates](../matrix_predicates.md).

API pages are generated or reviewed against exports and the upstream inventory.
Every public docstring must state signatures, mathematical
definition, input/output dimensions, conventions, keyword defaults, error
behavior, exact/numerical/heuristic/certified meaning, references and provenance,
a runnable example, and complexity for important operations.

Internal helpers are not public merely because Documenter can render them.

The generated reference is split by domain so each page remains
comfortably below Documenter's strict HTML-size limit:

- [Julia-native core and analysis](native.md)
- [Julia-native channels and general maps](native_channels.md)
- [Julia-native entanglement analysis](native_entanglement.md)
- [Julia-native optimization models](native_optimization.md)
- [MATLAB compatibility namespace](matlab_compat.md)
