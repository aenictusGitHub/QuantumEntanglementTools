# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is an independent Julia package for
quantum-information and entanglement calculations. At pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`, the strict static completion
checker records all 127/127 public rows with final `verified` status, all 36/36
internal helpers with terminal dispositions, an empty completion queue, and
458 exported bindings with matching provenance entries.

!!! warning "Scope of the completion evidence"
    The direct local full corpus passes 8,133/8,133 assertions, including 48
    executable-tutorial assertions, on Julia 1.12.6 and the installed Julia
    1.10.0. These static and local results do not establish MATLAB/QETLAB
    parity, remote supported-platform CI, comparative performance, API
    stability, release approval, or the required human review.

The architecture centers on ordinary Julia arrays, explicit subsystem
conventions, generic numeric types, sparse-aware algorithms, structured
diagnostic and certification results, and optional package extensions.

## Where to begin

- [Getting started](getting_started.md) explains the current local-development
  workflow.
- [Entanglement example code generator](code_generator.md) builds bounded,
  certificate-aware Julia scripts locally in the browser from curated state
  families and analysis routes.
- [Conventions](conventions.md) records the tested subsystem, indexing, and
  channel-representation choices.
- [Separability by example](separability_examples.md) constructs pure and mixed
  product-state examples and explains certified versus inconclusive outcomes.
- [Separability and local discrimination](separability_optimization.md)
  documents the structured, certificate-first `is_separable` interface,
  hierarchy boundaries, and local-measurement relaxations.
- [Symmetric SAPPT states and witnesses](paper_symmetric_separability.md)
  reproduces a research family with an explicit separable decomposition and
  constructive entanglement witnesses.
- [Product structure and separable-ball certificates](product_analysis.md)
  documents tolerance-aware product classifications and sufficient
  certificate semantics.
- [Matrix analysis](matrix_analysis.md) documents strong versus weak
  majorization, exact symmetric polynomials, and compound-matrix shape
  conventions.
- [Matrix predicates](matrix_predicates.md) documents structured
  positive-semidefinite and all-minor results, tolerance boundaries, witnesses,
  sparse policy, and combinatorial guards.
- [Migration from QETLAB](migration_from_qetlab.md) records the reviewed
  compatibility mappings, deliberate corrections, and result-semantics
  differences.
- [Entanglement backends](entanglement_backends.md) explains certificate and
  optional-backend boundaries.
- [Bell inequalities and nonlocal games](nonlocal_optimization.md) separates
  exact finite searches, numerical upper relaxations, and attained numerical
  candidates.
- [Minimum-error state discrimination](state_discrimination.md) demonstrates
  exact Helstrom certificates and status-aware POVM optimization.
- [Channel norms, discrimination, and output fidelity](channel_optimization.md)
  documents analytic certificates, solver-neutral SDPs, explicit bounds, and
  corrected pinned channel-probability and Kraus-rank defects.
- [S(k) operator norms and block positivity](sk_operator_norm_and_block_positivity.md)
  separates exact values, witnessed bounds, numerical relaxations, tri-state
  certificates, and bounded explicit-RNG searches.
- [Dual top-k p-norm model expressions](top_k_p_norm_dual_models.md) supplies
  the solver-neutral affine epigraph missing from the numeric dual-norm API.
- [Polynomial SOS hierarchy](polynomial_sos.md) separates numerical outer
  relaxations from attained sampled inner evidence.
- [Copositivity and clique-number bounds](copositivity_and_cliques.md)
  separates exact witnesses and graph certificates from numerical hierarchy
  evidence.
- [Executable tutorials](tutorials.md) runs the exact subsystem, channel,
  separability, symmetric-witness, and entanglement-analysis scripts exercised
  by the package tests.
- [EntanglementDetection.jl extension](entanglement_detection_extension.md)
  documents the optional child-process adapter and its uncertified-result
  boundary.
- [Reproducibility](reproducibility.md) lists the evidence expected with results.

Project-wide legal, milestone, build-environment, and session-handoff records
live at the top level of the `docs/` directory and are linked from the repository
README.

## Non-affiliation

This project is not affiliated with, endorsed by, or officially supported by
QETLAB, its maintainers, or their contributors.
