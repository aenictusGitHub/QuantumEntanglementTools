# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is an independent Julia package for
quantum-information and entanglement calculations. The experimental `v0.1.0`
release candidate exposes a scoped, certificate-aware API.

!!! warning "No completeness claim"
    The QETLAB inventory, validation matrix, and benchmark baseline are not
    complete. The release candidate covers only functions whose provenance,
    documentation, and tests are recorded; it does not claim QETLAB parity.

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
- [Migration from QETLAB](migration_from_qetlab.md) defines how the compatibility
  ledger will be generated.
- [Entanglement backends](entanglement_backends.md) explains certificate and
  optional-backend boundaries.
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
