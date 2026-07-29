# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is a pre-alpha, independent Julia package under
development for quantum-information and entanglement calculations.

!!! warning "No completeness claim"
    The public numerical API, QETLAB inventory, validation matrix, and benchmark
    baseline are incomplete. A page describing an intended capability is not
    evidence that the capability is implemented.

The architecture centers on ordinary Julia arrays, explicit subsystem
conventions, generic numeric types, sparse-aware algorithms, structured
diagnostic and certification results, and optional package extensions.

## Where to begin

- [Getting started](getting_started.md) explains the current local-development
  workflow.
- [Conventions](conventions.md) separates accepted Tier A choices from future
  choices that still need tests.
- [Separability by example](separability_examples.md) constructs pure and mixed
  product-state examples and explains certified versus inconclusive outcomes.
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
  separability, and entanglement-analysis scripts exercised by the package
  tests.
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
