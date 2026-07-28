# ADR 0005: solver-independent optimization layer

- Status: Proposed
- Date: 2026-07-28

## Context

Several upstream capabilities rely on CVX or other optimization toolboxes.
Julia offers multiple modeling layers and solvers with different cone support,
licenses, numerical behaviors, and status vocabularies. Choosing one without
implementing representative primal and dual prototypes would lock the package
to an untested abstraction.

## Proposed decision

Keep optimization out of the core dependency set. Define package-owned problem
configuration and result types, then use extensions for a selected modeling
layer and solvers. Every supported formulation must:

- be written mathematically in the docs before implementation;
- expose termination, primal and dual status, objectives, gaps, residuals,
  iterations, tolerances, and stopping reason;
- distinguish certificates, bounds, heuristics, and unknown outcomes;
- handle infeasible, unbounded, inaccurate, time-limited, and failed solves;
- have small primal/dual or independent cross-checks.

## Decision gate

Prototype representative SDP and non-SDP routines using public APIs of candidate
modeling layers. Compare cone coverage, extension ergonomics, type behavior,
result/status fidelity, maintenance, and licenses. Replace this ADR with an
accepted or superseding record only after those results are documented.
