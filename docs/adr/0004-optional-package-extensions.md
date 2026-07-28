# ADR 0004: optional package extension boundary

- Status: Accepted
- Date: 2026-07-28

## Context

Entanglement detection, optimization, symbolic work, plotting, and device
support can pull in large dependencies or global behaviors irrelevant to core
linear algebra. Users should not need a solver to compute a partial trace.

## Decision

Keep core result types and high-level operations package-owned. Integrate
third-party packages through Julia weak dependencies and package extensions.
Extensions:

- map stable configuration values to documented public backend APIs;
- report backend capabilities and versions;
- make representation conversions explicit;
- preserve witnesses, decompositions, statuses, residuals, and warnings;
- contain backend exceptions as structured failures rather than mathematical
  negatives;
- avoid type piracy, private APIs, re-exporting whole namespaces, and global
  RNG/logging/thread side effects;
- load correctly in both package orders and disappear cleanly when the optional
  dependency is absent.

## Consequences

- Dedicated environments and fresh-process load-order tests are required.
- A backend's raw result may be carried explicitly, but ordinary core results
  must not leak private backend types.
- Version-specific global-state hazards are release gates. In particular,
  EntanglementDetection.jl 0.2.2 is not considered safely integrated until its
  observed RNG, stdout/logging, and BLAS-thread mutations are isolated and
  regression-tested.
