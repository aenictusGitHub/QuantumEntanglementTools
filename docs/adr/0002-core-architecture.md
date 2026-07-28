# ADR 0002: one Julia-native core over standard arrays

- Status: Accepted
- Date: 2026-07-28

## Context

Quantum-information operations need clear subsystem metadata but should compose
with Julia's linear-algebra ecosystem. Requiring a custom state container would
add conversion and dispatch friction. A mechanical MATLAB transliteration would
also preserve accidental conventions and prevent idiomatic keyword and
multiple-dispatch design.

## Decision

Use one public top-level module, organized by source files rather than a forest
of user-visible nested modules. The native numerical API accepts standard Julia
abstractions:

- `AbstractVector{<:Number}` for pure states;
- `AbstractMatrix{<:Number}` for operators and density matrices;
- structured and sparse matrix wrappers where the operation supports them;
- tuples or vectors of positive subsystem dimensions.

Small immutable dimension/layout and reusable plan types may be introduced when
they clarify validation or have measured reuse benefits. Common operations must
remain callable without wrapping data in a custom quantum-state type.

The primary API uses lowercase `snake_case`, multiple dispatch, and keywords.
QETLAB/MATLAB call patterns belong in a separate compatibility namespace whose
wrappers delegate to verified native methods.

## Consequences

- Core algorithms need correct generic fallbacks and may add strided fast paths.
- Sparse preservation is an API obligation; implicit densification is a defect.
- Representation conversions must be explicit and tested.
- Source organization does not define separate type universes.
- Compatibility wrappers cannot become the only tested implementation path.
