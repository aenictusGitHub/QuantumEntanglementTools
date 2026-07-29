# ADR 0003: indexing and representation conventions

- Status: Accepted
- Date: 2026-07-28

## Context

Subsystem routines can be individually plausible while disagreeing about which
factor varies fastest, how pure states reduce, or whether partial transpose
conjugates. These choices are observable API behavior and must precede broad
porting.

## Decision

For dimensions `(d₁, …, dₙ)`:

- factors and dimensions are written left-to-right;
- subsystem labels are one-based;
- subsystem 1 is the most-significant/slowest-varying factor and subsystem `n`
  is the least-significant/fastest-varying factor;
- for one-based component indices `(i₁, …, iₙ)`, the one-based linear index is
  `1 + Σₖ (iₖ - 1) ∏_{j>k} dⱼ`;
- a pure-state reduction means reducing `ψ * ψ'` and returns an operator matrix;
- tracing all subsystems retains an explicit `1 × 1` matrix rather than changing
  the result type to a scalar;
- partial transpose swaps the selected subsystem's row and column indices. It is
  an ordinary transpose on those indices, with no conjugation;
- `vec` uses Julia/MATLAB column-major vectorization.

The all-subsystem `1 × 1` shape and the Tier A indexing decisions have local
executable coverage. Cross-platform and differential validation remain pending.

The channel layer uses the unnormalized Choi convention

```math
J(\Phi) = \sum_{i,j} |i\rangle\langle j| \otimes
          \Phi(|i\rangle\langle j|),
```

with input factor first and output factor second. Then trace preservation means
`tr_output(J) = I_input` and `tr(J) = d_input`. Choi/Kraus/superoperator
round-trip and physicality tests now cover this decision for the implemented
Tier C scope.

The default maximally entangled state is normalized. Normalizations and ordering
for the implemented Pauli, generalized Pauli, Gell-Mann, spin, Bell, and named
state families are documented and tested; families outside the reviewed
inventory scope remain unresolved.

## Consequences

- Reshape/permutation kernels must implement this product-basis order explicitly;
  Julia's column-major array dimension order must not be mistaken for subsystem
  significance.
- Unsorted subsystem selections may be accepted only with documented
  canonicalization; repeated and out-of-range labels are errors.
- Compatibility wrappers must translate a genuine QETLAB convention difference
  instead of silently changing the core contract.
- Executable tiny-system tests are required for every formula in this ADR.
