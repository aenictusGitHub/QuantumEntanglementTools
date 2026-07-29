# API reference

The Tier A subsystem kernel passes 304 local tests and the focused Tier B
operators/states/random suite passes 941 local tests on Julia 1.12.6. The Tier
C channels/maps suite passes 154 focused local assertions across 24 native
public bindings/types and 11 `MATLABCompat` wrappers. Tier D passes 162 focused
measures/criteria assertions across 22 native bindings and 11 wrappers. The
native entanglement pipeline passes 68 assertions across 10 project-native
orchestration exports, and the optional backend boundary adds two public
descriptor/configuration types with 125 focused extension assertions. Tier E
coherence passes 52. The integrated package corpus passes 2,318 assertions
(2,270 core plus 48 executable tutorials) on Julia 1.12.6 and 1.10.11. Tier E
product analysis passes 175 native and 52 compatibility assertions across 10
native bindings and six wrappers. Tier E
matrix analysis passes 127 native and 32 compatibility assertions across four
native bindings and four compatibility entry points. Tier E matrix predicates
pass 166 native and 37 compatibility assertions across nine native bindings
and four structured-result compatibility entry points. The reference below
covers all exported native operations and compatibility wrappers. This local
evidence does not establish full QETLAB parity, supported-platform coverage,
or API stability.

Randomized constructors require an explicit `rng::AbstractRNG`, including
their `MATLABCompat` spellings. See
[States, operators, and random objects](../states_operators_random.md) for
examples, physical parameter validation, and deferred capabilities.
Channel representation conventions, general-map limitations, raw compatibility
constructors, and the explicit-RNG `PauliChannel` form are documented in
[Migration from QETLAB](../migration_from_qetlab.md#tier-c-channels-and-maps).
Tier D tri-state criteria, certificate semantics, the intentionally partial
`Entropy` mapping, and the native pipeline are documented in
[Entanglement backends](../entanglement_backends.md) and
[Migration from QETLAB](../migration_from_qetlab.md#tier-d-measures-and-entanglement-criteria).
The exact-version optional adapter, its Julia 1.11 resolver floor, and its
uncertified child-process result semantics are documented in
[EntanglementDetection.jl extension](../entanglement_detection_extension.md).
Tier E operator Schmidt decompositions, structured product analysis,
closed-form entanglement of formation, and sufficient separable-ball
certificates are documented in
[Product structure and separable-ball certificates](../product_analysis.md).
Strong and weak majorization conventions, exact symmetric polynomials, and
compound-matrix boundary shapes are documented in
[Matrix analysis](../matrix_analysis.md).
Structured predicate outcomes, witnesses, tolerance boundaries, and
combinatorial guards are documented in
[Matrix predicates](../matrix_predicates.md).

API pages are generated or reviewed against exports and the upstream inventory.
Every public docstring must state signatures, mathematical
definition, input/output dimensions, conventions, keyword defaults, error
behavior, exact/numerical/heuristic/certified meaning, references and provenance,
a runnable example, and complexity for important operations.

Internal helpers are not public merely because Documenter can render them.

## Julia-native exports

```@autodocs
Modules = [QuantumEntanglementTools]
Private = false
```

## MATLAB compatibility namespace

```@autodocs
Modules = [QuantumEntanglementTools.MATLABCompat]
Private = false
```
