# QuantumEntanglementTools.jl

`QuantumEntanglementTools` is an independent Julia package for
quantum-information and entanglement calculations. It works with ordinary
Julia arrays, keeps subsystem dimensions explicit, preserves useful numeric
types and sparsity, and reports certificates separately from inconclusive
numerical evidence.

## Start with a certified result

```@example home-first-result
using QuantumEntanglementTools

ket0 = ComplexF64[1, 0]
ket1 = ComplexF64[0, 1]
psi = tensor_product(ket0, ket1)

analyze_entanglement(psi, (2, 2))
```

This report certifies separability from the pure-state product decomposition.
An `:unknown` report would instead mean that the selected methods established
neither separability nor entanglement.

Continue with the [five-minute quick start](getting_started.md) for
installation, result semantics, and a task-to-function chooser.

## Pick a learning path

| Goal | Start here | Continue with |
|---|---|---|
| Perform a first quantum-information calculation | [Five-minute quick start](getting_started.md) | [Mathematical conventions](conventions.md) |
| Learn subsystem and channel operations | [Executable tutorials](tutorials.md) | [States, operators, and random objects](states_operators_random.md) |
| Decide what a separability result proves | [Separability by example](separability_examples.md) | [Certificate-aware analysis](entanglement_backends.md) |
| Generate a complete example interactively | [Code generator](code_generator.md) | Download and run the generated Julia script |
| Move from QETLAB | [Migration from QETLAB](migration_from_qetlab.md) | [`MATLABCompat` reference](api/matlab_compat.md) |
| Configure an optional SDP or detection backend | [Optimization and solvers](optimization_and_solvers.md) | [External integrations](integrating_external_packages.md) |
| Find a specific function | [Task-based API overview](api/index.md) | Domain reference pages in the navigation |

## Result semantics in one minute

!!! tip "A criterion pass is not necessarily a separability proof"
    `CriterionEntanglementDetected` records positive entanglement evidence.
    `CriterionSatisfied` says only that the state passed that necessary test.
    `CriterionUnknown` records a numerical boundary. Composite reports retain
    `:unknown` whenever no requested route supplies a valid certificate.

The documentation groups the remaining material by user task:
entanglement and separability; states, channels, and maps; coherence and matrix
analysis; optimization and nonlocality; and package maintenance.

## Development status

At pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`, the strict static completion
checker records all 127/127 public rows with final `verified` status, all 36/36
internal helpers with terminal dispositions, an empty completion queue, and
477 exported bindings with matching provenance entries.

!!! warning "Scope of the completion evidence"
    The direct local full corpus passes 9,759/9,759 assertions: the core
    accounts for 9,675 assertions, and the seven standalone executable
    tutorials pass 84/84 assertions (48 existing plus 36 for the two new
    workflows), on Julia 1.12.6 and the installed Julia 1.10.11. These static and
    local results do not establish MATLAB/QETLAB parity, remote
    supported-platform CI, comparative performance, API stability, release
    approval, or the required human review.

Project-wide legal, milestone, build-environment, and session-handoff records
live at the top level of the `docs/` directory and are linked from the repository
README.

## Non-affiliation

This project is not affiliated with, endorsed by, or officially supported by
QETLAB, its maintainers, or their contributors.
