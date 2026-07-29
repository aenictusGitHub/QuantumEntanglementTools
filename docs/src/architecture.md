# Architecture

The accepted foundation is one Julia-native core over ordinary array
abstractions, plus isolated compatibility wrappers and optional package
extensions.

```text
standard vectors/matrices + subsystem dimensions
                      |
          Julia-native validated core
             /        |         \
     core results   MATLABCompat  optional extensions
                         |        /       |        \
                    QETLAB names backend  solver  future adapters
```

The core owns public inputs, configuration, result types, and certificate
semantics. Optional dependencies may implement methods but may not redefine the
meaning of a failed or inconclusive result.

## Layers

1. Dimension, indexing, validation, and tolerance policy.
2. Subsystem and representation primitives.
3. States, operators, channels, maps, measures, and criteria.
4. Structured analysis and certification orchestration.
5. Compatibility wrappers and optional backend/solver integrations.

Implementation follows dependency order. A later layer must not compensate for
an ambiguous convention in an earlier layer.

## Decisions

The architecture decision records are maintained in `docs/adr/`. In particular:

- ADR 0002 selects standard arrays and one top-level module.
- ADR 0003 fixes subsystem/indexing and implemented channel-representation
  conventions.
- ADR 0004 defines the optional-extension boundary.
- ADR 0005 deliberately leaves the optimization modeling layer proposed until
  representative prototypes are evaluated.
