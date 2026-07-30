# Architecture decision records

Architecture decision records (ADRs) capture choices that affect public
semantics, dependencies, or long-lived maintenance. `Proposed` decisions are not
implementation contracts. `Accepted` decisions remain revisable through a new
ADR that supersedes the old one; do not silently rewrite history after release.

| ADR | Title | Status |
|---|---|---|
| [0001](0001-package-name.md) | Provisional package name | Accepted, provisional |
| [0002](0002-core-architecture.md) | One Julia-native core over standard arrays | Accepted |
| [0003](0003-indexing-and-representation-conventions.md) | Indexing and representation conventions | Accepted |
| [0004](0004-optional-package-extensions.md) | Optional package extension boundary | Accepted |
| [0005](0005-optimization-layer.md) | Package-owned conic models with an optional JuMP/MOI translator | Accepted |
| [0006](0006-general-operator-space-maps.md) | Rectangular operator spaces and two-sided map sums | Accepted and implemented |
