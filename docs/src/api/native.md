# Julia-native core and analysis API

This page contains Julia-native subsystem, state, operator, measure,
entanglement, and analysis bindings. Channel and general-map bindings are on
the separate [native channels and maps](native_channels.md) page.
Certificate, tolerance, ordering, and sparse-storage contracts are defined in
each docstring and in the linked conceptual guides.

```@autodocs
Modules = [QuantumEntanglementTools]
Pages = [
    "QuantumEntanglementTools.jl",
    "dimensions.jl",
    "commutant.jl",
    "matrix_analysis.jl",
    "matrix_predicates.jl",
    "scalar_measures.jl",
    "parallel_repetition.jl",
    "operators.jl",
    "random_objects.jl",
    "entangled_subspace.jl",
    "states.jl",
    "symmetric_states.jl",
    "partial_trace.jl",
    "partial_transpose.jl",
    "permutation.jl",
    "projectors.jl",
    "realignment.jl",
    "tensor_products.jl",
    "coherence.jl",
    "pure_k_robustness.jl",
]
Private = false
```
