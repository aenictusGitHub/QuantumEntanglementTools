# Julia-native optimization models

This page contains the solver-neutral affine SDP representation, resource
limits, backend configuration, status-rich result types, and subsystem-aware
linear maps. JuMP and concrete solvers remain optional; see
[Optimization and solvers](../optimization_and_solvers.md) for the modeling,
certificate, and failure-semantics contract.

```@autodocs
Modules = [QuantumEntanglementTools]
Pages = [
    "optimization/interface.jl",
    "optimization/linear_maps.jl",
    "optimization/psd_constraints.jl",
    "optimization/matsumoto_fidelity_model.jl",
    "optimization/top_k_p_norm_epigraph.jl",
    "optimization/top_k_p_norm_dual_epigraph.jl",
    "optimization/sk_operator_norm.jl",
    "optimization/state_discrimination.jl",
    "optimization/channel_optimization.jl",
    "optimization/polynomial_sos.jl",
    "optimization/copositivity_clique.jl",
    "coherence/coherence_optimization.jl",
    "entanglement/separability_optimization.jl",
    "nonlocal_games/scenarios.jl",
    "nonlocal_games/npa.jl",
    "nonlocal_games/game_values.jl",
    "nonlocal_games/nonlocal_optimization.jl",
]
Private = false
```
