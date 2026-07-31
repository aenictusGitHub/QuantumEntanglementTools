# Julia-native entanglement analysis

This page contains certificate-aware entanglement criteria, product-structure
analysis, UPB construction and analysis, bounded local filtering, and optional
backend orchestration. Mathematical conclusions remain distinct from necessary tests,
heuristics, numerical boundaries, and resource-limited `unknown` outcomes.

```@autodocs
Modules = [QuantumEntanglementTools]
Pages = [
    "entanglement/backend_interface.jl",
    "result_interface.jl",
    "entanglement/criteria.jl",
    "entanglement/block_positivity.jl",
    "entanglement/entangling_gate.jl",
    "entanglement/absolute_ppt.jl",
    "entanglement/symmetric_extensions.jl",
    "entanglement/upb_size.jl",
    "entanglement/upb_catalog.jl",
    "entanglement/upb_analysis.jl",
    "entanglement/filter_normal_form.jl",
    "entanglement/operator_sinkhorn.jl",
    "entanglement/product_analysis.jl",
]
Private = false
```
