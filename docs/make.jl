using Documenter
using QuantumEntanglementTools

DocMeta.setdocmeta!(
    QuantumEntanglementTools,
    :DocTestSetup,
    :(using QuantumEntanglementTools);
    recursive=true,
)

makedocs(;
    sitename="QuantumEntanglementTools.jl",
    modules=[QuantumEntanglementTools],
    remotes=nothing,
    checkdocs=:exports,
    doctest=true,
    warnonly=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        collapselevel=1,
        edit_link=nothing,
        repolink=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Architecture" => "architecture.md",
        "Conventions" => "conventions.md",
        "States, operators, and random objects" => "states_operators_random.md",
        "Coherence" => "coherence.md",
        "Product structure and separable-ball certificates" => "product_analysis.md",
        "Matrix analysis" => "matrix_analysis.md",
        "Matrix predicates" => "matrix_predicates.md",
        "Migration from QETLAB" => "migration_from_qetlab.md",
        "Entanglement backends" => "entanglement_backends.md",
        "Optimization and solvers" => "optimization_and_solvers.md",
        "Sparse and large-scale work" => "sparse_and_large_scale.md",
        "Reproducibility" => "reproducibility.md",
        "Performance" => "performance.md",
        "External integrations" => "integrating_external_packages.md",
        "Maintaining upstream parity" => "maintaining_upstream_parity.md",
        "Tutorials" => "tutorials/index.md",
        "API reference" => "api/index.md",
    ],
)
