using Documenter
using QuantumEntanglementTools

include(joinpath(@__DIR__, "..", "scripts", "check_docs_math_compat.jl"))
DocsMathCompatibility.check_docs_math_compat() ||
    error("documentation math compatibility checks failed")

DocMeta.setdocmeta!(
    QuantumEntanglementTools,
    :DocTestSetup,
    :(using QuantumEntanglementTools);
    recursive=true,
)

makedocs(;
    sitename="QuantumEntanglementTools.jl",
    modules=[QuantumEntanglementTools],
    repo=Documenter.Remotes.GitHub("aenictusGitHub", "QuantumEntanglementTools"),
    checkdocs=:exports,
    doctest=true,
    warnonly=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        collapselevel=1,
        edit_link="main",
        assets=[
            Documenter.asset("assets/qet_code_generator.css"; class=:css, islocal=true),
            Documenter.asset(
                "assets/qet_code_generator_core.js";
                class=:js,
                islocal=true,
                attributes=Dict(:defer => ""),
            ),
            Documenter.asset(
                "assets/qet_code_generator_ui.js";
                class=:js,
                islocal=true,
                attributes=Dict(:defer => ""),
            ),
        ],
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Code generator" => "code_generator.md",
        "Separability by example" => "separability_examples.md",
        "Symmetric SAPPT states and witnesses" => "paper_symmetric_separability.md",
        "Architecture" => "architecture.md",
        "Conventions" => "conventions.md",
        "States, operators, and random objects" => "states_operators_random.md",
        "Coherence" => "coherence.md",
        "Product structure and separable-ball certificates" => "product_analysis.md",
        "Matrix analysis" => "matrix_analysis.md",
        "Matrix predicates" => "matrix_predicates.md",
        "Migration from QETLAB" => "migration_from_qetlab.md",
        "Entanglement backends" => "entanglement_backends.md",
        "EntanglementDetection extension" => "entanglement_detection_extension.md",
        "Optimization and solvers" => "optimization_and_solvers.md",
        "Sparse and large-scale work" => "sparse_and_large_scale.md",
        "Reproducibility" => "reproducibility.md",
        "Performance" => "performance.md",
        "External integrations" => "integrating_external_packages.md",
        "Maintaining upstream parity" => "maintaining_upstream_parity.md",
        "Executable tutorials" => "tutorials.md",
        "API reference" => [
            "Overview" => "api/index.md",
            "Julia-native exports" => "api/native.md",
            "MATLAB compatibility" => "api/matlab_compat.md",
        ],
    ],
)
