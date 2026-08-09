using Documenter
using QuantumEntanglementTools

edit_ref = if get(ENV, "GITHUB_REF_TYPE", "") == "tag"
    get(ENV, "GITHUB_REF_NAME", "main")
else
    "main"
end

include(joinpath(@__DIR__, "..", "scripts", "check_docs_math_compat.jl"))
DocsMathCompatibility.check_docs_math_compat() ||
    error("documentation math compatibility checks failed")

include(joinpath(@__DIR__, "..", "scripts", "check_optional_backend_docs.jl"))
OptionalBackendDocs.check_optional_backend_docs() ||
    error("optional-backend documentation checks failed")

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
        edit_link=edit_ref,
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
        "Start here" => [
            "Five-minute quick start" => "getting_started.md",
            "Executable tutorials" => "tutorials.md",
            "Interactive code generator" => "code_generator.md",
            "Mathematical conventions" => "conventions.md",
        ],
        "Entanglement and separability" => [
            "Separability by example" => "separability_examples.md",
            "Certificate-aware analysis" => "entanglement_backends.md",
            "Product structure and separable balls" => "product_analysis.md",
            "UPB construction catalog" => "upb_catalog.md",
            "UPB certificates and witnesses" => "is_upb.md",
            "Absolute PPT from a spectrum" => "absolute_ppt.md",
            "Symmetric-extension hierarchies" => "symmetric_extensions.md",
            "Bounded random PPT states" => "random_ppt_states.md",
            "Separability and local discrimination" => "separability_optimization.md",
            "Symmetric SAPPT states and witnesses" => "paper_symmetric_separability.md",
            "EntanglementDetection.jl extension" => "entanglement_detection_extension.md",
        ],
        "States, channels, and maps" => [
            "States, operators, and random objects" => "states_operators_random.md",
            "Symmetric multiqubit and multiqudit states" => "symmetric_states.md",
            "General and rectangular maps" => "general_maps.md",
            "Bounded random completely positive maps" => "random_superoperators.md",
            "Guarded group twirls" => "twirls.md",
            "Bounded operator Sinkhorn scaling" => "operator_sinkhorn.md",
            "Bounded filter normal form" => "filter_normal_form.md",
            "Channel norms and discrimination" => "channel_optimization.md",
        ],
        "Coherence and matrix analysis" => [
            "Coherence" => "coherence.md",
            "Pure-state robustness of k-coherence" => "pure_k_coherence_robustness.md",
            "Coherence criteria and optimization" => "coherence_optimization.md",
            "Matrix analysis" => "matrix_analysis.md",
            "Matrix predicates" => "matrix_predicates.md",
            "Induced Schatten norm lower bounds" => "induced_schatten_norm.md",
        ],
        "Optimization and nonlocality" => [
            "Optimization architecture" => "optimization_architecture.md",
            "Optimization and solvers" => "optimization_and_solvers.md",
            "Positive-semidefinite constraints" => "positive_semidefinite_constraints.md",
            "Matsumoto-fidelity SDP models" => "matsumoto_fidelity_model.md",
            "Top-k p-norm model atoms" => "top_k_p_norm_models.md",
            "Dual top-k p-norm model atoms" => "top_k_p_norm_dual_models.md",
            "S(k) norms and block positivity" => "sk_operator_norm_and_block_positivity.md",
            "Solver-free polynomial foundations" => "polynomial_foundations.md",
            "Polynomial SOS hierarchy" => "polynomial_sos.md",
            "Copositivity and clique-number bounds" => "copositivity_and_cliques.md",
            "Minimum-error state discrimination" => "state_discrimination.md",
            "Nonlocal games" => "nonlocal_games.md",
            "Nonlocal optimization and NPA" => "nonlocal_optimization.md",
        ],
        "Package use and maintenance" => [
            "Architecture" => "architecture.md",
            "Migration from QETLAB" => "migration_from_qetlab.md",
            "External integrations" => "integrating_external_packages.md",
            "Sparse and large-scale work" => "sparse_and_large_scale.md",
            "Reproducibility" => "reproducibility.md",
            "Performance" => "performance.md",
            "Maintaining upstream parity" => "maintaining_upstream_parity.md",
        ],
        "API reference" => [
            "Overview" => "api/index.md",
            "Julia-native core and analysis" => "api/native.md",
            "Julia-native channels and maps" => "api/native_channels.md",
            "Julia-native entanglement analysis" => "api/native_entanglement.md",
            "Julia-native optimization models" => "api/native_optimization.md",
            "MATLAB compatibility" => "api/matlab_compat.md",
        ],
    ],
)
