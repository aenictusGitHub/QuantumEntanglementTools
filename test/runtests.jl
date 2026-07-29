using Test
using QuantumEntanglementTools

include("tier_a_subsystem_kernel.jl")
include("tier_b_states_operators_random.jl")
include("tier_c_channels_maps.jl")
include("tier_d_measures_criteria.jl")
include("tier_d_entanglement_pipeline.jl")
include("tier_e_coherence.jl")
include("tier_e_product_analysis.jl")
include("tier_e_product_compat.jl")
include("tier_e_matrix_analysis.jl")
include("tier_e_matrix_analysis_compat.jl")
include("tier_e_matrix_predicates.jl")
include("tier_e_matrix_predicates_compat.jl")
include(joinpath(@__DIR__, "..", "tutorials", "runtests.jl"))
