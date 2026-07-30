using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const KpEpigraphOptQET = QuantumEntanglementTools
const KpEpigraphMOI = JuMP.MOI

if !isdefined(KpEpigraphOptQET, :TopKPNormEpigraph)
    Base.include(
        KpEpigraphOptQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "top_k_p_norm_epigraph.jl"
        ),
    )
end

const KpEpigraphExtension = Base.get_extension(
    KpEpigraphOptQET, :QuantumEntanglementToolsJuMPExt
)
isnothing(KpEpigraphExtension) && error("QuantumEntanglementToolsJuMPExt is not loaded")
if !isdefined(KpEpigraphExtension, :_top_k_p_hermitian_real_block)
    Base.include(
        KpEpigraphExtension,
        joinpath(@__DIR__, "..", "..", "..", "ext", "top_k_p_norm_epigraph.jl"),
    )
end

function fixed_kp_epigraph_model(
    matrix, k, p, optimizer_factory; optimizer_options=NamedTuple()
)
    rows, columns = size(matrix)
    affine = KpEpigraphOptQET.complex_affine_variable(
        :X, rows, columns; coefficient_type=Float64
    )
    atom = KpEpigraphOptQET.top_k_p_norm_epigraph(affine, k, p)
    model = JuMP.Model(optimizer_factory)
    JuMP.set_silent(model)
    for (name, value) in pairs(optimizer_options)
        JuMP.set_optimizer_attribute(model, String(name), value)
    end
    coordinates = JuMP.@variable(model, [1:(affine.variable_count)])
    handle = KpEpigraphOptQET.add_top_k_p_norm_epigraph!(
        model, atom, coordinates; allow_densify=true
    )
    target = vcat(vec(real.(matrix)), vec(imag.(matrix)))
    fixing = [
        JuMP.@constraint(model, coordinates[index] == target[index]) for
        index in eachindex(coordinates)
    ]
    JuMP.@objective(model, Min, handle.epigraph)
    return (; model, affine, atom, coordinates, handle, fixing)
end

@testset "WP3 top-k p-norm optional epigraph" begin
    matrix = ComplexF64[3 0; 0 1]
    expected_p3 = KpEpigraphOptQET.top_k_p_norm(matrix, 2, 3)
    hypatia = fixed_kp_epigraph_model(matrix, 2, 3, Hypatia.Optimizer)
    JuMP.optimize!(hypatia.model)
    @test JuMP.termination_status(hypatia.model) === KpEpigraphMOI.OPTIMAL
    @test JuMP.primal_status(hypatia.model) === KpEpigraphMOI.FEASIBLE_POINT
    @test JuMP.dual_status(hypatia.model) === KpEpigraphMOI.FEASIBLE_POINT
    @test JuMP.objective_value(hypatia.model) ≈ expected_p3 atol = 3e-7
    @test JuMP.value(hypatia.handle.epigraph) ≈ expected_p3 atol = 3e-7
    @test hypatia.handle.certified_epigraph
    @test hypatia.handle.formulation === :ky_fan_majorization_and_norm_cone
    @test length(hypatia.handle.constraints.psd) == 4
    @test length(hypatia.handle.constraints.majorization) == 2
    @test all(value -> value >= -3e-7, JuMP.value.(hypatia.handle.singular_value_majorant))

    nuclear = fixed_kp_epigraph_model(matrix, 2, 1, Hypatia.Optimizer)
    JuMP.optimize!(nuclear.model)
    @test JuMP.termination_status(nuclear.model) === KpEpigraphMOI.OPTIMAL
    @test JuMP.objective_value(nuclear.model) ≈ KpEpigraphOptQET.top_k_p_norm(matrix, 2, 1) atol =
        3e-7

    vector = reshape(ComplexF64[3 + 4im, -2, 1im], 1, :)
    vector_model = fixed_kp_epigraph_model(vector, 2, Inf, Hypatia.Optimizer)
    JuMP.optimize!(vector_model.model)
    @test JuMP.termination_status(vector_model.model) === KpEpigraphMOI.OPTIMAL
    @test JuMP.objective_value(vector_model.model) ≈
        KpEpigraphOptQET.top_k_p_norm(vec(vector), 2, Inf) atol = 4e-7

    scs = fixed_kp_epigraph_model(
        matrix, 2, 3, SCS.Optimizer; optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7)
    )
    JuMP.optimize!(scs.model)
    @test JuMP.termination_status(scs.model) in
        (KpEpigraphMOI.OPTIMAL, KpEpigraphMOI.ALMOST_OPTIMAL)
    @test JuMP.objective_value(scs.model) ≈ expected_p3 atol = 5e-5

    affine = KpEpigraphOptQET.complex_affine_variable(:small, 1, 2)
    atom = KpEpigraphOptQET.top_k_p_norm_epigraph(affine, 1, 2)
    empty_model = JuMP.Model(Hypatia.Optimizer)
    coordinates = JuMP.@variable(empty_model, [1:(affine.variable_count)])
    @test_throws ArgumentError KpEpigraphOptQET.add_top_k_p_norm_epigraph!(
        empty_model, atom, coordinates
    )
    @test_throws DimensionMismatch KpEpigraphOptQET.add_top_k_p_norm_epigraph!(
        empty_model, atom, coordinates[1:(end - 1)]; allow_densify=true
    )

    limited = fixed_kp_epigraph_model(
        matrix, 2, 3, Hypatia.Optimizer; optimizer_options=(iter_limit=0,)
    )
    JuMP.optimize!(limited.model)
    @test JuMP.termination_status(limited.model) === KpEpigraphMOI.ITERATION_LIMIT

    @test_throws ErrorException JuMP.Model(
        () -> error("intentional top-k p-norm optimizer factory failure")
    )
end
