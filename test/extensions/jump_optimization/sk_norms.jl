using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SCS
using Test

const SKNormOptQET = QuantumEntanglementTools
const SKNormMOI = JuMP.MOI

if !isdefined(SKNormOptQET, :TopKPNormDualEpigraph)
    Base.include(
        SKNormOptQET,
        joinpath(
            @__DIR__,
            "..",
            "..",
            "..",
            "src",
            "optimization",
            "top_k_p_norm_dual_epigraph.jl",
        ),
    )
end
if !isdefined(SKNormOptQET, :SKOperatorNormResult)
    Base.include(
        SKNormOptQET,
        joinpath(@__DIR__, "..", "..", "..", "src", "optimization", "sk_operator_norm.jl"),
    )
end
if !isdefined(SKNormOptQET, :BlockPositivityResult)
    Base.include(
        SKNormOptQET,
        joinpath(@__DIR__, "..", "..", "..", "src", "entanglement", "block_positivity.jl"),
    )
end

const SKNormExtension = Base.get_extension(SKNormOptQET, :QuantumEntanglementToolsJuMPExt)
isnothing(SKNormExtension) && error("QuantumEntanglementToolsJuMPExt is not loaded")
if !isdefined(SKNormExtension, :_top_k_p_dual_singular_majorant!)
    Base.include(
        SKNormExtension,
        joinpath(@__DIR__, "..", "..", "..", "ext", "top_k_p_norm_dual_epigraph.jl"),
    )
end
if !hasmethod(
    SKNormOptQET.sk_operator_norm, Tuple{Random.AbstractRNG,Any,SKNormOptQET.JuMPBackend}
)
    Base.include(
        SKNormExtension, joinpath(@__DIR__, "..", "..", "..", "ext", "sk_operator_norm.jl")
    )
end
if !hasmethod(
    SKNormOptQET.is_block_positive, Tuple{Random.AbstractRNG,Any,SKNormOptQET.JuMPBackend}
)
    Base.include(
        SKNormExtension, joinpath(@__DIR__, "..", "..", "..", "ext", "block_positivity.jl")
    )
end

function sknorm_hypatia_backend(; options=NamedTuple(), kwargs...)
    return SKNormOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_options=options,
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        kwargs...,
    )
end

function sknorm_scs_backend(; options=(eps_abs=1.0e-7, eps_rel=1.0e-7), kwargs...)
    return SKNormOptQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=options,
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=5.0e-6,
        rtol=5.0e-6,
        kwargs...,
    )
end

function fixed_dual_epigraph_model(
    matrix, k, p, optimizer_factory; optimizer_options=NamedTuple()
)
    rows, columns = size(matrix)
    affine = SKNormOptQET.complex_affine_variable(:X, rows, columns)
    atom = SKNormOptQET.top_k_p_norm_dual_epigraph(affine, k, p)
    model = JuMP.Model(optimizer_factory)
    JuMP.set_silent(model)
    for (name, value) in pairs(optimizer_options)
        JuMP.set_optimizer_attribute(model, String(name), value)
    end
    coordinates = JuMP.@variable(model, [1:(affine.variable_count)])
    handle = SKNormOptQET.add_top_k_p_norm_dual_epigraph!(
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

@testset "WP7 optional top-k p-norm dual epigraph" begin
    matrix = ComplexF64[3 0; 0 1]
    for p in (1, 2, 3, Inf)
        expected = SKNormOptQET.top_k_p_norm_dual(matrix, 2, p)
        built = fixed_dual_epigraph_model(matrix, 2, p, Hypatia.Optimizer)
        JuMP.optimize!(built.model)
        @test JuMP.termination_status(built.model) === SKNormMOI.OPTIMAL
        @test JuMP.primal_status(built.model) === SKNormMOI.FEASIBLE_POINT
        @test JuMP.dual_status(built.model) === SKNormMOI.FEASIBLE_POINT
        @test JuMP.objective_value(built.model) ≈ expected atol = 3e-7
        @test JuMP.value(built.handle.epigraph) ≈ expected atol = 3e-7
        @test built.handle.certified_epigraph
    end

    k_support = fixed_dual_epigraph_model(matrix, 1, 2, Hypatia.Optimizer)
    JuMP.optimize!(k_support.model)
    @test JuMP.objective_value(k_support.model) ≈ 4 atol = 4e-7
    @test k_support.handle.formulation === :k_support_perspective_power_cones
    @test length(k_support.handle.constraints.psd) == 4
    @test length(k_support.handle.constraints.majorization) == 2
    @test length(k_support.handle.constraints.gauge) == 6

    vector = reshape(ComplexF64[3 + 4im, -2, 1im], 1, :)
    vector_built = fixed_dual_epigraph_model(vector, 2, 2, Hypatia.Optimizer)
    JuMP.optimize!(vector_built.model)
    @test JuMP.objective_value(vector_built.model) ≈
        SKNormOptQET.top_k_p_norm_dual(vec(vector), 2, 2) atol = 5e-7

    scs = fixed_dual_epigraph_model(
        matrix, 2, 2, SCS.Optimizer; optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7)
    )
    JuMP.optimize!(scs.model)
    @test JuMP.termination_status(scs.model) in
        (SKNormMOI.OPTIMAL, SKNormMOI.ALMOST_OPTIMAL)
    @test JuMP.objective_value(scs.model) ≈ SKNormOptQET.top_k_p_norm_dual(matrix, 2, 2) atol =
        8e-5

    affine = SKNormOptQET.complex_affine_variable(:small, 1, 2)
    atom = SKNormOptQET.top_k_p_norm_dual_epigraph(affine, 1, 2)
    empty_model = JuMP.Model(Hypatia.Optimizer)
    coordinates = JuMP.@variable(empty_model, [1:(affine.variable_count)])
    @test_throws ArgumentError SKNormOptQET.add_top_k_p_norm_dual_epigraph!(
        empty_model, atom, coordinates
    )
    @test_throws DimensionMismatch SKNormOptQET.add_top_k_p_norm_dual_epigraph!(
        empty_model, atom, coordinates[1:(end - 1)]; allow_densify=true
    )
end

@testset "WP7 optional S(k) SDP relaxations" begin
    diagonal = Diagonal([4.0, 3.0, 2.0, 1.0])
    hypatia = sknorm_hypatia_backend()
    result = SKNormOptQET.sk_operator_norm(
        MersenneTwister(721),
        diagonal,
        hypatia;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test result.optimization_result.status === SKNormOptQET.OptimizationOptimal
    @test result.relaxation_upper_bound !== nothing
    @test result.relaxation_upper_bound ≈ 4 atol = 4e-6
    @test result.lower_bound == 4
    @test result.upper_witness.numerically_validated
    @test result.upper_witness.formulation === :ppt_relaxation
    @test result.optimization_result.primal_residual <= result.tolerance
    @test result.optimization_result.dual_residual <= result.tolerance

    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    positive = bell * bell' + 0.2Matrix{ComplexF64}(I, 4, 4)
    exact_ppt = SKNormOptQET.sk_operator_norm(
        MersenneTwister(722),
        positive,
        hypatia;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test exact_ppt.optimization_result.status === SKNormOptQET.OptimizationOptimal
    @test exact_ppt.lower_bound <= 0.7 + 5e-7
    @test exact_ppt.upper_bound ≈ 0.7 atol = 5e-6

    hierarchy_problem = SKNormOptQET.sk_operator_norm_problem(
        Matrix{Float64}(I, 4, 4); k=1, dims=(2, 2), hierarchy_level=2, allow_densify=true
    )
    hierarchy_solve = SKNormOptQET.solve_optimization(hierarchy_problem.program, hypatia)
    @test hierarchy_solve.status === SKNormOptQET.OptimizationOptimal
    @test hierarchy_solve.objective_value ≈ 1 atol = 3e-7
    @test hierarchy_solve.primal_residual <= 3e-7
    @test hierarchy_solve.dual_residual <= 3e-7

    scs = sknorm_scs_backend()
    scs_result = SKNormOptQET.sk_operator_norm(
        MersenneTwister(723),
        diagonal,
        scs;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        rtol=1e-5,
        allow_densify=true,
    )
    @test scs_result.optimization_result.status in
        (SKNormOptQET.OptimizationOptimal, SKNormOptQET.OptimizationFeasible)
    @test scs_result.optimization_result.objective_value ≈ 4 atol = 2e-4

    limited = sknorm_hypatia_backend(options=(iter_limit=0,))
    limited_result = SKNormOptQET.sk_operator_norm(
        MersenneTwister(724),
        positive,
        limited;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test limited_result.optimization_result.status === SKNormOptQET.OptimizationLimit
    @test limited_result.relaxation_upper_bound === nothing
    @test limited_result.upper_bound == opnorm(positive)

    broken = SKNormOptQET.JuMPBackend(
        () -> error("intentional S(k) optimizer factory failure");
        optimizer_name="broken",
        allow_densify=true,
    )
    broken_result = SKNormOptQET.sk_operator_norm(
        MersenneTwister(725),
        positive,
        broken;
        k=1,
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test broken_result.optimization_result.status ===
        SKNormOptQET.OptimizationMalformedBackend
    @test broken_result.relaxation_upper_bound === nothing
end

@testset "WP7 optional block-positivity upper certificate" begin
    swap = ComplexF64[
        1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1
    ]
    strictly_block_positive = swap + 0.1Matrix{ComplexF64}(I, 4, 4)
    result = SKNormOptQET.is_block_positive(
        MersenneTwister(726),
        strictly_block_positive,
        sknorm_hypatia_backend();
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test result.status === SKNormOptQET.BlockPositivityCertified
    @test result.verdict === true
    @test result.certified
    @test result.certificate_kind === :sk_sdp_upper_relaxation
    @test result.upper_witness.numerically_validated
    @test result.sk_norm_result.optimization_result.status ===
        SKNormOptQET.OptimizationOptimal

    inconclusive = SKNormOptQET.is_block_positive(
        MersenneTwister(727),
        swap,
        sknorm_hypatia_backend();
        dims=(2, 2),
        strength=1,
        max_restarts=0,
        allow_densify=true,
    )
    @test inconclusive.verdict !== false
    @test inconclusive.status in (
        SKNormOptQET.BlockPositivityUnknown,
        SKNormOptQET.BlockPositivityNumericalBoundary,
        SKNormOptQET.BlockPositivityCertified,
    )
end
