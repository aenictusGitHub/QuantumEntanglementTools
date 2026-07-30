using Hypatia
using JuMP
using QuantumEntanglementTools
using SCS
using Test

const AbsPPTOptQET = QuantumEntanglementTools

if !isdefined(AbsPPTOptQET, :is_abs_ppt)
    Base.include(
        AbsPPTOptQET,
        joinpath(@__DIR__, "..", "..", "..", "src", "entanglement", "absolute_ppt.jl"),
    )
end

function _abs_ppt_hypatia_backend()
    return AbsPPTOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=1.0e-7,
        rtol=1.0e-7,
    )
end

function _abs_ppt_scs_backend()
    return AbsPPTOptQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=5.0e-6,
        rtol=5.0e-6,
    )
end

@testset "absolute-PPT ordering-realization backends" begin
    spectrum = [1.0, 0.0, 0.0, 0.0]
    hypatia_result = AbsPPTOptQET.is_abs_ppt(
        spectrum;
        dims=(2, 2),
        max_realization_iterations=0,
        backend=_abs_ppt_hypatia_backend(),
    )
    @test hypatia_result.status === AbsPPTOptQET.AbsolutePPTCertifiedNot
    @test hypatia_result.verdict === false
    @test hypatia_result.backend_result.status === AbsPPTOptQET.OptimizationOptimal
    @test hypatia_result.backend_result.termination_status === :optimal
    @test hypatia_result.backend_result.primal_status === :feasible_point
    @test hypatia_result.backend_result.primal !== nothing
    @test hypatia_result.backend_result.primal_residual <= 1.0e-7
    @test hypatia_result.backend_result.optimizer.configured_optimizer_name == "Hypatia"
    @test hypatia_result.ordering_certificate.source === :verified_backend_primal
    @test hypatia_result.ordering_certificate.minimum_gap > 0
    @test hypatia_result.ordering_certificate.minimum_product_margin > 0

    scs_result = AbsPPTOptQET.is_abs_ppt(
        spectrum; dims=(2, 2), max_realization_iterations=0, backend=_abs_ppt_scs_backend()
    )
    @test scs_result.status === AbsPPTOptQET.AbsolutePPTCertifiedNot
    @test scs_result.verdict === false
    @test scs_result.backend_result.status in
        (AbsPPTOptQET.OptimizationOptimal, AbsPPTOptQET.OptimizationFeasible)
    @test scs_result.backend_result.primal !== nothing
    @test scs_result.backend_result.primal_residual <= 5.0e-6
    @test scs_result.ordering_certificate.minimum_product_margin > 0

    malformed_backend = AbsPPTOptQET.JuMPBackend(
        () -> error("intentional malformed optimizer factory");
        optimizer_name="intentional malformed factory",
    )
    malformed = AbsPPTOptQET.is_abs_ppt(
        spectrum; dims=(2, 2), max_realization_iterations=0, backend=malformed_backend
    )
    @test malformed.status === AbsPPTOptQET.AbsolutePPTBackendFailure
    @test malformed.verdict === nothing
    @test malformed.backend_result.status === AbsPPTOptQET.OptimizationMalformedBackend
    @test malformed.backend_result.termination_status === :backend_exception
    @test malformed.backend_result.primal === nothing
    @test malformed.violating_constraint !== nothing
end
