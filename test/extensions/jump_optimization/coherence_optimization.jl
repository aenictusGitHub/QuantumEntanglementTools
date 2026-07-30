#!/usr/bin/env julia

using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const QETCoherenceSDP = QuantumEntanglementTools
if !isdefined(QETCoherenceSDP, :is_k_incoherent)
    Base.include(
        QETCoherenceSDP,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "coherence", "coherence_optimization.jl"
        ),
    )
end

function _coherence_hypatia_backend(; options=NamedTuple())
    return QETCoherenceSDP.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_options=options,
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=1.0e-7,
        rtol=1.0e-7,
    )
end

function _coherence_scs_backend()
    return QETCoherenceSDP.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=1.0e-6, eps_rel=1.0e-6),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=2.0e-5,
        rtol=2.0e-5,
    )
end

@testset "coherence optimization JuMP extension" begin
    hypatia = _coherence_hypatia_backend()
    scs = _coherence_scs_backend()

    @testset "analytic certificates agree with both SDP backends" begin
        state = fill(0.5 + 0im, 4)
        rho = state * adjoint(state)
        analytic_robustness = QETCoherenceSDP.robustness_coherence(state)
        analytic_trace = QETCoherenceSDP.trace_distance_coherence(state)
        analytic_generalized = QETCoherenceSDP.generalized_robustness_k_coherence(state, 2)

        for backend in (hypatia, scs)
            robustness = QETCoherenceSDP.robustness_coherence(rho; strategy=:sdp, backend)
            trace = QETCoherenceSDP.trace_distance_coherence(rho; strategy=:sdp, backend)
            generalized = QETCoherenceSDP.generalized_robustness_k_coherence(
                rho, 2; strategy=:sdp, backend
            )
            for result in (robustness, trace, generalized)
                @test result.status === :optimal
                @test result.optimization.status === QETCoherenceSDP.OptimizationOptimal
                @test result.optimization.termination_status === :optimal
                @test result.optimization.primal_status === :feasible_point
                @test result.optimization.dual_status === :feasible_point
                @test result.optimization.primal_residual ≤ 3.0e-5
                @test result.optimization.dual_residual ≤ 3.0e-5
                @test !result.exact
            end
            @test robustness.value ≈ analytic_robustness.value atol = 4.0e-5
            @test trace.value ≈ analytic_trace.value atol = 4.0e-5
            @test generalized.value ≈ analytic_generalized.value atol = 4.0e-5
            @test real(tr(robustness.free_state)) ≈ 1 atol = 4.0e-5
            @test real(tr(generalized.free_state)) ≈ 1 atol = 4.0e-5
            @test minimum(eigvals(Hermitian(robustness.unnormalized_noise))) ≥ -4.0e-5
            @test minimum(eigvals(Hermitian(generalized.unnormalized_noise))) ≥ -4.0e-5

            incoherent = QETCoherenceSDP.is_k_incoherent(rho, 2; strategy=:sdp, backend)
            @test incoherent.status === :solver_infeasible
            @test incoherent.verdict === false
            @test incoherent.optimization.status === QETCoherenceSDP.OptimizationInfeasible
            @test !incoherent.exact
        end
    end

    @testset "general mixed-state Hypatia/SCS cross-check" begin
        rho = ComplexF64[
            0.40 0.05 0.02im
            0.05 0.35 0.03
            -0.02im 0.03 0.25
        ]
        hypatia_robustness = QETCoherenceSDP.robustness_coherence(
            rho; strategy=:sdp, backend=hypatia
        )
        scs_robustness = QETCoherenceSDP.robustness_coherence(
            rho; strategy=:sdp, backend=scs
        )
        hypatia_trace = QETCoherenceSDP.trace_distance_coherence(
            rho; strategy=:sdp, backend=hypatia
        )
        scs_trace = QETCoherenceSDP.trace_distance_coherence(
            rho; strategy=:sdp, backend=scs
        )
        @test hypatia_robustness.status === :optimal
        @test scs_robustness.status === :optimal
        @test hypatia_trace.status === :optimal
        @test scs_trace.status === :optimal
        @test hypatia_robustness.value ≈ scs_robustness.value atol = 4.0e-5
        @test hypatia_trace.value ≈ scs_trace.value atol = 4.0e-5
        @test hypatia_robustness.value ≈ 0.17689600744 atol = 2.0e-7
        @test hypatia_trace.value ≈ 0.12328828059 atol = 2.0e-7

        free_hypatia = QETCoherenceSDP.generalized_robustness_k_coherence(
            rho, 2; strategy=:sdp, backend=hypatia
        )
        free_scs = QETCoherenceSDP.generalized_robustness_k_coherence(
            rho, 2; strategy=:sdp, backend=scs
        )
        @test free_hypatia.value ≤ 3.0e-7
        @test free_scs.value ≤ 3.0e-5
        @test free_hypatia.noise_state === nothing
        @test free_scs.noise_state === nothing

        feasible_hypatia = QETCoherenceSDP.is_k_incoherent(
            rho, 2; strategy=:sdp, backend=hypatia
        )
        feasible_scs = QETCoherenceSDP.is_k_incoherent(rho, 2; strategy=:sdp, backend=scs)
        @test feasible_hypatia.verdict === true
        @test feasible_scs.verdict === true
        @test feasible_hypatia.status === :solver_feasible
        @test feasible_scs.status === :solver_feasible
        @test length(feasible_hypatia.decomposition.blocks) == 3
    end

    @testset "absolute (d-1)-incoherence theorem SDP" begin
        absolutely_free = Diagonal(ComplexF64[0.40, 0.35, 0.25]) |> Matrix
        not_absolutely_free = Diagonal(ComplexF64[1, 0, 0]) |> Matrix
        for backend in (hypatia, scs)
            positive = QETCoherenceSDP.is_absolutely_k_incoherent(
                absolutely_free, 2; strategy=:sdp, backend
            )
            negative = QETCoherenceSDP.is_absolutely_k_incoherent(
                not_absolutely_free, 2; strategy=:sdp, backend
            )
            @test positive.status === :solver_feasible
            @test positive.verdict === true
            @test positive.decomposition.theorem_matrix isa Matrix
            @test minimum(eigvals(Hermitian(positive.decomposition.theorem_matrix))) ≥
                -4.0e-5
            @test negative.status === :solver_infeasible
            @test negative.verdict === false
        end
    end

    @testset "backend limits and malformed factories stay explicit" begin
        rho = ComplexF64[
            0.40 0.05 0.02im
            0.05 0.35 0.03
            -0.02im 0.03 0.25
        ]
        limited = QETCoherenceSDP.robustness_coherence(
            rho; strategy=:sdp, backend=_coherence_hypatia_backend(options=(iter_limit=0,))
        )
        @test limited.status === :solver_limit
        @test limited.value === nothing
        @test limited.optimization.status === QETCoherenceSDP.OptimizationLimit

        malformed_backend = QETCoherenceSDP.JuMPBackend(
            () -> 42; optimizer_name="malformed test factory", allow_densify=true
        )
        malformed = QETCoherenceSDP.trace_distance_coherence(
            rho; strategy=:sdp, backend=malformed_backend
        )
        @test malformed.status === :malformed_backend
        @test malformed.value === nothing
        @test malformed.optimization.status === QETCoherenceSDP.OptimizationMalformedBackend
        @test occursin("failed", malformed.optimization.message)
    end
end
