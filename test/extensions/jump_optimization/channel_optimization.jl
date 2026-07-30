using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const ChannelOptimizationQET = QuantumEntanglementTools

if !isdefined(ChannelOptimizationQET, :ChannelNormResult)
    Base.include(
        ChannelOptimizationQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "channel_optimization.jl"
        ),
    )
end

function channel_optimization_hypatia_backend(; kwargs...)
    return ChannelOptimizationQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=3.0e-7,
        rtol=3.0e-7,
        kwargs...,
    )
end

function channel_optimization_scs_backend(; kwargs...)
    return ChannelOptimizationQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=1.0e-5,
        rtol=1.0e-5,
        kwargs...,
    )
end

function channel_optimization_transpose_map()
    swap = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    return ChannelOptimizationQET.ChoiRepresentation(swap, 2, 2)
end

function channel_optimization_identity()
    return ChannelOptimizationQET.KrausRepresentation([Matrix{ComplexF64}(I, 2, 2)])
end

function channel_optimization_rotation()
    return ChannelOptimizationQET.KrausRepresentation([ComplexF64[0.8 0.6; -0.6 0.8]])
end

function channel_optimization_depolarizing()
    identity = ComplexF64[1 0; 0 1]
    pauli_x = ComplexF64[0 1; 1 0]
    pauli_y = ComplexF64[0 -im; im 0]
    pauli_z = ComplexF64[1 0; 0 -1]
    return ChannelOptimizationQET.KrausRepresentation([
        identity / 2, pauli_x / 2, pauli_y / 2, pauli_z / 2
    ])
end

@testset "WP4 channel optional optimization" begin
    transpose_map = channel_optimization_transpose_map()
    identity_channel = channel_optimization_identity()
    rotation_channel = channel_optimization_rotation()
    depolarizing_channel = channel_optimization_depolarizing()

    @testset "positional extension methods" begin
        backend = channel_optimization_hypatia_backend()
        @test ChannelOptimizationQET.diamond_norm(identity_channel, backend).value == 1
        @test ChannelOptimizationQET.cb_norm(identity_channel, backend).value == 1
        @test ChannelOptimizationQET.channel_distinguishability(
            identity_channel, identity_channel, backend
        ).success_probability == 0.5
        @test ChannelOptimizationQET.maximum_output_fidelity(
            identity_channel, identity_channel, backend
        ).value == 1
    end

    @testset "Hypatia diamond and completely bounded norms" begin
        transpose_norm = ChannelOptimizationQET.diamond_norm(
            transpose_map; backend=channel_optimization_hypatia_backend()
        )
        @test transpose_norm.status ===
            ChannelOptimizationQET.ChannelOptimizationSolverOptimal
        @test transpose_norm.optimization_result.status ===
            ChannelOptimizationQET.OptimizationOptimal
        @test transpose_norm.optimization_result.termination_status === :optimal
        @test transpose_norm.optimization_result.primal_status === :feasible_point
        @test transpose_norm.optimization_result.dual_status === :feasible_point
        @test transpose_norm.value ≈ 2 atol = 5e-7
        @test transpose_norm.lower_bound ≈ 2 atol = 5e-7
        @test transpose_norm.upper_bound ≈ 2 atol = 5e-7
        @test transpose_norm.density_operators !== nothing
        @test transpose_norm.witness_operator !== nothing
        @test transpose_norm.residuals.trace_residual <= 5e-7
        @test transpose_norm.residuals.structure_residual <= 5e-7
        @test transpose_norm.optimization_result.primal_residual <= 5e-7
        @test transpose_norm.optimization_result.dual_residual <= 5e-7
        @test transpose_norm.optimization_result.absolute_gap <= 5e-7
        @test !transpose_norm.certified
        @test transpose_norm.certificate_kind === nothing

        transpose_cb = ChannelOptimizationQET.cb_norm(
            transpose_map; backend=channel_optimization_hypatia_backend()
        )
        @test transpose_cb.status ===
            ChannelOptimizationQET.ChannelOptimizationSolverOptimal
        @test transpose_cb.value ≈ 2 atol = 5e-7
        @test transpose_cb.quantity === :completely_bounded
        @test transpose_cb.certificate.transformation === :hilbert_schmidt_adjoint

        left = ComplexF64[2 0; 0 1; 0 0]
        right = ComplexF64[3 0; 0 1; 0 0]
        elementary = ChannelOptimizationQET.OperatorSumRepresentation([left], [right])
        elementary_norm = ChannelOptimizationQET.diamond_norm(
            elementary; backend=channel_optimization_hypatia_backend()
        )
        @test elementary_norm.status ===
            ChannelOptimizationQET.ChannelOptimizationSolverOptimal
        @test elementary_norm.value ≈ 6 atol = 2e-6
        @test elementary_norm.input_dimension == 2
        @test elementary_norm.output_dimension == 3
        @test elementary_norm.residuals.structure_residual <= 2e-6
    end

    @testset "channel discrimination with corrected probability" begin
        result = ChannelOptimizationQET.channel_distinguishability(
            identity_channel,
            depolarizing_channel;
            backend=channel_optimization_hypatia_backend(),
        )
        # ||identity - completely_depolarizing||_diamond = 3/2 for qubits.
        @test result.status === ChannelOptimizationQET.ChannelOptimizationSolverOptimal
        @test result.success_probability ≈ 7 / 8 atol = 8e-7
        @test result.lower_bound ≈ 7 / 8 atol = 8e-7
        @test result.upper_bound ≈ 7 / 8 atol = 8e-7
        @test result.diamond_norm_result.value ≈ 3 / 4 atol = 8e-7
        @test result.diamond_norm_result.optimization_result.status ===
            ChannelOptimizationQET.OptimizationOptimal
        @test result.diamond_norm_result.optimization_result.primal_status ===
            :feasible_point
        @test result.diamond_norm_result.optimization_result.dual_status === :feasible_point
        @test !result.certified
    end

    @testset "maximum output fidelity SDP" begin
        result = ChannelOptimizationQET.maximum_output_fidelity(
            identity_channel,
            rotation_channel;
            backend=channel_optimization_hypatia_backend(),
        )
        @test result.status === ChannelOptimizationQET.ChannelOptimizationSolverOptimal
        @test result.value ≈ 1 atol = 8e-7
        @test result.lower_bound ≈ 1 atol = 8e-7
        @test result.upper_bound ≈ 1 atol = 8e-7
        @test result.input_states !== nothing
        @test result.output_states !== nothing
        @test result.coupling !== nothing
        @test result.residuals.trace_residual <= 8e-7
        @test result.residuals.structure_residual <= 8e-7
        @test result.optimization_result.status ===
            ChannelOptimizationQET.OptimizationOptimal
        @test result.optimization_result.primal_residual <= 8e-7
        @test result.optimization_result.dual_residual <= 8e-7
        @test !result.certified
    end

    @testset "SCS cross-checks" begin
        transpose_norm = ChannelOptimizationQET.diamond_norm(
            transpose_map; backend=channel_optimization_scs_backend()
        )
        @test transpose_norm.status in (
            ChannelOptimizationQET.ChannelOptimizationSolverOptimal,
            ChannelOptimizationQET.ChannelOptimizationSolverFeasible,
        )
        @test transpose_norm.lower_bound ≈ 2 atol = 4e-5
        @test transpose_norm.witness_operator !== nothing
        @test transpose_norm.optimization_result.optimizer.reported_optimizer_name == "SCS"
        @test transpose_norm.optimization_result.termination_status in
            (:optimal, :almost_optimal)

        fidelity_result = ChannelOptimizationQET.maximum_output_fidelity(
            identity_channel, rotation_channel; backend=channel_optimization_scs_backend()
        )
        @test fidelity_result.status in (
            ChannelOptimizationQET.ChannelOptimizationSolverOptimal,
            ChannelOptimizationQET.ChannelOptimizationSolverFeasible,
        )
        @test fidelity_result.lower_bound ≈ 1 atol = 5e-5
        @test fidelity_result.output_states !== nothing
        @test fidelity_result.optimization_result.optimizer.reported_optimizer_name == "SCS"
    end

    @testset "limit and malformed backend statuses" begin
        limited = ChannelOptimizationQET.diamond_norm(
            transpose_map;
            backend=channel_optimization_hypatia_backend(optimizer_options=(iter_limit=0,)),
        )
        @test limited.status === ChannelOptimizationQET.ChannelOptimizationResourceLimit
        @test limited.optimization_result.status ===
            ChannelOptimizationQET.OptimizationLimit
        @test limited.optimization_result.termination_status === :iteration_limit
        @test limited.value === nothing

        malformed_backend = ChannelOptimizationQET.JuMPBackend(
            () -> error("intentional channel-optimization factory failure");
            optimizer_name="intentional malformed factory",
            allow_densify=true,
        )
        malformed_norm = ChannelOptimizationQET.diamond_norm(
            transpose_map; backend=malformed_backend
        )
        @test malformed_norm.status ===
            ChannelOptimizationQET.ChannelOptimizationBackendFailure
        @test malformed_norm.optimization_result.status ===
            ChannelOptimizationQET.OptimizationMalformedBackend
        @test malformed_norm.optimization_result.termination_status === :backend_exception
        @test malformed_norm.value === nothing
        @test malformed_norm.witness_operator === nothing

        malformed_fidelity = ChannelOptimizationQET.maximum_output_fidelity(
            identity_channel, rotation_channel; backend=malformed_backend
        )
        @test malformed_fidelity.status ===
            ChannelOptimizationQET.ChannelOptimizationBackendFailure
        @test malformed_fidelity.optimization_result.status ===
            ChannelOptimizationQET.OptimizationMalformedBackend
        @test malformed_fidelity.value === nothing
        @test malformed_fidelity.coupling === nothing
    end
end
