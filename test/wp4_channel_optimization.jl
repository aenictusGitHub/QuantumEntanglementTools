using LinearAlgebra
using SparseArrays
using Test

using QuantumEntanglementTools: QuantumEntanglementTools
const QETChannelOptimization = QuantumEntanglementTools

if !isdefined(QETChannelOptimization, :ChannelNormResult)
    Base.include(
        QETChannelOptimization,
        joinpath(@__DIR__, "..", "src", "optimization", "channel_optimization.jl"),
    )
end

function channel_opt_identity(::Type{T}=ComplexF64) where {T}
    return QETChannelOptimization.KrausRepresentation([Matrix{T}(I, 2, 2)])
end

function channel_opt_dephasing(::Type{T}=ComplexF64) where {T}
    return QETChannelOptimization.KrausRepresentation([T[1 0; 0 0], T[0 0; 0 1]])
end

function channel_opt_rotation(::Type{T}=ComplexF64) where {T}
    return QETChannelOptimization.KrausRepresentation([T[0.8 0.6; -0.6 0.8]])
end

function channel_opt_replacer_zero(::Type{T}=ComplexF64) where {T}
    return QETChannelOptimization.KrausRepresentation([T[1 0; 0 0], T[0 1; 0 0]])
end

function channel_opt_replacer_one(::Type{T}=ComplexF64) where {T}
    return QETChannelOptimization.KrausRepresentation([T[0 0; 1 0], T[0 0; 0 1]])
end

function channel_opt_transpose_map(::Type{T}=ComplexF64) where {T}
    swap = T[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    return QETChannelOptimization.ChoiRepresentation(swap, 2, 2)
end

@testset "WP4 channel optimization" begin
    identity_channel = channel_opt_identity()
    dephasing_channel = channel_opt_dephasing()
    rotation_channel = channel_opt_rotation()
    reset_zero = channel_opt_replacer_zero()
    reset_one = channel_opt_replacer_one()
    transpose_map = channel_opt_transpose_map()

    @testset "solver-free channel norms" begin
        identity_norm = QETChannelOptimization.diamond_norm(identity_channel)
        @test identity_norm isa QETChannelOptimization.ChannelNormResult{Float64}
        @test identity_norm.quantity === :diamond
        @test identity_norm.status ===
            QETChannelOptimization.ChannelOptimizationAnalyticOptimal
        @test identity_norm.value == 1
        @test identity_norm.lower_bound == identity_norm.upper_bound == 1
        @test identity_norm.certified
        @test identity_norm.certificate_kind === :completely_positive_adjoint_identity
        @test identity_norm.problem === nothing
        @test identity_norm.optimization_result === nothing

        scaled = QETChannelOptimization.KrausRepresentation([ComplexF64[2 0; 0 2]])
        scaled_norm = QETChannelOptimization.diamond_norm(scaled)
        @test scaled_norm.value == 4
        @test scaled_norm.certificate.adjoint_identity_spectrum == [4, 4]

        reset_cb = QETChannelOptimization.cb_norm(reset_zero)
        @test reset_cb.quantity === :completely_bounded
        @test reset_cb.value == 2
        @test reset_cb.input_dimension == 2
        @test reset_cb.output_dimension == 2
        @test reset_cb.certified
        @test reset_cb.certificate.transformation === :hilbert_schmidt_adjoint

        zero_map = QETChannelOptimization.ChoiRepresentation(
            zeros(ComplexF64, 6, 6), QETChannelOptimization.OperatorSpace(2, 2, 3, 3)
        )
        zero_norm = QETChannelOptimization.diamond_norm(zero_map)
        @test zero_norm.value == 0
        @test zero_norm.certificate_kind === :zero_map
        @test zero_norm.input_dimension == 2
        @test zero_norm.output_dimension == 3

        scalar_input = QETChannelOptimization.ChoiRepresentation(
            ComplexF64[3 0; 0 -2], QETChannelOptimization.OperatorSpace(1, 1, 2, 2)
        )
        scalar_input_norm = QETChannelOptimization.diamond_norm(scalar_input)
        @test scalar_input_norm.value == 5
        @test scalar_input_norm.certificate_kind === :scalar_input_trace_norm

        scalar_output = QETChannelOptimization.ChoiRepresentation(
            ComplexF64[1 2; 0 0], QETChannelOptimization.OperatorSpace(2, 2, 1, 1)
        )
        scalar_output_norm = QETChannelOptimization.diamond_norm(scalar_output)
        @test scalar_output_norm.value ≈ sqrt(5)
        @test scalar_output_norm.certificate_kind === :scalar_output_operator_norm

        identity32 = channel_opt_identity(ComplexF32)
        result32 = QETChannelOptimization.diamond_norm(identity32)
        @test result32 isa QETChannelOptimization.ChannelNormResult{Float32}
        @test result32.value isa Float32
        @test result32.value == 1

        sparse_identity = QETChannelOptimization.KrausRepresentation([
            sparse(Matrix{ComplexF64}(I, 2, 2))
        ])
        sparse_without_permission = QETChannelOptimization.diamond_norm(sparse_identity)
        @test sparse_without_permission.status ===
            QETChannelOptimization.ChannelOptimizationBackendUnavailable
        @test sparse_without_permission.value === nothing
        @test sparse_without_permission.problem !== nothing
        sparse_with_permission = QETChannelOptimization.diamond_norm(
            sparse_identity; allow_densify=true
        )
        @test sparse_with_permission.value == 1
    end

    @testset "diamond-norm affine model and no-backend status" begin
        problem = QETChannelOptimization.diamond_norm_problem(transpose_map)
        @test problem isa QETChannelOptimization.DiamondNormProblem
        @test problem.input_dimension == 2
        @test problem.output_dimension == 2
        @test problem.program.name === :diamond_norm_watrous_primal
        @test problem.program.sense === :maximize
        @test problem.program.variable_count == 72
        @test length(problem.program.equalities) == 34
        @test length(problem.program.psd_constraints) == 3
        @test length(problem.program.primal_views) == 3
        @test problem.program.known_feasible_point !== nothing
        @test QETChannelOptimization.primal_residual(
            problem.program, problem.program.known_feasible_point; allow_densify=true
        ) == 0
        @test QETChannelOptimization.evaluate_affine(
            problem.program.objective, problem.program.known_feasible_point
        ) == 0

        unavailable = QETChannelOptimization.diamond_norm(transpose_map)
        @test unavailable.status ===
            QETChannelOptimization.ChannelOptimizationBackendUnavailable
        @test unavailable.value === nothing
        @test unavailable.lower_bound == 0
        @test unavailable.upper_bound === nothing
        @test unavailable.problem isa QETChannelOptimization.DiamondNormProblem
        @test unavailable.optimization_result.status ===
            QETChannelOptimization.OptimizationBackendUnavailable
        @test !unavailable.certified

        left = ComplexF64[2 0; 0 1; 0 0]
        right = ComplexF64[3 0; 0 1; 0 0]
        rectangular_two_sided = QETChannelOptimization.OperatorSumRepresentation(
            [left], [right]
        )
        rectangular_problem = QETChannelOptimization.diamond_norm_problem(
            rectangular_two_sided
        )
        @test rectangular_problem.input_dimension == 2
        @test rectangular_problem.output_dimension == 3
        @test rectangular_problem.program.variable_count == 152
        @test length(rectangular_problem.program.equalities) == 74
        @test QETChannelOptimization.primal_residual(
            rectangular_problem.program,
            rectangular_problem.program.known_feasible_point;
            allow_densify=true,
        ) == 0

        @test_throws ArgumentError QETChannelOptimization.diamond_norm_problem(
            transpose_map; max_variables=71
        )
        @test_throws ArgumentError QETChannelOptimization.diamond_norm(
            identity_channel; max_choi_dimension=3
        )
        @test_throws ArgumentError QETChannelOptimization.diamond_norm_problem(
            transpose_map;
            limits=QETChannelOptimization.OptimizationLimits(max_psd_dimension=7),
        )

        nonsquare_space = QETChannelOptimization.OperatorSpace(2, 3, 2, 2)
        nonsquare_map = QETChannelOptimization.ChoiRepresentation(
            zeros(ComplexF64, 4, 6), nonsquare_space
        )
        @test_throws ArgumentError QETChannelOptimization.diamond_norm(nonsquare_map)
    end

    @testset "channel Holevo--Helstrom reduction" begin
        identical = QETChannelOptimization.channel_distinguishability(
            identity_channel, identity_channel; priors=[0.8, 0.2]
        )
        @test identical.status === QETChannelOptimization.ChannelOptimizationAnalyticOptimal
        @test identical.success_probability == 0.8
        @test identical.lower_bound == identical.upper_bound == 0.8
        @test identical.certified
        @test identical.certificate_kind === :identical_channels
        @test identical.diamond_norm_result === nothing
        # Pinned QETLAB returns |0.8 - 0.2| = 0.6 here and labels it a
        # probability; the channel Holevo--Helstrom theorem gives 0.8.
        @test identical.success_probability != 0.6

        deterministic = QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel; priors=[1.0, 0.0]
        )
        @test deterministic.success_probability == 1
        @test deterministic.certificate_kind === :deterministic_prior

        unavailable = QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel
        )
        @test unavailable.status ===
            QETChannelOptimization.ChannelOptimizationBackendUnavailable
        @test unavailable.success_probability === nothing
        @test unavailable.lower_bound == 0.5
        @test unavailable.upper_bound == 1
        @test unavailable.diamond_norm_result isa QETChannelOptimization.ChannelNormResult
        @test unavailable.diamond_norm_result.optimization_result.status ===
            QETChannelOptimization.OptimizationBackendUnavailable

        @test_throws DomainError QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel; priors=[1.1, -0.1]
        )
        @test_throws ArgumentError QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel; priors=[0.4, 0.4]
        )
        @test_throws DimensionMismatch QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel; priors=[1.0]
        )
        @test_throws DomainError QETChannelOptimization.channel_distinguishability(
            identity_channel, transpose_map
        )
        non_trace_preserving = QETChannelOptimization.KrausRepresentation([
            ComplexF64[2 0; 0 2]
        ])
        @test_throws DomainError QETChannelOptimization.channel_distinguishability(
            identity_channel, non_trace_preserving
        )
        one_dimensional = QETChannelOptimization.KrausRepresentation([
            ones(ComplexF64, 1, 1)
        ])
        @test_throws DimensionMismatch QETChannelOptimization.channel_distinguishability(
            identity_channel, one_dimensional
        )
    end

    @testset "maximum output fidelity" begin
        identical = QETChannelOptimization.maximum_output_fidelity(
            identity_channel, identity_channel
        )
        @test identical isa QETChannelOptimization.MaximumOutputFidelityResult{Float64}
        @test identical.status === QETChannelOptimization.ChannelOptimizationAnalyticOptimal
        @test identical.value == 1
        @test identical.lower_bound == identical.upper_bound == 1
        @test identical.certified
        @test identical.certificate_kind === :identical_channels
        @test identical.input_states !== nothing
        @test identical.output_states !== nothing

        orthogonal_replacers = QETChannelOptimization.maximum_output_fidelity(
            reset_zero, reset_one
        )
        @test orthogonal_replacers.value == 0
        @test orthogonal_replacers.lower_bound == orthogonal_replacers.upper_bound == 0
        @test orthogonal_replacers.certified
        @test orthogonal_replacers.certificate_kind === :replacer_channels
        @test orthogonal_replacers.certificate.output_fidelity_singular_values == [0, 0]

        phase_flip = QETChannelOptimization.KrausRepresentation([
            0.8Matrix{ComplexF64}(I, 2, 2), 0.6ComplexF64[1 0; 0 -1]
        ])
        common_output = QETChannelOptimization.maximum_output_fidelity(
            identity_channel, phase_flip
        )
        @test common_output.value == 1
        @test common_output.certificate_kind === :common_basis_output
        @test common_output.certificate.first_basis_index == 1
        @test common_output.output_states[1] == common_output.output_states[2]

        problem = QETChannelOptimization.maximum_output_fidelity_problem(
            identity_channel, rotation_channel
        )
        @test problem isa QETChannelOptimization.MaximumOutputFidelityProblem
        @test problem.program.name === :maximum_output_fidelity_primal
        @test problem.program.sense === :maximize
        @test problem.program.variable_count == 24
        @test length(problem.program.equalities) == 10
        @test length(problem.program.psd_constraints) == 3
        @test length(problem.program.primal_views) == 5
        @test QETChannelOptimization.primal_residual(
            problem.program, problem.program.known_feasible_point; allow_densify=true
        ) == 0

        unavailable = QETChannelOptimization.maximum_output_fidelity(
            identity_channel, rotation_channel
        )
        @test unavailable.status ===
            QETChannelOptimization.ChannelOptimizationBackendUnavailable
        @test unavailable.value === nothing
        @test unavailable.lower_bound == 0
        @test unavailable.upper_bound == 1
        @test unavailable.problem isa QETChannelOptimization.MaximumOutputFidelityProblem
        @test unavailable.optimization_result.status ===
            QETChannelOptimization.OptimizationBackendUnavailable

        identity32 = channel_opt_identity(ComplexF32)
        fidelity32 = QETChannelOptimization.maximum_output_fidelity(identity32, identity32)
        @test fidelity32 isa QETChannelOptimization.MaximumOutputFidelityResult{Float32}
        @test fidelity32.value isa Float32

        @test_throws DomainError QETChannelOptimization.maximum_output_fidelity(
            identity_channel, transpose_map
        )
        non_trace_preserving = QETChannelOptimization.KrausRepresentation([
            ComplexF64[2 0; 0 2]
        ])
        @test_throws DomainError QETChannelOptimization.maximum_output_fidelity(
            identity_channel, non_trace_preserving
        )
        @test QETChannelOptimization.maximum_output_fidelity(
            non_trace_preserving, non_trace_preserving; require_trace_preserving=false
        ).status === QETChannelOptimization.ChannelOptimizationBackendUnavailable
        @test_throws ArgumentError QETChannelOptimization.maximum_output_fidelity_problem(
            identity_channel, rotation_channel; max_variables=23
        )
        @test_throws DimensionMismatch QETChannelOptimization.maximum_output_fidelity(
            identity_channel,
            QETChannelOptimization.KrausRepresentation([ones(ComplexF64, 1, 1)]),
        )
    end

    @testset "nonmutation and display" begin
        first_choi = QETChannelOptimization.choi_matrix(identity_channel)
        second_choi = QETChannelOptimization.choi_matrix(dephasing_channel)
        first_copy = copy(first_choi)
        second_copy = copy(second_choi)
        result = QETChannelOptimization.channel_distinguishability(
            identity_channel, dephasing_channel
        )
        @test first_choi == first_copy
        @test second_choi == second_copy
        @test occursin("ChannelDistinguishabilityResult", sprint(show, result))
        @test occursin("ChannelNormResult", sprint(show, result.diamond_norm_result))
        @test occursin(
            "MaximumOutputFidelityResult",
            sprint(
                show,
                QETChannelOptimization.maximum_output_fidelity(
                    identity_channel, identity_channel
                ),
            ),
        )
    end
end
