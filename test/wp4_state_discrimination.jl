using Test
using LinearAlgebra
using SparseArrays

using QuantumEntanglementTools: QuantumEntanglementTools
const QETDiscrimination = QuantumEntanglementTools

if !isdefined(QETDiscrimination, :StateDiscriminationResult)
    Base.include(
        QETDiscrimination,
        joinpath(@__DIR__, "..", "src", "optimization", "state_discrimination.jl"),
    )
end

@testset "WP4 minimum-error state discrimination" begin
    ket0 = ComplexF64[1, 0]
    ket1 = ComplexF64[0, 1]
    ket_plus = ComplexF64[1, 1] / sqrt(2)
    rho0 = ket0 * ket0'
    rho1 = ket1 * ket1'
    rho_plus = ket_plus * ket_plus'

    @testset "trivial and Helstrom branches" begin
        single = QETDiscrimination.state_distinguishability((rho0,))
        @test single.status === QETDiscrimination.StateDiscriminationTrivialOptimal
        @test single.success_probability == 1
        @test single.lower_bound == single.upper_bound == 1
        @test single.certified
        @test single.certificate_kind === :deterministic_prior
        @test length(single.measurement) == 1
        @test single.measurement[1] == Matrix{ComplexF64}(I, 2, 2)
        @test single.residuals.valid

        deterministic = QETDiscrimination.state_distinguishability(
            (rho0, rho_plus, rho1); priors=[0.0, 1.0, 0.0]
        )
        @test deterministic.status === QETDiscrimination.StateDiscriminationTrivialOptimal
        @test deterministic.success_probability ≈ 1
        @test iszero(deterministic.measurement[1])
        @test deterministic.measurement[2] == Matrix{ComplexF64}(I, 2, 2)
        @test iszero(deterministic.measurement[3])

        orthogonal = QETDiscrimination.state_distinguishability(hcat(ket0, ket1))
        @test orthogonal.status === QETDiscrimination.StateDiscriminationHelstromOptimal
        @test orthogonal.success_probability == 1
        @test orthogonal.certified
        @test orthogonal.certificate_kind === :helstrom_theorem
        @test orthogonal.residuals.valid
        @test sum(orthogonal.measurement) ≈ Matrix{ComplexF64}(I, 2, 2)

        nonorthogonal = QETDiscrimination.state_distinguishability(hcat(ket0, ket_plus))
        expected = (1 + inv(sqrt(2))) / 2
        @test nonorthogonal.status === QETDiscrimination.StateDiscriminationHelstromOptimal
        @test nonorthogonal.success_probability ≈ expected
        @test nonorthogonal.lower_bound ≈ expected
        @test nonorthogonal.upper_bound ≈ expected
        @test nonorthogonal.certified
        @test nonorthogonal.residuals.valid
        @test nonorthogonal.residuals.helstrom_objective_residual ≤ 1e-14
        @test all(
            effect -> ishermitian(effect) && minimum(eigvals(Hermitian(effect))) ≥ -1e-14,
            nonorthogonal.measurement,
        )

        unequal = QETDiscrimination.state_distinguishability(
            (rho0, rho_plus); priors=[0.8, 0.2]
        )
        expected_unequal = (1 + sqrt((0.8 - 0.2)^2 + 4 * 0.8 * 0.2 * (1 - 0.5))) / 2
        @test unequal.success_probability ≈ expected_unequal
        @test unequal.residuals.valid

        identical = QETDiscrimination.state_distinguishability(
            (rho0, rho0); priors=[0.8, 0.2]
        )
        @test identical.success_probability ≈ 0.8
        @test identical.residuals.valid

        mixed_left = ComplexF64[0.75 0; 0 0.25]
        mixed_right = ComplexF64[0.25 0; 0 0.75]
        mixed = QETDiscrimination.state_distinguishability((mixed_left, mixed_right))
        @test mixed.success_probability ≈ 0.75
        @test mixed.residuals.valid
    end

    @testset "exact multi-state orthogonal shortcut" begin
        basis = Matrix{ComplexF64}(I, 3, 3)
        result = QETDiscrimination.state_distinguishability(basis)
        @test result.status === QETDiscrimination.StateDiscriminationOrthogonalOptimal
        @test result.success_probability == 1
        @test result.certified
        @test result.certificate_kind === :orthogonal_supports
        @test result.residuals.valid
        @test sum(result.measurement) ≈ Matrix{ComplexF64}(I, 3, 3)
        @test all(
            index ->
                real(dot(basis[:, index] * basis[:, index]', result.measurement[index])) ≈
                1,
            1:3,
        )

        block_states = (
            ComplexF64[0.5 0 0; 0 0.5 0; 0 0 0], ComplexF64[0 0 0; 0 0 0; 0 0 1]
        )
        block_result = QETDiscrimination.state_distinguishability(block_states)
        @test block_result.success_probability == 1
        @test block_result.residuals.valid
    end

    @testset "solver-neutral model and unavailable backend" begin
        omega = cis(2pi / 3)
        trine = hcat(
            ComplexF64[1, 1] / sqrt(2),
            ComplexF64[1, omega] / sqrt(2),
            ComplexF64[1, omega ^ 2] / sqrt(2),
        )
        problem = QETDiscrimination.state_discrimination_problem(trine)
        @test problem.dimension == 2
        @test problem.state_count == 3
        @test problem.input_kind === :pure_columns
        @test problem.program.name === :minimum_error_state_discrimination
        @test problem.program.sense === :maximize
        @test problem.program.variable_count == 12
        @test length(problem.program.equalities) == 4
        @test length(problem.program.psd_constraints) == 3
        @test length(problem.program.primal_views) == 3
        @test problem.program.known_feasible_point !== nothing
        @test QETDiscrimination.primal_residual(
            problem.program, problem.program.known_feasible_point; allow_densify=true
        ) ≤ 1e-14
        @test all(
            index ->
                QETDiscrimination.evaluate_affine(
                    problem.measurement_views[index], problem.program.known_feasible_point
                ) ≈ Matrix{ComplexF64}(I, 2, 2) / 3,
            1:3,
        )

        unavailable = QETDiscrimination.state_distinguishability(trine)
        @test unavailable.status === QETDiscrimination.StateDiscriminationBackendUnavailable
        @test unavailable.success_probability === nothing
        @test unavailable.lower_bound === nothing
        @test unavailable.upper_bound === nothing
        @test unavailable.measurement === nothing
        @test unavailable.problem isa QETDiscrimination.StateDiscriminationProblem
        @test unavailable.optimization_result.status ===
            QETDiscrimination.OptimizationBackendUnavailable
        @test !unavailable.certified
    end

    @testset "type behavior and nonmutation" begin
        pure32 = Float32[1 0 inv(sqrt(Float32(2))); 0 1 inv(sqrt(Float32(2)))]
        # The first two columns are orthogonal but the third is not, so the
        # no-backend route must retain a Float32 SDP problem.
        result32 = QETDiscrimination.state_distinguishability(pure32)
        @test result32 isa QETDiscrimination.StateDiscriminationResult{Float32}
        @test result32.problem.program.objective.constant isa Float32
        @test eltype(result32.priors) === Float32

        rho32 = Float32[1 0; 0 0]
        helstrom32 = QETDiscrimination.state_distinguishability((rho32, Float32[0 0; 0 1]))
        @test helstrom32.success_probability isa Float32
        @test all(effect -> eltype(effect) === ComplexF32, helstrom32.measurement)

        states = [copy(rho0), copy(rho_plus)]
        original_states = deepcopy(states)
        priors = [0.4, 0.6]
        original_priors = copy(priors)
        _ = QETDiscrimination.state_distinguishability(states; priors=priors)
        @test states == original_states
        @test priors == original_priors
    end

    @testset "validation and resource limits" begin
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(())
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            zeros(ComplexF64, 0, 2)
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            hcat(ket0, 2ket1)
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability((rho0, 2rho1))
        @test_throws DimensionMismatch QETDiscrimination.state_distinguishability((
            rho0, Matrix{ComplexF64}(I, 3, 3) / 3
        ))
        @test_throws DimensionMismatch QETDiscrimination.state_distinguishability((
            ones(ComplexF64, 2, 3),
        ))
        @test_throws ArgumentError QETDiscrimination.state_distinguishability((
            ComplexF64[1 im * eps(); 0 0],
        ))
        @test_throws DomainError QETDiscrimination.state_distinguishability((
            ComplexF64[1.1 0; 0 -0.1],
        ))
        @test_throws DimensionMismatch QETDiscrimination.state_distinguishability(
            (rho0, rho1); priors=[1.0]
        )
        @test_throws DomainError QETDiscrimination.state_distinguishability(
            (rho0, rho1); priors=[1.1, -0.1]
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); priors=[0.4, 0.4]
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); priors=[NaN, NaN]
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); max_dimension=1
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); max_states=1
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); max_variables=7
        )
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            (rho0, rho1); max_states=true
        )

        sparse_rho0 = sparse(rho0)
        sparse_rho1 = sparse(rho1)
        @test_throws ArgumentError QETDiscrimination.state_distinguishability((
            sparse_rho0, sparse_rho1
        ))
        sparse_result = QETDiscrimination.state_distinguishability(
            (sparse_rho0, sparse_rho1); allow_densify=true
        )
        @test sparse_result.success_probability == 1

        sparse_columns = sparse(hcat(ket0, ket1))
        @test_throws ArgumentError QETDiscrimination.state_distinguishability(
            sparse_columns
        )
        @test QETDiscrimination.state_distinguishability(
            sparse_columns; allow_densify=true
        ).success_probability == 1

        @test_throws ArgumentError QETDiscrimination.state_discrimination_problem(
            hcat(ket0, ket_plus, ket1);
            limits=QETDiscrimination.OptimizationLimits(max_variables=8),
        )
        @test_throws ArgumentError QETDiscrimination.state_discrimination_problem(
            hcat(ket0, ket_plus, ket1);
            limits=QETDiscrimination.OptimizationLimits(max_psd_blocks=2),
        )
        @test_throws ArgumentError QETDiscrimination.state_discrimination_problem(
            hcat(ket0, ket_plus, ket1);
            limits=QETDiscrimination.OptimizationLimits(max_psd_dimension=1),
        )
    end
end
