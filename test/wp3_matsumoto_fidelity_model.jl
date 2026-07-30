using LinearAlgebra
using QuantumEntanglementTools
using SparseArrays
using Test

const MatsumotoModelQET = QuantumEntanglementTools

if !isdefined(MatsumotoModelQET, :MatsumotoFidelityModel)
    Base.include(
        MatsumotoModelQET,
        joinpath(@__DIR__, "..", "src", "optimization", "matsumoto_fidelity_model.jl"),
    )
end

@testset "WP3 Matsumoto-fidelity SDP model" begin
    rho = Diagonal([0.7, 0.3])
    sigma = Diagonal([0.4, 0.6])
    expected = sum(sqrt.(diag(rho)) .* sqrt.(diag(sigma)))
    model = MatsumotoModelQET.matsumoto_fidelity_model(rho, sigma)

    @test model isa MatsumotoModelQET.MatsumotoFidelityModel{Float64}
    @test model.dimension == 2
    @test model.input_variable_count == 0
    @test model.variable_count == 4
    @test model.coupling_variables == 1:4
    @test model.coupling.name === :matsumoto_fidelity_coupling
    @test model.psd_constraint.name === :matsumoto_fidelity_block
    @test model.psd_constraint.dimension == 4
    traced_coupling = MatsumotoModelQET.trace_affine(model.coupling)
    @test model.objective.constant == traced_coupling.constant
    @test model.objective.coefficients == traced_coupling.coefficients

    zero_coordinates = zeros(4)
    zero_block = Matrix(
        MatsumotoModelQET.evaluate_affine(model.psd_constraint, zero_coordinates)
    )
    @test zero_block == [Matrix(rho) zeros(2, 2); zeros(2, 2) Matrix(sigma)]
    @test ishermitian(zero_block)
    @test minimum(eigvals(Hermitian(zero_block))) >= 0

    optimal_coordinates = [
        sqrt(rho[1, 1] * sigma[1, 1]), sqrt(rho[2, 2] * sigma[2, 2]), 0.0, 0.0
    ]
    @test MatsumotoModelQET.evaluate_affine(model.objective, optimal_coordinates) ≈ expected
    optimal_block = Matrix(
        MatsumotoModelQET.evaluate_affine(model.psd_constraint, optimal_coordinates)
    )
    @test minimum(eigvals(Hermitian(optimal_block))) >= -1e-14

    problem = MatsumotoModelQET.matsumoto_fidelity_problem(rho, sigma)
    @test problem.sense === :maximize
    @test problem.variable_count == 4
    @test length(problem.psd_constraints) == 1
    @test problem.known_feasible_point == zeros(4)
    @test MatsumotoModelQET.primal_residual(
        problem, problem.known_feasible_point; allow_densify=true
    ) == 0
    unavailable = MatsumotoModelQET.solve_optimization(problem)
    @test unavailable.status === MatsumotoModelQET.OptimizationBackendUnavailable
    @test unavailable.primal === nothing
    @test unavailable.dual === nothing
    @test !unavailable.certified

    rho32 = Diagonal(Float32[0.75, 0.25])
    sigma32 = Diagonal(Float32[0.5, 0.5])
    model32 = MatsumotoModelQET.matsumoto_fidelity_model(rho32, sigma32)
    @test model32 isa MatsumotoModelQET.MatsumotoFidelityModel{Float32}
    @test eltype(model32.psd_constraint.constant) === ComplexF32
    @test eltype(model32.objective.coefficients) === Float32

    sparse_rho = sparse(Matrix(rho))
    sparse_sigma = sparse(Matrix(sigma))
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(
        sparse_rho, sparse_sigma
    )
    sparse_model = MatsumotoModelQET.matsumoto_fidelity_model(
        sparse_rho, sparse_sigma; allow_densify=true
    )
    @test issparse(sparse_model.psd_constraint.constant)

    @test_throws DimensionMismatch MatsumotoModelQET.matsumoto_fidelity_model(
        rho, Matrix{Float64}(I, 3, 3) / 3
    )
    @test_throws DimensionMismatch MatsumotoModelQET.matsumoto_fidelity_model(
        ones(2, 3), ones(2, 3)
    )
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(
        ComplexF64[0.5 0.1im; 0.1im 0.5], sigma
    )
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(0.9rho, sigma)
    @test_throws DomainError MatsumotoModelQET.matsumoto_fidelity_model(
        [1.1 0.0; 0.0 -0.1], sigma
    )
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(
        rho, sigma; limits=MatsumotoModelQET.OptimizationLimits(max_variables=3)
    )
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(
        rho, sigma; limits=MatsumotoModelQET.OptimizationLimits(max_psd_dimension=3)
    )
    @test_throws ArgumentError MatsumotoModelQET.matsumoto_fidelity_model(
        rho, sigma; limits=MatsumotoModelQET.OptimizationLimits(max_model_entries=9)
    )

    affine_rho = MatsumotoModelQET.hermitian_variable(:rho, 2; coefficient_type=Float64)
    affine_sigma = MatsumotoModelQET.HermitianAffineMatrix(
        :sigma, Matrix(sigma), Int[], Matrix{Float64}[], affine_rho.variable_count
    )
    affine_model = MatsumotoModelQET.matsumoto_fidelity_model(
        affine_rho, affine_sigma; name=:affine_matsumoto
    )
    @test affine_model.input_variable_count == 4
    @test affine_model.variable_count == 8
    @test affine_model.coupling_variables == 5:8
    @test all(
        term.variable <= 4 for
        term in affine_model.psd_constraint.terms[1:length(affine_rho.terms)]
    )

    fixed_rho = Matrix(rho)
    rho_equalities = MatsumotoModelQET.hermitian_equalities(
        affine_rho; target=fixed_rho, name_prefix=:fixed_rho
    )
    affine_problem = MatsumotoModelQET.matsumoto_fidelity_problem(
        affine_rho,
        affine_sigma;
        equalities=rho_equalities,
        psd_constraints=[affine_rho],
        primal_views=[affine_rho],
        initial_point=[0.7, 0.3, 0.0, 0.0],
        known_feasible_point=[0.7, 0.3, 0.0, 0.0],
        name=:affine_matsumoto,
    )
    @test affine_problem.variable_count == 8
    @test length(affine_problem.equalities) == 4
    @test length(affine_problem.psd_constraints) == 2
    @test affine_problem.initial_point == [0.7, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test MatsumotoModelQET.primal_residual(
        affine_problem, affine_problem.known_feasible_point; allow_densify=true
    ) == 0

    mismatched_count = MatsumotoModelQET.HermitianAffineMatrix(
        :bad_sigma, Matrix(sigma), Int[], Matrix{Float64}[], 5
    )
    @test_throws DimensionMismatch MatsumotoModelQET.matsumoto_fidelity_model(
        affine_rho, mismatched_count
    )
end
