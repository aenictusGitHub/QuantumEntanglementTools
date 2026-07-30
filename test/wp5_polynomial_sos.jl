using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const PolynomialSOSQET = QuantumEntanglementTools

if !isdefined(PolynomialSOSQET, :PolynomialSOSResult)
    Base.include(
        PolynomialSOSQET,
        joinpath(@__DIR__, "..", "src", "optimization", "polynomial_sos.jl"),
    )
end

function polynomial_sos_hermitian_coordinates(matrix)
    dimension = size(matrix, 1)
    coordinates = Vector{real(eltype(matrix))}()
    append!(coordinates, real.(diag(matrix)))
    for column in 2:dimension, row in 1:(column - 1)
        push!(coordinates, real(matrix[row, column]))
    end
    for column in 2:dimension, row in 1:(column - 1)
        push!(coordinates, imag(matrix[row, column]))
    end
    return coordinates
end

@testset "WP5 polynomial SOS model and status contract" begin
    polynomial = PolynomialSOSQET.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
    problem = PolynomialSOSQET.polynomial_sos_problem(polynomial)
    @test problem.sense === :maximize
    @test problem.variable_count == 4
    @test length(problem.equalities) == 5
    @test length(problem.psd_constraints) == 1
    @test problem.psd_constraints[1].name === :polynomial_sos_moment
    @test problem.metadata.hierarchy_level == 0
    @test problem.metadata.tensor_copies == 1
    @test problem.metadata.full_tensor_dimension == 2
    @test problem.metadata.moment_dimension == 2
    @test problem.metadata.outer_semantics === :upper_bound

    maximally_mixed = Matrix{ComplexF64}(I, 2, 2) / 2
    mixed_coordinates = polynomial_sos_hermitian_coordinates(maximally_mixed)
    @test PolynomialSOSQET.primal_residual(
        problem, mixed_coordinates; allow_densify=true
    ) <= 1e-14
    @test PolynomialSOSQET.evaluate_affine(problem.objective, mixed_coordinates) ≈ 2

    minimum_problem = PolynomialSOSQET.polynomial_sos_problem(polynomial; sense=:min)
    @test minimum_problem.sense === :minimize
    @test minimum_problem.metadata.outer_semantics === :lower_bound

    level_one = PolynomialSOSQET.polynomial_sos_problem(polynomial; level=1)
    @test level_one.variable_count == 9
    @test level_one.metadata.tensor_copies == 2
    @test level_one.metadata.full_tensor_dimension == 4
    @test level_one.metadata.moment_dimension == 3
    @test length(level_one.equalities) == 17

    point = [3.0, 4.0] / 5
    tensor_point = kron(point, point)
    basis = PolynomialSOSQET.symmetric_subspace_basis(2, 2)
    symmetric_coordinates = transpose(basis) * tensor_point
    moment = symmetric_coordinates * transpose(symmetric_coordinates)
    level_one_coordinates = polynomial_sos_hermitian_coordinates(moment)
    @test PolynomialSOSQET.primal_residual(
        level_one, level_one_coordinates; allow_densify=true
    ) <= 2e-14
    @test PolynomialSOSQET.evaluate_affine(level_one.objective, level_one_coordinates) ≈
        PolynomialSOSQET.evaluate_polynomial(polynomial, point) atol = 2e-14

    rng = MersenneTwister(71)
    unavailable = PolynomialSOSQET.polynomial_sos_bounds(rng, polynomial; inner_samples=8)
    @test unavailable isa PolynomialSOSQET.PolynomialSOSResult
    @test unavailable.status === PolynomialSOSQET.OptimizationBackendUnavailable
    @test unavailable.outer_bound === nothing
    @test unavailable.inner_bound !== nothing
    @test unavailable.inner_kind === :sampled_feasible
    @test unavailable.samples_requested == 8
    @test unavailable.samples_evaluated == 8
    @test unavailable.best_point !== nothing
    @test norm(unavailable.best_point) ≈ 1 atol = 2e-15
    @test PolynomialSOSQET.evaluate_polynomial(polynomial, unavailable.best_point) ≈
        unavailable.inner_bound atol = 2e-15
    @test unavailable.optimization_result.status ===
        PolynomialSOSQET.OptimizationBackendUnavailable
    @test !unavailable.certified_outer

    same_one = PolynomialSOSQET.polynomial_sos_bounds(
        MersenneTwister(91), polynomial; inner_samples=5
    )
    same_two = PolynomialSOSQET.polynomial_sos_bounds(
        MersenneTwister(91), polynomial; inner_samples=5
    )
    @test same_one.inner_bound == same_two.inner_bound
    @test same_one.best_point == same_two.best_point

    Random.seed!(20260730)
    expected_global_draw = rand()
    Random.seed!(20260730)
    PolynomialSOSQET.polynomial_sos_bounds(MersenneTwister(3), polynomial; inner_samples=2)
    @test rand() == expected_global_draw

    constant = PolynomialSOSQET.HomogeneousPolynomial([2.5], 3, 0)
    constant_result = PolynomialSOSQET.polynomial_sos_bounds(
        MersenneTwister(1), constant; level=4, inner_samples=10
    )
    @test constant_result.status === PolynomialSOSQET.OptimizationOptimal
    @test constant_result.outer_bound == 2.5
    @test constant_result.inner_bound == 2.5
    @test constant_result.outer_kind === :exact_constant
    @test constant_result.inner_kind === :exact_constant
    @test constant_result.certified_outer
    @test constant_result.problem === nothing
    @test constant_result.optimization_result === nothing

    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        PolynomialSOSQET.HomogeneousPolynomial([1.0, 2.0], 2, 1)
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        PolynomialSOSQET.HomogeneousPolynomial(ComplexF64[1, 0, 1], 2, 2)
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(polynomial; level=-1)
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(constant)
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        polynomial; max_full_dimension=1
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        polynomial; limits=PolynomialSOSQET.OptimizationLimits(max_variables=3)
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        polynomial; limits=PolynomialSOSQET.OptimizationLimits(max_equalities=4)
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_problem(
        polynomial; max_work=15
    )
    @test_throws ArgumentError PolynomialSOSQET.polynomial_sos_bounds(
        MersenneTwister(1), polynomial; inner_samples=3, max_samples=2
    )
end
