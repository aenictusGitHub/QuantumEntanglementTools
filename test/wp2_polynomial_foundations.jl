using LinearAlgebra
using QuantumEntanglementTools
using Random
using SparseArrays
using Test

const PF = QuantumEntanglementTools
const PFCompat = QuantumEntanglementTools.MATLABCompat

struct _PolynomialZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_PolynomialZeroBasedVector) = size(vector.storage)
Base.axes(vector::_PolynomialZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_PolynomialZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_PolynomialZeroBasedVector, index::Int) = vector.storage[index + 1]

struct _PolynomialZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_PolynomialZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_PolynomialZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_PolynomialZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_PolynomialZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function _normalized_monomial_basis(exponents, point)
    degree = sum(@view exponents[1, :])
    values = Vector{Float64}(undef, size(exponents, 1))
    degree_factorial = factorial(big(degree))
    for row in axes(exponents, 1)
        denominator = prod(
            factorial(big(exponents[row, column])) for column in axes(exponents, 2)
        )
        normalization = sqrt(Float64(degree_factorial // denominator))
        values[row] =
            normalization *
            prod(point[column]^exponents[row, column] for column in axes(exponents, 2))
    end
    return values
end

@testset "WP2 solver-free polynomial foundations" begin
    @testset "owned homogeneous polynomials and QETLAB order" begin
        expected = [
            2 0 0
            1 1 0
            1 0 1
            0 2 0
            0 1 1
            0 0 2
        ]
        @test PF.monomial_exponents(3, 2) == expected

        coefficients = Rational{Int}[1, 2, 3, 4, 5, 6]
        polynomial = PF.HomogeneousPolynomial(coefficients, 3, 2)
        coefficients[1] = 99
        @test polynomial.coefficients[1] == 1
        @test eltype(polynomial) == Rational{Int}
        @test length(polynomial) == 6
        @test PF.monomial_exponents(polynomial) == expected
        @test occursin("dense Rational", sprint(show, polynomial))

        sparse_coefficients = sparsevec([1, 6], [2, 3], 6)
        sparse_polynomial = PF.HomogeneousPolynomial(sparse_coefficients, 3, 2)
        sparse_coefficients[1] = 8
        @test sparse_polynomial.coefficients isa SparseVector
        @test sparse_polynomial.coefficients[1] == 2
        @test nnz(sparse_polynomial.coefficients) == 2
        @test occursin("sparse Int", sprint(show, sparse_polynomial))
        @test copy(sparse_polynomial).coefficients !== sparse_polynomial.coefficients

        exact_point = Rational{Int}[1 // 2, 1 // 3, -1 // 4]
        exact_value = PF.evaluate_polynomial(polynomial, exact_point)
        direct_value = sum(
            polynomial.coefficients[row] *
            prod(exact_point[column]^expected[row, column] for column in axes(expected, 2))
            for row in axes(expected, 1)
        )
        @test exact_value == direct_value
        @test exact_value isa Rational{Int}

        @test_throws DimensionMismatch PF.HomogeneousPolynomial([1, 2], 2, 2)
        @test_throws ArgumentError PF.HomogeneousPolynomial(
            [1, 2, 3], 2, 2; coefficient_order=:graded_reverse_lexicographic
        )
        @test_throws ArgumentError PF.HomogeneousPolynomial(Bool[true, false, true], 2, 2)
        @test_throws ArgumentError PF.HomogeneousPolynomial([1.0, NaN, 2.0], 2, 2)
        @test_throws ArgumentError PF.HomogeneousPolynomial([1, 2, 3], 2, 2; max_terms=2)
        @test_throws ArgumentError PF.monomial_exponents(0, 2)
        @test_throws ArgumentError PF.monomial_exponents(2, -1)
        @test_throws ArgumentError PF.monomial_exponents(50, 0; max_exponent_entries=49)
        @test_throws ArgumentError PF.HomogeneousPolynomial(
            _PolynomialZeroBasedVector([1.0, 2.0, 3.0]), 2, 2
        )
        @test_throws DimensionMismatch PF.evaluate_polynomial(polynomial, [1, 2])
        @test_throws ArgumentError PF.evaluate_polynomial(polynomial, [1, 2, Inf])
        @test_throws ArgumentError PF.evaluate_polynomial(polynomial, [1, 2, 3]; max_work=1)
        @test_throws ArgumentError PF.evaluate_polynomial(
            polynomial, [1, 2, 3]; max_exponent_entries=17
        )
        @test_throws ArgumentError PF.evaluate_polynomial(
            polynomial, _PolynomialZeroBasedVector([1, 2, 3])
        )

        mutable_polynomial = PF.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
        mutable_polynomial.coefficients[2] = NaN
        @test_throws ArgumentError PF.evaluate_polynomial(mutable_polynomial, [1.0, 2.0])
    end

    @testset "copositive quartic construction" begin
        matrix = [1 2; 2 3]
        matrix_copy = copy(matrix)
        polynomial = PF.copositive_polynomial(matrix)
        @test matrix == matrix_copy
        @test polynomial.variables == 2
        @test polynomial.degree == 4
        @test polynomial.coefficients isa SparseVector{Int}
        @test collect(polynomial.coefficients) == [1, 0, 4, 0, 3]
        @test nnz(polynomial.coefficients) == 3

        point = Rational{Int}[2 // 3, -3 // 5]
        squared = point .^ 2
        @test PF.evaluate_polynomial(polynomial, point) == dot(squared, matrix * squared)

        sparse_matrix = sparse(matrix)
        sparse_result = PF.copositive_polynomial(sparse_matrix)
        @test sparse_result.coefficients == polynomial.coefficients
        compatibility = PFCompat.CopositivePolynomial(matrix)
        @test compatibility isa SparseVector{Int}
        @test compatibility == polynomial.coefficients
        @test PFCompat.CopositivePolynomial(matrix; dense_output=true) ==
            collect(polynomial.coefficients)

        result32 = PF.copositive_polynomial(Float32.(matrix))
        @test eltype(result32) == Float32
        @test collect(result32.coefficients) == Float32[1, 0, 4, 0, 3]

        @test_throws DimensionMismatch PF.copositive_polynomial(ones(2, 3))
        @test_throws ArgumentError PF.copositive_polynomial(zeros(0, 0))
        @test_throws ArgumentError PF.copositive_polynomial([1 2; 3 4])
        @test_throws ArgumentError PFCompat.CopositivePolynomial([1 2; 3 4])
        @test_throws ArgumentError PF.copositive_polynomial(ComplexF64[1 0; 0 1])
        @test_throws ArgumentError PF.copositive_polynomial([1.0 NaN; NaN 1.0])
        @test_throws ArgumentError PF.copositive_polynomial(Bool[true false; false true])
        @test_throws ArgumentError PF.copositive_polynomial(
            _PolynomialZeroBasedMatrix(matrix)
        )
        @test_throws ArgumentError PF.copositive_polynomial(matrix; max_terms=4)
        overflow_matrix = [0 typemax(Int); typemax(Int) 0]
        @test_throws OverflowError PF.copositive_polynomial(overflow_matrix)
    end

    @testset "compact symmetric polynomial matrices" begin
        quadratic = PF.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
        expected_level_zero = [1.0 1.0; 1.0 3.0]
        level_zero = PF.polynomial_as_matrix(quadratic)
        @test level_zero isa SparseMatrixCSC{Float64,Int}
        @test Matrix(level_zero) == expected_level_zero
        @test issymmetric(level_zero)
        @test PF.polynomial_as_matrix([1.0, 2.0, 3.0], 2, 1; sparse_output=false) ==
            expected_level_zero
        @test PFCompat.PolynomialAsMatrix([1.0 2.0 3.0], 2, 1) ==
            sparse(expected_level_zero)
        @test PFCompat.PolynomialAsMatrix(
            reshape([1.0, 2.0, 3.0], :, 1), 2, 1, 0; sparse_output=false
        ) == expected_level_zero

        # Source-free values captured from the pinned QETLAB matrix formula.
        expected_level_two = [
            1.0 1 / sqrt(3) 0.0 0.0
            1 / sqrt(3) 5 / 3 2 / 3 0.0
            0.0 2 / 3 7 / 3 1 / sqrt(3)
            0.0 0.0 1 / sqrt(3) 3.0
        ]
        @test Matrix(PF.polynomial_as_matrix(quadratic; level=2)) ≈ expected_level_two atol =
            2e-15 rtol = 2e-15

        quartic = PF.HomogeneousPolynomial([1.0, 0.0, 2.0, 0.0, 3.0], 2, 4)
        expected_quartic_level_one = [
            1.0 0.0 1 / (3sqrt(3)) 0.0
            0.0 7 / 9 0.0 1 / (3sqrt(3))
            1 / (3sqrt(3)) 0.0 13 / 9 0.0
            0.0 1 / (3sqrt(3)) 0.0 3.0
        ]
        @test Matrix(PF.polynomial_as_matrix(quartic; level=1)) ≈ expected_quartic_level_one atol =
            3e-15 rtol = 3e-15

        for variables in 1:4, half_degree in 0:3, level in 0:2
            coefficient_count = binomial(variables + 2 * half_degree - 1, 2 * half_degree)
            rng = MersenneTwister(10_000 * variables + 100 * half_degree + level)
            coefficients = randn(rng, coefficient_count)
            polynomial = PF.HomogeneousPolynomial(coefficients, variables, 2 * half_degree)
            matrix = Matrix(PF.polynomial_as_matrix(polynomial; level=level))
            exponents = PF.monomial_exponents(variables, half_degree + level)
            point = randn(rng, variables)
            basis = _normalized_monomial_basis(exponents, point)
            represented = dot(basis, matrix * basis)
            direct = PF.evaluate_polynomial(polynomial, point) * sum(abs2, point)^level
            @test represented ≈ direct atol = 2e-10 rtol = 2e-10
        end

        matrix32 = PF.polynomial_as_matrix(PF.HomogeneousPolynomial(Float32[1, 2, 3], 2, 2))
        @test eltype(matrix32) == Float32
        rational_matrix = PF.polynomial_as_matrix(
            PF.HomogeneousPolynomial(Rational{Int}[1, 2, 3], 2, 2)
        )
        @test eltype(rational_matrix) == Float64
        complex_matrix = PF.polynomial_as_matrix(
            PF.HomogeneousPolynomial(ComplexF32[1, 2im, 3], 2, 2)
        )
        @test eltype(complex_matrix) == ComplexF32
        @test transpose(complex_matrix) == complex_matrix
        @test !ishermitian(complex_matrix)

        @test_throws ArgumentError PF.polynomial_as_matrix(
            PF.HomogeneousPolynomial([1.0, 2.0], 2, 1)
        )
        @test_throws ArgumentError PF.polynomial_as_matrix(quadratic; level=-1)
        @test_throws ArgumentError PF.polynomial_as_matrix(quadratic; max_degree=1)
        @test_throws ArgumentError PF.polynomial_as_matrix(
            quadratic; level=2, max_dimension=3
        )
        @test_throws ArgumentError PF.polynomial_as_matrix(quadratic; max_work=3)
        @test_throws ArgumentError PF.polynomial_as_matrix(
            quadratic; level=2, max_exponent_entries=7
        )
        @test_throws ArgumentError PF.polynomial_as_matrix(
            quadratic; sparse_output=false, max_dense_entries=3
        )
        @test_throws ArgumentError PF.polynomial_as_matrix(quadratic; max_nonzeros=3)
        @test_throws DimensionMismatch PFCompat.PolynomialAsMatrix(ones(2, 2), 2, 1)

        corrupted = PF.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
        corrupted.coefficients[1] = Inf
        @test_throws ArgumentError PF.polynomial_as_matrix(corrupted)
    end

    @testset "deterministic hierarchy and sampled bounds" begin
        quadratic = PF.HomogeneousPolynomial([1.0, 2.0, 3.0], 2, 2)
        exact_minimum = 2 - sqrt(2)
        exact_maximum = 2 + sqrt(2)

        maximum_result = PF.polynomial_bounds(
            MersenneTwister(11),
            quadratic;
            level=1,
            sense=:max,
            inner_samples=32,
            allow_densify=true,
        )
        @test maximum_result isa PF.PolynomialOptimizationResult
        @test maximum_result.sense === :max
        @test maximum_result.outer_kind === :hierarchy_outer_bound
        @test maximum_result.outer_status === :validated_numerical
        @test maximum_result.inner_kind === :sampled_feasible
        @test maximum_result.outer_bound >= exact_maximum
        @test maximum_result.outer_bound ≈ exact_maximum atol = 5e-13
        @test maximum_result.inner_bound <= exact_maximum + 2e-15
        @test maximum_result.samples_requested == 32
        @test maximum_result.samples_evaluated == 32
        @test norm(maximum_result.best_point) ≈ 1 atol = 3e-15
        @test PF.evaluate_polynomial(quadratic, maximum_result.best_point) ≈
            maximum_result.inner_bound atol = 2e-15
        @test maximum_result.outer_uncertainty >= 0
        @test maximum_result.sphere_residual <= 3e-15

        minimum_result = PF.polynomial_bounds(
            MersenneTwister(12),
            quadratic;
            level=2,
            sense=:min,
            inner_samples=32,
            allow_densify=true,
        )
        @test minimum_result.outer_bound <= exact_minimum
        @test minimum_result.outer_bound ≈ exact_minimum atol = 5e-13
        @test minimum_result.inner_bound >= exact_minimum - 2e-15

        same_seed_one = PF.polynomial_bounds(
            MersenneTwister(99), quadratic; level=0, inner_samples=12, allow_densify=true
        )
        same_seed_two = PF.polynomial_bounds(
            MersenneTwister(99), quadratic; level=0, inner_samples=12, allow_densify=true
        )
        @test same_seed_one.inner_bound == same_seed_two.inner_bound
        @test same_seed_one.best_point == same_seed_two.best_point

        Random.seed!(20260730)
        expected_global_draw = rand()
        Random.seed!(20260730)
        PF.polynomial_bounds(
            MersenneTwister(7), quadratic; inner_samples=3, allow_densify=true
        )
        @test rand() == expected_global_draw

        no_sample_rng = MersenneTwister(123)
        untouched_rng = copy(no_sample_rng)
        no_sample = PF.polynomial_bounds(
            no_sample_rng, quadratic; inner_samples=0, allow_densify=true
        )
        @test no_sample.inner_kind === :not_requested
        @test no_sample.inner_bound === nothing
        @test no_sample.best_point === nothing
        @test rand(no_sample_rng) == rand(untouched_rng)

        # Regression for the pinned maximization recursion's TARGET-sign bug.
        negative_sphere = PF.HomogeneousPolynomial([-1.0, 0.0, -1.0], 2, 2)
        corrected_skip = PF.polynomial_bounds(
            MersenneTwister(31),
            negative_sphere;
            sense=:max,
            target=3.0,
            inner_samples=5,
            allow_densify=true,
        )
        @test corrected_skip.target_status === :outer_proves_at_most
        @test corrected_skip.inner_kind === :skipped_by_outer_target
        @test corrected_skip.samples_evaluated == 0
        @test !(-corrected_skip.outer_bound >= corrected_skip.target)
        @test -corrected_skip.outer_bound >= -corrected_skip.target

        positive_sphere = PF.HomogeneousPolynomial([1.0, 0.0, 1.0], 2, 2)
        corrected_no_skip = PF.polynomial_bounds(
            MersenneTwister(32),
            positive_sphere;
            sense=:max,
            target=-3.0,
            inner_samples=5,
            allow_densify=true,
        )
        @test corrected_no_skip.target_status === :not_proven
        @test corrected_no_skip.samples_evaluated == 5
        @test -corrected_no_skip.outer_bound >= corrected_no_skip.target
        @test !(-corrected_no_skip.outer_bound >= -corrected_no_skip.target)

        minimum_target = PF.polynomial_bounds(
            MersenneTwister(33),
            positive_sphere;
            sense=:min,
            target=0.0,
            inner_samples=5,
            allow_densify=true,
        )
        @test minimum_target.target_status === :outer_proves_at_least
        @test minimum_target.samples_evaluated == 0

        result32 = PF.polynomial_bounds(
            MersenneTwister(34),
            PF.HomogeneousPolynomial(Float32[1, 0, 1], 2, 2);
            inner_samples=2,
            allow_densify=true,
        )
        @test result32.outer_bound isa Float32
        @test result32.inner_bound isa Float32
        @test eltype(result32.best_point) == Float32

        compatibility = PFCompat.PolynomialOptimize(
            MersenneTwister(35),
            [1.0 2.0 3.0],
            2,
            1,
            1,
            "MAX",
            "NONE";
            inner_samples=4,
            allow_densify=true,
        )
        @test compatibility isa PF.PolynomialOptimizationResult
        @test compatibility.sense === :max
        @test compatibility.samples_evaluated == 4
        @test compatibility.outer_bound ≈ exact_maximum atol = 5e-13

        compatibility_target = PFCompat.PolynomialOptimize(
            MersenneTwister(36),
            reshape([-1.0, 0.0, -1.0], :, 1),
            2,
            1,
            0,
            :max,
            3.0;
            inner_samples=4,
            allow_densify=true,
        )
        @test compatibility_target.target_status === :outer_proves_at_most
        @test compatibility_target.samples_evaluated == 0

        @test_throws ArgumentError PF.polynomial_bounds(MersenneTwister(1), quadratic)
        @test_throws ArgumentError PFCompat.PolynomialOptimize(
            MersenneTwister(1), [1.0, 0.0, 1.0], 2, 1, 0
        )
        @test_throws ArgumentError PFCompat.PolynomialOptimize(
            MersenneTwister(1), [1.0, 0.0, 1.0], 2, 1, 0, "largest"; allow_densify=true
        )
        @test_throws ArgumentError PFCompat.PolynomialOptimize(
            MersenneTwister(1),
            [1.0, 0.0, 1.0],
            2,
            1,
            0,
            "max",
            "target";
            allow_densify=true,
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1),
            PF.HomogeneousPolynomial([1.0, 2.0], 2, 1);
            allow_densify=true,
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1),
            PF.HomogeneousPolynomial(ComplexF64[1, 0, 1], 2, 2);
            allow_densify=true,
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1),
            PF.HomogeneousPolynomial(BigFloat[1, 0, 1], 2, 2);
            allow_densify=true,
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; sense=:largest, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; target=Inf, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; inner_samples=-1, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1),
            quadratic;
            inner_samples=2,
            max_samples=1,
            allow_densify=true,
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; level=3, max_dimension=3, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; max_exponent_entries=3, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; max_dense_entries=3, allow_densify=true
        )
        @test_throws ArgumentError PF.polynomial_bounds(
            MersenneTwister(1), quadratic; max_work=1, allow_densify=true
        )
    end
end
