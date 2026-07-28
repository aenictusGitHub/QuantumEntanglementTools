using LinearAlgebra
using Random
using SparseArrays
using Test

function _tier_e_combination_products(values, order)
    order == 0 && return one(eltype(values))
    result = zero(eltype(values))
    selected = Vector{Int}(undef, order)
    function visit(position, next_index)
        if position > order
            term = one(eltype(values))
            for index in selected
                term *= values[index]
            end
            result += term
            return nothing
        end
        final_index = length(values) - (order - position)
        for index in next_index:final_index
            selected[position] = index
            visit(position + 1, index + 1)
        end
    end
    visit(1, 1)
    return result
end

function _tier_e_pair_sums(values)
    result = eltype(values)[]
    for second in 2:length(values), first in 1:(second - 1)
        push!(result, values[first] + values[second])
    end
    return sort(result)
end

@testset "Tier E matrix analysis" begin
    @testset "strong vector and singular-value majorization" begin
        @test majorizes([4, 1, 1], [3, 2, 1]; rtol=0)
        @test !majorizes([3, 2, 1], [4, 1, 1]; rtol=0)
        @test majorizes([1, -1], [1, 0, -1]; rtol=0)
        @test majorizes([1, 0, -1], [1, -1]; rtol=0)
        @test majorizes(Int[], [0, 0]; rtol=0)
        @test !majorizes(Int[], [0, 1]; rtol=0)
        @test !majorizes([2, 0], [1, 0]; rtol=0)
        @test majorizes([3 // 2, 1 // 2], [1, 1]; rtol=0)
        @test majorizes([1, 4, 1], [2, 1, 3]; rtol=0)
        largest = typemax(Int)
        @test majorizes([largest, largest, -largest], [largest, 0, 0]; rtol=0)

        first_matrix = Diagonal([3.0, 1.0])
        second_matrix = Diagonal([2.0, 2.0])
        @test majorizes(first_matrix, second_matrix; rtol=0)
        @test !majorizes(second_matrix, first_matrix; rtol=0)
        @test majorizes(first_matrix, [2.0, 2.0]; rtol=0)
        @test majorizes(Diagonal(ComplexF64[3im, -1]), second_matrix; rtol=0)
        @test !majorizes(Diagonal([3.0, 1.0]), Diagonal([2.0, 1.0]); rtol=0)
        @test majorizes([1.0 1.0], [sqrt(2.0)]; atol=1e-14, rtol=0)

        delta = 1e-12
        almost = [0.6 - delta, 0.3 + delta, 0.1]
        reference = [0.6, 0.3, 0.1]
        @test !majorizes(almost, reference; atol=0, rtol=0)
        @test majorizes(almost, reference; atol=2delta, rtol=0)
        @test majorizes(reference, almost; atol=0, rtol=0)

        @test_throws ArgumentError majorizes([1 + im, 2], [2, 1])
        @test_throws ArgumentError majorizes([NaN, 1.0], [1.0, 0.0])
        @test_throws ArgumentError majorizes([Inf 0.0; 0.0 1.0], Matrix{Float64}(I, 2, 2))
        @test_throws ArgumentError majorizes(sparsevec([1], [1.0], 2), [1, 0])
        @test_throws ArgumentError majorizes(
            sparse(Matrix{Float64}(I, 2, 2)), Matrix{Float64}(I, 2, 2)
        )
        @test_throws ArgumentError majorizes(BigFloat[1 0; 0 2], BigFloat[1 0; 0 2])
        @test_throws ArgumentError majorizes([1, 0], [1, 0]; atol=-1)
        @test_throws ArgumentError majorizes([1, 0], [1, 0]; rtol=Inf)
        @test_throws ArgumentError majorizes([1, 0], [1, 0]; atol=true)
    end

    @testset "elementary symmetric polynomials" begin
        values = [1, 2, 3, 4]
        @test elementary_symmetric_polynomial(values, 0) == 1
        @test elementary_symmetric_polynomial(values, 1) == 10
        @test elementary_symmetric_polynomial(values, 2) == 35
        @test elementary_symmetric_polynomial(values, 3) == 50
        @test elementary_symmetric_polynomial(values, 4) == 24
        @test elementary_symmetric_polynomial(Int[], 0) == 1

        rationals = Rational{Int}[1 // 2, 2 // 3, 3 // 4]
        rational_result = elementary_symmetric_polynomial(rationals, 2)
        @test rational_result == 29 // 24
        @test rational_result isa Rational{Int}
        @test elementary_symmetric_polynomial(rationals, 0) === 1 // 1
        @test elementary_symmetric_polynomial(Complex{Int}[1 + im, 2 - im, -1 + 2im], 3) ==
            prod(Complex{Int}[1 + im, 2 - im, -1 + 2im])

        sparse_values = sparsevec([1, 4], Rational{Int}[2 // 3, -3 // 2], 5)
        @test elementary_symmetric_polynomial(sparse_values, 1) == -5 // 6
        @test elementary_symmetric_polynomial(sparse_values, 2) == -1
        @test elementary_symmetric_polynomial(sparse_values, 3) == 0

        property_values = Rational{BigInt}.([2, -1, 3, 5, -4])
        for order in 0:length(property_values)
            @test elementary_symmetric_polynomial(property_values, order) ==
                _tier_e_combination_products(property_values, order)
        end
        for order in 0:length(values)
            @test elementary_symmetric_polynomial(3 .* values, order) ==
                3^order * elementary_symmetric_polynomial(values, order)
        end

        first = Rational{Int}[1 // 2, 2, -1]
        second = Rational{Int}[3 // 2, -2]
        for order in 0:(length(first) + length(second))
            convolution = sum(
                elementary_symmetric_polynomial(first, first_order) *
                elementary_symmetric_polynomial(second, order - first_order) for
                first_order in max(0, order - length(second)):min(order, length(first))
            )
            @test elementary_symmetric_polynomial(vcat(first, second), order) == convolution
        end
        @test elementary_symmetric_polynomial(@view(values[2:4]), 2) == 26
        @test elementary_symmetric_polynomial(Number[1, 2 // 3, 0.5], 2) ≈ 3 / 2

        @test_throws ArgumentError elementary_symmetric_polynomial(values, -1)
        @test_throws ArgumentError elementary_symmetric_polynomial(values, 5)
        @test_throws ArgumentError elementary_symmetric_polynomial(values, 1.5)
        @test_throws ArgumentError elementary_symmetric_polynomial(values, true)
        @test_throws ArgumentError elementary_symmetric_polynomial(Any[1, "x"], 1)
        @test_throws ArgumentError elementary_symmetric_polynomial([1.0, NaN], 1)
        @test_throws OverflowError elementary_symmetric_polynomial([typemax(Int), 2], 2)
    end

    @testset "multiplicative compound matrices" begin
        matrix = [1 2 3; 4 5 6; 7 8 10]
        @test compound_matrix(matrix, 0) == reshape([1], 1, 1)
        @test compound_matrix(matrix, 1) == matrix
        @test compound_matrix(matrix, 2) == [-3 -6 -3; -6 -11 -4; -3 -2 2]
        @test compound_matrix(matrix, 3) == reshape([-3], 1, 1)
        @test eltype(compound_matrix(matrix, 2)) == Int

        rectangular = [1 2; 3 5; 7 11]
        @test compound_matrix(rectangular, 2) == reshape([-1, -3, -2], 3, 1)
        @test size(compound_matrix(zeros(Int, 2, 3), 3)) == (0, 1)
        @test size(compound_matrix(zeros(Int, 2, 3), 4)) == (0, 0)
        @test size(compound_matrix(zeros(Int, 3, 2), 3)) == (1, 0)

        diagonal_values = Rational{Int}[2, -3, 5, 7]
        diagonal_compound = compound_matrix(Diagonal(diagonal_values), 2)
        @test diagonal_compound == Diagonal(
            Rational{Int}[
                diagonal_values[first] * diagonal_values[second] for first in 1:3 for
                second in (first + 1):4
            ],
        )
        @test eltype(diagonal_compound) == Rational{Int}

        first_factor = [
            1 2 0
            -1 3 4
            2 0 5
            3 -2 1
        ]
        second_factor = [
            2 1 0 -1
            0 3 2 1
            1 -2 4 2
        ]
        @test compound_matrix(first_factor * second_factor, 2) ==
            compound_matrix(first_factor, 2) * compound_matrix(second_factor, 2)
        @test compound_matrix(transpose(first_factor), 2) ==
            transpose(compound_matrix(first_factor, 2))

        rational_matrix = Rational{Int}[1//2 2//3; -3//4 5//7]
        @test compound_matrix(rational_matrix, 2)[1, 1] == det(rational_matrix)
        @test eltype(compound_matrix(rational_matrix, 2)) == Rational{Int}
        complex_matrix = Complex{Int}[1 + im 2; 3im 4 - im]
        @test compound_matrix(complex_matrix, 2)[1, 1] ==
            complex_matrix[1, 1] * complex_matrix[2, 2] -
              complex_matrix[1, 2] * complex_matrix[2, 1]
        @test eltype(compound_matrix(complex_matrix, 2)) == Complex{Int}
        @test compound_matrix(BigInt[typemax(Int) 0; 0 2], 2)[1, 1] ==
            BigInt(typemax(Int)) * 2
        @test eltype(compound_matrix(UInt[1 2; 3 4], 2)) == BigInt
        @test compound_matrix(Matrix{Number}([1 2; 3 5]), 2) == reshape([-1], 1, 1)

        @test_throws OverflowError compound_matrix([typemax(Int) 0; 0 2], 2)
        sparse_compound = compound_matrix(sparse(matrix), 2)
        @test issparse(sparse_compound)
        @test Matrix(sparse_compound) == compound_matrix(matrix, 2)
        @test !issparse(compound_matrix(sparse(matrix), 2; sparse_output=false))
        @test issparse(compound_matrix(matrix, 2; sparse_output=true))
        @test_throws ArgumentError compound_matrix(matrix, -1)
        @test_throws ArgumentError compound_matrix(matrix, 1.5)
        @test_throws ArgumentError compound_matrix(matrix, true)
        @test_throws ArgumentError compound_matrix(Any[1 "x"; 2 3], 1)
        @test_throws ArgumentError compound_matrix([1.0 Inf; 2.0 3.0], 1)
    end

    @testset "additive compound matrices" begin
        matrix = [1 2 3; 4 5 6; 7 8 10]
        expected_second = [6 6 -3; 8 11 2; -7 4 15]
        @test additive_compound_matrix(matrix, 0) == reshape([0], 1, 1)
        @test additive_compound_matrix(matrix, 1) == matrix
        @test additive_compound_matrix(matrix, 2) == expected_second
        @test additive_compound_matrix(matrix, 3) == reshape([tr(matrix)], 1, 1)
        @test size(additive_compound_matrix(matrix, 4)) == (0, 0)
        @test eltype(additive_compound_matrix(matrix, 2)) == Int

        diagonal_values = [-2.0, 0.5, 3.0, 7.0]
        diagonal_additive = additive_compound_matrix(Diagonal(diagonal_values), 2)
        lexicographic_pair_sums = [
            diagonal_values[first] + diagonal_values[second] for first in 1:3 for
            second in (first + 1):4
        ]
        @test diagonal_additive == Matrix(Diagonal(lexicographic_pair_sums))
        @test sort(eigvals(diagonal_additive)) == _tier_e_pair_sums(diagonal_values)

        for order in 1:4
            exact_matrix = reshape(collect(1:16), 4, 4)
            @test tr(additive_compound_matrix(exact_matrix, order)) ==
                binomial(3, order - 1) * tr(exact_matrix)
        end

        first = Rational{Int}.([1 2 0; -1 3 4; 2 0 5])
        second = Rational{Int}.([2 -1 3; 0 4 1; -2 5 2])
        @test additive_compound_matrix(first + second, 2) ==
            additive_compound_matrix(first, 2) + additive_compound_matrix(second, 2)
        @test eltype(additive_compound_matrix(first, 2)) == Rational{Int}

        rng = MersenneTwister(20260728)
        floating_matrix = randn(rng, 4, 4)
        step = 1e-5
        identity4 = Matrix{Float64}(I, 4, 4)
        central_difference =
            (
                compound_matrix(identity4 + step * floating_matrix, 3) -
                compound_matrix(identity4 - step * floating_matrix, 3)
            ) / (2step)
        @test central_difference ≈ additive_compound_matrix(floating_matrix, 3) atol = 3e-9 rtol =
            3e-9

        symmetric_matrix = Symmetric(randn(rng, 4, 4))
        source_eigenvalues = eigvals(symmetric_matrix)
        additive_eigenvalues = sort(
            eigvals(Symmetric(additive_compound_matrix(symmetric_matrix, 2)))
        )
        @test additive_eigenvalues ≈ _tier_e_pair_sums(source_eigenvalues) atol = 3e-12 rtol =
            3e-12

        complex_matrix = Complex{Rational{Int}}[
            1//2+im 2//3
            -3//4*im 5//7-im
        ]
        @test additive_compound_matrix(complex_matrix, 2) ==
            reshape([tr(complex_matrix)], 1, 1)
        @test eltype(additive_compound_matrix(complex_matrix, 2)) == Complex{Rational{Int}}
        @test eltype(additive_compound_matrix(UInt[1 2; 3 4], 1)) == BigInt

        @test_throws OverflowError additive_compound_matrix([typemax(Int) 0; 0 1], 2)
        @test_throws DimensionMismatch additive_compound_matrix(zeros(2, 3), 1)
        sparse_additive = additive_compound_matrix(sparse(matrix), 2)
        @test issparse(sparse_additive)
        @test Matrix(sparse_additive) == expected_second
        @test !issparse(additive_compound_matrix(sparse(matrix), 2; sparse_output=false))
        @test issparse(additive_compound_matrix(matrix, 2; sparse_output=true))
        @test_throws ArgumentError additive_compound_matrix(matrix, -1)
        @test_throws ArgumentError additive_compound_matrix(matrix, 1.5)
        @test_throws ArgumentError additive_compound_matrix(matrix, true)
        @test_throws ArgumentError additive_compound_matrix(Any[1 "x"; 2 3], 1)
        @test_throws ArgumentError additive_compound_matrix(ComplexF64[1 0; 0 Inf * im], 1)
    end
end
