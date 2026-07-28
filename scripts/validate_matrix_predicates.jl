#!/usr/bin/env julia

using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const QET = QuantumEntanglementTools

function cauchy_matrix(rows::Int, columns::Int, row_shift::Int, column_shift::Int)
    return Rational{BigInt}[
        1 // BigInt(row_shift + row + column_shift + column) for
        row in 1:rows, column in 1:columns
    ]
end

function compound_minor_values(matrix)
    return (vec(matrix), vec(QET.compound_matrix(matrix, 2)))
end

rng = MersenneTwister(0x5145544d50524544)

@testset "independent randomized matrix-predicate cross-checks" begin
    @testset "exact PSD versus eigenspectrum" begin
        for sample in 1:50
            factor = Rational{BigInt}.(rand(rng, -4:4, 3, 3))
            matrix = factor' * factor
            if iseven(sample)
                matrix[1, 1] = -one(eltype(matrix))
            end

            result = QET.is_positive_semidefinite(matrix)
            smallest = eigmin(Hermitian(Float64.(matrix)))
            expected = if smallest < -1.0e-10
                QET.MatrixPredicateViolated
            else
                QET.MatrixPredicateSatisfied
            end
            @test result.status === expected
        end
    end

    @testset "total positivity versus compound minors" begin
        for sample in 1:20
            row_shift = rand(rng, 0:5)
            column_shift = rand(rng, 5:10)
            matrix = cauchy_matrix(2, 3, row_shift, column_shift)
            if iseven(sample)
                matrix[1, 1] = zero(eltype(matrix))
            end

            minor_values = compound_minor_values(matrix)
            expected = all(value -> value > 0, Iterators.flatten(minor_values))
            result = QET.is_totally_positive(matrix; orders=(1, 2))
            @test result.status ===
                (expected ? QET.MatrixPredicateSatisfied : QET.MatrixPredicateViolated)
            @test result.planned == length(minor_values[1]) + length(minor_values[2])
        end
    end

    @testset "total nonsingularity versus compound minors" begin
        for sample in 1:20
            row_shift = rand(rng, 0:5)
            column_shift = rand(rng, 5:10)
            matrix = cauchy_matrix(2, 3, row_shift, column_shift)
            if iseven(sample)
                matrix[:, 3] = matrix[:, 2]
            end

            minor_values = compound_minor_values(matrix)
            expected = all(value -> !iszero(value), Iterators.flatten(minor_values))
            result = QET.is_totally_nonsingular(matrix; orders=(1, 2))
            @test result.status ===
                (expected ? QET.MatrixPredicateSatisfied : QET.MatrixPredicateViolated)
            @test result.planned == length(minor_values[1]) + length(minor_values[2])
        end
    end
end
