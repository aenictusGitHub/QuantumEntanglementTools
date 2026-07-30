using LinearAlgebra
using Random
using Test

using QuantumEntanglementTools
const SKNormCompat = QuantumEntanglementTools.MATLABCompat

@testset "WP7 S(k) compatibility mappings" begin
    variable = complex_affine_variable(:X_compat, 2, 2)
    atom = SKNormCompat.kpNormDual(variable, 2, 2)
    @test atom isa TopKPNormDualEpigraph
    @test atom.k == 2
    @test atom.p == 2

    diagonal = Diagonal([4.0, 3.0, 2.0, 1.0])
    exact_rng = MersenneTwister(0x5a01)
    exact_control = copy(exact_rng)
    exact = SKNormCompat.SkOperatorNorm(exact_rng, diagonal, 2, (2, 2), 0)
    @test exact isa SKOperatorNormResult
    @test exact.exact
    @test exact.lower_bound == exact.upper_bound == 4
    @test rand(exact_rng) == rand(exact_control)

    legacy = SKNormCompat.SkOperatorNorm(
        MersenneTwister(0x5a02), diagonal, 2, (2, 2), 0; structured=false
    )
    @test legacy.lb == legacy.ub == 4
    @test legacy.lwit.validated
    @test legacy.uwit.exact

    violation_rng = MersenneTwister(0x5a03)
    violation_control = copy(violation_rng)
    violation = SKNormCompat.IsBlockPositive(
        violation_rng, Diagonal([-1.0, 2.0, 2.0, 2.0]), 1, (2, 2), 0
    )
    @test violation isa BlockPositivityResult
    @test violation.verdict === false
    @test violation.witness.validated
    @test rand(violation_rng) == rand(violation_control)

    legacy_violation = SKNormCompat.IsBlockPositive(
        MersenneTwister(0x5a04),
        Diagonal([-1.0, 2.0, 2.0, 2.0]),
        1,
        (2, 2),
        0;
        structured=false,
    )
    @test legacy_violation.ibp == 0
    @test legacy_violation.wit !== nothing

    positive = SKNormCompat.IsBlockPositive(
        MersenneTwister(0x5a05), Matrix{Float64}(I, 4, 4), 1, (2, 2), 0; structured=false
    )
    @test positive == (ibp=1, wit=nothing)

    boundary = SKNormCompat.IsBlockPositive(
        MersenneTwister(0x5a06),
        Diagonal([-1.0e-10, 1.0, 1.0, 1.0]),
        2,
        (2, 2),
        0,
        1.0e-9;
        structured=false,
    )
    @test boundary == (ibp=-1, wit=nothing)

    @test_throws ArgumentError SKNormCompat.SkOperatorNorm(
        MersenneTwister(1), diagonal, 1, (2, 2), 0, -1, -1
    )
    @test_throws DimensionMismatch SKNormCompat.IsBlockPositive(
        MersenneTwister(1), diagonal, 1, (3, 2)
    )
end
