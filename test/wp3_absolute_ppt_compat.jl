using Test
using QuantumEntanglementTools

const AbsPPTCompatQET = QuantumEntanglementTools
const AbsPPTCompat = AbsPPTCompatQET.MATLABCompat

if !isdefined(AbsPPTCompatQET, :abs_ppt_constraints)
    Base.include(
        AbsPPTCompatQET, joinpath(@__DIR__, "..", "src", "entanglement", "absolute_ppt.jl")
    )
end

if !isdefined(AbsPPTCompat, :AbsPPTConstraints)
    Core.eval(
        AbsPPTCompat,
        quote
            using ..QuantumEntanglementTools:
                AbsPPTEnumerationConstraintLimit,
                AbsPPTEnumerationEarlyViolation,
                AffineScalar,
                abs_ppt_constraints,
                is_abs_ppt
        end,
    )
    Base.include(AbsPPTCompat, joinpath(@__DIR__, "..", "src", "compat", "absolute_ppt.jl"))
end

@testset "WP3 absolute-PPT MATLAB compatibility" begin
    matrices2 = AbsPPTCompat.AbsPPTConstraints([0.4, 0.3, 0.2, 0.1], [2, 2])
    @test length(matrices2) == 1
    @test matrices2[1] ≈ [0.2 -0.2; -0.2 0.6]

    spectrum4 = collect((16:-1:1) ./ sum(1:16))
    limited = AbsPPTCompat.AbsPPTConstraints(spectrum4, [4, 4], 0, 3)
    @test length(limited) == 3
    structured_limited = AbsPPTCompat.AbsPPTConstraints(
        spectrum4, [4, 4], 0, 3; structured=true
    )
    @test structured_limited.status === AbsPPTCompatQET.AbsPPTEnumerationConstraintLimit
    @test length(structured_limited.constraints) == 3

    early = AbsPPTCompat.AbsPPTConstraints([1.0; zeros(15)], [4, 4], 1)
    @test length(early) == 1
    structured_early = AbsPPTCompat.AbsPPTConstraints(
        [1.0; zeros(15)], [4, 4], 1; structured=true
    )
    @test structured_early.status === AbsPPTCompatQET.AbsPPTEnumerationEarlyViolation

    @test_throws DomainError AbsPPTCompat.AbsPPTConstraints(
        collect(9.0:-1:1), [3, 3]; max_work=10
    )
    work_limited = AbsPPTCompat.AbsPPTConstraints(
        collect(9.0:-1:1), [3, 3]; max_work=10, structured=true
    )
    @test work_limited.status === AbsPPTCompatQET.AbsPPTEnumerationWorkLimit

    # The two pinned routines give scalar DIM different meanings.
    @test_throws DimensionMismatch AbsPPTCompat.AbsPPTConstraints(ones(6), 2)
    @test AbsPPTCompat.IsAbsPPT(fill(1 / 6, 6), 2) == 1
    @test AbsPPTCompat.IsAbsPPT(fill(1 / 6, 6)) == 1

    @test AbsPPTCompat.IsAbsPPT(fill(0.25, 4), [2, 2]) == 1
    @test AbsPPTCompat.IsAbsPPT([1.0, 0.0, 0.0, 0.0], [2, 2]) == 0
    @test AbsPPTCompat.IsAbsPPT([0.5, 1 / 6, 1 / 6, 1 / 6], [2, 2]) == -1
    structured_result = AbsPPTCompat.IsAbsPPT([1.0, 0.0, 0.0, 0.0], [2, 2]; structured=true)
    @test structured_result.status === AbsPPTCompatQET.AbsolutePPTCertifiedNot
    @test structured_result.verdict === false

    affine_spectrum = [
        AbsPPTCompatQET.AffineScalar(0.0, [index == variable ? 1.0 : 0.0 for index in 1:4])
        for variable in 1:4
    ]
    affine_matrices = AbsPPTCompat.AbsPPTConstraints(affine_spectrum, [2, 2])
    @test length(affine_matrices) == 1
    @test affine_matrices[1] isa AbsPPTCompatQET.HermitianAffineMatrix
    @test_throws ArgumentError AbsPPTCompat.AbsPPTConstraints(affine_spectrum, [2, 2], 1)

    @test_throws ArgumentError AbsPPTCompat.AbsPPTConstraints(ones(4), [2, 2], 2)
    @test_throws ArgumentError AbsPPTCompat.AbsPPTConstraints(ones(4), [2, 2], 0, -1)
    @test_throws ArgumentError AbsPPTCompat.AbsPPTConstraints(
        ones(4), [2, 2], 0, BigInt(typemax(Int)) + 1
    )
    @test_throws ArgumentError AbsPPTCompat.IsAbsPPT(ones(4), [2, 2]; structured=2)
end
