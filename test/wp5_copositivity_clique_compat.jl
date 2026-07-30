using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const CopCliqueCompat = QuantumEntanglementTools.MATLABCompat

function wp5_compat_cycle(vertices::Int)
    adjacency = zeros(Int, vertices, vertices)
    for vertex in 1:vertices
        neighbor = mod1(vertex + 1, vertices)
        adjacency[vertex, neighbor] = 1
        adjacency[neighbor, vertex] = 1
    end
    return adjacency
end

@testset "WP5 copositivity and clique-number compatibility mappings" begin
    positive = CopCliqueCompat.IsCopositive(
        MersenneTwister(0xc001), Rational{Int}[1 -1; -1 1]
    )
    @test positive isa CopositivityResult
    @test positive.verdict === true
    @test positive.certified
    @test CopCliqueCompat.IsCopositive(
        MersenneTwister(0xc002), Rational{Int}[1 -1; -1 1], 0, :sos; structured=false
    ) == 1

    negative = CopCliqueCompat.IsCopositive(MersenneTwister(0xc003), [1 -2; -2 1])
    @test negative.verdict === false
    @test negative.witness.value == -1 // 2
    @test CopCliqueCompat.IsCopositive(
        MersenneTwister(0xc004), [1 -2; -2 1], 0, "SOS"; structured=false
    ) == 0

    boundary = [-5.0e-10 0.0; 0.0 1.0]
    @test_throws DomainError CopCliqueCompat.IsCopositive(
        MersenneTwister(0xc005), boundary; structured=false
    )
    @test_throws ArgumentError CopCliqueCompat.IsCopositive(
        MersenneTwister(0xc006), Matrix{Float64}(I, 2, 2), 0, "invalid"
    )
    @test !hasmethod(CopCliqueCompat.IsCopositive, Tuple{Matrix{Int}})

    cycle = wp5_compat_cycle(5)
    bounds = CopCliqueCompat.CliqueNumber(MersenneTwister(0xc007), cycle)
    @test bounds isa CliqueNumberResult
    @test bounds.bounds_certified
    @test (bounds.lower_bound, bounds.upper_bound) == (2, 3)
    @test bounds.uncertified_upper_candidate === nothing
    @test CopCliqueCompat.CliqueNumber(
        MersenneTwister(0xc008), cycle, 0, :sos; structured=false
    ) == (ub=3, lb=2)

    complete = ones(Int, 4, 4) - Matrix{Int}(I, 4, 4)
    @test CopCliqueCompat.CliqueNumber(
        MersenneTwister(0xc009), complete; structured=false
    ) == (ub=4, lb=4)
    @test_throws ArgumentError CopCliqueCompat.CliqueNumber(
        MersenneTwister(0xc00a), cycle, 0, "bad-mode"
    )
    @test_throws ArgumentError CopCliqueCompat.CliqueNumber(
        MersenneTwister(0xc00b), [0 1; 0 0]
    )
    @test !hasmethod(CopCliqueCompat.CliqueNumber, Tuple{Matrix{Int}})

    Random.seed!(0xc0ffee)
    expected_after_copositivity = rand(UInt64)
    expected_after_clique = rand(UInt64)
    Random.seed!(0xc0ffee)
    CopCliqueCompat.IsCopositive(MersenneTwister(0xc00c), Rational{Int}[1 -1; -1 1])
    @test rand(UInt64) == expected_after_copositivity
    CopCliqueCompat.CliqueNumber(MersenneTwister(0xc00d), cycle)
    @test rand(UInt64) == expected_after_clique
end
