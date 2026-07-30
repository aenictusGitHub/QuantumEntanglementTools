using LinearAlgebra
using QuantumEntanglementTools
using Random
using Test

const SeparabilityCompatQET = QuantumEntanglementTools
const SeparabilityCompat = QuantumEntanglementTools.MATLABCompat

if !isdefined(SeparabilityCompatQET, :StateDiscriminationResult)
    Base.include(
        SeparabilityCompatQET,
        joinpath(@__DIR__, "..", "src", "optimization", "state_discrimination.jl"),
    )
end
if !isdefined(SeparabilityCompatQET, :LocalDistinguishabilityResult)
    Base.include(
        SeparabilityCompatQET,
        joinpath(@__DIR__, "..", "src", "entanglement", "separability_optimization.jl"),
    )
end
if !isdefined(SeparabilityCompat, :_SEPARABILITY_OPTIMIZATION_COMPAT_LOADED)
    Core.eval(
        SeparabilityCompat,
        :(using ..QuantumEntanglementTools:
            local_distinguishability, is_separable, upb_sep_distinguishable),
    )
    Base.include(
        SeparabilityCompat,
        joinpath(@__DIR__, "..", "src", "compat", "separability_optimization.jl"),
    )
end

function wp6_compat_isotropic(weight)
    bell = zeros(Float64, 9, 9)
    indices = (1, 5, 9)
    bell[collect(indices), collect(indices)] .= 1 / 3
    state = (1 - weight) * Matrix{Float64}(I, 9, 9) / 9 + weight * bell
    state[end, end] += 1 - real(SeparabilityCompatQET._sepopt_represented_trace(state))
    return state
end

@testset "WP6 separability MATLAB compatibility" begin
    @testset "LocalDistinguishability safe legacy outputs" begin
        rho = Matrix{Float64}(I, 4, 4) / 4
        structured = SeparabilityCompat.LocalDistinguishability((rho,))
        @test structured isa SeparabilityCompatQET.LocalDistinguishabilityResult
        @test structured.status ===
            SeparabilityCompatQET.LocalDistinguishabilityTrivialOptimal

        legacy = SeparabilityCompat.LocalDistinguishability((rho,); structured=false)
        @test legacy.dist == 1
        @test length(legacy.meas) == 1
        @test legacy.meas[1] == Matrix{ComplexF64}(I, 4, 4)
        @test legacy.dual_sol == rho

        deterministic = SeparabilityCompat.LocalDistinguishability(
            (rho, rho), [0.0, 1.0], (2, 2); structured=false
        )
        @test deterministic.dist == 1
        @test iszero(deterministic.meas[1])
        @test deterministic.meas[2] == Matrix{ComplexF64}(I, 4, 4)
        @test deterministic.dual_sol == rho

        columns = hcat([1.0, 0, 0, 0], [0.0, 0, 0, 1])
        unavailable = SeparabilityCompat.LocalDistinguishability(columns, nothing, (2, 2))
        @test unavailable.status ===
            SeparabilityCompatQET.LocalDistinguishabilityBackendUnavailable
        @test_throws DomainError SeparabilityCompat.LocalDistinguishability(
            columns, nothing, (2, 2); structured=false
        )

        @test_throws ArgumentError SeparabilityCompat.LocalDistinguishability(
            hcat([2.0, 0, 0, 0], [0.0, 0, 0, 1]), nothing, (2, 2)
        )
        @test_throws ArgumentError SeparabilityCompat.LocalDistinguishability(
            (rho,), nothing, (2, 2), 2, 2
        )
        @test_throws ArgumentError SeparabilityCompat.LocalDistinguishability(
            (rho,), nothing, (2, 2), 2, 1, 1, -1
        )
    end

    @testset "IsSeparable explicit RNG and tri-state safety" begin
        bell = ComplexF64[
            0.5 0 0 0.5
            0 0 0 0
            0 0 0 0
            0.5 0 0 0.5
        ]
        entangled = SeparabilityCompat.IsSeparable(MersenneTwister(10), bell, (2, 2), 0, 0)
        @test entangled.status === :entangled
        @test entangled.certified
        @test SeparabilityCompat.IsSeparable(
            MersenneTwister(10), bell, (2, 2), 0, 0; structured=false
        ) == 0

        diagonal = Diagonal([0.4, 0.1, 0.2, 0.3])
        @test SeparabilityCompat.IsSeparable(
            MersenneTwister(11), diagonal, (2, 2), 0, 0; structured=false
        ) == 1

        rng = MersenneTwister(12)
        untouched = copy(rng)
        unknown = SeparabilityCompat.IsSeparable(
            rng, wp6_compat_isotropic(0.2), (3, 3), 0, 0
        )
        @test unknown.status === :unknown
        @test !unknown.certified
        @test rand(rng) == rand(untouched)
        @test_throws DomainError SeparabilityCompat.IsSeparable(
            MersenneTwister(12), wp6_compat_isotropic(0.2), (3, 3), 0, 0; structured=false
        )

        @test_throws MethodError SeparabilityCompat.IsSeparable(bell, (2, 2))
        @test_throws ArgumentError SeparabilityCompat.IsSeparable(
            MersenneTwister(13), bell, (2, 2), -1, 0
        )
        @test_throws ArgumentError SeparabilityCompat.IsSeparable(
            MersenneTwister(13), bell, (2, 2), 0, 2
        )
        @test_throws ArgumentError SeparabilityCompat.IsSeparable(
            MersenneTwister(13), bell, (2, 2), 0, 0, -1
        )
        @test_throws ArgumentError SeparabilityCompat.IsSeparable(
            MersenneTwister(13), bell, (2, 2), 7, 0
        )

        user_report = SeparabilityCompatQET.EntanglementReport(
            :separable,
            false,
            nothing,
            :user_claim,
            :user,
            (reconstruction=Matrix{Float64}(I, 4, 4) / 4,),
            SeparabilityCompatQET.EntanglementAttempt[],
            "unvalidated user claim",
        )
        @test_throws DomainError SeparabilityCompat._compat_is_separable_output(user_report)
    end

    @testset "UPBSepDistinguishable refuses numerical Boolean collapse" begin
        tiles = SeparabilityCompatQET.upb(:tiles)
        unavailable = SeparabilityCompat.UPBSepDistinguishable(tiles.local_factors...)
        @test unavailable isa SeparabilityCompatQET.UPBSeparableDiscriminationResult
        @test unavailable.status === :backend_unavailable
        @test unavailable.separably_distinguishable === nothing
        collapse_error = try
            SeparabilityCompat.UPBSepDistinguishable(tiles.local_factors...; structured=false)
            nothing
        catch error
            error
        end
        @test collapse_error isa DomainError
        @test occursin("floating-only implementation", sprint(showerror, collapse_error))
        @test occursin(
            "inconclusive floating-cone Farkas evidence", sprint(showerror, collapse_error)
        )

        state_limited = SeparabilityCompat.UPBSepDistinguishable(
            tiles.local_factors...; max_states=4
        )
        @test state_limited.status === :resource_limit
        @test state_limited.feasibility === :unknown
        @test state_limited.upb_analysis.reason === :state_limit
        @test_throws DomainError SeparabilityCompat.UPBSepDistinguishable(
            tiles.local_factors...; max_states=4, structured=false
        )

        invalid = SeparabilityCompat.UPBSepDistinguishable(
            Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 3, 3)
        )
        @test invalid.status === :invalid_input
        @test_throws DomainError SeparabilityCompat.UPBSepDistinguishable(
            Matrix{Float64}(I, 3, 3), Matrix{Float64}(I, 3, 3); structured=false
        )

        result_fields = ntuple(
            index -> getfield(unavailable, index), fieldcount(typeof(unavailable))
        )
        @test_throws MethodError SeparabilityCompatQET.UPBSeparableDiscriminationResult(
            result_fields...
        )
    end
end
