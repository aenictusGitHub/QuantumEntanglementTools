using Test
using QuantumEntanglementTools

const Robk = QuantumEntanglementTools
const RobkCompat = Robk.MATLABCompat

struct _RobkZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_RobkZeroBasedVector) = size(vector.storage)
Base.axes(vector::_RobkZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_RobkZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_RobkZeroBasedVector, index::Int) = vector.storage[index + 1]

function _robk_reference(sorted_magnitudes, k)
    n = length(sorted_magnitudes)
    branch = 1
    for candidate in k:-1:2
        tail = sum(@view sorted_magnitudes[candidate:n])
        if sorted_magnitudes[candidate - 1] >= tail / (k - candidate + 1)
            branch = candidate
            break
        end
    end
    tail = @view sorted_magnitudes[branch:n]
    return sum(tail)^2 / (k - branch + 1) - sum(abs2, tail), branch
end

@testset "RobkCohValue theorem-backed pure-state formula" begin
    @testset "closed-form values and branch selection" begin
        basis = ComplexF64[1, 0, 0, 0]
        basis_result = Robk.pure_k_coherence_robustness(basis, 2)
        @test basis_result.value == 0.0
        @test basis_result.branch_index == 2
        @test basis_result.branch_status == :stable
        @test basis_result.selected_gap == 1.0
        @test basis_result.next_gap === nothing

        maximally_coherent = fill(ComplexF64(1 / 2), 4)
        maximum_result = Robk.pure_k_coherence_robustness(maximally_coherent, 2)
        @test maximum_result.value ≈ 1.0
        @test maximum_result.branch_index == 1
        @test maximum_result.selected_gap === nothing
        @test maximum_result.next_gap > 0
        @test maximum_result.branch_status == :stable

        rank_two = ComplexF64[inv(sqrt(2)), im / sqrt(2), 0, 0]
        rank_two_result = Robk.pure_k_coherence_robustness(rank_two, 2)
        @test rank_two_result.value == 0.0
        @test rank_two_result.branch_index == 2

        state = ComplexF64[0.13 + 0.19im, -0.47im, 0.31, -0.28 + 0.08im]
        state ./= sqrt(sum(abs2, state))
        for k in 2:length(state)
            result = Robk.pure_k_coherence_robustness(state, k)
            expected, expected_branch = _robk_reference(sort(abs.(state); rev=true), k)
            @test result.value ≈ expected atol = 8eps(Float64)
            @test result.branch_index == expected_branch
            @test result.tail_average == result.tail_sum / (k - expected_branch + 1)
            @test occursin("branch_index=$(result.branch_index)", sprint(show, result))
        end
    end

    @testset "complex phases and coordinate order are immaterial" begin
        magnitudes = [0.7, 0.5, 0.4, sqrt(0.1)]
        magnitudes ./= sqrt(sum(abs2, magnitudes))
        phases = cis.([0.3, -1.1, 2.4, 0.8])
        complex_state = magnitudes .* phases
        permuted = complex_state[[3, 1, 4, 2]]

        for k in 2:4
            canonical = Robk.pure_k_coherence_robustness(magnitudes, k)
            phased = Robk.pure_k_coherence_robustness(complex_state, k)
            reordered = Robk.pure_k_coherence_robustness(permuted, k)
            @test phased.value ≈ canonical.value
            @test reordered.value ≈ canonical.value
            @test phased.branch_index == canonical.branch_index
            @test reordered.branch_index == canonical.branch_index
        end
    end

    @testset "equality and adjacent branch boundaries" begin
        # Before the common normalization factor, a1 = (a2 + a3 + a4) / 2,
        # while a2 < a3 + a4. The theorem therefore selects ell = 2 at the
        # equality rather than falling through to ell = 1.
        equality_state = Float64[0.5, 0.375, 0.375, 0.25]
        equality_state ./= sqrt(sum(abs2, equality_state))
        equality_result = Robk.pure_k_coherence_robustness(
            equality_state, 3; atol=4eps(Float64), rtol=0
        )
        @test equality_result.branch_index == 2
        @test equality_result.branch_status == :exact_equality
        @test equality_result.selected_gap == 0

        delta = 32eps(Float64)
        upper = copy(equality_state)
        upper[1] += delta
        upper ./= sqrt(sum(abs2, upper))
        lower = copy(equality_state)
        lower[1] -= delta
        lower ./= sqrt(sum(abs2, lower))
        upper_result = Robk.pure_k_coherence_robustness(upper, 3; atol=1e-12, rtol=0)
        lower_result = Robk.pure_k_coherence_robustness(lower, 3; atol=1e-12, rtol=0)
        @test upper_result.branch_index == 2
        @test lower_result.branch_index == 1
        @test upper_result.branch_status == :near_boundary
        @test lower_result.branch_status == :near_boundary
        @test upper_result.value ≈ equality_result.value atol = 1e-12
        @test lower_result.value ≈ equality_result.value atol = 1e-12

        exact_high_branch = Float64[0.75, 0.5, 0.375, 0.125]
        exact_high_branch ./= sqrt(sum(abs2, exact_high_branch))
        high_result = Robk.pure_k_coherence_robustness(
            exact_high_branch, 3; atol=4eps(Float64), rtol=0
        )
        @test high_result.branch_index == 3
        @test high_result.branch_status == :exact_equality
        @test high_result.selected_gap == 0
        @test high_result.next_gap === nothing
    end

    @testset "validation never repairs the input" begin
        normalized = ComplexF64[1, 1] / sqrt(2)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(ComplexF64[], 2)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(ComplexF64[1, NaN], 2)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness([1, 0], 2)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(
            _RobkZeroBasedVector(normalized), 2
        )
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(0.9normalized, 2)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(normalized, 2; atol=-1)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(normalized, 2; rtol=Inf)

        @test_throws ArgumentError Robk.pure_k_coherence_robustness(normalized, true)
        @test_throws ArgumentError Robk.pure_k_coherence_robustness(normalized, 2.0)
        @test_throws DomainError Robk.pure_k_coherence_robustness(normalized, 1)
        @test_throws DomainError Robk.pure_k_coherence_robustness(normalized, 3)

        almost_normalized = normalized .* (1 + 2eps(Float64))
        accepted = Robk.pure_k_coherence_robustness(
            almost_normalized, 2; atol=8eps(Float64), rtol=0
        )
        @test accepted.norm_squared != 1
        @test accepted.normalization_residual > 0
        @test accepted.value == 0
        @test almost_normalized == normalized .* (1 + 2eps(Float64))
    end

    @testset "meaningful floating precision is preserved" begin
        state32 = ComplexF32[1, im] / sqrt(2.0f0)
        result32 = Robk.pure_k_coherence_robustness(state32, 2)
        @test result32.value isa Float32
        @test result32.tail_sum isa Float32
        @test result32.branch_tolerance isa Float32

        setprecision(BigFloat, 192) do
            state_big = Complex{BigFloat}[1 + 2im, 3 - im, -2 + 4im]
            state_big ./= sqrt(sum(abs2, state_big))
            result_big = Robk.pure_k_coherence_robustness(
                state_big, 2; atol=big"1e-55", rtol=big"1e-50"
            )
            expected, branch = _robk_reference(sort(abs.(state_big); rev=true), 2)
            @test result_big.value isa BigFloat
            @test result_big.value ≈ expected rtol = big"1e-50"
            @test result_big.branch_index == branch
        end
    end

    @testset "two-output QETLAB compatibility" begin
        state = ComplexF64[0.13 + 0.19im, -0.47im, 0.31, -0.28 + 0.08im]
        state ./= sqrt(sum(abs2, state))
        native = Robk.pure_k_coherence_robustness(state, 3)

        value, branch = RobkCompat.RobkCohValue(state, 3)
        @test value ≈ native.value
        @test branch == native.branch_index
        @test RobkCompat.RobkCohValue(reshape(state, 1, :), 3) == (value, branch)
        @test RobkCompat.RobkCohValue(reshape(state, :, 1), 3) == (value, branch)

        unsorted = state[[4, 1, 3, 2]]
        phased = abs.(state) .* cis.([0.3, -1.1, 2.4, 0.8])
        phased ./= sqrt(sum(abs2, phased))
        @test RobkCompat.RobkCohValue(unsorted, 3) == (value, branch)
        phased_value, phased_branch = RobkCompat.RobkCohValue(phased, 3)
        magnitude_value, magnitude_branch = RobkCompat.RobkCohValue(abs.(state), 3)
        @test phased_value ≈ magnitude_value
        @test phased_branch == magnitude_branch

        @test_throws ArgumentError RobkCompat.RobkCohValue(0.9state, 3)
        @test_throws DomainError RobkCompat.RobkCohValue(state, 1)
        @test_throws DimensionMismatch RobkCompat.RobkCohValue(ones(2, 2), 2)
    end

    @testset "result constructor enforces diagnostic invariants" begin
        @test_throws ArgumentError Robk.PureKCoherenceRobustnessResult(
            -1.0, 2, 1, :stable, nothing, 0.1, 1e-8, 1.0, 0.5, 1.0, 0.0, 1e-8
        )
        @test_throws ArgumentError Robk.PureKCoherenceRobustnessResult(
            0.0, 2, 1, :invented, nothing, 0.1, 1e-8, 1.0, 0.5, 1.0, 0.0, 1e-8
        )
    end
end
