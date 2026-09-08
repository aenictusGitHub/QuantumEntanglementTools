using LinearAlgebra

struct _WP2ZeroBasedArray4{T,A<:Array{T,4}} <: AbstractArray{T,4}
    storage::A
end

Base.size(array::_WP2ZeroBasedArray4) = size(array.storage)
function Base.axes(array::_WP2ZeroBasedArray4)
    return ntuple(axis -> 0:(size(array.storage, axis) - 1), 4)
end
Base.IndexStyle(::Type{<:_WP2ZeroBasedArray4}) = IndexCartesian()
function Base.getindex(array::_WP2ZeroBasedArray4, indices::Vararg{Int,4})
    return array.storage[ntuple(axis -> indices[axis] + 1, 4)...]
end

@testset "WP2 parallel repetition" begin
    @testset "one-copy identity and input ownership" begin
        game = reshape(collect(1:24), 2, 3, 2, 2)
        repeated = parallel_repetition(game, 1)
        @test repeated == game
        @test repeated !== game
        game[1] = -100
        @test repeated[1] == 1

        compatibility = QuantumEntanglementTools.MATLABCompat.ParallelRepetition(
            repeated, 1
        )
        @test compatibility == repeated
        @test compatibility !== repeated
    end

    @testset "two-copy packed-index formula" begin
        game = reshape(collect(1:16), 2, 2, 2, 2)
        repeated = parallel_repetition(game, 2; max_entries=1_000)
        compatibility = QuantumEntanglementTools.MATLABCompat.ParallelRepetition(
            game, 2; max_entries=1_000
        )
        @test size(repeated) == (4, 4, 4, 4)
        @test compatibility == repeated

        pack(first, second, base) = (first - 1) * base + second
        for a1 in 1:2,
            a2 in 1:2, b1 in 1:2, b2 in 1:2, x1 in 1:2, x2 in 1:2, y1 in 1:2,
            y2 in 1:2

            @test repeated[
                pack(a1, a2, 2), pack(b1, b2, 2), pack(x1, x2, 2), pack(y1, y2, 2)
            ] == game[a1, b1, x1, y1] * game[a2, b2, x2, y2]
        end
    end

    @testset "rectangular, three-copy, and scalar games" begin
        rectangular = reshape(Rational{Int}.(1:12), 1, 2, 3, 2)
        repeated = parallel_repetition(
            rectangular, 3; max_entries=100_000, max_work=1_000_000
        )
        @test size(repeated) == (1, 8, 27, 8)
        @test eltype(repeated) == Rational{Int}

        pack3(first, second, third, base) =
            (first - 1) * base^2 + (second - 1) * base + third
        @test repeated[1, pack3(2, 1, 2, 2), pack3(3, 1, 2, 3), pack3(1, 2, 2, 2)] ==
            rectangular[1, 2, 3, 1] * rectangular[1, 1, 1, 2] * rectangular[1, 2, 2, 2]

        scalar_game = fill(Complex{BigFloat}(2, -1), 1, 1, 1, 1)
        scalar_repeated = parallel_repetition(scalar_game, 7; max_entries=1, max_work=7)
        @test size(scalar_repeated) == (1, 1, 1, 1)
        @test only(scalar_repeated) == only(scalar_game)^7
        @test eltype(scalar_repeated) == Complex{BigFloat}
    end

    @testset "checked fixed-width integer coefficients" begin
        safe_game = fill(Int8(3), 1, 1, 1, 1)
        safe_repeated = parallel_repetition(safe_game, 3; max_entries=1, max_work=3)
        @test only(safe_repeated) === Int8(27)
        @test eltype(safe_repeated) === Int8

        large_integer = typemax(Int)
        overflowing_game = fill(large_integer, 1, 1, 1, 1)
        overflowing_snapshot = copy(overflowing_game)
        @test only(parallel_repetition(overflowing_game, 1)) == large_integer
        @test_throws OverflowError parallel_repetition(
            overflowing_game, 2; max_entries=1, max_work=2
        )
        @test_throws OverflowError parallel_repetition(
            @view(overflowing_game[:, :, :, :]), 2; max_entries=1, max_work=2
        )
        @test_throws OverflowError QuantumEntanglementTools.MATLABCompat.ParallelRepetition(
            overflowing_game, 2; max_entries=1, max_work=2
        )
        @test overflowing_game == overflowing_snapshot

        exact_game = fill(big(large_integer), 1, 1, 1, 1)
        exact_repeated = parallel_repetition(exact_game, 2; max_entries=1, max_work=2)
        @test only(exact_repeated) == big(large_integer)^2
        @test eltype(exact_repeated) === BigInt
    end

    @testset "limits and validation" begin
        game = ones(Float32, 2, 2, 2, 2)
        @test_throws ArgumentError parallel_repetition(game, 0)
        @test_throws ArgumentError parallel_repetition(game, -1)
        @test_throws ArgumentError parallel_repetition(game, true)
        @test_throws MethodError parallel_repetition(game, 2.0)
        @test_throws ArgumentError QuantumEntanglementTools.MATLABCompat.ParallelRepetition(
            game, 0
        )
        @test_throws ArgumentError QuantumEntanglementTools.MATLABCompat.ParallelRepetition(
            game, 2.0
        )
        @test_throws ArgumentError parallel_repetition(game, big(typemax(Int)) + 1)
        @test_throws DimensionMismatch parallel_repetition(ones(2, 2, 2), 2)
        @test_throws ArgumentError parallel_repetition(zeros(0, 1, 1, 1), 1)
        @test_throws ArgumentError parallel_repetition(fill(NaN, 1, 1, 1, 1), 1)
        @test_throws ArgumentError parallel_repetition(fill(Inf, 1, 1, 1, 1), 1)
        @test_throws ArgumentError parallel_repetition(
            _WP2ZeroBasedArray4(ones(1, 1, 1, 1)), 1
        )

        @test_throws ArgumentError parallel_repetition(game, 2; max_entries=255)
        @test size(parallel_repetition(game, 2; max_entries=256)) == (4, 4, 4, 4)
        @test_throws ArgumentError parallel_repetition(game, 2; max_work=511)
        @test size(parallel_repetition(game, 2; max_work=512)) == (4, 4, 4, 4)
        @test_throws ArgumentError parallel_repetition(game, 2; max_entries=0)
        @test_throws ArgumentError parallel_repetition(game, 2; max_entries=true)
        @test_throws ArgumentError parallel_repetition(game, 2; max_entries=1.5)
        @test_throws ArgumentError parallel_repetition(game, 2; max_work=0)
        @test size(parallel_repetition(game, 2; max_entries=nothing, max_work=nothing)) ==
            (4, 4, 4, 4)

        @test_throws ArgumentError parallel_repetition(
            ones(Int, 2, 1, 1, 1), typemax(Int); max_entries=nothing, max_work=nothing
        )
        @test_throws ArgumentError parallel_repetition(
            ones(Int, 1, 1, 1, 1), typemax(Int); max_entries=nothing, max_work=1
        )
    end
end
