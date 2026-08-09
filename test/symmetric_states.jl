using LinearAlgebra
using SparseArrays

function _test_collective_ambient(local_operator, parties)
    local_dimension = size(local_operator, 1)
    value_type = eltype(local_operator)
    total_dimension = local_dimension^parties
    result = zeros(value_type, total_dimension, total_dimension)
    for active in 1:parties
        factors = [
            if party == active
                local_operator
            else
                Matrix{value_type}(I, local_dimension, local_dimension)
            end for party in 1:parties
        ]
        result += reduce(kron, factors)
    end
    return result
end

@testset "Symmetric-state occupation coordinates" begin
    @testset "dimension, ordering, and rank" begin
        @test symmetric_subspace_dimension(3, 4) == big(15)
        @test symmetric_subspace_dimension(5, 0) == big(1)
        @test symmetric_subspace_dimension(1, 1_000) == big(1)
        @test symmetric_occupations(3, 2) ==
            [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
        @test symmetric_occupations(4, 0) == [(0, 0, 0, 0)]
        @test symmetric_basis_index((1, 0, 1)) == big(3)
        @test symmetric_basis_index((0, 1, 1)) == big(5)
        @test symmetric_basis_index((2, 0, 0, 0)) == big(1)
        huge_middle_rank = symmetric_basis_index((0, typemax(Int), 0))
        @test huge_middle_rank == 1 + binomial(BigInt(typemax(Int)) + 1, 2)
        @test huge_middle_rank > typemax(Int)
        @test symmetric_basis_occupation(5, 3, 2) == (0, 1, 1)
        @test symmetric_basis_occupation(5, 3, 2; max_entries=3, max_work=100) == (0, 1, 1)

        large_vacuum = only(
            symmetric_occupations(1_000, 0; max_occupations=1, max_work=1_000_000)
        )
        @test length(large_vacuum) == 1_000
        @test all(iszero, large_vacuum)
        @test_throws ArgumentError symmetric_occupations(
            100_000, 0; max_occupations=1, max_entries=100_000, max_work=1_000_000
        )

        for local_dimension in 1:4, parties in 0:5
            occupations = symmetric_occupations(local_dimension, parties)
            @test length(occupations) ==
                symmetric_subspace_dimension(local_dimension, parties)
            for (index, occupation) in pairs(occupations)
                @test symmetric_basis_index(occupation) == index
                @test symmetric_basis_occupation(index, local_dimension, parties) ==
                    occupation
            end
        end

        @test_throws ArgumentError symmetric_subspace_dimension(true, 2)
        @test_throws ArgumentError symmetric_subspace_dimension(2, -1)
        @test_throws ArgumentError symmetric_occupations(0, 2)
        @test_throws ArgumentError symmetric_occupations(3, 3; max_occupations=9)
        @test_throws ArgumentError symmetric_occupations(446, 2)
        @test_throws ArgumentError symmetric_occupations(10, 2; max_entries=549)
        @test_throws ArgumentError symmetric_occupations(3, 3; max_work=1)
        @test_throws ArgumentError symmetric_occupations(3, 3; max_work=0)
        @test_throws ArgumentError symmetric_basis_index(())
        @test_throws ArgumentError symmetric_basis_index((1, -1))
        @test_throws ArgumentError symmetric_basis_index((1, true))
        @test_throws ArgumentError symmetric_basis_index((1, 0.5))
        @test_throws ArgumentError symmetric_basis_occupation(false, 2, 2)
        @test_throws ArgumentError symmetric_basis_occupation(0, 2, 2)
        @test_throws BoundsError symmetric_basis_occupation(4, 2, 2)
        @test_throws ArgumentError symmetric_basis_occupation(5, 3, 2; max_entries=2)
        @test_throws ArgumentError symmetric_basis_occupation(5, 3, 2; max_work=1)
        @test_throws ArgumentError symmetric_basis_occupation(1, typemax(Int), 0)
        @test_throws ArgumentError symmetric_occupations(typemax(Int), 0)
    end

    @testset "generalized Dicke states" begin
        for parties in 1:6, excitations in 0:parties
            generalized = generalized_dicke_state((parties - excitations, excitations))
            @test generalized == dicke_state(parties, excitations)
        end

        state = generalized_dicke_state((1, 2, 1))
        @test length(state) == 3^4
        @test nnz(state) == 12
        @test isapprox(norm(state), 1)
        @test swap_subsystems(state, ntuple(_ -> 3, 4), 1, 4) == state
        @test swap_subsystems(state, ntuple(_ -> 3, 4), 2, 3) == state

        for local_dimension in 1:4, parties in 1:4
            basis = symmetric_subspace_basis(local_dimension, parties)
            for (index, occupation) in
                pairs(symmetric_occupations(local_dimension, parties))
                @test generalized_dicke_state(occupation) == basis[:, index]
            end
        end

        raw = generalized_dicke_state((1, 1, 0); normalized=false, T=Int)
        @test eltype(raw) == Int
        @test nonzeros(raw) == [1, 1]
        @test generalized_dicke_state((0, 0, 0)) == sparsevec([1], [1.0], 1)
        @test generalized_dicke_state((2, 1); T=Float32) isa SparseVector{Float32}
        @test generalized_dicke_state((2, 1); T=BigFloat) isa SparseVector{BigFloat}
        @test generalized_dicke_state((2, 1); T=ComplexF32) isa SparseVector{ComplexF32}
        @test generalized_dicke_state((1, 1); sparse_output=false) isa Vector

        large_single_level = generalized_dicke_state(
            (100_000,); max_nonzeros=1, max_work=1_000_000
        )
        @test large_single_level == sparsevec([1], [1.0], 1)

        @test_throws ArgumentError generalized_dicke_state([])
        @test_throws ArgumentError generalized_dicke_state((1, -1))
        @test_throws ArgumentError generalized_dicke_state((2, 2); max_nonzeros=5)
        @test_throws ArgumentError generalized_dicke_state(
            (10, 0); sparse_output=false, max_dense_entries=1_000
        )
        @test_throws ArgumentError generalized_dicke_state((2, 2); max_work=5)
        @test_throws ArgumentError generalized_dicke_state((typemax(Int),))
    end

    @testset "product-state coordinates" begin
        local_state = ComplexF64[1, im, 2] / sqrt(6)
        for parties in 0:4
            coordinates = symmetric_product_coordinates(local_state, parties)
            @test length(coordinates) == symmetric_subspace_dimension(3, parties)
            @test isapprox(norm(coordinates), norm(local_state)^parties; atol=1e-13)
            if parties > 0
                basis = symmetric_subspace_basis(3, parties; T=ComplexF64)
                ambient = reduce(kron, ntuple(_ -> local_state, parties))
                @test isapprox(coordinates, adjoint(basis) * ambient; atol=1e-13)
            else
                @test coordinates == ComplexF64[1]
            end
        end

        sparse_local = sparsevec([1, 3], [1.0, 2.0], 3)
        sparse_coordinates = symmetric_product_coordinates(sparse_local, 3)
        @test issparse(sparse_coordinates)
        @test Vector(sparse_coordinates) ==
            symmetric_product_coordinates(Vector(sparse_local), 3)
        @test symmetric_product_coordinates(Float32[1, 1], 3) isa Vector{Float32}
        @test symmetric_product_coordinates(BigFloat[1, 1], 3) isa Vector{BigFloat}
        @test symmetric_product_coordinates(ComplexF32[1, im], 2) isa Vector{ComplexF32}

        for exact_local_state in (Int[2, -3, 5], Rational{Int}[1 // 2, -2 // 3, 7 // 5])
            vacuum = symmetric_product_coordinates(exact_local_state, 0)
            one_particle = symmetric_product_coordinates(exact_local_state, 1)
            @test eltype(vacuum) == eltype(exact_local_state)
            @test vacuum == [one(eltype(exact_local_state))]
            @test eltype(one_particle) == eltype(exact_local_state)
            @test one_particle == exact_local_state
        end

        abstract_sparse_local = sparsevec([1, 2, 3], Number[1, 2.0, 0], 3)
        abstract_sparse_product = symmetric_product_coordinates(abstract_sparse_local, 1)
        @test abstract_sparse_product isa SparseVector{Float64}
        @test nnz(abstract_sparse_product) == 2
        @test Vector(abstract_sparse_product) == [1.0, 2.0, 0.0]

        high_precision_local, high_precision_coordinates = setprecision(BigFloat, 1_024) do
            local_state = BigFloat[BigFloat(1) / 3, sqrt(BigFloat(2)), -BigFloat(7) / 11]
            return local_state, symmetric_product_coordinates(local_state, 1)
        end
        @test high_precision_coordinates == high_precision_local
        @test all(value -> precision(value) == 1_024, high_precision_coordinates)

        large = symmetric_product_coordinates(
            BigFloat[inv(sqrt(big(2))), inv(sqrt(big(2)))], 200
        )
        @test all(isfinite, large)
        @test isapprox(norm(large), big(1); rtol=big"1e-60")

        @test_throws ArgumentError symmetric_product_coordinates(Float64[], 2)
        @test_throws ArgumentError symmetric_product_coordinates([1.0, Inf], 2)
        @test_throws ArgumentError symmetric_product_coordinates(
            [1.0, 1.0], 10; max_occupations=10
        )
        @test_throws ArgumentError symmetric_product_coordinates([1.0, 1.0], 3; max_work=5)
        @test_throws OverflowError symmetric_product_coordinates(
            Float16[nextfloat(Float16(0)), 1], 2
        )

        exact_large = (BigInt(1) << 400) + 1
        mixed_precision_local = setprecision(BigFloat, 1_024) do
            Real[BigFloat(1) / 3, exact_large]
        end
        mixed_precision_product = symmetric_product_coordinates(
            mixed_precision_local, 1; sparse_output=false
        )
        exact_large_reference = setprecision(BigFloat, 1_024) do
            BigFloat(exact_large)
        end
        @test mixed_precision_product[2] == exact_large_reference
        @test precision(mixed_precision_product[1]) >= 1_024
    end

    @testset "collective operators" begin
        local_operator = ComplexF64[1 im 0; -im 2 1; 0 1 -1]
        for parties in 1:3
            basis = symmetric_subspace_basis(3, parties; T=ComplexF64)
            ambient = _test_collective_ambient(local_operator, parties)
            expected = adjoint(basis) * ambient * basis
            actual = symmetric_collective_operator(local_operator, parties)
            @test issparse(actual)
            @test isapprox(Matrix(actual), expected; atol=1e-12)
            @test ishermitian(actual)
        end

        nonhermitian_local = ComplexF64[1 2+im; 0 -im]
        nonhermitian_basis = symmetric_subspace_basis(2, 3; T=ComplexF64)
        nonhermitian_expected =
            adjoint(nonhermitian_basis) *
            _test_collective_ambient(nonhermitian_local, 3) *
            nonhermitian_basis
        nonhermitian_actual = symmetric_collective_operator(nonhermitian_local, 3)
        @test isapprox(Matrix(nonhermitian_actual), nonhermitian_expected; atol=1e-12)
        @test !ishermitian(nonhermitian_actual)

        integer_diagonal_collective = symmetric_collective_operator([2 0; 0 -3], 3)
        @test eltype(integer_diagonal_collective) == Int
        @test Matrix(integer_diagonal_collective) == [
            6 0 0 0
            0 1 0 0
            0 0 -4 0
            0 0 0 -9
        ]

        rational_diagonal = Rational{Int}[1//2 0; 0 -1//3]
        rational_diagonal_collective = symmetric_collective_operator(rational_diagonal, 2)
        @test eltype(rational_diagonal_collective) == Rational{Int}
        @test Matrix(rational_diagonal_collective) == Rational{Int}[
            1 0 0
            0 1//6 0
            0 0 -2//3
        ]

        for exact_local_operator in (Int[1 2; -3 4], Rational{Int}[1//2 2//3; -1//4 3//5])
            zero_particle_exact = symmetric_collective_operator(exact_local_operator, 0)
            one_particle_exact = symmetric_collective_operator(exact_local_operator, 1)
            @test eltype(zero_particle_exact) == eltype(exact_local_operator)
            @test size(zero_particle_exact) == (1, 1)
            @test iszero(zero_particle_exact)
            @test eltype(one_particle_exact) == eltype(exact_local_operator)
            @test Matrix(one_particle_exact) == exact_local_operator
        end

        zero_particle = symmetric_collective_operator([1 2; 3 4], 0)
        @test size(zero_particle) == (1, 1)
        @test iszero(zero_particle)
        @test Matrix(symmetric_collective_operator([1 0; 0 -1], 2)) == [
            2.0 0.0 0.0
            0.0 0.0 0.0
            0.0 0.0 -2.0
        ]
        @test symmetric_collective_operator(sparse([1.0 0.0; 0.0 -1.0]), 2) isa
            SparseMatrixCSC
        sparse_single_transition = sparse([1], [2], [3.0], 100, 100)
        @test symmetric_collective_operator(
            sparse_single_transition, 1; max_nonzeros=1, max_work=30_000
        ) == sparse_single_transition
        wrapped_diagonal = Hermitian(spdiagm(0 => [1.0, 2.0]))
        @test Matrix(symmetric_collective_operator(wrapped_diagonal, 2)) ==
            Matrix(symmetric_collective_operator(Matrix(wrapped_diagonal), 2))
        @test symmetric_collective_operator(adjoint(sparse_single_transition), 1) ==
            adjoint(sparse_single_transition)
        @test symmetric_collective_operator(Float32[1 0; 0 -1], 2; sparse_output=false) isa
            Matrix{Float32}

        float16_collective = symmetric_collective_operator(Float16[0 1; 1 0], 1_000)
        @test eltype(float16_collective) == Float16
        @test all(isfinite, nonzeros(float16_collective))
        @test maximum(abs, nonzeros(float16_collective)) > Float16(400)

        sx = ComplexF64[0 1; 1 0] / 2
        sy = ComplexF64[0 -im; im 0] / 2
        sz = ComplexF64[1 0; 0 -1] / 2
        jx = symmetric_collective_operator(sx, 5)
        jy = symmetric_collective_operator(sy, 5)
        jz = symmetric_collective_operator(sz, 5)
        @test isapprox(Matrix(jx * jy - jy * jx), Matrix(im * jz); atol=1e-12)

        @test_throws DimensionMismatch symmetric_collective_operator(ones(2, 3), 2)
        @test_throws ArgumentError symmetric_collective_operator([1.0 NaN; 0 1], 2)
        @test_throws ArgumentError symmetric_collective_operator(
            ones(2, 2), 3; max_nonzeros=1
        )
        @test_throws ArgumentError symmetric_collective_operator(
            ones(2, 2), 3; sparse_output=false, max_dense_entries=10
        )
        @test_throws ArgumentError symmetric_collective_operator(ones(2, 2), 3; max_work=1)

        tiny64 = nextfloat(0.0)
        cancellation_operator = Matrix(Diagonal([1.0, tiny64, -1.0]))
        cancellation_collective = symmetric_collective_operator(cancellation_operator, 3)
        cancellation_index = Int(symmetric_basis_index((1, 1, 1)))
        @test cancellation_collective[cancellation_index, cancellation_index] == tiny64

        exact_scale = Int(2)^54
        exact_cancellation_operator = Int[exact_scale 1; 0 -exact_scale+1]
        exact_cancellation_collective = symmetric_collective_operator(
            exact_cancellation_operator, 2
        )
        @test exact_cancellation_collective[2, 2] == 1

        hidden_nan_parent = sparse([1, 2], [1, 3], [1.0, NaN], 2, 3)
        finite_sparse_view = @view hidden_nan_parent[:, 1:2]
        @test symmetric_collective_operator(finite_sparse_view, 1) ==
            sparse([1], [1], [1.0], 2, 2)
    end

    @testset "split isometries and reduced states" begin
        for local_dimension in 1:4, parties in 0:4, keep in 0:parties
            isometry = symmetric_split_isometry(local_dimension, parties, keep)
            global_dimension = Int(symmetric_subspace_dimension(local_dimension, parties))
            kept_dimension = Int(symmetric_subspace_dimension(local_dimension, keep))
            traced_dimension = Int(
                symmetric_subspace_dimension(local_dimension, parties - keep)
            )
            @test size(isometry) == (kept_dimension * traced_dimension, global_dimension)
            @test nnz(isometry) == kept_dimension * traced_dimension
            @test isapprox(
                Matrix(adjoint(isometry) * isometry),
                Matrix{Float64}(I, global_dimension, global_dimension);
                atol=1e-13,
            )
        end

        dense_isometry = symmetric_split_isometry(2, 3, 1; sparse_output=false)
        @test dense_isometry isa Matrix
        @test symmetric_split_isometry(2, 3, 1; T=Float32) isa SparseMatrixCSC{Float32}
        @test_throws ArgumentError symmetric_split_isometry(2, 3, 4)
        @test_throws ArgumentError symmetric_split_isometry(3, 4, 2; max_nonzeros=10)
        @test_throws ArgumentError symmetric_split_isometry(
            3, 4, 2; sparse_output=false, max_dense_entries=10
        )
        @test_throws ArgumentError symmetric_split_isometry(3, 4, 2; max_work=10)
        @test_throws ArgumentError symmetric_split_isometry(typemax(Int), 0, 0)
        @test size(symmetric_split_isometry(2, 99, 0; max_nonzeros=100, max_work=1_000)) ==
            (100, 100)

        split = symmetric_split_isometry(3, 3, 1)
        global_basis = symmetric_subspace_basis(3, 3)
        kept_basis = symmetric_subspace_basis(3, 1)
        traced_basis = symmetric_subspace_basis(3, 2)
        @test isapprox(kron(kept_basis, traced_basis) * split, global_basis; atol=2e-13)

        local_state = ComplexF64[1, 1 + im, -im]
        local_state /= norm(local_state)
        parties = 3
        coordinates = symmetric_product_coordinates(local_state, parties)
        global_basis = symmetric_subspace_basis(3, parties; T=ComplexF64)
        ambient_state = global_basis * coordinates
        dimensions = ntuple(_ -> 3, parties)
        for keep in 0:parties
            reduced = symmetric_reduced_state(
                coordinates, 3, parties; keep=keep, sparse_output=false
            )
            traced_systems = Tuple((keep + 1):parties)
            ambient_reduced = partial_trace(
                ambient_state, dimensions; trace_out=traced_systems
            )
            expected = if keep == 0
                ambient_reduced
            else
                kept_basis = symmetric_subspace_basis(3, keep; T=ComplexF64)
                adjoint(kept_basis) * ambient_reduced * kept_basis
            end
            @test isapprox(reduced, expected; atol=2e-12)
            @test ishermitian(reduced)
            @test isapprox(tr(reduced), 1; atol=2e-12)

            mixed_reduced = symmetric_reduced_state(
                coordinates * adjoint(coordinates),
                3,
                parties;
                keep=keep,
                sparse_output=false,
            )
            @test isapprox(mixed_reduced, reduced; atol=2e-12)
        end

        dicke_coordinates = sparsevec([2], [1.0], 4)
        one_body = Matrix(symmetric_reduced_state(dicke_coordinates, 2, 3; keep=1))
        @test isapprox(one_body, [2 / 3 0; 0 1 / 3]; atol=1e-14)
        @test issparse(symmetric_reduced_state(dicke_coordinates, 2, 3; keep=1))
        @test symmetric_reduced_state(
            dicke_coordinates, 2, 3; keep=0, sparse_output=false
        ) == ones(1, 1)
        @test symmetric_reduced_state(
            dicke_coordinates * adjoint(dicke_coordinates),
            2,
            3;
            keep=3,
            sparse_output=false,
        ) == Matrix(dicke_coordinates * adjoint(dicke_coordinates))

        rational_coordinates = Rational{Int}[1 // 2, 1 // 3, 0]
        rational_zero_body = symmetric_reduced_state(
            rational_coordinates, 2, 2; keep=0, sparse_output=false
        )
        rational_full_body = symmetric_reduced_state(
            rational_coordinates, 2, 2; keep=2, sparse_output=false
        )
        @test rational_zero_body == reshape(Rational{Int}[13 // 36], 1, 1)
        @test eltype(rational_zero_body) == Rational{Int}
        @test rational_full_body == rational_coordinates * adjoint(rational_coordinates)
        @test eltype(rational_full_body) == Rational{Int}

        rational_operator = Rational{Int}[
            1//2 1//3 0
            0 1//4 1//5
            0 0 1//6
        ]
        @test symmetric_reduced_state(
            rational_operator, 2, 2; keep=0, sparse_output=false
        ) == reshape(Rational{Int}[11 // 12], 1, 1)
        @test symmetric_reduced_state(
            rational_operator, 2, 2; keep=2, sparse_output=false
        ) == rational_operator

        nonhermitian_operator = reshape(
            ComplexF64.(1:16) .+ im .* ComplexF64.(16:-1:1), 4, 4
        )
        nonhermitian_split = symmetric_split_isometry(2, 3, 1)
        lifted_operator =
            nonhermitian_split * nonhermitian_operator * adjoint(nonhermitian_split)
        expected_nonhermitian_reduction = partial_trace(
            lifted_operator, (2, 3); trace_out=(2,)
        )
        actual_nonhermitian_reduction = symmetric_reduced_state(
            nonhermitian_operator, 2, 3; keep=1, sparse_output=false
        )
        @test isapprox(
            actual_nonhermitian_reduction, expected_nonhermitian_reduction; atol=2e-13
        )
        @test !ishermitian(actual_nonhermitian_reduction)

        sparse_operator = sparse(
            [1, 1, 2, 3, 4], [1, 2, 3, 2, 4], ComplexF64[1, 2 + im, -im, 3, 4 - im], 4, 4
        )
        sparse_split = symmetric_split_isometry(2, 3, 1)
        expected_sparse_reduction = partial_trace(
            sparse_split * sparse_operator * adjoint(sparse_split), (2, 3); trace_out=(2,)
        )
        actual_sparse_reduction = symmetric_reduced_state(sparse_operator, 2, 3; keep=1)
        @test issparse(actual_sparse_reduction)
        @test isapprox(
            Matrix(actual_sparse_reduction), Matrix(expected_sparse_reduction); atol=2e-13
        )
        @test !iszero(actual_sparse_reduction[1, 2])

        large_rational = Rational{Int}(typemax(Int), 1)
        overflowing_rational_vector = Rational{Int}[large_rational, large_rational]
        widened_rational_norm = symmetric_reduced_state(
            overflowing_rational_vector, 2, 1; keep=0, sparse_output=false
        )
        @test eltype(widened_rational_norm) == Rational{BigInt}
        @test only(widened_rational_norm) == (2 * BigInt(typemax(Int))^2) // BigInt(1)
        overflowing_rational_matrix = Rational{Int}[
            large_rational 0
            0 large_rational
        ]
        widened_rational_trace = symmetric_reduced_state(
            overflowing_rational_matrix, 2, 1; keep=0, sparse_output=false
        )
        @test only(widened_rational_trace) == (2 * BigInt(typemax(Int))) // BigInt(1)

        for float_type in (Float16, Float32, Float64)
            tiny = nextfloat(zero(float_type))
            tiny_coordinates = float_type[tiny, 0]
            @test_throws OverflowError symmetric_reduced_state(
                tiny_coordinates, 2, 1; keep=0, sparse_output=false
            )
            @test_throws OverflowError symmetric_reduced_state(
                tiny_coordinates, 2, 1; keep=1, sparse_output=false
            )
        end
        tiny16 = nextfloat(Float16(0))
        tiny16_operator = spdiagm(0 => Float16[0, tiny16, 0])
        @test_throws OverflowError symmetric_reduced_state(tiny16_operator, 2, 2; keep=1)
        @test_throws OverflowError symmetric_reduced_state(
            tiny16_operator, 2, 2; keep=1, sparse_output=false
        )

        abstract_coordinates = Number[1, 2.0]
        @test symmetric_reduced_state(
            abstract_coordinates, 2, 1; keep=1, sparse_output=false
        ) isa Matrix{Float64}
        abstract_sparse_coordinates = sparsevec([1, 2], Number[1, 2.0], 2)
        abstract_sparse_reduction = symmetric_reduced_state(
            abstract_sparse_coordinates, 2, 1; keep=1
        )
        @test abstract_sparse_reduction isa SparseMatrixCSC{Float64}
        @test Matrix(abstract_sparse_reduction) == [1.0 2.0; 2.0 4.0]
        abstract_operator = Number[1 2.0; 3 4]
        @test symmetric_reduced_state(
            abstract_operator, 2, 1; keep=0, sparse_output=false
        ) == reshape([5.0], 1, 1)
        @test symmetric_reduced_state(
            abstract_operator, 2, 1; keep=1, sparse_output=false
        ) isa Matrix{Float64}

        wrapped_symmetric_operator = Hermitian(spdiagm(0 => ones(3)))
        @test symmetric_reduced_state(wrapped_symmetric_operator, 2, 2; keep=2) ==
            sparse(wrapped_symmetric_operator)
        @test symmetric_reduced_state(wrapped_symmetric_operator, 2, 2; keep=1) ==
            symmetric_reduced_state(sparse(wrapped_symmetric_operator), 2, 2; keep=1)

        huge_exact = BigInt(1) << 400
        exact_cancellation_state = Matrix(
            Diagonal(BigInt[huge_exact, -2huge_exact + 1, huge_exact])
        )
        exact_cancellation_reduction = symmetric_reduced_state(
            exact_cancellation_state, 2, 2; keep=1, sparse_output=false
        )
        expected_exact_cancellation = Matrix(Diagonal(fill(BigFloat(1) / 2, 2)))
        @test isapprox(
            exact_cancellation_reduction,
            expected_exact_cancellation;
            rtol=big"1e-70",
            atol=big"1e-70",
        )

        extreme_trace_state = Matrix(
            Diagonal([floatmax(Float64), nextfloat(0.0), -floatmax(Float64)])
        )
        @test only(
            symmetric_reduced_state(extreme_trace_state, 2, 2; keep=0, sparse_output=false)
        ) == nextfloat(0.0)

        zero_large_split = spzeros(Float64, 101)
        zero_large_reduction = symmetric_reduced_state(
            zero_large_split, 2, 100; keep=50, max_nonzeros=1, max_work=1_000_000
        )
        @test size(zero_large_reduction) == (51, 51)
        @test nnz(zero_large_reduction) == 0

        active_extreme_sector = sparsevec([1], [Float16(1)], 55)
        active_extreme_reduction = symmetric_reduced_state(
            active_extreme_sector, 2, 54; keep=27
        )
        @test active_extreme_reduction[1, 1] == Float16(1)
        @test nnz(active_extreme_reduction) == 1

        @test_throws DimensionMismatch symmetric_reduced_state([1.0, 0], 2, 2; keep=1)
        @test_throws DimensionMismatch symmetric_reduced_state(ones(2, 2), 2, 2; keep=1)
        @test_throws ArgumentError symmetric_reduced_state([1.0, NaN, 0], 2, 2; keep=1)
        @test_throws ArgumentError symmetric_reduced_state(
            ones(4), 2, 3; keep=1, max_nonzeros=1
        )
        @test_throws ArgumentError symmetric_reduced_state(
            ones(4), 2, 3; keep=2, sparse_output=false, max_dense_entries=2
        )
        @test_throws ArgumentError symmetric_reduced_state(
            ones(4), 2, 3; keep=2, max_work=1
        )
    end

    @testset "maximally mixed symmetric states" begin
        coordinates = symmetric_maximally_mixed_state(3, 2)
        @test coordinates == spdiagm(0 => fill(big(1) // big(6), 6))
        @test tr(coordinates) == 1
        @test symmetric_maximally_mixed_state(3, 2; sparse_output=false) isa
            Matrix{Rational{BigInt}}

        ambient = symmetric_maximally_mixed_state(
            2, 3; representation=:ambient, T=Rational{BigInt}
        )
        expected = symmetric_projector(2, 3; T=Rational{BigInt}) / 4
        @test ambient == expected
        @test tr(ambient) == 1
        @test symmetric_maximally_mixed_state(4, 0) == sparse([1], [1], [1 // 1])

        float_coordinates = symmetric_maximally_mixed_state(2, 3; T=Float32)
        @test float_coordinates isa SparseMatrixCSC{Float32}

        float16_coordinates = symmetric_maximally_mixed_state(
            2,
            70_000;
            T=Float16,
            max_occupations=70_001,
            max_nonzeros=70_001,
            max_work=1_000_000,
        )
        @test float16_coordinates isa SparseMatrixCSC{Float16}
        @test all(isfinite, nonzeros(float16_coordinates))
        @test all(!iszero, nonzeros(float16_coordinates))
        @test isapprox(sum(Float64, nonzeros(float16_coordinates)), 1.0; atol=2e-3, rtol=0)
        exact_complex_coordinates = symmetric_maximally_mixed_state(
            2, 2; T=Complex{Rational{BigInt}}
        )
        @test exact_complex_coordinates[1, 1] == (big(1) // big(3)) + 0im

        dense_ambient = symmetric_maximally_mixed_state(
            2, 3; representation=:ambient, sparse_output=false, max_dense_entries=64
        )
        @test size(dense_ambient) == (8, 8)
        @test_throws ArgumentError symmetric_maximally_mixed_state(
            2, 3; representation=:ambient, sparse_output=false, max_dense_entries=63
        )
        @test_throws ArgumentError symmetric_maximally_mixed_state(
            2, 2; representation=:full
        )
        @test_throws ArgumentError symmetric_maximally_mixed_state(3, 3; max_occupations=9)
        @test_throws ArgumentError symmetric_maximally_mixed_state(3, 3; max_work=1)
        @test_throws ArgumentError symmetric_maximally_mixed_state(
            3, 3; sparse_output=false, max_dense_entries=50
        )
    end

    @testset "overflow and retained precision regressions" begin
        @test_throws OverflowError symmetric_collective_operator([typemax(Int) 0; 0 0], 2)

        complex_integer_diagonal = Complex{Int}[
            complex(typemax(Int), typemax(Int)) 0
            0 0
        ]
        widened_collective = symmetric_collective_operator(complex_integer_diagonal, 2)
        twice_typemax = 2 * BigInt(typemax(Int))
        @test eltype(widened_collective) == Complex{BigInt}
        @test real(widened_collective[1, 1]) == twice_typemax
        @test imag(widened_collective[1, 1]) == twice_typemax

        integer_coordinates = Int[typemax(Int), -2, 3]
        wide_coordinates = BigInt.(integer_coordinates)
        expected_norm_squared = sum(abs2, wide_coordinates)
        integer_zero_body = symmetric_reduced_state(
            integer_coordinates, 2, 2; keep=0, sparse_output=false
        )
        integer_full_body = symmetric_reduced_state(
            integer_coordinates, 2, 2; keep=2, sparse_output=false
        )
        @test eltype(integer_zero_body) == BigInt
        @test integer_zero_body == reshape(BigInt[expected_norm_squared], 1, 1)
        @test eltype(integer_full_body) == BigInt
        @test integer_full_body == wide_coordinates * adjoint(wide_coordinates)

        large_float16_coordinates = fill(floatmax(Float16), 4)
        for keep in (0, 1, 3)
            @test_throws OverflowError symmetric_reduced_state(
                large_float16_coordinates, 2, 3; keep=keep, sparse_output=false
            )
        end

        cheap_sparse_coordinates = sparsevec([51], [2.0], 101)
        cheap_zero_body = symmetric_reduced_state(
            cheap_sparse_coordinates, 2, 100; keep=0, max_nonzeros=1
        )
        @test cheap_zero_body == sparse([1], [1], [4.0], 1, 1)

        high_precision_local_operator = setprecision(BigFloat, 1_024) do
            return BigFloat[
                BigFloat(1) / 3 sqrt(BigFloat(2))
                BigFloat(5) / 7 -BigFloat(11) / 13
            ]
        end
        high_precision_collective = symmetric_collective_operator(
            high_precision_local_operator, 3
        )
        @test eltype(high_precision_collective) == BigFloat
        @test all(value -> precision(value) >= 1_024, nonzeros(high_precision_collective))

        high_precision_coordinates = setprecision(BigFloat, 1_024) do
            return BigFloat[
                BigFloat(1) / 3, sqrt(BigFloat(2)) / 5, BigFloat(7) / 11, BigFloat(13) / 17
            ]
        end
        high_precision_reduction = symmetric_reduced_state(
            high_precision_coordinates, 2, 3; keep=1, sparse_output=false
        )
        @test eltype(high_precision_reduction) == BigFloat
        @test all(value -> precision(value) >= 1_024, high_precision_reduction)
    end
end
