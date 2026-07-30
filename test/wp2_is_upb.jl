using LinearAlgebra
using QuantumEntanglementTools
using SparseArrays
using Test

const IsolatedUPB = QuantumEntanglementTools
const CompatUPB = QuantumEntanglementTools.MATLABCompat

struct _UPBZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end
Base.size(matrix::_UPBZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_UPBZeroBasedMatrix)
    return ntuple(index -> 0:(size(matrix.storage, index) - 1), 2)
end
Base.IndexStyle(::Type{<:_UPBZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_UPBZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function _wp2_tiles(::Type{T}=Int) where {T}
    left = T[
        1 1 0 0 1
        0 -1 0 1 1
        0 0 1 -1 1
    ]
    right = T[
        1 0 0 1 1
        -1 0 1 0 1
        0 1 -1 0 1
    ]
    return left, right
end

function _wp2_normalized_tiles(::Type{T}) where {T<:AbstractFloat}
    left, right = _wp2_tiles(T)
    for factor in (left, right), state in axes(factor, 2)
        factor[:, state] ./= norm(@view factor[:, state])
    end
    return left, right
end

function _wp2_global_columns(factors)
    return hcat(
        [
            foldl(kron, (factor[:, state] for factor in factors)) for
            state in axes(factors[1], 2)
        ]...,
    )
end

function _wp2_check_extension(result, factors)
    result.witness_factors === nothing && return false
    state_count = size(factors[1], 2)
    for state in 1:state_count
        overlap = prod(
            dot(result.witness_factors[party], @view(factors[party][:, state])) for
            party in eachindex(factors)
        )
        iszero(overlap) || return false
    end
    partition_indices = sort!(vcat((collect(part) for part in result.witness_partition)...))
    return partition_indices == collect(1:state_count)
end

@testset "WP2 IsUPB certificate-aware analysis" begin
    @testset "exact UPB certificate and private partition supersession" begin
        left, right = _wp2_tiles()
        result = IsolatedUPB.is_upb((left, right); normalization=:allow)
        @test result isa IsolatedUPB.UPBAnalysisResult
        @test result.status === :upb
        @test result.reason === :unextendible
        @test result.is_upb === true
        @test result.certificate_kind === :exact
        @test result.witness_factors === nothing
        @test result.witness_vector === nothing
        @test result.dimensions == (3, 3)
        @test result.state_count == 5
        @test result.input_form === :local_factors
        @test result.input_normalized === false
        @test !result.analysis_rescaled
        @test !result.densified
        @test result.partitions_examined == 20
        @test result.work_used == 2640
        @test result.max_states == 256
        @test iszero(result.orthogonality_residual)
        @test occursin("status=upb", sprint(show, result))
        @test occursin("partitions_examined=20", sprint(show, result))

        # Ordered set partitions with both part sizes at least two:
        # C(5,2) + C(5,3) = 20. This directly exercises the private,
        # lazy replacement for QETLAB's vec_partitions helper.
        @test result.partitions_examined == binomial(5, 2) + binomial(5, 3)
    end

    @testset "floating and generic-precision certificates" begin
        left32, right32 = _wp2_normalized_tiles(Float32)
        result32 = IsolatedUPB.is_upb(left32, right32)
        @test result32.status === :upb
        @test result32.certificate_kind === :tolerance_robust
        @test result32.atol isa Float32
        @test result32.rtol isa Float32
        @test result32.orthogonality_residual isa Float32
        @test result32.analysis_rescaled

        setprecision(BigFloat, 192) do
            left_big, right_big = _wp2_normalized_tiles(BigFloat)
            result_big = IsolatedUPB.is_upb((left_big, right_big))
            @test result_big.status === :upb
            @test result_big.rtol isa BigFloat
            @test result_big.orthogonality_residual isa BigFloat
            @test result_big.partitions_examined == 20
        end

        rational_left, rational_right = _wp2_tiles(Rational{Int})
        rational_result = IsolatedUPB.is_upb(
            rational_left, rational_right; normalization=:allow
        )
        @test rational_result.status === :upb
        @test rational_result.certificate_kind === :exact
    end

    @testset "explicit product-vector extension witnesses" begin
        left = [1 0; 0 1]
        right = [1 1; 0 0]
        result = IsolatedUPB.is_upb((left, right); materialize_witness=true)
        @test result.status === :not_upb
        @test result.reason === :extension_witness
        @test result.is_upb === false
        @test result.certificate_kind === :exact
        @test result.witness_factors isa Tuple
        @test result.witness_vector == foldl(kron, result.witness_factors)
        @test _wp2_check_extension(result, (left, right))
        @test iszero(result.witness_residual)
        @test result.partitions_examined == 1

        # Complex Hermitian orthogonality is deliberate. The pinned MATLAB
        # implementation uses a nonconjugating transpose in this branch.
        complex_left = reshape(Complex{Int}[1, im], 2, 1)
        complex_right = reshape(Complex{Int}[1, 0], 2, 1)
        complex_result = IsolatedUPB.is_upb(
            complex_left, complex_right; normalization=:allow, materialize_witness=true
        )
        @test complex_result.status === :not_upb
        @test _wp2_check_extension(complex_result, (complex_left, complex_right))
        @test eltype(complex_result.witness_factors[1]) == Complex{Rational{BigInt}}
        @test dot(complex_result.witness_factors[1], @view(complex_left[:, 1])) == 0

        empty_left = reshape(Int[], 2, 0)
        empty_right = reshape(Int[], 3, 0)
        empty_result = IsolatedUPB.is_upb(empty_left, empty_right; materialize_witness=true)
        @test empty_result.status === :not_upb
        @test empty_result.reason === :extension_witness
        @test length(empty_result.witness_vector) == 6
        @test empty_result.witness_partition == ((), ())
    end

    @testset "global-vector input and product validation" begin
        local_factors = _wp2_normalized_tiles(Float64)
        global_vectors = _wp2_global_columns(local_factors)
        global_result = IsolatedUPB.is_upb(global_vectors, (3, 3))
        @test global_result.status === :upb
        @test global_result.input_form === :global_vectors
        @test global_result.dimensions == (3, 3)
        @test global_result.offending_state === nothing
        @test global_result.product_residual < 1e-14
        @test global_result.witness_residual === nothing

        exact_factors = _wp2_tiles()
        exact_global = _wp2_global_columns(exact_factors)
        exact_global_result = IsolatedUPB.is_upb(exact_global, [3, 3]; normalization=:allow)
        @test exact_global_result.status === :upb
        @test exact_global_result.certificate_kind === :exact

        bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
        not_product = IsolatedUPB.is_upb(bell, (2, 2))
        @test not_product.status === :not_upb
        @test not_product.reason === :not_product
        @test not_product.offending_state == 1
        @test not_product.product_residual > 0
        @test not_product.witness_residual === nothing

        product = kron([1.0, 0.0], [1.0, 0.0])
        one_state = IsolatedUPB.is_upb(product, (2, 2))
        @test one_state.status === :not_upb
        @test one_state.reason === :extension_witness
    end

    @testset "definition checks beyond pinned unextendibility" begin
        generic_left = [
            1 1 1 1 1
            0 1 2 3 4
            0 1 4 9 16
        ]
        generic_right = [
            1 1 1 1 1
            0 1 3 9 27
            0 1 9 81 243
        ]
        nonorthogonal = IsolatedUPB.is_upb(
            generic_left, generic_right; normalization=:allow
        )
        @test nonorthogonal.status === :not_upb
        @test nonorthogonal.reason === :not_orthogonal
        @test nonorthogonal.offending_state == (1, 2)
        @test nonorthogonal.orthogonality_residual > 0

        basis_left = [
            1 1 0 0
            0 0 1 1
        ]
        basis_right = [
            1 0 1 0
            0 1 0 1
        ]
        complete = IsolatedUPB.is_upb(basis_left, basis_right)
        @test complete.status === :not_upb
        @test complete.reason === :complete_basis
        @test complete.is_upb === false
        @test complete.witness_factors === nothing
    end

    @testset "tri-state boundaries and deterministic limits" begin
        delta = 1.5e-4
        first = [1.0 0.0; 0.0 1.0]
        second = [
            1.0 delta
            0.0 sqrt(1 - delta^2)
        ]
        repeated = [1.0 1.0; 0.0 0.0]
        boundary = IsolatedUPB.is_upb(
            second, repeated; atol=1e-4, rtol=0, boundary_factor=2
        )
        @test boundary.status === :unknown
        @test boundary.reason === :numerical_boundary
        @test boundary.is_upb === nothing
        @test boundary.offending_state == (1, 2)

        outside = IsolatedUPB.is_upb(
            [
                1.0 3e-4
                0.0 sqrt(1 - (3e-4)^2)
            ],
            repeated;
            atol=1e-4,
            rtol=0,
            boundary_factor=2,
        )
        @test outside.status === :not_upb
        @test outside.reason === :not_orthogonal

        near_product = kron([1.0, 0.0], [1.0, 0.0]) + delta * kron([0.0, 1.0], [0.0, 1.0])
        near_product ./= norm(near_product)
        product_boundary = IsolatedUPB.is_upb(
            near_product, (2, 2); atol=1e-4, rtol=0, boundary_factor=2
        )
        @test product_boundary.status === :unknown
        @test product_boundary.reason === :numerical_boundary
        @test product_boundary.offending_state == 1

        tiles = _wp2_normalized_tiles(Float64)
        rank_boundary = IsolatedUPB.is_upb(tiles; atol=0.2, rtol=0, boundary_factor=8)
        @test rank_boundary.status === :unknown
        @test rank_boundary.reason === :numerical_boundary
        @test rank_boundary.partitions_examined == 20

        exact_tiles = _wp2_tiles()
        partition_limited = IsolatedUPB.is_upb(
            exact_tiles; normalization=:allow, max_partitions=3
        )
        @test partition_limited.status === :unknown
        @test partition_limited.reason === :partition_limit
        @test partition_limited.partitions_examined == 3
        @test partition_limited.is_upb === nothing

        work_limited = IsolatedUPB.is_upb(exact_tiles; normalization=:allow, max_work=0)
        @test work_limited.status === :unknown
        @test work_limited.reason === :work_limit
        @test work_limited.partitions_examined == 1
        @test work_limited.work_used == 0

        state_limited = IsolatedUPB.is_upb(exact_tiles; normalization=:allow, max_states=4)
        @test state_limited.status === :unknown
        @test state_limited.reason === :state_limit
        @test state_limited.partitions_examined == 0
        global_state_limited = IsolatedUPB.is_upb(
            _wp2_global_columns(exact_tiles), (3, 3); normalization=:allow, max_states=4
        )
        @test global_state_limited.status === :unknown
        @test global_state_limited.reason === :state_limit
    end

    @testset "normalization and sparse policy" begin
        scaled_left = reshape([2.0, 0.0], 2, 1)
        unit_right = reshape([1.0, 0.0], 2, 1)
        @test_throws ArgumentError IsolatedUPB.is_upb(scaled_left, unit_right)
        allowed = IsolatedUPB.is_upb(scaled_left, unit_right; normalization=:allow)
        @test allowed.status === :not_upb
        @test allowed.input_normalized === false
        @test allowed.analysis_rescaled

        left, right = _wp2_tiles()
        @test_throws ArgumentError IsolatedUPB.is_upb(
            sparse(left), sparse(right); normalization=:allow
        )
        sparse_result = IsolatedUPB.is_upb(
            sparse(left), sparse(right); normalization=:allow, allow_densify=true
        )
        @test sparse_result.status === :upb
        @test sparse_result.densified
        @test_throws ArgumentError IsolatedUPB.is_upb(
            sparse(left),
            sparse(right);
            normalization=:allow,
            allow_densify=true,
            max_dense_entries=29,
        )

        global_sparse = sparse(_wp2_global_columns((left, right)))
        @test_throws ArgumentError IsolatedUPB.is_upb(
            global_sparse, (3, 3); normalization=:allow
        )
        sparse_global_result = IsolatedUPB.is_upb(
            global_sparse, (3, 3); normalization=:allow, allow_densify=true
        )
        @test sparse_global_result.status === :upb
        @test sparse_global_result.densified

        empty_left = reshape(Int[], 2, 0)
        empty_right = reshape(Int[], 2, 0)
        @test_throws ArgumentError IsolatedUPB.is_upb(
            empty_left, empty_right; materialize_witness=true, max_dense_entries=3
        )
    end

    @testset "MATLAB compatibility mapping" begin
        tiles = _wp2_tiles()
        @test CompatUPB.IsUPB(tiles...) === true
        structured = CompatUPB.IsUPB(tiles...; structured=true)
        @test structured isa UPBAnalysisResult
        @test structured.status === :upb
        @test structured.normalization === :allow

        extendible_left = [1 0; 0 1]
        extendible_right = [1 1; 0 0]
        flag, witness = CompatUPB.IsUPB(
            extendible_left, extendible_right; return_witness=true
        )
        @test flag === false
        @test witness isa Tuple
        @test all(
            state -> iszero(
                prod(
                    dot(
                        witness[party],
                        (extendible_left, extendible_right)[party][:, state],
                    ) for party in 1:2
                ),
            ),
            1:2,
        )

        complete_left = [1 1 0 0; 0 0 1 1]
        complete_right = [1 0 1 0; 0 1 0 1]
        @test CompatUPB.IsUPB(complete_left, complete_right) === false
        @test_throws DomainError CompatUPB.IsUPB(tiles...; max_partitions=1)
        @test_throws ArgumentError CompatUPB.IsUPB(tiles...; return_witness=2)
        @test_throws ArgumentError CompatUPB.IsUPB(
            tiles...; return_witness=true, structured=true
        )
    end

    @testset "validation and overflow" begin
        left, right = _wp2_tiles()
        @test_throws ArgumentError IsolatedUPB.is_upb(())
        @test_throws ArgumentError IsolatedUPB.is_upb((left,))
        @test_throws ArgumentError IsolatedUPB.is_upb(Any[left, 1])
        @test_throws DimensionMismatch IsolatedUPB.is_upb(left, right[:, 1:4])
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(1, 2), zeros(2, 2))
        @test_throws ArgumentError IsolatedUPB.is_upb(
            Bool[true false; false true], Bool[true false; false true]
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            [1.0 NaN; 0.0 1.0], [1.0 0.0; 0.0 1.0]
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            zeros(2, 1), reshape([1.0, 0.0], 2, 1)
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            _UPBZeroBasedMatrix(left), right; normalization=:allow
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(left, right; normalization=:normalize)
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, max_partitions=-1
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, max_work=true
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, max_dense_entries=-1
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, max_states=-1
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, boundary_factor=1
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, atol=1 // 10
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, allow_densify=1
        )
        @test_throws ArgumentError IsolatedUPB.is_upb(
            left, right; normalization=:allow, materialize_witness=1
        )

        @test_throws DimensionMismatch IsolatedUPB.is_upb(zeros(5), (2, 3))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (4,))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (2, 1, 2))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (2.0, 2))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (true, 4))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (typemax(Int), 2))
        @test_throws ArgumentError IsolatedUPB.is_upb([1.0, 0.0, Inf, 0.0], (2, 2))
        @test_throws ArgumentError IsolatedUPB.is_upb(zeros(4), (2, 2))
    end
end
