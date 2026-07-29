using LinearAlgebra
using Random
using SparseArrays

struct _ZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_ZeroBasedVector) = size(vector.storage)
Base.axes(vector::_ZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_ZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_ZeroBasedVector, index::Int) = vector.storage[index + 1]

struct _ZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_ZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_ZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_ZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_ZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

function _direct_basis(index::Int, dims::Tuple)
    zero_based = index - 1
    return ntuple(length(dims)) do system
        stride = prod(dims[(system + 1):end]; init=1)
        return mod(div(zero_based, stride), dims[system]) + 1
    end
end

function _direct_linear(indices, dims::Tuple)
    zero_based = 0
    for system in eachindex(dims)
        zero_based = zero_based * dims[system] + indices[system] - 1
    end
    return zero_based + 1
end

function _direct_partial_trace(matrix, dims::Tuple, trace_out::Tuple)
    keep = Tuple(system for system in eachindex(dims) if !(system in trace_out))
    kept_dims = Tuple(dims[system] for system in keep)
    traced_dims = Tuple(dims[system] for system in trace_out)
    kept_dimension = prod(kept_dims; init=1)
    traced_dimension = prod(traced_dims; init=1)
    result = zeros(eltype(matrix), kept_dimension, kept_dimension)
    for kept_column in 1:kept_dimension, kept_row in 1:kept_dimension
        row_keep_basis = _direct_basis(kept_row, kept_dims)
        column_keep_basis = _direct_basis(kept_column, kept_dims)
        for traced in 1:traced_dimension
            traced_basis = _direct_basis(traced, traced_dims)
            row_basis = Vector{Int}(undef, length(dims))
            column_basis = Vector{Int}(undef, length(dims))
            for (position, system) in pairs(keep)
                row_basis[system] = row_keep_basis[position]
                column_basis[system] = column_keep_basis[position]
            end
            for (position, system) in pairs(trace_out)
                row_basis[system] = traced_basis[position]
                column_basis[system] = traced_basis[position]
            end
            result[kept_row, kept_column] += matrix[
                _direct_linear(row_basis, dims), _direct_linear(column_basis, dims)
            ]
        end
    end
    return result
end

function _direct_partial_transpose(matrix, dims::Tuple, systems::Tuple)
    dimension = prod(dims)
    result = similar(matrix)
    for old_column in 1:dimension, old_row in 1:dimension
        new_row_basis = collect(_direct_basis(old_row, dims))
        new_column_basis = collect(_direct_basis(old_column, dims))
        for system in systems
            new_row_basis[system], new_column_basis[system] = new_column_basis[system],
            new_row_basis[system]
        end
        result[_direct_linear(new_row_basis, dims), _direct_linear(new_column_basis, dims)] = matrix[
            old_row, old_column
        ]
    end
    return result
end

@testset "Tier A subsystem kernel" begin
    @testset "layouts and basis indices" begin
        layout = SubsystemLayout((2, 3, 4))
        @test layout.dims == (2, 3, 4)
        @test layout.strides == (12, 4, 1)
        @test layout.total_dimension == 24
        @test Tuple(SubsystemLayout([2, 3, 4])) == (2, 3, 4)
        @test length(layout) == 3
        @test collect(layout) == [2, 3, 4]
        @test SubsystemLayout(()).total_dimension == 1

        @test basis_to_linear((1, 1), (2, 3)) == 1
        @test basis_to_linear((2, 1), (2, 3)) == 4
        @test basis_to_linear((2, 3), (2, 3)) == 6
        @test linear_to_basis(4, (2, 3)) == (2, 1)
        for index in 1:layout.total_dimension
            @test basis_to_linear(linear_to_basis(index, layout), layout) == index
        end

        @test_throws ArgumentError SubsystemLayout((2, 0))
        @test_throws ArgumentError SubsystemLayout((2, true))
        @test_throws ArgumentError SubsystemLayout((2, 1.5))
        @test_throws ArgumentError SubsystemLayout((typemax(Int), 2))
        @test_throws MethodError SubsystemLayout((2, 2), (2, 1), 4)
        @test_throws MethodError SubsystemLayout(Val(:validated), (2, 2), (2, 1), 4)
        forged_token = QuantumEntanglementTools._ValidatedConstructorToken()
        @test_throws ArgumentError SubsystemLayout(forged_token, (2, 2), (2, 1), 4)
        @test_throws DimensionMismatch basis_to_linear((1,), (2, 3))
        @test_throws ArgumentError basis_to_linear((3, 1), (2, 3))
        @test_throws ArgumentError linear_to_basis(0, (2, 3))
        @test_throws ArgumentError linear_to_basis(7, (2, 3))
    end

    @testset "tensor products, decomposition sums, and Kronecker sums" begin
        a = Float32[1, 2]
        b = Float32[3, 4, 5]
        @test tensor_product(a, b) == kron(a, b)
        @test eltype(tensor_product(a, b)) == Float32
        @test tensor_product(a) == a
        @test tensor_product(a) !== a
        @test tensor_product(a; copies=3) == kron(kron(a, a), a)
        @test tensor_power(a, 0) === 1.0f0
        @test length(methods(tensor_product, (Vector{Float32},))) == 1
        @test_throws ArgumentError tensor_product(a, b; copies=2)
        @test_throws ArgumentError tensor_power(a, -1)

        sparse_a = sparse([1 0; 0 2])
        sparse_b = sparse([0 3; 4 0])
        @test issparse(tensor_product(sparse_a, sparse_b))
        @test Matrix(tensor_product(sparse_a, sparse_b)) ==
            kron(Matrix(sparse_a), Matrix(sparse_b))

        columns = Matrix{Int}(I, 2, 2)
        @test tensor_sum(columns, columns) == [1, 0, 0, 1]
        @test tensor_sum(columns, columns; weights=[2, -1]) == [2, 0, 0, -1]
        matrix_terms = ([1 0; 0 0], [0 0; 0 1])
        @test tensor_sum(matrix_terms, matrix_terms) ==
            kron(matrix_terms[1], matrix_terms[1]) +
              kron(matrix_terms[2], matrix_terms[2])
        @test_throws DimensionMismatch tensor_sum(columns, columns[:, 1:1])
        @test_throws DimensionMismatch tensor_sum(columns, columns; weights=[1])

        A = [1 2; 3 4]
        B = [5 6 7; 8 9 10; 11 12 13]
        expected_sum = kron(A, Matrix{Int}(I, 3, 3)) + kron(Matrix{Int}(I, 2, 2), B)
        @test kronecker_sum(A, B) == expected_sum
        @test kronecker_sum(A; copies=2) ==
            kron(A, Matrix{Int}(I, 2, 2)) + kron(Matrix{Int}(I, 2, 2), A)
        @test issparse(kronecker_sum(sparse_a, sparse_b))
        @test_throws DimensionMismatch kronecker_sum(ones(2, 3), A)
        @test_throws ArgumentError kronecker_sum(A; copies=0)
    end

    @testset "subsystem permutations and swaps" begin
        dims = (2, 3, 2)
        permutation = (3, 1, 2)
        vector = collect(1:prod(dims))
        plan = SubsystemPermutationPlan(dims, permutation)
        permuted = permute_subsystems(vector, plan)
        @test plan.output_layout.dims == (2, 2, 3)
        @test_throws Base.CanonicalIndexError setindex!(plan.forward, 1, 1)
        @test permute_subsystems(vector, plan) == permuted
        @test_throws MethodError SubsystemPermutationPlan(
            plan.layout, plan.permutation, plan.output_layout, copy(plan.forward)
        )
        @test_throws MethodError SubsystemPermutationPlan(
            Val(:validated),
            plan.layout,
            plan.permutation,
            plan.output_layout,
            fill(1, length(plan.forward)),
        )
        @test_throws ArgumentError SubsystemPermutationPlan(
            QuantumEntanglementTools._ValidatedConstructorToken(),
            plan.layout,
            plan.permutation,
            plan.output_layout,
            fill(1, length(plan.forward)),
        )
        for old_index in eachindex(vector)
            old_basis = linear_to_basis(old_index, dims)
            new_basis = Tuple(old_basis[system] for system in permutation)
            new_index = basis_to_linear(new_basis, plan.output_layout)
            @test permuted[new_index] == vector[old_index]
        end

        inverse_plan = SubsystemPermutationPlan(
            plan.output_layout, permutation; inverse=true
        )
        @test permute_subsystems(permuted, inverse_plan) == vector
        next_permutation = (2, 3, 1)
        next_plan = SubsystemPermutationPlan(plan.output_layout, next_permutation)
        composed_permutation = ntuple(
            position -> permutation[next_permutation[position]], length(dims)
        )
        @test permute_subsystems(permuted, next_plan) ==
            permute_subsystems(vector, dims; permutation=composed_permutation)
        @test permute_subsystems(vector, dims; permutation=(1, 2, 3)) == vector

        matrix = reshape(collect(1:(prod(dims) ^ 2)), prod(dims), prod(dims))
        permutation_matrix = sparse(
            plan.forward, 1:prod(dims), ones(Int, prod(dims)), prod(dims), prod(dims)
        )
        @test permute_subsystems(matrix, plan) ==
            Matrix(permutation_matrix * matrix * transpose(permutation_matrix))
        @test permute_subsystems(matrix, plan; rows_only=true) ==
            Matrix(permutation_matrix * matrix)

        sparse_matrix = sparse(matrix .* (matrix .% 7 .== 0))
        @test issparse(permute_subsystems(sparse_matrix, plan))
        @test Matrix(permute_subsystems(sparse_matrix, plan)) ==
            permute_subsystems(Matrix(sparse_matrix), plan)
        sparse_vector = sparsevec([1, 8], [2.0, -1.0], prod(dims))
        @test issparse(permute_subsystems(sparse_vector, plan))
        @test Vector(permute_subsystems(sparse_vector, plan)) ==
            permute_subsystems(Vector(sparse_vector), plan)

        rectangular = reshape(1:48, 6, 8)
        row_plan = SubsystemPermutationPlan((2, 3), (2, 1))
        column_plan = SubsystemPermutationPlan((4, 2), (2, 1))
        rectangular_permuted = permute_subsystems(rectangular, row_plan, column_plan)
        @test size(rectangular_permuted) == (6, 8)
        @test permute_subsystems(
            rectangular_permuted,
            SubsystemPermutationPlan((3, 2), (2, 1)),
            SubsystemPermutationPlan((2, 4), (2, 1)),
        ) == rectangular

        @test swap_subsystems(vector, dims, 1, 3) ==
            permute_subsystems(vector, dims; permutation=(3, 2, 1))
        operator = permutation_operator(plan; T=Int)
        @test issparse(operator)
        @test nnz(operator) == prod(dims)
        @test operator * vector == permuted
        @test permutation_operator(dims, permutation; sparse_output=false) * vector ==
            permuted
        swap = swap_operator(dims; systems=(1, 3), T=Int)
        @test swap * vector == swap_subsystems(vector, dims, 1, 3)
        @test swap * transpose(swap) == sparse(I, prod(dims), prod(dims))
        @test_throws ArgumentError swap_subsystems(vector, dims, 1, 1)
        @test_throws ArgumentError SubsystemPermutationPlan(dims, (1, 1, 3))
        @test_throws ArgumentError SubsystemPermutationPlan(dims, (1, 2))
        @test_throws ArgumentError SubsystemPermutationPlan(dims, (1, 2, 4))
        @test_throws DimensionMismatch permute_subsystems(
            ones(5), SubsystemPermutationPlan((2, 3), (2, 1))
        )
    end

    @testset "partial trace" begin
        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * adjoint(bell)
        expected_reduction = Matrix{ComplexF64}(I, 2, 2) / 2
        @test partial_trace(bell, (2, 2); trace_out=(2,)) ≈ expected_reduction
        @test partial_trace(bell_density, (2, 2); trace_out=1) ≈ expected_reduction

        all_vector = partial_trace(bell, (2, 2); trace_out=(2, 1))
        all_matrix = partial_trace(bell_density, (2, 2); trace_out=(1, 2))
        @test size(all_vector) == (1, 1)
        @test size(all_matrix) == (1, 1)
        @test all_vector[1, 1] ≈ 1
        @test all_matrix[1, 1] ≈ 1
        @test partial_trace(bell_density, (2, 2); trace_out=()) == bell_density
        @test partial_trace(bell, (2, 2); trace_out=()) == bell_density

        dims = (2, 3, 2)
        rng = MersenneTwister(17)
        operator = randn(rng, ComplexF64, prod(dims), prod(dims))
        plan = PartialTracePlan(dims, (3, 1))
        reduced = partial_trace(operator, plan)
        @test size(reduced) == (3, 3)
        @test tr(reduced) ≈ tr(operator)
        @test_throws Base.CanonicalIndexError setindex!(plan.keep_index, 1, 1)
        @test partial_trace(operator, plan) == reduced
        @test_throws MethodError PartialTracePlan(
            plan.layout,
            plan.trace_out,
            plan.keep,
            plan.output_layout,
            copy(plan.keep_index),
            copy(plan.trace_index),
            plan.trace_dimension,
            copy(plan.source_indices),
        )
        @test_throws MethodError PartialTracePlan(
            Val(:validated),
            plan.layout,
            plan.trace_out,
            plan.keep,
            plan.output_layout,
            fill(0, length(plan.keep_index)),
            copy(plan.trace_index),
            plan.trace_dimension,
            copy(plan.source_indices),
        )
        @test_throws ArgumentError PartialTracePlan(
            QuantumEntanglementTools._ValidatedConstructorToken(),
            plan.layout,
            plan.trace_out,
            plan.keep,
            plan.output_layout,
            fill(0, length(plan.keep_index)),
            copy(plan.trace_index),
            plan.trace_dimension,
            copy(plan.source_indices),
        )
        exact_operator = reshape(collect(1:(prod(dims) ^ 2)), prod(dims), prod(dims))
        @test partial_trace(exact_operator, plan) ==
            _direct_partial_trace(exact_operator, dims, (1, 3))

        exact_vector = Rational{Int}[1 // 2, 0, 0, 1 // 3]
        exact_reduction = partial_trace(exact_vector, (2, 2); trace_out=(2,))
        @test exact_reduction == Rational{Int}[1 // 4 0; 0 1 // 9]
        @test eltype(exact_reduction) == Rational{Int}
        abstractly_typed = Matrix{Number}(I, 4, 4)
        abstractly_typed[4, 4] = 2.5
        abstract_reduction = partial_trace(abstractly_typed, (2, 2); trace_out=(2,))
        @test eltype(abstract_reduction) === Number
        @test abstract_reduction == Number[2 0; 0 3.5]

        big_vector = Complex{BigFloat}[
            BigFloat("0.5") + BigFloat("0.25")im, 0, 0, BigFloat("0.75")
        ]
        @test eltype(partial_trace(big_vector, (2, 2); trace_out=(1,))) == Complex{BigFloat}

        sparse_density = sparse(bell_density)
        sparse_reduction = partial_trace(sparse_density, (2, 2); trace_out=(2,))
        @test issparse(sparse_reduction)
        @test Matrix(sparse_reduction) ≈ expected_reduction
        sparse_all = partial_trace(sparse_density, (2, 2); trace_out=(1, 2))
        @test issparse(sparse_all)
        @test size(sparse_all) == (1, 1)

        @test_throws DimensionMismatch partial_trace(ones(5), (2, 3); trace_out=1)
        @test_throws DimensionMismatch partial_trace(ones(5, 5), (2, 3); trace_out=1)
        @test_throws ArgumentError PartialTracePlan(dims, (1, 1))
        @test_throws ArgumentError PartialTracePlan(dims, 4)
    end

    @testset "partial transpose" begin
        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * adjoint(bell)
        transposed = partial_transpose(bell_density, (2, 2); systems=(2,))
        @test eigvals(Hermitian(transposed)) ≈ [-0.5, 0.5, 0.5, 0.5]
        @test partial_transpose(transposed, (2, 2); systems=(2,)) == bell_density

        dims = (2, 3, 2)
        exact_matrix = reshape(collect(1:(prod(dims) ^ 2)), prod(dims), prod(dims))
        plan = PartialTransposePlan(dims, (3, 1))
        @test_throws Base.CanonicalIndexError setindex!(plan.row_to_row, 1, 1)
        @test_throws MethodError PartialTransposePlan(
            plan.row_layout,
            plan.column_layout,
            plan.systems,
            plan.output_row_layout,
            plan.output_column_layout,
            copy(plan.row_to_row),
            copy(plan.row_to_column),
            copy(plan.column_to_row),
            copy(plan.column_to_column),
        )
        @test_throws MethodError PartialTransposePlan(
            Val(:validated),
            plan.row_layout,
            plan.column_layout,
            plan.systems,
            plan.output_row_layout,
            plan.output_column_layout,
            fill(-1, length(plan.row_to_row)),
            copy(plan.row_to_column),
            copy(plan.column_to_row),
            copy(plan.column_to_column),
        )
        @test_throws ArgumentError PartialTransposePlan(
            QuantumEntanglementTools._ValidatedConstructorToken(),
            plan.row_layout,
            plan.column_layout,
            plan.systems,
            plan.output_row_layout,
            plan.output_column_layout,
            fill(-1, length(plan.row_to_row)),
            copy(plan.row_to_column),
            copy(plan.column_to_row),
            copy(plan.column_to_column),
        )
        twice = partial_transpose(partial_transpose(exact_matrix, plan), plan)
        @test twice == exact_matrix
        @test partial_transpose(exact_matrix, plan) ==
            _direct_partial_transpose(exact_matrix, dims, (1, 3))
        @test partial_transpose(exact_matrix, dims; systems=()) == exact_matrix
        @test eltype(partial_transpose(exact_matrix, plan)) == Int

        sparse_matrix = sparse(exact_matrix .* (exact_matrix .% 11 .== 0))
        @test issparse(partial_transpose(sparse_matrix, plan))
        @test Matrix(partial_transpose(sparse_matrix, plan)) ==
            partial_transpose(Matrix(sparse_matrix), plan)

        rectangular = reshape(collect(1:48), 6, 8)
        rectangular_plan = PartialTransposePlan((2, 3), (4, 2), (1,))
        rectangular_transposed = partial_transpose(rectangular, rectangular_plan)
        @test size(rectangular_transposed) == (12, 4)
        inverse_plan = PartialTransposePlan(
            rectangular_plan.output_row_layout, rectangular_plan.output_column_layout, (1,)
        )
        @test partial_transpose(rectangular_transposed, inverse_plan) == rectangular

        @test_throws ArgumentError PartialTransposePlan(dims, (2, 2))
        @test_throws ArgumentError PartialTransposePlan(dims, 0)
        @test_throws DimensionMismatch PartialTransposePlan((2, 3), (2, 2, 2), (1,))
        @test_throws DimensionMismatch partial_transpose(ones(5, 5), (2, 3); systems=1)
    end

    @testset "realignment and inverse reshuffling" begin
        basis_operator = zeros(Int, 6, 6)
        basis_operator[basis_to_linear((2, 3), (2, 3)), basis_to_linear((1, 2), (2, 3))] = 7
        realigned_basis = realign(basis_operator, (2, 3))
        @test size(realigned_basis) == (4, 9)
        @test realigned_basis[
            basis_to_linear((2, 1), (2, 2)), basis_to_linear((3, 2), (3, 3))
        ] == 7
        @test count(!iszero, realigned_basis) == 1

        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * adjoint(bell)
        @test realign(bell_density, (2, 2)) ≈ Matrix{ComplexF64}(I, 4, 4) / 2
        @test realignment(bell_density, (2, 2)) == realign(bell_density, (2, 2))

        dims = (2, 3, 2)
        matrix = reshape(Rational{Int}.(1:(prod(dims) ^ 2)), prod(dims), prod(dims))
        for systems in ((1,), (2, 1))
            plan = RealignmentPlan(dims; systems=systems)
            @test_throws Base.CanonicalIndexError setindex!(plan.row_a_index, 1, 1)
            @test_throws MethodError RealignmentPlan(
                plan.row_layout,
                plan.column_layout,
                plan.systems,
                plan.complement,
                plan.row_a_layout,
                plan.row_b_layout,
                plan.column_a_layout,
                plan.column_b_layout,
                copy(plan.row_a_index),
                copy(plan.row_b_index),
                copy(plan.column_a_index),
                copy(plan.column_b_index),
                copy(plan.row_from_groups),
                copy(plan.column_from_groups),
                plan.output_size,
            )
            @test_throws MethodError RealignmentPlan(
                Val(:validated),
                plan.row_layout,
                plan.column_layout,
                plan.systems,
                plan.complement,
                plan.row_a_layout,
                plan.row_b_layout,
                plan.column_a_layout,
                plan.column_b_layout,
                fill(0, length(plan.row_a_index)),
                copy(plan.row_b_index),
                copy(plan.column_a_index),
                copy(plan.column_b_index),
                copy(plan.row_from_groups),
                copy(plan.column_from_groups),
                plan.output_size,
            )
            @test_throws ArgumentError RealignmentPlan(
                QuantumEntanglementTools._ValidatedConstructorToken(),
                plan.row_layout,
                plan.column_layout,
                plan.systems,
                plan.complement,
                plan.row_a_layout,
                plan.row_b_layout,
                plan.column_a_layout,
                plan.column_b_layout,
                fill(0, length(plan.row_a_index)),
                copy(plan.row_b_index),
                copy(plan.column_a_index),
                copy(plan.column_b_index),
                copy(plan.row_from_groups),
                copy(plan.column_from_groups),
                plan.output_size,
            )
            aligned = realign(matrix, plan)
            @test inverse_realign(aligned, plan) == matrix
            @test inverse_realignment(aligned, plan) == matrix
            @test inverse_reshuffle(reshuffle(matrix, plan), plan) == matrix
        end
        @test realign(matrix, dims; systems=(1, 2)) != realign(matrix, dims; systems=(2, 1))

        sparse_matrix = sparse(matrix .* (matrix .% 13 .== 0))
        sparse_plan = RealignmentPlan(dims; systems=(2,))
        sparse_aligned = realign(sparse_matrix, sparse_plan)
        @test issparse(sparse_aligned)
        @test issparse(inverse_realign(sparse_aligned, sparse_plan))
        @test inverse_realign(sparse_aligned, sparse_plan) == sparse_matrix

        rectangular = reshape(1:48, 6, 8)
        rectangular_plan = RealignmentPlan((2, 3), (4, 2); systems=(1,))
        aligned = realign(rectangular, rectangular_plan)
        @test size(aligned) == (8, 6)
        @test inverse_realign(aligned, rectangular_plan) == rectangular

        @test_throws ArgumentError RealignmentPlan((6,); systems=(1,))
        @test_throws ArgumentError RealignmentPlan((2, 3); systems=())
        @test_throws ArgumentError RealignmentPlan((2, 3); systems=(1, 2))
        @test_throws ArgumentError RealignmentPlan((2, 3); systems=(1, 1))
        @test_throws DimensionMismatch inverse_realign(ones(3, 3), RealignmentPlan((2, 2)))
    end

    @testset "plan invariants and array axes" begin
        vector = _ZeroBasedVector(collect(1:4))
        matrix = _ZeroBasedMatrix(reshape(collect(1:16), 4, 4))
        permutation_plan = SubsystemPermutationPlan((2, 2), (2, 1))
        trace_plan = PartialTracePlan((2, 2), (2,))
        transpose_plan = PartialTransposePlan((2, 2), (2,))
        realignment_plan = RealignmentPlan((2, 2))

        @test_throws ArgumentError permute_subsystems(vector, permutation_plan)
        @test_throws ArgumentError permute_subsystems(matrix, permutation_plan)
        @test_throws ArgumentError permute_subsystems(
            matrix, permutation_plan; rows_only=true
        )
        @test_throws ArgumentError partial_trace(vector, trace_plan)
        @test_throws ArgumentError partial_trace(matrix, trace_plan)
        @test_throws ArgumentError partial_transpose(matrix, transpose_plan)
        @test_throws ArgumentError realign(matrix, realignment_plan)
        @test_throws ArgumentError inverse_realign(matrix, realignment_plan)
        @test_throws ArgumentError tensor_product(vector, ones(2))
        @test_throws ArgumentError tensor_power(vector, 2)
        @test_throws ArgumentError kronecker_sum(matrix)
        @test_throws ArgumentError tensor_sum(matrix, ones(2, 2))
        @test_throws ArgumentError tensor_sum(
            ones(2, 2), ones(2, 2); weights=_ZeroBasedVector([1, 1])
        )
        @test_throws ArgumentError tensor_sum(
            _ZeroBasedVector([[1, 0], [0, 1]]), ([1, 0], [0, 1])
        )
    end

    @testset "symmetric and antisymmetric projectors" begin
        for dimension in 1:3, copies in 1:3
            symmetric = symmetric_projector(dimension, copies)
            antisymmetric = antisymmetric_projector(dimension, copies)
            @test issparse(symmetric)
            @test issparse(antisymmetric)
            @test symmetric == adjoint(symmetric)
            @test antisymmetric == adjoint(antisymmetric)
            @test symmetric * symmetric ≈ symmetric
            @test antisymmetric * antisymmetric ≈ antisymmetric

            symmetric_basis = symmetric_subspace_basis(dimension, copies)
            antisymmetric_basis = antisymmetric_subspace_basis(dimension, copies)
            @test symmetric_basis * adjoint(symmetric_basis) ≈ symmetric
            @test antisymmetric_basis * adjoint(antisymmetric_basis) ≈ antisymmetric
            @test adjoint(symmetric_basis) * symmetric_basis ≈
                Matrix{Float64}(I, size(symmetric_basis, 2), size(symmetric_basis, 2))
            @test adjoint(antisymmetric_basis) * antisymmetric_basis ≈ Matrix{Float64}(
                I, size(antisymmetric_basis, 2), size(antisymmetric_basis, 2)
            )
            if copies > 1
                @test symmetric * antisymmetric ≈ spzeros(size(symmetric)...)
            end
        end

        exact_symmetric = symmetric_projector(2, 2; T=Rational{Int})
        exact_antisymmetric = antisymmetric_projector(2, 2; T=Rational{Int})
        @test symmetric_projection(2, 2; T=Rational{Int}) == exact_symmetric
        @test antisymmetric_projection(2, 2; T=Rational{Int}) == exact_antisymmetric
        @test eltype(exact_symmetric) == Rational{Int}
        @test exact_symmetric * exact_symmetric == exact_symmetric
        @test exact_antisymmetric * exact_antisymmetric == exact_antisymmetric
        @test !issparse(symmetric_projector(2; sparse_output=false))
        @test eltype(symmetric_projector(2; T=Float32)) == Float32
        @test eltype(antisymmetric_projector(2; T=BigFloat)) == BigFloat
        @test size(antisymmetric_subspace_basis(2, 3)) == (8, 0)
        @test iszero(antisymmetric_projector(2, 3))
        @test_throws ArgumentError symmetric_projector(0)
        @test_throws ArgumentError antisymmetric_projector(2, 0)
        @test_throws ArgumentError symmetric_projector(true)
    end

    @testset "MATLAB compatibility wrappers" begin
        A = [1 2; 3 4]
        B = [0 1; 1 0]
        @test MATLABCompat.Tensor(A, B) == tensor_product(A, B)
        @test MATLABCompat.Tensor(A, 3) == tensor_power(A, 3)
        @test MATLABCompat.Tensor((A, B)) == tensor_product(A, B)

        columns = Matrix{Int}(I, 2, 2)
        @test MATLABCompat.TensorSum(columns, columns) == tensor_sum(columns, columns)
        @test MATLABCompat.TensorSum([2, -1], columns, columns) ==
            tensor_sum(columns, columns; weights=[2, -1])
        @test_throws ArgumentError MATLABCompat.TensorSum(
            _ZeroBasedVector([2, -1]), columns, columns
        )
        @test MATLABCompat.KroneckerSum(A, B) == kronecker_sum(A, B)
        @test MATLABCompat.KroneckerSum(A, 2) == kronecker_sum(A; copies=2)

        vector = collect(1:6)
        compat_permuted = MATLABCompat.PermuteSystems(vector, (2, 1), (2, 3))
        @test compat_permuted == permute_subsystems(vector, (2, 3); permutation=(2, 1))
        @test MATLABCompat.PermuteSystems(compat_permuted, (2, 1), (3, 2), 0, 1) == vector

        rectangular = reshape(1:48, 6, 8)
        dim_matrix = [2 3; 4 2]
        @test MATLABCompat.PermuteSystems(rectangular, (2, 1), dim_matrix, 1) ==
            permute_subsystems(
            rectangular, SubsystemPermutationPlan((2, 3), (2, 1)); rows_only=true
        )
        @test_throws ArgumentError MATLABCompat.PermuteSystems(
            rectangular, (2, 1), _ZeroBasedMatrix(dim_matrix), 1
        )
        @test MATLABCompat.Swap(vector, (1, 2), (2, 3)) ==
            swap_subsystems(vector, (2, 3), 1, 2)
        compat_operator = MATLABCompat.PermutationOperator((2, 3), (2, 1), 0, 1)
        @test issparse(compat_operator)
        @test compat_operator * vector == compat_permuted
        @test !issparse(MATLABCompat.PermutationOperator((2, 3), (2, 1)))
        @test MATLABCompat.SwapOperator((2, 3), 1) * vector ==
            MATLABCompat.Swap(vector, (1, 2), (2, 3))

        row_vector = transpose(vector)
        row_permuted = MATLABCompat.PermuteSystems(row_vector, (2, 1), (2, 3))
        @test size(row_permuted) == (1, 6)
        @test vec(row_permuted) == compat_permuted
        @test size(MATLABCompat.Swap(row_vector, (1, 2), (2, 3))) == (1, 6)
        sparse_row = sparse(reshape(vector, 1, :))
        sparse_row_permuted = MATLABCompat.PermuteSystems(sparse_row, (2, 1), (2, 3))
        @test issparse(sparse_row_permuted)
        @test size(sparse_row_permuted) == (1, 6)
        @test vec(Matrix(sparse_row_permuted)) == compat_permuted

        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * adjoint(bell)
        @test MATLABCompat.PartialTrace(bell, 2, (2, 2), -1) ≈
            partial_trace(bell, (2, 2); trace_out=2)
        bell_row = transpose(bell)
        @test MATLABCompat.PartialTrace(bell_row, 2, (2, 2), -1) ≈
            partial_trace(bell, (2, 2); trace_out=2)
        @test_throws ArgumentError MATLABCompat.PartialTrace(
            _ZeroBasedMatrix(reshape(bell, 1, :)), 2, (2, 2), -1
        )
        @test MATLABCompat.PartialTranspose(bell_density, 2, (2, 2)) ==
            partial_transpose(bell_density, (2, 2); systems=2)
        rectangular_pt = MATLABCompat.PartialTranspose(rectangular, 1, dim_matrix)
        @test size(rectangular_pt) == (12, 4)

        realigned = MATLABCompat.Realignment(bell_density, (2, 2))
        @test realigned == realign(bell_density, (2, 2))
        @test MATLABCompat.InverseRealignment(realigned, (2, 2)) == bell_density
        @test MATLABCompat.SymmetricProjection(2, 2, 0) == symmetric_projector(2, 2)
        symmetric_basis = MATLABCompat.SymmetricProjection(2, 2, 1)
        @test symmetric_basis * adjoint(symmetric_basis) ≈ symmetric_projector(2, 2)
        @test MATLABCompat.AntisymmetricProjection(2, 2, 0) == antisymmetric_projector(2, 2)
        @test MATLABCompat.BasisToLinear((2, 1), (2, 3)) == 4
        @test MATLABCompat.LinearToBasis(4, (2, 3)) == (2, 1)

        @test_throws ArgumentError MATLABCompat.PartialTrace(bell_density, 2, (2, 2), 4)
        @test_throws ArgumentError MATLABCompat.PermuteSystems(vector, (2, 1), (2, 3), 2)
    end
end
