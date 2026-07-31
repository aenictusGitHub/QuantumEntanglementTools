using LinearAlgebra
using SparseArrays

const QETD = QuantumEntanglementTools
const CompatD = QuantumEntanglementTools.MATLABCompat

struct _TierDZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end

Base.size(vector::_TierDZeroBasedVector) = size(vector.storage)
Base.axes(vector::_TierDZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_TierDZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_TierDZeroBasedVector, index::Int) = vector.storage[index + 1]

struct _TierDZeroBasedMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    storage::M
end

Base.size(matrix::_TierDZeroBasedMatrix) = size(matrix.storage)
function Base.axes(matrix::_TierDZeroBasedMatrix)
    return (0:(size(matrix.storage, 1) - 1), 0:(size(matrix.storage, 2) - 1))
end
Base.IndexStyle(::Type{<:_TierDZeroBasedMatrix}) = IndexCartesian()
function Base.getindex(matrix::_TierDZeroBasedMatrix, row::Int, column::Int)
    return matrix.storage[row + 1, column + 1]
end

@testset "Tier D array axes validation" begin
    bell = _TierDZeroBasedVector([1.0, 0.0, 0.0, 1.0] / sqrt(2))
    density = _TierDZeroBasedMatrix(Matrix{Float64}(I, 4, 4) / 4)
    @test_throws ArgumentError QETD.trace_norm(density)
    @test_throws ArgumentError QETD.purity(density)
    @test_throws ArgumentError QETD.schmidt_coefficients(bell, (2, 2))
    @test_throws ArgumentError QETD.concurrence(bell)
    @test_throws ArgumentError CompatD.Purity(density)
    @test_throws ArgumentError CompatD.InSeparableBall(density, (2, 2))
end

@testset "Tier D matrix norms" begin
    matrix = Diagonal([3.0, 2.0, 1.0])
    @test QETD.trace_norm(matrix) == 6.0
    @test QETD.schatten_norm(matrix, 1) == 6.0
    @test QETD.schatten_norm(matrix, 2) ≈ sqrt(14.0)
    @test QETD.schatten_norm(matrix, Inf) == 3.0
    @test QETD.ky_fan_norm(matrix, 1) == 3.0
    @test QETD.ky_fan_norm(matrix, 2) == 5.0
    @test QETD.ky_fan_norm(matrix, 3) == 6.0
    @test QETD.schatten_norm(Float32[3 0; 0 1], 2) isa Float32

    rectangular = ComplexF64[1 im 0; 0 2 -im]
    singular_values = svdvals(rectangular)
    @test QETD.trace_norm(rectangular) ≈ sum(singular_values)
    @test QETD.schatten_norm(rectangular, 3) ≈ sum(singular_values .^ 3)^(1 / 3)

    @test QETD.trace_norm(zeros(0, 3)) == 0.0
    @test QETD.schatten_norm(zeros(0, 3), 2) == 0.0
    @test_throws ArgumentError QETD.ky_fan_norm(zeros(0, 3), 1)
    @test_throws ArgumentError QETD.schatten_norm(matrix, 0.5)
    @test_throws ArgumentError QETD.schatten_norm(matrix, NaN)
    @test_throws ArgumentError QETD.schatten_norm(matrix, true)
    @test_throws ArgumentError QETD.ky_fan_norm(matrix, 0)
    @test_throws ArgumentError QETD.ky_fan_norm(matrix, 4)
    @test_throws ArgumentError QETD.ky_fan_norm(matrix, 1.0)
    @test_throws ArgumentError QETD.trace_norm([1.0 NaN; 0.0 1.0])
    @test_throws ArgumentError QETD.trace_norm(BigFloat[1 0; 0 1])

    sparse_matrix = sparse(matrix)
    @test_throws ArgumentError QETD.trace_norm(sparse_matrix)
    @test QETD.trace_norm(sparse_matrix; allow_densify=true) == 6.0
    @test_throws ArgumentError QETD.schatten_norm(sparse_matrix, 2)
    @test QETD.schatten_norm(sparse_matrix, 2; allow_densify=true) ≈ sqrt(14.0)
end

@testset "Tier D density-matrix scalar measures" begin
    pure_zero = [1.0 0.0; 0.0 0.0]
    pure_one = [0.0 0.0; 0.0 1.0]
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4

    @test QETD.purity(pure_zero) == 1.0
    @test QETD.purity(maximally_mixed) == 0.25
    @test QETD.von_neumann_entropy(pure_zero; base=2) == 0.0
    @test QETD.von_neumann_entropy(maximally_mixed; base=2) == 2.0
    @test QETD.von_neumann_entropy(maximally_mixed; base=exp(1)) ≈ log(4)
    @test_throws UndefKeywordError QETD.von_neumann_entropy(maximally_mixed)
    @test_throws ArgumentError QETD.von_neumann_entropy(maximally_mixed; base=1)
    @test_throws ArgumentError QETD.von_neumann_entropy(maximally_mixed; base=-2)
    @test_throws ArgumentError QETD.von_neumann_entropy(maximally_mixed; base=0.5)

    @test QETD.fidelity(pure_zero, pure_zero) == 1.0
    @test QETD.fidelity(pure_zero, pure_one) == 0.0
    pure_plus = [0.5 0.5; 0.5 0.5]
    @test QETD.fidelity(pure_zero, pure_plus) ≈ inv(sqrt(2))
    @test QETD.fidelity(pure_zero, pure_plus; squared=true) ≈ 0.5
    rho = ComplexF64[0.7 0.1im; -0.1im 0.3]
    sigma = ComplexF64[0.4 0.05; 0.05 0.6]
    @test QETD.fidelity(rho, sigma) ≈ QETD.fidelity(sigma, rho)

    @test QETD.trace_distance(pure_zero, pure_zero) == 0.0
    @test QETD.trace_distance(pure_zero, pure_one) == 1.0
    @test QETD.trace_distance(rho, sigma) ≈ QETD.trace_distance(sigma, rho)

    @test_throws DimensionMismatch QETD.fidelity(pure_zero, maximally_mixed)
    @test_throws DimensionMismatch QETD.trace_distance(pure_zero, maximally_mixed)
    @test_throws DimensionMismatch QETD.purity(zeros(2, 3))
    @test_throws ArgumentError QETD.purity(zeros(0, 0))
    @test_throws ArgumentError QETD.purity([1.0 0.1; 0.0 0.0])
    @test_throws ArgumentError QETD.purity([0.9 0.0; 0.0 0.0])
    @test_throws DomainError QETD.purity([1.1 0.0; 0.0 -0.1])
    @test_throws ArgumentError QETD.purity([1.0 Inf; Inf 0.0])
    @test_throws ArgumentError QETD.purity(pure_zero; atol=-1)
    @test_throws ArgumentError QETD.purity(pure_zero; rtol=Inf)

    sparse_state = sparse(maximally_mixed)
    @test_throws ArgumentError QETD.purity(sparse_state)
    @test QETD.purity(sparse_state; allow_densify=true) == 0.25

    # A tiny negative eigenvalue is not silently clipped even when it lies
    # within a deliberately coarse state-validation tolerance.
    boundary_invalid = Diagonal([-1e-10, 0.5, 0.3, 0.2000000001])
    @test_throws DomainError QETD.von_neumann_entropy(
        boundary_invalid; base=2, atol=1e-9, rtol=0
    )
end

@testset "Tier D density-matrix validation diagnostics" begin
    bell = [0.5 0.0 0.0 0.5; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.5]
    original = copy(bell)
    report = QETD.validate_density_matrix(bell, (2, 2))
    @test report isa QETD.DensityMatrixValidationReport
    @test report.status === :valid
    @test report.valid
    @test report.complete
    @test report.shape == (4, 4)
    @test report.dimensions == (2, 2)
    @test report.expected_dimension == 4
    @test report.dimension_match
    @test report.finite === true
    @test report.exactly_hermitian === true
    @test report.normalized === true
    @test report.positive_semidefinite === true
    @test report.minimum_eigenvalue >= -report.spectral_tolerance
    @test report.spectral_analysis === :eigendecomposition
    @test !report.spectral_densification_required
    @test bell == original
    @test occursin("status=valid", sprint(show, report))
    @test occursin("passes", only(report.messages))
    rich_display = sprint(show, MIME("text/plain"), report)
    @test occursin("Density-matrix validation: VALID", rich_display)
    @test occursin("Shape and dimensions: pass", rich_display)
    @test occursin("Positive semidefinite: pass", rich_display)
    @test occursin("Guidance:", rich_display)

    bad_trace = QETD.validate_density_matrix(Diagonal([0.6, 0.5, 0.0, 0.0]), (2, 2))
    @test bad_trace.status === :invalid
    @test !bad_trace.valid
    @test bad_trace.complete
    @test bad_trace.normalized === false
    @test bad_trace.positive_semidefinite === true
    @test occursin("not one", join(bad_trace.messages, " "))

    negative = QETD.validate_density_matrix(Diagonal([1.1, -0.1]), (2,))
    @test negative.status === :invalid
    @test negative.positive_semidefinite === false
    @test negative.minimum_eigenvalue == -0.1
    @test occursin("No eigenvalue was clipped", join(negative.messages, " "))

    wrong_dims = QETD.validate_density_matrix(Matrix{Float64}(I, 2, 2) / 2, (2, 2))
    @test wrong_dims.status === :invalid
    @test !wrong_dims.dimension_match
    @test wrong_dims.spectral_analysis === :skipped_invalid_structure
    @test occursin("prod(dims)=4", join(wrong_dims.messages, " "))

    nonhermitian = QETD.validate_density_matrix([0.5 0.1; 0.0 0.5], (2,))
    @test nonhermitian.status === :invalid
    @test nonhermitian.hermitian_within_tolerance === false
    @test nonhermitian.minimum_eigenvalue === nothing
    @test occursin("No Hermitian part", join(nonhermitian.messages, " "))

    near_hermitian = ComplexF64[0.5 1.0e-10; 0.0 0.5]
    conservative = QETD.validate_density_matrix(near_hermitian, (2,); atol=2e-10, rtol=0)
    @test conservative.status === :incomplete
    @test conservative.hermitian_within_tolerance === true
    @test conservative.exactly_hermitian === false
    @test conservative.spectral_analysis === :skipped_nonexact_hermitian
    @test near_hermitian == ComplexF64[0.5 1.0e-10; 0.0 0.5]

    sparse_state = sparse([0.5 0.1; 0.1 0.5])
    sparse_snapshot = copy(sparse_state)
    sparse_report = QETD.validate_density_matrix(sparse_state, (2,))
    @test sparse_report.status === :incomplete
    @test sparse_report.sparse
    @test sparse_report.spectral_densification_required
    @test sparse_report.positive_semidefinite === nothing
    @test sparse_report.spectral_analysis === :requires_densification_opt_in
    @test occursin("allow_densify=true", join(sparse_report.messages, " "))
    @test sparse_state == sparse_snapshot
    sparse_display = sprint(show, MIME("text/plain"), sparse_report)
    @test occursin("INCOMPLETE", sparse_display)
    @test occursin("Positive semidefinite: not checked", sparse_display)
    @test occursin("allow_densify=true", sparse_display)

    sparse_opt_in = QETD.validate_density_matrix(sparse_state, (2,); allow_densify=true)
    @test sparse_opt_in.valid
    @test sparse_opt_in.spectral_analysis === :eigendecomposition
    @test sparse_opt_in.densification_permitted
    sparse_limited = QETD.validate_density_matrix(
        sparse_state, (2,); allow_densify=true, max_dense_entries=3
    )
    @test sparse_limited.status === :incomplete
    @test sparse_limited.spectral_analysis === :resource_limit

    sparse_diagonal = sparse(Diagonal([0.5, 0.5]))
    diagonal_report = QETD.validate_density_matrix(sparse_diagonal, (2,))
    @test diagonal_report.valid
    @test diagonal_report.spectral_analysis === :diagonal
    @test !diagonal_report.spectral_densification_required

    exact_diagonal = QETD.validate_density_matrix(
        Diagonal(Rational{Int}[1 // 3, 2 // 3]), (2,)
    )
    @test exact_diagonal.valid
    @test exact_diagonal.minimum_eigenvalue == 1 // 3
    exact_nondiagonal = QETD.validate_density_matrix(
        Rational{Int}[1//2 1//4; 1//4 1//2], (2,)
    )
    @test exact_nondiagonal.status === :incomplete
    @test exact_nondiagonal.spectral_analysis === :unsupported_eltype

    tolerance_boundary = QETD.validate_density_matrix(
        Diagonal([-1.0e-10, 0.5, 0.3, 0.2000000001]), (2, 2); atol=1e-9, rtol=0
    )
    @test tolerance_boundary.status === :valid_within_tolerance
    @test tolerance_boundary.valid
    @test tolerance_boundary.boundary_uncertain
    @test tolerance_boundary.minimum_eigenvalue < 0

    nonfinite = QETD.validate_density_matrix([1.0 NaN; NaN 0.0], (2,))
    @test nonfinite.status === :invalid
    @test nonfinite.finite === false
    @test nonfinite.trace_value === nothing

    zero_based = QETD.validate_density_matrix(
        _TierDZeroBasedMatrix(Matrix{Float64}(I, 2, 2) / 2), (2,)
    )
    @test zero_based.status === :invalid
    @test !zero_based.one_based_indexing
    @test_throws ArgumentError QETD.validate_density_matrix(
        Matrix{Float64}(I, 2, 2) / 2, (2,); atol=-1
    )
    @test_throws ArgumentError QETD.validate_density_matrix(
        Matrix{Float64}(I, 2, 2) / 2, (2,); max_dense_entries=-1
    )
end

@testset "Tier D negativity and logarithmic negativity" begin
    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    product = [1.0, 0.0, 0.0, 0.0]
    bell_density = bell * adjoint(bell)

    @test QETD.negativity(bell, (2, 2)) ≈ 0.5
    @test QETD.negativity(bell_density, [2, 2]) ≈ 0.5
    @test QETD.negativity(product, (2, 2)) == 0.0
    @test QETD.logarithmic_negativity(bell, (2, 2); base=2) ≈ 1.0
    @test QETD.logarithmic_negativity(product, (2, 2); base=2) == 0.0
    @test_throws UndefKeywordError QETD.logarithmic_negativity(bell, (2, 2))
    @test_throws DimensionMismatch QETD.negativity(bell, (4, 2))
    @test_throws ArgumentError QETD.negativity(bell, (4,))
    @test_throws ArgumentError QETD.negativity(bell, (2, 2); systems=())
    @test_throws ArgumentError QETD.negativity(bell, (2, 2); systems=(1, 2))
    @test_throws ArgumentError QETD.negativity(0.9bell, (2, 2))
end

@testset "Tier D Schmidt analysis" begin
    bell = ComplexF64[1, 0, 0, im] / sqrt(2)
    decomposition = QETD.schmidt_decomposition(bell, (2, 2))
    @test decomposition.coefficients ≈ [inv(sqrt(2)), inv(sqrt(2))]
    reconstructed = sum(
        decomposition.coefficients[index] *
        kron(decomposition.left_vectors[:, index], decomposition.right_vectors[:, index])
        for index in eachindex(decomposition.coefficients)
    )
    @test reconstructed ≈ bell
    @test QETD.schmidt_coefficients(bell, [2, 2]) ≈ decomposition.coefficients
    @test QETD.schmidt_rank(bell, (2, 2)) == 2
    @test QETD.schmidt_rank(ComplexF64[1, 0, 0, 0], (2, 2)) == 1
    @test QETD.schmidt_rank(zeros(ComplexF64, 4), (2, 2)) == 0

    rectangular = Float64[1, 2, 3, 4, 5, 6]
    rectangular_decomposition = QETD.schmidt_decomposition(rectangular, (2, 3))
    rectangular_reconstruction = sum(
        rectangular_decomposition.coefficients[index] * kron(
            rectangular_decomposition.left_vectors[:, index],
            rectangular_decomposition.right_vectors[:, index],
        ) for index in eachindex(rectangular_decomposition.coefficients)
    )
    @test rectangular_reconstruction ≈ rectangular

    near_product = Float64[1, 0, 0, 1e-7]
    @test QETD.schmidt_rank(near_product, (2, 2); atol=1e-6, rtol=0) == 1
    @test QETD.schmidt_rank(near_product, (2, 2); atol=1e-8, rtol=0) == 2

    sparse_bell = sparsevec([1, 4], [inv(sqrt(2)), inv(sqrt(2))], 4)
    @test_throws ArgumentError QETD.schmidt_coefficients(sparse_bell, (2, 2))
    @test QETD.schmidt_coefficients(sparse_bell, (2, 2); allow_densify=true) ≈
        [inv(sqrt(2)), inv(sqrt(2))]
    @test_throws DimensionMismatch QETD.schmidt_coefficients(bell, (2, 3))
    @test_throws ArgumentError QETD.schmidt_coefficients(bell, (4,))
    @test_throws ArgumentError QETD.schmidt_coefficients(ComplexF64[1, NaN, 0, 0], (2, 2))
    @test_throws ArgumentError QETD.schmidt_coefficients(BigFloat[1, 0, 0, 1], (2, 2))
end

@testset "Tier D two-qubit concurrence" begin
    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    product = [1.0, 0.0, 0.0, 0.0]
    theta = 0.31
    partially_entangled = [cos(theta), 0.0, 0.0, sin(theta)]

    @test QETD.concurrence(bell) ≈ 1.0
    @test QETD.concurrence(bell * adjoint(bell)) ≈ 1.0
    @test QETD.concurrence(product) == 0.0
    @test QETD.concurrence(partially_entangled) ≈ sin(2theta)
    @test QETD.concurrence(Matrix{Float64}(I, 4, 4) / 4) == 0.0
    big_bell = BigFloat[1, 0, 0, 1] / sqrt(big(2))
    @test QETD.concurrence(big_bell) isa BigFloat
    @test QETD.concurrence(big_bell) ≈ one(BigFloat)
    @test_throws DimensionMismatch QETD.concurrence(ones(3) / sqrt(3))
    @test_throws DimensionMismatch QETD.concurrence(Matrix{Float64}(I, 2, 2) / 2)
    @test_throws ArgumentError QETD.concurrence(0.9bell)

    sparse_bell = sparsevec([1, 4], [inv(sqrt(2)), inv(sqrt(2))], 4)
    @test QETD.concurrence(sparse_bell) ≈ 1.0
end

@testset "Tier D structured necessary criteria" begin
    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    bell_density = bell * adjoint(bell)
    product = [1.0, 0.0, 0.0, 0.0]
    product_density = product * adjoint(product)
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4

    ppt_bell = QETD.ppt_criterion(bell_density, (2, 2))
    @test ppt_bell isa QETD.CriterionResult
    @test ppt_bell.criterion === :ppt
    @test ppt_bell.status === QETD.CriterionEntanglementDetected
    @test ppt_bell.value ≈ -0.5
    @test ppt_bell.threshold == 0.0
    @test ppt_bell.witness !== nothing
    partial = QETD.partial_transpose(bell_density, (2, 2); systems=(2,))
    @test real(dot(ppt_bell.witness, partial * ppt_bell.witness)) < 0
    @test occursin("certifies entanglement", ppt_bell.message)

    ppt_mixed = QETD.ppt_criterion(maximally_mixed, (2, 2))
    @test ppt_mixed.status === QETD.CriterionSatisfied
    @test ppt_mixed.witness === nothing
    @test !occursin("is separable", lowercase(ppt_mixed.message))

    ppt_product = QETD.ppt_criterion(product_density, (2, 2))
    @test ppt_product.status === QETD.CriterionUnknown
    @test ppt_product.witness !== nothing
    @test QETD.ppt_criterion(maximally_mixed, (2, 2); atol=0.3, rtol=0).status ===
        QETD.CriterionUnknown

    realignment_bell = QETD.realignment_criterion(bell_density, (2, 2))
    @test realignment_bell.status === QETD.CriterionEntanglementDetected
    @test realignment_bell.value ≈ 2.0
    @test realignment_bell.threshold ≈ 1.0
    @test QETD.realignment_criterion(maximally_mixed, (2, 2)).status ===
        QETD.CriterionSatisfied
    @test QETD.realignment_criterion(product_density, (2, 2)).status ===
        QETD.CriterionUnknown

    reduction_bell = QETD.reduction_criterion(bell_density, (2, 2))
    @test reduction_bell.status === QETD.CriterionEntanglementDetected
    @test reduction_bell.value ≈ -0.5
    @test reduction_bell.witness.side in (:a, :b)
    @test QETD.reduction_criterion(maximally_mixed, (2, 2)).status ===
        QETD.CriterionSatisfied
    @test QETD.reduction_criterion(product_density, (2, 2)).status === QETD.CriterionUnknown
    @test QETD.reduction_criterion(bell_density, (2, 2); side=:a).status ===
        QETD.CriterionEntanglementDetected
    @test QETD.reduction_criterion(bell_density, (2, 2); side=:b).status ===
        QETD.CriterionEntanglementDetected

    boundary_invalid = Diagonal([-1e-10, 0.5, 0.3, 0.2000000001])
    boundary_result = QETD.ppt_criterion(boundary_invalid, (2, 2); atol=1e-9, rtol=0)
    @test boundary_result.status === QETD.CriterionUnknown
    @test occursin("input", boundary_result.message)
    near_hermitian = ComplexF64.(maximally_mixed)
    near_hermitian[1, 2] = 1e-9
    @test QETD.ppt_criterion(near_hermitian, (2, 2)).status === QETD.CriterionUnknown

    grossly_invalid = Diagonal([-0.1, 0.4, 0.3, 0.4])
    @test_throws DomainError QETD.ppt_criterion(grossly_invalid, (2, 2); atol=1e-9, rtol=0)
    @test_throws ArgumentError QETD.ppt_criterion(bell_density, (2, 2); systems=())
    @test_throws ArgumentError QETD.realignment_criterion(
        bell_density, (2, 2); systems=(1, 2)
    )
    @test_throws ArgumentError QETD.reduction_criterion(bell_density, (2, 2); side=:invalid)
    @test_throws ArgumentError QETD.reduction_criterion(bell_density, (2, 2, 1))

    sparse_mixed = sparse(maximally_mixed)
    @test_throws ArgumentError QETD.ppt_criterion(sparse_mixed, (2, 2))
    @test QETD.ppt_criterion(sparse_mixed, (2, 2); allow_densify=true).status ===
        QETD.CriterionSatisfied

    shown = sprint(show, ppt_bell)
    @test occursin("CriterionResult", shown)
    @test occursin("ppt", shown)
end

@testset "Tier D MATLAB/QETLAB compatibility wrappers" begin
    diagonal = Diagonal([3.0, 2.0, 1.0])
    @test CompatD.TraceNorm(diagonal) == 6.0
    @test CompatD.SchattenNorm(diagonal, 2) ≈ sqrt(14.0)
    @test CompatD.KyFanNorm(diagonal, 2) == 5.0
    @test_throws ArgumentError CompatD.TraceNorm(sparse(diagonal))
    @test CompatD.TraceNorm(sparse(diagonal); allow_densify=true) == 6.0

    # QETLAB Purity is an unchecked algebraic scalar operation. The
    # compatibility wrapper preserves that behavior while the native
    # `purity` API performs density-matrix validation.
    @test CompatD.Purity([2.0 0.0; 0.0 0.0]) == 4.0
    @test_throws DimensionMismatch CompatD.Purity(ones(2, 3))

    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    @test CompatD.Entropy(maximally_mixed) == 2.0
    @test CompatD.Entropy(maximally_mixed, exp(1), 1) ≈ log(4)
    @test CompatD.Entropy(maximally_mixed, 2, 2) == 2.0

    pure_zero = [1.0 0.0; 0.0 0.0]
    pure_plus = [0.5 0.5; 0.5 0.5]
    @test CompatD.Fidelity(pure_zero, pure_plus) ≈ inv(sqrt(2))

    bell = ComplexF64[1, 0, 0, im] / sqrt(2)
    bell_density = bell * bell'
    @test CompatD.Negativity(bell) ≈ 0.5
    @test CompatD.Negativity(bell_density, 2) ≈ 0.5

    decomposition = CompatD.SchmidtDecomposition(bell, (2, 2), -1)
    @test decomposition.coefficients ≈ [inv(sqrt(2)), inv(sqrt(2))]
    reconstructed = sum(
        decomposition.coefficients[index] *
        kron(decomposition.left_vectors[:, index], decomposition.right_vectors[:, index])
        for index in eachindex(decomposition.coefficients)
    )
    @test reconstructed ≈ bell
    @test length(
        CompatD.SchmidtDecomposition(ComplexF64[1, 0, 0, 0], (2, 2), 0).coefficients
    ) == 1
    @test length(CompatD.SchmidtDecomposition(bell, (2, 2), 1).coefficients) == 1
    @test_throws ArgumentError CompatD.SchmidtDecomposition(bell, (2, 2), 3)
    @test CompatD.SchmidtRank(bell, (2, 2)) == 2
    @test CompatD.SchmidtRank(ComplexF64[1, 0, 0, 1e-7], (2, 2), 1e-6) == 1

    @test CompatD.Concurrence(bell) ≈ 1.0
    @test CompatD.Concurrence(maximally_mixed) == 0.0

    ppt_bell = CompatD.IsPPT(bell_density, 2, (2, 2), 1e-12)
    @test ppt_bell.status === QETD.CriterionEntanglementDetected
    @test ppt_bell.value ≈ -0.5
    @test ppt_bell.witness !== nothing

    ppt_mixed = CompatD.IsPPT(maximally_mixed, 2, (2, 2), 1e-12)
    @test ppt_mixed.status === QETD.CriterionSatisfied
    @test ppt_mixed.witness === nothing

    near_hermitian = copy(maximally_mixed)
    near_hermitian[1, 2] = 1e-10
    saved_near_hermitian = copy(near_hermitian)
    ppt_near_hermitian = CompatD.IsPPT(near_hermitian, 2, (2, 2), 1e-8)
    @test ppt_near_hermitian.status === QETD.CriterionUnknown
    @test ppt_near_hermitian.witness.kind === :hermiticity_boundary
    @test ppt_near_hermitian.witness.residual == 1e-10
    @test abs(ppt_near_hermitian.witness.difference) == 1e-10
    @test ppt_near_hermitian.value == ppt_near_hermitian.witness.residual
    @test near_hermitian == saved_near_hermitian

    product_density = Diagonal([1.0, 0.0, 0.0, 0.0])
    @test CompatD.IsPPT(Matrix(product_density), 2, (2, 2), 1e-12).status ===
        QETD.CriterionUnknown

    # Unlike the native density criterion, compatibility IsPPT accepts
    # unnormalized Hermitian operators, but still returns a tri-state result.
    @test CompatD.IsPPT(Matrix{Float64}(I, 4, 4), 2, (2, 2), 1e-12).status ===
        QETD.CriterionSatisfied
    @test_throws ArgumentError CompatD.IsPPT(sparse(maximally_mixed))
    @test CompatD.IsPPT(sparse(maximally_mixed); allow_densify=true).status ===
        QETD.CriterionSatisfied
    @test_throws ArgumentError CompatD.IsPPT(
        ComplexF64[
            1 1 0 0
            0 0 0 0
            0 0 0 0
            0 0 0 0
        ],
        2,
        (2, 2),
        1e-12,
    )
end
