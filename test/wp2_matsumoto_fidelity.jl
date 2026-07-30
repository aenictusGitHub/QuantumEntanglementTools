using LinearAlgebra
using SparseArrays

const QETMatsumoto = QuantumEntanglementTools
const CompatMatsumoto = QuantumEntanglementTools.MATLABCompat

function _matsumoto_test_pd_formula(left, right)
    decomposition = eigen(Hermitian(left))
    square_root =
        decomposition.vectors *
        Diagonal(sqrt.(decomposition.values)) *
        adjoint(decomposition.vectors)
    inverse_square_root =
        decomposition.vectors *
        Diagonal(inv.(sqrt.(decomposition.values))) *
        adjoint(decomposition.vectors)
    relative = inverse_square_root * right * inverse_square_root
    relative_decomposition = eigen(Hermitian((relative + adjoint(relative)) / 2))
    relative_root =
        relative_decomposition.vectors *
        Diagonal(sqrt.(relative_decomposition.values)) *
        adjoint(relative_decomposition.vectors)
    return real(tr(square_root * relative_root * square_root))
end

@testset "WP2 Matsumoto fidelity commuting and generic paths" begin
    rho = Diagonal([0.7, 0.2, 0.1])
    sigma = Diagonal([0.4, 0.5, 0.1])
    expected = sum(sqrt.(diag(rho)) .* sqrt.(diag(sigma)))
    @test QETMatsumoto.matsumoto_fidelity(rho, sigma) ≈ expected
    @test QETMatsumoto.matsumoto_fidelity(sigma, rho) ≈ expected
    @test QETMatsumoto.matsumoto_fidelity(rho, rho) ≈ 1
    @test CompatMatsumoto.MatsumotoFidelity(rho, sigma) ≈ expected

    pure_zero = Diagonal([1.0, 0.0])
    pure_one = Diagonal([0.0, 1.0])
    @test QETMatsumoto.matsumoto_fidelity(pure_zero, pure_zero) == 1
    @test QETMatsumoto.matsumoto_fidelity(pure_zero, pure_one) == 0

    rho32 = Diagonal(Float32[0.7, 0.2, 0.1])
    sigma32 = Diagonal(Float32[0.4, 0.5, 0.1])
    value32 = QETMatsumoto.matsumoto_fidelity(rho32, sigma32)
    @test value32 isa Float32
    @test value32 ≈ Float32(expected) atol = 8eps(Float32)

    big_rho = Diagonal(BigFloat[big"0.7", big"0.2", big"0.1"])
    big_sigma = Diagonal(BigFloat[big"0.4", big"0.5", big"0.1"])
    big_value = QETMatsumoto.matsumoto_fidelity(big_rho, big_sigma)
    @test big_value isa BigFloat
    @test big_value ≈
        sqrt(big"0.7" * big"0.4") + sqrt(big"0.2" * big"0.5") + sqrt(big"0.1" * big"0.1")
    @test QETMatsumoto.matsumoto_fidelity(
        big_rho, big_sigma; support_boundary_policy=:project
    ) == big_value

    @test_throws DimensionMismatch QETMatsumoto.matsumoto_fidelity(
        Diagonal([1.0, 0.0]), Diagonal([1.0, 0.0, 0.0])
    )
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        Diagonal(Float64[]), Diagonal(Float64[])
    )
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity(
        Diagonal([1.1, -0.1]), pure_zero
    )
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        Diagonal([0.9, 0.0]), pure_zero
    )
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        Diagonal(ComplexF64[1 + 1e-3im, 0]), pure_zero
    )
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        Diagonal([1, 0]), Diagonal([1, 0])
    )
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        rho, sigma; support_boundary_policy=:unsupported
    )
end

@testset "WP2 Matsumoto fidelity positive-definite path" begin
    rho = [0.7 0.1; 0.1 0.3]
    sigma = [0.4 0.05; 0.05 0.6]
    value = QETMatsumoto.matsumoto_fidelity(rho, sigma)
    ratio = sqrt(det(sigma) / det(rho))
    independent_2x2 = real(tr(sigma + ratio * rho)) / sqrt(real(tr(rho \ sigma)) + 2ratio)
    @test value ≈ independent_2x2
    @test value ≈ _matsumoto_test_pd_formula(rho, sigma)
    @test value ≈ QETMatsumoto.matsumoto_fidelity(sigma, rho)
    @test 0 < value < QETMatsumoto.fidelity(rho, sigma) < 1
    @test QETMatsumoto.matsumoto_fidelity(rho, rho) ≈ 1

    common_rotation = [0.6 -0.8; 0.8 0.6]
    commuting_rho = Matrix(
        Hermitian(common_rotation * Diagonal([0.8, 0.2]) * transpose(common_rotation))
    )
    commuting_sigma = Matrix(
        Hermitian(common_rotation * Diagonal([0.3, 0.7]) * transpose(common_rotation))
    )
    @test QETMatsumoto.matsumoto_fidelity(commuting_rho, commuting_sigma) ≈
        sqrt(0.8 * 0.3) + sqrt(0.2 * 0.7)

    rotation = [cos(0.37) -sin(0.37); sin(0.37) cos(0.37)]
    rotated_rho = Matrix(Hermitian(rotation * rho * transpose(rotation)))
    rotated_sigma = Matrix(Hermitian(rotation * sigma * transpose(rotation)))
    @test QETMatsumoto.matsumoto_fidelity(rotated_rho, rotated_sigma) ≈ value

    complex_unitary = ComplexF64[1 im; im 1] / sqrt(2)
    complex_rho = Matrix(Hermitian(complex_unitary * rho * adjoint(complex_unitary)))
    complex_sigma = Matrix(Hermitian(complex_unitary * sigma * adjoint(complex_unitary)))
    @test QETMatsumoto.matsumoto_fidelity(complex_rho, complex_sigma) ≈ value

    rho32 = Float32.(rho)
    sigma32 = Float32.(sigma)
    value32 = QETMatsumoto.matsumoto_fidelity(rho32, sigma32)
    @test value32 isa Float32
    @test value32 ≈ Float32(value) atol = 32eps(Float32)

    sparse_rho = sparse(rho)
    sparse_sigma = sparse(sigma)
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(sparse_rho, sparse_sigma)
    @test QETMatsumoto.matsumoto_fidelity(sparse_rho, sparse_sigma; allow_densify=true) ≈
        value
    @test CompatMatsumoto.MatsumotoFidelity(sparse_rho, sparse_sigma; allow_densify=true) ≈
        value

    near_hermitian = copy(rho)
    near_hermitian[1, 2] += 1e-10
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        near_hermitian, sigma; atol=1e-9, rtol=0
    )
    @test_throws DimensionMismatch QETMatsumoto.matsumoto_fidelity(ones(2, 3), ones(2, 3))
    @test_throws DimensionMismatch QETMatsumoto.matsumoto_fidelity(rho, Matrix(I, 3, 3) / 3)
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(0.9rho, sigma)
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity([1.1 0.0; 0.0 -0.1], sigma)
    invalid = copy(rho)
    invalid[1, 1] = Inf
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(invalid, sigma)
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(
        BigFloat.(rho), BigFloat.(sigma)
    )
end

@testset "WP2 Matsumoto fidelity singular supports" begin
    full_rank = [0.7 0.1; 0.1 0.3]
    vector = [3.0, 4.0] / 5
    pure = vector * transpose(vector)
    expected_pure = inv(sqrt(dot(vector, full_rank \ vector)))
    @test pure * pure == pure
    @test QETMatsumoto.matsumoto_fidelity(full_rank, pure) ≈ expected_pure
    @test QETMatsumoto.matsumoto_fidelity(pure, full_rank) ≈ expected_pure
    @test QETMatsumoto.matsumoto_fidelity(pure, pure) ≈ 1

    distinct_pure = [1.0 0.0; 0.0 0.0]
    @test dot(vector, [1.0, 0.0]) != 0
    @test QETMatsumoto.fidelity(pure, distinct_pure) ≈ abs(vector[1])
    @test QETMatsumoto.matsumoto_fidelity(pure, distinct_pure) == 0

    rank_two_left = [0.6 0.1 0.0; 0.1 0.4 0.0; 0.0 0.0 0.0]
    rank_two_right = [0.55 0.0 0.12; 0.0 0.0 0.0; 0.12 0.0 0.45]
    left_block = rank_two_left[1:2, 1:2]
    right_block = rank_two_right[[1, 3], [1, 3]]
    left_shorted = inv((left_block \ [1.0, 0.0])[1])
    right_shorted = inv((right_block \ [1.0, 0.0])[1])
    expected_intersection = sqrt(left_shorted * right_shorted)
    intersection_value = QETMatsumoto.matsumoto_fidelity(rank_two_left, rank_two_right)
    @test intersection_value ≈ expected_intersection
    @test intersection_value ≈
        QETMatsumoto.matsumoto_fidelity(rank_two_right, rank_two_left)

    witness = zeros(3, 3)
    witness[1, 1] = expected_intersection
    block_constraint = [rank_two_left witness; witness rank_two_right]
    @test minimum(eigvals(Hermitian(block_constraint))) ≥ -2e-15
    @test tr(witness) ≈ intersection_value

    disjoint_left = Diagonal([0.6, 0.4, 0.0, 0.0])
    disjoint_right = Diagonal([0.0, 0.0, 0.3, 0.7])
    @test QETMatsumoto.matsumoto_fidelity(Matrix(disjoint_left), Matrix(disjoint_right)) ==
        0

    boundary_size = 1e-10
    spectral_boundary = Matrix(Diagonal([1 - boundary_size, boundary_size]))
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity(
        spectral_boundary, distinct_pure; atol=1e-9, rtol=0
    )
    @test QETMatsumoto.matsumoto_fidelity(
        spectral_boundary,
        distinct_pure;
        atol=1e-9,
        rtol=0,
        support_boundary_policy=:project,
    ) ≈ sqrt(1 - boundary_size)
    @test CompatMatsumoto.MatsumotoFidelity(
        spectral_boundary,
        distinct_pure;
        atol=1e-9,
        rtol=0,
        support_boundary_policy=:project,
    ) ≈ sqrt(1 - boundary_size)

    negative_boundary = Matrix(Diagonal([1 + boundary_size, -boundary_size]))
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity(
        negative_boundary, distinct_pure; atol=1e-9, rtol=0
    )
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity(
        negative_boundary,
        distinct_pure;
        atol=1e-9,
        rtol=0,
        support_boundary_policy=:project,
    )

    pythagorean_m = 20_000
    near_cosine = (pythagorean_m^2 - 1) / (pythagorean_m^2 + 1)
    near_sine = 2pythagorean_m / (pythagorean_m^2 + 1)
    near_vector = [near_cosine, near_sine]
    near_pure = near_vector * transpose(near_vector)
    @test near_pure * near_pure == near_pure
    @test_throws DomainError QETMatsumoto.matsumoto_fidelity(distinct_pure, near_pure)
    @test QETMatsumoto.matsumoto_fidelity(
        distinct_pure, near_pure; support_boundary_policy=:project
    ) ≈ 1

    sparse_pure = sparse(pure)
    @test_throws ArgumentError QETMatsumoto.matsumoto_fidelity(sparse_pure, sparse_pure)
    @test QETMatsumoto.matsumoto_fidelity(sparse_pure, sparse_pure; allow_densify=true) ≈ 1
end

@testset "WP2 Matsumoto fidelity independent singular limit" begin
    left = [0.6 0.1 0.0; 0.1 0.4 0.0; 0.0 0.0 0.0]
    right = [0.55 0.0 0.12; 0.0 0.0 0.0; 0.12 0.0 0.45]
    target = QETMatsumoto.matsumoto_fidelity(left, right)
    regularizers = (1e-4, 1e-6, 1e-8, 1e-10)
    approximations = [
        _matsumoto_test_pd_formula(left + regularizer * I, right + regularizer * I) for
        regularizer in regularizers
    ]
    errors = abs.(approximations .- target)
    @test all(diff(errors) .< 0)
    @test last(errors) < 2e-5
end
