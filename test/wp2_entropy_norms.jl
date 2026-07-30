using LinearAlgebra
using SparseArrays

const QETWP2 = QuantumEntanglementTools
const CompatWP2 = QuantumEntanglementTools.MATLABCompat

@testset "WP2 Rényi entropy" begin
    pure = Diagonal([1.0, 0.0, 0.0])
    mixed = Diagonal([0.5, 0.25, 0.25])
    maximally_mixed = Diagonal(fill(0.25, 4))

    for alpha in (0, 0.5, 1, 2, 10, Inf)
        @test QETWP2.renyi_entropy(pure, alpha; base=2) == 0.0
        @test QETWP2.renyi_entropy(maximally_mixed, alpha; base=2) ≈ 2.0
    end
    @test QETWP2.renyi_entropy(mixed, 0; base=2) ≈ log2(3)
    @test QETWP2.renyi_entropy(mixed, 1; base=2) ≈ 1.5
    @test QETWP2.renyi_entropy(mixed, 2; base=2) ≈ -log2(3 / 8)
    @test QETWP2.renyi_entropy(mixed, Inf; base=2) == 1.0
    @test QETWP2.renyi_entropy(mixed, 1e300; base=2) ≈ 1.0
    @test QETWP2.renyi_entropy(mixed, 1; base=exp(1)) ≈ 1.5 * log(2)
    @test QETWP2.von_neumann_entropy(mixed; base=2) ==
        QETWP2.renyi_entropy(mixed, 1; base=2)

    entropy_at_one = QETWP2.renyi_entropy(mixed, 1; base=2)
    @test QETWP2.renyi_entropy(mixed, 1 - 1e-6; base=2) ≈ entropy_at_one atol = 1e-6
    @test QETWP2.renyi_entropy(mixed, 1 + 1e-6; base=2) ≈ entropy_at_one atol = 1e-6

    float32_entropy = QETWP2.renyi_entropy(
        Diagonal(Float32[0.5, 0.25, 0.25]), Float32(2); base=Float32(2)
    )
    @test float32_entropy isa Float32
    @test float32_entropy ≈ Float32(-log2(3 / 8))

    big_mixed = Diagonal(BigFloat[big"0.5", big"0.25", big"0.25"])
    big_entropy = QETWP2.renyi_entropy(big_mixed, big(2); base=big(2))
    @test big_entropy isa BigFloat
    @test big_entropy ≈ -log2(big(3) / big(8))
    @test_throws ArgumentError QETWP2.renyi_entropy(Matrix(big_mixed), 2; base=2)

    sparse_mixed = sparse(Matrix(mixed))
    @test_throws ArgumentError QETWP2.renyi_entropy(sparse_mixed, 2; base=2)
    @test QETWP2.renyi_entropy(sparse_mixed, 2; base=2, allow_densify=true) ≈ -log2(3 / 8)

    @test CompatWP2.Entropy(mixed) == 1.5
    @test CompatWP2.Entropy(mixed, 2, 0) ≈ log2(3)
    @test CompatWP2.Entropy(mixed, 2, 2) ≈ -log2(3 / 8)
    @test CompatWP2.Entropy(mixed, 2, Inf) == 1.0

    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, -1; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, -Inf; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, NaN; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, true; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, 2; base=1)
    @test_throws ArgumentError QETWP2.renyi_entropy(mixed, 2; base=Inf)
    @test_throws ArgumentError QETWP2.renyi_entropy(Diagonal([0.6, 0.3]), 2; base=2)
    @test_throws DomainError QETWP2.renyi_entropy(Diagonal([1.1, -0.1]), 2; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(ComplexF64[0.5 0.1; 0.0 0.5], 2; base=2)
    @test_throws ArgumentError QETWP2.renyi_entropy(
        Diagonal(ComplexF64[0.5 + 0.1im, 0.5 - 0.1im]), 2; base=2
    )
    @test_throws ArgumentError QETWP2.renyi_entropy(Diagonal(Float64[]), 2; base=2)
    @test_throws DomainError QETWP2.renyi_entropy(
        Diagonal([-1e-10, 0.5, 0.3, 0.2000000001]), 2; base=2, atol=1e-9, rtol=0
    )
end

@testset "WP2 top-k p-norm" begin
    signed = [-5.0, 4.0, 3.0]
    @test QETWP2.top_k_p_norm(signed, 2, 1) == 9.0
    @test QETWP2.top_k_p_norm(signed, 2, 2) == sqrt(41.0)
    @test QETWP2.top_k_p_norm(signed, 2, Inf) == 5.0
    @test QETWP2.top_k_p_norm(signed, 99, 2) == norm(signed)
    @test QETWP2.top_k_p_norm(ComplexF64[3im, -4, 0], 1, 3) ≈ 4.0
    @test QETWP2.top_k_p_norm(Float64[], 1, 2) == 0.0

    sparse_vector = sparsevec([1, 4], BigFloat[big(3), big(4)], 10)
    sparse_norm = QETWP2.top_k_p_norm(sparse_vector, 2, big(2))
    @test sparse_norm isa BigFloat
    @test sparse_norm == big(5)

    diagonal = Diagonal([3.0, 2.0, 1.0])
    @test QETWP2.top_k_p_norm(diagonal, 2, 1) == 5.0
    @test QETWP2.top_k_p_norm(diagonal, 2, 2) == sqrt(13.0)
    @test QETWP2.top_k_p_norm(diagonal, 20, 2) == sqrt(14.0)
    @test QETWP2.top_k_p_norm(Float32[3 0; 0 1], 2, 2) isa Float32

    big_diagonal = Diagonal(BigFloat[big(3), big(2), big(1)])
    @test QETWP2.top_k_p_norm(big_diagonal, 2, big(2)) isa BigFloat
    @test QETWP2.top_k_p_norm(big_diagonal, 2, big(2)) == sqrt(big(13))
    @test_throws ArgumentError QETWP2.top_k_p_norm(Matrix(big_diagonal), 2, 2)

    rectangular = ComplexF64[1 im 0; 0 2 -im]
    singular_values = svdvals(rectangular)
    @test QETWP2.top_k_p_norm(rectangular, 1, 3) == singular_values[1]
    @test QETWP2.top_k_p_norm(rectangular, 2, 3) ≈ norm(singular_values, 3)

    sparse_diagonal = sparse(diagonal)
    @test_throws ArgumentError QETWP2.top_k_p_norm(sparse_diagonal, 2, 2)
    @test QETWP2.top_k_p_norm(sparse_diagonal, 2, 2; allow_densify=true) ≈ sqrt(13.0)

    @test CompatWP2.kpNorm(reshape(signed, 1, :), 2, 1) == 9.0
    @test CompatWP2.kpNorm(reshape(signed, :, 1), 2, 1) == 9.0
    @test CompatWP2.kpNorm(sparse(reshape(signed, 1, :)), 2, 1) == 9.0
    @test CompatWP2.kpNorm(diagonal, 2, 1) == 5.0

    @test_throws ArgumentError QETWP2.top_k_p_norm(signed, 0, 2)
    @test_throws ArgumentError QETWP2.top_k_p_norm(signed, true, 2)
    @test_throws ArgumentError QETWP2.top_k_p_norm(signed, 2, 0.5)
    @test_throws ArgumentError QETWP2.top_k_p_norm(signed, 2, NaN)
    @test_throws ArgumentError QETWP2.top_k_p_norm(signed, 2, true)
    @test_throws ArgumentError QETWP2.top_k_p_norm([1.0, NaN], 1, 2)
end

@testset "WP2 dual top-k p-norm" begin
    values = [3.0, 2.0, 1.0]
    @test QETWP2.top_k_p_norm_dual(values, 2, 1) == 3.0
    @test QETWP2.top_k_p_norm_dual(values, 2, 2) ≈ sqrt(18.0)
    @test QETWP2.top_k_p_norm_dual(values, 2, 3) ≈ 3 * 2^(2 / 3)
    @test QETWP2.top_k_p_norm_dual(values, 2, Inf) == 6.0
    @test QETWP2.top_k_p_norm_dual(values, 1, 3) == 6.0
    @test QETWP2.top_k_p_norm_dual(values, 3, 1) == 3.0
    @test QETWP2.top_k_p_norm_dual(values, 3, 2) ≈ sqrt(14.0)
    @test QETWP2.top_k_p_norm_dual(values, 3, 3) ≈ norm(values, 3 / 2)
    @test QETWP2.top_k_p_norm_dual(values, 30, 2) ≈ sqrt(14.0)
    @test QETWP2.top_k_p_norm_dual(Float64[], 1, 2) == 0.0
    big_dual = QETWP2.top_k_p_norm_dual(BigFloat[3, 2, 1], 2, big(2))
    @test big_dual isa BigFloat
    @test big_dual ≈ sqrt(big(18))

    diagonal = Diagonal(values)
    @test QETWP2.top_k_p_norm_dual(diagonal, 2, 2) ≈ sqrt(18.0)
    @test CompatWP2.kpNormDual(diagonal, 2, 2) ≈ sqrt(18.0)
    @test CompatWP2.kpNormDual(reshape(values, 1, :), 2, 2) ≈ sqrt(18.0)

    sparse_values = sparsevec([1], [2.0], 5)
    @test QETWP2.top_k_p_norm_dual(sparse_values, 3, 1) == 2.0
    @test QETWP2.top_k_p_norm_dual(sparse_values, 3, 2) == 2.0
    huge_ambient_sparse = sparsevec([1], [2.0], 100_000_000)
    @test QETWP2.top_k_p_norm_dual(huge_ambient_sparse, 100_000_000, 2) == 2.0
    @test isfinite(QETWP2.top_k_p_norm_dual(fill(typemax(Int), 2), 2, Inf))

    # Independent finite search over nonnegative directions. The exact
    # optimizer [1, 1, 1] is included, and every sampled direction obeys
    # the defining duality inequality.
    dual_value = QETWP2.top_k_p_norm_dual(values, 2, 2)
    brute_force = 0.0
    sampled_inequality_holds = true
    for first_entry in 0:12, second_entry in 0:12, third_entry in 0:12
        direction = Float64[first_entry, second_entry, third_entry]
        all(iszero, direction) && continue
        primal_value = QETWP2.top_k_p_norm(direction, 2, 2)
        ratio = dot(values, direction) / primal_value
        brute_force = max(brute_force, ratio)
        sampled_inequality_holds &= ratio <= dual_value + 16eps(dual_value)
    end
    @test sampled_inequality_holds
    @test brute_force ≈ dual_value

    for direction in ([1.0, -2.0, 0.5], [-3.0, 1.0, 4.0], [0.25, 0.75, -1.5])
        @test abs(dot(values, direction)) <=
            dual_value * QETWP2.top_k_p_norm(direction, 2, 2) + 16eps(dual_value)
    end

    @test_throws ArgumentError QETWP2.top_k_p_norm_dual(values, 0, 2)
    @test_throws ArgumentError QETWP2.top_k_p_norm_dual(values, 2, 0.5)
end

@testset "WP2 Schmidt k-norm" begin
    rectangular = [3.0, 0.0, 0.0, 0.0, 0.0, 4.0]
    @test QETWP2.schmidt_coefficients(rectangular, (2, 3)) ≈ [4.0, 3.0]
    @test QETWP2.schmidt_k_norm(rectangular, (2, 3), 1) == 4.0
    @test QETWP2.schmidt_k_norm(rectangular, (2, 3), 2) == 5.0
    @test QETWP2.schmidt_k_norm(rectangular, (2, 3), 20) == 5.0
    @test QETWP2.schmidt_k_norm(zeros(6), (2, 3), 1) == 0.0

    sparse_rectangular = sparsevec([1, 6], [3.0, 4.0], 6)
    @test_throws ArgumentError QETWP2.schmidt_k_norm(sparse_rectangular, (2, 3), 1)
    @test QETWP2.schmidt_k_norm(sparse_rectangular, (2, 3), 1; allow_densify=true) == 4.0
    @test QETWP2.schmidt_k_norm(sparse_rectangular, (2, 3), 2) == 5.0

    big_rectangular = BigFloat[3, 0, 0, 0, 0, 4]
    @test QETWP2.schmidt_k_norm(big_rectangular, (2, 3), 2) isa BigFloat
    @test QETWP2.schmidt_k_norm(big_rectangular, (2, 3), 2) == big(5)
    @test_throws ArgumentError QETWP2.schmidt_k_norm(big_rectangular, (2, 3), 1)

    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    @test CompatWP2.SkVectorNorm(bell) ≈ inv(sqrt(2))
    @test CompatWP2.SkVectorNorm(bell, 2) ≈ 1.0
    @test CompatWP2.SkVectorNorm(rectangular) == 4.0
    @test CompatWP2.SkVectorNorm(rectangular, 1, 2) == 4.0
    @test CompatWP2.SkVectorNorm(rectangular, 2, (2, 3)) == 5.0

    @test_throws ArgumentError QETWP2.schmidt_k_norm(rectangular, (2, 3), 0)
    @test_throws ArgumentError QETWP2.schmidt_k_norm(rectangular, (2, 3), true)
    @test_throws DimensionMismatch QETWP2.schmidt_k_norm(rectangular, (2, 2), 1)
    @test_throws ArgumentError CompatWP2.SkVectorNorm(ones(3))
    @test_throws ArgumentError CompatWP2.SkVectorNorm(rectangular, 1, 4)
end
