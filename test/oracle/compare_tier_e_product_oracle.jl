using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QuantumEntanglementTools.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"

function fixture_matrix(fixture)
    rows = Int(fixture.rows)
    columns = Int(fixture.columns)
    real_parts = fixture.real isa Number ? [Float64(fixture.real)] : Float64.(fixture.real)
    imaginary_parts = if fixture.imaginary isa Number
        [Float64(fixture.imaginary)]
    else
        Float64.(fixture.imaginary)
    end
    length(real_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, rows, columns)
end

as_matrix(value::AbstractVector) = reshape(collect(value), :, 1)
as_matrix(value::Number) = reshape([value], 1, 1)

function oracle_inputs()
    e11 = [1.0 0.0; 0.0 0.0]
    e22 = [0.0 0.0; 0.0 1.0]
    basis_12 = zeros(3, 3)
    basis_12[1, 2] = 1
    basis_23 = zeros(3, 3)
    basis_23[2, 3] = 1
    operator = 3 * QET.tensor_product(e11, basis_12) + 2 * QET.tensor_product(e22, basis_23)

    rectangular_left = [1.0 0.0 0.0; 0.0 0.0 0.0]
    rectangular_operator = QET.tensor_product(rectangular_left, e11)
    hermitian_operator = Matrix(Diagonal(collect(1.0:6.0)))

    product_vector = QET.tensor_product(
        ComplexF64[1, 2im], ComplexF64[3, -1], ComplexF64[2, im]
    )
    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)

    first_operator = [1.0 2.0; 3.0 4.0]
    second_operator = [0.0 1.0 2.0; 3.0 4.0 5.0]
    third_operator = reshape([1.0, 2.0], 2, 1)
    product_operator = QET.tensor_product(first_operator, second_operator, third_operator)
    rank_two_operator = QET.tensor_product(e11, e11) + QET.tensor_product(e22, e22)

    angle = 0.31
    rectangular_pure = [cos(angle), 0.0, 0.0, 0.0, sin(angle), 0.0]
    bell_density = bell * adjoint(bell)
    mixed_bell = 0.7 * bell_density + 0.3 * Matrix{Float64}(I, 4, 4) / 4
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    unnormalized_identity = 2 * Matrix{Float64}(I, 4, 4)

    return (;
        operator,
        rectangular_operator,
        hermitian_operator,
        product_vector,
        bell,
        product_operator,
        rank_two_operator,
        rectangular_pure,
        mixed_bell,
        maximally_mixed,
        unnormalized_identity,
    )
end

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function normwise_match(result, expected, fixture)
    atol = Float64(fixture.atol)
    rtol = Float64(fixture.rtol)
    return norm(result - expected) <= atol + rtol * norm(expected)
end

fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "tier_e_product_oracle.json")
    committed = joinpath(
        @__DIR__, "fixtures", "tier_e_product_octave_11_3_qetlab_d858961.json"
    )
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier E product-analysis oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_e_product_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "E-product-analysis" ||
    error("expected a Tier E product-analysis fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")
occursin("rectangular operators", String(metadata.reviewed_operator_schmidt_discrepancy)) ||
    error("fixture does not record the reviewed operator-Schmidt discrepancy")
occursin("NaN", String(metadata.reviewed_entformation_discrepancy)) ||
    error("fixture does not record the reviewed EntFormation discrepancy")
occursin("structured boundary", String(metadata.product_classification_semantics)) ||
    error("fixture does not record the product-classification semantic difference")
occursin("not entangled", String(metadata.separable_ball_semantics)) ||
    error("fixture does not record the sufficient-certificate semantics")

inputs = oracle_inputs()

@testset "QETLAB Tier E product-analysis differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected = fixture_matrix(fixture)

        if name == "operator_schmidt_coefficients_2x3"
            native = as_matrix(QET.operator_schmidt_coefficients(inputs.operator, (2, 3)))
            @testset "$name/native" begin
                @test size(native) == size(expected)
                @test normwise_match(native, expected, fixture)
                compatibility = as_matrix(
                    Compat.OperatorSchmidtDecomposition(inputs.operator, (2, 3), 4).coefficients,
                )
                @test normwise_match(compatibility, expected, fixture)
            end
        elseif name == "operator_schmidt_rank_2x3"
            @testset "$name/native" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 2
                @test QET.operator_schmidt_rank(inputs.operator, (2, 3)) == only(expected)
                @test Compat.OperatorSchmidtRank(inputs.operator, (2, 3)) == only(expected)
            end
        elseif name == "operator_schmidt_rectangular_qetlab_error"
            native = QET.operator_schmidt_decomposition(
                inputs.rectangular_operator, (2, 2), (3, 2)
            )
            reconstruction = QET.tensor_sum(
                native.left_factors, native.right_factors; weights=native.coefficients
            )
            @testset "$name/reviewed_discrepancy" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native.row_dims == (2, 2)
                @test native.column_dims == (3, 2)
                @test reconstruction ≈ inputs.rectangular_operator
                compatibility = Compat.OperatorSchmidtDecomposition(
                    inputs.rectangular_operator, [2 2; 3 2], -1
                )
                @test QET.tensor_sum(
                    compatibility.left_factors,
                    compatibility.right_factors;
                    weights=compatibility.coefficients,
                ) ≈ inputs.rectangular_operator
            end
        elseif name == "operator_schmidt_hermitian_2x3_qetlab_error"
            native = QET.operator_schmidt_decomposition(inputs.hermitian_operator, (2, 3))
            reconstruction = QET.tensor_sum(
                native.left_factors, native.right_factors; weights=native.coefficients
            )
            @testset "$name/reviewed_discrepancy" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test reconstruction ≈ inputs.hermitian_operator
                @test all(isfinite, native.coefficients)
                compatibility = Compat.OperatorSchmidtDecomposition(
                    inputs.hermitian_operator, (2, 3), -1
                )
                @test QET.tensor_sum(
                    compatibility.left_factors,
                    compatibility.right_factors;
                    weights=compatibility.coefficients,
                ) ≈ inputs.hermitian_operator
            end
        elseif name == "is_product_vector_three_party_product"
            native = QET.is_product_vector(inputs.product_vector, (2, 2, 2))
            @testset "$name/structured_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native.status === :within_tolerance
                @test QET.tensor_product(native.factors...) ≈ inputs.product_vector
                @test Compat.IsProductVector(inputs.product_vector, (2, 2, 2)).status ===
                    :within_tolerance
            end
        elseif name == "is_product_vector_bell_nonproduct"
            native = QET.is_product_vector(inputs.bell, (2, 2))
            @testset "$name/structured_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 0
                @test native.status === :outside_tolerance
                @test any(native.cut_residuals .> native.cut_thresholds)
                @test Compat.IsProductVector(inputs.bell, (2, 2)).status ===
                    :outside_tolerance
            end
        elseif name == "is_product_operator_rectangular_three_party_product"
            native = QET.is_product_operator(inputs.product_operator, (2, 2, 2), (2, 3, 1))
            @testset "$name/structured_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native.status === :within_tolerance
                @test QET.tensor_product(native.factors...) ≈ inputs.product_operator
                @test Compat.IsProductOperator(inputs.product_operator, [2 2 2; 2 3 1]).status ===
                    :within_tolerance
            end
        elseif name == "is_product_operator_rank_two_nonproduct"
            native = QET.is_product_operator(inputs.rank_two_operator, (2, 2))
            @testset "$name/structured_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 0
                @test native.status === :outside_tolerance
                @test any(native.cut_residuals .> native.cut_thresholds)
                @test Compat.IsProductOperator(inputs.rank_two_operator, (2, 2)).status ===
                    :outside_tolerance
            end
        elseif name == "entformation_pure_2x3"
            native = as_matrix(
                QET.entanglement_of_formation(inputs.rectangular_pure, (2, 3))
            )
            @testset "$name/native" begin
                @test size(native) == size(expected)
                @test normwise_match(native, expected, fixture)
                @test normwise_match(
                    as_matrix(Compat.EntFormation(inputs.rectangular_pure, (2, 3))),
                    expected,
                    fixture,
                )
            end
        elseif name == "entformation_mixed_bell"
            native = as_matrix(QET.entanglement_of_formation(inputs.mixed_bell, (2, 2)))
            @testset "$name/native" begin
                @test size(native) == size(expected)
                @test normwise_match(native, expected, fixture)
                @test normwise_match(
                    as_matrix(Compat.EntFormation(inputs.mixed_bell, (2, 2))),
                    expected,
                    fixture,
                )
            end
        elseif name == "entformation_zero_concurrence_qetlab_nan"
            native = QET.entanglement_of_formation(inputs.maximally_mixed, (2, 2))
            @testset "$name/reviewed_discrepancy" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native == 0
                @test isfinite(native)
                @test Compat.EntFormation(inputs.maximally_mixed, (2, 2)) == 0
            end
        elseif name == "in_separable_ball_maximally_mixed"
            native = QET.in_separable_ball(inputs.maximally_mixed, (2, 2))
            @testset "$name/certificate_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native.status === :separable_certified
                @test native.purity < native.boundary - native.tolerance
                @test Compat.InSeparableBall(inputs.maximally_mixed).status ===
                    :separable_certified
            end
        elseif name == "in_separable_ball_pure_spectrum_outside"
            native = QET.in_separable_ball([1.0, 0.0, 0.0, 0.0], (2, 2))
            @testset "$name/certificate_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 0
                @test native.status === :outside_ball
                @test occursin("no entanglement conclusion", lowercase(native.message))
                @test !occursin("entangled", lowercase(native.message))
                @test Compat.InSeparableBall([1.0, 0.0, 0.0, 0.0],).status === :outside_ball
            end
        elseif name == "in_separable_ball_unnormalized_qetlab_normalizes"
            @testset "$name/reviewed_semantic_difference" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test_throws ArgumentError QET.in_separable_ball(
                    inputs.unnormalized_identity, (2, 2)
                )
                @test Compat.InSeparableBall(inputs.unnormalized_identity).status ===
                    :separable_certified
            end
        else
            error("unknown Tier E product-analysis oracle fixture: $name")
        end
    end

    @testset "phase-independent operator Schmidt properties" begin
        decomposition = QET.operator_schmidt_decomposition(inputs.operator, (2, 3))
        reconstruction = QET.tensor_sum(
            decomposition.left_factors,
            decomposition.right_factors;
            weights=decomposition.coefficients,
        )
        @test reconstruction ≈ inputs.operator
        @test decomposition.coefficients ≈ [3.0, 2.0, 0.0, 0.0]
    end

    @testset "native tolerance boundary remains tri-state" begin
        boundary = QET.in_separable_ball([1 / 3, 1 / 3, 1 / 3, 0.0], (2, 2))
        @test boundary.status === :unknown
        @test occursin("tolerance band", lowercase(boundary.message))
    end
end

println(
    "Verified $(length(payload.fixtures)) Tier E product-analysis fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
println(
    "Product booleans were compared only away from tolerance boundaries; ",
    "native boundary behavior remains structured.",
)
println(
    "Separable-ball failure is preserved as failure of a sufficient ",
    "certificate, never as an entanglement conclusion.",
)
if String(metadata.engine) == "Octave"
    println(
        "Octave evidence is function-specific and supplemental; ",
        "it is not a general MATLAB-equivalence claim.",
    )
end
