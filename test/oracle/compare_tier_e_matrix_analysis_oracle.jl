using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QET.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const AGREEMENT_FIXTURE_NAMES = Set([
    "majorizes_vector_true",
    "majorizes_vector_false",
    "majorizes_matrix_singular_values_true",
    "majorizes_complex_matrix_singular_values_true",
    "elem_sym_poly_order_zero",
    "elem_sym_poly_order_two",
    "elem_sym_poly_order_full",
    "elem_sym_poly_complex_order_two",
    "compound_order_zero",
    "compound_square_order_two",
    "compound_square_order_three",
    "compound_rectangular_order_two",
    "compound_complex_order_two",
    "additive_order_one",
    "additive_order_two",
    "additive_order_full_trace",
    "additive_diagonal_order_two",
])
const DISCREPANCY_FIXTURE_NAMES = Set([
    "majorizes_weak_total_qetlab_true",
    "majorizes_negative_padding_qetlab_false",
    "majorizes_row_matrix_qetlab_vector_semantics",
    "compound_rectangular_high_order_qetlab_zero_by_zero",
    "additive_order_zero_qetlab_error",
])

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

as_matrix(value::AbstractMatrix) = Matrix(value)
as_matrix(value::Number) = reshape([value], 1, 1)

function oracle_inputs()
    square_matrix = Float64[1 2 3; 4 5 6; 7 8 10]
    rectangular_matrix = Float64[1 2; 3 5; 7 11]
    complex_matrix = ComplexF64[1 + im 2; 3im 4 - im]
    polynomial_values = Float64[1, 2, 3, 4]
    complex_polynomial_values = ComplexF64[1 + 2im, -3 + im, 2 - im]
    diagonal_matrix = Diagonal(Float64[-2, 0.5, 3, 7])
    return (;
        square_matrix,
        rectangular_matrix,
        complex_matrix,
        polynomial_values,
        complex_polynomial_values,
        diagonal_matrix,
    )
end

function agreement_result(name::AbstractString, inputs)
    if name == "majorizes_vector_true"
        return QET.majorizes([4.0, 1.0, 1.0], [3.0, 2.0, 1.0])
    elseif name == "majorizes_vector_false"
        return QET.majorizes([3.0, 2.0, 1.0], [4.0, 1.0, 1.0])
    elseif name == "majorizes_matrix_singular_values_true"
        return QET.majorizes(Diagonal([3.0, 1.0]), Diagonal([2.0, 2.0]))
    elseif name == "majorizes_complex_matrix_singular_values_true"
        return QET.majorizes(Diagonal(ComplexF64[3im, -1]), Diagonal([2.0, 2.0]))
    elseif name == "elem_sym_poly_order_zero"
        return QET.elementary_symmetric_polynomial(inputs.polynomial_values, 0)
    elseif name == "elem_sym_poly_order_two"
        return QET.elementary_symmetric_polynomial(inputs.polynomial_values, 2)
    elseif name == "elem_sym_poly_order_full"
        return QET.elementary_symmetric_polynomial(inputs.polynomial_values, 4)
    elseif name == "elem_sym_poly_complex_order_two"
        return QET.elementary_symmetric_polynomial(inputs.complex_polynomial_values, 2)
    elseif name == "compound_order_zero"
        return QET.compound_matrix(inputs.square_matrix, 0)
    elseif name == "compound_square_order_two"
        return QET.compound_matrix(inputs.square_matrix, 2)
    elseif name == "compound_square_order_three"
        return QET.compound_matrix(inputs.square_matrix, 3)
    elseif name == "compound_rectangular_order_two"
        return QET.compound_matrix(inputs.rectangular_matrix, 2)
    elseif name == "compound_complex_order_two"
        return QET.compound_matrix(inputs.complex_matrix, 2)
    elseif name == "additive_order_one"
        return QET.additive_compound_matrix(inputs.square_matrix, 1)
    elseif name == "additive_order_two"
        return QET.additive_compound_matrix(inputs.square_matrix, 2)
    elseif name == "additive_order_full_trace"
        return QET.additive_compound_matrix(inputs.square_matrix, 3)
    elseif name == "additive_diagonal_order_two"
        return QET.additive_compound_matrix(inputs.diagonal_matrix, 2)
    end
    return error("unknown agreement fixture: $name")
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

committed_fixture = joinpath(
    @__DIR__, "fixtures", "tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "tier_e_matrix_analysis_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier E matrix-analysis oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_e_matrix_analysis_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "E-matrix-analysis" ||
    error("expected a Tier E matrix-analysis fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")
Int(metadata.agreement_fixture_count) == length(AGREEMENT_FIXTURE_NAMES) ||
    error("agreement fixture count does not match the reviewed comparator set")
Int(metadata.reviewed_discrepancy_fixture_count) == length(DISCREPANCY_FIXTURE_NAMES) ||
    error("discrepancy fixture count does not match the reviewed comparator set")
Bool(metadata.source_free_fixture) || error("fixture is not marked as source-free")
occursin("strong majorization", String(metadata.reviewed_majorization_total_semantics)) ||
    error("fixture does not record strong-versus-weak majorization")
occursin(
    "sorts before appending zeros", String(metadata.reviewed_majorization_padding_semantics)
) || error("fixture does not record the negative-padding discrepancy")
occursin("one-row matrices", String(metadata.reviewed_row_matrix_semantics)) ||
    error("fixture does not record the row-matrix discrepancy")
occursin("0-by-0", String(metadata.reviewed_rectangular_compound_semantics)) ||
    error("fixture does not record the rectangular-compound discrepancy")
occursin("errors at order zero", String(metadata.reviewed_additive_zero_semantics)) ||
    error("fixture does not record the additive order-zero discrepancy")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixture_names = Set(String(fixture.name) for fixture in payload.fixtures)
fixture_names == union(AGREEMENT_FIXTURE_NAMES, DISCREPANCY_FIXTURE_NAMES) ||
    error("fixture names do not match the reviewed comparator set")

inputs = oracle_inputs()

@testset "QETLAB Tier E matrix-analysis differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        classification = String(fixture.classification)
        expected = fixture_matrix(fixture)

        if name in AGREEMENT_FIXTURE_NAMES
            classification == "agreement" ||
                error("$name is not classified as an agreement fixture")
            result = as_matrix(agreement_result(name, inputs))
            @testset "$name/agreement" begin
                @test size(result) == size(expected)
                if String(fixture.comparison) == "exact"
                    @test result == expected
                elseif String(fixture.comparison) == "normwise"
                    @test normwise_match(result, expected, fixture)
                else
                    error("unknown comparison mode for $name")
                end
            end
        elseif name in DISCREPANCY_FIXTURE_NAMES
            classification == "reviewed_discrepancy" ||
                error("$name is not classified as a reviewed discrepancy")
            if name == "majorizes_weak_total_qetlab_true"
                native_result = QET.majorizes([2, 0], [1, 0]; rtol=0)
                @testset "$name/reviewed_discrepancy" begin
                    @test size(expected) == (1, 1)
                    @test only(expected) == 1
                    @test native_result === false
                    @test native_result != Bool(only(expected))
                    @test Compat.Majorizes([2, 0], [1, 0]; rtol=0)
                end
            elseif name == "majorizes_negative_padding_qetlab_false"
                native_result = QET.majorizes([1, -1], [1, 0, -1]; rtol=0)
                @testset "$name/reviewed_discrepancy" begin
                    @test size(expected) == (1, 1)
                    @test only(expected) == 0
                    @test native_result === true
                    @test native_result != Bool(only(expected))
                    @test !Compat.Majorizes([1, -1], [1, 0, -1]; rtol=0)
                end
            elseif name == "majorizes_row_matrix_qetlab_vector_semantics"
                native_result = QET.majorizes(
                    [1.0 1.0], [sqrt(2.0) 0.0]; atol=1e-14, rtol=0
                )
                @testset "$name/reviewed_discrepancy" begin
                    @test size(expected) == (1, 1)
                    @test only(expected) == 0
                    @test native_result === true
                    @test native_result != Bool(only(expected))
                    @test !Compat.Majorizes([1.0 1.0], [sqrt(2.0) 0.0]; atol=1e-14, rtol=0)
                end
            elseif name == "compound_rectangular_high_order_qetlab_zero_by_zero"
                native_result = QET.compound_matrix(zeros(2, 3), 3)
                @testset "$name/reviewed_discrepancy" begin
                    @test size(expected) == (0, 0)
                    @test size(native_result) == (0, 1)
                    @test isempty(native_result)
                    @test size(native_result) != size(expected)
                    @test size(Compat.CompoundMatrix(zeros(2, 3), 3)) == size(expected)
                end
            elseif name == "additive_order_zero_qetlab_error"
                native_result = QET.additive_compound_matrix(Matrix{Float64}(I, 2, 2), 0)
                @testset "$name/reviewed_discrepancy" begin
                    @test size(expected) == (1, 1)
                    @test only(expected) == 1
                    @test size(native_result) == (1, 1)
                    @test native_result == zeros(1, 1)
                    @test_throws ArgumentError Compat.AdditiveCompoundMatrix(
                        Matrix{Float64}(I, 2, 2), 0
                    )
                end
            end
        end
    end
end

println(
    "Verified $(length(payload.fixtures)) Tier E matrix-analysis fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
println(
    "$(metadata.agreement_fixture_count) fixtures were compared normally; ",
    "$(metadata.reviewed_discrepancy_fixture_count) reviewed QETLAB behaviors ",
    "were checked without becoming native expected values.",
)
if String(metadata.engine) == "Octave"
    println(
        "Octave evidence is function-specific and supplemental; ",
        "it is not a general MATLAB-equivalence claim.",
    )
end
