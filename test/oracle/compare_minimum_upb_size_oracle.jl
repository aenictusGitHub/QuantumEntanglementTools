using JSON3
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QET.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const FIXTURE_ARGUMENTS = Dict(
    "minimum_upb_2x3" => [2, 3],
    "minimum_upb_3x3x3" => [3, 3, 3],
    "minimum_upb_2x3x3" => [2, 3, 3],
    "minimum_upb_4x6" => [4, 6],
    "minimum_upb_four_qubits" => fill(2, 4),
    "minimum_upb_six_qubits" => fill(2, 6),
    "minimum_upb_eight_qubits" => fill(2, 8),
    "minimum_upb_twelve_qubits" => fill(2, 12),
    "minimum_upb_2x2x3" => [2, 2, 3],
    "minimum_upb_2x2x5" => [2, 2, 5],
    "minimum_upb_2x2x9" => [2, 2, 9],
)

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function scalar_fixture(fixture)
    Tuple(Int.(fixture.dims)) == (1, 1) || error("fixture $(fixture.name) is not scalar")
    real_value = fixture.real isa Number ? Float64(fixture.real) : only(fixture.real)
    imaginary_value =
        fixture.imaginary isa Number ? Float64(fixture.imaginary) : only(fixture.imaginary)
    iszero(imaginary_value) ||
        error("fixture $(fixture.name) has an unexpected imaginary value")
    return real_value
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "minimum_upb_size_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "minimum_upb_size_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "MinUPBSize oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_minimum_upb_size_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-minimum-upb-size" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Int.(metadata.unknown_dimensions) == [2, 3, 4] ||
    error("fixture does not record the reviewed unknown case")
occursin("MinSizeUnknown", String(metadata.unknown_identifier)) ||
    error("pinned unknown case did not report the reviewed identifier")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixture_names = Set(String(fixture.name) for fixture in payload.fixtures)
fixture_names == Set(keys(FIXTURE_ARGUMENTS)) ||
    error("fixture names do not match the reviewed comparator set")

@testset "QETLAB MinUPBSize differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        dimensions = FIXTURE_ARGUMENTS[name]
        expected = scalar_fixture(fixture)
        native = QET.minimum_upb_size(dimensions)
        @test String(fixture.comparison) == "exact"
        @test native.status === :known
        @test native.size == expected
        @test Compat.MinUPBSize(dimensions, 0) == expected
    end
    unknown = QET.minimum_upb_size([2, 3, 4])
    @test unknown.status === :unknown
    @test unknown.size === nothing
    @test_throws DomainError Compat.MinUPBSize([2, 3, 4], 0)
end

println(
    "MinUPBSize oracle comparison passed for ",
    length(payload.fixtures),
    " exact fixtures and one unknown case from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
