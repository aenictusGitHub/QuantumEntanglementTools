using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QET.MATLABCompat
const TWIRL_PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const TWIRL_PINNED_SOURCE_SHA = "a2ed4de3c937cb3e451d0b8e84508cdbe4b3bc23013558015d0e2b972853cd33"
const TWIRL_FIXTURE_NAMES = Set([
    "werner_p2", "isotropic_p2", "real_p2", "pauli_p2", "werner_p3", "real_p3"
])

function twirl_verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function twirl_fixture_matrix(fixture)
    dimensions = Tuple(Int.(fixture.dims))
    length(dimensions) == 2 ||
        error("fixture $(fixture.name) does not record two dimensions")
    real_parts = Float64.(fixture.real)
    imaginary_parts = Float64.(fixture.imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    return reshape(complex.(real_parts, imaginary_parts), dimensions)
end

committed_fixture = joinpath(@__DIR__, "fixtures", "twirl_octave_11_3_qetlab_d858961.json")
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "twirl_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Twirl oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_twirl_oracle.sh first",
)

digest = twirl_verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-twirl" || error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == TWIRL_PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.source_sha256) == TWIRL_PINNED_SOURCE_SHA ||
    error("fixture does not record the pinned Twirl.m SHA-256")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
!isempty(String(metadata.invalid_type_identifier)) ||
    error("fixture does not record invalid TYPE behavior")
!isempty(String(metadata.invalid_isotropic_copies_identifier)) ||
    error("fixture does not record the pinned isotropic P > 2 rejection")
!isempty(String(metadata.invalid_pauli_dimension_identifier)) ||
    error("fixture does not record invalid Pauli dimension behavior")
occursin("requires P = 2", String(metadata.corrected_copy_validation)) ||
    error("fixture does not record the corrected bipartite-copy validation")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(
    String(fixture.name) => twirl_fixture_matrix(fixture) for fixture in payload.fixtures
)
Set(keys(fixtures)) == TWIRL_FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

X2 = reshape(Float64.(1:16), 4, 4) + im * reshape(Float64.([0:7; -8:-1]), 4, 4)
X3 = reshape(Float64.(1:64), 8, 8) + im * reshape(Float64.([0:31; -32:-1]), 8, 8)

@testset "QETLAB Twirl differential fixture" begin
    cases = (
        ("werner_p2", X2, :werner, 2, "werner"),
        ("isotropic_p2", X2, :isotropic, 2, "isotropic"),
        ("real_p2", X2, :real, 2, "real"),
        ("pauli_p2", X2, :pauli, 2, "pauli"),
        ("werner_p3", X3, :werner, 3, "werner"),
        ("real_p3", X3, :real, 3, "real"),
    )
    for (name, input, kind, copies, compat_kind) in cases
        upstream = fixtures[name]
        native = QET.twirl(input; kind, copies)
        compatibility = Compat.Twirl(input, compat_kind, copies)
        @test isapprox(native, upstream; atol=8e-13, rtol=8e-13)
        @test isapprox(compatibility, upstream; atol=8e-13, rtol=8e-13)
        @test isapprox(QET.twirl(native; kind, copies), native; atol=8e-13, rtol=8e-13)
    end
end

println(
    "Twirl oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
