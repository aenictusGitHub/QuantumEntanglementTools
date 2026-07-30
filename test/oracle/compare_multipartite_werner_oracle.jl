using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using SparseArrays
using Test

const QET = QuantumEntanglementTools
const Compat = QET.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const FIXTURE_NAMES = Set([
    "multipartite_werner_pinned_overwrite", "werner_one_entry_vector"
])

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function fixture_matrix(fixture)
    dimensions = Tuple(Int.(fixture.dims))
    length(dimensions) == 2 ||
        error("fixture $(fixture.name) does not record two dimensions")
    real_parts = Float64.(fixture.real)
    imaginary_parts = Float64.(fixture.imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "multipartite_werner_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "multipartite_werner_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "WernerState oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_multipartite_werner_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-multipartite-werner" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
occursin("retains only", String(metadata.pinned_defect)) ||
    error("fixture does not record the reviewed loop-overwrite defect")
!isempty(String(metadata.invalid_identifier)) ||
    error("fixture does not record the invalid parameter-length route")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(
    String(fixture.name) => fixture_matrix(fixture) for fixture in payload.fixtures
)
Set(keys(fixtures)) == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

parameters = Float64.(metadata.parameters)
permutations = [(1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]
identity_matrix = Matrix{Float64}(I, 27, 27)
last_permutation = Matrix(QET.permutation_operator((3, 3, 3), permutations[end]))
pinned_overwrite = identity_matrix - parameters[end] * last_permutation
pinned_overwrite ./= tr(pinned_overwrite)

intended_raw = copy(identity_matrix)
for (coefficient, permutation) in zip(parameters, permutations)
    intended_raw .-= coefficient .* Matrix(QET.permutation_operator((3, 3, 3), permutation))
end
intended = intended_raw / tr(intended_raw)
native = QET.werner_state(3, parameters; sparse_output=false)
compatibility = Matrix(Compat.WernerState(3, parameters))

@testset "QETLAB multipartite Werner differential fixture" begin
    upstream = fixtures["multipartite_werner_pinned_overwrite"]
    @test isapprox(upstream, pinned_overwrite; atol=3e-14, rtol=3e-14)
    @test !isapprox(upstream, intended; atol=3e-14, rtol=3e-14)
    @test isapprox(native, intended; atol=3e-14, rtol=3e-14)
    @test isapprox(compatibility, intended; atol=3e-14, rtol=3e-14)
    @test ishermitian(native)
    @test tr(native) ≈ 1
    @test minimum(eigvals(Hermitian(native))) >= -1e-14

    one_entry = fixtures["werner_one_entry_vector"]
    @test isapprox(one_entry, QET.werner_state(3, 0.2); atol=3e-14, rtol=3e-14)
    @test Compat.WernerState(3, [0.2]) == QET.werner_state(3, 0.2)
    @test_throws ArgumentError Compat.WernerState(2, [0.1, 0.2])
end

println(
    "Multipartite Werner oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
