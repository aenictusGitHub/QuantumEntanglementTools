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
    "is_entangling_gate_identity_flag",
    "is_entangling_gate_swap_flag",
    "is_entangling_gate_cnot_flag",
    "is_entangling_gate_cnot_witness",
    "is_entangling_gate_controlled_z_flag",
    "is_entangling_gate_controlled_z_last_candidate",
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

function fixture_array(fixture)
    dimensions = Tuple(Int.(fixture.dims))
    real_parts = fixture.real isa Number ? Float64[fixture.real] : Float64.(fixture.real)
    imaginary_parts = if fixture.imaginary isa Number
        Float64[fixture.imaginary]
    else
        Float64.(fixture.imaginary)
    end
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "entangling_gate_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "entangling_gate_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "IsEntanglingGate oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_entangling_gate_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-entangling-gate" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
occursin("not a valid witness", String(metadata.witness_review)) ||
    error("fixture does not record the reviewed controlled-Z witness defect")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(
    String(fixture.name) => fixture_array(fixture) for fixture in payload.fixtures
)
Set(keys(fixtures)) == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

identity_gate = Matrix{Float64}(I, 4, 4)
swap_gate = Matrix(QET.swap_operator((2, 2)))
cnot_gate = [
    1.0 0 0 0
    0 1 0 0
    0 0 0 1
    0 0 1 0
]
controlled_z = Matrix(Diagonal([1.0, 1, 1, -1]))

native_results = (
    QET.is_entangling_gate(identity_gate),
    QET.is_entangling_gate(swap_gate),
    QET.is_entangling_gate(cnot_gate),
    QET.is_entangling_gate(controlled_z),
)
compatibility_results = (
    Compat.IsEntanglingGate(identity_gate),
    Compat.IsEntanglingGate(swap_gate),
    Compat.IsEntanglingGate(cnot_gate),
    Compat.IsEntanglingGate(controlled_z),
)

@testset "QETLAB IsEntanglingGate differential fixture" begin
    expected_flags = (
        fixtures["is_entangling_gate_identity_flag"][],
        fixtures["is_entangling_gate_swap_flag"][],
        fixtures["is_entangling_gate_cnot_flag"][],
        fixtures["is_entangling_gate_controlled_z_flag"][],
    )
    @test expected_flags == (0.0, 0.0, 1.0, 1.0)
    @test map(result -> result.status, native_results) ==
        (:not_entangling, :not_entangling, :entangling, :entangling)
    @test map(result -> result.status, compatibility_results) ==
        map(result -> result.status, native_results)

    cnot_witness = vec(fixtures["is_entangling_gate_cnot_witness"])
    @test count(!iszero, cnot_witness) <= 2
    @test norm(cnot_witness) ≈ 1
    @test QET.is_product_vector(cnot_witness, (2, 2)).status === :within_tolerance
    @test QET.is_product_vector(cnot_gate * cnot_witness, (2, 2)).status ===
        :outside_tolerance

    controlled_z_candidate = vec(fixtures["is_entangling_gate_controlled_z_last_candidate"])
    @test count(!iszero, controlled_z_candidate) <= 2
    @test norm(controlled_z_candidate) ≈ sqrt(2)
    @test QET.is_product_vector(controlled_z_candidate, (2, 2)).status === :within_tolerance
    @test QET.is_product_vector(controlled_z * controlled_z_candidate, (2, 2)).status ===
        :within_tolerance
    @test nnz(native_results[4].product_witness) <= 4
    @test QET.is_product_vector(
        controlled_z * native_results[4].product_witness, (2, 2)
    ).status === :outside_tolerance
end

println(
    "IsEntanglingGate oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
