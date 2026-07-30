using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SHA
using Test

const QETInducedSchattenOracle = QuantumEntanglementTools
const CompatInducedSchattenOracle = QuantumEntanglementTools.MATLABCompat
const INDUCED_SCHATTEN_PINNED_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const INDUCED_SCHATTEN_FIXTURE_NAMES = Set([
    "induced_schatten_identity_3to2_value",
    "induced_schatten_identity_3to2_witness",
    "induced_schatten_diagonal_3to2_value",
    "induced_schatten_diagonal_3to2_witness",
    "induced_schatten_amplitude_4to2_value",
    "induced_schatten_amplitude_4to2_witness",
    "induced_schatten_amplitude_2to2_value",
    "induced_schatten_amplitude_2to2_witness",
])

function induced_schatten_verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function induced_schatten_fixture_array(fixture)
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

committed_induced_schatten_fixture = joinpath(
    @__DIR__, "fixtures", "induced_schatten_octave_11_3_qetlab_d858961.json"
)
induced_schatten_fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "induced_schatten_oracle.json")
    isfile(generated) ? generated : committed_induced_schatten_fixture
else
    abspath(ARGS[1])
end
isfile(induced_schatten_fixture_path) || error(
    "InducedSchattenNorm oracle fixture not found at " *
    "$induced_schatten_fixture_path; run " *
    "scripts/matlab_oracle/run_induced_schatten_oracle.sh first",
)

induced_schatten_digest = induced_schatten_verify_digest(induced_schatten_fixture_path)
induced_schatten_payload = JSON3.read(read(induced_schatten_fixture_path, String))
induced_schatten_metadata = induced_schatten_payload.metadata
String(induced_schatten_metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(induced_schatten_metadata.schema)")
String(induced_schatten_metadata.tier) == "WP2-induced-schatten" ||
    error("unexpected fixture tier: $(induced_schatten_metadata.tier)")
String(induced_schatten_metadata.qetlab_commit) == INDUCED_SCHATTEN_PINNED_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(induced_schatten_metadata.source_free_fixture) ||
    error("fixture is not marked source-free")
Bool(induced_schatten_metadata.seeded_random_start) ||
    error("fixture does not record its seeded random-start policy")
occursin("unreachable", String(induced_schatten_metadata.evidence_scope)) ||
    error("fixture does not record the pinned TOL/X0 forwarding defect")
occursin("TooManyArguments", String(induced_schatten_metadata.tol_x0_error_identifier)) ||
    error("fixture does not retain the pinned TOL/X0 failure identifier")
Int(induced_schatten_metadata.fixture_count) == length(induced_schatten_payload.fixtures) ||
    error("fixture count does not match metadata")

if normpath(induced_schatten_fixture_path) == normpath(committed_induced_schatten_fixture)
    String(induced_schatten_metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(induced_schatten_metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

induced_schatten_fixtures = Dict(
    String(fixture.name) => induced_schatten_fixture_array(fixture) for
    fixture in induced_schatten_payload.fixtures
)
Set(keys(induced_schatten_fixtures)) == INDUCED_SCHATTEN_FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

induced_schatten_identity = QETInducedSchattenOracle.KrausRepresentation([
    Matrix{Float64}(I, 2, 2)
])
induced_schatten_diagonal = QETInducedSchattenOracle.KrausRepresentation([
    Diagonal([1.0, 2.0])
])
induced_schatten_amplitude_operators = [
    Diagonal([1.0, sqrt(0.7)]),
    [
        0.0 sqrt(0.3)
        0.0 0.0
    ],
]
induced_schatten_amplitude = QETInducedSchattenOracle.KrausRepresentation(
    induced_schatten_amplitude_operators
)

function induced_schatten_fixture_value(name)
    return only(induced_schatten_fixtures[name])
end

function induced_schatten_witness_ratio(map, witness, p, q)
    input_norm = QETInducedSchattenOracle.schatten_norm(witness, p)
    output_norm = QETInducedSchattenOracle.schatten_norm(
        QETInducedSchattenOracle.apply_channel(witness, map), q
    )
    return output_norm / input_norm
end

@testset "QETLAB InducedSchattenNorm differential fixture" begin
    cases = (
        (
            "identity_3to2",
            induced_schatten_identity,
            3,
            2,
            "induced_schatten_identity_3to2_value",
            "induced_schatten_identity_3to2_witness",
        ),
        (
            "diagonal_3to2",
            induced_schatten_diagonal,
            3,
            2,
            "induced_schatten_diagonal_3to2_value",
            "induced_schatten_diagonal_3to2_witness",
        ),
        (
            "amplitude_4to2",
            induced_schatten_amplitude,
            4,
            2,
            "induced_schatten_amplitude_4to2_value",
            "induced_schatten_amplitude_4to2_witness",
        ),
    )

    for (_, map, p, q, value_name, witness_name) in cases
        expected = induced_schatten_fixture_value(value_name)
        witness = induced_schatten_fixtures[witness_name]
        @test QETInducedSchattenOracle.schatten_norm(witness, p) ≈ 1 atol = 3e-7
        @test induced_schatten_witness_ratio(map, witness, p, q) ≈ expected atol = 3e-10 rtol =
            3e-10

        native = QETInducedSchattenOracle.induced_schatten_lower_bound(
            MersenneTwister(1),
            map,
            p;
            q=q,
            initial_matrix=witness,
            tolerance=sqrt(eps(Float64)),
            max_iterations=1_000,
        )
        @test native.bound_kind === :lower_bound
        @test native.status === :converged_lower_bound
        @test native.value + 3e-10 >= expected
        @test native.value ≈ induced_schatten_witness_ratio(map, native.witness, p, q) atol =
            3e-10 rtol = 3e-10
    end

    identity_value = induced_schatten_fixture_value("induced_schatten_identity_3to2_value")
    @test identity_value ≈ 2.0^(1 / 6) atol = 3e-9

    exact_expected = induced_schatten_fixture_value("induced_schatten_amplitude_2to2_value")
    exact_witness = induced_schatten_fixtures["induced_schatten_amplitude_2to2_witness"]
    exact_native = QETInducedSchattenOracle.induced_schatten_lower_bound(
        MersenneTwister(2), induced_schatten_amplitude, 2
    )
    @test exact_native.status === :exact
    @test exact_native.value ≈ exact_expected atol = 3e-14 rtol = 3e-14
    @test induced_schatten_witness_ratio(induced_schatten_amplitude, exact_witness, 2, 2) ≈
        exact_expected atol = 3e-14 rtol = 3e-14

    compat_exact = CompatInducedSchattenOracle.InducedSchattenNorm(
        MersenneTwister(3), induced_schatten_amplitude_operators, "fro"
    )
    @test compat_exact.value == exact_native.value
    @test compat_exact.bound_kind === :exact
end

println(
    "InducedSchattenNorm oracle comparison passed for ",
    length(induced_schatten_payload.fixtures),
    " fixtures from ",
    induced_schatten_metadata.engine,
    " ",
    induced_schatten_metadata.engine_version,
    " (SHA-256 ",
    induced_schatten_digest,
    ").",
)
