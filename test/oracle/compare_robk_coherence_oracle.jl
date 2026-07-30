using JSON3
using QuantumEntanglementTools
using SHA
using Test

const RobkOracle = QuantumEntanglementTools
const RobkCompatOracle = RobkOracle.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const PINNED_SOURCE_SHA256 = "99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5"

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function fixture_state(fixture)
    real_values = Float64.(fixture.state_real)
    imaginary_values = Float64.(fixture.state_imaginary)
    return complex.(real_values, imaginary_values)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "robk_coherence_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "robk_coherence_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "RobkCohValue oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_robk_coherence_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-robk-coherence" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.qetlab_source_sha256) == PINNED_SOURCE_SHA256 ||
    error("fixture records the wrong RobkCohValue.m digest")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Int(metadata.correction_case_count) == length(payload.correction_cases) ||
    error("correction-case count does not match metadata")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB RobkCohValue differential fixture" begin
    for fixture in payload.fixtures
        state = fixture_state(fixture)
        result = RobkOracle.pure_k_coherence_robustness(
            state, Int(fixture.k); atol=Float64(fixture.atol), rtol=Float64(fixture.rtol)
        )
        @test String(fixture.comparison) == "approximate"
        @test iszero(Float64(fixture.robustness_imaginary))
        @test result.value ≈ Float64(fixture.robustness_real) atol = Float64(fixture.atol) rtol = Float64(
            fixture.rtol
        )
        @test result.branch_index == Int(fixture.branch_index)
        compatibility_value, compatibility_branch = RobkCompatOracle.RobkCohValue(
            state, Int(fixture.k); atol=Float64(fixture.atol), rtol=Float64(fixture.rtol)
        )
        @test compatibility_value ≈ Float64(fixture.robustness_real) atol = Float64(
            fixture.atol
        ) rtol = Float64(fixture.rtol)
        @test compatibility_branch == Int(fixture.branch_index)
    end

    corrections = Dict(
        String(fixture.name) => fixture for fixture in payload.correction_cases
    )
    @test Set(keys(corrections)) ==
        Set(("unsorted_coefficients", "complex_phases", "nonnormalized_coefficients"))

    unsorted = corrections["unsorted_coefficients"]
    unsorted_result = RobkOracle.pure_k_coherence_robustness(
        fixture_state(unsorted), Int(unsorted.k)
    )
    canonical_result = RobkOracle.pure_k_coherence_robustness(
        Float64.(unsorted.canonical_magnitudes), Int(unsorted.k)
    )
    @test unsorted_result.value ≈ canonical_result.value
    @test unsorted_result.branch_index == canonical_result.branch_index
    @test RobkCompatOracle.RobkCohValue(fixture_state(unsorted), Int(unsorted.k)) ==
        (unsorted_result.value, unsorted_result.branch_index)
    @test !isapprox(
        unsorted_result.value,
        Float64(unsorted.qetlab_robustness_real);
        atol=1e-12,
        rtol=1e-12,
    )

    complex_case = corrections["complex_phases"]
    complex_result = RobkOracle.pure_k_coherence_robustness(
        fixture_state(complex_case), Int(complex_case.k)
    )
    complex_canonical = RobkOracle.pure_k_coherence_robustness(
        Float64.(complex_case.canonical_magnitudes), Int(complex_case.k)
    )
    @test complex_result.value ≈ complex_canonical.value
    @test RobkCompatOracle.RobkCohValue(fixture_state(complex_case), Int(complex_case.k)) ==
        (complex_result.value, complex_result.branch_index)
    @test !iszero(Float64(complex_case.qetlab_robustness_imaginary))

    nonnormalized = corrections["nonnormalized_coefficients"]
    @test_throws ArgumentError RobkOracle.pure_k_coherence_robustness(
        fixture_state(nonnormalized), Int(nonnormalized.k)
    )
    @test_throws ArgumentError RobkCompatOracle.RobkCohValue(
        fixture_state(nonnormalized), Int(nonnormalized.k)
    )
end

println(
    "RobkCohValue oracle comparison passed for ",
    length(payload.fixtures),
    " theorem-domain fixtures and ",
    length(payload.correction_cases),
    " reviewed correction cases from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
