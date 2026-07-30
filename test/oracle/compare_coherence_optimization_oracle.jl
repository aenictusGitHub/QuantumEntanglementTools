using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const CoherenceOptimizationOracle = QuantumEntanglementTools
if !isdefined(CoherenceOptimizationOracle, :is_k_incoherent)
    Base.include(
        CoherenceOptimizationOracle,
        joinpath(@__DIR__, "..", "..", "src", "coherence", "coherence_optimization.jl"),
    )
end

const COHERENCE_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const COHERENCE_SOURCE_HASHES = (
    isk="d924a12b239a1ea9a784ed3252ff0f289f2281aeae6022b645b169c4d2779f6a",
    absolute="6bf68ba8b87e441dd6b2ab2b718f621e2fa6a8a67a92f016b3629266496658ef",
    robustness="b8af3c8f056ba49ad703e361b1b6f8554776747b0790efefed34a359a28b7986",
    trace="f74de030eb833da8c1bd2352961e1080c652746fa11572812000f50a0cd2048b",
    generalized="aa3d429fc0fb9f715391620de67d9407f638bd58307922e49bf7ab5b57e37e78",
    band="d4b6b6be5793515b3b98f98ba585f04041e6e7875f22437a36b413f6308be8d7",
    robk="99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5",
)

function coherence_verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function coherence_decode_array(encoded)
    dimensions = Tuple(Int.(encoded.dims))
    values = complex.(Float64.(encoded.real), Float64.(encoded.imaginary))
    return reshape(values, dimensions)
end

function coherence_decode_state(encoded)
    array = coherence_decode_array(encoded)
    return ndims(array) == 2 && min(size(array)...) == 1 ? vec(array) : array
end

committed_coherence_fixture = joinpath(
    @__DIR__, "fixtures", "coherence_optimization_octave_11_3_qetlab_d858961.json"
)
coherence_fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "coherence_optimization_oracle.json")
    isfile(generated) ? generated : committed_coherence_fixture
else
    abspath(ARGS[1])
end
isfile(coherence_fixture_path) || error(
    "coherence optimization oracle fixture not found at $coherence_fixture_path; " *
    "run scripts/matlab_oracle/run_coherence_optimization_oracle.sh first",
)

coherence_digest = coherence_verify_digest(coherence_fixture_path)
coherence_payload = JSON3.read(read(coherence_fixture_path, String))
coherence_metadata = coherence_payload.metadata
String(coherence_metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(coherence_metadata.schema)")
String(coherence_metadata.tier) == "WP4-coherence-optimization" ||
    error("unexpected oracle tier: $(coherence_metadata.tier)")
String(coherence_metadata.qetlab_commit) == COHERENCE_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(coherence_metadata.source_free_fixture) || error("fixture is not marked source-free")
Bool(coherence_metadata.solver_free_branches_only) ||
    error("fixture unexpectedly records solver-dependent QETLAB outputs")
String(coherence_metadata.isk_source_sha256) == COHERENCE_SOURCE_HASHES.isk ||
    error("fixture records the wrong IskIncoherent.m digest")
String(coherence_metadata.is_absk_source_sha256) == COHERENCE_SOURCE_HASHES.absolute ||
    error("fixture records the wrong IsAbskIncoh.m digest")
String(coherence_metadata.robustness_source_sha256) == COHERENCE_SOURCE_HASHES.robustness ||
    error("fixture records the wrong RobustnessCoherence.m digest")
String(coherence_metadata.trace_distance_source_sha256) == COHERENCE_SOURCE_HASHES.trace ||
    error("fixture records the wrong TraceDistanceCoherence.m digest")
String(coherence_metadata.generalized_source_sha256) ==
COHERENCE_SOURCE_HASHES.generalized ||
    error("fixture records the wrong GenRobustnesskCoherence.m digest")
String(coherence_metadata.band_helper_source_sha256) == COHERENCE_SOURCE_HASHES.band ||
    error("fixture records the wrong has_band_k_ordering.m digest")
String(coherence_metadata.robk_source_sha256) == COHERENCE_SOURCE_HASHES.robk ||
    error("fixture records the wrong RobkCohValue.m digest")
Int(coherence_metadata.isk_fixture_count) == length(coherence_payload.isk_fixtures) ||
    error("k-incoherence fixture count does not match metadata")
Int(coherence_metadata.absolute_fixture_count) ==
length(coherence_payload.absolute_fixtures) ||
    error("absolute fixture count does not match metadata")
Int(coherence_metadata.robustness_fixture_count) ==
length(coherence_payload.robustness_fixtures) ||
    error("robustness fixture count does not match metadata")
Int(coherence_metadata.trace_fixture_count) == length(coherence_payload.trace_fixtures) ||
    error("trace fixture count does not match metadata")
Int(coherence_metadata.generalized_fixture_count) ==
length(coherence_payload.generalized_fixtures) ||
    error("generalized fixture count does not match metadata")

if normpath(coherence_fixture_path) == normpath(committed_coherence_fixture)
    String(coherence_metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(coherence_metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB coherence optimization analytic fixtures" begin
    for fixture in coherence_payload.isk_fixtures
        result = CoherenceOptimizationOracle.is_k_incoherent(
            coherence_decode_state(fixture.state), Int(fixture.k)
        )
        expected = Int(fixture.qetlab_value)
        @test expected in (0, 1)
        @test result.verdict === isone(expected)
    end

    for fixture in coherence_payload.absolute_fixtures
        result = CoherenceOptimizationOracle.is_absolutely_k_incoherent(
            coherence_decode_state(fixture.state), Int(fixture.k)
        )
        expected = Int(fixture.qetlab_value)
        @test expected in (-1, 0, 1)
        if expected == -1
            @test result.verdict === nothing
            @test result.status === :unknown
        else
            @test result.verdict === isone(expected)
        end
    end

    for fixture in coherence_payload.robustness_fixtures
        result = CoherenceOptimizationOracle.robustness_coherence(
            coherence_decode_state(fixture.state)
        )
        @test result.value ≈ Float64(fixture.qetlab_value) atol = Float64(fixture.atol) rtol = Float64(
            fixture.rtol
        )
        @test result.free_state !== nothing
        if iszero(result.value)
            @test result.noise_state === nothing
        end
    end

    for fixture in coherence_payload.trace_fixtures
        result = CoherenceOptimizationOracle.trace_distance_coherence(
            coherence_decode_state(fixture.state)
        )
        @test result.value ≈ Float64(fixture.qetlab_value) atol = Float64(fixture.atol) rtol = Float64(
            fixture.rtol
        )
        @test real.(diag(result.free_state)) ≈ Float64.(fixture.closest_diagonal) atol = Float64(
            fixture.atol
        ) rtol = Float64(fixture.rtol)
    end

    for fixture in coherence_payload.generalized_fixtures
        result = CoherenceOptimizationOracle.generalized_robustness_k_coherence(
            coherence_decode_state(fixture.state), Int(fixture.k)
        )
        @test result.value ≈ Float64(fixture.qetlab_value) atol = Float64(fixture.atol) rtol = Float64(
            fixture.rtol
        )
        @test result.diagnostics.branch_index == Int(fixture.qetlab_branch_index)
        @test result.decomposition !== nothing
    end

    defect = coherence_payload.band_defect
    @test !Bool(defect.pinned_has_ordering)
    @test Bool(defect.bruteforce_has_ordering)
    @test isempty(defect.pinned_ordering)
    matrix = coherence_decode_array(defect.matrix)
    corrected = CoherenceOptimizationOracle._cohopt_band_ordering(
        matrix, Int(defect.k); max_search_nodes=100
    )
    @test corrected.status === :found
    adjacency = .!iszero.(matrix)
    @test CoherenceOptimizationOracle._cohopt_bandwidth(adjacency, corrected.ordering) <=
        Int(defect.k)
    valid_ordering = Int.(defect.valid_ordering)
    @test CoherenceOptimizationOracle._cohopt_bandwidth(adjacency, valid_ordering) <=
        Int(defect.k)
end

println(
    "Coherence optimization oracle comparison passed for ",
    length(coherence_payload.isk_fixtures),
    " k-incoherence, ",
    length(coherence_payload.absolute_fixtures),
    " absolute, ",
    length(coherence_payload.robustness_fixtures),
    " robustness, ",
    length(coherence_payload.trace_fixtures),
    " trace, and ",
    length(coherence_payload.generalized_fixtures),
    " pure generalized fixtures from ",
    coherence_metadata.engine,
    " ",
    coherence_metadata.engine_version,
    " (SHA-256 ",
    coherence_digest,
    ").",
)
