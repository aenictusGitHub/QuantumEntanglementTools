using JSON3
using QuantumEntanglementTools
using Random
using SHA
using Test

const CopCliqueOracleQET = QuantumEntanglementTools

if !isdefined(CopCliqueOracleQET, :CopositivityResult)
    Base.include(
        CopCliqueOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "optimization", "copositivity_clique.jl"),
    )
end

const COPCLIQUE_PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const COPCLIQUE_PINNED_IS_COPOSITIVE_SHA = "6637ea6c9d4e2f7647520087cb16b2608d5e32587bc3f87e62340ca102a93d22"
const COPCLIQUE_PINNED_CLIQUE_NUMBER_SHA = "dc99fa9c91a6e2be985881a3598013347acf63b42fb59303b48d8e0aac6a1eff"
const COPCLIQUE_PINNED_POLYNOMIAL_SHA = "4814f9bfbfe4379b725035c1a6c9d60652216278b20e74bdc745c148ca0544bb"
const COPCLIQUE_PINNED_POLYNOMIAL_OPTIMIZE_SHA = "b70a02cc12c9b341a0aaf11e22a4e0382d3e5c18f9b345ff7e08791ad59a24d7"
const COPCLIQUE_COPOSITIVITY_NAMES = Set([
    "positive_diagonal",
    "negative_diagonal",
    "fixed_threshold_negative",
    "negative_pair_ray",
    "horn_boundary",
])
const COPCLIQUE_CLIQUE_NAMES = Set([
    "edgeless_three", "complete_three", "cycle_five", "triangle_plus_isolate", "path_four"
])

function verify_copclique_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function copclique_fixture_matrix(fixture, field::Symbol; integer::Bool)
    dimension = Int(fixture.dimension)
    values = getproperty(fixture, field)
    length(values) == dimension^2 ||
        error("fixture $(fixture.name) has an inconsistent matrix length")
    converted = integer ? Int.(values) : Float64.(values)
    return reshape(converted, dimension, dimension)
end

committed_copclique_fixture = joinpath(
    @__DIR__, "fixtures", "copositivity_clique_octave_11_3_qetlab_d858961.json"
)
copclique_fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "copositivity_clique_oracle.json")
    isfile(generated) ? generated : committed_copclique_fixture
else
    abspath(ARGS[1])
end
isfile(copclique_fixture_path) || error(
    "copositivity/clique oracle fixture not found at $copclique_fixture_path; " *
    "run scripts/matlab_oracle/run_copositivity_clique_oracle.sh first",
)

copclique_digest = verify_copclique_digest(copclique_fixture_path)
copclique_payload = JSON3.read(read(copclique_fixture_path, String))
copclique_metadata = copclique_payload.metadata
String(copclique_metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(copclique_metadata.schema)")
String(copclique_metadata.tier) == "WP5-copositivity-clique" ||
    error("unexpected fixture tier: $(copclique_metadata.tier)")
String(copclique_metadata.qetlab_commit) == COPCLIQUE_PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(copclique_metadata.source_free_fixture) || error("fixture is not marked source-free")
String(copclique_metadata.is_copositive_sha256) == COPCLIQUE_PINNED_IS_COPOSITIVE_SHA ||
    error("fixture records the wrong IsCopositive.m hash")
String(copclique_metadata.clique_number_sha256) == COPCLIQUE_PINNED_CLIQUE_NUMBER_SHA ||
    error("fixture records the wrong CliqueNumber.m hash")
String(copclique_metadata.copositive_polynomial_sha256) ==
COPCLIQUE_PINNED_POLYNOMIAL_SHA ||
    error("fixture records the wrong CopositivePolynomial.m hash")
String(copclique_metadata.polynomial_optimize_sha256) ==
COPCLIQUE_PINNED_POLYNOMIAL_OPTIMIZE_SHA ||
    error("fixture records the wrong PolynomialOptimize.m hash")
Int(copclique_metadata.copositivity_fixture_count) ==
length(copclique_payload.copositivity_fixtures) ||
    error("copositivity fixture count does not match metadata")
Int(copclique_metadata.clique_fixture_count) == length(copclique_payload.clique_fixtures) ||
    error("clique fixture count does not match metadata")
Set(String(fixture.name) for fixture in copclique_payload.copositivity_fixtures) ==
COPCLIQUE_COPOSITIVITY_NAMES ||
    error("copositivity fixture names do not match the reviewed set")
Set(String(fixture.name) for fixture in copclique_payload.clique_fixtures) ==
COPCLIQUE_CLIQUE_NAMES || error("clique fixture names do not match the reviewed set")
occursin("-1e-9", String(copclique_metadata.fixed_threshold_defect)) ||
    error("the pinned fixed-threshold discrepancy is not recorded")
occursin("elapsed wall time", String(copclique_metadata.evidence_scope)) ||
    error("the pinned randomized-oracle limitation is not recorded")

if normpath(copclique_fixture_path) == normpath(committed_copclique_fixture)
    String(copclique_metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(copclique_metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB copositivity and clique-number source-free fixture" begin
    for fixture in copclique_payload.copositivity_fixtures
        name = String(fixture.name)
        matrix = copclique_fixture_matrix(fixture, :matrix; integer=false)
        expected = Int(fixture.expected)
        result = CopCliqueOracleQET.copositivity_criterion(
            MersenneTwister(100 + Int(fixture.dimension)),
            matrix;
            hierarchy=:nosdp,
            allow_densify=true,
            inner_samples=0,
        )
        @test result.samples_evaluated == 0
        if name == "positive_diagonal"
            @test expected == 1
            @test result.verdict === true
            @test result.status === CopCliqueOracleQET.CopositivityCertifiedTrue
            @test result.certificate_kind === :entrywise_nonnegative
        elseif name in ("negative_diagonal", "negative_pair_ray")
            @test expected == 0
            @test result.verdict === false
            @test result.status === CopCliqueOracleQET.CopositivityCertifiedFalse
            @test result.witness.value < 0
        elseif name == "fixed_threshold_negative"
            @test expected == 1
            @test result.verdict === nothing
            @test result.status === CopCliqueOracleQET.CopositivityNumericalBoundary
            @test result.tolerance > abs(matrix[1, 1])
        else
            @test name == "horn_boundary"
            @test expected == -1
            @test result.verdict === nothing
            @test result.status === CopCliqueOracleQET.CopositivityHierarchyUnknown
            @test result.lower_bound < 0
            @test result.upper_bound === nothing
        end
    end

    oracle_rng = MersenneTwister(901)
    untouched_oracle_rng = copy(oracle_rng)
    for fixture in copclique_payload.clique_fixtures
        name = String(fixture.name)
        adjacency = copclique_fixture_matrix(fixture, :adjacency; integer=true)
        expected_upper = Int(fixture.expected_upper)
        expected_lower = Int(fixture.expected_lower)
        result = CopCliqueOracleQET.clique_number_bounds(
            oracle_rng, adjacency; hierarchy=:nosdp, allow_densify=true, inner_samples=0
        )
        @test result.bounds_certified
        @test result.samples_evaluated == 0
        @test result.lower_bound >= expected_lower
        @test result.upper_bound <= expected_upper
        @test result.lower_bound <= result.upper_bound
        @test length(result.best_clique) <= result.lower_bound
        @test all(
            adjacency[left, right] == 1 for
            (position, left) in enumerate(result.best_clique) for
            right in result.best_clique[(position + 1):end]
        )
        if name in ("edgeless_three", "complete_three", "triangle_plus_isolate")
            @test result.lower_bound == expected_lower
            @test result.upper_bound == expected_upper
            @test result.exact
        elseif name == "cycle_five"
            @test (expected_lower, expected_upper) == (2, 3)
            @test (result.lower_bound, result.upper_bound) == (2, 3)
            @test !result.exact
        else
            @test name == "path_four"
            @test (expected_lower, expected_upper) == (2, 3)
            @test (result.lower_bound, result.upper_bound) == (2, 2)
            @test result.exact
            @test :greedy_coloring in result.upper_certificate_kinds
        end
    end
    @test rand(oracle_rng) == rand(untouched_oracle_rng)
end

println(
    "Copositivity/clique oracle comparison passed for ",
    length(copclique_payload.copositivity_fixtures),
    " copositivity and ",
    length(copclique_payload.clique_fixtures),
    " clique fixtures from ",
    copclique_metadata.engine,
    " ",
    copclique_metadata.engine_version,
    " (SHA-256 ",
    copclique_digest,
    ").",
)
