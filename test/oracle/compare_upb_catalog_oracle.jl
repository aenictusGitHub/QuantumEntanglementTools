using JSON3
using LinearAlgebra
using Random
using SHA
using Test

module OracleUPBCatalog
using LinearAlgebra
using Random

include("../../src/dimensions.jl")
include("../../src/entanglement/upb_size.jl")
include("../../src/entanglement/upb_catalog.jl")
end

const QETUPB = OracleUPBCatalog
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const PINNED_UPB_SHA256 = "d60d9d2375d773c7e59930991ddf73427c9ffb0c77ded640b7f635f99ec13544"
const PINNED_ONE_FACTORIZATION_SHA256 = "7b042b3729b96dc4646c57adc24dcd5e4f59cf710520b0a6f98e7f3fd0c531ca"
const FIXTURE_NAMES = Set([
    "pyramid",
    "tiles",
    "generalized_tiles_1_d4",
    "generalized_tiles_2_3x4",
    "minimum_4x4",
    "quad_residue_d3",
    "six_parameter",
    "generalized_shifts_p5",
    "feng_2x2x2x2",
    "johnston_2_power_8",
    "chen_johnston_4x6",
    "feng_2x2x3",
    "feng_2x2x5",
    "feng_4x4",
    "feng_2x2x2x4",
    "feng_2x2x2x2x5",
    "dimension_dispatch_3x2",
])

function verify_upb_catalog_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function upb_catalog_fixture_factors(fixture)
    dimensions = Tuple(Int.(fixture.dimensions))
    cardinality = Int(fixture.cardinality)
    total_rows = sum(dimensions)
    Int(fixture.factor_rows) == total_rows ||
        error("fixture $(fixture.name) has inconsistent factor_rows")
    real_parts = Float64.(fixture.real)
    imaginary_parts = Float64.(fixture.imaginary)
    length(real_parts) == total_rows * cardinality ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == total_rows * cardinality ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    concatenated = reshape(values, total_rows, cardinality)
    offsets = cumsum((0, dimensions...))
    return ntuple(
        party -> concatenated[(offsets[party] + 1):offsets[party + 1], :],
        length(dimensions),
    )
end

function upb_catalog_phase_residual(actual, expected)
    size(actual) == size(expected) || return Inf
    residual = 0.0
    for column in axes(expected, 2)
        actual_column = @view actual[:, column]
        expected_column = @view expected[:, column]
        overlap = dot(actual_column, expected_column)
        iszero(overlap) && return Inf
        aligned = actual_column .* (overlap / abs(overlap))
        residual = max(residual, norm(aligned - expected_column))
    end
    return residual
end

function upb_catalog_fixture_construction(name::String)
    name == "pyramid" && return QETUPB.upb(:pyramid)
    name == "tiles" && return QETUPB.upb(:tiles)
    name == "generalized_tiles_1_d4" && return QETUPB.upb(:generalized_tiles_1, 4)
    name == "generalized_tiles_2_3x4" && return QETUPB.upb(:generalized_tiles_2, 3, 4)
    name == "minimum_4x4" && return QETUPB.upb(:minimum_4x4)
    name == "quad_residue_d3" && return QETUPB.upb(:quad_residue, 3)
    name == "six_parameter" &&
        return QETUPB.upb(:six_parameter, [0.31, 0.47, -0.2, 0.53, 0.61, 0.4])
    name == "generalized_shifts_p5" && return QETUPB.upb(:generalized_shifts, 5)
    name == "feng_2x2x2x2" && return QETUPB.upb(:feng_2x2x2x2)
    name == "johnston_2_power_8" && return QETUPB.upb(:johnston_2_power_8)
    name == "chen_johnston_4x6" && return QETUPB.upb(:chen_johnston_4x6)
    name == "feng_2x2x3" && return QETUPB.upb(:feng_2x2x3)
    name == "feng_2x2x5" && return QETUPB.upb(:feng_2x2x5)
    name == "feng_4x4" && return QETUPB.upb(:feng_4x4)
    name == "feng_2x2x2x4" && return QETUPB.upb(:feng_2x2x2x4)
    name == "feng_2x2x2x2x5" && return QETUPB.upb(:feng_2x2x2x2x5)
    name == "dimension_dispatch_3x2" && return QETUPB.upb((3, 2))
    return error("unreviewed UPB fixture: $name")
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "upb_catalog_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "upb_catalog_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "UPB catalog oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_upb_catalog_oracle.sh first",
)

digest = verify_upb_catalog_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-UPB-catalog" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.upb_sha256) == PINNED_UPB_SHA256 ||
    error("fixture does not identify the reviewed UPB.m source")
String(metadata.one_factorization_sha256) == PINNED_ONE_FACTORIZATION_SHA256 ||
    error("fixture does not identify the reviewed one_factorization.m source")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Int(metadata.upstream_defect_count) == 1 ||
    error("fixture does not retain the reviewed upstream defect")
occursin("independent local-vector phases", String(metadata.comparison_scope)) ||
    error("fixture does not record its gauge-aware comparison scope")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixture_names = Set(String(fixture.name) for fixture in payload.fixtures)
fixture_names == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

@testset "QETLAB UPB deterministic catalog differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected_factors = upb_catalog_fixture_factors(fixture)
        construction = upb_catalog_fixture_construction(name)
        @test String(fixture.comparison) == "local_vectors_up_to_independent_phase"
        @test construction.dimensions == Tuple(Int.(fixture.dimensions))
        @test construction.cardinality == Int(fixture.cardinality)
        @test length(construction.local_factors) == length(expected_factors)
        atol = Float64(fixture.atol)
        rtol = Float64(fixture.rtol)
        for party in eachindex(expected_factors)
            residual = upb_catalog_phase_residual(
                construction.local_factors[party], expected_factors[party]
            )
            @testset "$name party $party phase gauge" begin
                @test isapprox(residual, 0; atol, rtol)
            end
        end
        @test construction.orthogonality_residual < 2e-11
    end
end

@testset "QETLAB UPB pinned 4k reshape-order defect" begin
    defect = payload.upstream_defects
    @test String(defect.name) == "johnston_2_power_4k_reshape_order"
    @test String(defect.family) == "John2^4k"
    @test Int(defect.parameter) == 8
    @test Int(defect.cardinality) == 12
    @test Float64(defect.pinned_orthogonality_residual) > 0.1
    @test occursin("corrected", String(defect.native_disposition))
    corrected = QETUPB.upb(:johnston_2_power_4k, 8)
    @test corrected.cardinality == 12
    @test corrected.orthogonality_residual < 1e-12
end

println(
    "UPB catalog oracle comparison passed for ",
    length(payload.fixtures),
    " deterministic fixtures and one pinned-defect fixture from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
