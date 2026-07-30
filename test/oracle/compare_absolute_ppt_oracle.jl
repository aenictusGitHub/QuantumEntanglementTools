using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const AbsPPTOracleQET = QuantumEntanglementTools
const AbsPPTOracleCompat = AbsPPTOracleQET.MATLABCompat

if !isdefined(AbsPPTOracleQET, :abs_ppt_constraints)
    Base.include(
        AbsPPTOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "entanglement", "absolute_ppt.jl"),
    )
end

const ABS_PPT_PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const ABS_PPT_CONSTRAINTS_SOURCE_SHA = "e7cdc38ce5b16a68165f0b983a6033496147fd8a6c1df4a33743c15bd11a5b75"
const IS_ABS_PPT_SOURCE_SHA = "f6272ff80fa53789ea92485f32b9155bf55327670d8e24c58f51f59f1286a5d2"
const ABS_PPT_FIXTURE_NAMES = Set([
    "p2_constraint_1",
    "p3_constraint_1",
    "p3_constraint_2",
    "p4_limited_constraint_1",
    "p4_limited_constraint_2",
    "p4_limited_constraint_3",
    "p4_early_violation",
])

function _abs_ppt_verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function _abs_ppt_fixture_matrix(fixture)
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

committed_fixture = joinpath(
    @__DIR__, "fixtures", "absolute_ppt_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "absolute_ppt_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "absolute-PPT oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_absolute_ppt_oracle.sh first",
)

digest = _abs_ppt_verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP3-absolute-ppt" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == ABS_PPT_PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.abs_constraints_source_sha256) == ABS_PPT_CONSTRAINTS_SOURCE_SHA ||
    error("fixture does not record the pinned AbsPPTConstraints.m SHA-256")
String(metadata.is_abs_ppt_source_sha256) == IS_ABS_PPT_SOURCE_SHA ||
    error("fixture does not record the pinned IsAbsPPT.m SHA-256")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Int(metadata.p2_constraint_count) == 1 ||
    error("fixture does not retain the pinned p=2 count")
Int(metadata.p3_constraint_count) == 2 ||
    error("fixture does not retain the pinned p=3 count")
Int(metadata.p4_full_constraint_count) == 10 ||
    error("fixture does not retain the pinned p=4 count")
Int(metadata.p4_limited_constraint_count) == 3 ||
    error("fixture does not retain LIM behavior")
Int(metadata.p4_early_constraint_count) == 1 ||
    error("fixture does not retain ESC_IF_NPOS behavior")
occursin("while pinned", lowercase(String(metadata.scalar_dim_inconsistency))) ||
    error("fixture does not record the pinned scalar-DIM inconsistency")
occursin("2612", String(metadata.p6_redundant_constraint_note)) ||
    error("fixture does not record the pinned p=6 redundant-family evidence")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(
    String(fixture.name) => _abs_ppt_fixture_matrix(fixture) for fixture in payload.fixtures
)
Set(keys(fixtures)) == ABS_PPT_FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

@testset "QETLAB absolute-PPT source-free fixture" begin
    family2 = AbsPPTOracleQET.abs_ppt_constraints([0.4, 0.3, 0.2, 0.1]; dims=(2, 2))
    @test family2.constraints[1].matrix ≈ fixtures["p2_constraint_1"] atol = 8e-14 rtol =
        8e-14

    family3 = AbsPPTOracleQET.abs_ppt_constraints(
        [0.30, 0.19, 0.14, 0.11, 0.09, 0.07, 0.05, 0.03, 0.02]; dims=(3, 3)
    )
    for index in 1:2
        @test family3.constraints[index].matrix ≈ fixtures["p3_constraint_$index"] atol =
            8e-14 rtol = 8e-14
    end

    spectrum4 = collect((16:-1:1) ./ sum(1:16))
    limited4 = AbsPPTOracleQET.abs_ppt_constraints(
        spectrum4; dims=(4, 4), max_constraints=3
    )
    @test limited4.status === AbsPPTOracleQET.AbsPPTEnumerationConstraintLimit
    for index in 1:3
        @test limited4.constraints[index].matrix ≈ fixtures["p4_limited_constraint_$index"] atol =
            8e-14 rtol = 8e-14
    end

    early4 = AbsPPTOracleQET.abs_ppt_constraints(
        [1.0; zeros(15)]; dims=(4, 4), stop_on_violation=true, max_constraints=nothing
    )
    @test early4.status === AbsPPTOracleQET.AbsPPTEnumerationEarlyViolation
    @test length(early4.constraints) == 1
    @test early4.constraints[1].matrix ≈ fixtures["p4_early_violation"] atol = 8e-14 rtol =
        8e-14

    outputs = payload.is_abs_ppt_outputs
    @test Int(outputs.maximally_mixed_2x2) == 1
    @test AbsPPTOracleQET.is_abs_ppt(fill(0.25, 4); dims=(2, 2)).verdict === true
    @test Int(outputs.exhaustive_yes_2x2) == 1
    @test AbsPPTOracleQET.is_abs_ppt([0.45, 0.35, 0.1, 0.1]; dims=(2, 2)).status ===
        AbsPPTOracleQET.AbsolutePPTExhaustiveCertified
    @test Int(outputs.pure_no_2x2) == 0
    @test AbsPPTOracleQET.is_abs_ppt([1.0, 0.0, 0.0, 0.0]; dims=(2, 2)).verdict === false
    @test Int(outputs.boundary_2x2) == 1
    @test AbsPPTOracleQET.is_abs_ppt([0.5, 1 / 6, 1 / 6, 1 / 6]; dims=(2, 2)).status ===
        AbsPPTOracleQET.AbsolutePPTNumericalBoundary
    @test Int(outputs.default_rectangular_length6) == 1
    @test AbsPPTOracleQET.is_abs_ppt(fill(1 / 6, 6); dims=(2, 3)).verdict === true
    @test_throws ArgumentError AbsPPTOracleQET.is_abs_ppt(fill(1 / 6, 6))

    if isdefined(AbsPPTOracleCompat, :AbsPPTConstraints)
        compat2 = AbsPPTOracleCompat.AbsPPTConstraints([0.4, 0.3, 0.2, 0.1], [2, 2])
        @test compat2[1] ≈ fixtures["p2_constraint_1"] atol = 8e-14 rtol = 8e-14
        compat_limited = AbsPPTOracleCompat.AbsPPTConstraints(spectrum4, [4, 4], 0, 3)
        @test length(compat_limited) == 3
        @test all(
            isapprox(
                compat_limited[index],
                fixtures["p4_limited_constraint_$index"];
                atol=8e-14,
                rtol=8e-14,
            ) for index in 1:3
        )
    end
    if isdefined(AbsPPTOracleCompat, :IsAbsPPT)
        @test AbsPPTOracleCompat.IsAbsPPT(fill(0.25, 4), [2, 2]) == 1
        @test AbsPPTOracleCompat.IsAbsPPT([1.0, 0.0, 0.0, 0.0], [2, 2]) == 0
        @test AbsPPTOracleCompat.IsAbsPPT(fill(1 / 6, 6)) == 1
        # The native numerical-boundary correction is intentionally tri-state.
        @test AbsPPTOracleCompat.IsAbsPPT([0.5, 1 / 6, 1 / 6, 1 / 6], [2, 2]) == -1
    end
end

println(
    "Absolute-PPT oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
