using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SHA
using Test

const SymExtOracleQET = QuantumEntanglementTools
const SYMEXT_PINNED_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const SYMEXT_SOURCE_SHA = "ca8e1aaf9b3766a3e29bc1375eecf473cdb77bc71159a033411c77eacb6e7e59"
const SYMINNER_SOURCE_SHA = "f73e38499197fb18c1ff9cff554fd77e13ab7ccfbe48590a25437cd70d7f39da"
const RANDOM_PPT_SOURCE_SHA = "e5208eac2c04dac95051d8c5850bb6fdaaca6430117416d0675b2c9a26ef04fd"
const JACOBI_SOURCE_SHA = "3001148f6fb136306a23a48ce6176b2b2ed65b278e4f1f998cb180486124d513"

if !isdefined(SymExtOracleQET, :SymmetricExtensionStatus)
    Base.include(
        SymExtOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "entanglement", "symmetric_extensions.jl"),
    )
end

function symext_oracle_digest(path)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function symext_oracle_matrix(name)
    name == "mixed" && return Matrix{Float64}(I, 4, 4) / 4
    if name == "bell"
        vector = [inv(sqrt(2.0)), 0, 0, inv(sqrt(2.0))]
        return vector * adjoint(vector)
    end
    return error("unreviewed fixture matrix $name")
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "symmetric_extensions_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "symmetric_extensions_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "symmetric-extension oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_symmetric_extensions_oracle.sh first",
)

digest = symext_oracle_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema")
String(metadata.tier) == "WP3-symmetric-extensions-random-ppt" ||
    error("unexpected fixture tier")
String(metadata.qetlab_commit) == SYMEXT_PINNED_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.symmetric_extension_source_sha256) == SYMEXT_SOURCE_SHA ||
    error("fixture has the wrong SymmetricExtension.m hash")
String(metadata.symmetric_inner_extension_source_sha256) == SYMINNER_SOURCE_SHA ||
    error("fixture has the wrong SymmetricInnerExtension.m hash")
String(metadata.random_ppt_state_source_sha256) == RANDOM_PPT_SOURCE_SHA ||
    error("fixture has the wrong RandomPPTState.m hash")
String(metadata.jacobi_poly_source_sha256) == JACOBI_SOURCE_SHA ||
    error("fixture has the wrong jacobi_poly.m hash")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.extension_fixture_count) == length(payload.extension_fixtures) ||
    error("extension fixture count mismatch")
Int(metadata.jacobi_fixture_count) == length(payload.jacobi_fixtures) ||
    error("Jacobi fixture count mismatch")
occursin("solver-free", String(metadata.solver_scope)) ||
    error("fixture does not delimit its solver scope")
occursin("not compared", String(metadata.random_scope)) ||
    error("fixture does not state its random-stream comparison boundary")
occursin("not automatically", String(metadata.inner_warning)) ||
    error("fixture omits the inner dual warning")
occursin("bounded", String(metadata.native_low_rank_deviation)) ||
    error("fixture omits the bounded low-rank deviation")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" || error("the committed fixture must record Octave")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB symmetric-extension solver-free fixture" begin
    @test length(payload.extension_fixtures) == 5
    for fixture in payload.extension_fixtures
        state = symext_oracle_matrix(String(fixture.matrix_name))
        result = SymExtOracleQET.symmetric_extension(
            state;
            order=Int(fixture.order),
            dims=Tuple(Int.(fixture.dimensions)),
            ppt=Bool(fixture.ppt),
            bosonic=Bool(fixture.bosonic),
        )
        @test result.verdict === Bool(fixture.expected)
        @test result.status in (
            SymExtOracleQET.SymmetricExtensionExactPresent,
            SymExtOracleQET.SymmetricExtensionAnalyticPresent,
            SymExtOracleQET.SymmetricExtensionAnalyticAbsent,
        )
    end
end

@testset "QETLAB Jacobi helper fixtures" begin
    @test length(payload.jacobi_fixtures) == 3
    for fixture in payload.jacobi_fixtures
        actual = SymExtOracleQET._jacobi_polynomial_coefficients(
            Float64(fixture.alpha), Float64(fixture.beta), Int(fixture.degree)
        )
        expected = Float64.(fixture.coefficients)
        @test actual ≈ expected atol = 2e-14 rtol = 2e-14
    end
end

@testset "QETLAB and Julia random-PPT property evidence" begin
    fixture = payload.random_fixture
    @test String(fixture.name) == "full_rank_shifted_seeded_properties"
    @test Float64(fixture.hermiticity_residual) <= 2e-12
    @test Float64(fixture.minimum_eigenvalue) >= -2e-12
    @test Float64(fixture.minimum_partial_transpose_eigenvalue) >= -2e-12
    @test Float64(fixture.trace) ≈ 1 atol = 2e-14
    @test Int(fixture.numerical_rank) <= Int(fixture.requested_rank)
    @test Int(fixture.partial_transpose_numerical_rank) <= Int(fixture.requested_rank)

    native = SymExtOracleQET.random_ppt_state(
        Xoshiro(Int(fixture.seed)), Tuple(Int.(fixture.dimensions))
    )
    @test native.status === SymExtOracleQET.RandomPPTConstructed
    @test native.verified
    @test native.trace_residual <= native.tolerance
    @test native.minimum_eigenvalue >= -native.tolerance
    @test native.minimum_partial_transpose_eigenvalue >= -native.tolerance

    low_rank = SymExtOracleQET.random_ppt_state(
        Xoshiro(7302), Tuple(Int.(fixture.dimensions)); ranks=(3, 4)
    )
    @test low_rank.status === SymExtOracleQET.RandomPPTConstructed
    @test low_rank.construction === :separable_mixture
    @test low_rank.numerical_ranks[1] <= 3
    @test low_rank.numerical_ranks[2] <= 4
end

println(
    "Symmetric-extension/Jacobi/random-PPT oracle comparison passed for ",
    length(payload.extension_fixtures),
    " extension and ",
    length(payload.jacobi_fixtures),
    " Jacobi fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
