using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SHA
using Test

const SKNormOracleQET = QuantumEntanglementTools
const SKNORM_PINNED_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const KPNORMDUAL_SOURCE_SHA = "ecedafbbb1a58977b9447006942c23f729e79b8445139733aa6fb3ee2da2730b"
const SKOPERATORNORM_SOURCE_SHA = "67fc32b39abb01fb1c098750b19de1ef0cadcf622d610ec9446d9d2a0d87eaac"
const ISBLOCKPOSITIVE_SOURCE_SHA = "cfb38c34eacb699e8eb5760f1d1e6257c4f0b4a574f3713f0c3d6a2472d575a1"

if !isdefined(SKNormOracleQET, :TopKPNormDualEpigraph)
    Base.include(
        SKNormOracleQET,
        joinpath(
            @__DIR__, "..", "..", "src", "optimization", "top_k_p_norm_dual_epigraph.jl"
        ),
    )
end
if !isdefined(SKNormOracleQET, :SKOperatorNormResult)
    Base.include(
        SKNormOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "optimization", "sk_operator_norm.jl"),
    )
end
if !isdefined(SKNormOracleQET, :BlockPositivityResult)
    Base.include(
        SKNormOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "entanglement", "block_positivity.jl"),
    )
end

function sknorm_oracle_digest(path)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function sknorm_fixture_array(fixture)
    dimensions = Tuple(Int.(fixture.dimensions))
    real_parts = fixture.real isa Number ? Float64[fixture.real] : Float64.(fixture.real)
    imaginary_parts = if fixture.imaginary isa Number
        Float64[fixture.imaginary]
    else
        Float64.(fixture.imaginary)
    end
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

function sknorm_named_matrix(name)
    name == "diagonal" && return Diagonal([4.0, 3.0, 2.0, 1.0])
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    product = ComplexF64[1, 0, 0, 0]
    name == "rank_one" && return 2bell * product'
    name == "identity" && return Matrix{Float64}(I, 4, 4)
    name == "negative_diagonal" && return Diagonal([-1.0, 2, 2, 2])
    name == "bell_boundary" && return 0.5Matrix{ComplexF64}(I, 4, 4) - bell * bell'
    name == "negative_bell" && return -bell * bell'
    return error("unreviewed fixture matrix $name")
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "sk_norms_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "sk_norms_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "S(k)-norm oracle fixture not found at $fixture_path; run " *
    "scripts/matlab_oracle/run_sk_norms_oracle.sh first",
)

digest = sknorm_oracle_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema")
String(metadata.tier) == "WP7-sk-norms-block-positivity" || error("unexpected fixture tier")
String(metadata.qetlab_commit) == SKNORM_PINNED_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.kp_norm_dual_source_sha256) == KPNORMDUAL_SOURCE_SHA ||
    error("fixture has the wrong kpNormDual.m hash")
String(metadata.sk_operator_norm_source_sha256) == SKOPERATORNORM_SOURCE_SHA ||
    error("fixture has the wrong SkOperatorNorm.m hash")
String(metadata.is_block_positive_source_sha256) == ISBLOCKPOSITIVE_SOURCE_SHA ||
    error("fixture has the wrong IsBlockPositive.m hash")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.dual_fixture_count) == length(payload.dual_fixtures) ||
    error("dual fixture count mismatch")
Int(metadata.sk_fixture_count) == length(payload.sk_fixtures) ||
    error("S(k) fixture count mismatch")
Int(metadata.block_fixture_count) == length(payload.block_fixtures) ||
    error("block fixture count mismatch")
occursin("CVX is not invoked", String(metadata.solver_scope)) ||
    error("fixture does not delimit its solver scope")
occursin("no random", String(metadata.random_scope)) ||
    error("fixture does not state its random-stream boundary")
occursin("inconclusive", String(metadata.boundary_scope)) ||
    error("fixture omits its numerical-boundary scope")
occursin("Hypatia and SCS", String(metadata.model_expression_scope)) ||
    error("fixture omits its model-expression validation boundary")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" || error("the committed fixture must record Octave")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB kpNormDual numeric fixtures" begin
    @test length(payload.dual_fixtures) == 4
    for fixture in payload.dual_fixtures
        input = sknorm_fixture_array(fixture)
        value = String(fixture.kind) == "vector" ? vec(input) : input
        order = if fixture.p isa Number
            Float64(fixture.p)
        elseif String(fixture.p) == "Inf"
            Inf
        else
            error("unreviewed norm order $(fixture.p)")
        end
        actual = SKNormOracleQET.top_k_p_norm_dual(value, Int(fixture.k), order)
        @test actual ≈ Float64(fixture.expected) atol = 3e-14 rtol = 3e-14
    end
end

@testset "QETLAB SkOperatorNorm solver-free exact fixtures" begin
    @test length(payload.sk_fixtures) == 2
    for fixture in payload.sk_fixtures
        result = SKNormOracleQET.sk_operator_norm(
            MersenneTwister(731),
            sknorm_named_matrix(String(fixture.matrix_name));
            k=Int(fixture.k),
            dims=Tuple(Int.(fixture.dimensions)),
            strength=Int(fixture.strength),
        )
        @test result.status === SKNormOracleQET.SKOperatorNormExact
        @test result.exact
        @test result.lower_bound ≈ Float64(fixture.lower) atol = 3e-14 rtol = 3e-14
        @test result.upper_bound ≈ Float64(fixture.upper) atol = 3e-14 rtol = 3e-14
        @test result.lower_witness.validated
    end
end

@testset "QETLAB IsBlockPositive solver-free fixtures" begin
    @test length(payload.block_fixtures) == 4
    for fixture in payload.block_fixtures
        result = SKNormOracleQET.is_block_positive(
            MersenneTwister(732),
            sknorm_named_matrix(String(fixture.matrix_name));
            k=Int(fixture.k),
            dims=Tuple(Int.(fixture.dimensions)),
            strength=Int(fixture.strength),
        )
        expected = Int(fixture.expected)
        if expected == 1
            @test result.verdict === true
            @test result.certified
        elseif expected == 0
            @test result.verdict === false
            @test result.witness.validated
        else
            @test expected == -1
            @test result.verdict === nothing
            @test !result.certified
        end
    end
end

println(
    "S(k)-norm/block-positivity oracle comparison passed for ",
    length(payload.dual_fixtures) +
    length(payload.sk_fixtures) +
    length(payload.block_fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
