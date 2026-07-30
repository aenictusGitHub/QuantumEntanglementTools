using JSON3
using QuantumEntanglementTools
using SHA
using Test

const NonlocalOracleQET = QuantumEntanglementTools
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const XOR_FIXTURE_NAMES = Set(["xor_chsh_classical", "xor_rectangular_classical"])
const BELL_FIXTURE_NAMES = Set([
    "bell_chsh_full_correlator_classical",
    "bell_affine_full_correlator_classical",
    "bell_ternary_binary_full_probability_classical",
])

for source_file in ("scenarios.jl", "npa.jl", "game_values.jl", "nonlocal_optimization.jl")
    symbol = if source_file == "scenarios.jl"
        :BellScenario
    elseif source_file == "npa.jl"
        :NPAWord
    elseif source_file == "game_values.jl"
        :NonlocalValueResult
    else
        :NonlocalLowerBoundResult
    end
    if !isdefined(NonlocalOracleQET, symbol)
        Base.include(
            NonlocalOracleQET,
            joinpath(@__DIR__, "..", "..", "src", "nonlocal_games", source_file),
        )
    end
end

module _NonlocalOracleCompatHarness
using QuantumEntanglementTools:
    BellScenario,
    CollinsGisinBehavior,
    NoOptimizationBackend,
    bcs_game_lower_bound,
    bcs_game_value,
    bell_inequality_bound,
    bell_inequality_qubit_bound,
    nonlocal_game_lower_bound,
    npa_membership,
    xor_game_value

include(joinpath(@__DIR__, "..", "..", "src", "compat", "nonlocal_games.jl"))
end

const NonlocalOracleCompat = _NonlocalOracleCompatHarness

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function fixture_array(values, dimensions, ::Type{T}) where {T}
    shape = Tuple(Int.(dimensions))
    length(values) == prod(shape) ||
        error("fixture array data are inconsistent with dimensions $shape")
    return reshape(T.(values), shape)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "nonlocal_games_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "nonlocal_games_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Nonlocal-games oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_nonlocal_games_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP6-nonlocal-games" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) ==
length(payload.xor_fixtures) + length(payload.bell_fixtures) + 1 ||
    error("fixture count does not match metadata")
occursin("NonlocalGameValue.m", String(metadata.upstream_blocker)) ||
    error("fixture does not record the pinned BCSGameValue dependency blocker")
occursin("CVX-dependent", String(metadata.excluded_scope)) ||
    error("fixture does not delimit its solver-free evidence")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

Set(String(fixture.name) for fixture in payload.xor_fixtures) == XOR_FIXTURE_NAMES ||
    error("XOR fixture names do not match the reviewed comparator set")
Set(String(fixture.name) for fixture in payload.bell_fixtures) == BELL_FIXTURE_NAMES ||
    error("Bell fixture names do not match the reviewed comparator set")
String(payload.bcs_fixture.name) == "bcs_to_nonlocal_active_and_inactive_variables" ||
    error("unexpected BCS helper fixture")

@testset "QETLAB nonlocal-game differential fixtures" begin
    @testset "XORGameValue classical" begin
        for fixture in payload.xor_fixtures
            probabilities = fixture_array(
                fixture.probabilities, fixture.probability_dims, Float64
            )
            parity = fixture_array(fixture.parity, fixture.parity_dims, Int)
            native = NonlocalOracleQET.xor_game_value(
                probabilities, parity; regime=:classical
            )
            compatibility = NonlocalOracleCompat.XORGameValue(
                probabilities, parity, "classical"; structured=false
            )
            @test native.status === NonlocalOracleQET.NonlocalValueExact
            @test native.exact
            @test native.certificate_kind === :exhaustive_xor_strategy_reduction
            @test native.value ≈ Float64(fixture.expected) atol = Float64(fixture.atol) rtol = Float64(
                fixture.rtol
            )
            @test compatibility ≈ Float64(fixture.expected) atol = Float64(fixture.atol) rtol = Float64(
                fixture.rtol
            )
        end
    end

    @testset "BellInequalityMax classical" begin
        for fixture in payload.bell_fixtures
            coefficients = fixture_array(
                fixture.coefficients, fixture.coefficient_dims, Int
            )
            desc = Int.(fixture.desc)
            notation = Symbol(String(fixture.notation))
            native = NonlocalOracleQET.bell_inequality_bound(
                coefficients,
                NonlocalOracleQET.BellScenario(desc);
                notation=notation,
                regime=:classical,
            )
            compatibility = NonlocalOracleCompat.BellInequalityMax(
                coefficients, desc, notation, "classical", 1; structured=false
            )
            @test String(fixture.comparison) == "exact"
            @test native.status === NonlocalOracleQET.NonlocalValueExact
            @test native.exact
            @test native.value == Int(fixture.expected)
            @test compatibility == Int(fixture.expected)
        end
    end

    @testset "bcs_to_nonlocal helper" begin
        fixture = payload.bcs_fixture
        constraint_1 = fixture_array(fixture.constraint_1, fixture.constraint_1_dims, Int)
        constraint_2 = fixture_array(fixture.constraint_2, fixture.constraint_2_dims, Int)
        game = NonlocalOracleQET.nonlocal_game(
            NonlocalOracleQET.BCSGame([constraint_1, constraint_2])
        )
        expected_probabilities = fixture_array(
            fixture.probabilities, fixture.probability_dims, Float64
        )
        expected_payoff = fixture_array(fixture.payoff, fixture.payoff_dims, Int)
        @test String(fixture.comparison) == "exact"
        @test size(game.probabilities) == size(expected_probabilities)
        @test Float64.(game.probabilities) == expected_probabilities
        @test size(game.payoff) == size(expected_payoff)
        @test game.payoff == expected_payoff
    end
end

println(
    "Nonlocal-game oracle comparison passed for ",
    metadata.fixture_count,
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
