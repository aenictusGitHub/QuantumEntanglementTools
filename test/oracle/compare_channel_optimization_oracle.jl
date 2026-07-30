using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const ChannelOptimizationOracleQET = QuantumEntanglementTools

if !isdefined(ChannelOptimizationOracleQET, :ChannelNormResult)
    Base.include(
        ChannelOptimizationOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "optimization", "channel_optimization.jl"),
    )
end

const CHANNEL_OPTIMIZATION_PINNED_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const CHANNEL_OPTIMIZATION_SOURCE_HASHES = Dict(
    "diamond_norm_sha256" => "a26804a1acf940b90a90d0b1fc67874cc30c0b075285ceb51849b3a0ffeccbb3",
    "cb_norm_sha256" => "04c2dd5cfe2fb6d1864b51b58fd8028018eb3a4c53fe8a6d9b9933f05d931d2c",
    "channel_distinguishability_sha256" => "dc579aefa18c9cfb11c2975e0191b93dfaaf930664b76b05b54021e38c8155f5",
    "maximum_output_fidelity_sha256" => "8a8b0900640db599558582d53683012e7c35be8f226b2864ac146ec91d1fc1b0",
)
const CHANNEL_OPTIMIZATION_FIXTURE_NAMES = Set([
    "diamond_identity",
    "diamond_rectangular_trace",
    "diamond_scaled_identity",
    "cb_identity",
    "cb_reset",
    "channel_identical_unequal",
    "channel_deterministic",
    "maximum_identical",
    "maximum_orthogonal_replacers",
    "maximum_identity_dephasing",
    "maximum_rank_truncation_defect",
])

function verify_channel_optimization_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function channel_optimization_oracle_kraus(fixture, prefix::String)
    dims = Tuple(Int.(getproperty(fixture, Symbol(prefix, "_kraus_dims"))))
    length(dims) == 3 || error("fixture $(fixture.name) has invalid Kraus dimensions")
    dims == (0, 0, 0) && return nothing
    output_dimension, input_dimension, count = dims
    real_values = Float64.(getproperty(fixture, Symbol(prefix, "_kraus_real")))
    imaginary_values = Float64.(getproperty(fixture, Symbol(prefix, "_kraus_imaginary")))
    length(real_values) == output_dimension * input_dimension * count ||
        error("fixture $(fixture.name) has inconsistent real Kraus data")
    length(imaginary_values) == length(real_values) ||
        error("fixture $(fixture.name) has inconsistent imaginary Kraus data")
    values = complex.(real_values, imaginary_values)
    stack = reshape(values, output_dimension, input_dimension, count)
    operators = Tuple(copy(@view(stack[:, :, index])) for index in 1:count)
    @assert Int(getproperty(fixture, Symbol(prefix, "_input_dimension"))) == input_dimension
    @assert Int(getproperty(fixture, Symbol(prefix, "_output_dimension"))) ==
        output_dimension
    return ChannelOptimizationOracleQET.KrausRepresentation(operators)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "channel_optimization_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "channel_optimization_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "channel-optimization oracle fixture not found at $fixture_path; run " *
    "scripts/matlab_oracle/run_channel_optimization_oracle.sh first",
)

digest = verify_channel_optimization_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP4-channel-optimization" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == CHANNEL_OPTIMIZATION_PINNED_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
!Bool(metadata.cvx_used) ||
    error("the committed fixture must not claim unavailable CVX evidence")
for (field, expected) in CHANNEL_OPTIMIZATION_SOURCE_HASHES
    String(getproperty(metadata, Symbol(field))) == expected ||
        error("fixture records the wrong $field")
end
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Set(String(fixture.name) for fixture in payload.fixtures) ==
CHANNEL_OPTIMIZATION_FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")
occursin("(1+that norm)/2", String(metadata.channel_probability_defect)) ||
    error("fixture does not record the channel-probability defect")
occursin("returns 0.8", String(metadata.maximum_output_fidelity_defect)) ||
    error("fixture does not record the Kraus-rank truncation defect")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB channel-optimization solver-free fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        routine = String(fixture.routine)
        first_map = channel_optimization_oracle_kraus(fixture, "first")
        second_map = channel_optimization_oracle_kraus(fixture, "second")
        expected_qetlab = Float64(fixture.expected_value)
        expected_native = Float64(fixture.native_expected_value)
        atol = Float64(fixture.atol)
        rtol = Float64(fixture.rtol)
        deviation = String(fixture.deviation)

        result, native_value = if routine == "diamond_norm"
            local norm_result = ChannelOptimizationOracleQET.diamond_norm(first_map)
            norm_result, norm_result.value
        elseif routine == "cb_norm"
            local norm_result = ChannelOptimizationOracleQET.cb_norm(first_map)
            norm_result, norm_result.value
        elseif routine == "channel_distinguishability"
            priors = Float64.(fixture.priors)
            local discrimination = ChannelOptimizationOracleQET.channel_distinguishability(
                first_map, second_map; priors=priors
            )
            discrimination, discrimination.success_probability
        elseif routine == "maximum_output_fidelity"
            local fidelity_result = ChannelOptimizationOracleQET.maximum_output_fidelity(
                first_map, second_map
            )
            fidelity_result, fidelity_result.value
        else
            error("unreviewed fixture routine $routine")
        end

        @test result.status ===
            ChannelOptimizationOracleQET.ChannelOptimizationAnalyticOptimal
        @test native_value !== nothing
        @test isapprox(native_value, expected_native; atol=atol, rtol=rtol)
        @test result.lower_bound !== nothing
        @test result.upper_bound !== nothing
        @test isapprox(result.lower_bound, expected_native; atol=atol, rtol=rtol)
        @test isapprox(result.upper_bound, expected_native; atol=atol, rtol=rtol)
        @test result.certified
        @test result.certificate_kind !== nothing
        if isempty(deviation)
            @test isapprox(native_value, expected_qetlab; atol=atol, rtol=rtol)
        elseif name == "channel_identical_unequal"
            @test deviation == "pinned_missing_channel_Helstrom_affine_conversion"
            @test isapprox(expected_qetlab, 0.6; atol=atol, rtol=rtol)
            @test isapprox(native_value, 0.8; atol=atol, rtol=rtol)
            @test !isapprox(native_value, expected_qetlab; atol=atol, rtol=rtol)
            @test result.certificate_kind === :identical_channels
        elseif name == "maximum_rank_truncation_defect"
            @test deviation == "pinned_minimum_Kraus_rank_truncation"
            @test isapprox(expected_qetlab, 0.8; atol=atol, rtol=rtol)
            @test native_value == 1
            @test !isapprox(native_value, expected_qetlab; atol=atol, rtol=rtol)
            @test result.certificate_kind === :common_basis_output
            @test result.output_states[1] == result.output_states[2]
            @test real(tr(result.output_states[1])) == 1
        else
            error("unreviewed deviation $deviation for fixture $name")
        end
    end

    @test Float64(metadata.negative_prior_accepted_value) == 1
    identity_channel = ChannelOptimizationOracleQET.KrausRepresentation([
        Matrix{ComplexF64}(I, 2, 2)
    ])
    dephasing_channel = ChannelOptimizationOracleQET.KrausRepresentation([
        ComplexF64[1 0; 0 0], ComplexF64[0 0; 0 1]
    ])
    @test_throws DomainError ChannelOptimizationOracleQET.channel_distinguishability(
        identity_channel, dephasing_channel; priors=[1.1, -0.1]
    )
end

println(
    "Channel-optimization oracle comparison passed for ",
    length(payload.fixtures),
    " deterministic solver-free fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
