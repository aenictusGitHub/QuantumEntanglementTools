using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const SeparabilityOracleQET = QuantumEntanglementTools

if !isdefined(SeparabilityOracleQET, :StateDiscriminationResult)
    Base.include(
        SeparabilityOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "optimization", "state_discrimination.jl"),
    )
end
if !isdefined(SeparabilityOracleQET, :LocalDistinguishabilityResult)
    Base.include(
        SeparabilityOracleQET,
        joinpath(
            @__DIR__, "..", "..", "src", "entanglement", "separability_optimization.jl"
        ),
    )
end

const SEPARABILITY_PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const PINNED_IS_SEPARABLE_SHA256 = "be078f7087c1d64ca8c1e29778fbb78690d9d06b569dc6a617cbb1b5a3656f97"
const PINNED_LOCAL_DISTINGUISHABILITY_SHA256 = "a83329f5ab1d4e751ad55748a56907dfa25bb7241909dad5c0222a4f55d47853"
const PINNED_UPB_SEP_DISTINGUISHABLE_SHA256 = "651c0ff6d8a59d8da6bf1d0632e7e118362b02e9742262e576c39e7aa4491035"
const SEPARABILITY_FIXTURE_NAMES = Set([
    "one_dimensional_local_factor",
    "two_qubit_product",
    "two_qubit_bell",
    "three_by_three_maximally_mixed",
    "three_by_three_isotropic_q_0_2",
    "three_by_three_isotropic_q_0_3",
])
const LOCAL_FIXTURE_NAMES = Set([
    "single_scaled_pure_state", "deterministic_prior_scaled_pure_states"
])

function verify_separability_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function separability_oracle_vector(value)
    value isa AbstractVector && return Float64.(value)
    return Float64[Float64(value)]
end

function separability_oracle_matrix(fixture)
    dimensions = Tuple(Int.(fixture.matrix_dims))
    real_parts = Float64.(fixture.matrix_real)
    imaginary_parts = Float64.(fixture.matrix_imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real matrix data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary matrix data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

function local_oracle_states(fixture)
    dimensions = Tuple(Int.(fixture.input_dims))
    real_parts = Float64.(fixture.state_real)
    imaginary_parts = Float64.(fixture.state_imaginary)
    length(real_parts) == prod(dimensions) ||
        error("local fixture $(fixture.name) has inconsistent state data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "separability_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "separability_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "separability oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_separability_oracle.sh first",
)

digest = verify_separability_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP6-separability" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == SEPARABILITY_PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
String(metadata.is_separable_sha256) == PINNED_IS_SEPARABLE_SHA256 ||
    error("fixture records the wrong IsSeparable.m hash")
String(metadata.local_distinguishability_sha256) ==
PINNED_LOCAL_DISTINGUISHABILITY_SHA256 ||
    error("fixture records the wrong LocalDistinguishability.m hash")
String(metadata.upb_sep_distinguishable_sha256) == PINNED_UPB_SEP_DISTINGUISHABLE_SHA256 ||
    error("fixture records the wrong UPBSepDistinguishable.m hash")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Int(metadata.local_fixture_count) == length(payload.local_fixtures) ||
    error("local fixture count does not match metadata")
Set(String(fixture.name) for fixture in payload.fixtures) == SEPARABILITY_FIXTURE_NAMES ||
    error("separability fixture names do not match the reviewed comparator set")
Set(String(fixture.name) for fixture in payload.local_fixtures) == LOCAL_FIXTURE_NAMES ||
    error("local fixture names do not match the reviewed comparator set")
occursin("Xsep2", String(metadata.is_separable_stale_update_defect)) ||
    error("fixture does not record the randomized stale-update defect")
occursin("floating boundary unknown", String(metadata.rank_one_tolerance_promotion)) ||
    error("fixture does not record the rank-one tolerance-promotion difference")
occursin("nonconjugating", String(metadata.upb_complex_transpose_defect)) ||
    error("fixture does not record the complex-transpose defect")
occursin("2-by-2", String(metadata.local_trivial_measurement_shape_defect)) ||
    error("fixture does not record the trivial-measurement shape defect")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
    !Bool(metadata.cvx_available) ||
        error("the committed solver-free fixture unexpectedly records CVX")
end

@testset "QETLAB separability solver-free differential fixture" begin
    julia_certificates = Dict(
        "one_dimensional_local_factor" => :one_dimensional_local_factor_theorem,
        "two_qubit_product" => :computational_basis_product_decomposition,
        "two_qubit_bell" => :negative_partial_transpose_witness,
        "three_by_three_maximally_mixed" => :computational_basis_product_decomposition,
        "three_by_three_isotropic_q_0_3" => :negative_partial_transpose_witness,
    )
    qetlab_to_status = Dict(-1 => :unknown, 0 => :entangled, 1 => :separable)
    for fixture in payload.fixtures
        name = String(fixture.name)
        matrix = separability_oracle_matrix(fixture)
        dims = Tuple(Int.(fixture.dims))
        qetlab_status = qetlab_to_status[Int(fixture.qetlab_sep)]
        @test Symbol(String(fixture.expected_status)) === qetlab_status
        represented_trace = SeparabilityOracleQET._sepopt_represented_trace(matrix)
        normalization_residual = abs(real(represented_trace) - 1)
        @test normalization_residual <= Float64(fixture.atol)
        if !iszero(normalization_residual)
            # QETLAB accepted these decimal fixtures within its tolerance. The
            # native certificate API intentionally requires the represented
            # trace to be exactly one and does not repair the fixture.
            @test_throws ArgumentError SeparabilityOracleQET.is_separable(
                matrix, dims; strategies=:qetlab_deterministic
            )
            continue
        end
        result = SeparabilityOracleQET.is_separable(
            matrix, dims; strategies=:qetlab_deterministic
        )
        if name == "three_by_three_isotropic_q_0_2"
            @test qetlab_status === :separable
            @test result.status === :unknown
            @test !result.certified
            @test result.certificate_kind === nothing
            rank_one_attempt = only(
                attempt for
                attempt in result.attempts if attempt.method === :rank_one_identity
            )
            @test rank_one_attempt.status === :unknown
            @test !rank_one_attempt.certified
        else
            @test result.status === qetlab_status
            @test result.certified
            @test result.certificate_kind === julia_certificates[name]
        end
        @test represented_trace == 1
    end

    @test Int(metadata.trace_two_bell_accepted_sep) == 0
    bell_fixture = only(
        fixture for fixture in payload.fixtures if String(fixture.name) == "two_qubit_bell"
    )
    trace_two_bell = 2 * separability_oracle_matrix(bell_fixture)
    @test_throws ArgumentError SeparabilityOracleQET.is_separable(trace_two_bell, (2, 2))
end

@testset "QETLAB local-discrimination trivial differential fixture" begin
    for fixture in payload.local_fixtures
        scaled_states = local_oracle_states(fixture)
        priors = separability_oracle_vector(fixture.priors)
        @test Float64(fixture.qetlab_distance) == 1
        @test Tuple(Int.(fixture.qetlab_measurement_dims)) == (2, 2)
        @test Tuple(Int.(fixture.qetlab_dual_dims)) == (4, 4)
        @test_throws ArgumentError SeparabilityOracleQET.local_distinguishability(
            scaled_states, (2, 2); priors=priors
        )

        normalized_states = copy(scaled_states)
        for column in axes(normalized_states, 2)
            normalized_states[:, column] ./= norm(normalized_states[:, column])
        end
        result = SeparabilityOracleQET.local_distinguishability(
            normalized_states, (2, 2); priors=priors
        )
        @test result.status === SeparabilityOracleQET.LocalDistinguishabilityTrivialOptimal
        @test result.relaxation_value == Float64(fixture.qetlab_distance)
        @test result.certified
        @test all(
            size(effect) == (
                Int(fixture.native_measurement_dimension),
                Int(fixture.native_measurement_dimension),
            ) for effect in result.measurement
        )
        @test result.residuals.valid
    end
end

println(
    "Separability oracle comparison passed for ",
    length(payload.fixtures),
    " IsSeparable fixtures and ",
    length(payload.local_fixtures),
    " LocalDistinguishability fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
