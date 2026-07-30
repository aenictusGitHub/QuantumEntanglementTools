using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const DiscriminationOracleQET = QuantumEntanglementTools

if !isdefined(DiscriminationOracleQET, :StateDiscriminationResult)
    Base.include(
        DiscriminationOracleQET,
        joinpath(@__DIR__, "..", "..", "src", "optimization", "state_discrimination.jl"),
    )
end

const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const PINNED_DISTINGUISHABILITY_SHA256 = "714b3fe26124fa08520f19ed363f7cf4c856d752f6ce8b746dfc282715f111f6"
const PINNED_NORMALIZE_COLS_SHA256 = "9c1f28f57a263233e1f947a27e7ffde896a50e4c2b5294f7aff6795b68c2310d"
const FIXTURE_NAMES = Set([
    "pure_orthogonal_equal",
    "pure_nonorthogonal_equal",
    "pure_nonorthogonal_unequal",
    "density_mixed_equal",
    "density_identical_unequal",
    "pure_orthogonal_three",
    "density_orthogonal_three",
    "single_pure_state",
])

function verify_state_discrimination_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function oracle_vector(value)
    value isa AbstractVector && return Float64.(value)
    return Float64[Float64(value)]
end

function oracle_state_data(fixture)
    dimensions = Tuple(Int.(fixture.state_dims))
    real_parts = Float64.(fixture.state_real)
    imaginary_parts = Float64.(fixture.state_imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real state data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary state data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    array = reshape(values, dimensions)
    input_kind = Symbol(String(fixture.input_kind))
    if input_kind === :pure_columns
        return array
    elseif input_kind === :density_stack
        dimension = Int(fixture.dimension)
        state_count = Int(fixture.state_count)
        size(array) == (dimension, dimension, state_count) ||
            error("fixture $(fixture.name) has inconsistent density-stack dimensions")
        return Tuple(copy(@view(array[:, :, index])) for index in 1:state_count)
    end
    return error("unreviewed input kind $(fixture.input_kind)")
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "state_discrimination_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "state_discrimination_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "state-discrimination oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_state_discrimination_oracle.sh first",
)

digest = verify_state_discrimination_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP4-state-discrimination" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
String(metadata.distinguishability_sha256) == PINNED_DISTINGUISHABILITY_SHA256 ||
    error("fixture records the wrong Distinguishability.m hash")
String(metadata.normalize_cols_sha256) == PINNED_NORMALIZE_COLS_SHA256 ||
    error("fixture records the wrong normalize_cols.m hash")
occursin("can exceed 1", String(metadata.unequal_pure_prior_formula_defect)) ||
    error("fixture does not record the unequal-prior pure-state formula defect")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
Set(String(fixture.name) for fixture in payload.fixtures) == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

@testset "QETLAB state-discrimination solver-free differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        states = oracle_state_data(fixture)
        priors = oracle_vector(fixture.priors)
        result = DiscriminationOracleQET.state_distinguishability(states; priors=priors)
        expected = Float64(fixture.expected_value)
        measurement_objective = Float64(fixture.measurement_objective)
        atol = Float64(fixture.atol)
        rtol = Float64(fixture.rtol)
        @test result.success_probability !== nothing
        if name == "pure_nonorthogonal_unequal"
            @test expected > 1
            @test !isapprox(result.success_probability, expected; atol=atol, rtol=rtol)
            @test isapprox(
                result.success_probability, measurement_objective; atol=atol, rtol=rtol
            )
            @test isapprox(result.lower_bound, measurement_objective; atol=atol, rtol=rtol)
            @test isapprox(result.upper_bound, measurement_objective; atol=atol, rtol=rtol)
            @test Float64(fixture.measurement_objective_residual) > 0.09
        else
            @test isapprox(result.success_probability, expected; atol=atol, rtol=rtol)
            @test isapprox(result.lower_bound, expected; atol=atol, rtol=rtol)
            @test isapprox(result.upper_bound, expected; atol=atol, rtol=rtol)
            @test isapprox(measurement_objective, expected; atol=atol, rtol=rtol)
            @test Float64(fixture.measurement_objective_residual) <= 2e-10
        end
        @test result.certified
        @test result.measurement !== nothing
        @test result.residuals.valid
        @test result.residuals.completeness_residual <= 2e-10
        @test result.residuals.positivity_violation <= 2e-10
        @test Float64(fixture.measurement_completeness_residual) <= 2e-10
        @test Float64(fixture.measurement_positivity_violation) <= 2e-10
        if occursin("three", name)
            @test result.status ===
                DiscriminationOracleQET.StateDiscriminationOrthogonalOptimal
        elseif name == "single_pure_state"
            @test result.status ===
                DiscriminationOracleQET.StateDiscriminationTrivialOptimal
        else
            @test result.status ===
                DiscriminationOracleQET.StateDiscriminationHelstromOptimal
        end
    end

    @test Float64(metadata.pure_silent_normalization_value) == 1
    @test Float64(metadata.density_silent_normalization_value) == 1
    @test Float64(metadata.negative_prior_accepted_value) == 1

    ket0 = ComplexF64[1, 0]
    ket1 = ComplexF64[0, 1]
    rho0 = ket0 * ket0'
    rho1 = ket1 * ket1'
    @test_throws ArgumentError DiscriminationOracleQET.state_distinguishability(
        hcat(2ket0, 3ket1)
    )
    @test_throws ArgumentError DiscriminationOracleQET.state_distinguishability((
        2rho0, 3rho1
    ))
    @test_throws DomainError DiscriminationOracleQET.state_distinguishability(
        (rho0, rho1); priors=[1.1, -0.1]
    )
end

println(
    "State-discrimination oracle comparison passed for ",
    length(payload.fixtures),
    " deterministic fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
