using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using Random
using SHA
using Test

const RandomSuperoperatorQET = QuantumEntanglementTools
const RandomSuperoperatorCompat = QuantumEntanglementTools.MATLABCompat
const RANDOM_SUPEROPERATOR_PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const RANDOM_SUPEROPERATOR_PINNED_SOURCE_SHA = "ac2862092569e0878c0f5dac68f4cd58670992ae44adea26d2d7dfb36c07d96f"
const RANDOM_SUPEROPERATOR_FIXTURE_NAMES = Set([
    "unconstrained_complex",
    "trace_preserving_real",
    "unital_real",
    "bistochastic_rank_one_real",
    "unequal_balanced_is_not_unital",
])

if !isdefined(RandomSuperoperatorQET, :RandomSuperoperatorResult)
    Base.include(
        RandomSuperoperatorQET,
        joinpath(@__DIR__, "..", "..", "src", "channels", "random_superoperator.jl"),
    )
end

function random_superoperator_verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function random_superoperator_fixture_matrix(fixture)
    dimensions = Tuple(Int.(fixture.matrix_dimensions))
    dimensions == (prod(Int.(fixture.dimensions)), prod(Int.(fixture.dimensions))) ||
        error("fixture $(fixture.name) has inconsistent Choi dimensions")
    real_parts = Float64.(fixture.real)
    imaginary_parts = Float64.(fixture.imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, dimensions)
end

function random_superoperator_choi_marginals(choi, input_dimension, output_dimension)
    input_marginal = zeros(eltype(choi), input_dimension, input_dimension)
    output_marginal = zeros(eltype(choi), output_dimension, output_dimension)
    for input_column in 1:input_dimension,
        input_row in 1:input_dimension,
        output_index in 1:output_dimension

        input_marginal[input_row, input_column] += choi[
            output_index + (input_row - 1) * output_dimension,
            output_index + (input_column - 1) * output_dimension,
        ]
    end
    for output_column in 1:output_dimension,
        output_row in 1:output_dimension,
        input_index in 1:input_dimension

        output_marginal[output_row, output_column] += choi[
            output_row + (input_index - 1) * output_dimension,
            output_column + (input_index - 1) * output_dimension,
        ]
    end
    return input_marginal, output_marginal
end

committed_fixture = joinpath(
    @__DIR__, "fixtures", "random_superoperator_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "random_superoperator_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "RandomSuperoperator oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_random_superoperator_oracle.sh first",
)

digest = random_superoperator_verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-random-superoperator" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == RANDOM_SUPEROPERATOR_PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
String(metadata.source_sha256) == RANDOM_SUPEROPERATOR_PINNED_SOURCE_SHA ||
    error("fixture does not record the pinned RandomSuperoperator.m SHA-256")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
occursin("not compared entrywise", String(metadata.comparison_scope)) ||
    error("fixture does not state the cross-engine random-stream limitation")
occursin("cannot be unital", String(metadata.unequal_dimension_defect)) ||
    error("fixture does not record the unequal-dimensional unitality defect")
occursin("TP=2", String(metadata.permissive_flag_behavior)) ||
    error("fixture does not record permissive upstream flag handling")
occursin("larger than prod(DIM)", String(metadata.oversized_rank_behavior)) ||
    error("fixture does not record permissive oversized-rank handling")
Int(metadata.oversized_rank_numerical_rank) == 4 ||
    error("fixture does not demonstrate the impossible KR=5 exact-rank claim")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(String(fixture.name) => fixture for fixture in payload.fixtures)
Set(keys(fixtures)) == RANDOM_SUPEROPERATOR_FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

@testset "QETLAB RandomSuperoperator seeded property fixture" begin
    for name in sort!(collect(keys(fixtures)))
        fixture = fixtures[name]
        choi = random_superoperator_fixture_matrix(fixture)
        input_dimension, output_dimension = Int.(fixture.dimensions)
        input_marginal, output_marginal = random_superoperator_choi_marginals(
            choi, input_dimension, output_dimension
        )
        proportional_factor = input_dimension / output_dimension
        input_identity = Matrix{eltype(choi)}(I, input_dimension, input_dimension)
        output_identity = Matrix{eltype(choi)}(I, output_dimension, output_dimension)
        trace_residual = norm(input_marginal - input_identity)
        unital_residual = norm(output_marginal - output_identity)
        proportional_residual = norm(
            output_marginal - proportional_factor * output_identity
        )

        @test norm(choi - adjoint(choi)) <= 2e-12
        @test eigmin(Hermitian(choi)) >= -2e-12
        @test rank(choi; atol=1e-10) == Int(fixture.numerical_rank)
        @test Int(fixture.numerical_rank) <= Int(fixture.kraus_rank)
        @test tr(choi) ≈ Float64(fixture.choi_trace) atol = 3e-14
        @test trace_residual ≈ Float64(fixture.trace_preservation_residual) atol = 3e-14
        @test unital_residual ≈ Float64(fixture.unitality_residual) atol = 3e-14
        @test proportional_residual ≈ Float64(fixture.proportional_unitality_residual) atol =
            3e-14
        if Int(fixture.real_output) == 1
            @test eltype(choi) == Float64
        end
        if Int(fixture.trace_preserving) == 1
            @test trace_residual <= 5e-8
        end
        if Int(fixture.unital) == 1 && input_dimension == output_dimension
            @test unital_residual <= 5e-8
        elseif Int(fixture.unital) == 1
            @test unital_residual > 0.5
            @test proportional_residual <= 5e-8
        end
    end
end

@testset "Julia random-superoperator branch properties" begin
    cases = (
        (
            name=:unconstrained,
            seed=4101,
            dims=(2, 3),
            trace_preserving=false,
            unital=false,
            proportional=false,
            real=false,
            rank=3,
        ),
        (
            name=:trace_preserving,
            seed=4102,
            dims=(3, 3),
            trace_preserving=true,
            unital=false,
            proportional=false,
            real=true,
            rank=2,
        ),
        (
            name=:unital,
            seed=4103,
            dims=(3, 3),
            trace_preserving=false,
            unital=true,
            proportional=false,
            real=true,
            rank=2,
        ),
        (
            name=:bistochastic,
            seed=4104,
            dims=(2, 2),
            trace_preserving=true,
            unital=true,
            proportional=false,
            real=true,
            rank=1,
        ),
        (
            name=:proportional,
            seed=4105,
            dims=(2, 3),
            trace_preserving=true,
            unital=false,
            proportional=true,
            real=true,
            rank=2,
        ),
    )
    for case in cases
        result = RandomSuperoperatorQET.random_superoperator(
            Xoshiro(case.seed),
            case.dims;
            trace_preserving=case.trace_preserving,
            unital=case.unital,
            proportional_unital=case.proportional,
            real=case.real,
            kraus_rank=case.rank,
            representation=:choi,
        )
        choi = RandomSuperoperatorQET.choi_matrix(result.representation)
        @test result.status === :success
        @test result.complete_positivity_guaranteed
        @test ishermitian(choi)
        @test eigmin(Hermitian(choi)) >= -3e-12
        @test result.numerical_kraus_rank <= case.rank
        @test case.real ? eltype(choi) == Float64 : eltype(choi) == ComplexF64
        if case.trace_preserving
            @test result.trace_preservation_guaranteed
            @test result.trace_preservation_residual <= result.trace_preservation_tolerance
        end
        if case.unital
            @test result.unitality_guaranteed
            @test result.unitality_residual <= result.unitality_tolerance
        end
        if case.proportional
            @test result.proportional_unitality_guaranteed
            @test !result.unitality_guaranteed
            @test result.proportional_unitality_residual <= result.unitality_tolerance
        end
    end
    @test_throws ArgumentError RandomSuperoperatorQET.random_superoperator(
        Xoshiro(4199), (2, 3); trace_preserving=2
    )
    @test_throws ArgumentError RandomSuperoperatorQET.random_superoperator(
        Xoshiro(4200), (2, 2); trace_preserving=false, kraus_rank=5
    )
end

@testset "RandomSuperoperator compatibility surface" begin
    raw = RandomSuperoperatorCompat.RandomSuperoperator(Xoshiro(4301), (2, 3), 1, 0, 1, 2)
    diagnostic = RandomSuperoperatorCompat.RandomSuperoperator(
        Xoshiro(4301), (2, 3), 1, 0, 1, 2; diagnostics=true
    )
    @test raw == RandomSuperoperatorQET.choi_matrix(diagnostic.representation)
    @test size(raw) == (6, 6)
    @test eltype(raw) == Float64
    @test diagnostic.status === :success
    @test diagnostic.trace_preservation_guaranteed

    proportional = RandomSuperoperatorCompat.RandomSuperoperator(
        Xoshiro(4302), (2, 3), 1, 1, 1, 2; diagnostics=true, allow_proportional_unital=true
    )
    @test proportional.status === :success
    @test proportional.proportional_unitality_guaranteed
    @test !proportional.unitality_guaranteed
    @test proportional.proportional_unitality_residual <= proportional.unitality_tolerance
    @test_throws ArgumentError RandomSuperoperatorCompat.RandomSuperoperator(
        Xoshiro(4303), (2, 3), 1, 1
    )
end

println(
    "RandomSuperoperator oracle/property comparison passed for ",
    length(payload.fixtures),
    " pinned fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
