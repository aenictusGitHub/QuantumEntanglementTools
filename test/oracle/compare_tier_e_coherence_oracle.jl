using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"

function fixture_matrix(fixture)
    rows = Int(fixture.rows)
    columns = Int(fixture.columns)
    real_parts = fixture.real isa Number ? [Float64(fixture.real)] : Float64.(fixture.real)
    imaginary_parts = if fixture.imaginary isa Number
        [Float64(fixture.imaginary)]
    else
        Float64.(fixture.imaginary)
    end
    length(real_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, rows, columns)
end

as_matrix(value::Number) = reshape([value], 1, 1)

function oracle_inputs()
    plus_ququart = fill(0.5, 4)
    mixed_qubit = ComplexF64[0.6 0.2im; -0.2im 0.4]
    plus_qubit = fill(inv(sqrt(2.0)), 2)
    entropy_state = ComplexF64[0.7 0.1im; -0.1im 0.3]
    basis_qutrit = [1.0, 0.0, 0.0]
    return (; plus_ququart, mixed_qubit, plus_qubit, entropy_state, basis_qutrit)
end

function local_result(name::AbstractString, inputs)
    if name == "l1_plus_ququart"
        return QET.l1_coherence(inputs.plus_ququart)
    elseif name == "l1_mixed_qubit"
        return QET.l1_coherence(inputs.mixed_qubit)
    elseif name == "relative_entropy_plus_qubit_base2"
        return QET.relative_entropy_coherence(inputs.plus_qubit; base=2)
    elseif name == "relative_entropy_mixed_qubit_base2"
        return QET.relative_entropy_coherence(inputs.entropy_state; base=2)
    end
    return error("unknown numeric Tier E coherence oracle fixture: $name")
end

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "tier_e_coherence_oracle.json")
    committed = joinpath(
        @__DIR__, "fixtures", "tier_e_coherence_octave_11_3_qetlab_d858961.json"
    )
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier E coherence oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_e_coherence_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "E-coherence" ||
    error("expected a Tier E coherence fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")
occursin("counts entries", String(metadata.reviewed_upstream_discrepancy)) ||
    error("fixture does not record the reviewed CoherenceRank discrepancy")

inputs = oracle_inputs()

@testset "QETLAB Tier E coherence differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected = fixture_matrix(fixture)
        if name == "coherence_rank_qetlab_plus_bug"
            @testset "$name/reviewed_discrepancy" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 0
                @test QET.coherence_rank(inputs.plus_qubit) == 2
                @test QET.coherence_rank(inputs.plus_qubit) != only(expected)
            end
        elseif name == "coherence_rank_qetlab_basis_bug"
            @testset "$name/reviewed_discrepancy" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 2
                @test QET.coherence_rank(inputs.basis_qutrit) == 1
                @test QET.coherence_rank(inputs.basis_qutrit) != only(expected)
            end
        else
            result = as_matrix(local_result(name, inputs))
            @testset "$name/native" begin
                @test size(result) == size(expected)
                atol = Float64(fixture.atol)
                rtol = Float64(fixture.rtol)
                @test norm(result - expected) <= atol + rtol * norm(expected)
            end
            compatibility_result = if startswith(name, "l1_")
                QET.MATLABCompat.L1NormCoherence(
                    name == "l1_plus_ququart" ? inputs.plus_ququart : inputs.mixed_qubit
                )
            else
                QET.MATLABCompat.RelEntCoherence(
                    if name == "relative_entropy_plus_qubit_base2"
                        inputs.plus_qubit
                    else
                        inputs.entropy_state
                    end,
                )
            end
            @testset "$name/compatibility" begin
                @test norm(as_matrix(compatibility_result) - expected) <=
                    Float64(fixture.atol) + Float64(fixture.rtol) * norm(expected)
            end
        end
    end

    @testset "independent coherence identities" begin
        @test QET.l1_coherence(inputs.plus_ququart) ≈ 3
        @test QET.relative_entropy_coherence(inputs.plus_qubit; base=2) ≈ 1
        @test QET.coherence_rank(inputs.basis_qutrit) == 1
        @test QET.MATLABCompat.CoherenceRank(inputs.plus_qubit) == 2
        @test QET.MATLABCompat.CoherenceRank(inputs.basis_qutrit) == 1
    end
end

println(
    "Verified $(length(payload.fixtures)) Tier E coherence fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
println(
    "The CoherenceRank fixtures preserve a reviewed upstream bug; ",
    "the Julia API follows the documented nonzero-coefficient definition.",
)
if String(metadata.engine) == "Octave"
    println(
        "Octave evidence is function-specific and supplemental; ",
        "it is not a general MATLAB-equivalence claim.",
    )
end
