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

as_matrix(value::AbstractMatrix) = Matrix(value)
as_matrix(value::AbstractVector) = reshape(collect(value), :, 1)
as_matrix(value::Number) = reshape([value], 1, 1)

function oracle_inputs()
    operator = ComplexF64[
        (1+2im) (-2) (0.5im)
        3 (-im) (4-2im)
        0.25 (2+im) (-3)
    ]
    diagonal_state = Matrix(Diagonal([0.5, 0.3, 0.2]))
    rho = ComplexF64[0.7 0.1im; -0.1im 0.3]
    sigma = ComplexF64[0.4 0.05; 0.05 0.6]
    bell = [1.0, 0.0, 0.0, 1.0] / sqrt(2)
    bell_density = bell * adjoint(bell)
    schmidt_vector = ComplexF64[1, 2im, 3, 4im, 5, 6im]
    mixed_bell = 0.7 * bell_density + 0.3 * Matrix{Float64}(I, 4, 4) / 4
    maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
    return (;
        operator,
        diagonal_state,
        rho,
        sigma,
        bell_density,
        schmidt_vector,
        mixed_bell,
        maximally_mixed,
    )
end

function local_result(name::AbstractString, inputs)
    if name == "trace_norm"
        return QET.trace_norm(inputs.operator)
    elseif name == "schatten_norm_p3"
        return QET.schatten_norm(inputs.operator, 3)
    elseif name == "ky_fan_norm_k2"
        return QET.ky_fan_norm(inputs.operator, 2)
    elseif name == "purity"
        return QET.purity(inputs.diagonal_state)
    elseif name == "entropy_base3_alpha1"
        return QET.von_neumann_entropy(inputs.diagonal_state; base=3)
    elseif name == "fidelity_root"
        return QET.fidelity(inputs.rho, inputs.sigma)
    elseif name == "negativity_bell"
        return QET.negativity(inputs.bell_density, (2, 2))
    elseif name == "schmidt_coefficients_2x3"
        return QET.schmidt_coefficients(inputs.schmidt_vector, (2, 3))
    elseif name == "schmidt_rank_2x3"
        return QET.schmidt_rank(inputs.schmidt_vector, (2, 3); atol=1e-10, rtol=0)
    elseif name == "concurrence_mixed_bell"
        return QET.concurrence(inputs.mixed_bell)
    elseif name == "realignment_trace_norm_bell"
        return QET.realignment_criterion(inputs.bell_density, (2, 2)).value
    end
    return error("unknown numeric Tier D oracle fixture: $name")
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
    generated = joinpath(@__DIR__, "generated", "tier_d_oracle.json")
    committed = joinpath(@__DIR__, "fixtures", "tier_d_octave_11_3_qetlab_d858961.json")
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier D oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_d_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "D" ||
    error("expected a Tier D fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")

inputs = oracle_inputs()
fixture_lookup = Dict(String(fixture.name) => fixture for fixture in payload.fixtures)

@testset "QETLAB Tier D differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected = fixture_matrix(fixture)

        if name == "is_ppt_bell_detected"
            # QETLAB's false result and the native structured certificate have
            # different types. Check their matching mathematical meaning
            # explicitly rather than coercing CriterionResult to a boolean.
            native = QET.ppt_criterion(inputs.bell_density, (2, 2); atol=1e-10, rtol=0)
            @testset "$name/native_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 0
                @test native.status === QET.CriterionEntanglementDetected
                @test native.value < native.threshold - native.tolerance
            end
        elseif name == "is_ppt_maximally_mixed"
            # CriterionSatisfied means only that this necessary condition
            # holds with margin; it is not treated as a separability result.
            native = QET.ppt_criterion(inputs.maximally_mixed, (2, 2); atol=1e-10, rtol=0)
            @testset "$name/native_semantics" begin
                @test size(expected) == (1, 1)
                @test only(expected) == 1
                @test native.status === QET.CriterionSatisfied
                @test native.value > native.threshold + native.tolerance
            end
        else
            result = as_matrix(local_result(name, inputs))
            @testset "$name/native" begin
                @test size(result) == size(expected)
                if String(fixture.comparison) == "exact"
                    @test result == expected
                elseif String(fixture.comparison) == "normwise"
                    atol = Float64(fixture.atol)
                    rtol = Float64(fixture.rtol)
                    @test norm(result - expected) <= atol + rtol * norm(expected)
                else
                    error(
                        "unsupported comparison $(fixture.comparison) " *
                        "for fixture $name",
                    )
                end
            end
        end
    end

    @testset "phase-independent Schmidt property" begin
        decomposition = QET.schmidt_decomposition(inputs.schmidt_vector, (2, 3))
        reconstructed = sum(
            decomposition.coefficients[index] * kron(
                decomposition.left_vectors[:, index], decomposition.right_vectors[:, index]
            ) for index in eachindex(decomposition.coefficients)
        )
        @test reconstructed ≈ inputs.schmidt_vector
    end

    @testset "fidelity convention and criterion semantics" begin
        root_fidelity = only(fixture_matrix(fixture_lookup["fidelity_root"]))
        @test QET.fidelity(inputs.rho, inputs.sigma; squared=true) ≈ root_fidelity^2

        realignment = QET.realignment_criterion(inputs.bell_density, (2, 2))
        @test realignment.status === QET.CriterionEntanglementDetected
        @test realignment.value > realignment.threshold + realignment.tolerance
    end
end

println(
    "Verified $(length(payload.fixtures)) native Tier D fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
if String(metadata.engine) == "Octave"
    println(
        "Octave evidence is function-specific and supplemental; ",
        "it is not a general MATLAB-equivalence claim.",
    )
end
