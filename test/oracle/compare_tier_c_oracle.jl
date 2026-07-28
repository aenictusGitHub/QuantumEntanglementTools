using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QuantumEntanglementTools.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"

function fixture_matrix(fixture)
    rows = Int(fixture.rows)
    columns = Int(fixture.columns)
    real_parts = Float64.(fixture.real)
    imaginary_parts = Float64.(fixture.imaginary)
    length(real_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == rows * columns ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values =
        all(iszero, imaginary_parts) ? real_parts : complex.(real_parts, imaginary_parts)
    return reshape(values, rows, columns)
end

function local_result(name::AbstractString)
    if name == "depolarizing_channel"
        return QET.choi_matrix(QET.depolarizing_channel(3, 0.2))
    elseif name == "dephasing_channel"
        return QET.choi_matrix(QET.dephasing_channel(3, 0.25))
    elseif name == "pauli_channel"
        return QET.choi_matrix(QET.pauli_channel([0.1, 0.2, 0.3, 0.4]))
    elseif name == "choi_map"
        return QET.choi_matrix(QET.choi_map(1, 1, 0))
    elseif name == "reduction_map"
        return QET.choi_matrix(QET.reduction_map(3, 2))
    elseif name == "apply_map"
        input = ComplexF64[
            0.6 0.1+0.2im 0.0
            0.1-0.2im 0.3 0.05im
            0.0 -0.05im 0.1
        ]
        return QET.apply_channel(input, QET.depolarizing_channel(3, 0.2))
    elseif name == "partial_map"
        product_input = kron(
            ComplexF64[0.7 0.1im; -0.1im 0.3], ComplexF64[0.4 0.05; 0.05 0.6]
        )
        return QET.partial_map(product_input, QET.dephasing_channel(2), 2, (2, 2))
    end
    return error("unknown Tier C oracle fixture: $name")
end

function compat_result(name::AbstractString)
    if name == "depolarizing_channel"
        return Compat.DepolarizingChannel(3, 0.2)
    elseif name == "dephasing_channel"
        return Compat.DephasingChannel(3, 0.25)
    elseif name == "pauli_channel"
        return Compat.PauliChannel([0.1, 0.2, 0.3, 0.4])
    elseif name == "choi_map"
        return Compat.ChoiMap(1, 1, 0)
    elseif name == "reduction_map"
        return Compat.ReductionMap(3, 2)
    elseif name == "apply_map"
        input = ComplexF64[
            0.6 0.1+0.2im 0.0
            0.1-0.2im 0.3 0.05im
            0.0 -0.05im 0.1
        ]
        return Compat.ApplyMap(input, Compat.DepolarizingChannel(3, 0.2))
    elseif name == "partial_map"
        product_input = kron(
            ComplexF64[0.7 0.1im; -0.1im 0.3], ComplexF64[0.4 0.05; 0.05 0.6]
        )
        return Compat.PartialMap(product_input, Compat.DephasingChannel(2), 2, (2, 2))
    end
    return error("unknown Tier C compatibility oracle fixture: $name")
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
    generated = joinpath(@__DIR__, "generated", "tier_c_oracle.json")
    committed = joinpath(@__DIR__, "fixtures", "tier_c_octave_11_3_qetlab_d858961.json")
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier C oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_c_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "C" ||
    error("expected a Tier C fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")

@testset "QETLAB Tier C differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected = fixture_matrix(fixture)
        for (surface, result) in (
            ("native", Matrix(local_result(name))), ("compat", Matrix(compat_result(name)))
        )
            @testset "$name/$surface" begin
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
end

println(
    "Verified $(length(payload.fixtures)) Tier C fixtures from ",
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
