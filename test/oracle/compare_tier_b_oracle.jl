using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const MC = QuantumEntanglementTools.MATLABCompat
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

as_matrix(value::AbstractMatrix) = Matrix(value)
as_matrix(value::AbstractVector) = reshape(collect(value), :, 1)

function local_results(name::AbstractString)
    if name == "pauli"
        return QET.pauli([1, 2, 3]), MC.Pauli([1, 2, 3], 0)
    elseif name == "generalized_pauli"
        return QET.generalized_pauli(1, 2, 3), MC.GenPauli(1, 2, 3)
    elseif name == "gell_mann"
        return QET.gell_mann(8), MC.GellMann(8)
    elseif name == "generalized_gell_mann"
        return (QET.generalized_gell_mann(2, 1, 4), MC.GenGellMann(2, 1, 4))
    elseif name == "fourier_matrix"
        return QET.fourier_matrix(4), MC.FourierMatrix(4)
    elseif name == "maximally_entangled"
        return QET.maximally_entangled(3), MC.MaxEntangled(3)
    elseif name == "bell_state"
        return QET.bell_state(1), MC.Bell(5)
    elseif name == "ghz_state"
        coefficients = [1, 2im, -3]
        return (
            QET.ghz_state(3, 2; coefficients, sparse_output=false),
            MC.GHZState(3, 2, coefficients),
        )
    elseif name == "w_state"
        coefficients = [1, 2, 3, 4]
        return (
            QET.w_state(4; coefficients, sparse_output=false), MC.WState(4, coefficients)
        )
    elseif name == "dicke_state"
        return QET.dicke_state(5, 2), MC.DickeState(5, 2)
    elseif name == "isotropic_state"
        return QET.isotropic_state(3, 0.25), MC.IsotropicState(3, 0.25)
    elseif name == "werner_state"
        return QET.werner_state(3, 0.2), MC.WernerState(3, 0.2)
    elseif name == "horodecki_3x3"
        return (QET.horodecki_state(0.3; dims=(3, 3)), MC.HorodeckiState(0.3, (3, 3)))
    elseif name == "horodecki_2x4"
        return (QET.horodecki_state(0.3; dims=(2, 4)), MC.HorodeckiState(0.3, (2, 4)))
    elseif name == "gisin_state"
        return QET.gisin_state(0.4, 0.7), MC.GisinState(0.4, 0.7)
    elseif name == "breuer_state"
        return QET.breuer_state(4, 0.35), MC.BreuerState(4, 0.35)
    elseif name == "brauer_states"
        return QET.brauer_states(2, 2), MC.BrauerStates(2, 2)
    elseif name == "chessboard_state"
        return (
            QET.chessboard_state(1, 2, 3, 4, 5, 6), MC.ChessboardState(1, 2, 3, 4, 5, 6)
        )
    end
    return error("unknown Tier B oracle fixture: $name")
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
    generated = joinpath(@__DIR__, "generated", "tier_b_oracle.json")
    committed = joinpath(@__DIR__, "fixtures", "tier_b_octave_11_3_qetlab_d858961.json")
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "Tier B oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_b_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "B" ||
    error("expected a Tier B fixture, got tier $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")

@testset "QETLAB Tier B differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        expected = fixture_matrix(fixture)
        native, compatibility = local_results(name)
        for (api, result) in
            ("native" => as_matrix(native), "MATLABCompat" => as_matrix(compatibility))
            @testset "$name/$api" begin
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
    "Verified $(length(payload.fixtures)) Tier B fixtures from ",
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
