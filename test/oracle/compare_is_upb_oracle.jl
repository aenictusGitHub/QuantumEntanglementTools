using JSON3
using LinearAlgebra
using QuantumEntanglementTools
using SHA
using SparseArrays
using Test

const IsUPBOracleNative = QuantumEntanglementTools
const IsUPBOracleCompat = QuantumEntanglementTools.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const FIXTURE_NAMES = Set([
    "is_upb_tiles_flag",
    "is_upb_shifts_flag",
    "is_upb_extendible_flag",
    "is_upb_extendible_witness_left",
    "is_upb_extendible_witness_right",
    "is_upb_complete_basis_pinned_flag",
    "is_upb_nonorthogonal_pinned_flag",
    "is_upb_complex_flag",
    "is_upb_complex_witness_left",
    "is_upb_complex_witness_right",
])

function verify_digest(path::AbstractString)
    digest_path = path * ".sha256"
    isfile(digest_path) || error("missing fixture digest file: $digest_path")
    recorded = first(split(strip(read(digest_path, String))))
    actual = bytes2hex(sha256(read(path)))
    recorded == actual ||
        error("fixture digest mismatch: recorded $recorded, computed $actual")
    return actual
end

function fixture_array(fixture)
    dimensions = Tuple(Int.(fixture.dims))
    real_parts = fixture.real isa Number ? [fixture.real] : collect(fixture.real)
    imaginary_parts =
        fixture.imaginary isa Number ? [fixture.imaginary] : collect(fixture.imaginary)
    length(real_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent real data")
    length(imaginary_parts) == prod(dimensions) ||
        error("fixture $(fixture.name) has inconsistent imaginary data")
    values = if all(iszero, imaginary_parts)
        Float64.(real_parts)
    else
        complex.(Float64.(real_parts), Float64.(imaginary_parts))
    end
    return reshape(values, dimensions)
end

function tiles_factors()
    left = [
        1 1 0 0 1
        0 -1 0 1 1
        0 0 1 -1 1
    ]
    right = [
        1 0 0 1 1
        -1 0 1 0 1
        0 1 -1 0 1
    ]
    return left, right
end

function shifts_factors()
    first = [1 0 1 1; 0 1 -1 1]
    second = [1 1 0 1; 0 1 1 -1]
    third = [1 1 1 0; 0 -1 1 1]
    return first, second, third
end

committed_fixture = joinpath(@__DIR__, "fixtures", "is_upb_octave_11_3_qetlab_d858961.json")
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "is_upb_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "IsUPB oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_is_upb_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-is-upb" || error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
occursin("nonconjugating-transpose", String(metadata.reviewed_deviations)) ||
    error("fixture does not record the reviewed complex-witness deviation")
occursin("incompleteness", String(metadata.reviewed_deviations)) ||
    error("fixture does not record the reviewed complete-basis deviation")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixtures = Dict(
    String(fixture.name) => fixture_array(fixture) for fixture in payload.fixtures
)
Set(keys(fixtures)) == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

@testset "QETLAB IsUPB differential and corrected semantics" begin
    @test fixtures["is_upb_tiles_flag"][] == 1
    @test fixtures["is_upb_shifts_flag"][] == 1
    @test fixtures["is_upb_extendible_flag"][] == 0
    @test fixtures["is_upb_complete_basis_pinned_flag"][] == 1
    @test fixtures["is_upb_nonorthogonal_pinned_flag"][] == 1
    @test fixtures["is_upb_complex_flag"][] == 0

    tiles = IsUPBOracleNative.is_upb(tiles_factors(); normalization=:allow)
    shifts = IsUPBOracleNative.is_upb(shifts_factors(); normalization=:allow)
    @test tiles.status === :upb
    @test tiles.certificate_kind === :exact
    @test tiles.partitions_examined == 20
    @test shifts.status === :upb
    @test shifts.certificate_kind === :exact
    @test shifts.partitions_examined == 36
    @test IsUPBOracleCompat.IsUPB(tiles_factors()...) === true
    @test IsUPBOracleCompat.IsUPB(shifts_factors()...) === true

    extendible_left = [1 0; 0 1]
    extendible_right = [1 1; 0 0]
    extendible = IsUPBOracleNative.is_upb(extendible_left, extendible_right)
    @test extendible.status === :not_upb
    @test extendible.reason === :extension_witness
    @test IsUPBOracleCompat.IsUPB(extendible_left, extendible_right) === false
    @test all(
        state -> iszero(
            prod(
                dot(
                    extendible.witness_factors[party],
                    (extendible_left, extendible_right)[party][:, state],
                ) for party in 1:2
            ),
        ),
        1:2,
    )

    pinned_extendible_witness = (
        vec(fixtures["is_upb_extendible_witness_left"]),
        vec(fixtures["is_upb_extendible_witness_right"]),
    )
    @test all(
        state -> iszero(
            prod(
                dot(
                    pinned_extendible_witness[party],
                    (extendible_left, extendible_right)[party][:, state],
                ) for party in 1:2
            ),
        ),
        1:2,
    )

    complete_left = [1 1 0 0; 0 0 1 1]
    complete_right = [1 0 1 0; 0 1 0 1]
    complete = IsUPBOracleNative.is_upb(complete_left, complete_right)
    @test complete.status === :not_upb
    @test complete.reason === :complete_basis
    @test IsUPBOracleCompat.IsUPB(complete_left, complete_right) === false

    generic_left = [
        1 1 1 1 1
        0 1 2 3 4
        0 1 4 9 16
    ]
    generic_right = [
        1 1 1 1 1
        0 2 4 6 8
        0 4 16 36 64
    ]
    nonorthogonal = IsUPBOracleNative.is_upb(
        generic_left, generic_right; normalization=:allow
    )
    @test nonorthogonal.status === :not_upb
    @test nonorthogonal.reason === :not_orthogonal
    @test IsUPBOracleCompat.IsUPB(generic_left, generic_right) === false
    @test Float64(metadata.nonorthogonal_first_overlap) > 0

    complex_left = ComplexF64[1, 0]
    complex_right = ComplexF64[1, im]
    pinned_complex_witness = (
        vec(fixtures["is_upb_complex_witness_left"]),
        vec(fixtures["is_upb_complex_witness_right"]),
    )
    pinned_hilbert_residual = abs(
        dot(pinned_complex_witness[1], complex_left) *
        dot(pinned_complex_witness[2], complex_right),
    )
    pinned_bilinear_residual = abs(
        transpose(pinned_complex_witness[1]) *
        complex_left *
        (transpose(pinned_complex_witness[2]) * complex_right),
    )
    @test pinned_hilbert_residual ≈ Float64(metadata.complex_witness_hilbert_residual)
    @test pinned_hilbert_residual > 1
    @test pinned_bilinear_residual ≈ Float64(metadata.complex_witness_bilinear_residual) atol =
        1e-14
    @test abs(
        dot(conj(pinned_complex_witness[1]), complex_left) *
        dot(conj(pinned_complex_witness[2]), complex_right),
    ) ≈ Float64(metadata.conjugated_witness_hilbert_residual) atol = 1e-14

    corrected = IsUPBOracleNative.is_upb(
        reshape(Complex{Int}[1, 0], 2, 1),
        reshape(Complex{Int}[1, im], 2, 1);
        normalization=:allow,
    )
    @test corrected.status === :not_upb
    @test corrected.reason === :extension_witness
    @test iszero(
        dot(corrected.witness_factors[1], complex_left) *
        dot(corrected.witness_factors[2], complex_right),
    )
end

println(
    "IsUPB oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
