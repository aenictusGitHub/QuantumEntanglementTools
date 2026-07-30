using JSON3
using QuantumEntanglementTools
using SHA
using Test

const QET = QuantumEntanglementTools
const Compat = QET.MATLABCompat
const PINNED_QETLAB_COMMIT = "d8589610f00cff106537268dee2e2a1153f3a601"
const FIXTURE_NAMES = Set([
    "entangled_subspace_equal_default",
    "entangled_subspace_rectangular_wide",
    "entangled_subspace_rectangular_r2",
    "entangled_subspace_prefix",
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

function fixture_matrix(fixture)
    dimensions = Tuple(Int.(fixture.dims))
    length(dimensions) == 2 ||
        error("fixture $(fixture.name) does not record two dimensions")
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

committed_fixture = joinpath(
    @__DIR__, "fixtures", "entangled_subspace_octave_11_3_qetlab_d858961.json"
)
fixture_path = if isempty(ARGS)
    generated = joinpath(@__DIR__, "generated", "entangled_subspace_oracle.json")
    isfile(generated) ? generated : committed_fixture
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "EntangledSubspace oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_entangled_subspace_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.tier) == "WP2-entangled-subspace" ||
    error("unexpected fixture tier: $(metadata.tier)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Bool(metadata.source_free_fixture) || error("fixture is not marked source-free")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("fixture count does not match metadata")
occursin("slowest-varying", String(metadata.vectorization_order)) ||
    error("fixture does not record subsystem vectorization order")

if normpath(fixture_path) == normpath(committed_fixture)
    String(metadata.engine) == "Octave" ||
        error("the committed fixture must record Octave as its engine")
    startswith(String(metadata.engine_version), "11.3") ||
        error("the committed fixture must come from Octave 11.3")
end

fixture_names = Set(String(fixture.name) for fixture in payload.fixtures)
fixture_names == FIXTURE_NAMES ||
    error("fixture names do not match the reviewed comparator set")

function fixture_arguments(name::String)
    name == "entangled_subspace_equal_default" && return (4, 3, 1)
    name == "entangled_subspace_rectangular_wide" && return (3, [2, 4], 1)
    name == "entangled_subspace_rectangular_r2" && return (2, [3, 4], 2)
    name == "entangled_subspace_prefix" && return (2, [3, 3], 1)
    return error("unreviewed EntangledSubspace fixture: $name")
end

@testset "QETLAB EntangledSubspace differential fixture" begin
    for fixture in payload.fixtures
        name = String(fixture.name)
        subspace_dimension, local_dims, r = fixture_arguments(name)
        expected = fixture_matrix(fixture)
        native = QET.entangled_subspace(subspace_dimension, local_dims; r)
        compatibility = Compat.EntangledSubspace(subspace_dimension, local_dims, r)
        @test String(fixture.comparison) == "exact"
        @test size(native) == size(expected)
        @test Matrix(native) == expected
        @test Matrix(compatibility) == expected
    end
end

println(
    "EntangledSubspace oracle comparison passed for ",
    length(payload.fixtures),
    " fixtures from ",
    metadata.engine,
    " ",
    metadata.engine_version,
    " (SHA-256 ",
    digest,
    ").",
)
