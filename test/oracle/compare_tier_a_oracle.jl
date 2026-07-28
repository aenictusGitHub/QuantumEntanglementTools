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
    values = if all(iszero, imaginary_parts)
        real_parts
    else
        complex.(real_parts, imaginary_parts)
    end
    return reshape(values, rows, columns)
end

as_matrix(value::AbstractMatrix) = Matrix(value)
as_matrix(value::AbstractVector) = reshape(collect(value), :, 1)
as_matrix(value::Number) = reshape([value], 1, 1)

function local_results(name::AbstractString)
    a2 = [1 2; 3 4]
    b2 = [0 1; 1 0]
    a3 = [0 1 0; 1 0 1; 0 1 0]
    dims = (2, 3, 2)
    permutation = (3, 1, 2)
    vector12 = collect(1:12)
    matrix12 = reshape(collect(1:144), 12, 12)

    if name == "tensor"
        return (QET.tensor_product(a2, b2), MC.Tensor(a2, b2))
    elseif name == "tensor_sum"
        identity2 = Matrix{Int}(I, 2, 2)
        return (
            QET.tensor_sum(identity2, a2; weights=[2, -1]),
            MC.TensorSum([2, -1], identity2, a2),
        )
    elseif name == "kronecker_sum"
        return (QET.kronecker_sum(a2, a3), MC.KroneckerSum(a2, a3))
    elseif name == "permute_systems"
        return (
            QET.permute_subsystems(vector12, dims; permutation=permutation),
            MC.PermuteSystems(vector12, permutation, dims),
        )
    elseif name == "permutation_operator"
        return (
            QET.permutation_operator(dims, permutation; sparse_output=false, T=Int),
            MC.PermutationOperator(dims, permutation),
        )
    elseif name == "swap"
        return (QET.swap_subsystems(vector12, dims, 1, 3), MC.Swap(vector12, (1, 3), dims))
    elseif name == "swap_operator"
        return (
            QET.swap_operator((2, 3); sparse_output=false, T=Int), MC.SwapOperator((2, 3))
        )
    elseif name == "partial_trace"
        return (
            QET.partial_trace(matrix12, dims; trace_out=(1, 3)),
            MC.PartialTrace(matrix12, (1, 3), dims),
        )
    elseif name == "partial_transpose"
        return (
            QET.partial_transpose(matrix12, dims; systems=(1, 3)),
            MC.PartialTranspose(matrix12, (1, 3), dims),
        )
    elseif name == "realignment"
        matrix6 = reshape(collect(1:36), 6, 6)
        return (QET.realign(matrix6, (2, 3)), MC.Realignment(matrix6, (2, 3)))
    elseif name == "complex_partial_transpose"
        complex_matrix = reshape(ComplexF64.(1:16) .+ im .* ComplexF64.(16:-1:1), 4, 4)
        return (
            QET.partial_transpose(complex_matrix, (2, 2); systems=(1,)),
            MC.PartialTranspose(complex_matrix, 1, (2, 2)),
        )
    elseif name == "symmetric_projection"
        return (
            QET.symmetric_projector(3, 2; sparse_output=false),
            MC.SymmetricProjection(3, 2, 0),
        )
    elseif name == "antisymmetric_projection"
        return (
            QET.antisymmetric_projector(3, 2; sparse_output=false),
            MC.AntisymmetricProjection(3, 2, 0),
        )
    end
    return error("unknown oracle fixture: $name")
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
    generated = joinpath(@__DIR__, "generated", "tier_a_oracle.json")
    committed = joinpath(@__DIR__, "fixtures", "tier_a_octave_11_3_qetlab_d858961.json")
    isfile(generated) ? generated : committed
else
    abspath(ARGS[1])
end
isfile(fixture_path) || error(
    "oracle fixture not found at $fixture_path; " *
    "run scripts/matlab_oracle/run_tier_a_oracle.sh first",
)

digest = verify_digest(fixture_path)
payload = JSON3.read(read(fixture_path, String))
metadata = payload.metadata
String(metadata.schema) == "quantum-entanglement-tools-oracle-v1" ||
    error("unsupported oracle schema: $(metadata.schema)")
String(metadata.qetlab_commit) == PINNED_QETLAB_COMMIT ||
    error("fixture was not generated from the pinned QETLAB revision")
Int(metadata.fixture_count) == length(payload.fixtures) ||
    error("oracle fixture count does not match metadata")

@testset "QETLAB Tier A differential fixture" begin
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
    "Verified $(length(payload.fixtures)) fixtures from ",
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
