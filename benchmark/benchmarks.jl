#!/usr/bin/env julia

using BenchmarkTools
using Dates
using InteractiveUtils
using LinearAlgebra
using Random
using SparseArrays
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)

using QuantumEntanglementTools

function parse_options(args)
    quick = false
    save_results = true
    output_prefix = nothing
    for arg in args
        if arg == "--quick"
            quick = true
        elseif arg == "--no-save"
            save_results = false
        elseif startswith(arg, "--output=")
            output_prefix = split(arg, "="; limit=2)[2]
        elseif arg in ("-h", "--help")
            println(
                """
                Usage: julia --project=benchmark benchmark/benchmarks.jl [options]

                  --quick          short smoke benchmark
                  --no-save        print results without writing raw artifacts
                  --output=PREFIX  artifact prefix (default: benchmark/results/local/<timestamp>)
                """,
            )
            exit(0)
        else
            error("unknown option: $arg")
        end
    end
    return (; quick, save_results, output_prefix)
end

function benchmark_suite()
    rng = MersenneTwister(0x5145544c41424a55)
    suite = BenchmarkGroup()

    dense_dims = (2, 3, 4, 2)
    dense_dimension = prod(dense_dims)
    dense_state = randn(rng, ComplexF64, dense_dimension)
    dense_operator = randn(rng, ComplexF64, dense_dimension, dense_dimension)
    trace_plan = PartialTracePlan(dense_dims, (2, 4))
    transpose_plan = PartialTransposePlan(dense_dims, (2, 4))
    permutation_plan = SubsystemPermutationPlan(dense_dims, (4, 2, 1, 3))

    suite["partial_trace/dense_pure_plan_reuse"] = @benchmarkable partial_trace(
        $dense_state, $trace_plan
    )
    suite["partial_trace/dense_matrix_plan_reuse"] = @benchmarkable partial_trace(
        $dense_operator, $trace_plan
    )
    suite["partial_trace/plan_construction"] = @benchmarkable PartialTracePlan(
        $dense_dims, (2, 4)
    )
    suite["partial_transpose/dense_plan_reuse"] = @benchmarkable partial_transpose(
        $dense_operator, $transpose_plan
    )
    suite["permutation/dense_matrix_plan_reuse"] = @benchmarkable permute_subsystems(
        $dense_operator, $permutation_plan
    )

    sparse_dims = (4, 4, 4, 4)
    sparse_dimension = prod(sparse_dims)
    sparse_operator = sprand(rng, ComplexF64, sparse_dimension, sparse_dimension, 0.005)
    sparse_trace_plan = PartialTracePlan(sparse_dims, (2, 4))
    sparse_transpose_plan = PartialTransposePlan(sparse_dims, (2, 4))
    suite["partial_trace/sparse_matrix_plan_reuse"] = @benchmarkable partial_trace(
        $sparse_operator, $sparse_trace_plan
    )
    suite["partial_transpose/sparse_plan_reuse"] = @benchmarkable partial_transpose(
        $sparse_operator, $sparse_transpose_plan
    )

    bipartite_dims = (8, 8)
    bipartite_operator = randn(rng, ComplexF64, 64, 64)
    realignment_plan = RealignmentPlan(bipartite_dims)
    suite["realignment/dense_plan_reuse"] = @benchmarkable realign(
        $bipartite_operator, $realignment_plan
    )

    factor_a = randn(rng, ComplexF64, 8, 8)
    factor_b = randn(rng, ComplexF64, 6, 6)
    suite["tensor_product/dense_two_factor"] = @benchmarkable tensor_product(
        $factor_a, $factor_b
    )
    suite["projector/symmetric_d4_p4_sparse"] = @benchmarkable symmetric_projector(4, 4)

    suite["operators/fourier_d32"] = @benchmarkable fourier_matrix(32)
    suite["states/werner_d8_sparse"] = @benchmarkable werner_state(8, 0.2)

    density_rng = MersenneTwister(0x5145544c41424431)
    unitary_rng = MersenneTwister(0x5145544c41424432)
    povm_rng = MersenneTwister(0x5145544c41424433)
    suite["random/density_d32_rank16"] = @benchmarkable random_density_matrix(
        $density_rng, 32; rank=16
    )
    suite["random/unitary_d32"] = @benchmarkable random_unitary($unitary_rng, 32)
    suite["random/povm_d8_outcomes4"] = @benchmarkable random_povm($povm_rng, 8, 4)

    channel_input = randn(rng, ComplexF64, 16, 16)
    channel = dephasing_channel(16, 0.25)
    channel_kraus = kraus_representation(channel)
    suite["channels/apply_dephasing_d16"] = @benchmarkable apply_channel(
        $channel_input, $channel
    )
    suite["channels/kraus_to_choi_dephasing_d16"] = @benchmarkable choi_representation(
        $channel_kraus
    )
    suite["channels/kraus_to_superoperator_dephasing_d16"] = @benchmarkable superoperator_representation(
        $channel_kraus
    )

    tier_d_operator = randn(rng, ComplexF64, 64, 48)
    tier_d_factor = randn(rng, ComplexF64, 64, 64)
    tier_d_density = tier_d_factor * adjoint(tier_d_factor)
    tier_d_density ./= real(tr(tier_d_density))
    tier_d_other_factor = randn(rng, ComplexF64, 64, 64)
    tier_d_other_density = tier_d_other_factor * adjoint(tier_d_other_factor)
    tier_d_other_density ./= real(tr(tier_d_other_density))
    fidelity_factor = randn(rng, ComplexF64, 32, 32)
    fidelity_density = fidelity_factor * adjoint(fidelity_factor)
    fidelity_density ./= real(tr(fidelity_density))
    fidelity_other_factor = randn(rng, ComplexF64, 32, 32)
    fidelity_other_density = fidelity_other_factor * adjoint(fidelity_other_factor)
    fidelity_other_density ./= real(tr(fidelity_other_density))
    schmidt_vector = randn(rng, ComplexF64, 32 * 32)
    bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
    concurrence_density =
        0.7 * (bell * adjoint(bell)) + 0.3 * Matrix{ComplexF64}(I, 4, 4) / 4

    suite["measures/trace_norm_64x48"] = @benchmarkable trace_norm($tier_d_operator)
    suite["measures/schatten_p3_64x48"] = @benchmarkable schatten_norm($tier_d_operator, 3)
    suite["measures/ky_fan_k16_64x48"] = @benchmarkable ky_fan_norm($tier_d_operator, 16)
    suite["measures/purity_density_d64"] = @benchmarkable purity($tier_d_density)
    suite["measures/entropy_density_d64"] = @benchmarkable von_neumann_entropy(
        $tier_d_density; base=2
    )
    suite["measures/fidelity_density_d32"] = @benchmarkable fidelity(
        $fidelity_density, $fidelity_other_density
    )
    suite["measures/trace_distance_density_d64"] = @benchmarkable trace_distance(
        $tier_d_density, $tier_d_other_density
    )
    suite["measures/concurrence_mixed_two_qubit"] = @benchmarkable concurrence(
        $concurrence_density
    )
    suite["entanglement/negativity_8x8"] = @benchmarkable negativity(
        $tier_d_density, (8, 8)
    )
    suite["entanglement/logarithmic_negativity_8x8"] = @benchmarkable logarithmic_negativity(
        $tier_d_density, (8, 8); base=2
    )
    suite["entanglement/schmidt_decomposition_32x32"] = @benchmarkable schmidt_decomposition(
        $schmidt_vector, (32, 32)
    )
    suite["entanglement/schmidt_rank_32x32"] = @benchmarkable schmidt_rank(
        $schmidt_vector, (32, 32)
    )
    suite["criteria/ppt_8x8"] = @benchmarkable ppt_criterion($tier_d_density, (8, 8))
    suite["criteria/realignment_8x8"] = @benchmarkable realignment_criterion(
        $tier_d_density, (8, 8)
    )
    suite["criteria/reduction_8x8"] = @benchmarkable reduction_criterion(
        $tier_d_density, (8, 8)
    )

    coherence_plus = fill(ComplexF64(inv(sqrt(4096))), 4096)
    coherence_basis = fourier_matrix(64)
    coherence_state = randn(rng, ComplexF64, 64)
    coherence_state ./= norm(coherence_state)
    suite["coherence/l1_pure_d4096"] = @benchmarkable l1_coherence($coherence_plus)
    suite["coherence/relative_entropy_pure_d4096"] = @benchmarkable relative_entropy_coherence(
        $coherence_plus; base=2
    )
    suite["coherence/rank_basis_transform_d64"] = @benchmarkable coherence_rank(
        $coherence_state; basis=($coherence_basis)
    )

    product_operator = tensor_product(
        randn(rng, ComplexF64, 4, 4),
        randn(rng, ComplexF64, 8, 8),
        randn(rng, ComplexF64, 2, 2),
    )
    suite["product/operator_schmidt_decomposition_8x8"] = @benchmarkable operator_schmidt_decomposition(
        $tier_d_density, (8, 8)
    )
    suite["product/operator_analysis_4x8x2"] = @benchmarkable is_product_operator(
        $product_operator, (4, 8, 2)
    )
    suite["product/formation_mixed_two_qubit"] = @benchmarkable entanglement_of_formation(
        $concurrence_density, (2, 2)
    )

    majorization_first = randn(rng, 32, 32)
    majorization_second = randn(rng, 32, 32)
    compound_input = randn(rng, 8, 8)
    additive_compound_input = randn(rng, 16, 16)
    suite["matrix_analysis/majorizes_singular_values_32x32"] = @benchmarkable majorizes(
        $majorization_first, $majorization_second
    )
    suite["matrix_analysis/compound_dense_8_order3"] = @benchmarkable compound_matrix(
        $compound_input, 3
    )
    suite["matrix_analysis/additive_compound_dense_16_order2"] = @benchmarkable additive_compound_matrix(
        $additive_compound_input, 2
    )

    return suite
end

function command_output(args...)
    try
        stdout = IOBuffer()
        process = run(pipeline(ignorestatus(Cmd(collect(args))); stdout, stderr=devnull))
        success(process) || return "unavailable"
        return chomp(String(take!(stdout)))
    catch
        return "unavailable"
    end
end

function result_rows(results::BenchmarkGroup)
    rows = Vector{Dict{String,Any}}()
    for name in sort!(collect(keys(results)); by=string)
        estimate = minimum(results[name])
        push!(
            rows,
            Dict(
                "name" => string(name),
                "minimum_time_ns" => estimate.time,
                "minimum_gctime_ns" => estimate.gctime,
                "memory_bytes" => estimate.memory,
                "allocations" => estimate.allocs,
                "samples" => length(results[name].times),
            ),
        )
    end
    return rows
end

function print_results(rows)
    println("case\tminimum time (ns)\tmemory (bytes)\tallocations\tsamples")
    for row in rows
        println(
            row["name"],
            '\t',
            row["minimum_time_ns"],
            '\t',
            row["memory_bytes"],
            '\t',
            row["allocations"],
            '\t',
            row["samples"],
        )
    end
end

function metadata(rows, options)
    return Dict(
        "schema_version" => 1,
        "recorded_at_utc" => string(now(UTC)),
        "quick" => options.quick,
        "repository_commit" =>
            command_output("git", "-C", REPOSITORY_ROOT, "rev-parse", "HEAD"),
        "repository_status" =>
            command_output("git", "-C", REPOSITORY_ROOT, "status", "--short"),
        "julia_version" => string(VERSION),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "cpu_name" => Sys.CPU_NAME,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_configuration" => sprint(show, BLAS.get_config()),
        "benchmarktools_version" => string(Base.pkgversion(BenchmarkTools)),
        "rng" => "MersenneTwister",
        "rng_seed_hex" => "0x5145544c41424a55",
        "cases" => rows,
    )
end

function main(args)
    options = parse_options(args)
    seconds = options.quick ? 0.2 : 2.0
    samples = options.quick ? 20 : 10_000
    suite = benchmark_suite()
    results = run(suite; verbose=(!options.quick), seconds, samples, evals=1)
    rows = result_rows(results)
    print_results(rows)

    if options.save_results
        default_stamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS")
        prefix = something(
            options.output_prefix, joinpath(@__DIR__, "results", "local", default_stamp)
        )
        mkpath(dirname(prefix))
        raw_path = prefix * ".json"
        metadata_path = prefix * ".toml"
        BenchmarkTools.save(raw_path, results)
        open(metadata_path, "w") do io
            return TOML.print(io, metadata(rows, options); sorted=true)
        end
        println("raw results: ", relpath(raw_path, REPOSITORY_ROOT))
        println("metadata: ", relpath(metadata_path, REPOSITORY_ROOT))
    end
end

main(ARGS)
