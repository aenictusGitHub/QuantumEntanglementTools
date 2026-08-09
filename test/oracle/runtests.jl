#!/usr/bin/env julia

using SHA

const ORACLE_ROOT = @__DIR__
const FIXTURE_ROOT = joinpath(ORACLE_ROOT, "fixtures")

# Keep this map explicit. The manifest audit below fails when a comparator or
# committed fixture is added without a reviewed pairing here.
const ORACLE_CASES = (
    (
        script="compare_absolute_ppt_oracle.jl",
        fixture="absolute_ppt_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_channel_optimization_oracle.jl",
        fixture="channel_optimization_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_coherence_optimization_oracle.jl",
        fixture="coherence_optimization_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_copositivity_clique_oracle.jl",
        fixture="copositivity_clique_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_entangled_subspace_oracle.jl",
        fixture="entangled_subspace_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_entangling_gate_oracle.jl",
        fixture="entangling_gate_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_induced_schatten_oracle.jl",
        fixture="induced_schatten_octave_11_3_qetlab_d858961.json",
    ),
    (script="compare_is_upb_oracle.jl", fixture="is_upb_octave_11_3_qetlab_d858961.json"),
    (
        script="compare_minimum_upb_size_oracle.jl",
        fixture="minimum_upb_size_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_multipartite_werner_oracle.jl",
        fixture="multipartite_werner_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_nonlocal_games_oracle.jl",
        fixture="nonlocal_games_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_parallel_repetition_oracle.jl",
        fixture="parallel_repetition_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_random_superoperator_oracle.jl",
        fixture="random_superoperator_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_robk_coherence_oracle.jl",
        fixture="robk_coherence_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_separability_oracle.jl",
        fixture="separability_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_sk_norms_oracle.jl",
        fixture="sk_norms_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_state_discrimination_oracle.jl",
        fixture="state_discrimination_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_symmetric_extensions_oracle.jl",
        fixture="symmetric_extensions_octave_11_3_qetlab_d858961.json",
    ),
    (script="compare_tier_a_oracle.jl", fixture="tier_a_octave_11_3_qetlab_d858961.json"),
    (script="compare_tier_b_oracle.jl", fixture="tier_b_octave_11_3_qetlab_d858961.json"),
    (script="compare_tier_c_oracle.jl", fixture="tier_c_octave_11_3_qetlab_d858961.json"),
    (script="compare_tier_d_oracle.jl", fixture="tier_d_octave_11_3_qetlab_d858961.json"),
    (
        script="compare_tier_e_coherence_oracle.jl",
        fixture="tier_e_coherence_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_tier_e_matrix_analysis_oracle.jl",
        fixture="tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json",
    ),
    (
        script="compare_tier_e_product_oracle.jl",
        fixture="tier_e_product_octave_11_3_qetlab_d858961.json",
    ),
    (script="compare_twirl_oracle.jl", fixture="twirl_octave_11_3_qetlab_d858961.json"),
    (
        script="compare_upb_catalog_oracle.jl",
        fixture="upb_catalog_octave_11_3_qetlab_d858961.json",
    ),
)

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia --project=test/oracle test/oracle/runtests.jl [--check-manifest]

With no option, validate the explicit comparator/fixture manifest and every
fixture digest, run all comparator processes in deterministic order with an
explicit committed fixture path, then verify the fixture and sidecar bytes
again. --check-manifest performs only the manifest and digest checks.
""",
    )
end

file_digest(path::AbstractString) = bytes2hex(sha256(read(path)))

function verify_fixture_digest(fixture_name::AbstractString)
    fixture_path = joinpath(FIXTURE_ROOT, fixture_name)
    digest_path = fixture_path * ".sha256"
    isfile(fixture_path) || error("missing committed oracle fixture: $fixture_path")
    isfile(digest_path) || error("missing fixture digest sidecar: $digest_path")

    fields = split(strip(read(digest_path, String)))
    length(fields) == 2 ||
        error("fixture digest sidecar must contain a SHA-256 and filename: $digest_path")
    recorded = lowercase(fields[1])
    occursin(r"^[0-9a-f]{64}$", recorded) ||
        error("fixture digest sidecar contains a malformed SHA-256: $digest_path")
    recorded_name = lstrip(fields[2], '*')
    recorded_name == fixture_name || error(
        "fixture digest sidecar names $(repr(recorded_name)), expected $(repr(fixture_name))",
    )
    actual = file_digest(fixture_path)
    recorded == actual || error(
        "fixture digest mismatch for $fixture_name: recorded $recorded, computed $actual",
    )
    return (fixture=actual, sidecar=file_digest(digest_path))
end

function validate_oracle_manifest()
    mapped_scripts = Set(case.script for case in ORACLE_CASES)
    mapped_fixtures = Set(case.fixture for case in ORACLE_CASES)
    length(mapped_scripts) == length(ORACLE_CASES) ||
        error("oracle manifest maps a comparator more than once")
    length(mapped_fixtures) == length(ORACLE_CASES) ||
        error("oracle manifest maps a fixture more than once")

    actual_scripts = Set(
        file for file in readdir(ORACLE_ROOT) if
        startswith(file, "compare_") && endswith(file, "_oracle.jl")
    )
    actual_fixtures = Set(file for file in readdir(FIXTURE_ROOT) if endswith(file, ".json"))
    actual_sidecars = Set(
        chop(file; tail=7) for
        file in readdir(FIXTURE_ROOT) if endswith(file, ".json.sha256")
    )

    actual_scripts == mapped_scripts || error(
        "oracle comparator manifest drift: missing mappings=$(sort!(collect(setdiff(actual_scripts, mapped_scripts)))); " *
        "stale mappings=$(sort!(collect(setdiff(mapped_scripts, actual_scripts))))",
    )
    actual_fixtures == mapped_fixtures || error(
        "oracle fixture manifest drift: missing mappings=$(sort!(collect(setdiff(actual_fixtures, mapped_fixtures)))); " *
        "stale mappings=$(sort!(collect(setdiff(mapped_fixtures, actual_fixtures))))",
    )
    actual_sidecars == mapped_fixtures || error(
        "oracle digest-sidecar manifest drift: missing sidecars=$(sort!(collect(setdiff(mapped_fixtures, actual_sidecars)))); " *
        "unexpected sidecars=$(sort!(collect(setdiff(actual_sidecars, mapped_fixtures))))",
    )

    return Dict(
        case.fixture => verify_fixture_digest(case.fixture) for case in ORACLE_CASES
    )
end

function run_comparators()
    before = validate_oracle_manifest()
    for (index, case) in enumerate(ORACLE_CASES)
        script_path = joinpath(ORACLE_ROOT, case.script)
        fixture_path = joinpath(FIXTURE_ROOT, case.fixture)
        println("[", index, "/", length(ORACLE_CASES), "] ", case.script)
        command = `$(Base.julia_cmd()) --startup-file=no --project=$ORACLE_ROOT $script_path $fixture_path`
        run(command)
    end
    after = validate_oracle_manifest()
    before == after ||
        error("an oracle fixture or digest sidecar changed during comparison")
    println(
        "all ",
        length(ORACLE_CASES),
        " committed oracle comparators passed with unchanged verified fixtures",
    )
    return nothing
end

function main(args)
    if isempty(args)
        run_comparators()
    elseif args == ["--check-manifest"]
        validate_oracle_manifest()
        println(
            "oracle manifest and SHA-256 sidecars passed for ",
            length(ORACLE_CASES),
            " comparator/fixture pairs",
        )
    elseif args in (["-h"], ["--help"])
        usage()
    else
        usage(stderr)
        error("unknown arguments: $(join(args, ' '))")
    end
    return nothing
end

main(ARGS)
