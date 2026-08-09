module ReleaseGateChecks

using Test

include(joinpath(@__DIR__, "..", "scripts", "check_release.jl"))

const DEVELOPMENT_CHANGELOG = """
# Changelog

## [Unreleased]
"""
const DATED_CHANGELOG = DEVELOPMENT_CHANGELOG * """

## [0.1.0] - 2026-07-31
"""
const DEVELOPMENT_CITATION = """
cff-version: 1.2.0
version: 0.1.0
repository-code: https://github.com/aenictusGitHub/QuantumEntanglementTools
"""
const DATED_CITATION = DEVELOPMENT_CITATION * "date-released: 2026-07-31\n"

@testset "Release-check option modes" begin
    withenv("GITHUB_REF_TYPE" => "", "GITHUB_REF_NAME" => "") do
        development = parse_options(String[])
        @test !development.release_candidate
        @test isnothing(development.tag)

        candidate = parse_options(["--release-candidate", "--treeish", "candidate"])
        @test candidate.release_candidate
        @test candidate.treeish == "candidate"
        @test isnothing(candidate.tag)

        tagged = parse_options(["--tag=v0.1.0"])
        @test tagged.tag == "v0.1.0"
        @test !tagged.release_candidate
        @test_throws ErrorException parse_options([
            "--release-candidate", "--tag", "v0.1.0"
        ])
    end
end

@testset "Development and release-candidate metadata" begin
    failures = String[]
    validate_release_metadata!(
        failures, v"0.1.0", DEVELOPMENT_CITATION, DEVELOPMENT_CHANGELOG; dated=false
    )
    @test isempty(failures)

    failures = String[]
    validate_release_metadata!(
        failures, v"0.1.0", DATED_CITATION, DATED_CHANGELOG; dated=true
    )
    @test isempty(failures)

    failures = String[]
    validate_release_metadata!(
        failures, v"0.1.0", DEVELOPMENT_CITATION, DEVELOPMENT_CHANGELOG; dated=true
    )
    @test "CHANGELOG.md has no dated release heading for 0.1.0" in failures
    @test "CITATION.cff has no release date" in failures

    failures = String[]
    validate_release_metadata!(
        failures, v"0.1.0", DATED_CITATION, DATED_CHANGELOG; dated=false
    )
    @test any(contains("must not contain a dated 0.1.0 release heading"), failures)
    @test "unreleased CITATION.cff must not contain date-released" in failures

    failures = String[]
    invalid_citation = replace(DATED_CITATION, "2026-07-31" => "2026-02-31")
    validate_release_metadata!(
        failures, v"0.1.0", invalid_citation, DATED_CHANGELOG; dated=true
    )
    @test "CITATION.cff release date is not a valid ISO calendar date" in failures
    @test "CHANGELOG.md and CITATION.cff release dates disagree" in failures
end

@testset "Tagged release copy" begin
    clean_documents = Dict(
        "AGENTS.md" => "Release candidate maintenance instructions.\n",
        "CONTRIBUTING.md" => "Release candidate contribution instructions.\n",
        "README.md" => "Install with rev=\"v0.1.0\".\n",
        "SECURITY.md" => "Security support starts with version 0.1.0.\n",
        "CITATION.cff" => "version: 0.1.0\ndate-released: 2026-07-31\n",
        "CITATION.bib" => "note = {QuantumEntanglementTools.jl 0.1.0}\n",
        "docs/CONVERGENCE_AUDIT.md" => "Release candidate audit.\n",
        "docs/PORTING_STATUS.md" => "Release candidate status.\n",
        "docs/SESSION_HANDOFF.md" => "Release candidate handoff.\n",
        "docs/VALIDATION_REPORT.md" => "Release candidate validation.\n",
        "docs/src/getting_started.md" => "Pkg.add(; rev=\"v0.1.0\")\n",
    )
    failures = String[]
    validate_release_copy!(failures, v"0.1.0", clean_documents)
    @test isempty(failures)

    stale_documents = copy(clean_documents)
    stale_documents["README.md"] = "This unreleased version uses rev=\"main\".\n"
    stale_documents["SECURITY.md"] = "No version has been published.\n"
    failures = String[]
    validate_release_copy!(failures, v"0.1.0", stale_documents)
    @test any(contains("README.md retains pre-release wording"), failures)
    @test any(contains("SECURITY.md retains pre-release wording"), failures)
    @test any(contains("README.md must contain the stable installation"), failures)

    stale_documents = copy(clean_documents)
    stale_documents["AGENTS.md"] = "This package is unreleased pre-1.0 software.\n"
    failures = String[]
    validate_release_copy!(failures, v"0.1.0", stale_documents)
    @test any(contains("AGENTS.md retains pre-release wording"), failures)

    historical_documents = copy(clean_documents)
    historical_documents["CONTRIBUTING.md"] =
        join(fill("current release guidance", 120), '\n') *
        "\nUpdate CHANGELOG.md under Unreleased for later work.\n"
    failures = String[]
    validate_release_copy!(failures, v"0.1.0", historical_documents)
    @test isempty(failures)

    missing_documents = copy(clean_documents)
    delete!(missing_documents, "CITATION.bib")
    failures = String[]
    validate_release_copy!(failures, v"0.1.0", missing_documents)
    @test "release-copy check has no CITATION.bib" in failures
end

@testset "Generated completion-artifact subprocess" begin
    mktempdir() do candidate_root
        missing = check_completion_artifacts(candidate_root)
        @test !missing.success
        @test occursin("has no scripts/build_qetlab_completion_plan.jl", missing.output)

        scripts_root = joinpath(candidate_root, "scripts")
        mkpath(scripts_root)
        generator = joinpath(scripts_root, "build_qetlab_completion_plan.jl")
        open(generator, "w") do io
            return write(
                io, "ARGS == [\"--check\"] || exit(2)\nprintln(\"artifacts current\")\n"
            )
        end
        current = check_completion_artifacts(candidate_root)
        @test current.success
        @test current.output == "artifacts current"

        open(generator, "w") do io
            return write(io, "println(stderr, \"stale generated artifact\")\nexit(1)\n")
        end
        stale = check_completion_artifacts(candidate_root)
        @test !stale.success
        @test stale.output == "stale generated artifact"
    end
end

end
