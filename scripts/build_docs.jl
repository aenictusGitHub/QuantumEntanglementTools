#!/usr/bin/env julia

using Pkg

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DOCS_PROJECT = joinpath(REPOSITORY_ROOT, "docs", "Project.toml")
const DOCS_BUILD_SCRIPT = joinpath(REPOSITORY_ROOT, "docs", "make.jl")

"""
    build_docs()

Build the documentation with a fresh, Julia-version-appropriate dependency
resolution. The disposable environment is removed afterward, while Documenter
continues to write the rendered site to `docs/build`.
"""
function build_docs()
    isfile(DOCS_PROJECT) || error("documentation project not found at $DOCS_PROJECT")
    isfile(DOCS_BUILD_SCRIPT) ||
        error("documentation build script not found at $DOCS_BUILD_SCRIPT")

    mktempdir(; prefix="qet-docs-") do environment
        cp(DOCS_PROJECT, joinpath(environment, "Project.toml"))
        Pkg.activate(environment)
        Pkg.develop(PackageSpec(; path=REPOSITORY_ROOT))
        Pkg.instantiate()
        return include(DOCS_BUILD_SCRIPT)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ||
        throw(ArgumentError("usage: julia --startup-file=no scripts/build_docs.jl"))
    build_docs()
end
