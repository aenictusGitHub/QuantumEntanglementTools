using JuliaFormatter

repository_root = normpath(joinpath(@__DIR__, ".."))
paths = [
    joinpath(repository_root, path) for
    path in ("src", "ext", "test", "tutorials", "benchmark", "scripts", "docs", "quality")
]
already_formatted = JuliaFormatter.format(
    paths; overwrite=false, throw_on_error=true, verbose=true
)
already_formatted ||
    error("Julia source is not formatted; run JuliaFormatter.format on the listed paths")
