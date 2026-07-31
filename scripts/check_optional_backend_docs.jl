module OptionalBackendDocs

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DOCS_SOURCE = joinpath(REPOSITORY_ROOT, "docs", "src")

function _julia_blocks(source::AbstractString)
    return (
        match.captures[1] for
        match in eachmatch(r"```(?:julia|@example[^\n]*)\n(.*?)```"s, source)
    )
end

"""
    check_optional_backend_docs(; io=stdout) -> Bool

Check that every documentation block constructing a `JuMPBackend` explicitly
loads JuMP in that same copyable block. This keeps optional-solver examples
independent of execution order and prevents an apparently configured backend
from leaving the package extension inactive.
"""
function check_optional_backend_docs(; io::IO=stdout)
    failures = String[]
    checked = 0
    for path in sort!(filter(endswith(".md"), readdir(DOCS_SOURCE; join=true)))
        source = read(path, String)
        for block in _julia_blocks(source)
            occursin("JuMPBackend(", block) || continue
            checked += 1
            occursin(r"(?m)^\s*using\s+(?:[^#\n]*,\s*)?JuMP(?:\s*,|$)", block) || push!(
                failures,
                "$(relpath(path, REPOSITORY_ROOT)) has a JuMPBackend block " *
                "without an explicit `using JuMP`",
            )
        end
    end
    if isempty(failures)
        println(io, "Optional-backend documentation check passed: $checked JuMP blocks.")
        return true
    end
    foreach(message -> println(io, "ERROR: ", message), failures)
    return false
end

end # module OptionalBackendDocs

if abspath(PROGRAM_FILE) == @__FILE__
    OptionalBackendDocs.check_optional_backend_docs() || exit(1)
end
