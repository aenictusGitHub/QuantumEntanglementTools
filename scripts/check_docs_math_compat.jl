#!/usr/bin/env julia

# Keep documentation equations readable in both GitHub's repository preview
# and Documenter. GitHub does not interpret Documenter's double-backtick
# inline-math extension, while both render `$...$` inline math and fenced
# `math` display blocks.

module DocsMathCompatibility

export check_docs_math_compat

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DOCS_SOURCE = joinpath(REPOSITORY_ROOT, "docs", "src")
const FORBIDDEN_MATH_COMMANDS = ("\\operatorname",)

function markdown_files(root::AbstractString)
    files = String[]
    for (directory, _, names) in walkdir(root)
        for name in names
            endswith(name, ".md") || continue
            push!(files, joinpath(directory, name))
        end
    end
    return sort!(files)
end

function unescaped_dollar_indices(text::AbstractString)
    indices = Int[]
    for index in eachindex(text)
        text[index] == '$' || continue
        backslashes = 0
        cursor = prevind(text, index)
        while cursor >= firstindex(text) && text[cursor] == '\\'
            backslashes += 1
            cursor == firstindex(text) && break
            cursor = prevind(text, cursor)
        end
        iseven(backslashes) && push!(indices, index)
    end
    return indices
end

function without_inline_code(text::AbstractString)
    return replace(text, r"`[^`]*`" => "")
end

function line_math_fragments(text::AbstractString, dollar_indices::Vector{Int})
    fragments = SubString{String}[]
    for pair_start in 1:2:length(dollar_indices)
        opening = dollar_indices[pair_start]
        closing = dollar_indices[pair_start + 1]
        first = nextind(text, opening)
        last = prevind(text, closing)
        first <= last && push!(fragments, SubString(text, first, last))
    end
    return fragments
end

function forbidden_command(fragment::AbstractString)
    return findfirst(command -> occursin(command, fragment), FORBIDDEN_MATH_COMMANDS)
end

function check_file(path::AbstractString)
    errors = String[]
    inline_count = 0
    display_count = 0
    fence_character = nothing
    fence_length = 0
    math_fence = false

    for (line_number, line) in enumerate(eachline(path))
        opening = match(r"^\s*(`{3,}|~{3,})([^`]*)$", line)
        if isnothing(fence_character)
            if !isnothing(opening)
                delimiter = opening.captures[1]
                fence_character = first(delimiter)
                fence_length = length(delimiter)
                info = strip(opening.captures[2])
                math_fence = info == "math"
                display_count += math_fence
                continue
            end
        else
            closing_pattern = Regex(
                "^\\s*" * string(fence_character) * "{" * string(fence_length) * ",}\\s*\$"
            )
            if occursin(closing_pattern, line)
                fence_character = nothing
                fence_length = 0
                math_fence = false
                continue
            end
            if math_fence
                command_index = forbidden_command(line)
                if !isnothing(command_index)
                    command = FORBIDDEN_MATH_COMMANDS[command_index]
                    push!(
                        errors,
                        "$(path):$(line_number): unsupported GitHub math command `$command`",
                    )
                end
            end
            continue
        end

        if occursin(r"(?<!`)``(?!`)", line)
            push!(
                errors,
                "$(path):$(line_number): use GitHub-compatible `\$...\$` inline math, not Documenter-only double backticks",
            )
        end

        prose = without_inline_code(line)
        dollar_indices = unescaped_dollar_indices(prose)
        if isodd(length(dollar_indices))
            push!(errors, "$(path):$(line_number): unbalanced inline-math dollar delimiter")
            continue
        end

        fragments = line_math_fragments(prose, dollar_indices)
        inline_count += length(fragments)
        for fragment in fragments
            command_index = forbidden_command(fragment)
            if !isnothing(command_index)
                command = FORBIDDEN_MATH_COMMANDS[command_index]
                push!(
                    errors,
                    "$(path):$(line_number): unsupported GitHub math command `$command`",
                )
            end
            if startswith(strip(prose), "|") && occursin('|', fragment)
                push!(
                    errors,
                    "$(path):$(line_number): use `\\vert` instead of a literal `|` inside table-cell math",
                )
            end
        end
    end

    if !isnothing(fence_character)
        push!(errors, "$(path): unclosed Markdown fence")
    end
    return (; errors, inline_count, display_count)
end

function check_docs_math_compat(; output::IO=stdout, error_output::IO=stderr)
    files = markdown_files(DOCS_SOURCE)
    all_errors = String[]
    inline_count = 0
    display_count = 0
    for path in files
        result = check_file(path)
        append!(all_errors, result.errors)
        inline_count += result.inline_count
        display_count += result.display_count
    end

    if !isempty(all_errors)
        foreach(error -> println(error_output, error), all_errors)
        return false
    end

    println(
        output,
        "Documentation math compatibility passed: ",
        length(files),
        " Markdown files, ",
        inline_count,
        " inline spans, ",
        display_count,
        " fenced display blocks.",
    )
    return true
end

end # module DocsMathCompatibility

if abspath(PROGRAM_FILE) == @__FILE__
    exit(DocsMathCompatibility.check_docs_math_compat() ? 0 : 1)
end
