#!/usr/bin/env julia

include(joinpath(@__DIR__, "qetlab_completion_common.jl"))
using .QETLABCompletion

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/build_qetlab_completion_plan.jl [--check]

Generate the all-public-row QETLAB completion ledger, the incomplete-row queue,
and the human-readable completion plan. With --check, fail when an artifact is
missing or stale instead of writing it.
""",
    )
end

function main(args)
    if args == ["--check"]
        QETLABCompletion.write_or_check_outputs(; check=true) || exit(1)
    elseif isempty(args)
        QETLABCompletion.write_or_check_outputs(; check=false) || exit(1)
    elseif args in (["-h"], ["--help"])
        usage()
    else
        usage(stderr)
        error("unknown arguments: $(join(args, ' '))")
    end
    return nothing
end

main(ARGS)
