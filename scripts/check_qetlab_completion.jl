#!/usr/bin/env julia

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)
using QuantumEntanglementTools: QuantumEntanglementTools

include(joinpath(@__DIR__, "qetlab_completion_common.jl"))
using .QETLABCompletion

const RUNTIME_PUBLIC_MODULES = Dict(
    "QuantumEntanglementTools" => QuantumEntanglementTools,
    "QuantumEntanglementTools.MATLABCompat" => QuantumEntanglementTools.MATLABCompat,
)

function usage(io::IO=stdout)
    return print(
        io,
        """
Usage: julia scripts/check_qetlab_completion.jl [--check | --strict]

No argument validates the source ledgers and prints the current completion
report. --check additionally requires every generated completion artifact to
be current. --strict applies the implementation-completeness gate after the
same generated-artifact check; it is expected to fail while public rows remain
partial, deferred, or blocked. Every mode also checks runtime exports,
provenance/file evidence, compatibility delegation, source placeholders and
global-RNG candidates, and row-specific RNG/solver evidence.
""",
    )
end

function audited_report()
    model = QETLABCompletion.build_completion_model()
    audit = QETLABCompletion.repository_quality_audit(
        model; runtime_modules=RUNTIME_PUBLIC_MODULES
    )
    QETLABCompletion.completion_report(model; audit)
    return (; model, audit)
end

function main(args)
    if isempty(args)
        result = audited_report()
        isempty(result.audit.failures) || exit(1)
    elseif args == ["--check"]
        QETLABCompletion.write_or_check_outputs(; check=true) || exit(1)
        result = audited_report()
        isempty(result.audit.failures) || exit(1)
    elseif args == ["--strict"]
        QETLABCompletion.write_or_check_outputs(; check=true) || exit(1)
        result = audited_report()
        QETLABCompletion.strict_completion_check(result.model; audit=result.audit) ||
            exit(1)
    elseif args in (["-h"], ["--help"])
        usage()
    else
        usage(stderr)
        error("unknown arguments: $(join(args, ' '))")
    end
    return nothing
end

main(ARGS)
