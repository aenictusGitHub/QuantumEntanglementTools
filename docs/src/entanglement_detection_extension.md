# EntanglementDetection.jl extension

EntanglementDetection.jl 0.2.2 is an optional heuristic backend. It is not a
core dependency, and the adapter deliberately does not turn the backend's
`true`, `false`, or `nothing` result into a mathematical certificate.

## Loading and discovery

On Julia 1.11 or later, install EntanglementDetection.jl 0.2.2 in the active
environment, then load both packages in either order:

```julia
using Pkg
Pkg.add(name="EntanglementDetection", version="0.2.2")

using QuantumEntanglementTools
using EntanglementDetection

backends = available_entanglement_backends()
backend = only(filter(b -> b isa EntanglementDetectionBackend, backends))
backend_capabilities(backend)
```

The capability record reports `side_effect_free=false` for the upstream
backend, `caller_state_isolated=true`, `isolation=:child_process`,
`source_integrity_enforced=false`, and `resource_sandboxed=false`. Its
`transport_trust=:same_version_local_worker` value makes the IPC assumption
machine-readable. Discovery does not load EntanglementDetection.jl implicitly.
Before the weak dependency is loaded, only the native backend appears.

## Running a search

```julia
using LinearAlgebra

rho = Matrix{Float64}(I, 4, 4) / 4
method = EntanglementDetectionSearch(
    timeout_seconds=120,
    max_iteration=10_000,
    epsilon=1e-6,
)
report = detect_entanglement(rho, (2, 2), method)

@assert report.status === :unknown
@assert !report.certified
report.evidence.backend_conclusion
```

The final expression is one of `:entangled`, `:separable`, or
`:inconclusive`. It is candidate evidence only. In particular, a backend
`:separable` candidate is not a validated decomposition, and an
`:entangled` candidate is not a witness certificate independently verified by
this package. The report therefore remains `:unknown` with
`certified_by_adapter=false`.

Density matrices must be finite, normalized, positive semidefinite, and
compatible with the supplied dimensions. Sparse matrices are rejected unless
`allow_densify=true` is set explicitly. The adapter never normalizes or clips
the input. In the child only, a real matrix is represented as a complex matrix
with the same real component type. This explicit representation conversion
works around EntanglementDetection.jl 0.2.2's same-element-type dispatch
requirement; it does not change numerical precision or matrix entries. The
chosen representation is recorded in
`report.evidence.backend_input_representation`.

## Isolation and failures

EntanglementDetection.jl 0.2.2 seeds its default RNG during the public search.
Its optional logfile path can redirect global stdout, and its parallel LMO
option can change the process-wide BLAS thread count. The adapter does not try
to patch or reverse those behaviors. Every search starts a fresh Julia child
process in the current active project, fixes verbosity to zero, supplies no
logfile, captures child output, and discards the process after the result.

The child-process boundary adds startup and serialization overhead. It also
requires EntanglementDetection.jl to be resolvable from the active project.
Package-owned density-matrix validation and request serialization happen
before the child starts. `timeout_seconds` covers child startup and
computation; bounded termination cleanup can add a small amount of time after
that budget. A timeout, launch failure, nonzero child exit, malformed response,
or backend exception returns an uncertified `:unknown` report whose
`evidence.execution.error_kind` records the failure. Invalid package-owned
configuration or invalid density-matrix input still throws before a child is
started.

On timeout or an abnormal wait interruption, the adapter sends `SIGTERM`, waits
at most 0.1 seconds, then uses `SIGKILL` on Unix if needed. On Windows, libuv's
termination operation is used again. Forced-process reaping is bounded by two
seconds; an unreaped child is reported as `:termination_failed`, and its
temporary directory is left for Julia's process-exit cleanup instead of
starting synchronous cleanup against files that may still be open. User
interrupts are rethrown only after this bounded cleanup attempt.

Child stdout and stderr evidence is read only up to 8 KiB per stream. The files
are not size-limited while the child writes them; the audited call fixes
verbosity to zero, and a finite timeout bounds ordinary execution. Before
decoding, `response.bin` is limited to
`min(256 MiB, 1 MiB + 64 * length(rho))`. An oversized response is
`:invalid_response`.

The request and response use Julia's `Serialization` between the parent and a
worker launched from the same Julia executable and active project. This is
trusted, same-version local transport, not a security sandbox. The size cap and
post-decode schema/type validation contain ordinary corrupt or incompatible
worker output, but Julia serialization is not safe for hostile bytes, and a
maliciously altered payload may not be recoverable as a structured report.

The adapter enforces loaded package version 0.2.2. It does not hash the
installed source tree at runtime: the registered release or the pinned checkout
in `UpstreamManifest.toml` must be selected by the controlled environment. A
modified path dependency that retains version 0.2.2 is therefore trusted and
passes the version gate. The adapter exposes no in-process escape hatch and no
logfile or backend-parallelism option.

## Julia-version availability

The core package still supports Julia 1.10. The optional backend environment
currently resolves only on Julia 1.11 or later: EntanglementDetection.jl 0.2.2
requires Ket 0.9, and the compatible Ket 0.9 registry releases require Julia
1.11. This transitive constraint is stricter than EntanglementDetection.jl's
own declared Julia 1.9 compatibility. `backend_capabilities` records
`minimum_resolvable_julia=v"1.11.0"`, and the dedicated extension environment
therefore starts at Julia 1.11.

## Validation scope

The dedicated suite passes 125/125 assertions locally on Julia 1.12.6. It
covers both load orders, repeated loading, method ambiguities, configuration
and density validation, real-to-complex representation conversion, a live
heuristic search, preservation of caller RNG, stdout, logger, and BLAS state,
child RNG mutation, bounded timeout escalation, forced termination of a
TERM-resistant process on Unix, interrupt/error cleanup, response schema and
size checks, output truncation, and cleanup when output-stream setup fails
after launch. The core pipeline suite separately checks dependency absence and
the actionable error path. Independent smoke calls also covered `Float32`,
`Float64`, `ComplexF32`, `ComplexF64`, and multipartite input.

A Julia 1.11/1.12 Linux/macOS/Windows workflow is configured but has not run
remotely. The local evidence does not make backend candidates certificates or
make `Serialization` an adversarial isolation boundary.
