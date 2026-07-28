# ADR 0001: provisional package name

- Status: Accepted, provisional
- Date: 2026-07-28
- Decision owners: project maintainers

## Context

The initial brief used `QETLAB.jl` as a default but required collision and
non-affiliation review. Reusing an upstream name can imply endorsement and makes
a later identity dispute expensive. Development must still continue before a
final public URL or registry submission exists.

Search evidence captured on 2026-07-28:

- GitHub's public repository search API returned zero exact-name repositories
  for `QuantumEntanglementTools.jl` and `QETLAB.jl`.
- A local Julia General registry snapshot with tree SHA
  `8bf2aa07751110d39d9833c3d4150b8655913ece` contained neither
  `Q/QuantumEntanglementTools/Package.toml` nor `Q/QETLAB/Package.toml`.
- GitHub General-registry content lookups for both names returned 404.
- JuliaHub returned a client-rendered shell, so an independent JuliaHub search
  result remains unverified.

Searches are point-in-time evidence, not a trademark clearance or a reservation.

## Decision

Use the neutral provisional repository and module name
`QuantumEntanglementTools` with package UUID
`45675e5b-5c8b-4983-b92d-4c3725d56c4e`.

Describe the project as an independent Julia implementation. Do not use upstream
logos or imply that QETLAB, QUBIT4MATLAB, their maintainers, or contributors
endorse it.

## Consequences

- The Julia-native API can evolve without claiming to be official QETLAB.
- A separate compatibility namespace can use documented upstream function names
  where legally and technically appropriate.
- Filenames, module references, docs configuration, and CI use one isolated
  package name so a pre-release rename remains mechanical.
- The name must be rechecked against the authoritative General registry,
  JuliaHub, GitHub, and relevant naming/trademark concerns immediately before
  publication. Configure a canonical repository URL only after ownership is
  known.
