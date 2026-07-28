# QUBIT4MATLAB v6.5 port report

Status: **not started; blocked before implementation-source inspection**

Report date: 2026-07-28

This is the required future branch report skeleton. It is not evidence that a
branch, mapping, port, test, or benchmark exists.

## Branch and source

| Field | Current value |
|---|---|
| Feature branch | Not created |
| Base main SHA | Not applicable; repository had no base commit at scaffold time |
| Archive version | Advertised v6.5 |
| Archive SHA-256 | `282628dad2b8e134c7d88834770c1293abde3baa05664b161ce698ae30375547` |
| Archive license | Materially conflicts with separately supplied license |
| Human clarification | Required; not available |
| `.m` implementation files inspected | No |

See `QUBIT4MATLAB_LICENSE_AND_PROVENANCE.md` for hashes, the exact discrepancy,
and the separately supplied text.

## Function mapping

Not generated. After legal clearance, every archive function must be classified
as exact reuse of the verified native core, a thin compatibility wrapper, an
extension, a genuinely new feature, deprecated/superseded, or explicitly
blocked. Do not maintain duplicate production partial-trace or
partial-transpose algorithms solely to mirror two MATLAB projects.

## Validation

No QUBIT4MATLAB-derived validation fixtures or results exist.

## Benchmarks

No QUBIT4MATLAB-derived benchmark results exist.

## API additions and compatibility risk

Not assessed because the source inventory is blocked.

## Human review checklist

- [ ] Authoritative clarification resolves which terms govern v6.5.
- [ ] Clarification and applicable verbatim terms are preserved as evidence.
- [ ] Isolated branch is created from a tested main SHA.
- [ ] Archive contents receive a per-file authorship/license audit.
- [ ] Inventory and overlap mapping are complete.
- [ ] Every derived function/file has provenance and applicable attribution.
- [ ] Analytic, property, independent/differential, error, sparse/generic, and
      RNG/global-state tests pass as applicable.
- [ ] Documentation and benchmarks support every claim.
- [ ] Notices ship in source and binary/package materials.
- [ ] Wording and assets do not imply endorsement.
- [ ] Human technical and legal review explicitly approves merge.
