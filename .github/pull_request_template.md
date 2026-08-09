## Purpose

Describe the user-visible problem and the scoped change.

## Evidence

List the focused commands run, Julia versions, platforms, fixtures, and relevant numerical or analytic cross-checks. A passing source-free fixture is supplemental evidence, not a general MATLAB/QETLAB parity claim.

## Review checklist

- [ ] Public names and documented conventions remain compatible, or the breaking change is explicitly justified for a pre-1.0 release.
- [ ] Shapes, subsystem order, vectorization, normalization, tolerances, failure behavior, and complexity/resource limits are documented.
- [ ] Sparse inputs remain sparse unless guarded by an explicit, tested `allow_densify=true` path.
- [ ] Randomized APIs accept an explicit `rng::AbstractRNG` and do not mutate the caller's global random stream.
- [ ] Necessary tests, heuristics, numerical evidence, and certificates remain distinguishable; inconclusive results remain `unknown`.
- [ ] Tests cover normal, boundary, invalid, generic-numeric, and sparse/resource-limit behavior as applicable.
- [ ] Every public/source-informed change has synchronized specification, provenance, source attribution, tests, and documentation.
- [ ] Optional solvers/integrations remain in extensions and use public dependency APIs.
- [ ] Generated ledgers and maintained claim counts were regenerated and checked.
- [ ] No credentials, private data, proprietary solver files, ignored manifests, or unrelated workspace changes are included.

## Maintainer review

For release-bound changes, identify the human reviewer and durable review record. Do not check this box on an AI agent's behalf.

- [ ] Human mathematical/API/provenance review completed and recorded.
