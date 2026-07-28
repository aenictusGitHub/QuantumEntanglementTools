# Legal and provenance status

Last evidence update: 2026-07-28.

This file documents an engineering compliance review, not legal advice. A source
being publicly accessible does not by itself authorize copying, adaptation, or
redistribution. Function-level implementation origin must also be recorded in
`PROVENANCE.toml`.

## New project work

New, independently written project code and documentation are licensed under
the BSD 3-Clause License in the repository root. The provisional copyright
notice is:

```text
Copyright (c) 2026, John Martin
```

Contributors must have the right to submit their changes. Files adapted from
other projects may carry additional compatible notices that the root license
does not replace.

## QETLAB

The inspected upstream source is:

- repository: <https://github.com/nathanieljohnston/QETLAB>
- commit: `d8589610f00cff106537268dee2e2a1153f3a601`
- retrieval/inspection date: 2026-07-28
- root `license.txt` SHA-256:
  `f861a171af77f8377e973b22ed08361b0069e90b846e4b4f49be58217a08792f`
- root license: BSD 2-Clause, `Copyright 2014 Nathaniel Johnston`

Two helper notices were also identified:

- `helpers/license_spnull.txt`
- `helpers/license_sporth.txt`
- identical SHA-256:
  `ff2b568397db2682226129f05a9425c07afca47e0d22d302fc5577fa49373993`
- BSD 2-Clause, `Copyright 2010 Bruno Luong`

Those helper notices must be preserved if `spnull`, `sporth`, or source derived
from them is redistributed. Their presence does not mean those helpers have been
ported. Every direct or source-informed adaptation must name its source file,
revision, implementation kind, tests, and docs in `PROVENANCE.toml`. Upstream
MATLAB code, CVX, solver binaries, and documentation assets must not be vendored
without a specific redistribution review.

Verbatim retained copies are:

- `licenses/QETLAB-LICENSE.txt`;
- `licenses/QETLAB-BRUNO-LUONG-HELPERS-LICENSE.txt`.

For a Julia file whose implementation is informed nontrivially by QETLAB source,
use a concise header like the following, with the exact filename filled in:

```julia
# Source-informed independent Julia implementation based on the specification
# and QETLAB <UpstreamFile.m> at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
```

The header is not a substitute for a `PROVENANCE.toml` entry. If code is a direct
translation or adaptation rather than an independently structured
source-informed implementation, label that implementation kind accurately. For
`spnull`/`sporth`-informed work, name Bruno Luong and point to the separate helper
license instead.

## EntanglementDetection.jl

The inspected reference is:

- repository: <https://github.com/ZIB-IOL/EntanglementDetection.jl>
- version: 0.2.2
- commit/tag SHA: `5f60da1ceef6442acb669e10acc2fa47670bab06`
- MIT license SHA-256:
  `8d7174972190a50508c86dcc936c31d23cccea7ab4b381aa49ab8735ae326824`

The intended integration is an optional Julia package dependency, not vendored
source. Version 0.2.2 has caller-visible global-state hazards:
`separable_distance` calls `Random.seed!(0)`, may redirect global stdout when a
log file is supplied, and performs logging/flush operations;
`AlternatingSeparableLMO(..., parallelism=true)` globally sets BLAS threads to
one. The extension is not safe or complete until those effects are isolated and
covered by regression tests.

## QUBIT4MATLAB v6.5: blocking conflict

The maintainer supplied a BSD 3-Clause–style text with
`Copyright (c) 2005-2024, Geza Toth`. It is preserved byte-for-byte, including
single-space blank lines and the malformed phrase `Neither the name of  nor`, at
[`licenses/QUBIT4MATLAB-LICENSE.txt`](../licenses/QUBIT4MATLAB-LICENSE.txt).
The malformed clause is conservatively understood as a non-endorsement
restriction covering the project and contributor names; it has not been
silently repaired.

Audit evidence for the author-hosted v6.5 archive:

- retrieval date: 2026-07-28
- archive SHA-256:
  `282628dad2b8e134c7d88834770c1293abde3baa05664b161ce698ae30375547`
- bundled license path: `QUBIT4MATLAB/bsd.txt`
- bundled license SHA-256:
  `b0bd73519dd6d16963c8aeb21ef9659e09dd3beb2536d9f736ec1d151f51ca81`
- bundled README SHA-256:
  `0fecb2c1ea22f352efd59c586a3a09f22c96a5f5dfc8525b8bec80c017b788e2`

The bundled license materially conflicts with the supplied text: it identifies
copyright years 2005–2015, contains only two redistribution clauses, omits the
non-endorsement clause present in the supplied text, and contains a trailing
quote. It is neither byte-identical nor materially identical.

Consequently:

1. The preserved supplied text is evidence only; this project does **not** claim
   that it governs the v6.5 archive.
2. No `.m` implementation files from the archive were inspected during the
   audit.
3. Direct or source-informed QUBIT4MATLAB porting is blocked pending human
   clarification of the conflicting terms.
4. If clarified, every archive file still requires an individual authorship,
   generated-code, and third-party-license audit before adaptation.
5. QUBIT4MATLAB work must remain on an isolated branch and requires human review
   before merge.

See `docs/QUBIT4MATLAB_LICENSE_AND_PROVENANCE.md` for the evidence ledger and
future review checklist.

## Distribution checklist

Before any release or source archive:

- verify that `LICENSE`, `NOTICE`, and `THIRD_PARTY_LICENSES.md` ship;
- ensure `UpstreamManifest.toml` pins every inspected source and license hash;
- ensure every public/source-derived function is in `PROVENANCE.toml`;
- reproduce all applicable upstream notices in source and binary documentation;
- check every vendored fixture or asset independently;
- ensure the README and package metadata do not imply upstream endorsement;
- resolve or exclude every `blocked` or `unknown` license classification.
