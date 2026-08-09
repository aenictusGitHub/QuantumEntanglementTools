# Legal and provenance status

Last evidence update: 2026-07-31.

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

## Browser code generator

The interaction pattern of the documentation's browser code generator was
inspired by the public
[PermutationalInvariantDynamics.jl model code generator](https://aenictusgithub.github.io/PermutationalInvariantDynamics.jl/dev/model_code_generator/),
inspected in a clean disposable checkout on 2026-07-31 at commit
[`49c64b1c0fc5b301531582d470144c5b6b3d4030`](https://github.com/aenictusGitHub/PermutationalInvariantDynamics.jl/tree/49c64b1c0fc5b301531582d470144c5b6b3d4030)
(tree `39407ccd1d38f1332b987b3365f5d8ef0b317600`). The inspected paths were
`docs/src/model_code_generator.md` and its three
`docs/src/assets/model_code_generator_{core.js,ui.js,css}` assets. Exact path
hashes are recorded in `UpstreamManifest.toml`.

The reference repository's `REUSE.toml` applies GPL-3.0-only to the complete
tree. Its root `LICENSE` and `LICENSES/GPL-3.0-only.txt` both have SHA-256
`8ceb4b9ee5adedde47b31e975c1d90c73ad27b6b165a1dcd80c7c545eb65b903`.
The reference generator documentation and emitted headers explicitly classify
its generated Julia templates as GPL-3.0-only without a generated-output
exception.

No reference HTML, CSS, JavaScript, prose, identifiers, or Julia templates were
copied or adapted. The QuantumEntanglementTools generator is an independently
written implementation based only on this package's public API contracts and
the general idea of a deterministic browser form with copy/download controls.
Its original source assets and generated template text are covered by this
repository's BSD 3-Clause License. Emitted Julia files carry SPDX copyright and
license identifiers plus a link to the complete project license. Redistribution
of substantial generated template text must retain that notice and the
applicable BSD conditions. Configuration JSON contains normalized user choices,
not the Julia template; users remain responsible for rights in data or other
material they supply.

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

The integration is an exact-version optional Julia package dependency, not
vendored source. Version 0.2.2 has caller-visible global-state hazards:
`separable_distance` calls `Random.seed!(0)`, may redirect global stdout when a
log file is supplied, and performs logging/flush operations;
`AlternatingSeparableLMO(..., parallelism=true)` globally sets BLAS threads to
one.

The package-owned extension calls the documented public backend API only in a
fresh Julia child process. It supplies no logfile or backend-parallelism option,
contains failures as uncertified `unknown` reports, and never promotes the
backend's heuristic Boolean to a certificate. The local 130-assertion extension
suite covers caller RNG/stdout/logger/BLAS preservation, load order, lifecycle,
timeout escalation, and structured response handling on Julia 1.12.6. This is
engineering evidence for the stated isolation boundary, not a legal conclusion,
a security sandbox, or remote supported-platform validation. The optional
dependency currently resolves only on Julia 1.11 or later because registered
Ket 0.9 releases have that compatibility floor. Runtime compatibility checks
enforce version 0.2.2 but do not authenticate a source tree; controlled
validation and release environments must use the registered release or the
recorded pinned checkout.

## Published SAPPT example

The executable tutorial `tutorials/symmetric_sappt_witnesses.jl` cites and
reproduces selected calculations from:

- J. Louvet, E. Serrano-Ensástiga, T. Bastin, and J. Martin,
  “Nonequivalence between absolute separability and positive partial
  transposition in the symmetric subspace,” *Physical Review A* **111**,
  042418 (2025), <https://doi.org/10.1103/PhysRevA.111.042418>.

The inspected author-supplied publisher PDF has SHA-256
`c9a32248a9f0f730409dc690db0dfc8ab97063bd600041cb7e567d8b0e9aaf7a`
and states `Copyright 2025 American Physical Society`. The PDF, paper prose,
and figures are not vendored or redistributed. The tutorial is independently
written Julia code under this repository's BSD 3-Clause License. It cites the
analytic SAPPT threshold and reconstructs rounded witness coefficients printed
in the paper; it does not copy or rerun the authors' source SDP. The full source
relationship and artifact evidence are recorded in `UpstreamManifest.toml`.

## Distribution checklist

Before any release or source archive:

- verify that `LICENSE`, `NOTICE`, and `THIRD_PARTY_LICENSES.md` ship;
- ensure `UpstreamManifest.toml` pins every inspected source and license hash;
- ensure every public/source-derived function is in `PROVENANCE.toml`;
- reproduce all applicable upstream notices in source and binary documentation;
- check every vendored fixture or asset independently;
- ensure the README and package metadata do not imply upstream endorsement;
- resolve or exclude every `blocked` or `unknown` license classification.
