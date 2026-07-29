# Third-party licenses

This document is a compliance ledger, not a statement that every project listed
below is a runtime dependency or that its source is redistributed here. Exact
file/function provenance belongs in `PROVENANCE.toml`.

## QETLAB

- Project: QETLAB
- Maintainer: Nathaniel Johnston
- Canonical source: <https://github.com/nathanieljohnston/QETLAB>
- Inspected commit: `d8589610f00cff106537268dee2e2a1153f3a601`
- License: BSD 2-Clause
- Root license SHA-256:
  `f861a171af77f8377e973b22ed08361b0069e90b846e4b4f49be58217a08792f`
- Verbatim copy: [`licenses/QETLAB-LICENSE.txt`](licenses/QETLAB-LICENSE.txt)

Do not copy MATLAB source, documentation assets, solver components, or helper
routines merely because they are present upstream. Preserve the exact upstream
notice for every adapted portion and audit separately licensed helpers.

Full QETLAB terms:

```text
Copyright (c) 2014, Nathaniel Johnston
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

### QETLAB helpers by Bruno Luong

The inspected QETLAB tree contains separate, byte-identical BSD 2-Clause notices
for `helpers/spnull.m` and `helpers/sporth.m`:

- notice paths: `helpers/license_spnull.txt`,
  `helpers/license_sporth.txt`
- SHA-256:
  `ff2b568397db2682226129f05a9425c07afca47e0d22d302fc5577fa49373993`
- verbatim copy:
  [`licenses/QETLAB-BRUNO-LUONG-HELPERS-LICENSE.txt`](licenses/QETLAB-BRUNO-LUONG-HELPERS-LICENSE.txt)

Those helpers are not claimed as ported. If their implementations inform future
work, retain these complete terms:

```text
Copyright (c) 2010, Bruno Luong
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

## EntanglementDetection.jl

- Project: EntanglementDetection.jl
- Canonical source:
  <https://github.com/ZIB-IOL/EntanglementDetection.jl>
- Inspected version and commit: 0.2.2,
  `5f60da1ceef6442acb669e10acc2fa47670bab06`
- License: MIT; inspected SHA-256:
  `8d7174972190a50508c86dcc936c31d23cccea7ab4b381aa49ab8735ae326824`
- Intended use: exact-version optional dependency through a Julia extension and
  a package-owned child-process adapter; no upstream implementation is
  vendored.

The audited version mutates its process-local RNG and has stdout/logging and
BLAS-thread options with global effects. The adapter therefore invokes its
public API only in a fresh child process and treats every backend conclusion as
uncertified candidate evidence. Its package-manager installation supplies the
applicable MIT text. The exact version and isolation assumptions must be
reassessed before widening compatibility.

## System and development tools

Julia, Git, Graphviz, GNU Octave, GLPK, documentation tools, CI actions, and
other developer-installed programs are not redistributed merely because they
are used to build or test this repository. Their own licenses apply in their
respective distributions.
