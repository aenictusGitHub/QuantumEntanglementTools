# QUBIT4MATLAB license and provenance

Status: **blocked by a material license conflict**

Evidence date: 2026-07-28

Target: QUBIT4MATLAB v6.5

## Decision

Do not inspect, translate, or adapt implementation `.m` files from the v6.5
archive until a human resolves the conflict between the license bundled in the
archive and the license separately supplied to the maintainer. No permission
claim for archive-derived work is made here.

## Archive evidence

| Evidence | Value |
|---|---|
| Canonical project page | <https://www.gtoth.eu/qubit4matlab.html> |
| Author-hosted archive | `QUBIT4MATLAB6.5.zip` |
| Retrieval date | 2026-07-28 |
| Archive SHA-256 | `282628dad2b8e134c7d88834770c1293abde3baa05664b161ce698ae30375547` |
| Bundled license path | `QUBIT4MATLAB/bsd.txt` |
| Bundled license SHA-256 | `b0bd73519dd6d16963c8aeb21ef9659e09dd3beb2536d9f736ec1d151f51ca81` |
| Bundled README SHA-256 | `0fecb2c1ea22f352efd59c586a3a09f22c96a5f5dfc8525b8bec80c017b788e2` |
| Implementation files inspected | No |
| Per-file audit | Not started; blocked |

The archive's bundled license identifies copyright years 2005–2015, has two
redistribution clauses, omits the supplied non-endorsement clause, and contains
a trailing quote. The separately supplied text identifies 2005–2024 and has
three clauses. The difference is material, not formatting-only.

## Required human action

Obtain authoritative clarification from the rightsholder or maintainer about
which terms apply to the v6.5 archive and whether redistribution of a
source-informed Julia translation is permitted. Preserve the response as review
evidence. Do not overwrite either existing license artifact. If clarification
changes the applicable terms, add a new verbatim file and record its provenance.

## Supplied text

The following text was supplied separately and is preserved verbatim in
`licenses/QUBIT4MATLAB-LICENSE.txt`. The blank after `Neither the name of` is an
anomaly in the supplied text and must not be filled in. This reproduction is not
a conclusion that the text governs the downloaded archive.

```text
Copyright (c) 2005-2024, Geza Toth
All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
 
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
 
* Neither the name of  nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Policy if the conflict is resolved

Before source inspection:

1. Record the clarification, archive checksum, applicable license path and
   checksum, retrieval date, and comparison result in `UpstreamManifest.toml`.
2. Create `feature/qubit4matlab-v6.5` from a tested main commit and record its
   base SHA.
3. Inventory every archive file without assuming that the project-level license
   overrides file-level notices.

For each adapted function, record the upstream filename, archive checksum,
implementation kind, Julia file, tests, and docs in `PROVENANCE.toml`. Add a
concise header to nontrivially adapted Julia files pointing to the applicable
full license, and keep the independent-project/non-endorsement wording in all
promotional material.

## File audit ledger

No implementation-file ledger exists because source inspection correctly stopped
at the conflicting archive license. When unblocked, replace this paragraph with
a reviewed per-file table covering authorship, notices, generated or third-party
code, overlap with the existing core, adaptation decision, and provenance entry.
