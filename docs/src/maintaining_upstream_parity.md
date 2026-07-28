# Maintaining upstream parity

Upstream synchronization is a reviewed comparison, never an automatic overwrite
of Julia implementations.

## Update procedure

1. Fetch upstream into an ignored development checkout.
2. Record the exact new revision, retrieval date, and license hashes.
3. Run `scripts/check_upstream_updates.jl --candidate qetlab=PATH` against the
   candidate checkout. The default `PROVENANCE.toml` supplies mapped Julia
   functions and tests in the reverse-impact report; use
   `--provenance OTHER_PATH` only to select another ledger.
4. Review added, removed, and changed MATLAB files, signatures, help text, and
   license files.
5. Map affected Julia functions through `PROVENANCE.toml`.
6. Update inventory classifications deliberately.
7. Rerun affected analytic, property, differential/independent, sparse/generic,
   documentation, and benchmark checks.

Do not replace a pin with a moving branch name and do not describe a checkout as
reproducible without its commit or archive checksum.

The checker is implemented and read-only. It compares `.m` membership and
SHA-256 hashes, detects signature/help/license changes, reports reverse
dependency impact from the generated graph, and can print provenance mappings.
It does not fetch, merge, or rewrite either checkout. Its offline self-check
against the pinned checkout and a historical QETLAB tag passed locally on
2026-07-28; this is tooling evidence, not a claim that the inventory's manual
classification is complete.

The companion consistency commands are:

```sh
julia --project=. scripts/build_upstream_inventory.jl --check
julia --project=. scripts/check_public_api.jl
```

The second command is intentionally strict: an export without provenance, or an
implemented inventory row without native and compatibility mappings, fails the
check.
