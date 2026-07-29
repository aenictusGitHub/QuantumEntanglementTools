# Release checklist

This checklist prepares the experimental `v0.1.0` scope. An unchecked item is
not release evidence. Record the exact commit, command output, reviewer, date,
and any accepted limitation in `docs/SESSION_HANDOFF.md` before tagging.

The two release targets are intentionally separate:

1. a private GitHub release, visible only to repository collaborators; and
2. registration in Julia's public General registry.

A private release does not imply General-registry readiness. Do not change
repository visibility, rename the repository, create a tag, publish a release,
or invoke Registrator without the maintainer's explicit approval.

## Non-delegable human review

These items must be completed personally by the maintainer. An AI agent may
prepare evidence but may not check them on the maintainer's behalf.

- [ ] Read and understand every substantive Codex-assisted source and extension
      file included in the release.
- [ ] Review the mathematical contracts, failure modes, tolerance boundaries,
      sparse-densification gates, explicit-RNG behavior, and certificate versus
      heuristic semantics.
- [ ] Review the tests and executable examples critically rather than treating a
      passing generated test as proof that its expected result is correct.
- [ ] Confirm each QETLAB-informed implementation against its specification,
      provenance entry, source header, and applicable BSD-2 notice.
- [ ] Record the reviewed commit, reviewer name, review date, discovered issues,
      and their resolution in a durable release issue or signed review record.
- [ ] Confirm that the OpenAI Codex disclosure in `README.md` accurately
      describes the assistance received.

General states that an LLM-assisted package is suitable only when its human
maintainer understands the generated code; unreviewed “vibe-coded” packages are
not suitable for registration.

## Scope and metadata freeze

- [ ] Confirm that `Project.toml` contains the intended name, UUID, version
      `0.1.0`, authors, Julia floor, weak dependency, extension, and compat
      bounds.
- [ ] State that `v0.1.0` covers only the API recorded in `PROVENANCE.toml`.
- [ ] Keep all 64 partial, deferred, or explicitly blocked public QETLAB rows
      visible; zero pending classifications does not imply complete QETLAB
      parity.
- [ ] Run `julia --startup-file=no --project=. scripts/check_public_api.jl` and
      verify every exported binding has specification, provenance, tests, and
      documentation.
- [ ] Confirm that `CHANGELOG.md`, `CITATION.cff`, `CITATION.bib`,
      `SECURITY.md`, `README.md`, and `docs/PORTING_STATUS.md` describe the same
      version and release scope.
- [ ] Immediately before the approved tag, confirm that the dated changelog
      heading and `CITATION.cff` `date-released` equal the actual publication
      date. If publication is not 2026-07-29, update both on a new candidate
      commit and rerun every gate; the preparation date is not a release date.
- [ ] While iterating, run the explicitly non-evidentiary worktree preflight
      after all intended release files are tracked:

  ```sh
  julia --startup-file=no --project=. scripts/check_release.jl --allow-dirty
  ```

- [ ] Validate `CITATION.cff` against CFF schema 1.2.0 and record the validator
      version:

  ```sh
  cffconvert --version
  cffconvert --validate --infile CITATION.cff
  ```

- [ ] Confirm that all canonical repository and documentation links match the
      chosen release target.

## Local release gates

Run each command from the repository root on the exact candidate commit:

- [ ] Core package:

  ```sh
  julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
  julia +1.10 --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
  ```

- [ ] Executable tutorials:

  ```sh
  julia --startup-file=no --project=. tutorials/runtests.jl
  ```

- [ ] Strict documentation and browser-generator checks:

  ```sh
  julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
  julia --startup-file=no --project=docs docs/make.jl
  generator_smoke="$(mktemp -d)"
  node docs/test/code_generator_node.js "$generator_smoke"
  julia --startup-file=no --project=docs "$generator_smoke/generated_smoke.jl"
  julia +1.10 --startup-file=no --project=docs "$generator_smoke/generated_smoke.jl"
  ```

  If Node is unavailable, replace its command with the exact macOS
  JavaScriptCore fallback below and record which engine was used:

  ```sh
  /usr/bin/osascript -l JavaScript docs/test/code_generator_jxa.js "$generator_smoke"
  ```

- [ ] Quality, formatting, inventory, and independent predicate validation:

  ```sh
  julia --startup-file=no --project=quality quality/run_quality.jl
  julia --startup-file=no --project=quality quality/format.jl
  julia --startup-file=no --project=. scripts/build_upstream_inventory.jl --check
  julia --startup-file=no --project=. scripts/check_public_api.jl
  julia --startup-file=no --project=. scripts/validate_matrix_predicates.jl
  ```

- [ ] Optional EntanglementDetection.jl environment on Julia 1.11 or later:

  ```sh
  julia --project=test/extensions/entanglement_detection -e '
      using Pkg
      Pkg.develop(PackageSpec(path=pwd()))
      Pkg.instantiate()
  '
  julia --startup-file=no --project=test/extensions/entanglement_detection \
    test/extensions/entanglement_detection/runtests.jl
  ```

- [ ] Quick benchmark smoke, without presenting it as a comparative performance
      result:

  ```sh
  julia --startup-file=no --project=benchmark \
    benchmark/benchmarks.jl --quick --no-save
  ```

- [ ] Run every committed source-free oracle comparator. Each comparator
      verifies its fixture's committed `.sha256` sidecar before comparing
      native and compatibility results:

  ```sh
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_a_oracle.jl test/oracle/fixtures/tier_a_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_b_oracle.jl test/oracle/fixtures/tier_b_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_c_oracle.jl test/oracle/fixtures/tier_c_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_d_oracle.jl test/oracle/fixtures/tier_d_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_e_coherence_oracle.jl test/oracle/fixtures/tier_e_coherence_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_e_product_oracle.jl test/oracle/fixtures/tier_e_product_octave_11_3_qetlab_d858961.json
  julia --startup-file=no --project=test/oracle test/oracle/compare_tier_e_matrix_analysis_oracle.jl test/oracle/fixtures/tier_e_matrix_analysis_octave_11_3_qetlab_d858961.json
  ```

  Record that Octave evidence is function-specific and does not establish
  general MATLAB equivalence.

- [ ] On a clean exact candidate commit, run the release archive checks on both
      the current Julia and the minimum supported Julia:

  ```sh
  julia --startup-file=no --project=. scripts/check_release.jl --archive-smoke
  julia +1.10 --startup-file=no --project=. scripts/check_release.jl --archive-smoke
  ```

## Remote evidence

- [ ] Core CI passes on every declared Julia version and supported operating
      system.
- [ ] Documentation CI and the rendered-site artifact pass.
- [ ] Quality, formatting, inventory, and API/provenance workflows pass.
- [ ] Coverage tests pass and Codecov confirms successful ingestion of the
      report; a green job that skipped or failed upload is not sufficient.
- [ ] The EntanglementDetection Julia 1.11/1.12 Linux, macOS, and Windows
      matrix passes.
- [ ] Nightly failures are either fixed or recorded as upstream/nightly-only;
      they are not silently ignored.
- [ ] Required status checks and branch protection refer to the current workflow
      job names.
- [ ] A GitHub ruleset protects release tags matching `v*` from update or
      deletion after publication.

## Legal and archive inspection

- [ ] Re-run the distribution checklist in `docs/LEGAL.md`.
- [ ] Verify the hashes and classifications in `UpstreamManifest.toml`.
- [ ] Confirm that no `blocked` or `unknown` license classification is included
      in the release.
- [ ] Confirm that `LICENSE`, `NOTICE`, `THIRD_PARTY_LICENSES.md`,
      `PROVENANCE.toml`, `UpstreamManifest.toml`, and both files under
      `licenses/` are tracked and present in the archive.
- [ ] Review every committed fixture, paper-derived number, documentation asset,
      and generated-template notice for redistribution rights and attribution.
- [ ] Search the complete history and candidate tree for credentials, private
      data, proprietary solver files, unapproved upstream assets, and accidental
      development checkouts.
- [ ] Resolve the reachable-history finding recorded in
      `docs/SESSION_HANDOFF.md`: either obtain explicit approval for a
      coordinated history rewrite/clean migration, or record the maintainer's
      reviewed decision that removal applies to the release tree only. Do not
      make the repository public while this decision is unresolved.
- [ ] Confirm that ignored manifests, generated docs, benchmark output, and local
      workspace files are absent from the committed archive.
- [ ] Inspect and smoke-test the exact candidate commit before tagging:

  ```sh
  candidate_commit="$(git rev-parse HEAD)"
  julia --startup-file=no --project=. scripts/check_release.jl \
    --treeish "$candidate_commit" --archive-smoke
  git archive --format=tar "$candidate_commit" | tar -tf -
  ```

- [ ] After explicit approval creates the tag, verify that the tag targets the
      tested commit, re-run the archive smoke test from that exact tag, and
      publish the source archive's SHA-256 checksum:

  ```sh
  git rev-parse 'v0.1.0^{commit}'
  julia --startup-file=no --project=. scripts/check_release.jl \
    --treeish v0.1.0 --tag v0.1.0 --archive-smoke
  release_archive="/tmp/QuantumEntanglementTools-v0.1.0.tar"
  git archive --format=tar v0.1.0 > "$release_archive"
  shasum -a 256 "$release_archive"
  ```

## Private GitHub release

Use this section only if the repository is to remain private.

- [ ] Obtain explicit approval to publish a private GitHub release.
- [ ] Confirm that collaborators understand the release and its documentation
      are not publicly accessible.
- [ ] Confirm the release commit is pushed, branch protection is active, and all
      required checks pass on that exact commit.
- [ ] Decide whether GitHub should mark `v0.1.0` as a pre-release; do not call
      the experimental API stable.
- [ ] Create an annotated `v0.1.0` tag on the exact tested commit.
- [ ] Push that tag only after explicit approval, then wait for every
      tag-triggered Core, Documentation, Coverage, Quality, and Optional
      integration workflow to pass on the tag's exact commit before publishing
      the GitHub release.
- [ ] Create release notes from the `0.1.0` changelog section, including the
      no-parity statement, Julia requirements, optional-backend floor, and known
      limitations.
- [ ] Attach only inspected artifacts and checksums.
- [ ] In a fresh depot with authenticated private-repository access, install the
      exact tag and run a package-load plus small certificate smoke test.
- [ ] Verify that the tag, GitHub release, changelog, citation metadata, and
      source archive identify the same tree.

If a later General registration reuses `v0.1.0`, its registered tree must match
the existing tag exactly and the TagBot strategy must avoid conflicting tags.

## Public General-registry registration

Complete this section in addition to every applicable section above.

- [ ] Obtain explicit approval before changing the repository from private to
      public.
- [ ] Before making it public, complete the history-wide secret, privacy,
      proprietary-material, and license review.
- [ ] Obtain explicit approval before renaming the repository to
      `QuantumEntanglementTools.jl`. General's automatic-merge policy expects
      the repository URL to end in `/QuantumEntanglementTools.jl.git`.
- [ ] After any rename or visibility change, update remotes, documentation,
      citation metadata, changelog links, badges, branch rules, and deployment
      settings; re-run all remote checks.
- [ ] Run General's current AutoMerge.jl exact-name and similarity checks
      against an up-to-date registry snapshot.
- [ ] After personally completing and recording the non-delegable review, run
      the local partial registry preflight. The environment variable is only a
      maintainer attestation; it is not independent proof of review or
      RegistryCI acceptance:

  ```sh
  QET_HUMAN_REVIEW_CONFIRMED=true julia --startup-file=no --project=. \
    scripts/check_release.jl --registry
  ```

- [ ] Confirm the package name and non-affiliation wording with the maintainer;
      registration makes the name and UUID effectively permanent.
- [ ] Confirm every root dependency and weak dependency version is registered
      and every compat entry satisfies current RegistryCI upper-bound rules.
- [ ] Confirm the public repository exposes an OSI-approved top-level license
      and all required third-party notices.
- [ ] Keep the detailed Codex-assistance disclosure in the public README and
      link the maintainer's completed human-review record.
- [ ] Keep the README concise and user-oriented; move detailed validation data
      to the evidence documents.
- [ ] Deploy public, versioned documentation with working repository and source
      links.
- [ ] Configure Registrator and choose either TagBot or a documented manual
      tagging process. Configure CompatHelper after the canonical URL is final.
- [ ] Trigger Registrator only from the exact approved release commit and
      include scoped release notes.
- [ ] Review the General pull request, respond personally to maintainer feedback,
      and observe the new-package waiting period.
- [ ] After merge, use a fresh depot to run:

  ```sh
  release_depot="$(mktemp -d)"
  JULIA_DEPOT_PATH="$release_depot" julia --startup-file=no -e '
      using Pkg
      Pkg.add("QuantumEntanglementTools")
      using QuantumEntanglementTools
      @assert bell_state() isa AbstractVector
  '
  ```

## Post-release

- [ ] Verify the published tag and source archive resolve to the tested commit.
- [ ] Verify installation and documentation links from a machine or account
      without maintainer-local state.
- [ ] Verify GitHub displays the intended license, citation, security policy, and
      release notes.
- [ ] Publish versioned documentation and retain prior versions.
- [ ] Leave an empty `Unreleased` changelog section for subsequent changes.
- [ ] Update `docs/PORTING_STATUS.md` and `docs/SESSION_HANDOFF.md` with commands
      actually run, remote URLs, failures, remaining external gates, and the
      final dirty state.
- [ ] Announce only the verified scope; do not imply full QETLAB parity,
      general MATLAB equivalence, solver certification, or comparative speed.

## Authoritative external guidance

- [Julia General registry README](https://github.com/JuliaRegistries/General)
- [RegistryCI automatic-merge guidelines](https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/)
- [Registrator instructions and private-package policy](https://github.com/JuliaRegistries/Registrator.jl)
- [Pkg project and manifest documentation](https://pkgdocs.julialang.org/v1/toml-files/)
- [Pkg compatibility documentation](https://pkgdocs.julialang.org/v1/compatibility/)
- [TagBot documentation](https://github.com/JuliaRegistries/TagBot)
- [Citation File Format 1.2.0 guide](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md)
