# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/) together with
Julia's pre-`1.0` compatibility convention: breaking public-API changes require
a minor-version increment.

## [Unreleased]

The changes below are part of the unreleased `0.1.0` development milestone
toward behavioral coverage of the public API at the pinned QETLAB revision. No
version has been tagged or published.

<!-- qetlab-current-claims: begin -->

The inventory is pinned to QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`. Its strict static ledger reports
127/127 public rows are verified with the required final status, 36/36 internal
helpers have terminal dispositions, the completion queue contains 0 public
rows, 0 required internal helpers remain, and 0 static completion failures.
The package has 477 public bindings (347 native/module and 130
`MATLABCompat`) with matching provenance entries.

The 9,484-assertion full package suite passed 9,484/9,484: the core accounts
for 9,400 assertions, and the seven executable tutorials account for 84/84
assertions (48 existing plus 36 for the two new workflows), on Julia
1.12.6 and the installed Julia 1.10.11. The full optional JuMP suite passed
847/847 on both Julia lines. The
EntanglementDetection.jl extension passed 141/141 focused assertions on the
current compatible Julia. All 114 declared quick benchmark cases completed
without failure on the current uncommitted release-hardening worktree based on
`8b2fcbaf`. This is local smoke evidence only; targeted paired observations
remain diagnostics, not a stable comparative-performance baseline.

These static and local results are not QETLAB/MATLAB parity,
supported-platform remote CI, comparative performance, API stability, release
approval, or non-delegable human review. No release, tag, or publication has
been made.

<!-- qetlab-current-claims: end -->

### Added

- A project-native symmetric-state toolkit for qubits and qudits: exact
  occupation dimensions and rank/unrank, generalized Dicke states, stable
  symmetric product coordinates, compressed collective one-body operators,
  sparse bipartition isometries, direct reduced states, and exact
  coordinate/ambient maximally mixed symmetric states. All combinatorial
  allocation paths have explicit `BigInt` resource preflights.
- Initial Julia package shell.
- Documentation, governance, legal, citation, and CI scaffolding.
- Explicit experimental status and milestone ledgers.
- Tier A subsystem/indexing, tensor, permutation, partial trace/transpose,
  realignment, projector/basis, reusable-plan, and MATLAB-compatibility APIs.
- Tier B Pauli, generalized Pauli, Gell-Mann, generalized Gell-Mann, Fourier,
  named-state, and explicit-RNG random-object constructors with compatibility
  wrappers.
- Tier B analytic/property/error/generic/sparse/global-RNG-isolation tests and
  18 source-free deterministic Octave/QETLAB fixtures covering 72
  native/wrapper assertions.
- Tier C Kraus, Choi, and superoperator representations; channel conversion,
  application, physicality diagnostics, constructors, positive maps, and
  reviewed compatibility wrappers.
- Tier D scalar measures, structured entanglement criteria, and a
  certificate-first native entanglement pipeline that preserves inconclusive
  outcomes.
- Additive `conclusion`, `is_conclusive`, `is_certified`, and `explain`
  helpers, with readable rich displays for status-bearing result types.
- Non-mutating `validate_density_matrix` diagnostics that report trace,
  Hermiticity, positivity, subsystem, and sparse-storage issues without
  repairing caller input.
- Public separability-strategy metadata and side-effect-free backend readiness
  discovery through `available_separability_strategies`, `describe_strategy`,
  and `backend_status`.
- Tier E coherence measures and rank analysis, including a documented correction
  of the pinned QETLAB coherence-rank implementation discrepancy.
- Tier E operator Schmidt, structured product-analysis,
  entanglement-of-formation, and sufficient separable-ball APIs, with six
  reviewed compatibility wrappers and 14 source-free Octave/QETLAB fixtures
  covering 68 assertions.
- Tier E strong majorization, elementary symmetric polynomial, compound-matrix,
  and additive-compound APIs with reviewed compatibility behavior. The
  source-free matrix-analysis artifact covers 22 fixtures and passes 59
  assertions.
- Structured positive-semidefinite, locally positive-semidefinite, totally
  positive, and totally nonsingular matrix predicates, including witnesses,
  three-valued numerical-boundary results, exact arithmetic, explicit sparse
  densification, combinatorial guards, and four compatibility entry points.
  Their package-owned mappings use explicit structured results and documented
  numerical boundaries; function-specific MATLAB-family oracle coverage is not
  implied.
- A 9,484-assertion full package suite passing 9,484/9,484 locally on Julia
  1.12.6 and the installed Julia 1.10.11, comprising 9,400 core
  assertions and 84/84 executable-tutorial assertions (48 existing plus 36 for
  the two new workflows), plus consistency checks over 477 public bindings and
  all 163 source-reviewed inventory rows.
- Seven deterministic executable tutorials for subsystem reductions, local
  channel noise, separability certificates, symmetric SAPPT states and
  constructive witnesses, certificate-aware entanglement analysis, seeded
  Schmidt decomposition and reconstruction, and Tiles-UPB bound entanglement.
  The exact scripts run standalone, in `Pkg.test()`, and during the strict
  documentation build.
- A research example for the symmetric state family in Phys. Rev. A 111,
  042418 (2025), including an exact 19-term separable decomposition at the
  five-qubit SAPPT threshold, same-spectrum separable/entangled
  representatives, published witness reconstruction, a decomposable NPT
  witness, and explicit GHZ phase handling.
- A benchmark harness containing 114 declared quick cases, all of which
  completed for the clean `f32dd233` baseline, the historical performance
  candidate `a50f516`, and the current dirty release-hardening tree.
  Targeted paired observations remain local diagnostics, not a stable
  comparative-performance claim or regression baseline.
- An optional JuMP extension suite passing 847/847 locally on Julia 1.12.6 and
  the installed Julia 1.10.11.
- An optional EntanglementDetection.jl 0.2.2 extension that runs heuristic
  searches only in an isolated child Julia process, preserves candidate results
  as uncertified evidence, and contains backend failures as `:unknown` reports.
  Its 141 focused assertions cover load order, caller-state isolation, bounded
  termination, malformed IPC, bounded live output, inherited descriptors, and
  post-launch output-capture failure cleanup.
- Task-first onboarding, grouped documentation navigation, progressive
  tutorials, copyable examples, and a GitHub Pages deployment for the rendered
  equations and browser-local code generator.
- A schema-2 browser code generator with 11 bounded state families and nine
  curated workflows, including exact Tiles-UPB and explicit-seed random-state
  construction; density validation, marginals, Schmidt diagnostics, bounded
  core separability profiles, backend readiness, conservative result helpers,
  structural assertions, live resource estimates, stale-output protection,
  and strict 16 KiB versioned JSON import/export. Schema-1 configurations
  migrate without enabling new analyses.
- A development disclosure recording substantial OpenAI Codex assistance
  without treating that disclosure as evidence of the maintainer's required
  human review.
- Release-safety guards for non-one-based caller arrays, immutable subsystem
  plan lookups and Kraus collections, validated internal construction paths,
  and combinatorial Brauer-state generation.
- Exact-archive release preflights, immutable GitHub Action pins, strict
  Codecov failure handling, and a documented private-versus-General release
  checklist.

### Changed

- Dense homogeneous BLAS-float Kraus and operator-sum Choi construction now
  uses compact factor-column products, while mixed, exact,
  arbitrary-precision, and sparse inputs retain the prior path.
- Second- and third-order compound matrices now avoid per-minor matrix
  allocation for standard BLAS floating-point and exact types. Floating-point
  kernels retain partial pivoting; other numeric types and larger minors retain
  the standard-library determinant path.
- Randomized native and `MATLABCompat` APIs require a leading explicit
  `rng::AbstractRNG`; no public random constructor draws from Julia's global
  stream.
- Documentation equations use GitHub-compatible roman-text notation instead of
  the unsupported operator-name macro.

### Security

- Added a private vulnerability-reporting policy.
