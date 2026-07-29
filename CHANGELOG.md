# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/) together with
Julia's pre-`1.0` compatibility convention: breaking public-API changes require
a minor-version increment.

## [Unreleased]

No user-visible changes yet.

## [0.1.0] - 2026-07-29

### Added

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
  Their focused suites pass 170 native and 37 compatibility assertions; no
  MATLAB-family predicate oracle has been run. The `IsPSD` mapping remains
  partial because the pinned CVX symbolic branch is omitted.
- A 2,417-assertion full package suite passing locally on Julia 1.12.6 and
  Julia 1.10.11, comprising 2,369 core assertions and 48 executable-tutorial
  assertions, plus consistency checks over 214 public bindings and all 163
  source-reviewed inventory rows.
- Five deterministic executable tutorials for subsystem reductions, local
  channel noise, separability certificates, symmetric SAPPT states and
  constructive witnesses, and certificate-aware entanglement analysis. The
  exact scripts run standalone, in `Pkg.test()`, and during the strict
  documentation build.
- A research example for the symmetric state family in Phys. Rev. A 111,
  042418 (2025), including an exact 19-term separable decomposition at the
  five-qubit SAPPT threshold, same-spectrum separable/entangled
  representatives, published witness reconstruction, a decomposable NPT
  witness, and explicit GHZ phase handling.
- A 42-case quick benchmark smoke suite and quality checks in which Aqua passes
  11 assertions and 25 representative JET probes pass. The JET set does not
  cover the matrix predicates.
- An optional EntanglementDetection.jl 0.2.2 extension that runs heuristic
  searches only in an isolated child Julia process, preserves candidate results
  as uncertified evidence, and contains backend failures as `:unknown` reports.
  Its 125 focused assertions cover load order, caller-state isolation, bounded
  termination, malformed IPC, and post-launch output-capture failure cleanup.
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

- Randomized native and `MATLABCompat` APIs require a leading explicit
  `rng::AbstractRNG`; no public random constructor draws from Julia's global
  stream.
- Documentation equations use GitHub-compatible roman-text notation instead of
  the unsupported operator-name macro.

### Security

- Added a private vulnerability-reporting policy.

[Unreleased]: https://github.com/aenictusGitHub/QuantumEntanglementTools/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/aenictusGitHub/QuantumEntanglementTools/releases/tag/v0.1.0
