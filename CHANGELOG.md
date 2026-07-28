# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and released versions will follow [Semantic Versioning](https://semver.org/) once
a stable public API exists.

## [Unreleased]

### Added

- Initial Julia package shell.
- Documentation, governance, legal, citation, and CI scaffolding.
- Explicit pre-alpha status and milestone ledgers.
- Preserved copy of the separately supplied QUBIT4MATLAB license text and
  documentation of its material conflict with the v6.5 archive license.
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
  Their focused suites pass 166 native and 37 compatibility assertions; no
  MATLAB-family predicate oracle has been run. The `IsPSD` mapping remains
  partial because the pinned CVX symbolic branch is omitted.
- A 2,270-assertion full package suite passing locally on Julia 1.12.6 and
  Julia 1.10.11, plus consistency checks over 212 public bindings and 76
  manually reviewed inventory rows.
- A 42-case quick benchmark smoke suite and quality checks in which Aqua passes
  11 assertions and 24 representative JET probes pass. The JET set does not
  cover the matrix predicates.

### Changed

- Randomized native and `MATLABCompat` APIs require a leading explicit
  `rng::AbstractRNG`; no public random constructor draws from Julia's global
  stream.

### Security

- Added a private vulnerability-reporting policy.
