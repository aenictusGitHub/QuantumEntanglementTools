# QETLAB inventory source review

Last updated: 2026-07-29.

This report records the source-level disposition of every QETLAB row that was
previously awaiting review. The comparison target is the pinned QETLAB revision
`d8589610f00cff106537268dee2e2a1153f3a601`. The detailed, machine-readable
reason for each row is in [`porting/qetlab_status.toml`](../porting/qetlab_status.toml);
the generated result is in
[`porting/qetlab_inventory.toml`](../porting/qetlab_inventory.toml).

“Source-reviewed” means that the pinned MATLAB body, dependencies, current
Julia API, provenance, tests, and documentation were compared and a deliberate
disposition was recorded. It is not a claim that a human maintainer has
completed the non-delegable release review, and it is not synonymous with
implemented.

## Result

All 163 rows now have reviewed dispositions; none remains in the generated
`automated_parse_pending_manual_review` state.

The 127 public QETLAB rows are classified as:

| Disposition | Rows | Meaning |
|---|---:|---|
| Implemented | 63 | Existing reviewed mappings with recorded implementation evidence |
| Partial | 15 | A named, documented subset exists; omitted inputs, branches, or compatibility behavior remain explicit |
| Deferred | 19 | Implementable work is outside the current release scope; two were already deferred and 17 were classified in this sweep |
| Blocked with explicit reason | 30 | A solver architecture, missing pinned dependency, prerequisite API, or unresolved result/certificate contract is required |
| Pending | 0 | No public row remains automatically classified |

The 36 private helper rows are tracked separately:

| Helper disposition | Rows |
|---|---:|
| Replaced by Julia/Base or a tested private parent-specific algorithm | 13 |
| Intentionally excluded because the pinned private helper has no caller | 2 |
| Partially covered in current parents but still needed by deferred parents | 4 |
| Deferred with the blocked or deferred parent scope | 17 |
| Pending | 0 |

Private helpers are not public API promises. They remain in the inventory for
dependency, provenance, and license review.

## Newly reviewed public rows

The 51 public rows reviewed in this sweep have the following dispositions.

### Partial

| Function | Implemented subset and explicit remainder |
|---|---|
| `IsCP` | `is_completely_positive` covers typed Kraus, Choi, and superoperator representations; no `MATLABCompat.IsCP`, raw two-sided Kraus input, or QETLAB tolerance parity is claimed |
| `IsSeparable` | The certificate-first native pipeline covers PPT, realignment, reduction, exact low-dimensional PPT sufficiency, pure-state Schmidt certificates, and a separate separable-ball test; the composite QETLAB heuristic/solver pipeline is not ported |
| `kpNorm` | Verified Schatten, trace, and Ky Fan APIs cover the full-spectrum and top-`k`, `p=1` special cases; no general top-`k`, `p` binding or CVX branch exists |
| `SkVectorNorm` | The value can be composed from the largest `k` values returned by `schmidt_coefficients`; no dedicated entry point or compatibility alias exists |

### Blocked with an explicit reason

| Area | Functions | Common blocker |
|---|---|---|
| Absolute PPT and S(`k`) criteria | `AbsPPTConstraints`, `IsAbsPPT`, `IsBlockPositive`, `SkOperatorNorm` | Status-aware SDP/constraint infrastructure, bounded combinatorics, and certified bound/witness semantics |
| Filter and norm dependencies | `FilterNormalForm`, `kpNormDual` | Unimplemented prerequisite algorithms and stable convergence/duality contracts |
| UPB optimization | `UPBSepDistinguishable` | UPB enumeration plus an optional solver with feasibility and unknown-result semantics |
| Channel and state optimization | `CBNorm`, `ChannelDistinguishability`, `DiamondNorm`, `Distinguishability`, `LocalDistinguishability`, `MaximumOutputFidelity` | Diamond-norm or state-discrimination SDPs, solver statuses, and reviewed map/measurement result types |
| Extension hierarchies | `SymmetricExtension`, `SymmetricInnerExtension` | Optional SDP backend and explicit extension, separating-witness, dual, and inconclusive outcomes |
| Nonlocal games | `BCSGameLB`, `BCSGameValue`, `BellInequalityMax`, `BellInequalityMaxQubits`, `NonlocalGameLB`, `NPAHierarchy`, `XORGameValue` | CVX/NPA formulations, exponential work guards, explicit RNG for heuristics, and one dependency absent from the pinned tree |
| Polynomial optimization | `CliqueNumber`, `IsCopositive`, `PolynomialSOS` | SOS/SDP backend, reproducible randomized bounds, and tri-state bound semantics |
| Coherence optimization | `GenRobustnesskCoherence`, `IsAbskIncoh`, `IskIncoherent`, `RobustnessCoherence`, `TraceDistanceCoherence` | CVX branches, upstream formula/helper defects, and status-aware optimizer or closest-state results |

Blocked does not mean impossible. It means that marking the row implemented
before the named prerequisite exists would violate the package's solver,
certificate, RNG, or validation rules.

### Deferred

| Area | Functions | Required work |
|---|---|---|
| Entanglement and UPBs | `EntangledSubspace`, `IsEntanglingGate`, `IsUPB`, `MinUPBSize`, `UPB` | Guarded constructions, theorem/certificate specifications, structured witnesses, and family-by-family validation |
| Matrix analysis | `InducedMatrixNorm`, `MatsumotoFidelity` | Exact-versus-bound results, explicit RNG and convergence for iterations, and stable support-aware matrix geometry |
| Channels and general information | `InducedSchattenNorm`, `IsHermPreserving`, `Commutant`, `OperatorSinkhorn`, `ParallelRepetition`, `Twirl` | Bounded iterative/combinatorial designs, sparse policies, standalone diagnostics, and invariant tests |
| Polynomial and coherence formulas | `CopositivePolynomial`, `PolynomialAsMatrix`, `PolynomialOptimize`, `RobkCohValue` | A typed polynomial representation, exact coefficient conventions, deterministic budgets, and corrected theorem/input specifications |

The two previously reviewed deferred public rows are `RandomSuperoperator` and
`RandomPPTState`, bringing the public deferred total to 19.

## Internal helper dispositions

The following helpers are replaced by Base/stdlib operations or tested private
algorithms inside their reviewed parent implementations:

`asum_vector`, `asymind`, `dec_to_bin`, `glob_ind`, `iden`, `opt_args`,
`opt_disp`, `pad_array`, `perfect_matchings`, `perm_inv`, `perm_sign`,
`pure_to_mixed`, and `sporth`.

`chshd` and `ffl` have no caller in the pinned tree and are intentionally not
made package APIs.

The currently covered part is recorded, but deferred parents still require
additional semantics for `sum_vector`, `superoperator_dims`, `symind`, and
`update_odometer`.

The following helpers remain coupled to a deferred or blocked parent:

`bcs_to_nonlocal`, `CG2FC`, `CG2FP`, `exp2ind`, `FC2CG`, `FC2FP`, `FP2CG`,
`FP2FC`, `has_band_k_ordering`, `jacobi_poly`, `normalize_cols`,
`one_factorization`, `poly_rand_input`, `sk_iterate`, `spnull`, `symindfind`,
and `vec_partitions`.

The pinned `spnull` and `sporth` files contain separate Bruno Luong attribution;
that attribution must be preserved if source-derived code is ever ported.

## What this closes—and what it does not

This sweep closes the inventory-classification gap and replaces the misleading
“87 pending” number with an explicit roadmap. It does not complete QETLAB
parity: 64 public rows are still partial, deferred, or blocked. Advancing any
one of them requires the function-level specification, provenance,
implementation, tests, and documentation required by `AGENTS.md`.
