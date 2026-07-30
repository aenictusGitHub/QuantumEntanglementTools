# Bell inequalities and nonlocal games

This interface separates exact finite enumeration, numerical upper
relaxations, and attained numerical candidates. A solver result is never silently
rounded to a theorem:

- exact classical searches put the optimum in `value` and set `exact=true`;
- no-signalling and finite-NPA solves put a numerical relaxation value in
  `upper_bound`, with the solver record and witness retained;
- the see-saw puts its independently re-evaluated, tolerance-feasible
  candidate value in `lower_bound` without marking it certified;
- unavailable backends, resource limits, and invalid certificates have
  distinct statuses and do not masquerade as a negative result.

The full-probability axis order is always
`P[a, b, x, y]`: Alice outcome, Bob outcome, Alice setting, Bob setting.
Question distributions use `p[x, y]`. Constructors reject inconsistent
dimensions, nonfinite data, and probabilities outside the documented
tolerance. They do not normalize, clip, symmetrize, or otherwise repair
inputs.

## A solver-free CHSH calculation

The classical XOR route uses a guarded exhaustive search and returns a
deterministic strategy as its certificate.

```@example nonlocal-core
using QuantumEntanglementTools
const QET = QuantumEntanglementTools

question_probability = fill(1 // 4, 2, 2)
winning_parity = [0 0; 0 1]
classical = QET.xor_game_value(
    question_probability,
    winning_parity;
    regime=:classical,
    max_strategies=16,
)

(
    status=classical.status === QET.NonlocalValueExact,
    value=classical.value,
    exact=classical.exact,
    certificate=classical.certificate_kind,
)
```

For a general Bell functional, construct a scenario and state the notation
explicitly. The familiar full-correlator CHSH coefficients have the constant
and marginal terms in row and column one.

```@example nonlocal-core
scenario = QET.BellScenario(2, 2, 2, 2)
chsh = [
    0  0  0
    0  1  1
    0  1 -1
]
bell = QET.bell_inequality_bound(
    chsh,
    scenario;
    notation=:full_correlator,
    regime=:classical,
    max_strategies=16,
)

(value=bell.value, strategy=bell.strategy)
```

`BellFunctional` also accepts `:full_probability`/`:fp` and
`:collins_gisin`/`:cg`. Behavior conversions are represented by the distinct
`FullProbabilityBehavior` and `CollinsGisinBehavior` types, preventing a
coefficient tensor from being mistaken for a probability behavior.

## NPA membership

`npa_membership(behavior; level=0)` checks positivity, normalization, and
no-signalling without a solver. Only this exact branch has a Boolean verdict.

```@example nonlocal-core
probabilities = zeros(Rational{Int}, 2, 2, 2, 2)
probabilities[1, 1, :, :] .= 1
behavior = QET.FullProbabilityBehavior(
    probabilities, scenario; atol=0, rtol=0
)
basic = QET.npa_membership(behavior; level=0)

(
    status=basic.status === QET.NPABasicConditionsSatisfied,
    verdict=basic.verdict,
    certified=basic.certified,
)
```

Positive levels build a package-owned `SemidefiniteProgram`. Integer levels
include every canonical projector word up to that length. Strings such as
`"1+ab"` add the indicated intermediate words. Alice and Bob letters commute
with one another; same-setting projectors are idempotent or orthogonal.
`max_words` and `max_word_generation_work` are mandatory practical guards
against the exponential word catalog.

For a numerical membership solve, inspect `status`, `optimization_result`,
`residuals`, and `moment_matrix`. A finite-level feasible moment matrix is not
proof that a behavior is quantum, so `verdict` remains `nothing`. Likewise,
backend-reported infeasibility is retained as numerical evidence rather than
being converted automatically into an exact exclusion certificate.

## Optional solver-backed bounds

Install an optimization backend in the optional JuMP environment. For
example, with Hypatia:

```julia
using Hypatia
using QuantumEntanglementTools

backend = JuMPBackend(
    Hypatia.Optimizer;
    optimizer_name="Hypatia",
    optimizer_version=Base.pkgversion(Hypatia),
    allow_densify=true,
)

quantum = xor_game_value(
    fill(0.25, 2, 2),
    [0 0; 0 1];
    regime=:quantum,
    backend=backend,
)
quantum.upper_bound
quantum.optimization_result
```

The XOR quantum formulation is the Tsirelson Gram SDP. General quantum Bell
bounds use the requested finite NPA level. No-signalling bounds use a linear
program with positivity, normalization, and both parties' marginal
equalities. Solver-derived upper values are deliberately marked
`certified_upper=false`: they are useful numerical bounds, but not
interval-arithmetic or exact-dual certificates.

All programs accept `OptimizationLimits`. Backend status, primal and dual
status, residuals, objective values, and objective bounds remain available in
the structured result.

## Explicit-RNG attained numerical candidates

`nonlocal_game_lower_bound` is a bounded alternating SDP see-saw. It requires
an `AbstractRNG`; it never reads or mutates Julia's global random stream.

```julia
using Random

payoff = zeros(Float64, 2, 2, 2, 2)
for a in 1:2, b in 1:2, x in 1:2, y in 1:2
    payoff[a, b, x, y] = a == b
end
game = NonlocalGame(fill(0.25, 2, 2), payoff)

lower = nonlocal_game_lower_bound(
    MersenneTwister(41),
    2,
    game;
    backend=backend,
    max_iterations=20,
    max_initialization_attempts=4,
)
lower.lower_bound
```

Each candidate is reconstructed from the returned state and POVMs and then
re-evaluated independently. The conic residual checks establish only numerical
feasibility within the selected tolerances, so `certified_lower=false`: the
result is useful attained candidate evidence, but it is not an exact
mathematical lower-bound certificate. If the iteration limit is reached, the
best candidate is preserved with a resource-limit status.

## Binary constraint-system games

A BCS constraint is a binary array with one length-two axis per variable. The
conversion is exact: question probabilities use rational arithmetic, and
Alice's answer labels follow the pinned helper's most-significant-variable
bit order.

```@example nonlocal-core
constraints = [
    [1 0; 0 1],
    [0 1; 1 0],
]
bcs = QET.BCSGame(constraints)
converted = QET.nonlocal_game(bcs)
bcs_classical = QET.bcs_game_value(
    bcs; regime=:classical, max_strategies=256
)

(
    scenario=Tuple(converted.scenario),
    question_mass=sum(converted.probabilities),
    value=bcs_classical.value,
    status=bcs_classical.status === QET.NonlocalValueExact,
)
```

The pinned `BCSGameValue.m` calls `NonlocalGameValue.m`, but that dependency is
absent from the pinned QETLAB tree. `bcs_game_value` therefore uses an
independently specified dispatcher over the reviewed classical,
no-signalling, and NPA implementations. This is explicit missing-support
reconstruction, not a claim that the unavailable upstream routine was
executed.

`bcs_game_lower_bound` converts the exact rational game to a numerical solver
type through the public `coefficient_type` keyword. Only `Float32` and
`Float64` are accepted; no hidden precision conversion is performed.

## Fixed-qubit PPT relaxation

`bell_inequality_qubit_bound` accepts a rectangular joint coefficient matrix,
separate Alice and Bob marginal coefficients, and exactly two outcome values
for each party. The pinned MATLAB routine assigns Alice's setting count to
Bob's loops and dimensions. The Julia model corrects this defect and records
`upstream_rectangular_loop_corrected=true` in its diagnostics.

The distinguished Alice/Bob clone pair is grouped as a four-dimensional
subsystem, preserving their shared entanglement. Remaining clones are
two-dimensional, with one PPT constraint per complementary bipartition.
`max_dimension` and `max_ppt_constraints` guard the exponential relaxation.
The returned `relaxation_state` is the PPT model variable, not a physical
two-qubit state or an attained strategy, and its numerical objective is not
reported as an exact value.

## MATLAB-compatible entry points

The compatibility namespace preserves QETLAB's argument order while making
unsafe output conversion opt-in:

```julia
MATLABCompat.NPAHierarchy(cg, desc, k=1; structured=true)
MATLABCompat.NonlocalGameLB(rng, d, p, V, verbose=1; structured=true)
MATLABCompat.XORGameValue(p, f, vtype="classical"; structured=true)
MATLABCompat.BellInequalityMax(coefficients, desc, notation,
                              mtype="classical", k=1; structured=true)
MATLABCompat.BellInequalityMaxQubits(joint, alice, bob,
                                    alice_values, bob_values; structured=true)
MATLABCompat.BCSGameLB(rng, d, constraints, verbose=1; structured=true)
MATLABCompat.BCSGameValue(constraints, mtype="classical", k=1; structured=true)
```

The RNG argument added to both lower-bound wrappers is mandatory. With
`structured=false`, exact values may be returned as scalars. Numerical
see-saw candidates, inconclusive or unavailable outcomes, and uncertified
numerical upper bounds throw instead of discarding their status. The
fixed-qubit wrapper returns its legacy tuple only if a certified upper bound
exists.

## Cost, storage, and evidence boundary

General deterministic searches are exponential in the number of settings;
the XOR reduction enumerates the smaller party and takes a best response for
the other. BCS conversion creates up to `2^k` Alice outputs when a constraint
contains `k` active variables. NPA word catalogs and fixed-qubit PPT matrices
also grow exponentially. The `max_strategies`, `max_work`, `max_entries`,
`max_words`, `max_word_generation_work`, `max_dimension`,
`max_ppt_constraints`, `max_iterations`, `max_initialization_attempts`,
`max_dense_entries`, and `max_initialization_work` controls should be chosen
deliberately. See-saw dense-storage, initialization-work, and conic-model
limits are checked before the RNG is consumed.

Sparse question distributions remain sparse. Four-dimensional payoff tensors
are not silently converted to a two-dimensional sparse format. Solver-side
densification is allowed only when the selected backend explicitly permits it
and is guarded by optimization limits. Public behaviors, games, strategies,
solver witnesses, and returned operator families own read-only array views;
`copy` produces a mutable caller-owned array when needed.

The committed source-free Octave fixture exercises two classical XOR games,
three classical Bell functionals, and the BCS conversion helper at QETLAB commit
`d8589610f00cff106537268dee2e2a1153f3a601`. It does not claim CVX evidence for
NPA, no-signalling, see-saw, or fixed-qubit models. The optional Julia tests
exercise those package-owned models with Hypatia and SCS.

The implementation follows the finite-level moment-relaxation construction of
[Navascués, Pironio, and Acín](https://arxiv.org/abs/0803.4290), the XOR-game
SDP characterization discussed by
[Cleve, Høyer, Toner, and Watrous](https://arxiv.org/abs/quant-ph/0404076),
and the finite-dimensional hierarchy referenced by the pinned qubit routine,
[Navascués, de la Torre, and Vértesi](https://arxiv.org/abs/1308.3410).
