# Random PPT states

`random_ppt_state` constructs a bipartite density matrix whose requested
partial transpose is positive semidefinite. The random-number generator and
local dimensions are mandatory.

```julia
using QuantumEntanglementTools
using Random

result = random_ppt_state(
    Xoshiro(20260730),
    (2, 3);
    ranks=(3, 4),
)

@assert result.status == RandomPPTConstructed
@assert result.verified
@assert result.numerical_ranks[1] <= 3
@assert result.numerical_ranks[2] <= 4
@assert abs(tr(result.state) - 1) <= result.tolerance
```

The global random stream is never read or mutated. Repeating the call with a
new generator initialized from the same seed reproduces the state.

## Constructions

`construction=:auto` selects one of two bounded mathematical constructions:

- `:shifted_induced` for full requested ranks. A Ginibre Gram state is shifted
  by an explicit identity multiple sufficient to give its partial transpose a
  positive margin, then normalized once as part of the construction.
- `:separable_mixture` for low requested ranks. It draws
  `min(ranks...)` normalized product vectors and a simplex weight vector. Their
  convex mixture is separable, hence PPT, and both its rank and its
  partial-transpose rank are at most the requested bounds.

The pinned QETLAB contract promises neither a named distribution nor
entanglement, and specifies ranks as upper bounds. The guaranteed
separable-mixture route therefore covers the documented low-rank capability
without QETLAB's potentially unbounded pseudoinverse iteration.

The result identifies the selected construction and records the requested and
numerical ranks, iteration count, convergence history, tolerance, work
estimate, allocation limit, and whether normalization was part of the
construction.

## Verification and failure behavior

A candidate is promoted to `result.state` only after checking:

- Hermiticity;
- positive semidefiniteness;
- unit trace;
- Hermiticity and positive semidefiniteness after partial transpose;
- both requested numerical rank bounds.

If any check fails, `state === nothing`; the unverified matrix is retained only
as `candidate`, with `RandomPPTVerificationFailed`. There is no fallback
spectral clipping, PSD projection, Hermitian averaging, or post-failure
normalization.

`max_dense_entries` and `max_work` are checked before consuming the supplied
RNG. `max_iterations=0` similarly returns `RandomPPTIterationLimit` without a
draw. The current direct constructions need exactly one bounded step.
Float32, Float64, real, and complex output paths are supported.

Malformed dimensions, ranks larger than the total dimension, unsupported
scalar types, and incompatible construction choices throw an input error. A
rank-limited request cannot select `:shifted_induced`, because that route is
full rank by construction.

## Upstream relationship

The full-rank shift is source-informed by QETLAB `RandomPPTState` at pinned
revision `d8589610f00cff106537268dee2e2a1153f3a601`. The bounded low-rank route
is an intentional Julia-native replacement for the upstream unbounded
iteration, whose source cites [Leinaas, Myrheim, and
Sollid](https://arxiv.org/abs/1002.1949). Source-free cross-language evidence
compares properties rather than random entries, since Julia and MATLAB/Octave
do not share an RNG algorithm.
