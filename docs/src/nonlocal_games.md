# Nonlocal games

This page currently documents the solver-free tensor operation used to form
parallel repetitions. The optimization-dependent Bell, XOR, BCS, and NPA
rows remain tracked separately by the generated completion ledger; a tensor
constructor is not evidence that any classical or quantum game value has been
computed.

## Full-probability coefficient order

A bipartite game coefficient tensor has axes

```math
V[a,b,x,y],
```

in the order Alice output, Bob output, Alice setting, Bob setting. For `r`
parallel copies, [`parallel_repetition`](@ref) returns a dense tensor of shape

```math
(o_A^r,o_B^r,m_A^r,m_B^r)
```

whose entries are products of the corresponding one-copy coefficients. Each
packed axis follows the package's tensor convention: copy one is the
most-significant digit and the last copy varies fastest.

```@example parallel-repetition
using QuantumEntanglementTools

game = reshape(1:16, 2, 2, 2, 2)
repeated = parallel_repetition(game, 2; max_entries=1_000)

packed(i, j, base) = (i - 1) * base + j
entry = repeated[
    packed(1, 2, 2),
    packed(2, 1, 2),
    packed(1, 1, 2),
    packed(2, 2, 2),
]
expected = game[1, 2, 1, 2] * game[2, 1, 1, 2]

(size=size(repeated), entry=entry, agrees=entry == expected)
```

The operation is exact for exact coefficient types and does not normalize or
repair its input. A one-copy request returns an owned copy. Nonpositive copy
counts are rejected, correcting the pinned MATLAB routine's accidental
unchanged-input behavior for such values.

## Work and storage limits

The result is necessarily dense in the standard-library implementation because
Julia's standard sparse arrays are two-dimensional. Before allocation,
`max_entries` checks the output element count and `max_work` checks the
conservative estimate `r * max_entries_required`. Both limits can be disabled
with `nothing`, but doing so is an explicit acceptance of exponential storage
and work.

[`MATLABCompat.ParallelRepetition`](@ref) preserves QETLAB's argument order and
delegates to the same native kernel and guards. It deliberately retains the
native positive-integer validation rather than reproducing the invalid
nonpositive-copy behavior.

A committed source-free Octave 11.3.0 fixture from the exact pinned QETLAB
revision covers both the one-copy and two-copy branches. It passes eight
native/compatibility comparisons and records the full four-dimensional output
shape and flattened coefficients. This is function-specific supplemental
evidence; MATLAB was not run.

The implementation uses $O(r\prod_i d_i^r)$ scalar work and one dense output
allocation of $\prod_i d_i^r$ elements.
