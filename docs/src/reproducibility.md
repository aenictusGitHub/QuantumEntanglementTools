# Reproducibility

A reproducible numerical result records enough context to repeat the same
algorithm and interpret differences.

## Required metadata

- package commit/version and dirty-worktree state;
- Julia version, operating system, architecture, BLAS vendor, and thread counts;
- dependency manifest;
- numeric element type, dimensions, normalization, and tolerance policy;
- algorithm, backend, solver, versions, options, statuses, and residuals;
- RNG type and explicit seed for randomized work;
- source fixture provenance and comparison rule.

Randomized package APIs must accept an explicit `rng::AbstractRNG` and must not
change the process-global random stream.

For example:

```julia
using QuantumEntanglementTools
using Random

seed = 0x514554
first_run = random_density_matrix(Xoshiro(seed), 4; rank = 2)
second_run = random_density_matrix(Xoshiro(seed), 4; rank = 2)
@assert first_run == second_run
```

Passing the same mutable RNG through several calls intentionally advances that
stream. Reconstruct the seeded RNG to replay a call. The local Tier B
regression also verifies that all six native constructors and their
`MATLABCompat` wrappers leave Julia's global stream unchanged.

## Comparisons

Use exact comparison for combinatorial data, scale-aware numeric tolerances for
floating-point data, projectors/subspaces for phase- or basis-ambiguous
eigenvectors, sorted spectra where order is unspecified, and objectives,
residuals, and certificates rather than bitwise solver vectors.

MATLAB is the intended behavioral oracle where available, but it is not assumed
in public CI. Octave may supplement a check only after function-specific
compatibility is demonstrated. Neither replaces analytic fixtures,
property-based tests, or independent formulations.

The committed Tier B artifact
`test/oracle/fixtures/tier_b_octave_11_3_qetlab_d858961.json` contains 18
deterministic operator/state fixtures and drives 72 native/wrapper assertions.
It does not contain randomized fixtures: equal seeds do not define the same
stream across Julia and MATLAB-family engines. Octave 11.3.0 is recorded as the
generating engine and the artifact explicitly marks itself as not generally
MATLAB-compatible.
