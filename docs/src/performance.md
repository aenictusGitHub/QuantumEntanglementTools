# Performance

The project currently makes no performance claim.

## Benchmark contract

Performance work should compare a clear Julia reference, the optimized Julia
path, and QETLAB/MATLAB where available. Record:

- repository and dependency revisions;
- hardware, OS, Julia/MATLAB and BLAS versions;
- Julia, BLAS, and solver thread counts;
- input dimensions, density/sparsity, structure, and numeric type;
- warm runtime separately from first-call latency;
- samples, allocations, and peak memory where measurable;
- numerical error relative to the chosen reference.

Any statement that Julia is faster must be limited to the measured operation,
sizes, environment, and revision. Correctness tolerances may not be weakened to
improve a benchmark.

## Optimization priorities

Measure before selecting reshape/permutation, contraction, sparse accumulation,
matrix-free, or explicit loop kernels. Use views, `mul!`, factorizations, and
preallocation only where they preserve semantics and show a material benefit.
Plans and caches need thread-safety and bounded-memory review.

See the top-level `docs/BENCHMARK_REPORT.md` for the current evidence ledger.
