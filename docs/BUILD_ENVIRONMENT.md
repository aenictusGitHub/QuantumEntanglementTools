# Build environment

Observed locally on 2026-07-28. This is environment evidence, not a statement
that every tool below has successfully built, tested, or validated the package.

## Host

| Item | Observed value |
|---|---|
| Operating system | macOS 26.5.2, build 25F84 |
| Kernel | Darwin 25.5.0 |
| Architecture | Apple arm64/aarch64 |
| CPU | Apple M4 Pro |
| Memory | 24 GiB |
| Shell | Bash |

## Julia

| Item | Observed value |
|---|---|
| Julia | 1.12.6, official release |
| Commit | `15346901f00` |
| LLVM | 18.1.7 |
| Word size | 64 bit |
| Julia threads | 1 default, 1 interactive, 1 GC |
| BLAS configuration | `LBTConfig([ILP64] libopenblas64_.dylib)` |

The project policy minimum is Julia 1.10. The full local package corpus also
passes on Julia 1.10.11, but these two local versions do not substitute for
remote CI on supported platforms or multiple thread counts.

## Development and oracle tools

| Tool | Location/version | Status |
|---|---|---|
| Git | 2.52.0 | Present |
| GNU Make | 3.81 | Present |
| Graphviz `dot` | 12.2.1 | Present; docs diagrams not yet validated |
| GNU Octave | 11.3.0, x86_64 build | Present; not accepted as a QETLAB oracle by default |
| MATLAB | — | Not found on `PATH` |
| GLPK/`glpsol` | 5.0 | Present; not an SDP solver and not validated for package APIs |
| Other solver executables checked | SCS, CSDP, SDPA, Mosek, Gurobi, CBC, HiGHS, Ipopt | Not found on `PATH` |

No MATLAB, CVX, or solver-backed certification results are available from this
environment snapshot. Exact Octave/QETLAB Tier A fixtures have passed for
permutation, trace, partial transpose, realignment, and
symmetric/antisymmetric projectors; this function-specific evidence does not
make Octave a generally equivalent MATLAB oracle. Package-managed Julia solver
libraries may differ from command-line executables and must be reported by each
test environment.

## Local checks recorded

- Core: the integrated 2,270-assertion package corpus passed under Julia 1.12.6
  and Julia 1.10.11.
- Documentation: strict Documenter build and doctests passed with Documenter
  1.17.0; the generated API page emitted only a non-failing size warning.
- Quality: Aqua passed 11/11 and 24 representative JET probes passed. Those
  probes do not yet include the matrix-predicate slice.
- Benchmark smoke: all 42 cases ran with BenchmarkTools 1.8.0; see
  `BENCHMARK_REPORT.md`. This is not a comparative or release baseline.

## Reproduction metadata to capture

Validation and benchmark reports must additionally record:

- exact repository commit and dirty-worktree state;
- `Manifest.toml` or equivalent dependency resolution;
- `versioninfo()` and `LinearAlgebra.BLAS.get_config()`;
- Julia and BLAS thread counts;
- RNG type and seed for randomized inputs;
- backend/solver versions, options, statuses, and residuals;
- benchmark warmup, samples, allocations, and problem dimensions.
