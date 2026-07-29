# Build environment

Observed locally on 2026-07-29. This is environment evidence, not a statement
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

The core project policy minimum is Julia 1.10. The full 2,318-assertion local
package corpus passes on Julia 1.12.6 and Julia 1.10.11, but these two local
versions do not substitute for remote CI on supported platforms or multiple
thread counts.

The isolated EntanglementDetection.jl 0.2.2 test environment has an effective
Julia 1.11 resolver floor because compatible Ket 0.9 releases require Julia
1.11. Its 125-assertion focused suite was recorded locally on Julia 1.12.6.
Julia 1.11 is therefore a dependency-resolution floor and configured CI target,
not a locally cited pass in this report. The core package continues to support
Julia 1.10.

## Development and oracle tools

| Tool | Location/version | Status |
|---|---|---|
| Git | 2.52.0 | Present |
| GNU Make | 3.81 | Present |
| Graphviz `dot` | 12.2.1 | Present; docs diagrams not yet validated |
| GNU Octave | 11.3.0, x86_64 build | Present; not accepted as a QETLAB oracle by default |
| MATLAB | — | Not found on `PATH` |
| EntanglementDetection.jl | Exact 0.2.2 in an isolated test environment | 125/125 focused assertions pass on Julia 1.12.6; remote Julia 1.11/1.12 platform matrix pending |
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

- Package: the integrated 2,318-assertion corpus, comprising 2,270 core
  assertions plus 48 executable-tutorial assertions, passed under Julia 1.12.6
  and Julia 1.10.11.
- Executable tutorials: the standalone 48-assertion runner passed under Julia
  1.12.6 and Julia 1.10.11; the same scripts run from `Pkg.test()` and as live
  documentation examples.
- Optional integration: the exact EntanglementDetection.jl 0.2.2 focused suite
  passed 125/125 on Julia 1.12.6. Searches use child-process isolation, bounded
  reads, and explicit cleanup; backend candidates always remain `unknown` and
  uncertified. The configured remote platform matrix has not run.
- Documentation: strict Documenter build, doctests, and live tutorial examples
  passed on Julia 1.12.6 and Julia 1.10.11 with Documenter 1.17.0; the generated
  API page emitted only a non-failing size warning.
- Quality: Aqua passed 11/11 and 25 representative JET probes passed. Those
  probes do not yet include the matrix-predicate slice.
- API consistency: the public-API/provenance gate passed over 214 public
  bindings.
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
- optional-backend active project/manifest, child-process launch and timeout
  policy, read-size bounds, termination/reaping outcome, and captured failure
  metadata;
- benchmark warmup, samples, allocations, and problem dimensions.

The EntanglementDetection.jl adapter's Julia `Serialization` channel is trusted
local worker IPC, not a security boundary for untrusted peers or payloads.
