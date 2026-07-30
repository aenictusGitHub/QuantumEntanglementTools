# Build environment

Current convergence evidence observed locally on 2026-07-30. Environment
details and predecessor results dated 2026-07-29 are retained below as a
historical baseline.

## Current convergence environment

<!-- qetlab-current-claims: begin -->

| Item | Current local evidence |
|---|---|
| Repository | Dirty `main` worktree based on `f32dd233e478dd6e2642f11fab088f6c8febc420`; `origin/main` is the same base; no commit or publication action |
| Julia | 1.12.6 and installed minimum-line binary 1.10.0 |
| Completion ledger | 127/127 public rows are verified with final status; completion queue contains 0 public rows; 36/36 internal helpers have terminal dispositions; 0 required helpers remain; 0 static completion failures |
| API/provenance | 458 exports and 458 matching provenance records |
| Package corpus | 8,133/8,133 assertions including 48 tutorials on both installed Julia lines |
| JuMP optimization | 836/836 assertions with package-managed Hypatia/SCS on both installed Julia lines |
| EntanglementDetection.jl | Exact 0.2.2 environment; 125/125 assertions on Julia 1.12.6; effective resolver floor Julia 1.11 |
| Source-free oracles | 27/27 comparators and 1,347/1,347 assertions on both installed Julia lines; all fixture digests verified |
| Quality | Aqua 11/11, JET 25/25, formatter gate passing, randomized matrix predicates 130/130 on both lines |
| Benchmarks | All 114 declared quick benchmark cases completed for a clean `f32dd233` baseline and the candidate with one Julia and one BLAS thread; targeted paired results are diagnostic only |
| Distribution preflight | Isolated dirty-worktree archive and fresh-depot load smoke pass on Julia 1.12.6 and 1.10.0; not committed-tree release evidence |

`Pkg.test()` initially encountered sandbox-only permission failures when Julia
attempted to write `~/.julia/logs/manifest_usage.toml.pid`; the approved
outside-sandbox reruns passed. Julia 1.10 warned that the ignored local
manifest had been resolved by Julia 1.12 and that project requirements had
changed; `Pkg.test()` generated an isolated temporary test environment and
passed.

The optimization evidence uses Julia packages rather than command-line solver
executables. A successful numerical solve is not automatically a mathematical
certificate. These are local dirty-worktree results, not MATLAB parity,
supported-platform CI, comparative performance, API stability, release
approval, or human review.

<!-- qetlab-current-claims: end -->

## Historical environment and evidence (2026-07-29)

## Evidence boundaries

At the start of the convergence audit, local `main`, `origin/main`, and `HEAD`
were clean and equal to
`2d965bfbbcb6f3af350e293f1034b1fdcd4f937e`. The repository now contains
uncommitted Phase 0 improvements. The results below are the recorded local
comparison baseline; they do not certify the current prospective worktree or a
future release commit. Exact local and remote reruns remain required after the
changes are committed.

The repository is private, and no tag, GitHub release, or GitHub Pages
deployment exists. See [`CONVERGENCE_AUDIT.md`](CONVERGENCE_AUDIT.md) for the
separation between exact baseline evidence, prospective local work, and
external release gates.

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

The core project policy minimum is Julia 1.10. The full 2,417-assertion local
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
| EntanglementDetection.jl | Exact 0.2.2 in an isolated test environment | 125/125 focused assertions pass on Julia 1.12.6; an exact rerun remains required for the commit produced by the current convergence work |
| GLPK/`glpsol` | 5.0 | Present; not an SDP solver and not validated for package APIs |
| Other solver executables checked | SCS, CSDP, SDPA, Mosek, Gurobi, CBC, HiGHS, Ipopt | Not found on `PATH` |

No MATLAB, CVX, or solver-backed certification results are available from this
environment snapshot. The source-free differential fixture corpus contains 338
passing assertions against QETLAB commit
`d8589610f00cff106537268dee2e2a1153f3a601`. This function-specific evidence
does not make Octave a generally equivalent MATLAB oracle. Package-managed
Julia solver libraries may differ from command-line executables and must be
reported by each test environment.

## Local checks recorded

- Package: the integrated 2,417-assertion corpus, comprising 2,369 core
  assertions plus 48 executable-tutorial assertions, passed under Julia 1.12.6
  and Julia 1.10.11.
- Executable tutorials: the standalone 48-assertion runner passed under Julia
  1.12.6 and Julia 1.10.11; the same scripts run from `Pkg.test()` and as live
  documentation examples.
- Optional integration: the exact EntanglementDetection.jl 0.2.2 focused suite
  passed 125/125 on Julia 1.12.6. Searches use child-process isolation, bounded
  reads, and explicit cleanup; backend candidates always remain `unknown` and
  uncertified. An exact remote rerun remains required after the current
  prospective changes are committed.
- Inventory: the generator passed with 163 source rows, 503 dependency edges,
  and zero automatically detected cycles against QETLAB commit
  `d8589610f00cff106537268dee2e2a1153f3a601`.
- Documentation: strict Documenter build, doctests, and live tutorial examples
  passed on Julia 1.12.6 and Julia 1.10.11 with Documenter 1.17.0. The API
  reference is split into native and compatibility pages; the 157 KiB native
  page remains below the 200 KiB hard limit and emits only a non-failing size
  warning. A Julia 1.12.6 CI-mode build also produced the expected pretty-URL
  generator page and local assets.
- Code generator: JavaScriptCore passed 71 deterministic assertions and emitted
  a generated Julia smoke program whose nine state-family branches passed on
  Julia 1.12.6 and Julia 1.10.11.
- Quality: Aqua passed 11/11, 25 representative JET probes passed, and the
  independent matrix-predicate validation passed 130/130 assertions.
- API consistency: the public-API/provenance gate passed over 214 exports and
  214 provenance records.
- Differential fixtures: 338/338 assertions passed against the exact QETLAB
  source pin. This is supplemental, function-specific evidence.
- Benchmark smoke: all 42 cases ran with BenchmarkTools 1.8.0; see
  `BENCHMARK_REPORT.md`. This is not a comparative or release baseline.
- Release integrity: CFFConvert 2.0.0 validates `CITATION.cff` against schema
  1.2.0, and the offline upstream audit reports 10/10 pin/license checks clean.
  The current prospective worktree still requires an exact committed-tree
  rerun.

## Exact remote baseline and current blockers

At `2d965bfbbcb6f3af350e293f1034b1fdcd4f937e`, Core and Documentation passed.
Quality failed only because a tracked editor-workspace file made the
release-integrity check reject the tree. Coverage's package-test steps passed,
but Codecov rejected the upload with `Repository not found`.

The effective branch controls observed protect against deletion and
non-fast-forward updates. Required status checks and required pull-request
reviews are absent. The repository remains private and has no tag, GitHub
release, or Pages deployment. These are repository or service observations,
not properties of the local host.

A source-scope scan found no excluded-source identifier in the current tracked
tree or current commit archive. A history-only occurrence remains, and no
history rewrite was performed.

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
