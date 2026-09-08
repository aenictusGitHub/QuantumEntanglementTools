# Benchmark report

Evidence date: 2026-09-08.

Status: paired local quick-run diagnostic; no stable comparative-performance
claim or regression threshold.

## Current 120-case quick smoke

<!-- qetlab-current-claims: begin -->

The complete set of 120 declared quick benchmark cases completed with exit
status zero on Julia 1.12.6 for the uncommitted performance-and-stability
worktree based on `1d611e4f61f2d740602018dce05afd47f2ecf620`:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --compiled-modules=no --startup-file=no --project=benchmark \
  benchmark/benchmarks.jl --quick \
  --output=/private/tmp/qet-performance-final-20260908
```

Each case produced 20 samples with one evaluation. The environment was an
Apple M4 (`apple-m4`) on arm64 macOS, Julia 1.12.6, BenchmarkTools 1.8.0,
and ILP64 OpenBLAS through libblastrampoline, with one Julia and one BLAS
thread. `Manifest.toml` and `benchmark/Manifest.toml` had SHA-256 digests
`6faecb9f1c414a769c09ad82d5cbcb5188ec1c36fba0a465688485631bf1245f` and
`73b333138b13ff43b3d8470bedaaed547c855a51ab1d9cc49dab19f02500aa86`.

The raw JSON/TOML SHA-256 pair is
`ceda9894530b72eddd76df145189a234a8d7368d129f5781ffb018a26c10347a` /
`78b2e4efa676a6097598c95954b28c4e8c39c2e115acd9633569fb2d86dcc067`.
The files are temporary local artifacts under `/private/tmp`; their hashes
detect accidental substitution but do not make the evidence durable.

### Targeted paired diagnostics

The following minima compare the specialized candidate with the retained
generic/vector-producing path on the same worktree, except for the additive
compound rows, whose reference is a clean archive of `1d611e4`. Diagonal rows
used 100 samples, Schmidt rows used 300, and additive rows used 1,000 (`n=16`)
or 300 (`n=32`), always with one evaluation and the same one-thread runtime.

| Operation and input | Reference minimum | Candidate minimum | Reference → candidate memory | Allocations |
|---|---:|---:|---:|---:|
| diagonal purity, `n=512`, `Float64` | 17,935,750 ns | 375 ns | 14,890,384 → 16 bytes | 49 → 1 |
| diagonal fidelity, `n=512`, `Float64` | 82,447,709 ns | 1,166 ns | 42,687,872 → 16 bytes | 125 → 1 |
| diagonal trace distance, `n=512`, `Float64` | 67,037,958 ns | 958 ns | 34,290,624 → 16 bytes | 107 → 1 |
| diagonal `l1` coherence, `n=512`, `Float64` | 17,871,875 ns | 292 ns | 14,894,544 → 16 bytes | 52 → 1 |
| diagonal relative-entropy coherence, `n=512`, `Float64` | 17,965,667 ns | 333 ns | 14,898,736 → 16 bytes | 55 → 1 |
| Schmidt coefficients, `32 × 32`, `ComplexF64` | 110,458 ns | 41,000 ns | 167,552 → 87,536 bytes | 46 → 40 |
| dense order-two additive compound, `n=16`, `Float64` | 69,917 ns | 3,666 ns | 429,016 → 132,488 bytes | 7,297 → 68 |
| dense order-two additive compound, `n=32`, `Float64` | 633,334 ns | 30,333 ns | 4,476,296 → 1,983,880 bytes | 61,607 → 68 |

The diagonal formulas agreed with the retained dense methods in 800 seeded
Float32/Float64 comparisons; the largest absolute difference was
`6.88e-7`, from Float32 relative-entropy cancellation. Values-only Schmidt
results were compared with the full decomposition over real and complex,
square and rectangular inputs. Additive-compound results were bitwise equal to
the retained generic construction for floating inputs and exact for integer and
rational inputs in the focused suite. Sparse and higher-order additive paths
are unchanged.

These quick-run and paired minima are diagnostic observations, not stable
speedup guarantees. The additive rows were rerun after final review with seed
`0x4144445045524632`, 1,000/300 samples for dimensions 16/32, one evaluation,
and a three-second limit in both this tree and an extracted archive of
`1d611e4`. The other targeted commands and all raw trial objects were not
retained, so those rows remain ephemeral observations rather than durable
benchmark evidence. There is no repeated-session noise study, acceptance
threshold, cross-platform comparison, or regression gate.

<!-- qetlab-current-claims: end -->

## Previous 114-case paired smoke (2026-07-30)

The complete set of 114 declared quick benchmark cases completed with exit
status zero on Julia 1.12.6 for both the clean baseline
`f32dd233e478dd6e2642f11fab088f6c8febc420` and the candidate worktree. Both
runs fixed Julia and BLAS to one thread:

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --compiled-modules=no --startup-file=no --project=benchmark \
  benchmark/benchmarks.jl --quick \
  --output=/private/tmp/qet-quick-f32dd233-julia1.12.6-t1

JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --compiled-modules=no --startup-file=no --project=benchmark \
  benchmark/benchmarks.jl --quick \
  --output=/private/tmp/qet-quick-candidate-final-julia1.12.6-t1
```

Each case produced 20 samples with one evaluation. The environment was Apple
arm64 (`apple-m4`), Julia 1.12.6, BenchmarkTools 1.8.0, and ILP64 OpenBLAS
through libblastrampoline. The candidate differed from the clean baseline only
in the two optimized source files and their two regression-test files when the
benchmark metadata was captured.

The targeted observations were:

| Case | Baseline minimum | Candidate minimum | Baseline → candidate memory | Allocations |
|---|---:|---:|---:|---:|
| `channels/kraus_to_choi_dephasing_d16` | 789,125 ns | 106,167 ns (7.4× lower) | 16,779,968 → 1,081,584 bytes | 98 → 9 |
| `channels/operator_sum_to_choi_16x12_to_20x10_terms4` | 108,292 ns | 84,208 ns (1.29× lower) | 4,981,728 → 1,274,432 bytes | 34 → 20 |
| `matrix_analysis/compound_dense_8_order3` | 2,947,542 ns | 21,792 ns (135× lower) | 1,738,768 → 32,784 bytes | 53,367 → 55 |

Independent same-size checks measured maximum absolute differences of
`1.11e-16` for the Kraus Choi result, `1.84e-15` for the paired operator-sum
Choi result, and `3.56e-15` for the order-three compound matrix against the
previous termwise or standard-library constructions. Exact rational results
remained exact, `BigFloat` retained the standard-library determinant path, an
ill-conditioned three-by-three probe agreed exactly with `det`, and sparse
outputs remained sparse.

The raw JSON/TOML SHA-256 pairs are
`8f74f49a36fba76a9ba6e53519ac88c8b36cb17c2fe80c47706f0c880080a1d8` /
`4fc49c1cf77adfbd56311120adeabbaed0c3499e89f5e5dfa34b28f57955d85a`
for the baseline and
`54098ff810148f5fa9acd84580c61a66b1c3c98f5b9dfb88e4f408948723c9bd` /
`b7e2f298f90c5221b688b41a285b699aa66705f7af236b014c50bf5fd8b25db3`
for the candidate. The files are local temporary artifacts under
`/private/tmp`; the hashes make accidental substitution detectable but do not
make them durable repository evidence.

These quick-run minima are diagnostic observations, not stable speedup
guarantees. There is no repeated-session noise study, acceptance threshold,
cross-platform comparison, or regression gate.

## Historical 42-case quick smoke

The organized 42-case suite completed across subsystem kernels, states,
operators, random objects, channels, scalar measures, entanglement criteria,
coherence, product analysis, and matrix analysis. These single-worktree minima
are diagnostic evidence only: there is no committed baseline,
repeated-environment noise study, MATLAB comparison, allocation budget, or
sparse-scaling study. Matrix predicates do not yet have benchmark cases.

## Local quick-smoke evidence

Command:

```sh
julia --startup-file=no --project=benchmark benchmark/benchmarks.jl --quick --no-save
```

Environment: commit
`9b0d0d3b8177eded5b8743d250044d5428625b1b`; clean tracked worktree; Julia
1.12.6; Apple M4 Pro; one Julia thread; ten BLAS threads; ILP64 OpenBLAS
through libblastrampoline; BenchmarkTools 1.8.0. Each quick case used 20
samples and one evaluation.

| Case | Minimum time (ns) | Memory (bytes) | Allocations |
|---|---:|---:|---:|
| `channels/apply_dephasing_d16` | 4,000 | 39,616 | 54 |
| `channels/kraus_to_choi_dephasing_d16` | 273,833 | 16,780,480 | 114 |
| `channels/kraus_to_superoperator_dephasing_d16` | 288,875 | 16,779,968 | 98 |
| `coherence/l1_pure_d4096` | 7,875 | 0 | 0 |
| `coherence/rank_basis_transform_d64` | 76,458 | 264,160 | 16 |
| `coherence/relative_entropy_pure_d4096` | 15,291 | 32,832 | 3 |
| `criteria/ppt_8x8` | 808,834 | 931,408 | 159 |
| `criteria/realignment_8x8` | 604,292 | 756,736 | 156 |
| `criteria/reduction_8x8` | 1,152,833 | 1,467,920 | 226 |
| `entanglement/logarithmic_negativity_8x8` | 572,542 | 864,592 | 153 |
| `entanglement/negativity_8x8` | 573,250 | 864,592 | 153 |
| `entanglement/schmidt_decomposition_32x32` | 91,667 | 167,520 | 45 |
| `entanglement/schmidt_rank_32x32` | 91,875 | 167,520 | 45 |
| `matrix_analysis/additive_compound_dense_16_order2` | 63,500 | 427,872 | 7,242 |
| `matrix_analysis/compound_dense_8_order3` | 3,270,208 | 1,738,768 | 53,367 |
| `matrix_analysis/majorizes_singular_values_32x32` | 56,333 | 72,720 | 187 |
| `measures/concurrence_mixed_two_qubit` | 2,709 | 15,696 | 69 |
| `measures/entropy_density_d64` | 376,125 | 527,872 | 37 |
| `measures/fidelity_density_d32` | 198,916 | 444,512 | 108 |
| `measures/ky_fan_k16_64x48` | 119,042 | 170,928 | 19 |
| `measures/purity_density_d64` | 420,417 | 527,872 | 37 |
| `measures/schatten_p3_64x48` | 118,667 | 170,928 | 19 |
| `measures/trace_distance_density_d64` | 988,875 | 1,277,456 | 93 |
| `measures/trace_norm_64x48` | 118,375 | 170,928 | 19 |
| `operators/fourier_d32` | 5,250 | 16,464 | 3 |
| `partial_trace/dense_matrix_plan_reuse` | 208 | 1,136 | 2 |
| `partial_trace/dense_pure_plan_reuse` | 166 | 2,000 | 4 |
| `partial_trace/plan_construction` | 958 | 5,184 | 99 |
| `partial_trace/sparse_matrix_plan_reuse` | 1,375 | 22,256 | 43 |
| `partial_transpose/dense_plan_reuse` | 1,875 | 49,232 | 3 |
| `partial_transpose/sparse_plan_reuse` | 3,583 | 63,600 | 52 |
| `permutation/dense_matrix_plan_reuse` | 1,333 | 49,232 | 3 |
| `product/formation_mixed_two_qubit` | 3,083 | 16,528 | 92 |
| `product/operator_analysis_4x8x2` | 163,375 | 1,010,064 | 727 |
| `product/operator_schmidt_decomposition_8x8` | 594,917 | 1,169,856 | 941 |
| `projector/symmetric_d4_p4_sparse` | 68,625 | 335,840 | 1,627 |
| `random/density_d32_rank16` | 5,333 | 49,552 | 15 |
| `random/povm_d8_outcomes4` | 6,917 | 33,520 | 42 |
| `random/unitary_d32` | 44,417 | 149,392 | 35 |
| `realignment/dense_plan_reuse` | 5,250 | 65,616 | 3 |
| `states/werner_d8_sparse` | 37,708 | 29,584 | 348 |
| `tensor_product/dense_two_factor` | 1,167 | 50,368 | 5 |

Do not use this table to claim a speedup or stable regression threshold.
`--no-save` intentionally produced no raw artifact.

## Required result schema

Future entries must record:

- exact repository commit and dependency manifest;
- hardware, operating system, Julia/MATLAB, and BLAS;
- Julia, BLAS, and solver threads;
- operation, algorithm, problem dimensions, numeric type, and sparsity;
- setup/warmup exclusion, samples, estimator, time, allocations, and peak memory
  where measurable;
- numerical error against the validation reference;
- raw machine-readable result location.

## Covered smoke categories

- subsystem permutation;
- partial trace of vectors and matrices;
- partial transpose;
- realignment/reshuffling;
- repeated operations with and without reusable plans;
- dense, structured, sparse, and representative generic numeric inputs.

These categories have smoke-suite cases, but not all have scaling studies or
comparative baselines. Performance gates will be set only after correctness
baselines and CI-noise measurements exist.
