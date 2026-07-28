# Benchmark report

Evidence date: 2026-07-28.

Status: local quick smoke only; no comparative performance claim.

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

Environment: uncommitted pre-alpha worktree; Julia 1.12.6; Apple M4 Pro; one
Julia thread; ten BLAS threads; ILP64 OpenBLAS through libblastrampoline;
BenchmarkTools 1.8.0. Each quick case used 20 samples and one evaluation.

| Case | Minimum time (ns) | Memory (bytes) | Allocations |
|---|---:|---:|---:|
| `channels/apply_dephasing_d16` | 3,667 | 39,616 | 54 |
| `channels/kraus_to_choi_dephasing_d16` | 293,167 | 16,780,480 | 114 |
| `channels/kraus_to_superoperator_dephasing_d16` | 342,584 | 16,779,968 | 98 |
| `coherence/l1_pure_d4096` | 7,458 | 0 | 0 |
| `coherence/rank_basis_transform_d64` | 80,958 | 264,160 | 16 |
| `coherence/relative_entropy_pure_d4096` | 17,708 | 32,832 | 3 |
| `criteria/ppt_8x8` | 806,750 | 931,408 | 159 |
| `criteria/realignment_8x8` | 577,625 | 756,736 | 156 |
| `criteria/reduction_8x8` | 1,157,750 | 1,467,920 | 226 |
| `entanglement/logarithmic_negativity_8x8` | 563,458 | 864,592 | 153 |
| `entanglement/negativity_8x8` | 573,292 | 864,592 | 153 |
| `entanglement/schmidt_decomposition_32x32` | 93,875 | 167,520 | 45 |
| `entanglement/schmidt_rank_32x32` | 91,875 | 167,520 | 45 |
| `matrix_analysis/additive_compound_dense_16_order2` | 67,208 | 427,872 | 7,242 |
| `matrix_analysis/compound_dense_8_order3` | 3,203,791 | 1,738,768 | 53,367 |
| `matrix_analysis/majorizes_singular_values_32x32` | 55,000 | 72,720 | 187 |
| `measures/concurrence_mixed_two_qubit` | 2,625 | 15,696 | 69 |
| `measures/entropy_density_d64` | 383,208 | 527,872 | 37 |
| `measures/fidelity_density_d32` | 208,875 | 444,512 | 108 |
| `measures/ky_fan_k16_64x48` | 121,334 | 170,928 | 19 |
| `measures/purity_density_d64` | 421,042 | 527,872 | 37 |
| `measures/schatten_p3_64x48` | 123,750 | 170,928 | 19 |
| `measures/trace_distance_density_d64` | 986,625 | 1,277,456 | 93 |
| `measures/trace_norm_64x48` | 119,042 | 170,928 | 19 |
| `operators/fourier_d32` | 5,167 | 16,464 | 3 |
| `partial_trace/dense_matrix_plan_reuse` | 208 | 1,136 | 2 |
| `partial_trace/dense_pure_plan_reuse` | 167 | 2,000 | 4 |
| `partial_trace/plan_construction` | 916 | 5,184 | 99 |
| `partial_trace/sparse_matrix_plan_reuse` | 1,250 | 22,256 | 43 |
| `partial_transpose/dense_plan_reuse` | 1,916 | 49,232 | 3 |
| `partial_transpose/sparse_plan_reuse` | 3,750 | 63,600 | 52 |
| `permutation/dense_matrix_plan_reuse` | 1,250 | 49,232 | 3 |
| `product/formation_mixed_two_qubit` | 3,125 | 16,528 | 92 |
| `product/operator_analysis_4x8x2` | 164,583 | 1,010,064 | 727 |
| `product/operator_schmidt_decomposition_8x8` | 609,417 | 1,169,856 | 941 |
| `projector/symmetric_d4_p4_sparse` | 64,792 | 335,840 | 1,627 |
| `random/density_d32_rank16` | 5,458 | 49,552 | 15 |
| `random/povm_d8_outcomes4` | 6,958 | 33,520 | 42 |
| `random/unitary_d32` | 44,208 | 149,392 | 35 |
| `realignment/dense_plan_reuse` | 5,291 | 65,616 | 3 |
| `states/werner_d8_sparse` | 37,125 | 29,584 | 348 |
| `tensor_product/dense_two_factor` | 1,167 | 50,368 | 5 |

Do not use this table to claim a speedup or stable regression threshold. The
worktree had no commit identity and `--no-save` intentionally produced no raw
artifact.

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
