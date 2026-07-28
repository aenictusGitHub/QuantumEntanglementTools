# Sparse and large-scale computations

Sparse preservation is part of the API contract, not an afterthought.

## Policy

- Never silently call `Matrix(x)` on a sparse or structured input.
- Prefer index plans, contractions, sparse accumulation, or matrix-free actions
  over explicit enormous permutation or superoperator matrices.
- If densification is mathematically unavoidable, make it explicit in the API,
  estimate output size first, guard it, and document peak-memory implications.
- Reusable plans and workspaces should be caller-owned unless a bounded,
  thread-safe cache is justified by measurements.
- Mutating methods require clear aliasing rules and an actual allocation benefit.

## Evidence required

An operation advertised as sparse-aware needs tests for correctness, output
structure where promised, unusual sparse index patterns, dimension failures, and
allocation behavior. Large-scale guidance must name dimensions, sparsity,
numeric type, Julia/BLAS threads, and the exact algorithm; broad memory or speed
claims are not acceptable.

GPU and tensor-network support are not implied by accepting `AbstractArray`.
Device-specific work belongs in optional extensions with scalar-indexing and
transfer tests.
