# Independently written Julia rendition of the pure-state workflow on QETLAB's
# "Getting started" homepage: https://qetlab.com/
#
# The package APIs exercised here are independent Julia reimplementations of
# QETLAB RandomStateVector.m, SchmidtDecomposition.m, Tensor.m, and TensorSum.m
# at revision d8589610f00cff106537268dee2e2a1153f3a601. QETLAB is copyright
# Nathaniel Johnston and distributed under BSD-2-Clause; see
# licenses/QETLAB-LICENSE.txt and PROVENANCE.toml.

module TutorialQETLABIntroSchmidt

using LinearAlgebra
using QuantumEntanglementTools
using Random: Xoshiro

"""
    run(; io=stdout)

Generate a reproducible random pure state on `3 × 3`, compute its Schmidt
decomposition, and reconstruct the state both term by term and with
[`tensor_sum`](@ref).
"""
function run(; io::IO=stdout)
    seed = UInt64(0x5145544c41424a55)
    dimensions = (3, 3)
    state = random_state_vector(Xoshiro(seed), dimensions)
    decomposition = schmidt_decomposition(state, dimensions)

    coefficients = decomposition.coefficients
    left_vectors = decomposition.left_vectors
    right_vectors = decomposition.right_vectors

    manual_reconstruction = sum(
        coefficients[k] * tensor_product(left_vectors[:, k], right_vectors[:, k]) for
        k in eachindex(coefficients)
    )
    tensor_sum_reconstruction = tensor_sum(
        left_vectors, right_vectors; weights=coefficients
    )

    normalization_error = abs(norm(state) - 1)
    coefficient_normalization_error = abs(sum(abs2, coefficients) - 1)
    left_orthogonality_error = norm(adjoint(left_vectors) * left_vectors - I, Inf)
    right_orthogonality_error = norm(adjoint(right_vectors) * right_vectors - I, Inf)
    manual_reconstruction_error = norm(manual_reconstruction - state)
    tensor_sum_reconstruction_error = norm(tensor_sum_reconstruction - state)
    reconstruction_agreement_error = norm(manual_reconstruction - tensor_sum_reconstruction)
    numerical_schmidt_rank = schmidt_rank(state, dimensions; atol=1e-12, rtol=0)
    replay_exact = random_state_vector(Xoshiro(seed), dimensions) == state
    coefficients_nonnegative = all(coefficient -> coefficient >= 0, coefficients)
    coefficients_descending = issorted(coefficients; rev=true)

    @assert length(state) == prod(dimensions)
    @assert size(left_vectors) == (3, 3)
    @assert size(right_vectors) == (3, 3)
    @assert length(coefficients) == 3
    @assert coefficients_nonnegative
    @assert coefficients_descending
    @assert numerical_schmidt_rank == 3
    @assert normalization_error <= 1e-12
    @assert coefficient_normalization_error <= 1e-12
    @assert left_orthogonality_error <= 1e-12
    @assert right_orthogonality_error <= 1e-12
    @assert manual_reconstruction_error <= 1e-12
    @assert tensor_sum_reconstruction_error <= 1e-12
    @assert reconstruction_agreement_error <= 1e-12
    @assert replay_exact

    println(io, "Seeded random state on 3×3")
    println(io, "Schmidt coefficients: ", coefficients)
    println(io, "Numerical Schmidt rank: ", numerical_schmidt_rank)
    println(io, "Manual reconstruction error: ", manual_reconstruction_error)
    println(io, "tensor_sum reconstruction error: ", tensor_sum_reconstruction_error)

    return (
        seed=seed,
        local_dimensions=dimensions,
        state_dimension=length(state),
        term_count=length(coefficients),
        numerical_schmidt_rank=numerical_schmidt_rank,
        coefficients_nonnegative=coefficients_nonnegative,
        coefficients_descending=coefficients_descending,
        normalization_error=normalization_error,
        coefficient_normalization_error=coefficient_normalization_error,
        left_orthogonality_error=left_orthogonality_error,
        right_orthogonality_error=right_orthogonality_error,
        manual_reconstruction_error=manual_reconstruction_error,
        tensor_sum_reconstruction_error=tensor_sum_reconstruction_error,
        reconstruction_agreement_error=reconstruction_agreement_error,
        replay_exact=replay_exact,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
