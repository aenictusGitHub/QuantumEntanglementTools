module TutorialSubsystemReductions

using LinearAlgebra
using QuantumEntanglementTools

"""
    run(; io=stdout)

Build a three-qubit state whose last two qubits form a Bell pair, then verify
reductions, a subsystem permutation, and a partial-transpose spectrum.
"""
function run(; io::IO=stdout)
    zero_state = [1.0, 0.0]
    bell = bell_state()
    state = tensor_product(zero_state, bell)

    reduced_pair = partial_trace(state, (2, 2, 2); trace_out=(1,))
    reduced_qubit = partial_trace(state, (2, 2, 2); trace_out=(1, 3))

    permuted_state = permute_subsystems(state, (2, 2, 2); permutation=(2, 3, 1))
    pair_after_permutation = partial_trace(permuted_state, (2, 2, 2); trace_out=(3,))

    bell_density = bell * bell'
    maximally_mixed_qubit = Matrix{Float64}(I, 2, 2) / 2
    transposed_pair = partial_transpose(reduced_pair, (2, 2); systems=(2,))
    transposed_spectrum = eigvals(Hermitian(transposed_pair))

    @assert isapprox(reduced_pair, bell_density; atol=1e-12, rtol=0)
    @assert isapprox(reduced_qubit, maximally_mixed_qubit; atol=1e-12, rtol=0)
    @assert isapprox(pair_after_permutation, reduced_pair; atol=1e-12, rtol=0)
    @assert isapprox(tr(reduced_pair), 1; atol=1e-12, rtol=0)
    @assert isapprox(minimum(transposed_spectrum), -0.5; atol=1e-12, rtol=0)

    println(io, "Three-qubit state dimension: ", length(state))
    println(io, "Retained Bell-pair purity: ", purity(reduced_pair))
    println(io, "Single-qubit reduced purity: ", purity(reduced_qubit))
    println(io, "Minimum partial-transpose eigenvalue: ", minimum(transposed_spectrum))

    return (
        state_dimension=length(state),
        reduced_pair_size=size(reduced_pair),
        reduced_pair_trace=real(tr(reduced_pair)),
        reduced_pair_purity=purity(reduced_pair),
        reduced_qubit_purity=purity(reduced_qubit),
        minimum_partial_transpose_eigenvalue=minimum(transposed_spectrum),
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
