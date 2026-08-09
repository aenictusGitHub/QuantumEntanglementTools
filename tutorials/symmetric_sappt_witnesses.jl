# SPDX-FileCopyrightText: 2026 John Martin
# SPDX-License-Identifier: BSD-3-Clause
#
# Independently written executable tutorial citing selected analytic results
# and rounded witness coefficients from J. Louvet, E. Serrano-Ensástiga,
# T. Bastin, and J. Martin, "Nonequivalence between absolute separability and
# positive partial transposition in the symmetric subspace," Phys. Rev. A 111,
# 042418 (2025), https://doi.org/10.1103/PhysRevA.111.042418.
# Inspected author-supplied publisher PDF SHA-256:
# c9a32248a9f0f730409dc690db0dfc8ab97063bd600041cb7e567d8b0e9aaf7a.
# The paper PDF, prose, figures, and source SDP are not copied or redistributed;
# this tutorial's original Julia implementation is covered by the repository's
# BSD 3-Clause License.

module TutorialSymmetricSAPPTWitnesses

using LinearAlgebra
using QuantumEntanglementTools

const PUBLISHED_WITNESS_COEFFICIENTS = Dict(
    5 => (diagonal=[0.0366656, -0.134595, 1.0, 1.0, -0.134595, 0.0366656], corner=-9.31947),
    7 => (
        diagonal=[
            0.00197514, 0.0643064, -0.189017, 1.0, 1.0, -0.189017, 0.0643064, 0.00197514
        ],
        corner=-31.2405,
    ),
    9 => (
        diagonal=[
            0.00235791,
            -0.013747,
            0.0621661,
            -0.1636915,
            1.0,
            1.0,
            -0.1636915,
            0.0621661,
            -0.013747,
            0.00235791,
        ],
        corner=-114.305,
    ),
)

"""
    sappt_threshold(qubits)

Return the exact spectral threshold above which the paper's state family is
symmetric absolutely positive under partial transpose (SAPPT).
"""
function sappt_threshold(qubits::Integer)
    qubits >= 2 || throw(ArgumentError("qubits must be at least two"))
    qubit_count = BigInt(qubits)
    scale = (qubit_count + 1) * binomial(qubit_count, qubits ÷ 2)
    if scale <= typemax(Int) - 2
        machine_scale = Int(scale)
        return machine_scale // (machine_scale + 2)
    end
    return scale // (scale + 2)
end

function _ghz_dicke_vector(qubits::Integer; phase::Real=0)
    state = zeros(ComplexF64, qubits + 1)
    state[1] = inv(sqrt(2))
    state[end] = cis(phase) / sqrt(2)
    return state
end

function symmetric_family_state(qubits::Integer, p::Real, pure_state::AbstractVector)
    length(pure_state) == qubits + 1 ||
        throw(DimensionMismatch("pure_state must have $(qubits + 1) Dicke coefficients"))
    0 <= p <= 1 || throw(ArgumentError("p must lie in [0, 1]"))
    state_norm = norm(pure_state)
    isfinite(state_norm) || throw(ArgumentError("pure_state must have a finite norm"))
    isapprox(state_norm, one(state_norm); atol=1e-12, rtol=1e-12) ||
        throw(ArgumentError("pure_state must be normalized; received norm $state_norm"))
    scalar_type = promote_type(eltype(pure_state), typeof(p))
    identity_state = Matrix{scalar_type}(I, qubits + 1, qubits + 1) / (qubits + 1)
    return p * identity_state + (one(p) - p) * (pure_state * pure_state')
end

function published_symmetric_witness(qubits::Integer; phase::Real=0)
    haskey(PUBLISHED_WITNESS_COEFFICIENTS, qubits) || throw(
        ArgumentError("published coefficients are available only for 5, 7, and 9 qubits"),
    )
    coefficients = PUBLISHED_WITNESS_COEFFICIENTS[qubits]
    witness = Matrix(Diagonal(ComplexF64.(coefficients.diagonal)))

    # The displayed negative corner is phase-matched to GHZ+. A collective
    # local phase rotation preserves block positivity and treats any other GHZ
    # phase without silently changing conventions.
    witness[1, end] = coefficients.corner * cis(-phase)
    witness[end, 1] = conj(witness[1, end])
    return witness
end

function _coherent_dicke_vector(qubits::Integer, z::Real, phase::Real)
    amplitude_zero = sqrt((1 + z) / 2)
    amplitude_one = sqrt((1 - z) / 2) * cis(phase)
    return ComplexF64[
        sqrt(binomial(qubits, excitation)) *
        amplitude_zero^(qubits - excitation) *
        amplitude_one^excitation for excitation in 0:qubits
    ]
end

function _identity_decomposition_n5()
    # Three-point Gauss--Legendre quadrature integrates every diagonal
    # degree-five polynomial exactly. Six equally spaced phases cancel all
    # off-diagonal Dicke-basis terms.
    z_nodes = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
    z_weights = (5 / 18, 4 / 9, 5 / 18)
    phases = ntuple(index -> 2pi * (index - 1) / 6, 6)

    states = Vector{Vector{ComplexF64}}()
    weights = Float64[]
    for (z, z_weight) in zip(z_nodes, z_weights)
        for phase in phases
            push!(states, _coherent_dicke_vector(5, z, phase))
            push!(weights, z_weight / 6)
        end
    end
    return states, weights
end

function _separable_state_and_decomposition_n5(p::Real)
    states, weights = _identity_decomposition_n5()
    decomposition_states = copy(states)
    decomposition_weights = Float64(p) .* weights

    product_zero = zeros(ComplexF64, 6)
    product_zero[1] = 1
    push!(decomposition_states, product_zero)
    push!(decomposition_weights, 1 - Float64(p))

    state = sum(
        weight * (vector * vector') for
        (weight, vector) in zip(decomposition_weights, decomposition_states)
    )
    return state, decomposition_states, decomposition_weights
end

function _restricted_partial_transpose_minimum(
    state_dicke::AbstractMatrix, qubits::Integer, first_group::Integer
)
    basis = symmetric_subspace_basis(2, qubits; sparse_output=false)
    state_full = basis * state_dicke * basis'
    split_basis = tensor_product(
        symmetric_subspace_basis(2, first_group; sparse_output=false),
        symmetric_subspace_basis(2, qubits - first_group; sparse_output=false),
    )
    state_split = split_basis' * state_full * split_basis
    transposed = partial_transpose(
        state_split, (first_group + 1, qubits - first_group + 1); systems=(1,)
    )
    return minimum(eigvals(Hermitian(Matrix(transposed))))
end

function _published_witness_summary(qubits::Integer)
    threshold = sappt_threshold(qubits)
    ghz = _ghz_dicke_vector(qubits)
    state = symmetric_family_state(qubits, threshold, ghz)
    witness = published_symmetric_witness(qubits)
    expectation = real(tr(witness * state))

    mixed_expectation = real(tr(witness)) / (qubits + 1)
    ghz_expectation = real(dot(ghz, witness * ghz))
    detection_limit = -ghz_expectation / (mixed_expectation - ghz_expectation)
    return (
        qubits=qubits,
        threshold=threshold,
        expectation=expectation,
        detection_limit=detection_limit,
    )
end

function _power_to_bernstein(power_coefficients::AbstractVector{<:Rational})
    degree = length(power_coefficients) - 1
    return [
        sum(
            power_coefficients[power + 1] *
            (binomial(index, power) // binomial(degree, power)) for power in 0:index
        ) for index in 0:degree
    ]
end

function _elevate_bernstein(coefficients::AbstractVector{<:Rational})
    degree = length(coefficients) - 1
    elevated = similar(coefficients, length(coefficients) + 1)
    elevated[1] = first(coefficients)
    for index in 1:degree
        fraction = index // (degree + 1)
        elevated[index + 1] =
            fraction * coefficients[index] +
            (one(fraction) - fraction) * coefficients[index + 1]
    end
    elevated[end] = last(coefficients)
    return elevated
end

function _witness_product_minimum_n5()
    # Treat each published decimal as the exact rational it prints. With
    # u=2t, factor f(u/2)-f(1/2)=(1-u)q(u). Exact positivity of q on [0,1]
    # follows when a degree-elevated Bernstein representation has only positive
    # coefficients.
    a = 366656 // 10_000_000
    b = -134595 // 1_000_000
    corner = -931947 // 100_000
    exact_minimum = (a + 5b + 10 + corner) / 16

    difference_power = Rational{Int}[
        a - exact_minimum, 0, 5(b - a) / 4, 0, (5a - 15b + 10) / 16, corner / 16
    ]
    quotient_power = similar(difference_power, 5)
    quotient_power[1] = difference_power[1]
    for index in 2:5
        quotient_power[index] = difference_power[index] + quotient_power[index - 1]
    end
    @assert difference_power[6] == -quotient_power[5]

    bernstein = _power_to_bernstein(quotient_power)
    for _ in 5:8
        bernstein = _elevate_bernstein(bernstein)
    end
    @assert all(>(0), bernstein)

    return (
        minimum=Float64(exact_minimum),
        exact_minimum=exact_minimum,
        t=0.5,
        bernstein_coefficients=bernstein,
    )
end

function _npt_witness_example(p::Real)
    qubits = 5
    ghz = _ghz_dicke_vector(qubits)
    state_dicke = symmetric_family_state(qubits, p, ghz)
    basis = symmetric_subspace_basis(2, qubits; sparse_output=false)
    state_full = basis * state_dicke * basis'

    # This is the negative-eigenvalue vector of rho^(T_A) for the 2|3 split.
    eta =
        (
            tensor_product(dicke_state(2, 0), dicke_state(3, 3)) -
            tensor_product(dicke_state(2, 2), dicke_state(3, 0))
        ) / sqrt(2)
    witness = partial_transpose(eta * eta', ntuple(_ -> 2, qubits); systems=(1, 2))
    expectation = real(tr(witness * state_full))
    minimum_eigenvalue = _restricted_partial_transpose_minimum(state_dicke, qubits, 2)
    return (expectation=expectation, minimum_eigenvalue=minimum_eigenvalue, witness=witness)
end

"""
    run(; io=stdout)

Reproduce the five-qubit separability comparison and witness calculations from
Louvet et al., Phys. Rev. A 111, 042418 (2025). The rounded published
coefficients are reconstructed; their symmetric-extension SDP derivation is
not rerun.
"""
function run(; io::IO=stdout)
    qubits = 5
    threshold = sappt_threshold(qubits)
    product_zero = zeros(ComplexF64, qubits + 1)
    product_zero[1] = 1
    ghz_plus = _ghz_dicke_vector(qubits)

    separable_state, decomposition_states, decomposition_weights = _separable_state_and_decomposition_n5(
        threshold
    )
    separable_state_direct = symmetric_family_state(qubits, threshold, product_zero)
    ghz_state_at_threshold = symmetric_family_state(qubits, threshold, ghz_plus)

    decomposition_error = norm(separable_state - separable_state_direct)
    normalization_error = maximum(abs(norm(state) - 1) for state in decomposition_states)
    same_spectrum_error = norm(
        eigvals(Hermitian(separable_state_direct)) -
        eigvals(Hermitian(ghz_state_at_threshold)),
    )

    witness = published_symmetric_witness(qubits)
    witness_at_threshold = real(tr(witness * ghz_state_at_threshold))
    witness_summary = _published_witness_summary(qubits)
    product_minimum = _witness_product_minimum_n5()

    # This rational point lies strictly inside the SAPPT interval detected by
    # the rounded W5 witness, rather than on either numerical boundary.
    p_demo = 121 // 125
    ghz_demo = symmetric_family_state(qubits, p_demo, ghz_plus)
    restricted_pt_minimum_1_4 = _restricted_partial_transpose_minimum(ghz_demo, qubits, 1)
    restricted_pt_minimum_2_3 = _restricted_partial_transpose_minimum(ghz_demo, qubits, 2)
    witness_at_demo = real(tr(witness * ghz_demo))

    p_npt = 29 // 30
    npt_example = _npt_witness_example(p_npt)

    # The journal PDF writes GHZ- while the displayed negative witness corner
    # matches GHZ+. Rotating the witness corner makes the minus convention
    # explicit and yields the same physical expectation.
    minus_phase = pi
    ghz_minus = _ghz_dicke_vector(qubits; phase=minus_phase)
    minus_state = symmetric_family_state(qubits, threshold, ghz_minus)
    minus_witness = published_symmetric_witness(qubits; phase=minus_phase)
    minus_expectation = real(tr(minus_witness * minus_state))
    unmatched_minus_expectation = real(tr(witness * minus_state))

    published = Tuple(_published_witness_summary(n) for n in (5, 7, 9))

    @assert threshold == 30 // 31
    @assert isapprox(sum(decomposition_weights), 1; atol=1e-14, rtol=0)
    @assert normalization_error <= 1e-14
    @assert decomposition_error <= 1e-14
    @assert same_spectrum_error <= 1e-14
    @assert threshold < p_demo < witness_summary.detection_limit
    @assert min(restricted_pt_minimum_1_4, restricted_pt_minimum_2_3) > 0
    @assert witness_at_demo < 0
    @assert isapprox(
        witness_summary.detection_limit, 0.9686241592915386; atol=1e-14, rtol=0
    )
    @assert isapprox(product_minimum.minimum, 0.0027637875; atol=1e-13, rtol=0)
    @assert product_minimum.exact_minimum == 221103 // 80_000_000
    @assert npt_example.minimum_eigenvalue < 0
    @assert isapprox(
        npt_example.expectation, npt_example.minimum_eigenvalue; atol=1e-14, rtol=0
    )
    @assert unmatched_minus_expectation > 0
    @assert isapprox(minus_expectation, witness_at_threshold; atol=1e-13, rtol=0)
    @assert isapprox(published[2].detection_limit, 0.9930282521603; atol=1e-12, rtol=0)
    @assert isapprox(published[3].detection_limit, 0.9984502357594; atol=1e-12, rtol=0)

    println(io, "Five-qubit family rho(p) in the six-dimensional Dicke basis")
    println(io, "  SAPPT threshold: p_min = ", threshold, " = ", Float64(threshold))
    println(
        io,
        "  Explicit separable representative: ",
        length(decomposition_weights),
        " product terms; reconstruction error = ",
        decomposition_error,
    )
    println(io, "  Same spectrum as GHZ representative: error = ", same_spectrum_error)
    println(
        io,
        "  Strict SAPPT example at p = ",
        p_demo,
        ": restricted-sector PT minima (1|4, 2|3) = ",
        (restricted_pt_minimum_1_4, restricted_pt_minimum_2_3),
    )
    println(
        io,
        "  Published W5 expectation = ",
        witness_at_demo,
        " < 0; rounded-coefficient endpoint = ",
        witness_summary.detection_limit,
    )
    println(
        io,
        "  Exact minimum of the rounded W5 over symmetric products = ",
        product_minimum.exact_minimum,
        " at cos(theta/2)sin(theta/2) = ",
        product_minimum.t,
    )
    println(
        io,
        "  Below p_min at p = ",
        p_npt,
        ": PT minimum = witness expectation = ",
        npt_example.expectation,
    )
    println(
        io,
        "  Printed GHZ- with unmatched corner gives ",
        unmatched_minus_expectation,
        "; phase-matched expectation = ",
        minus_expectation,
    )
    println(io, "  Published rounded witnesses:")
    for entry in published
        println(
            io,
            "    N=",
            entry.qubits,
            ": Tr(W rho(p_min))=",
            entry.expectation,
            ", endpoint=",
            entry.detection_limit,
        )
    end

    return (
        qubits=qubits,
        p_sappt_min=threshold,
        p_demo=p_demo,
        same_spectrum_error=same_spectrum_error,
        separable_decomposition_terms=length(decomposition_weights),
        separable_decomposition_weight_sum=sum(decomposition_weights),
        separable_decomposition_error=decomposition_error,
        restricted_pt_minimum_1_4=restricted_pt_minimum_1_4,
        restricted_pt_minimum_2_3=restricted_pt_minimum_2_3,
        witness_expectation_at_pmin=witness_at_threshold,
        witness_expectation_demo=witness_at_demo,
        witness_detection_limit=witness_summary.detection_limit,
        product_witness_minimum=product_minimum.minimum,
        npt_parameter=p_npt,
        npt_minimum=npt_example.minimum_eigenvalue,
        npt_witness_expectation=npt_example.expectation,
        unmatched_minus_witness_expectation=unmatched_minus_expectation,
        minus_phase_witness_expectation=minus_expectation,
        published=published,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run()
end

end
