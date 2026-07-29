# Source-informed independent Julia implementation based on the specifications
# of QETLAB ComplementaryMap.m, DualMap.m, PartialMap.m,
# DepolarizingChannel.m, DephasingChannel.m, PauliChannel.m, ChoiMap.m, and
# ReductionMap.m at d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

"""
    dual_channel(map)

Return the Hilbert--Schmidt adjoint map, characterized by
`tr(Y' * Φ(X)) == tr(dual_channel(Φ)(Y)' * X)`.

The returned representation has the same representation kind as the input.
For Kraus operators this replaces every `K` by `K'`; for a superoperator it
uses the adjoint transfer matrix.
"""
function dual_channel(map::KrausRepresentation)
    operators = _validated_kraus_operators(map)
    return KrausRepresentation([copy(adjoint(operator)) for operator in operators])
end

function dual_channel(map::SuperoperatorRepresentation)
    return SuperoperatorRepresentation(
        copy(adjoint(_validated_representation_matrix(map))), map.output_dim, map.input_dim
    )
end

function dual_channel(map::ChoiRepresentation)
    return choi_representation(dual_channel(superoperator_representation(map)))
end

"""
    complementary_channel(map; atol, rtol)

Construct a complementary completely positive map from the supplied Kraus
dilation. If `K_r` has shape `(d_out, d_in)`, the complement has environment
dimension equal to the number of Kraus operators and Kraus matrices
`L_a[r,i] = K_r[a,i]` for `a = 1:d_out`.

A Kraus input therefore fixes the dilation exactly. Choi and superoperator
inputs first use the tolerance-controlled canonical eigendecomposition from
`kraus_representation`; their complement is consequently defined only up to an
environment isometry. The output representation kind matches the input.
"""
function complementary_channel(
    map::KrausRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    _checked_tolerances(map, atol, rtol)
    operators = _validated_kraus_operators(map)
    environment_dim = length(operators)
    scalar_type = eltype(map)
    complementary_operators = Vector{Matrix{scalar_type}}(undef, map.output_dim)
    @inbounds for output_index in 1:map.output_dim
        operator = Matrix{scalar_type}(undef, environment_dim, map.input_dim)
        for kraus_index in 1:environment_dim, input_index in 1:map.input_dim
            operator[kraus_index, input_index] = operators[kraus_index][
                output_index, input_index
            ]
        end
        complementary_operators[output_index] = operator
    end
    return KrausRepresentation(complementary_operators)
end

function complementary_channel(
    map::ChoiRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    kraus = kraus_representation(map; atol=atol, rtol=rtol)
    return choi_representation(complementary_channel(kraus))
end

function complementary_channel(
    map::SuperoperatorRepresentation; atol=zero(_default_rtol(map)), rtol=_default_rtol(map)
)
    kraus = kraus_representation(map; atol=atol, rtol=rtol)
    return superoperator_representation(complementary_channel(kraus))
end

function _map_has_sparse_storage(map::KrausRepresentation)
    return all(issparse, _validated_kraus_operators(map))
end

function _map_has_sparse_storage(map::ChoiRepresentation)
    return issparse(_validated_representation_matrix(map))
end

function _map_has_sparse_storage(map::SuperoperatorRepresentation)
    return issparse(_validated_representation_matrix(map))
end

"""
    partial_map(input, map, subsystem, dims)

Apply `map` to one subsystem of a multipartite square operator. `dims` lists
subsystems in tensor-product order (subsystem 1 slowest-varying), and
`dims[subsystem]` must equal `input_dimension(map)`. The output subsystem
dimension is replaced by `output_dimension(map)`.

The implementation permutes the selected subsystem to the final tensor
position, applies the map independently to every operator block, and restores
the original subsystem order. Sparse storage is preserved when both the input
and the representation data are sparse; otherwise the result is dense.
"""
function partial_map(input::AbstractMatrix, map::AbstractMapRepresentation, subsystem, dims)
    layout = _as_layout(dims)
    _validate_matrix_dimension(input, layout)
    _validate_numeric_matrix(input, "partial-map input")
    selected = _normalize_systems(
        subsystem, length(layout); name="subsystem", allow_empty=false
    )
    length(selected) == 1 ||
        throw(ArgumentError("partial_map applies to exactly one subsystem"))
    system = only(selected)
    layout[system] == map.input_dim || throw(
        DimensionMismatch(
            "dims[$system]=$(layout[system]) does not match map input dimension " *
            "$(map.input_dim)",
        ),
    )

    permutation = Tuple(
        vcat([index for index in 1:length(layout) if index != system], [system])
    )
    permuted = permute_subsystems(input, layout; permutation=permutation)
    remainder_dims = Tuple(layout[index] for index in permutation[1:(end - 1)])
    remainder = _checked_product(remainder_dims, "unmapped subsystem dimensions")
    output_total = Base.checked_mul(remainder, map.output_dim)
    output_type = promote_type(eltype(input), eltype(map))
    sparse_output = issparse(input) && _map_has_sparse_storage(map)
    permuted_output = if sparse_output
        spzeros(output_type, output_total, output_total)
    else
        zeros(output_type, output_total, output_total)
    end

    input_dim = map.input_dim
    output_dim = map.output_dim
    application_map = map isa ChoiRepresentation ? superoperator_representation(map) : map
    @inbounds for block_column in 1:remainder, block_row in 1:remainder
        input_rows = ((block_row - 1) * input_dim + 1):(block_row * input_dim)
        input_columns = ((block_column - 1) * input_dim + 1):(block_column * input_dim)
        output_rows = ((block_row - 1) * output_dim + 1):(block_row * output_dim)
        output_columns = ((block_column - 1) * output_dim + 1):(block_column * output_dim)
        permuted_output[output_rows, output_columns] = apply_channel(
            permuted[input_rows, input_columns], application_map
        )
    end

    permuted_output_dims = (remainder_dims..., output_dim)
    return permute_subsystems(
        permuted_output,
        permuted_output_dims;
        permutation=Tuple(invperm(collect(permutation))),
    )
end

function _channel_finite_real_parameter(value, name::AbstractString)
    value isa Real ||
        throw(ArgumentError("$name must be a real number; got $(repr(value))"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    return value
end

function _unnormalized_maximally_entangled(dim::Int, ::Type{T}) where {T}
    total = _squared_dimension(dim, "dim")
    return sparsevec(1:(dim + 1):total, fill(one(T), dim), total)
end

"""
    depolarizing_channel(dim, p=0)

Return the Choi representation of
``Φ(X) = p*X + (1-p)*tr(X)*I/dim``.

The full completely-positive trace-preserving range
`-1/(dim^2-1) <= p <= 1` is accepted for `dim > 1`. In dimension one every
finite `p` gives the identity map. The returned Choi matrix is sparse.
"""
function depolarizing_channel(dim, p=0)
    dimension = _map_dimension(dim, "dim")
    parameter = _channel_finite_real_parameter(p, "p")
    choi_dimension = _squared_dimension(dimension, "dim")
    if dimension > 1
        lower = -one(parameter) / (choi_dimension - 1)
        lower <= parameter <= one(parameter) || throw(
            DomainError(
                parameter, "p must lie in the CPTP range [$lower, 1] for dim=$dimension"
            ),
        )
    end
    scalar_type = promote_type(
        typeof((one(parameter) - parameter) / dimension), typeof(parameter)
    )
    omega = _unnormalized_maximally_entangled(dimension, scalar_type)
    diagonal_weight = (one(parameter) - parameter) / dimension
    matrix =
        spdiagm(0 => fill(convert(scalar_type, diagonal_weight), choi_dimension)) +
        parameter * (omega * adjoint(omega))
    return ChoiRepresentation(matrix, dimension, dimension)
end

"""
    dephasing_channel(dim, p=0)

Return the sparse Choi representation of the channel
``Φ(E_ii)=E_ii`` and ``Φ(E_ij)=p*E_ij`` for `i != j`.

The full completely-positive range `-1/(dim-1) <= p <= 1` is accepted for
`dim > 1`. In dimension one every finite `p` gives the identity map.
"""
function dephasing_channel(dim, p=0)
    dimension = _map_dimension(dim, "dim")
    parameter = _channel_finite_real_parameter(p, "p")
    if dimension > 1
        lower = -one(parameter) / (dimension - 1)
        lower <= parameter <= one(parameter) || throw(
            DomainError(
                parameter, "p must lie in the CP range [$lower, 1] for dim=$dimension"
            ),
        )
    end
    scalar_type = promote_type(typeof(one(parameter) - parameter), typeof(parameter))
    omega = _unnormalized_maximally_entangled(dimension, scalar_type)
    coherent = omega * adjoint(omega)
    matrix =
        (one(parameter) - parameter) * spdiagm(0 => diag(coherent)) + parameter * coherent
    return ChoiRepresentation(matrix, dimension, dimension)
end

function _number_of_qubits(probability_count::Int)
    probability_count >= 4 || throw(
        ArgumentError("a Pauli channel needs 4^q probabilities for a positive integer q"),
    )
    qubits = 0
    count = 1
    while count < probability_count
        count = try
            Base.checked_mul(count, 4)
        catch error
            error isa OverflowError || rethrow()
            break
        end
        qubits += 1
    end
    count == probability_count || throw(
        ArgumentError("probability count $probability_count is not 4^q for an integer q"),
    )
    return qubits
end

function _base_four_labels(index::Int, qubits::Int)
    labels = Vector{Int}(undef, qubits)
    remainder = index
    @inbounds for position in qubits:-1:1
        labels[position] = mod(remainder, 4)
        remainder = div(remainder, 4)
    end
    return labels
end

"""
    pauli_channel(probabilities; atol=0, rtol=sqrt(eps(T)))

Return the sparse Choi representation of a multi-qubit Pauli channel.
`probabilities` must contain exactly `4^q` finite, nonnegative real entries
for a positive integer `q`, in lexicographic base-four order with labels
`I, X, Y, Z`. Their sum must equal one within the explicit tolerances.

The probabilities are neither clipped nor normalized.
"""
function pauli_channel(probabilities::AbstractArray; atol=nothing, rtol=nothing)
    Base.require_one_based_indexing(probabilities)
    probability_vector = vec(probabilities)
    qubits = _number_of_qubits(length(probability_vector))
    all(probability -> probability isa Real, probability_vector) ||
        throw(ArgumentError("probabilities must be real numbers"))
    all(isfinite, probability_vector) ||
        throw(ArgumentError("probabilities must contain only finite values"))
    all(probability -> probability >= 0, probability_vector) ||
        throw(DomainError(probabilities, "probabilities must be nonnegative"))
    probability_type = foldl(
        promote_type, (typeof(probability) for probability in probability_vector)
    )
    default_tolerance = sqrt(eps(_real_float_type(probability_type)))
    checked_atol = atol === nothing ? zero(default_tolerance) : atol
    checked_rtol = rtol === nothing ? default_tolerance : rtol
    checked_atol, checked_rtol = _checked_tolerances(nothing, checked_atol, checked_rtol)
    isapprox(
        sum(probability_vector), one(probability_type); atol=checked_atol, rtol=checked_rtol
    ) || throw(
        DomainError(
            sum(probability_vector),
            "probabilities must sum to one at the requested tolerances",
        ),
    )

    dimension = _checked_power(2, qubits, "Pauli channel dimension")
    choi_dimension = _squared_dimension(dimension, "Pauli channel dimension")
    scalar_type = promote_type(probability_type, Complex{Int})
    matrix = spzeros(scalar_type, choi_dimension, choi_dimension)
    for (zero_based_index, probability) in enumerate(probability_vector)
        iszero(probability) && continue
        labels = _base_four_labels(zero_based_index - 1, qubits)
        operator = pauli(labels; sparse_output=true)
        vector = _column_vector(operator)
        matrix += probability * (vector * adjoint(vector))
    end
    return ChoiRepresentation(matrix, dimension, dimension)
end

"""
    choi_map(a=1, b=1, c=0)

Return the sparse Choi representation of the generalized three-dimensional
Choi map from S. J. Cho, S.-H. Kye, and S. G. Lee (1992), using the parameter
ordering and normalization of QETLAB's `ChoiMap`.

This constructor represents a general linear map: it deliberately does not
claim complete positivity or trace preservation.
"""
function choi_map(a=1, b=1, c=0)
    first = _channel_finite_real_parameter(a, "a")
    second = _channel_finite_real_parameter(b, "b")
    third = _channel_finite_real_parameter(c, "c")
    scalar_type = promote_type(typeof(first), typeof(second), typeof(third))
    omega = _unnormalized_maximally_entangled(3, scalar_type)
    diagonal = scalar_type[
        first + one(first),
        third,
        second,
        second,
        first + one(first),
        third,
        third,
        second,
        first + one(first),
    ]
    matrix = spdiagm(0 => diagonal) - omega * adjoint(omega)
    return ChoiRepresentation(matrix, 3, 3)
end

"""
    reduction_map(dim, k=1)

Return the sparse Choi representation of
``R_k(X) = k*tr(X)*I - X``. The default `k=1` is the reduction map.

This is a general linear-map constructor. Positivity depends on `k`; the
constructor does not infer or claim positivity, complete positivity, or trace
preservation.
"""
function reduction_map(dim, k=1)
    dimension = _map_dimension(dim, "dim")
    parameter = _channel_finite_real_parameter(k, "k")
    scalar_type = typeof(parameter)
    omega = _unnormalized_maximally_entangled(dimension, scalar_type)
    matrix =
        parameter * spdiagm(0 => fill(one(parameter), dimension^2)) - omega * adjoint(omega)
    return ChoiRepresentation(matrix, dimension, dimension)
end
