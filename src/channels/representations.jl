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

function dual_channel(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return OperatorSumRepresentation(
        [copy(adjoint(operator)) for operator in left_operators],
        [copy(adjoint(operator)) for operator in right_operators],
    )
end

function dual_channel(map::SuperoperatorRepresentation)
    space = operator_space(map)
    dual_space = OperatorSpace(output_size(space), input_size(space))
    return SuperoperatorRepresentation(
        copy(adjoint(_validated_representation_matrix(map))), dual_space
    )
end

function dual_channel(map::ChoiRepresentation)
    return choi_representation(dual_channel(superoperator_representation(map)))
end

function _complementary_factor_collection(
    factors, output_dimension::Int, input_dimension::Int, ::Type{T}
) where {T}
    environment_dimension = length(factors)
    if all(issparse, factors)
        output_rows = [Int[] for _ in 1:output_dimension]
        input_columns = [Int[] for _ in 1:output_dimension]
        stored_values = [T[] for _ in 1:output_dimension]
        for (factor_index, factor) in enumerate(factors)
            rows, columns, values = findnz(factor)
            for index in eachindex(values)
                output_index = rows[index]
                push!(output_rows[output_index], factor_index)
                push!(input_columns[output_index], columns[index])
                push!(stored_values[output_index], convert(T, values[index]))
            end
        end
        return [
            sparse(
                output_rows[index],
                input_columns[index],
                stored_values[index],
                environment_dimension,
                input_dimension,
            ) for index in 1:output_dimension
        ]
    end

    complementary_factors = [
        Matrix{T}(undef, environment_dimension, input_dimension) for _ in 1:output_dimension
    ]
    @inbounds for output_index in 1:output_dimension,
        factor_index in 1:environment_dimension,
        input_index in 1:input_dimension

        complementary_factors[output_index][factor_index, input_index] = factors[factor_index][
            output_index, input_index
        ]
    end
    return complementary_factors
end

"""
    complementary_channel(map; atol=0, rtol=sqrt(eps(T)),
                          allow_densify=false)

Construct a complementary map from a supplied dilation. For Kraus matrices
`K_r` of shape `(d_out, d_in)`, the complement has environment dimension
equal to the number of factors and matrices `L_a[r,i] = K_r[a,i]`.

A general paired map `sum(A_r * X * B_r')` has a complementary paired map
only when its output operator space is square. Its factors are
`L_a[r,i] = A_r[a,i]` and `R_a[r,j] = B_r[a,j]`. This supports a rectangular
input operator space but rejects unequal output row/column dimensions, where
there is no canonical output trace pairing. This is the explicit correction
to the pinned routine's mismatched `cell2mat`/`mat2cell` dimensions.

Kraus and operator-sum inputs preserve the supplied dilation exactly. Choi
and superoperator inputs first use [`canonical_map_decomposition`](@ref);
their complement therefore depends on that canonical numerical
factorization. Sparse spectral input requires `allow_densify=true`. The
output representation kind matches the input.
"""
function complementary_channel(
    map::KrausRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    _checked_tolerances(map, atol, rtol)
    operators = _validated_kraus_operators(map)
    complementary_operators = _complementary_factor_collection(
        operators, map.output_dim, map.input_dim, eltype(map)
    )
    return KrausRepresentation(complementary_operators)
end

function complementary_channel(
    map::OperatorSumRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    _checked_tolerances(map, atol, rtol)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    space = operator_space(map)
    space.output_rows == space.output_columns || throw(
        ArgumentError(
            "a paired complementary map requires a square output operator " *
            "space; got output_size=$(output_size(space))",
        ),
    )
    complementary_left = _complementary_factor_collection(
        left_operators, space.output_rows, space.input_rows, eltype(map)
    )
    complementary_right = _complementary_factor_collection(
        right_operators, space.output_columns, space.input_columns, eltype(map)
    )
    return OperatorSumRepresentation(complementary_left, complementary_right)
end

function complementary_channel(
    map::ChoiRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    decomposition = canonical_map_decomposition(
        map; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    return choi_representation(complementary_channel(decomposition.representation))
end

function complementary_channel(
    map::SuperoperatorRepresentation;
    atol=zero(_default_rtol(map)),
    rtol=_default_rtol(map),
    allow_densify::Bool=false,
)
    decomposition = canonical_map_decomposition(
        map; atol=atol, rtol=rtol, allow_densify=allow_densify
    )
    return superoperator_representation(complementary_channel(decomposition.representation))
end

function _map_has_sparse_storage(map::KrausRepresentation)
    return all(issparse, _validated_kraus_operators(map))
end

function _map_has_sparse_storage(map::OperatorSumRepresentation)
    left_operators, right_operators = _validated_operator_sum_factors(map)
    return all(issparse, left_operators) && all(issparse, right_operators)
end

function _map_has_sparse_storage(map::ChoiRepresentation)
    return issparse(_validated_representation_matrix(map))
end

function _map_has_sparse_storage(map::SuperoperatorRepresentation)
    return issparse(_validated_representation_matrix(map))
end

"""
    partial_map(input, map, subsystem, dims)
    partial_map(input, map, subsystem, row_dims, column_dims)

Apply `map` to one subsystem of a multipartite operator. Subsystems use
tensor-product order with subsystem 1 slowest-varying. The four-argument form
uses the same dimensions for rows and columns. The five-argument form accepts
independent layouts and requires
`(row_dims[subsystem], column_dims[subsystem]) == input_size(map)`. The
selected output dimensions are replaced by `output_size(map)`.

The implementation permutes the selected subsystem to the final tensor
position, applies the map directly to each operator block, and restores the
original subsystem order. It does not construct a full tensor-product
superoperator. Sparse storage is preserved when both the input and the
representation data are sparse; otherwise the result is dense.
"""
function partial_map(input::AbstractMatrix, map::AbstractMapRepresentation, subsystem, dims)
    layout = _as_layout(dims)
    return _partial_map_rectangular(input, map, subsystem, layout, layout)
end

function partial_map(
    input::AbstractMatrix, map::AbstractMapRepresentation, subsystem, row_dims, column_dims
)
    row_layout = _as_layout(row_dims)
    column_layout = _as_layout(column_dims)
    return _partial_map_rectangular(input, map, subsystem, row_layout, column_layout)
end

function _partial_map_rectangular(
    input::AbstractMatrix,
    map::AbstractMapRepresentation,
    subsystem,
    row_layout::SubsystemLayout,
    column_layout::SubsystemLayout,
)
    length(row_layout) == length(column_layout) || throw(
        DimensionMismatch(
            "row_dims describe $(length(row_layout)) subsystems but column_dims " *
            "describe $(length(column_layout))",
        ),
    )
    _validate_matrix_dimensions(input, row_layout, column_layout)
    _validate_numeric_matrix(input, "partial-map input")
    selected = _normalize_systems(
        subsystem, length(row_layout); name="subsystem", allow_empty=false
    )
    length(selected) == 1 ||
        throw(ArgumentError("partial_map applies to exactly one subsystem"))
    system = only(selected)
    map_input_rows, map_input_columns = input_size(map)
    map_output_rows, map_output_columns = output_size(map)
    row_layout[system] == map_input_rows || throw(
        DimensionMismatch(
            "row_dims[$system]=$(row_layout[system]) does not match map input " *
            "row dimension $map_input_rows",
        ),
    )
    column_layout[system] == map_input_columns || throw(
        DimensionMismatch(
            "column_dims[$system]=$(column_layout[system]) does not match map " *
            "input column dimension $map_input_columns",
        ),
    )

    permutation = Tuple(
        vcat([index for index in 1:length(row_layout) if index != system], [system])
    )
    row_plan = SubsystemPermutationPlan(row_layout.dims, permutation)
    column_plan = SubsystemPermutationPlan(column_layout.dims, permutation)
    permuted = permute_subsystems(input, row_plan, column_plan)
    remainder_row_dims = Tuple(row_layout[index] for index in permutation[1:(end - 1)])
    remainder_column_dims = Tuple(
        column_layout[index] for index in permutation[1:(end - 1)]
    )
    remainder_rows = _checked_product(
        remainder_row_dims, "unmapped row subsystem dimensions"
    )
    remainder_columns = _checked_product(
        remainder_column_dims, "unmapped column subsystem dimensions"
    )
    output_row_count = _checked_product(
        (remainder_rows, map_output_rows), "partial-map output row dimensions"
    )
    output_column_count = _checked_product(
        (remainder_columns, map_output_columns), "partial-map output column dimensions"
    )
    output_type = promote_type(eltype(input), eltype(map))
    sparse_output = issparse(input) && _map_has_sparse_storage(map)
    permuted_output = if sparse_output
        spzeros(output_type, output_row_count, output_column_count)
    else
        zeros(output_type, output_row_count, output_column_count)
    end

    application_map = map isa ChoiRepresentation ? superoperator_representation(map) : map
    @inbounds for block_column in 1:remainder_columns, block_row in 1:remainder_rows
        input_rows = ((block_row - 1) * map_input_rows + 1):(block_row * map_input_rows)
        input_columns =
            ((block_column - 1) * map_input_columns + 1):(block_column * map_input_columns)
        output_rows = ((block_row - 1) * map_output_rows + 1):(block_row * map_output_rows)
        output_columns =
            ((block_column - 1) * map_output_columns + 1):(block_column * map_output_columns)
        permuted_output[output_rows, output_columns] = apply_channel(
            permuted[input_rows, input_columns], application_map
        )
    end

    permuted_output_row_dims = (remainder_row_dims..., map_output_rows)
    permuted_output_column_dims = (remainder_column_dims..., map_output_columns)
    inverse_permutation = Tuple(invperm(collect(permutation)))
    output_row_plan = SubsystemPermutationPlan(
        permuted_output_row_dims, inverse_permutation
    )
    output_column_plan = SubsystemPermutationPlan(
        permuted_output_column_dims, inverse_permutation
    )
    return permute_subsystems(permuted_output, output_row_plan, output_column_plan)
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
