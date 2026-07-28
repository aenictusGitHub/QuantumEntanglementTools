# Source-informed independent Julia implementations based on the specifications
# and QETLAB MaxEntangled.m, Bell.m, GHZState.m, WState.m, DickeState.m,
# IsotropicState.m, WernerState.m, HorodeckiState.m, GisinState.m,
# BreuerState.m, BrauerStates.m, and ChessboardState.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and named coauthors,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.

export maximally_entangled,
    bell_state,
    ghz_state,
    w_state,
    dicke_state,
    isotropic_state,
    werner_state,
    horodecki_state,
    gisin_state,
    breuer_state,
    brauer_states,
    chessboard_state

function _state_coefficients(coefficients, expected::Int, default_value)
    if coefficients === nothing
        return fill(default_value, expected)
    elseif coefficients isa AbstractVector || coefficients isa Tuple
        length(coefficients) == expected || throw(
            DimensionMismatch(
                "coefficients has length $(length(coefficients)); expected $expected"
            ),
        )
        all(value -> value isa Number, coefficients) ||
            throw(ArgumentError("coefficients must contain only numbers"))
        return collect(coefficients)
    end
    return throw(
        ArgumentError(
            "coefficients must be `nothing` or a tuple/vector of numbers; got $(typeof(coefficients))",
        ),
    )
end

function _finite_real_parameter(value, name::AbstractString)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be real; got $(repr(value))"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    return value
end

function _unit_interval_parameter(value, name::AbstractString)
    checked = _finite_real_parameter(value, name)
    zero(checked) <= checked <= one(checked) ||
        throw(ArgumentError("$name must lie in [0, 1]; got $checked"))
    return checked
end

"""
    maximally_entangled(dim; normalized=true, sparse_output=false, T=Float64)

Return the standard bipartite vector
``\\sum_{j=0}^{dim-1}|j,j\\rangle``, normalized by `sqrt(dim)` by default.

Subsystem 1 is the slowest-varying tensor factor, so nonzeros occur at
indices `1, dim+2, 2dim+3, ...`.  `T` controls the real element type.
"""
function maximally_entangled(
    dim; normalized::Bool=true, sparse_output::Bool=false, T::Type{<:AbstractFloat}=Float64
)
    dimension = _positive_int(dim, "dim")
    total = _checked_power(dimension, 2, "dim")
    indices = [1 + (position - 1) * (dimension + 1) for position in 1:dimension]
    value = normalized ? inv(sqrt(T(dimension))) : one(T)
    result = sparsevec(indices, fill(value, dimension), total)
    return sparse_output ? result : Vector(result)
end

"""
    bell_state(index=0; normalized=true, sparse_output=false, T=Float64)

Return Bell state `0` through `3` in the order
``|00\\rangle+|11\\rangle``, ``|00\\rangle-|11\\rangle``,
``|01\\rangle+|10\\rangle``, and ``|01\\rangle-|10\\rangle``.
"""
function bell_state(
    index=0;
    normalized::Bool=true,
    sparse_output::Bool=false,
    T::Type{<:AbstractFloat}=Float64,
)
    checked_index = _operator_index(index, 3, "index")
    scale = normalized ? inv(sqrt(T(2))) : one(T)
    indices = checked_index < 2 ? [1, 4] : [2, 3]
    values = T[scale, isodd(checked_index) ? -scale : scale]
    result = sparsevec(indices, values, 4)
    return sparse_output ? result : Vector(result)
end

"""
    ghz_state(dim, parties; coefficients=nothing, sparse_output=true, T=Float64)

Construct ``\\sum_j c_j |j\\rangle^{\\otimes parties}`` in local dimension
`dim`.  The default coefficients are all `1/sqrt(dim)`.  Supplied
coefficients are used exactly and are never silently normalized.
"""
function ghz_state(
    dim,
    parties;
    coefficients=nothing,
    sparse_output::Bool=true,
    T::Type{<:AbstractFloat}=Float64,
)
    dimension = _positive_int(dim, "dim")
    party_count = _positive_int(parties, "parties")
    total = _checked_power(dimension, party_count, "dim")
    default_value = inv(sqrt(T(dimension)))
    values = _state_coefficients(coefficients, dimension, default_value)

    repeated_digit_stride = 0
    power = 1
    for _ in 1:party_count
        repeated_digit_stride += power
        power = Base.checked_mul(power, dimension)
    end
    indices = [1 + (digit - 1) * repeated_digit_stride for digit in 1:dimension]
    result = sparsevec(indices, values, total)
    return sparse_output ? result : Vector(result)
end

"""
    w_state(parties; coefficients=nothing, sparse_output=true, T=Float64)

Construct the `parties`-qubit one-excitation W state.  Coefficient `j`
multiplies the term whose excitation is on subsystem `j`.  The default is
the normalized equal superposition; supplied coefficients are not altered.
"""
function w_state(
    parties;
    coefficients=nothing,
    sparse_output::Bool=true,
    T::Type{<:AbstractFloat}=Float64,
)
    party_count = _positive_int(parties, "parties")
    party_count >= 2 || throw(ArgumentError("parties must be at least 2 for a W state"))
    total = _checked_power(2, party_count, "parties")
    default_value = inv(sqrt(T(party_count)))
    values = _state_coefficients(coefficients, party_count, default_value)
    indices = [
        1 + _checked_power(2, party_count - party, "parties") for party in 1:party_count
    ]
    result = sparsevec(indices, values, total)
    return sparse_output ? result : Vector(result)
end

function _dicke_indices(parties::Int, excitations::Int)
    indices = Int[]
    selected = Int[]
    function visit(first_position::Int, remaining::Int)
        if remaining == 0
            zero_based = 0
            for position in selected
                zero_based += _checked_power(2, parties - position, "parties")
            end
            push!(indices, zero_based + 1)
            return nothing
        end
        final_start = parties - remaining + 1
        for position in first_position:final_start
            push!(selected, position)
            visit(position + 1, remaining - 1)
            pop!(selected)
        end
    end
    visit(1, excitations)
    return indices
end

"""
    dicke_state(parties, excitations=1;
                normalized=true, sparse_output=true, T=Float64)

Return the equal superposition of all computational-basis states containing
exactly `excitations` ones among `parties` qubits.  The vector contains
`binomial(parties, excitations)` nonzeros.
"""
function dicke_state(
    parties,
    excitations=1;
    normalized::Bool=true,
    sparse_output::Bool=true,
    T::Type{<:AbstractFloat}=Float64,
)
    party_count = _positive_int(parties, "parties")
    excitation_count = _nonnegative_int(excitations, "excitations")
    excitation_count <= party_count || throw(
        ArgumentError("excitations=$excitation_count must not exceed parties=$party_count"),
    )
    total = _checked_power(2, party_count, "parties")
    indices = _dicke_indices(party_count, excitation_count)
    value = normalized ? inv(sqrt(T(length(indices)))) : one(T)
    result = sparsevec(indices, fill(value, length(indices)), total)
    return sparse_output ? result : Vector(result)
end

"""
    isotropic_state(dim, alpha; sparse_output=true)

Construct
``(1-alpha)I/dim^2 + alpha |Phi\\rangle\\langle Phi|``.

For `dim > 1`, the physical positivity range
``-1/(dim^2-1) <= alpha <= 1`` is enforced.  In the degenerate
one-dimensional case every finite `alpha` produces the same unit state.
No clipping or repair is performed.
"""
function isotropic_state(dim, alpha; sparse_output::Bool=true)
    dimension = _positive_int(dim, "dim")
    parameter = _finite_real_parameter(alpha, "alpha")
    total = _checked_power(dimension, 2, "dim")
    if dimension > 1
        lower = -one(parameter) / (total - 1)
        lower <= parameter <= one(parameter) || throw(
            ArgumentError(
                "alpha must lie in [$lower, 1] for dim=$dimension; got $parameter"
            ),
        )
    end

    scalar = (one(parameter) - parameter) / total
    identity_part = spdiagm(0 => fill(scalar, total))
    entangled_indices = [1 + (position - 1) * (dimension + 1) for position in 1:dimension]
    unnormalized_phi = sparsevec(entangled_indices, fill(one(parameter), dimension), total)
    result =
        identity_part + (parameter / dimension) * (unnormalized_phi * unnormalized_phi')
    return sparse_output ? result : Matrix(result)
end

"""
    werner_state(dim, alpha; sparse_output=true)

Construct the normalized bipartite Werner state
``(I-alpha*S)/(dim*(dim-alpha))``, where `S` swaps two `dim`-dimensional
systems.  Positivity requires and this constructor enforces
`-1 <= alpha <= 1`.
"""
function werner_state(dim, alpha; sparse_output::Bool=true)
    dimension = _positive_int(dim, "dim")
    dimension >= 2 || throw(ArgumentError("dim must be at least 2 for a Werner state"))
    parameter = _finite_real_parameter(alpha, "alpha")
    -one(parameter) <= parameter <= one(parameter) ||
        throw(ArgumentError("alpha must lie in [-1, 1]; got $parameter"))
    total = _checked_power(dimension, 2, "dim")
    identity_matrix = spdiagm(0 => fill(one(parameter), total))
    raw_swap = swap_operator((dimension, dimension); sparse_output=true)
    swap = convert(SparseMatrixCSC{typeof(one(parameter)),Int}, raw_swap)
    result = (identity_matrix - parameter * swap) / (dimension * (dimension - parameter))
    return sparse_output ? result : Matrix(result)
end

"""
    horodecki_state(a; dims=(3, 3))

Return the seminal PPT-entangled Horodecki family in either `3 x 3` or
`2 x 4` local dimensions.  `a` must lie in `[0,1]`; the endpoints are
separable.  The basis ordering follows the package subsystem convention.
"""
function horodecki_state(a; dims=(3, 3))
    parameter = _unit_interval_parameter(a, "a")
    local_dims = if dims isa Tuple || dims isa AbstractVector
        Tuple(dims)
    else
        throw(ArgumentError("dims must be a tuple or vector"))
    end
    local_dims == (3, 3) ||
        local_dims == (2, 4) ||
        throw(ArgumentError("dims must be (3, 3) or (2, 4); got $local_dims"))

    b = (one(parameter) + parameter) / 2
    c = sqrt(one(parameter) - parameter^2) / 2
    if local_dims == (3, 3)
        normalization = inv(8 * parameter + one(parameter))
        return normalization * [
            parameter 0 0 0 parameter 0 0 0 parameter
            0 parameter 0 0 0 0 0 0 0
            0 0 parameter 0 0 0 0 0 0
            0 0 0 parameter 0 0 0 0 0
            parameter 0 0 0 parameter 0 0 0 parameter
            0 0 0 0 0 parameter 0 0 0
            0 0 0 0 0 0 b 0 c
            0 0 0 0 0 0 0 parameter 0
            parameter 0 0 0 parameter 0 c 0 b
        ]
    end

    normalization = inv(7 * parameter + one(parameter))
    return normalization * [
        parameter 0 0 0 0 parameter 0 0
        0 parameter 0 0 0 0 parameter 0
        0 0 parameter 0 0 0 0 parameter
        0 0 0 parameter 0 0 0 0
        0 0 0 0 b 0 0 c
        parameter 0 0 0 0 parameter 0 0
        0 parameter 0 0 0 0 parameter 0
        0 0 parameter 0 c 0 0 b
    ]
end

"""
    gisin_state(lambda, theta)

Construct Gisin's two-qubit mixed-state family.  `lambda` is a mixing
probability in `[0,1]`; `theta` is any finite real angle in radians.
"""
function gisin_state(lambda, theta)
    weight = _unit_interval_parameter(lambda, "lambda")
    angle = _finite_real_parameter(theta, "theta")
    sine = sin(angle)
    cosine = cos(angle)
    cross = -sin(2 * angle) / 2
    rho_theta = [
        zero(cross) 0 0 0
        0 sine^2 cross 0
        0 cross cosine^2 0
        0 0 0 zero(cross)
    ]
    endpoint_mix = Diagonal([one(weight) / 2, zero(weight), zero(weight), one(weight) / 2])
    return weight * rho_theta + (one(weight) - weight) * endpoint_mix
end

"""
    breuer_state(dim, lambda; sparse_output=true)

Construct the even-local-dimension Breuer family as a convex mixture of a
specific maximally entangled projector and the normalized symmetric
projector.  `lambda` must lie in `[0,1]`.
"""
function breuer_state(dim, lambda; sparse_output::Bool=true)
    dimension = _positive_int(dim, "dim")
    iseven(dimension) || throw(ArgumentError("dim must be even; got $dimension"))
    dimension >= 2 || throw(ArgumentError("dim must be at least 2"))
    weight = _unit_interval_parameter(lambda, "lambda")
    total = _checked_power(dimension, 2, "dim")

    coefficient_type = typeof(float(weight))
    rows = Int[]
    values = coefficient_type[]
    for column in 1:dimension
        row = dimension - column + 1
        push!(rows, (column - 1) * dimension + row)
        sign = isodd(row) ? -one(coefficient_type) : one(coefficient_type)
        push!(values, sign / sqrt(coefficient_type(dimension)))
    end
    psi = sparsevec(rows, values, total)
    symmetric = symmetric_projector(dimension, 2; sparse_output=true)
    result =
        weight * (psi * psi') +
        (one(weight) - weight) *
        (coefficient_type(2) / (dimension * (dimension + 1))) *
        symmetric
    return sparse_output ? result : Matrix(result)
end

function _perfect_matchings(vertex_count::Int)
    vertex_count >= 0 && iseven(vertex_count) ||
        throw(ArgumentError("vertex_count must be a nonnegative even integer"))
    function recurse(vertices::Vector{Int})
        isempty(vertices) && return [Tuple{Int,Int}[]]
        first_vertex = first(vertices)
        matchings = Vector{Vector{Tuple{Int,Int}}}()
        for partner_position in 2:length(vertices)
            partner = vertices[partner_position]
            remaining = [
                vertices[position] for position in eachindex(vertices) if
                position != 1 && position != partner_position
            ]
            for tail in recurse(remaining)
                push!(matchings, vcat([(first_vertex, partner)], tail))
            end
        end
        return matchings
    end
    return recurse(collect(1:vertex_count))
end

"""
    brauer_states(dim, pairs; T=Float64)

Return all unnormalized Brauer vectors as columns of a sparse matrix.
There are `(2*pairs-1)!!` columns; each is the tensor product of `pairs`
unnormalized maximally entangled pairs arranged according to one perfect
matching of the `2*pairs` labeled subsystems.
"""
function brauer_states(dim, pairs; T::Type{<:Number}=Float64)
    dimension = _positive_int(dim, "dim")
    pair_count = _positive_int(pairs, "pairs")
    party_count = _checked_product((2, pair_count), "pairs")
    total = _checked_power(dimension, party_count, "dim")
    assignments = _checked_power(dimension, pair_count, "dim")
    matchings = _perfect_matchings(party_count)
    nonzeros = _checked_product((assignments, length(matchings)), "Brauer nonzeros")

    rows = Vector{Int}(undef, nonzeros)
    columns = Vector{Int}(undef, nonzeros)
    values = fill(one(T), nonzeros)
    cursor = 1
    dims = ntuple(_ -> dimension, party_count)
    digits = Vector{Int}(undef, party_count)
    assignment_digits = Vector{Int}(undef, pair_count)
    for (column, matching) in Base.pairs(matchings)
        for assignment in 0:(assignments - 1)
            residual = assignment
            for pair in pair_count:-1:1
                assignment_digits[pair] = mod(residual, dimension) + 1
                residual = div(residual, dimension)
            end
            for (pair, (first_party, second_party)) in Base.pairs(matching)
                digits[first_party] = assignment_digits[pair]
                digits[second_party] = assignment_digits[pair]
            end
            rows[cursor] = basis_to_linear(Tuple(digits), dims)
            columns[cursor] = column
            cursor += 1
        end
    end
    return sparse(rows, columns, values, total, length(matchings))
end

"""
    chessboard_state(a, b, c, d, m, n; s=nothing, t=nothing)

Construct the normalized `3 x 3` chessboard state from four rank-one
vectors.  Defaults are `s = a*conj(c)/conj(n)` and `t = a*d/m`; zero
denominators are rejected when those defaults are requested.

The result is positive semidefinite by construction.  This routine does not
claim or test the PPT property for arbitrary explicit `s` and `t`.
"""
function chessboard_state(a, b, c, d, m, n; s=nothing, t=nothing)
    all(value -> value isa Number, (a, b, c, d, m, n)) ||
        throw(ArgumentError("all chessboard parameters must be numbers"))
    all(isfinite, (a, b, c, d, m, n)) ||
        throw(ArgumentError("all chessboard parameters must be finite"))
    s_value = if s === nothing
        iszero(n) && throw(ArgumentError("n must be nonzero when s is omitted"))
        a * conj(c) / conj(n)
    else
        s isa Number || throw(ArgumentError("s must be a number; got $(repr(s))"))
        isfinite(s) || throw(ArgumentError("s must be finite; got $(repr(s))"))
        s
    end
    t_value = if t === nothing
        iszero(m) && throw(ArgumentError("m must be nonzero when t is omitted"))
        a * d / m
    else
        t isa Number || throw(ArgumentError("t must be a number; got $(repr(t))"))
        isfinite(t) || throw(ArgumentError("t must be finite; got $(repr(t))"))
        t
    end

    vectors = (
        [m, 0, s_value, 0, n, 0, 0, 0, 0],
        [0, a, 0, b, 0, c, 0, 0, 0],
        [conj(n), 0, 0, 0, -conj(m), 0, t_value, 0, 0],
        [0, conj(b), 0, -conj(a), 0, 0, 0, d, 0],
    )
    rho = sum(conj.(vector) * transpose(vector) for vector in vectors)
    normalization = real(tr(rho))
    iszero(normalization) &&
        throw(ArgumentError("the supplied parameters produce the zero matrix"))
    return rho / normalization
end
