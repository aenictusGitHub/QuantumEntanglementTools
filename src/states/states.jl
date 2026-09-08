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

const _DICKE_DEFAULT_MAX_NONZEROS = 1_000_000
const _DICKE_DEFAULT_MAX_DENSE_ENTRIES = 10_000_000
const _DICKE_DEFAULT_MAX_WORK = 100_000_000

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

function _equal_superposition_amplitude(count::Int, ::Type{T}) where {T<:AbstractFloat}
    count > 0 || throw(ArgumentError("count must be positive; got $count"))
    # Float16 overflows on ordinary dimensions above 65_504.  Evaluate its
    # normalization in Float32, while retaining the requested arithmetic for
    # Float32, Float64, BigFloat, and custom AbstractFloat implementations.
    work_type = T === Float16 ? Float32 : T
    work_count = convert(work_type, count)
    isfinite(work_count) || throw(
        ArgumentError(
            "count=$count is not representable in the normalization work type $work_type",
        ),
    )
    value = convert(T, inv(sqrt(work_count)))
    isfinite(value) && !iszero(value) || throw(
        ArgumentError(
            "the normalized amplitude for count=$count is not representable as a nonzero $T",
        ),
    )
    return value
end

"""
    maximally_entangled(dim; normalized=true, sparse_output=false, T=Float64)

Return the standard bipartite vector
``\\sum_{j=0}^{dim-1}|j,j\\rangle``, normalized by `sqrt(dim)` by default.

Subsystem 1 is the slowest-varying tensor factor, so nonzeros occur at
indices `1, dim+2, 2dim+3, ...`.  `T` controls the real element type.
Normalization throws `ArgumentError` if `T` cannot represent the required
amplitude as a finite, nonzero value.
"""
function maximally_entangled(
    dim; normalized::Bool=true, sparse_output::Bool=false, T::Type{<:AbstractFloat}=Float64
)
    dimension = _positive_int(dim, "dim")
    total = _checked_power(dimension, 2, "dim")
    indices = [1 + (position - 1) * (dimension + 1) for position in 1:dimension]
    value = normalized ? _equal_superposition_amplitude(dimension, T) : one(T)
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
    scale = normalized ? _equal_superposition_amplitude(2, T) : one(T)
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
Default normalization throws `ArgumentError` if `T` cannot represent the
required amplitude as a finite, nonzero value.
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
    default_value = _equal_superposition_amplitude(dimension, T)
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
Default normalization throws `ArgumentError` if `T` cannot represent the
required amplitude as a finite, nonzero value.
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
    default_value = _equal_superposition_amplitude(party_count, T)
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

function _dicke_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("$name must be a positive integer or nothing, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer or nothing"))
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _dicke_check_resource(
    planned::BigInt, limit, name::AbstractString, resource::AbstractString
)
    limit === nothing && return nothing
    planned <= limit || throw(
        ArgumentError(
            "dicke_state requires $planned $resource, exceeding $name=$limit; " *
            "raise the explicit guard only after reviewing the combinatorial cost",
        ),
    )
    return nothing
end

"""
    dicke_state(
        parties, excitations=1;
        normalized=true, sparse_output=true, T=Float64,
        max_nonzeros=1_000_000, max_dense_entries=10_000_000,
        max_work=100_000_000,
    )

Return the equal superposition of all computational-basis states containing
exactly `excitations` ones among `parties` qubits.  The vector contains
`binomial(parties, excitations)` nonzeros.

`max_nonzeros` and `max_work` reject excessive combination enumeration before
index allocation. `max_dense_entries` bounds the ambient vector length when
`sparse_output=false`. Set an individual guard to `nothing` only after
reviewing the requested memory and work; the ambient dimension must always fit
`Int`. Normalization throws `ArgumentError` if `T` cannot represent the
required amplitude as a finite, nonzero value.
"""
function dicke_state(
    parties,
    excitations=1;
    normalized::Bool=true,
    sparse_output::Bool=true,
    T::Type{<:AbstractFloat}=Float64,
    max_nonzeros=_DICKE_DEFAULT_MAX_NONZEROS,
    max_dense_entries=_DICKE_DEFAULT_MAX_DENSE_ENTRIES,
    max_work=_DICKE_DEFAULT_MAX_WORK,
)
    party_count = _positive_int(parties, "parties")
    excitation_count = _nonnegative_int(excitations, "excitations")
    excitation_count <= party_count || throw(
        ArgumentError("excitations=$excitation_count must not exceed parties=$party_count"),
    )
    nonzero_limit = _dicke_limit(max_nonzeros, "max_nonzeros")
    dense_limit = _dicke_limit(max_dense_entries, "max_dense_entries")
    work_limit = _dicke_limit(max_work, "max_work")
    total = _checked_power(2, party_count, "parties")
    nonzero_count = binomial(BigInt(party_count), excitation_count)
    nonzero_count <= typemax(Int) || throw(
        ArgumentError(
            "dicke_state requires $nonzero_count nonzeros, which cannot be " *
            "represented as an array length",
        ),
    )
    _dicke_check_resource(nonzero_count, nonzero_limit, "max_nonzeros", "nonzeros")
    if !sparse_output
        _dicke_check_resource(
            BigInt(total), dense_limit, "max_dense_entries", "dense output entries"
        )
    end
    work =
        BigInt(party_count) +
        nonzero_count * (BigInt(1) + BigInt(excitation_count) * BigInt(party_count))
    _dicke_check_resource(work, work_limit, "max_work", "estimated scalar operations")
    indices = _dicke_indices(party_count, excitation_count)
    value = normalized ? _equal_superposition_amplitude(length(indices), T) : one(T)
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

function _werner_party_count(permutation_count::Int)
    permutation_count >= 2 ||
        throw(ArgumentError("multipartite alpha must contain at least one parameter"))
    factorial_value = 1
    parties = 1
    while factorial_value < permutation_count
        parties += 1
        factorial_value = try
            Base.checked_mul(factorial_value, parties)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("length(alpha) + 1 exceeds representable factorials"))
        end
    end
    factorial_value == permutation_count || throw(
        ArgumentError(
            "multipartite alpha must contain p! - 1 parameters for some p >= 2; " *
            "got $(permutation_count - 1)",
        ),
    )
    return parties
end

function _werner_next_permutation!(permutation::Vector{Int})
    pivot = length(permutation) - 1
    while pivot >= 1 && permutation[pivot] >= permutation[pivot + 1]
        pivot -= 1
    end
    pivot == 0 && return false
    successor = length(permutation)
    while permutation[successor] <= permutation[pivot]
        successor -= 1
    end
    permutation[pivot], permutation[successor] = permutation[successor], permutation[pivot]
    reverse!(permutation, pivot + 1, length(permutation))
    return true
end

function _werner_permutations(parties::Int)
    permutation = collect(1:parties)
    permutations = NTuple{parties,Int}[]
    while _werner_next_permutation!(permutation)
        push!(permutations, Tuple(permutation))
    end
    return permutations
end

function _werner_limit(value, name::AbstractString)
    value === nothing && return nothing
    return _nonnegative_int(value, name)
end

@doc raw"""
    werner_state(
        dim,
        alpha::AbstractVector;
        sparse_output=true,
        max_permutations=40_320,
        max_nonzeros=1_000_000,
        max_dense_entries=1_000_000,
        max_work=100_000_000,
        atol=nothing,
        rtol=nothing,
    )

Construct the normalized multipartite Werner operator

```math
\rho =
\frac{I-\sum_{j=2}^{p!}\alpha_{j-1}P_j}
     {\operatorname{tr}\left(I-\sum_{j=2}^{p!}\alpha_{j-1}P_j\right)},
```

where `length(alpha) == p! - 1` and `P_j` follows lexicographic permutation
order. Every local subsystem has dimension `dim`, and subsystem `1` is the
slowest-varying tensor factor. A one-entry vector delegates to the verified
bipartite scalar family.

The coefficient of each permutation must be the conjugate of the coefficient
of its inverse, making the unnormalized operator exactly Hermitian. Positivity
is certified without densification when `sum(abs, alpha) <= 1`. Otherwise a
dense Hermitian eigendecomposition is performed only for BLAS floating types
and only within `max_dense_entries` and `max_work`; a negative or numerical
boundary eigenvalue is rejected rather than clipped. The trace must be finite,
real, and strictly positive. No normalization other than the displayed trace
division, symmetrization, projection, or coefficient repair is performed.

The implementation intentionally corrects the pinned multipartite loop, which
overwrites the accumulator at every permutation and therefore retains only the
last parameter. Permutation count, construction nonzeros, dense validation,
and work are checked before their corresponding allocations.
""" function werner_state(
    dim,
    alpha::AbstractVector;
    sparse_output::Bool=true,
    max_permutations=40_320,
    max_nonzeros=1_000_000,
    max_dense_entries=1_000_000,
    max_work=100_000_000,
    atol=nothing,
    rtol=nothing,
)
    Base.require_one_based_indexing(alpha)
    length(alpha) == 1 && return werner_state(dim, only(alpha); sparse_output=sparse_output)
    dimension = _positive_int(dim, "dim")
    dimension >= 2 || throw(ArgumentError("dim must be at least 2 for a Werner state"))
    isempty(alpha) &&
        throw(ArgumentError("multipartite alpha must contain p! - 1 parameters"))
    all(parameter -> parameter isa Number && !(parameter isa Bool), alpha) ||
        throw(ArgumentError("alpha must contain only numeric parameters other than Bool"))
    all(isfinite, alpha) ||
        throw(ArgumentError("alpha must contain only finite parameters"))

    permutation_count = try
        Base.checked_add(length(alpha), 1)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("length(alpha) + 1 exceeds typemax(Int)"))
    end
    parties = _werner_party_count(permutation_count)
    permutation_limit = _werner_limit(max_permutations, "max_permutations")
    if permutation_limit !== nothing && permutation_count > permutation_limit
        throw(
            ArgumentError(
                "multipartite Werner construction requires $permutation_count " *
                "permutations including identity, exceeding " *
                "max_permutations=$permutation_limit",
            ),
        )
    end

    total_dimension = _checked_power(dimension, parties, "dim")
    stored_entry_estimate = try
        Base.checked_mul(permutation_count, total_dimension)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the Werner construction nonzero estimate exceeds Int"))
    end
    nonzero_limit = _werner_limit(max_nonzeros, "max_nonzeros")
    if nonzero_limit !== nothing && stored_entry_estimate > nonzero_limit
        throw(
            ArgumentError(
                "multipartite Werner construction may require " *
                "$stored_entry_estimate stored entries, exceeding " *
                "max_nonzeros=$nonzero_limit",
            ),
        )
    end
    work_limit = _werner_limit(max_work, "max_work")
    if work_limit !== nothing && stored_entry_estimate > work_limit
        throw(
            ArgumentError(
                "multipartite Werner construction requires estimated work " *
                "$stored_entry_estimate, exceeding max_work=$work_limit",
            ),
        )
    end
    dense_entries = try
        Base.checked_mul(total_dimension, total_dimension)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("the Werner dense entry count exceeds Int"))
    end
    dense_limit = _werner_limit(max_dense_entries, "max_dense_entries")
    if !sparse_output && dense_limit !== nothing && dense_entries > dense_limit
        throw(
            ArgumentError(
                "dense Werner output requires $dense_entries entries, exceeding " *
                "max_dense_entries=$dense_limit",
            ),
        )
    end

    promoted_parameter_type = promote_type(map(typeof, alpha)...)
    coefficient_type = typeof(one(promoted_parameter_type) / one(promoted_parameter_type))
    coefficients = coefficient_type.(alpha)
    permutations = _werner_permutations(parties)
    length(permutations) == length(coefficients) ||
        error("internal multipartite Werner permutation enumeration is inconsistent")
    permutation_indices = Dict(
        permutation => index for (index, permutation) in pairs(permutations)
    )
    for (index, permutation) in pairs(permutations)
        inverse_permutation = Tuple(invperm(collect(permutation)))
        inverse_index = permutation_indices[inverse_permutation]
        coefficients[index] == conj(coefficients[inverse_index]) || throw(
            ArgumentError(
                "alpha[$index] for permutation $permutation must equal the " *
                "conjugate coefficient of inverse permutation " *
                "$inverse_permutation; coefficients are never symmetrized",
            ),
        )
    end

    local_dims = ntuple(_ -> dimension, parties)
    raw = spdiagm(0 => fill(one(coefficient_type), total_dimension))
    for (coefficient, permutation) in zip(coefficients, permutations)
        raw -=
            coefficient * permutation_operator(
                local_dims, permutation; T=coefficient_type, sparse_output=true
            )
    end
    normalization = tr(raw)
    isfinite(normalization) ||
        throw(DomainError(normalization, "multipartite Werner trace must be finite"))
    isreal(normalization) || throw(
        DomainError(
            normalization,
            "multipartite Werner trace must be exactly real; coefficients are " *
            "never repaired",
        ),
    )
    real_normalization = real(normalization)
    real_normalization > zero(real_normalization) || throw(
        DomainError(
            real_normalization,
            "multipartite Werner unnormalized operator must have positive trace",
        ),
    )

    coefficient_norm = sum(abs, coefficients)
    if coefficient_norm > one(coefficient_norm)
        coefficient_type <: LinearAlgebra.BlasFloat || throw(
            ArgumentError(
                "parameters with sum(abs, alpha) > 1 require a dense Hermitian " *
                "eigendecomposition and therefore Float32, Float64, ComplexF32, " *
                "or ComplexF64 coefficients",
            ),
        )
        if dense_limit !== nothing && dense_entries > dense_limit
            throw(
                ArgumentError(
                    "PSD validation requires $dense_entries dense entries, " *
                    "exceeding max_dense_entries=$dense_limit",
                ),
            )
        end
        spectral_work = try
            Base.checked_mul(dense_entries, total_dimension)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("the Werner spectral work estimate exceeds Int"))
        end
        total_work = try
            Base.checked_add(stored_entry_estimate, spectral_work)
        catch err
            err isa OverflowError || rethrow()
            throw(ArgumentError("the Werner total work estimate exceeds Int"))
        end
        if work_limit !== nothing && total_work > work_limit
            throw(
                ArgumentError(
                    "PSD validation requires estimated total work $total_work, " *
                    "exceeding max_work=$work_limit",
                ),
            )
        end
        eigenvalues = eigvals(Hermitian(Matrix(raw)))
        minimum_eigenvalue = minimum(eigenvalues)
        absolute, relative = _tierd_tolerances(
            typeof(real(zero(coefficient_type))), atol, rtol
        )
        scale = maximum(abs, eigenvalues; init=zero(eltype(eigenvalues)))
        threshold = absolute + relative * scale
        minimum_eigenvalue >= zero(minimum_eigenvalue) || throw(
            DomainError(
                minimum_eigenvalue,
                if minimum_eigenvalue >= -threshold
                    "multipartite Werner positivity lies on a numerical boundary; " *
                    "the operator is not projected or clipped"
                else
                    "multipartite Werner parameters produce a non-positive operator"
                end,
            ),
        )
    end

    result = raw / real_normalization
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

function _brauer_matching_count(pair_count::Int, max_matchings, max_nonzeros)
    matching_limit = if max_matchings === nothing
        nothing
    else
        BigInt(_positive_int(max_matchings, "max_matchings"))
    end
    nonzero_limit = if max_nonzeros === nothing
        nothing
    else
        BigInt(_positive_int(max_nonzeros, "max_nonzeros"))
    end
    count = BigInt(1)
    for factor in 1:2:(2 * pair_count - 1)
        count *= factor
        matching_limit !== nothing &&
            count > matching_limit &&
            _brauer_complexity_guard(
                count, matching_limit, "max_matchings", "perfect matchings"
            )
        nonzero_limit !== nothing &&
            count > nonzero_limit &&
            _brauer_complexity_guard(count, nonzero_limit, "max_nonzeros", "stored entries")
    end
    return count
end

function _brauer_complexity_guard(
    planned::BigInt, limit, keyword::AbstractString, resource::AbstractString
)
    limit === nothing && return nothing
    checked_limit = _positive_int(limit, keyword)
    planned <= checked_limit || throw(
        ArgumentError(
            "brauer_states would generate $planned $resource, exceeding " *
            "$keyword=$checked_limit; raise the guard only after reviewing " *
            "the combinatorial memory cost",
        ),
    )
    return nothing
end

"""
    brauer_states(
        dim, pairs;
        T=Float64, max_matchings=100_000, max_nonzeros=1_000_000
    )

Return all unnormalized Brauer vectors as columns of a sparse matrix.
There are `(2*pairs-1)!!` columns; each is the tensor product of `pairs`
unnormalized maximally entangled pairs arranged according to one perfect
matching of the `2*pairs` labeled subsystems. The output has `dim^pairs`
stored entries per column.

`max_matchings` and `max_nonzeros` reject excessive work before matching
generation or allocation. Setting either guard to `nothing` explicitly disables
that guard; overflow and allocation failure remain possible for sufficiently
large requests.
"""
function brauer_states(
    dim, pairs; T::Type{<:Number}=Float64, max_matchings=100_000, max_nonzeros=1_000_000
)
    dimension = _positive_int(dim, "dim")
    pair_count = _positive_int(pairs, "pairs")
    party_count = _checked_product((2, pair_count), "pairs")
    matching_count_big = _brauer_matching_count(pair_count, max_matchings, max_nonzeros)
    _brauer_complexity_guard(
        matching_count_big, max_matchings, "max_matchings", "perfect matchings"
    )
    matching_count_big <= typemax(Int) ||
        throw(ArgumentError("Brauer matching count exceeds typemax(Int)"))
    matching_count = Int(matching_count_big)
    assignments = _checked_power(dimension, pair_count, "dim")
    nonzero_count_big = BigInt(assignments) * matching_count_big
    _brauer_complexity_guard(
        nonzero_count_big, max_nonzeros, "max_nonzeros", "stored entries"
    )
    nonzero_count_big <= typemax(Int) ||
        throw(ArgumentError("Brauer stored-entry count exceeds typemax(Int)"))
    nonzeros = Int(nonzero_count_big)
    total = _checked_power(dimension, party_count, "dim")
    matchings = _perfect_matchings(party_count)
    length(matchings) == matching_count ||
        error("internal Brauer matching enumeration count is inconsistent")

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
    return sparse(rows, columns, values, total, matching_count)
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
