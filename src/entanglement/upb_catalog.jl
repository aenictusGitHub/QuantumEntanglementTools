# Source-informed independent Julia implementation based on the specification
# and QETLAB UPB.m and helpers/one_factorization.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2022 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The family specifications and construction invariants were also checked
# against the primary sources cited by each result.

using Random

const _UPBC_SOURCE_REVISION = "d8589610f00cff106537268dee2e2a1153f3a601"
const _UPBC_DEFAULT_MAX_LOCAL_ENTRIES = 1_000_000
const _UPBC_DEFAULT_MAX_GLOBAL_ENTRIES = 2_000_000
const _UPBC_DEFAULT_MAX_WORK = 50_000_000
const _UPBC_DEFAULT_MAX_MINORS = 250_000
const _UPBC_DEFAULT_MAX_ATTEMPTS = 64

"""
    UPBConstruction

Auditable construction returned by [`upb`](@ref).

`local_factors[p][:, j]` is the local vector of party `p` in product state
`j`, while `global_vectors[:, j]` is their tensor product in the package's
subsystem order. `dimensions`, `cardinality`, and `family` are canonical
construction metadata. The reference fields identify the primary source for
the selected family.

`construction_kind` is `:closed_form`, `:randomized_full_spark`, or
`:complete_product_basis`. It describes how the vectors were obtained; it is
not a claim that floating-point entries are exact. `verification_kind` is
`:exact_structure` for a computational product basis and
`:tolerance_robust` for closed-form or randomized floating constructions.
The three residuals retain the normalization, pairwise-orthogonality, and
tensor-product reconstruction checks performed before the result is returned.

Randomized constructions record the caller-supplied RNG use, attempts, draw
count, and full-spark minors checked. All constructions retain the explicit
entry, work, minor, and retry limits used for the call.
"""
struct UPBConstruction{L,G,D,R,N,T}
    status::Symbol
    family::Symbol
    requested::R
    dimensions::D
    cardinality::Int
    local_factors::L
    global_vectors::G
    reference_key::Symbol
    reference::String
    reference_url::String
    construction_kind::Symbol
    verification_kind::Symbol
    deterministic::Bool
    rng_used::Bool
    normalized::Bool
    pairwise_orthogonal::Bool
    normalization_residual::T
    orthogonality_residual::T
    product_residual::T
    attempts::Int
    random_draws::Int
    minors_checked::Int
    work_used::Int
    max_attempts::Int
    max_local_entries::Union{Nothing,Int}
    max_global_entries::Union{Nothing,Int}
    max_work::Union{Nothing,Int}
    max_minors::Union{Nothing,Int}
    arithmetic::Symbol
    symbolic_simplification::Symbol
    source_revision::String
    notes::N
end

function Base.show(io::IO, construction::UPBConstruction)
    return print(
        io,
        "UPBConstruction(family=",
        construction.family,
        ", dimensions=",
        construction.dimensions,
        ", cardinality=",
        construction.cardinality,
        ", construction_kind=",
        construction.construction_kind,
        ")",
    )
end

"""
    UPBConstructionUnavailable

Exception used when the pinned public contract intentionally reports that the
minimum cardinality is unknown or that a known minimum has no construction in
the catalog. `known_minimum` is retained when the reviewed theorem table
determines it. It is also used with reason `:randomized_search_exhausted` when
a bounded full-spark search uses every permitted attempt. No vectors are
fabricated in any case.
"""
struct UPBConstructionUnavailable{D,S} <: Exception
    reason::Symbol
    dimensions::D
    known_minimum::S
    message::String
end

function Base.showerror(io::IO, err::UPBConstructionUnavailable)
    return print(
        io,
        "UPB construction unavailable (",
        err.reason,
        ") for dimensions ",
        err.dimensions,
        ": ",
        err.message,
    )
end

"""
    UPBResourceLimitError

Exception raised before a UPB construction exceeds an explicit deterministic
entry, work, full-spark-minor, or retry budget.
"""
struct UPBResourceLimitError <: Exception
    resource::Symbol
    required::Union{Nothing,BigInt}
    limit::Int
    message::String
end

function Base.showerror(io::IO, err::UPBResourceLimitError)
    required = err.required === nothing ? "" : " (required $(err.required))"
    return print(
        io,
        "UPB resource limit ",
        err.resource,
        " exceeded",
        required,
        "; limit ",
        err.limit,
        ": ",
        err.message,
    )
end

mutable struct _UPBCBudget
    max_work::Union{Nothing,Int}
    max_minors::Union{Nothing,Int}
    work_used::Int
    minors_checked::Int
    random_draws::Int
end

struct _UPBCLocalConstruction{F,N}
    family::Symbol
    factors::F
    reference_key::Symbol
    construction_kind::Symbol
    attempts::Int
    notes::N
end

function _upbc_limit(value, name::AbstractString; positive::Bool=false)
    value === nothing && return nothing
    value isa Bool && throw(ArgumentError("$name must be an integer or nothing, not Bool"))
    value isa Integer ||
        throw(ArgumentError("$name must be an integer or nothing; got $(typeof(value))"))
    valid = positive ? value > 0 : value >= 0
    valid || throw(ArgumentError("$name must be $(positive ? "positive" : "nonnegative")"))
    value <= typemax(Int) || throw(ArgumentError("$name exceeds typemax(Int)"))
    return Int(value)
end

function _upbc_real_type(real_type)
    real_type isa Type ||
        throw(ArgumentError("real_type must be a concrete floating-point type"))
    isconcretetype(real_type) && real_type <: AbstractFloat || throw(
        ArgumentError(
            "real_type must be a concrete floating-point type such as Float32, " *
            "Float64, or BigFloat; got $real_type",
        ),
    )
    return real_type
end

function _upbc_consume!(budget::_UPBCBudget, amount::Integer, label::AbstractString)
    amount >= 0 || error("internal UPB work estimate is negative")
    required = big(budget.work_used) + big(amount)
    if budget.max_work !== nothing && required > budget.max_work
        throw(
            UPBResourceLimitError(
                :work, required, budget.max_work, "$label would exceed max_work"
            ),
        )
    end
    required <= typemax(Int) ||
        throw(ArgumentError("UPB work accounting exceeds typemax(Int)"))
    budget.work_used = Int(required)
    return nothing
end

function _upbc_record_draws!(budget::_UPBCBudget, amount::Integer)
    required = big(budget.random_draws) + big(amount)
    required <= typemax(Int) ||
        throw(ArgumentError("UPB random-draw accounting exceeds typemax(Int)"))
    budget.random_draws = Int(required)
    return nothing
end

function _upbc_guard_entries(
    amount::Integer, limit::Union{Nothing,Int}, resource::Symbol, label::AbstractString
)
    amount >= 0 || error("internal UPB entry estimate is negative")
    if limit !== nothing && big(amount) > limit
        throw(UPBResourceLimitError(resource, big(amount), limit, label))
    end
    amount <= typemax(Int) || throw(ArgumentError("$label exceeds typemax(Int) entries"))
    return Int(amount)
end

function _upbc_checked_product(values, label::AbstractString)
    value = big(1)
    for factor in values
        value *= factor
    end
    value <= typemax(Int) || throw(ArgumentError("$label exceeds typemax(Int)"))
    return Int(value)
end

function _upbc_dimensions(dimensions; allow_one_party::Bool=true)
    dimensions isa Tuple ||
        dimensions isa AbstractVector ||
        throw(ArgumentError("dimensions must be a tuple or one-based vector"))
    dimensions isa AbstractVector && Base.require_one_based_indexing(dimensions)
    minimum_parties = allow_one_party ? 1 : 2
    length(dimensions) >= minimum_parties ||
        throw(ArgumentError("dimensions must contain at least $minimum_parties entries"))
    checked = Vector{Int}(undef, length(dimensions))
    for (index, dimension) in pairs(dimensions)
        dimension isa Bool &&
            throw(ArgumentError("dimensions[$index] must be an integer, not Bool"))
        dimension isa Integer ||
            throw(ArgumentError("dimensions[$index] must be an integer"))
        dimension >= 2 || throw(ArgumentError("dimensions[$index] must be at least 2"))
        dimension <= typemax(Int) ||
            throw(ArgumentError("dimensions[$index] exceeds typemax(Int)"))
        checked[index] = Int(dimension)
    end
    return Tuple(checked)
end

function _upbc_positive_int(value, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a positive integer, not Bool"))
    value isa Integer || throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be positive"))
    value <= typemax(Int) || throw(ArgumentError("$name exceeds typemax(Int)"))
    return Int(value)
end

function _upbc_isprime(value::Int)
    value >= 2 || return false
    value == 2 && return true
    iseven(value) && return false
    divisor = 3
    while divisor <= value ÷ divisor
        value % divisor == 0 && return false
        divisor += 2
    end
    return true
end

function _upbc_reference(key::Symbol)
    key === :bennett_divincenzo_mor_shor_smolin_terhal_1999 && return (
        "C. H. Bennett, D. P. DiVincenzo, T. Mor, P. W. Shor, " *
        "J. A. Smolin, and B. M. Terhal, Phys. Rev. Lett. 82, 5385–5388 (1999).",
        "https://arxiv.org/abs/quant-ph/9808030",
    )
    key === :divincenzo_mor_shor_smolin_terhal_2003 && return (
        "D. P. DiVincenzo, T. Mor, P. W. Shor, J. A. Smolin, and " *
        "B. M. Terhal, Commun. Math. Phys. 238, 379–410 (2003).",
        "https://arxiv.org/abs/quant-ph/9908070",
    )
    key === :feng_2006 && return (
        "K. Feng, Discrete Appl. Math. 154, 942–949 (2006).",
        "https://doi.org/10.1016/j.dam.2005.10.011",
    )
    key === :alon_lovasz_2001 && return (
        "N. Alon and L. Lovász, J. Combin. Theory Ser. A 95, 169–179 (2001).",
        "https://doi.org/10.1006/jcta.2000.3122",
    )
    key === :chen_johnston_2015 && return (
        "J. Chen and N. Johnston, Commun. Math. Phys. 333, 351–365 (2015).",
        "https://arxiv.org/abs/1301.1406",
    )
    key === :johnston_2013 && return (
        "N. Johnston, Proc. TQC 2013, LIPIcs 22, 93–105 (2013).",
        "https://arxiv.org/abs/1302.1604",
    )
    key === :pedersen_2002 && return (
        "T. B. Pedersen, Characteristics of Unextendible Product Bases, " *
        "Aarhus University thesis (2002).",
        "https://www.cs.au.dk/~tomp/UPB/",
    )
    key === :elementary_product_basis && return (
        "Complete product-basis construction (elementary finite-dimensional " *
        "linear algebra).",
        "",
    )
    return error("internal UPB reference key is unknown: $key")
end

function _upbc_normalize_columns(matrix::AbstractMatrix)
    Base.require_one_based_indexing(matrix)
    result = copy(matrix)
    for column in axes(result, 2)
        column_norm = norm(@view result[:, column])
        isfinite(column_norm) && !iszero(column_norm) || throw(
            ArgumentError("UPB construction produced a zero or nonfinite local vector")
        )
        @views result[:, column] ./= column_norm
    end
    return result
end

function _upbc_rotation(::Type{T}, denominator::Int) where {T<:AbstractFloat}
    angle = convert(T, 2) * T(pi) / convert(T, denominator)
    cosine = cos(angle)
    sine = sin(angle)
    return T[cosine sine; -sine cosine]
end

function _upbc_fourier(::Type{T}, dimension::Int) where {T<:AbstractFloat}
    root_angle = convert(T, 2) * T(pi) / convert(T, dimension)
    root = cis(root_angle)
    scale = inv(sqrt(convert(T, dimension)))
    matrix = Matrix{Complex{T}}(undef, dimension, dimension)
    for column in 0:(dimension - 1), row in 0:(dimension - 1)
        matrix[row + 1, column + 1] = scale * root^(row * column)
    end
    return matrix
end

function _upbc_one_factorization(objects)
    labels = if objects isa Integer
        count = _upbc_positive_int(objects, "objects")
        collect(1:count)
    elseif objects isa Tuple
        collect(objects)
    elseif objects isa AbstractVector
        Base.require_one_based_indexing(objects)
        collect(objects)
    else
        throw(ArgumentError("objects must be an even integer or a collection"))
    end
    count = length(labels)
    count > 0 || throw(ArgumentError("objects must not be empty"))
    iseven(count) ||
        throw(ArgumentError("a one-factorization requires an even number of objects"))
    length(unique(labels)) == count ||
        throw(ArgumentError("one-factorization object labels must be distinct"))

    factorization = Matrix{eltype(labels)}(undef, count - 1, count)
    for round in 1:(count - 1)
        factorization[round, 1] = labels[round]
        factorization[round, 2] = labels[count]
        for pair in 2:(count ÷ 2)
            left = mod(round - (pair - 1) - 1, count - 1) + 1
            right = mod(round + (pair - 1) - 1, count - 1) + 1
            factorization[round, 2pair - 1] = labels[left]
            factorization[round, 2pair] = labels[right]
        end
    end
    return factorization
end

function _upbc_tensor_column(factors::Tuple, column::Int)
    value = copy(@view factors[1][:, column])
    for party in 2:length(factors)
        value = kron(value, @view(factors[party][:, column]))
    end
    return value
end

function _upbc_closed_form_metadata(family::Symbol)
    family in (:pyramid, :tiles) && return :bennett_divincenzo_mor_shor_smolin_terhal_1999
    family in (
        :quad_residue,
        :six_parameter,
        :generalized_shifts,
        :generalized_tiles_1,
        :generalized_tiles_2,
    ) && return :divincenzo_mor_shor_smolin_terhal_2003
    family in
    (:feng_2x2x3, :feng_2x2x5, :feng_2x2x2x2, :feng_4x4, :feng_2x2x2x4, :feng_2x2x2x2x5) &&
        return :feng_2006
    family in (:johnston_2_power_8, :johnston_2_power_4k) && return :johnston_2013
    family in (:chen_johnston_4k1, :chen_johnston_bipartite, :chen_johnston_4x6) &&
        return :chen_johnston_2015
    family === :minimum_4x4 && return :pedersen_2002
    family === :alon_lovasz && return :alon_lovasz_2001
    return error("internal UPB family has no reference: $family")
end

function _upbc_pyramid(::Type{T}) where {T<:AbstractFloat}
    root_five = sqrt(convert(T, 5))
    height = sqrt(one(T) + root_five) / convert(T, 2)
    scale = convert(T, 2) / sqrt(convert(T, 5) + root_five)
    first = Matrix{T}(undef, 3, 5)
    for index in 0:4
        angle = convert(T, 2index) * T(pi) / convert(T, 5)
        first[:, index + 1] = scale .* T[cos(angle), sin(angle), height]
    end
    second = first[:, [1, 3, 5, 2, 4]]
    return (first, second)
end

function _upbc_tiles(::Type{T}) where {T<:AbstractFloat}
    inverse_sqrt_two = inv(sqrt(convert(T, 2)))
    inverse_sqrt_three = inv(sqrt(convert(T, 3)))
    first = zeros(T, 3, 5)
    second = zeros(T, 3, 5)
    first[:, 1] = T[1, 0, 0]
    first[:, 2] = inverse_sqrt_two .* T[1, -1, 0]
    first[:, 3] = T[0, 0, 1]
    first[:, 4] = inverse_sqrt_two .* T[0, 1, -1]
    first[:, 5] = fill(inverse_sqrt_three, 3)
    second[:, 1] = inverse_sqrt_two .* T[1, -1, 0]
    second[:, 2] = T[0, 0, 1]
    second[:, 3] = inverse_sqrt_two .* T[0, 1, -1]
    second[:, 4] = T[1, 0, 0]
    second[:, 5] = fill(inverse_sqrt_three, 3)
    return (first, second)
end

function _upbc_generalized_tiles_1(::Type{T}, dimension::Int) where {T<:AbstractFloat}
    dimension >= 4 && iseven(dimension) || throw(
        ArgumentError(
            "generalized_tiles_1 requires an even local dimension of at least 4; " *
            "the pinned n=2 acceptance produces an extendible singleton and is rejected",
        ),
    )
    cardinality = dimension^2 - 2dimension + 1
    first = zeros(Complex{T}, dimension, cardinality)
    second = zeros(Complex{T}, dimension, cardinality)
    identity_matrix = Matrix{Complex{T}}(I, dimension, dimension)
    root = cis(convert(T, 4) * T(pi) / convert(T, dimension))
    half = dimension ÷ 2
    column = 1
    for mode in 1:(half - 1)
        wave = zeros(Complex{T}, dimension, dimension)
        for offset in 0:(dimension - 1), shift in 0:(half - 1)
            row = mod(shift + offset, dimension) + 1
            wave[row, offset + 1] += root^(shift * mode)
        end
        for offset in 0:(dimension - 1)
            first[:, column] = identity_matrix[:, offset + 1]
            second[:, column] =
                wave[:, mod(offset + 1, dimension) + 1] / sqrt(convert(T, half))
            first[:, column + 1] = wave[:, offset + 1] / sqrt(convert(T, half))
            second[:, column + 1] = identity_matrix[:, offset + 1]
            column += 2
        end
    end
    stopper = inv(sqrt(convert(T, dimension)))
    first[:, cardinality] .= stopper
    second[:, cardinality] .= stopper
    column == cardinality || error("internal generalized_tiles_1 cardinality mismatch")
    return (first, second)
end

function _upbc_generalized_tiles_2(
    ::Type{T}, first_dimension::Int, second_dimension::Int
) where {T<:AbstractFloat}
    first_dimension >= 3 ||
        throw(ArgumentError("generalized_tiles_2 requires first dimension at least 3"))
    second_dimension >= 4 ||
        throw(ArgumentError("generalized_tiles_2 requires second dimension at least 4"))
    second_dimension >= first_dimension || throw(
        ArgumentError(
            "generalized_tiles_2 requires second dimension not smaller than the first"
        ),
    )
    m = first_dimension
    n = second_dimension
    cardinality = m * n - 2m + 1
    first = zeros(Complex{T}, m, cardinality)
    second = zeros(Complex{T}, n, cardinality)
    identity_m = Matrix{Complex{T}}(I, m, m)
    identity_n = Matrix{Complex{T}}(I, n, n)
    inverse_sqrt_two = inv(sqrt(convert(T, 2)))
    root = cis(convert(T, 2) * T(pi) / convert(T, n - 2))

    for column in 1:m
        next_column = mod(column, m) + 1
        first[:, column] =
            inverse_sqrt_two .* (identity_m[:, column] - identity_m[:, next_column])
        second[:, column] = identity_n[:, column]
    end

    column = m + 1
    for offset in 0:(m - 1), mode in 1:(n - 3)
        first[:, column] = identity_m[:, offset + 1]
        for exponent in 0:(m - 3)
            row = mod(exponent + offset + 1, m) + 1
            second[row, column] = root^(exponent * mode)
        end
        for exponent in (m - 2):(n - 3)
            second[exponent + 3, column] = root^(exponent * mode)
        end
        second[:, column] ./= sqrt(convert(T, n - 2))
        column += 1
    end
    stopper_first = inv(sqrt(convert(T, m)))
    stopper_second = inv(sqrt(convert(T, n)))
    first[:, cardinality] .= stopper_first
    second[:, cardinality] .= stopper_second
    column == cardinality || error("internal generalized_tiles_2 cardinality mismatch")
    return (first, second)
end

function _upbc_minimum_4x4(::Type{T}) where {T<:AbstractFloat}
    root_two = sqrt(convert(T, 2))
    first = zeros(T, 4, 8)
    second = zeros(T, 4, 8)
    first[:, 1] = T[1, -3, 1, 1] / sqrt(convert(T, 12))
    first[:, 2] = T[1, 0, 0, 0]
    first[:, 3] = T[0, 1, 2, 1] / sqrt(convert(T, 6))
    first[:, 4] = T[1, 0, 0, -1] / sqrt(convert(T, 2))
    first[:, 5] = T[0, 1, 0, 0]
    first[:, 6] = T[3, 1, -1, 1] / sqrt(convert(T, 12))
    first[:, 7] = T[0, 1, 1, 0] / sqrt(convert(T, 2))
    first[:, 8] = T[0, 0, 1, 0]
    second[:, 1] = T[0, 1, -3 - root_two, -1 - root_two] / sqrt(convert(T, 15) + 8root_two)
    second[:, 2] = T[1, 0, 0, 0]
    second[:, 3] = T[1, 0, root_two - 1, 1] / sqrt(convert(T, 5) - 2root_two)
    second[:, 4] = T[0, 1, 0, 0]
    second[:, 5] = T[-1, 1 + root_two, 0, 1] / sqrt(convert(T, 5) + 2root_two)
    second[:, 6] = T[0, 0, 1, 0]
    second[:, 7] = T[1, 1, 1, -root_two] / sqrt(convert(T, 5))
    second[:, 8] = T[-1, 1 + root_two, 0, 1] / sqrt(convert(T, 5) + 2root_two)
    return (first, second)
end

function _upbc_quad_residue(::Type{T}, dimension::Int) where {T<:AbstractFloat}
    dimension >= 3 && isodd(dimension) ||
        throw(ArgumentError("quad_residue requires an odd dimension of at least 3"))
    prime = 2dimension - 1
    _upbc_isprime(prime) ||
        throw(ArgumentError("quad_residue requires 2*dimension-1 to be prime; got $prime"))
    residues = sort!(unique(mod(index^2, prime) for index in 1:(prime ÷ 2)))
    residue_set = Set(residues)
    nonresidue = findfirst(index -> !(index in residue_set), 1:(prime - 1))
    nonresidue === nothing && error("internal quadratic nonresidue search failed")
    root = cis(convert(T, 2) * T(pi) / convert(T, prime))
    gauss_sum = sum(root^residue for residue in residues)
    gauss_imaginary = abs(imag(gauss_sum))
    gauss_tolerance = convert(T, 128prime) * eps(T)
    gauss_imaginary <= gauss_tolerance || throw(
        ArgumentError(
            "quadratic-residue Gauss sum has unexpected imaginary residual " *
            "$gauss_imaginary",
        ),
    )
    gauss_real = real(gauss_sum)
    normalization = max(-gauss_real, one(T) + gauss_real)
    normalization > zero(T) ||
        error("internal quadratic-residue normalization is not positive")
    fourier = _upbc_fourier(T, prime)
    fourier[1, :] .*= sqrt(normalization)
    first_rows = (1, (residues .+ 1)...)
    second_residues = mod.(nonresidue .* residues, prime)
    second_rows = (1, (second_residues .+ 1)...)
    first = _upbc_normalize_columns(fourier[collect(first_rows), :])
    second = _upbc_normalize_columns(fourier[collect(second_rows), :])
    return (first, second)
end

function _upbc_six_parameter(::Type{T}, parameters) where {T<:AbstractFloat}
    parameters isa Tuple ||
        parameters isa AbstractVector ||
        throw(ArgumentError("six_parameter requires a six-entry tuple or vector"))
    parameters isa AbstractVector && Base.require_one_based_indexing(parameters)
    length(parameters) == 6 ||
        throw(ArgumentError("six_parameter requires exactly six parameters"))
    values = ntuple(
        index -> begin
            value = try
                convert(T, parameters[index])
            catch
                throw(ArgumentError("parameter $index must be representable as $T"))
            end
            isfinite(value) || throw(ArgumentError("parameter $index must be finite"))
            value
        end,
        6,
    )
    gamma_a, theta_a, phi_a, gamma_b, theta_b, phi_b = values
    threshold = convert(T, 10) * eps(T)
    for (name, angle) in (
        ("gamma_a", gamma_a),
        ("theta_a", theta_a),
        ("gamma_b", gamma_b),
        ("theta_b", theta_b),
    )
        abs(sin(angle)) > threshold && abs(cos(angle)) > threshold ||
            throw(ArgumentError("$name must not be within 10 eps of a multiple of pi/2"))
    end
    normalizer_a = sqrt(cos(gamma_a)^2 + sin(gamma_a)^2 * cos(theta_a)^2)
    normalizer_b = sqrt(cos(gamma_b)^2 + sin(gamma_b)^2 * cos(theta_b)^2)
    first = zeros(Complex{T}, 3, 5)
    second = zeros(Complex{T}, 3, 5)
    first[:, 1] = T[1, 0, 0]
    first[:, 2] = T[0, 1, 0]
    first[:, 3] = T[cos(theta_a), 0, sin(theta_a)]
    first[:, 4] = Complex{T}[
        sin(gamma_a) * sin(theta_a), cos(gamma_a) * cis(phi_a), -sin(gamma_a) * cos(theta_a)
    ]
    first[:, 5] =
        Complex{T}[0, sin(gamma_a) * cos(theta_a) * cis(phi_a), cos(gamma_a)] / normalizer_a
    second[:, 1] = T[0, 1, 0]
    second[:, 2] = Complex{T}[
        sin(gamma_b) * sin(theta_b), cos(gamma_b) * cis(phi_b), -sin(gamma_b) * cos(theta_b)
    ]
    second[:, 3] = T[1, 0, 0]
    second[:, 4] = T[cos(theta_b), 0, sin(theta_b)]
    second[:, 5] =
        Complex{T}[0, sin(gamma_b) * cos(theta_b) * cis(phi_b), cos(gamma_b)] / normalizer_b
    return (first, second)
end

function _upbc_generalized_shifts(::Type{T}, parties::Int) where {T<:AbstractFloat}
    parties >= 3 && isodd(parties) ||
        throw(ArgumentError("generalized_shifts requires an odd party count of at least 3"))
    half_cardinality = (parties + 1) ÷ 2
    cardinality = 2half_cardinality
    base = Matrix{T}(undef, 2, cardinality)
    for index in 0:(cardinality - 1)
        angle = convert(T, index) * T(pi) / convert(T, 2half_cardinality)
        base[:, index + 1] = T[cos(angle), sin(angle)]
    end
    order = vcat(
        1, collect((half_cardinality + 1):-1:2), collect((half_cardinality + 2):cardinality)
    )
    base = base[:, order]
    factors = Vector{Matrix{T}}(undef, parties)
    factors[1] = base
    tail = collect(2:cardinality)
    for party in 2:parties
        factors[party] = base[:, vcat(1, circshift(tail, party - 1))]
    end
    return Tuple(factors)
end

function _upbc_feng_2x2x2x2(::Type{T}) where {T<:AbstractFloat}
    basis_1 = Matrix{T}(I, 2, 2)
    basis_2 = T[1 1; 1 -1] / sqrt(convert(T, 2))
    angle = T(pi) / convert(T, 3)
    basis_3 = T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    return (
        hcat(
            basis_1[:, 1],
            basis_1[:, 2],
            basis_1[:, 1],
            basis_2[:, 1],
            basis_2[:, 2],
            basis_2[:, 1],
        ),
        hcat(
            basis_1[:, 1],
            basis_2[:, 1],
            basis_1[:, 2],
            basis_1[:, 2],
            basis_2[:, 2],
            basis_1[:, 1],
        ),
        hcat(
            basis_1[:, 1],
            basis_2[:, 1],
            basis_3[:, 1],
            basis_2[:, 2],
            basis_3[:, 2],
            basis_1[:, 2],
        ),
        hcat(
            basis_1[:, 1],
            basis_2[:, 1],
            basis_3[:, 1],
            basis_3[:, 2],
            basis_1[:, 2],
            basis_2[:, 2],
        ),
    )
end

function _upbc_johnston_2_power_8(::Type{T}) where {T<:AbstractFloat}
    b1 = Matrix{T}(I, 2, 2)
    b2 = T[1 1; 1 -1] / sqrt(convert(T, 2))
    b3 = let angle = T(pi) / convert(T, 3)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    b4 = let angle = T(pi) / convert(T, 5)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    b5 = let angle = T(pi) / convert(T, 7)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    vectors = (
        (
            (b1, 1),
            (b2, 1),
            (b2, 2),
            (b1, 2),
            (b3, 1),
            (b4, 1),
            (b4, 1),
            (b3, 1),
            (b1, 2),
            (b3, 2),
            (b4, 2),
        ),
        (
            (b1, 1),
            (b1, 2),
            (b2, 1),
            (b3, 1),
            (b4, 1),
            (b4, 1),
            (b4, 2),
            (b4, 2),
            (b3, 1),
            (b2, 2),
            (b3, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b3, 1),
            (b3, 2),
            (b1, 2),
            (b3, 2),
            (b1, 2),
            (b4, 1),
            (b3, 1),
            (b2, 2),
            (b4, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b3, 1),
            (b3, 1),
            (b3, 2),
            (b4, 1),
            (b3, 2),
            (b2, 2),
            (b4, 1),
            (b4, 2),
            (b1, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b1, 2),
            (b3, 1),
            (b4, 1),
            (b2, 2),
            (b5, 1),
            (b3, 2),
            (b3, 1),
            (b5, 2),
            (b4, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b3, 1),
            (b2, 2),
            (b4, 1),
            (b4, 2),
            (b5, 1),
            (b5, 2),
            (b2, 2),
            (b1, 2),
            (b3, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b3, 1),
            (b4, 1),
            (b2, 2),
            (b4, 2),
            (b2, 2),
            (b1, 2),
            (b3, 2),
            (b5, 1),
            (b5, 2),
        ),
        (
            (b1, 1),
            (b2, 1),
            (b3, 1),
            (b4, 1),
            (b5, 1),
            (b1, 2),
            (b5, 1),
            (b3, 2),
            (b5, 2),
            (b4, 2),
            (b2, 2),
        ),
    )
    return ntuple(
        party -> hcat((matrix[:, column] for (matrix, column) in vectors[party])...), 8
    )
end

function _upbc_johnston_2_power_4k(::Type{T}, parties::Int) where {T<:AbstractFloat}
    parties >= 8 && parties % 4 == 0 || throw(
        ArgumentError(
            "johnston_2_power_4k requires a party count divisible by 4 and at least 8"
        ),
    )
    k = parties ÷ 4
    block_size = k + 1
    vertex_count = 2block_size
    cardinality = 2vertex_count
    bases = Array{T}(undef, 2, 2, vertex_count)
    for index in 1:vertex_count
        bases[:, :, index] = _upbc_rotation(T, 2index - 1)
    end

    factors = Vector{Matrix{T}}(undef, parties)
    # State order is v₀,w₀,…,vₖ,wₖ,x₀,y₀,…,xₖ,yₖ. The first three
    # orthogonality graphs are B_{0,k}, B_{1,k}, and B_{2,k}.
    for graph in 0:2
        factor = Matrix{T}(undef, 2, cardinality)
        for index in 0:k
            v = 2index + 1
            w = v + 1
            x = vertex_count + 2index + 1
            y = x + 1
            shifted = mod(index - graph, block_size) + 1
            factor[:, v] = bases[:, 1, index + 1]
            factor[:, w] = bases[:, 1, index + 1]
            factor[:, x] = bases[:, 2, shifted]
            factor[:, y] = bases[:, 2, shifted]
        end
        factors[graph + 1] = factor
    end

    # Each remaining B_{j,k}, 3 ≤ j ≤ k, is represented on two parties.
    next_party = 4
    for graph in 3:k
        first_factor = Matrix{T}(undef, 2, cardinality)
        second_factor = Matrix{T}(undef, 2, cardinality)
        for index in 0:k
            v = 2index + 1
            w = v + 1
            x = vertex_count + 2index + 1
            y = x + 1
            shifted = mod(index - graph, block_size) + 1
            first_factor[:, v] = bases[:, 1, index + 1]
            first_factor[:, w] = bases[:, 1, block_size + index + 1]
            first_factor[:, x] = bases[:, 2, shifted]
            first_factor[:, y] = bases[:, 2, block_size + shifted]
            second_factor[:, v] = bases[:, 1, index + 1]
            second_factor[:, w] = bases[:, 1, block_size + index + 1]
            second_factor[:, x] = bases[:, 2, block_size + shifted]
            second_factor[:, y] = bases[:, 2, shifted]
        end
        factors[next_party] = first_factor
        factors[next_party + 1] = second_factor
        next_party += 2
    end

    # A one-factorization of each of the two K_{2k+2} components covers the
    # complement of the complete bipartite graph. Distinct local bases are
    # used for the two components so every local vector is repeated only
    # within its matched edge, as required by Johnston's unextendibility
    # count.
    factorization = _upbc_one_factorization(vertex_count)
    for round in 1:(vertex_count - 1)
        factor = Matrix{T}(undef, 2, cardinality)
        for pair in 1:block_size
            first_left = factorization[round, 2pair - 1]
            first_right = factorization[round, 2pair]
            factor[:, first_left] = bases[:, 1, pair]
            factor[:, first_right] = bases[:, 2, pair]
            second_left = vertex_count + first_left
            second_right = vertex_count + first_right
            factor[:, second_left] = bases[:, 1, block_size + pair]
            factor[:, second_right] = bases[:, 2, block_size + pair]
        end
        factors[next_party] = factor
        next_party += 1
    end
    next_party == parties + 1 || error("internal Johnston party count is inconsistent")
    all(isassigned(factors, party) for party in eachindex(factors)) ||
        error("internal johnston_2_power_4k party construction is incomplete")
    return Tuple(factors)
end

function _upbc_chen_johnston_4x6(::Type{T}) where {T<:AbstractFloat}
    root = convert(T, 1.64451358502312496885542269243)
    first = _upbc_normalize_columns(
        T[
            3 3 3 1 1 -1 -1 2 -2 0
            2 1 1 2 3 -2 0 0 -2 -1
            2 1 1 1 1 5 -3 -4 2 1
            2 2 3 2 2 0 2 1 3 0
        ],
    )
    tail = T[
        0 0 (root-2)/(1-root) 1 (root-1)/(root-2)
        root*(root-2)/(2root-3) 0 0 (3-2root)/(root-2) 1
        1 -1-root 0 0 1/(1+root)
        1 1 -root 0 0
        0 root-1 1 1/(1-root) 0
        root 1 1 1 1
    ]
    second = hcat(Matrix{T}(I, 6, 6)[:, 1:5], _upbc_normalize_columns(tail))
    return (first, second)
end

function _upbc_feng_2x2x3(::Type{T}) where {T<:AbstractFloat}
    b1 = Matrix{T}(I, 2, 2)
    b2 = T[1 1; 1 -1] / sqrt(convert(T, 2))
    third = _upbc_normalize_columns(T[1 1 2 1 0 0; 0 0 -3 1 1 1; 0 1 -1 -1 -3 0])
    first = hcat(b1[:, 1], b1[:, 2], b1[:, 1], b2[:, 1], b2[:, 2], b2[:, 1])
    second = hcat(b1[:, 1], b2[:, 1], b1[:, 2], b1[:, 2], b2[:, 2], b1[:, 1])
    return (first, second, third)
end

function _upbc_feng_2x2x5(::Type{T}) where {T<:AbstractFloat}
    root = cis(convert(T, 2) * T(pi) / convert(T, 3))
    b1 = Matrix{T}(I, 2, 2)
    b2 = T[1 1; 1 -1] / sqrt(convert(T, 2))
    b3 = let angle = T(pi) / convert(T, 3)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    b4 = let angle = T(pi) / convert(T, 5)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    third = _upbc_normalize_columns(
        Complex{T}[
            1 conj(root) root 0 0 0 1 0
            0 root conj(root) 1 0 1 0 0
            0 conj(root) 0 0 0 1 root 1
            0 0 root 0 1 1 conj(root) 0
            0 1 1 0 0 1 1 0
        ],
    )
    first = hcat(
        b2[:, 2], b2[:, 1], b1[:, 2], b1[:, 1], b1[:, 1], b1[:, 2], b2[:, 1], b2[:, 2]
    )
    second = hcat(
        b4[:, 2], b1[:, 2], b4[:, 1], b1[:, 1], b3[:, 1], b2[:, 1], b3[:, 2], b2[:, 2]
    )
    return (first, second, third)
end

function _upbc_feng_4x4(::Type{T}) where {T<:AbstractFloat}
    inverse_root_three = inv(sqrt(convert(T, 3)))
    first = zeros(T, 4, 8)
    second = zeros(T, 4, 8)
    first[:, 1:4] = Matrix{T}(I, 4, 4)
    first[:, 5] = inverse_root_three .* T[0, 1, 1, 1]
    first[:, 6] = inverse_root_three .* T[1, 0, -1, 1]
    first[:, 7] = inverse_root_three .* T[1, 1, 0, -1]
    first[:, 8] = inverse_root_three .* T[1, -1, 1, 0]
    second[:, [1, 7, 6, 4]] = Matrix{T}(I, 4, 4)
    second[:, 2] = inverse_root_three .* T[1, 0, -1, 1]
    second[:, 3] = inverse_root_three .* T[1, 1, 0, -1]
    second[:, 5] = inverse_root_three .* T[1, -1, 1, 0]
    second[:, 8] = inverse_root_three .* T[0, 1, 1, 1]
    return (first, second)
end

function _upbc_feng_2x2x2x4(::Type{T}) where {T<:AbstractFloat}
    b1 = Matrix{T}(I, 2, 2)
    b2 = T[1 1; 1 -1] / sqrt(convert(T, 2))
    b3 = let angle = T(pi) / convert(T, 3)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    b4 = let angle = T(pi) / convert(T, 5)
        T[cos(angle) sin(angle); -sin(angle) cos(angle)]
    end
    bases = hcat(b1, b2, b3, b4)
    fourth = zeros(T, 4, 8)
    fourth[:, 1:4] = Matrix{T}(I, 4, 4)
    inverse_root_three = inv(sqrt(convert(T, 3)))
    fourth[:, 5] = inverse_root_three .* T[0, 1, 1, 1]
    fourth[:, 6] = inverse_root_three .* T[1, 0, -1, 1]
    fourth[:, 7] = inverse_root_three .* T[1, 1, 0, -1]
    fourth[:, 8] = inverse_root_three .* T[1, -1, 1, 0]
    return (
        bases[:, [1, 3, 7, 5, 8, 6, 2, 4]],
        bases[:, [1, 3, 7, 5, 4, 2, 6, 8]],
        bases[:, [1, 5, 7, 3, 4, 8, 6, 2]],
        fourth,
    )
end

function _upbc_feng_2x2x2x2x5(::Type{T}) where {T<:AbstractFloat}
    root = cis(convert(T, 2) * T(pi) / convert(T, 3))
    fifth = zeros(Complex{T}, 5, 10)
    fifth[:, 1:5] = Matrix{Complex{T}}(I, 5, 5)
    fifth[:, 6] = Complex{T}[0, 1, 1, 1, 1] / convert(T, 2)
    fifth[:, 7] = Complex{T}[1, 0, 1, root, conj(root)] / convert(T, 2)
    fifth[:, 8] = Complex{T}[1, 1, 0, conj(root), root] / convert(T, 2)
    fifth[:, 9] = Complex{T}[1, root, conj(root), 0, 1] / convert(T, 2)
    fifth[:, 10] = Complex{T}[1, conj(root), root, 1, 0] / convert(T, 2)
    bases = hcat(
        Matrix{T}(I, 2, 2),
        T[1 1; 1 -1] / sqrt(convert(T, 2)),
        let angle = T(pi) / convert(T, 3)
            T[cos(angle) sin(angle); -sin(angle) cos(angle)]
        end,
        let angle = T(pi) / convert(T, 5)
            T[cos(angle) sin(angle); -sin(angle) cos(angle)]
        end,
        let angle = T(pi) / convert(T, 7)
            T[cos(angle) sin(angle); -sin(angle) cos(angle)]
        end,
    )
    return (
        bases[:, [1, 3, 5, 7, 9, 4, 6, 8, 10, 2]],
        bases[:, [1, 3, 5, 7, 9, 10, 2, 4, 6, 8]],
        bases[:, [1, 3, 5, 7, 9, 8, 10, 2, 4, 6]],
        bases[:, [1, 3, 5, 7, 9, 6, 8, 10, 2, 4]],
        fifth,
    )
end

function _upbc_require_randomized_type(::Type{T}) where {T<:AbstractFloat}
    T in (Float32, Float64) || throw(
        ArgumentError(
            "randomized UPB constructions currently require real_type=Float32 or " *
            "Float64 because their rank decisions use LAPACK SVD; deterministic " *
            "families support BigFloat",
        ),
    )
    return nothing
end

function _upbc_random_normal(
    rng::AbstractRNG, ::Type{T}, rows::Int, columns::Int, budget::_UPBCBudget
) where {T<:AbstractFloat}
    _upbc_record_draws!(budget, rows * columns)
    _upbc_consume!(budget, rows * columns, "Gaussian sampling")
    return randn(rng, T, rows, columns)
end

function _upbc_numerical_nullspace(
    matrix::AbstractMatrix{T}, rtol::T, budget::_UPBCBudget
) where {T<:Union{Float32,Float64}}
    rows, columns = size(matrix)
    _upbc_consume!(
        budget, max(rows, 1) * max(columns, 1) * max(min(rows, columns), 1), "nullspace SVD"
    )
    decomposition = svd(matrix; full=true)
    largest = isempty(decomposition.S) ? zero(T) : maximum(decomposition.S)
    threshold = rtol * max(largest, one(T)) * convert(T, max(rows, columns, 1))
    numerical_rank = count(value -> value > threshold, decomposition.S)
    return decomposition.V[:, (numerical_rank + 1):end]
end

function _upbc_advance_combination!(indices::Vector{Int}, count::Int)
    width = length(indices)
    position = width
    while position >= 1 && indices[position] == count - width + position
        position -= 1
    end
    position == 0 && return false
    indices[position] += 1
    for next_position in (position + 1):width
        indices[next_position] = indices[next_position - 1] + 1
    end
    return true
end

function _upbc_full_spark(
    matrix::AbstractMatrix{T},
    subset_size::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    rows, columns = size(matrix)
    subset_size <= rows && subset_size <= columns ||
        throw(ArgumentError("full-spark subset size exceeds matrix dimensions"))
    row_combinations = binomial(big(rows), subset_size)
    column_combinations = binomial(big(columns), subset_size)
    required_minors = row_combinations * column_combinations
    if budget.max_minors !== nothing &&
        big(budget.minors_checked) + required_minors > budget.max_minors
        throw(
            UPBResourceLimitError(
                :minors,
                big(budget.minors_checked) + required_minors,
                budget.max_minors,
                "full-spark verification would exceed max_minors",
            ),
        )
    end
    required_minors <= typemax(Int) ||
        throw(ArgumentError("full-spark minor count exceeds typemax(Int)"))

    row_indices = collect(1:subset_size)
    while true
        column_indices = collect(1:subset_size)
        while true
            _upbc_consume!(budget, subset_size^3, "full-spark singular-value verification")
            singular_values = svdvals(Matrix(@view matrix[row_indices, column_indices]))
            largest = maximum(singular_values)
            threshold = rtol * max(largest, one(T)) * convert(T, max(subset_size, 1))
            budget.minors_checked += 1
            minimum(singular_values) > boundary_factor * threshold || return false
            _upbc_advance_combination!(column_indices, columns) || break
        end
        _upbc_advance_combination!(row_indices, rows) || break
    end
    return true
end

function _upbc_alon_lovasz_odd(
    rng::AbstractRNG,
    ::Type{T},
    rows::Int,
    columns::Int,
    offset::Int,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    half = (rows - 1) ÷ 2
    for attempt in 1:max_attempts
        matrix = _upbc_random_normal(rng, T, rows, columns, budget)
        viable = true
        for column in 2:columns
            candidates = vcat(
                collect((column - half - offset):(column - 1 - offset)),
                collect((column + 1 + offset):(column + half + offset)),
            )
            wrapped = unique(mod.(candidates .- 1, columns) .+ 1)
            indices = sort!(intersect(collect(1:(column - 1)), wrapped))
            null_basis = _upbc_numerical_nullspace(
                transpose(@view(matrix[:, indices])), rtol, budget
            )
            if size(null_basis, 2) == 0
                viable = false
                break
            end
            coefficients = _upbc_random_normal(rng, T, size(null_basis, 2), 1, budget)
            matrix[:, column] = null_basis * vec(coefficients)
        end
        viable || continue
        matrix = _upbc_normalize_columns(matrix)
        _upbc_full_spark(matrix, rows, rtol, boundary_factor, budget) &&
            return matrix, attempt
    end
    return throw(
        UPBConstructionUnavailable(
            :randomized_search_exhausted,
            (rows, columns),
            nothing,
            "bounded Alon–Lovász odd-dimensional full-spark search exhausted " *
            "$max_attempts attempts",
        ),
    )
end

function _upbc_alon_lovasz_even(
    rng::AbstractRNG,
    ::Type{T},
    rows::Int,
    columns::Int,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    complement_width = columns - rows + 1
    isodd(complement_width) ||
        error("internal Alon–Lovász even construction has an invalid graph width")
    half = (complement_width - 1) ÷ 2
    for attempt in 1:max_attempts
        matrix = _upbc_random_normal(rng, T, rows, columns, budget)
        viable = true
        for column in 2:columns
            window = unique(
                mod.(collect((column - half):(column + half)) .- 1, columns) .+ 1
            )
            indices = sort!(setdiff(collect(1:columns), window))
            null_basis = _upbc_numerical_nullspace(
                transpose(@view(matrix[:, indices])), rtol, budget
            )
            if size(null_basis, 2) == 0
                viable = false
                break
            end
            coefficients = _upbc_random_normal(rng, T, size(null_basis, 2), 1, budget)
            matrix[:, column] = null_basis * vec(coefficients)
        end
        viable || continue
        matrix = _upbc_normalize_columns(matrix)
        _upbc_full_spark(matrix, rows, rtol, boundary_factor, budget) &&
            return matrix, attempt
    end
    return throw(
        UPBConstructionUnavailable(
            :randomized_search_exhausted,
            (rows, columns),
            nothing,
            "bounded Alon–Lovász even-dimensional full-spark search exhausted " *
            "$max_attempts attempts",
        ),
    )
end

function _upbc_alon_lovasz(
    rng::AbstractRNG,
    ::Type{T},
    dimensions::Tuple,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    parties = length(dimensions)
    cardinality = sum(dimensions) - parties + 1
    even_count = count(iseven, dimensions)
    ((isodd(sum(dimensions) - parties) || even_count == 0) && even_count <= 1) || throw(
        ArgumentError(
            "dimensions $dimensions do not satisfy the Alon–Lovász " *
            "minimum-cardinality construction condition",
        ),
    )
    factors = Vector{Matrix{T}}(undef, parties)
    offset = 0
    total_attempts = 0

    first_dimension = dimensions[1]
    if iseven(first_dimension)
        factors[1], attempts = _upbc_alon_lovasz_even(
            rng,
            T,
            first_dimension,
            cardinality,
            max_attempts,
            rtol,
            boundary_factor,
            budget,
        )
    else
        factors[1], attempts = _upbc_alon_lovasz_odd(
            rng,
            T,
            first_dimension,
            cardinality,
            offset,
            max_attempts,
            rtol,
            boundary_factor,
            budget,
        )
        offset += (first_dimension - 1) ÷ 2
    end
    total_attempts += attempts

    for party in parties:-1:2
        dimension = dimensions[party]
        if iseven(dimension)
            factors[party], attempts = _upbc_alon_lovasz_even(
                rng, T, dimension, cardinality, max_attempts, rtol, boundary_factor, budget
            )
        else
            factors[party], attempts = _upbc_alon_lovasz_odd(
                rng,
                T,
                dimension,
                cardinality,
                offset,
                max_attempts,
                rtol,
                boundary_factor,
                budget,
            )
            offset += (dimension - 1) ÷ 2
        end
        total_attempts += attempts
    end
    return Tuple(factors), total_attempts
end

function _upbc_hollow_unitary(::Type{T}, dimension::Int) where {T<:AbstractFloat}
    dimension == 2 ||
        dimension >= 4 ||
        throw(ArgumentError("a hollow unitary exists only in dimension 2 or at least 4"))
    root = cis(convert(T, 2) * T(pi) / convert(T, dimension - 2))
    fourier = _upbc_fourier(T, dimension - 1)
    lower = zeros(Complex{T}, dimension - 1, dimension - 1)
    for index in 1:(dimension - 2)
        column = @view fourier[:, index + 1]
        lower .+= root^index .* (column * adjoint(column))
    end
    unitary = zeros(Complex{T}, dimension, dimension)
    edge = inv(sqrt(convert(T, dimension - 1)))
    unitary[1, 2:end] .= edge
    unitary[2:end, 1] .= edge
    unitary[2:end, 2:end] = lower
    return unitary
end

function _upbc_chen_johnston_lemma5(
    rng::AbstractRNG,
    ::Type{T},
    q::Int,
    r::Int,
    s::Int,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    q >= r + s || throw(ArgumentError("Chen–Johnston Lemma 5 requires q >= r+s"))
    rows = r + 1
    for attempt in 1:max_attempts
        matrix = _upbc_random_normal(rng, T, rows, 2q, budget)
        matrix[:, (q + 1):(2q)] = _upbc_normalize_columns(matrix[:, (q + 1):(2q)])
        viable = true
        for index in 0:(q - 1)
            selected = q .+ (mod.(index .+ collect(s:(s + r - 1)), q) .+ 1)
            null_basis = _upbc_numerical_nullspace(
                transpose(@view(matrix[:, selected])), rtol, budget
            )
            if size(null_basis, 2) != 1
                viable = false
                break
            end
            matrix[:, index + 1] = null_basis[:, 1]
        end
        viable || continue
        matrix = reverse(matrix; dims=2)
        _upbc_full_spark(matrix, rows, rtol, boundary_factor, budget) &&
            return matrix, attempt
    end
    return throw(
        UPBConstructionUnavailable(
            :randomized_search_exhausted,
            (q, r + 1),
            2q,
            "bounded Chen–Johnston Lemma 5 full-spark search exhausted " *
            "$max_attempts attempts",
        ),
    )
end

function _upbc_chen_johnston_bipartite(
    rng::AbstractRNG,
    ::Type{T},
    dimensions::Tuple,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    parties = length(dimensions)
    largest = dimensions[end]
    largest - 1 == sum(dimension - 1 for dimension in dimensions[1:(end - 1)]) || throw(
        ArgumentError(
            "Chen–Johnston bipartite-style construction requires " *
            "d_max-1 == sum(d_j-1) over the other parties",
        ),
    )
    sum(dimension - 1 for dimension in dimensions) >= 3 ||
        throw(ArgumentError("Chen–Johnston construction requires lower-bound excess >= 3"))
    iseven(sum(dimensions) - parties) ||
        throw(ArgumentError("Chen–Johnston construction requires even sum(d_j-1)"))

    factors = Vector{AbstractMatrix}(undef, parties)
    hollow = _upbc_hollow_unitary(T, largest)
    factors[parties] = hcat(Matrix{Complex{T}}(I, largest, largest), hollow)
    offset = 0
    total_attempts = 0
    for party in 1:(parties - 1)
        factor, attempts = _upbc_chen_johnston_lemma5(
            rng,
            T,
            largest,
            dimensions[party] - 1,
            1 + offset,
            max_attempts,
            rtol,
            boundary_factor,
            budget,
        )
        factors[party] = factor
        offset += dimensions[party] - 1
        total_attempts += attempts
    end
    return Tuple(factors), total_attempts
end

function _upbc_chen_johnston_lemma6(
    rng::AbstractRNG,
    ::Type{T},
    k::Int,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    dimension = 4k + 1
    cardinality = dimension + 3
    half = cardinality ÷ 2
    for attempt in 1:max_attempts
        matrix = _upbc_random_normal(rng, T, dimension, cardinality, budget)
        matrix[:, 1:2] = _upbc_normalize_columns(matrix[:, 1:2])
        viable = true
        for index in 2:(dimension + 2)
            zero_based = if index < half
                setdiff(collect(0:(index - 2)), (mod(index + 1, half),))
            else
                setdiff(
                    collect(0:(index - 1)),
                    (
                        index - half,
                        half + mod(index - 1, half),
                        half + mod(index + 1, half),
                    ),
                )
            end
            selected = zero_based .+ 1
            null_basis = _upbc_numerical_nullspace(
                transpose(@view(matrix[:, selected])), rtol, budget
            )
            if size(null_basis, 2) == 0
                viable = false
                break
            end
            matrix[:, index + 1] = null_basis[:, 1]
        end
        viable || continue
        matrix = _upbc_normalize_columns(matrix)
        _upbc_full_spark(matrix, dimension, rtol, boundary_factor, budget) &&
            return matrix, attempt
    end
    return throw(
        UPBConstructionUnavailable(
            :randomized_search_exhausted,
            (2, 2, dimension),
            cardinality,
            "bounded Chen–Johnston Lemma 6 full-spark search exhausted " *
            "$max_attempts attempts",
        ),
    )
end

function _upbc_chen_johnston_4k1(
    rng::AbstractRNG,
    ::Type{T},
    dimension::Int,
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:Union{Float32,Float64}}
    dimension >= 5 && dimension % 4 == 1 ||
        throw(ArgumentError("chen_johnston_4k1 requires dimension == 1 mod 4 and >= 5"))
    basis_count = (dimension + 3) ÷ 2
    cardinality = dimension + 3
    bases = Array{T}(undef, 2, 2, basis_count)
    for index in 1:basis_count
        bases[:, :, index] = _upbc_rotation(T, 2index - 1)
    end
    selected_basis_count = (dimension + 3) ÷ 4
    flattened_selected = reshape(bases[:, :, 1:selected_basis_count], 2, basis_count)
    first = hcat(flattened_selected, flattened_selected)
    pair_swap = collect(1:basis_count)
    for index in 1:2:basis_count
        pair_swap[index], pair_swap[index + 1] = pair_swap[index + 1], pair_swap[index]
    end
    first[:, (basis_count + 1):cardinality] = first[:, basis_count .+ pair_swap]

    flattened = reshape(bases, 2, cardinality)
    second = circshift(flattened, (0, 1))
    second[:, [basis_count, cardinality]] = second[:, [cardinality, basis_count]]
    third, attempts = _upbc_chen_johnston_lemma6(
        rng, T, (dimension - 1) ÷ 4, max_attempts, rtol, boundary_factor, budget
    )
    return (first, second, third), attempts
end

struct _UPBCRequest{A,D,P,R,N}
    family::Symbol
    arguments::A
    sorted_dimensions::D
    output_permutation::P
    cardinality::Int
    randomized::Bool
    requested::R
    notes::N
end

function _upbc_canonical_name(name)
    name isa Symbol ||
        name isa AbstractString ||
        throw(ArgumentError("UPB family name must be a Symbol or string"))
    normalized = lowercase(String(name))
    normalized = replace(normalized, "_" => "", "-" => "", " " => "")
    aliases = Dict(
        "pyramid" => :pyramid,
        "tiles" => :tiles,
        "gentiles1" => :generalized_tiles_1,
        "generalizedtiles1" => :generalized_tiles_1,
        "gentiles2" => :generalized_tiles_2,
        "generalizedtiles2" => :generalized_tiles_2,
        "min4x4" => :minimum_4x4,
        "minimum4x4" => :minimum_4x4,
        "quadres" => :quad_residue,
        "quadresidue" => :quad_residue,
        "quadraticresidue" => :quad_residue,
        "sixparam" => :six_parameter,
        "sixparameter" => :six_parameter,
        "shifts" => :shifts,
        "genshifts" => :generalized_shifts,
        "generalizedshifts" => :generalized_shifts,
        "feng2x2x2x2" => :feng_2x2x2x2,
        "john2^8" => :johnston_2_power_8,
        "john2power8" => :johnston_2_power_8,
        "johnston2power8" => :johnston_2_power_8,
        "john2^4k" => :johnston_2_power_4k,
        "john2power4k" => :johnston_2_power_4k,
        "johnston2power4k" => :johnston_2_power_4k,
        "cj4k1" => :chen_johnston_4k1,
        "chenjohnston4k1" => :chen_johnston_4k1,
        "cjbip" => :chen_johnston_bipartite,
        "chenjohnstonbipartite" => :chen_johnston_bipartite,
        "cjbip46" => :chen_johnston_4x6,
        "chenjohnston4x6" => :chen_johnston_4x6,
        "feng2x2x3" => :feng_2x2x3,
        "feng2x2x5" => :feng_2x2x5,
        "feng4x4" => :feng_4x4,
        "feng2x2x2x4" => :feng_2x2x2x4,
        "feng2x2x2x2x5" => :feng_2x2x2x2x5,
        "alonlovasz" => :alon_lovasz,
        "alonlovász" => :alon_lovasz,
    )
    haskey(aliases, normalized) || throw(
        ArgumentError(
            "unknown UPB family $(repr(name)); see the `upb` documentation for " *
            "the executable catalog",
        ),
    )
    return aliases[normalized]
end

function _upbc_no_arguments(arguments, family::Symbol)
    isempty(arguments) ||
        throw(ArgumentError("$family does not accept positional family parameters"))
    return nothing
end

function _upbc_named_request(name, arguments...)
    family = _upbc_canonical_name(name)
    identity_permutation(dimensions) = ntuple(identity, length(dimensions))
    if family === :pyramid
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (3, 3), 5, ()
    elseif family === :tiles
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (3, 3), 5, ()
    elseif family === :generalized_tiles_1
        length(arguments) == 1 ||
            throw(ArgumentError("generalized_tiles_1 requires one local dimension"))
        dimension = _upbc_positive_int(arguments[1], "dimension")
        dimension >= 4 && iseven(dimension) || throw(
            ArgumentError(
                "generalized_tiles_1 requires an even local dimension of at least 4"
            ),
        )
        dimensions = (dimension, dimension)
        cardinality = _upbc_checked_product(
            (dimension - 1, dimension - 1), "generalized_tiles_1 cardinality"
        )
        parsed = (dimension,)
    elseif family === :generalized_tiles_2
        if length(arguments) == 1 &&
            (arguments[1] isa Tuple || arguments[1] isa AbstractVector)
            dimensions = _upbc_dimensions(arguments[1])
            length(dimensions) == 2 ||
                throw(ArgumentError("generalized_tiles_2 dimensions must have length 2"))
        elseif length(arguments) == 2
            dimensions = (
                _upbc_positive_int(arguments[1], "first_dimension"),
                _upbc_positive_int(arguments[2], "second_dimension"),
            )
        else
            throw(
                ArgumentError(
                    "generalized_tiles_2 requires two dimensions or one two-entry collection",
                ),
            )
        end
        dimensions[1] >= 3 && dimensions[2] >= 4 && dimensions[2] >= dimensions[1] ||
            throw(ArgumentError("generalized_tiles_2 requires n >= m, m >= 3, and n >= 4"))
        cardinality_big = big(dimensions[1]) * dimensions[2] - 2dimensions[1] + 1
        cardinality_big <= typemax(Int) ||
            throw(ArgumentError("generalized_tiles_2 cardinality exceeds typemax(Int)"))
        cardinality = Int(cardinality_big)
        parsed = dimensions
    elseif family === :minimum_4x4
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (4, 4), 8, ()
    elseif family === :quad_residue
        length(arguments) == 1 ||
            throw(ArgumentError("quad_residue requires one local dimension"))
        dimension = _upbc_positive_int(arguments[1], "dimension")
        dimension >= 3 && isodd(dimension) ||
            throw(ArgumentError("quad_residue requires an odd dimension of at least 3"))
        prime_big = 2big(dimension) - 1
        prime_big <= typemax(Int) ||
            throw(ArgumentError("2*dimension-1 exceeds typemax(Int)"))
        _upbc_isprime(Int(prime_big)) ||
            throw(ArgumentError("quad_residue requires 2*dimension-1 to be prime"))
        dimensions, cardinality, parsed = (dimension, dimension),
        Int(prime_big),
        (dimension,)
    elseif family === :six_parameter
        length(arguments) == 1 || throw(
            ArgumentError("six_parameter requires one six-entry parameter collection")
        )
        dimensions, cardinality, parsed = (3, 3), 5, (arguments[1],)
    elseif family === :shifts
        _upbc_no_arguments(arguments, family)
        family = :generalized_shifts
        dimensions, cardinality, parsed = (2, 2, 2), 4, (3,)
    elseif family === :generalized_shifts
        length(arguments) == 1 ||
            throw(ArgumentError("generalized_shifts requires one party count"))
        parties = _upbc_positive_int(arguments[1], "parties")
        parties >= 3 && isodd(parties) ||
            throw(ArgumentError("generalized_shifts requires odd parties >= 3"))
        dimensions = ntuple(_ -> 2, parties)
        cardinality, parsed = parties + 1, (parties,)
    elseif family === :feng_2x2x2x2
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (2, 2, 2, 2), 6, ()
    elseif family === :johnston_2_power_8
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = ntuple(_ -> 2, 8), 11, ()
    elseif family === :johnston_2_power_4k
        length(arguments) == 1 ||
            throw(ArgumentError("johnston_2_power_4k requires one party count"))
        parties = _upbc_positive_int(arguments[1], "parties")
        parties >= 8 && parties % 4 == 0 ||
            throw(ArgumentError("johnston_2_power_4k requires parties >= 8 and 0 mod 4"))
        dimensions = ntuple(_ -> 2, parties)
        cardinality, parsed = parties + 4, (parties,)
    elseif family === :chen_johnston_4k1
        length(arguments) == 1 ||
            throw(ArgumentError("chen_johnston_4k1 requires one odd local dimension"))
        dimension = _upbc_positive_int(arguments[1], "dimension")
        dimension >= 5 && dimension % 4 == 1 ||
            throw(ArgumentError("chen_johnston_4k1 requires dimension >= 5 and 1 mod 4"))
        dimensions, cardinality, parsed = (2, 2, dimension), dimension + 3, (dimension,)
    elseif family === :chen_johnston_bipartite
        length(arguments) == 1 || throw(
            ArgumentError("chen_johnston_bipartite requires one dimension collection")
        )
        dimensions = _upbc_dimensions(arguments[1])
        sorted = Tuple(sort(collect(dimensions)))
        sorted[end] - 1 == sum(dimension - 1 for dimension in sorted[1:(end - 1)]) || throw(
            ArgumentError(
                "chen_johnston_bipartite dimensions do not satisfy the " *
                "largest-dimension equality",
            ),
        )
        dimensions = sorted
        cardinality = 2dimensions[end]
        parsed = (dimensions,)
    elseif family === :chen_johnston_4x6
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (4, 6), 10, ()
    elseif family === :feng_2x2x3
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (2, 2, 3), 6, ()
    elseif family === :feng_2x2x5
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (2, 2, 5), 8, ()
    elseif family === :feng_4x4
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (4, 4), 8, ()
    elseif family === :feng_2x2x2x4
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (2, 2, 2, 4), 8, ()
    elseif family === :feng_2x2x2x2x5
        _upbc_no_arguments(arguments, family)
        dimensions, cardinality, parsed = (2, 2, 2, 2, 5), 10, ()
    elseif family === :alon_lovasz
        length(arguments) == 1 ||
            throw(ArgumentError("alon_lovasz requires one dimension collection"))
        dimensions = Tuple(sort(collect(_upbc_dimensions(arguments[1]))))
        parties = length(dimensions)
        cardinality = sum(dimensions) - parties + 1
        even_count = count(iseven, dimensions)
        ((isodd(sum(dimensions) - parties) || even_count == 0) && even_count <= 1) || throw(
            ArgumentError(
                "dimensions do not satisfy the Alon–Lovász construction condition"
            ),
        )
        parsed = (dimensions,)
    else
        error("internal named UPB dispatch is incomplete for $family")
    end
    randomized = family in (:alon_lovasz, :chen_johnston_bipartite, :chen_johnston_4k1)
    notes = if family === :generalized_tiles_1
        ("The mathematically invalid n=2 branch accepted by the pinned source is rejected.",)
    elseif randomized
        (
            "The upstream unbounded/global-RNG search is replaced by an explicit-RNG bounded full-spark search.",
        )
    else
        ()
    end
    return _UPBCRequest(
        family,
        parsed,
        dimensions,
        identity_permutation(dimensions),
        cardinality,
        randomized,
        name,
        notes,
    )
end

function _upbc_dimension_request(dimensions)
    requested_dimensions = _upbc_dimensions(dimensions)
    permutation = sortperm(collect(requested_dimensions); alg=MergeSort)
    sorted_dimensions = Tuple(requested_dimensions[index] for index in permutation)
    output_permutation = Tuple(invperm(permutation))
    parties = length(sorted_dimensions)
    total_excess = sum(dimension - 1 for dimension in sorted_dimensions)

    family = nothing
    arguments = ()
    cardinality = 0
    randomized = false
    notes = ()
    if parties == 1
        family = :complete_product_basis
        cardinality = sorted_dimensions[1]
    elseif parties == 2 && sorted_dimensions[1] <= 2
        family = :complete_product_basis
        cardinality = _upbc_checked_product(sorted_dimensions, "complete basis cardinality")
    elseif parties == 2 && sorted_dimensions == (3, 3)
        family, cardinality = :tiles, 5
    elseif parties == 2 && sorted_dimensions == (4, 4)
        family, cardinality = :feng_4x4, 8
    elseif parties == 2 &&
        sorted_dimensions[1] == sorted_dimensions[2] &&
        isodd(sorted_dimensions[1]) &&
        2big(sorted_dimensions[1]) - 1 <= typemax(Int) &&
        _upbc_isprime(2sorted_dimensions[1] - 1)
        family = :quad_residue
        arguments = (sorted_dimensions[1],)
        cardinality = 2sorted_dimensions[1] - 1
    elseif parties == 3 && sorted_dimensions == (2, 2, 2)
        family, arguments, cardinality = :generalized_shifts, (3,), 4
    elseif parties == 3 && sorted_dimensions == (2, 2, 3)
        family, cardinality = :feng_2x2x3, 6
    elseif parties == 3 && sorted_dimensions == (2, 2, 5)
        family, cardinality = :feng_2x2x5, 8
    elseif parties == 4 && sorted_dimensions == (2, 2, 2, 2)
        family, cardinality = :feng_2x2x2x2, 6
    elseif parties == 4 && sorted_dimensions == (2, 2, 2, 4)
        family, cardinality = :feng_2x2x2x4, 8
    elseif parties == 5 && sorted_dimensions == (2, 2, 2, 2, 5)
        family, cardinality = :feng_2x2x2x2x5, 10
    elseif parties == 8 && all(==(2), sorted_dimensions)
        family, cardinality = :johnston_2_power_8, 11
    elseif parties % 4 == 0 && all(==(2), sorted_dimensions)
        family = :johnston_2_power_4k
        arguments = (parties,)
        cardinality = parties + 4
    elseif isodd(parties) && all(==(2), sorted_dimensions)
        family = :generalized_shifts
        arguments = (parties,)
        cardinality = parties + 1
    else
        even_count = count(iseven, sorted_dimensions)
        if ((isodd(total_excess) || even_count == 0) && even_count <= 1)
            family = :alon_lovasz
            arguments = (sorted_dimensions,)
            cardinality = total_excess + 1
            randomized = true
        elseif sorted_dimensions[end] - 1 ==
               sum(dimension - 1 for dimension in sorted_dimensions[1:(end - 1)]) &&
            total_excess >= 3 &&
            iseven(total_excess)
            family = :chen_johnston_bipartite
            arguments = (sorted_dimensions,)
            cardinality = 2sorted_dimensions[end]
            randomized = true
        elseif parties == 2 && sorted_dimensions == (4, 6)
            family, cardinality = :chen_johnston_4x6, 10
        elseif parties == 3 &&
            sorted_dimensions[1:2] == (2, 2) &&
            sorted_dimensions[3] % 4 == 1
            family = :chen_johnston_4k1
            arguments = (sorted_dimensions[3],)
            cardinality = sorted_dimensions[3] + 3
            randomized = true
        end
    end

    if family === nothing
        minimum = minimum_upb_size(sorted_dimensions)
        if minimum.status === :unknown
            throw(
                UPBConstructionUnavailable(
                    :minimum_unknown,
                    requested_dimensions,
                    nothing,
                    "the reviewed theorem table does not determine a minimum UPB size",
                ),
            )
        end
        throw(
            UPBConstructionUnavailable(
                :known_but_not_in_catalog,
                requested_dimensions,
                minimum.size,
                "a minimum of size $(minimum.size) is known, but the pinned public " *
                "catalog provides no construction for this routed dimension family",
            ),
        )
    end
    randomized =
        randomized || family in (:alon_lovasz, :chen_johnston_bipartite, :chen_johnston_4k1)
    randomized && (
        notes = (
            "The upstream global-RNG search is replaced by an explicit-RNG bounded full-spark search.",
        )
    )
    return _UPBCRequest(
        family,
        arguments,
        sorted_dimensions,
        output_permutation,
        cardinality,
        randomized,
        requested_dimensions,
        notes,
    )
end

function _upbc_complete_product_basis(::Type{T}, dimensions::Tuple) where {T}
    if length(dimensions) == 1
        return (Matrix{T}(I, dimensions[1], dimensions[1]),)
    end
    length(dimensions) == 2 ||
        error("internal complete-product-basis request has more than two parties")
    first_dimension, second_dimension = dimensions
    cardinality = first_dimension * second_dimension
    first = zeros(T, first_dimension, cardinality)
    second = zeros(T, second_dimension, cardinality)
    column = 1
    for first_index in 1:first_dimension, second_index in 1:second_dimension
        first[first_index, column] = one(T)
        second[second_index, column] = one(T)
        column += 1
    end
    return (first, second)
end

function _upbc_build_local(
    request::_UPBCRequest,
    rng::Union{Nothing,AbstractRNG},
    ::Type{T},
    max_attempts::Int,
    rtol::T,
    boundary_factor::T,
    budget::_UPBCBudget,
) where {T<:AbstractFloat}
    family = request.family
    arguments = request.arguments
    if request.randomized
        rng === nothing && throw(
            ArgumentError(
                "the $family UPB construction is randomized; call " *
                "upb(rng, ...) with an explicit AbstractRNG",
            ),
        )
        _upbc_require_randomized_type(T)
    end

    attempts = 0
    factors = if family === :complete_product_basis
        _upbc_complete_product_basis(T, request.sorted_dimensions)
    elseif family === :pyramid
        _upbc_pyramid(T)
    elseif family === :tiles
        _upbc_tiles(T)
    elseif family === :generalized_tiles_1
        _upbc_generalized_tiles_1(T, arguments[1])
    elseif family === :generalized_tiles_2
        _upbc_generalized_tiles_2(T, arguments[1], arguments[2])
    elseif family === :minimum_4x4
        _upbc_minimum_4x4(T)
    elseif family === :quad_residue
        _upbc_quad_residue(T, arguments[1])
    elseif family === :six_parameter
        _upbc_six_parameter(T, arguments[1])
    elseif family === :generalized_shifts
        _upbc_generalized_shifts(T, arguments[1])
    elseif family === :feng_2x2x2x2
        _upbc_feng_2x2x2x2(T)
    elseif family === :johnston_2_power_8
        _upbc_johnston_2_power_8(T)
    elseif family === :johnston_2_power_4k
        _upbc_johnston_2_power_4k(T, arguments[1])
    elseif family === :chen_johnston_4k1
        built, attempts = _upbc_chen_johnston_4k1(
            rng, T, arguments[1], max_attempts, rtol, boundary_factor, budget
        )
        built
    elseif family === :chen_johnston_bipartite
        built, attempts = _upbc_chen_johnston_bipartite(
            rng,
            T,
            request.sorted_dimensions,
            max_attempts,
            rtol,
            boundary_factor,
            budget,
        )
        built
    elseif family === :chen_johnston_4x6
        _upbc_chen_johnston_4x6(T)
    elseif family === :feng_2x2x3
        _upbc_feng_2x2x3(T)
    elseif family === :feng_2x2x5
        _upbc_feng_2x2x5(T)
    elseif family === :feng_4x4
        _upbc_feng_4x4(T)
    elseif family === :feng_2x2x2x4
        _upbc_feng_2x2x2x4(T)
    elseif family === :feng_2x2x2x2x5
        _upbc_feng_2x2x2x2x5(T)
    elseif family === :alon_lovasz
        built, attempts = _upbc_alon_lovasz(
            rng,
            T,
            request.sorted_dimensions,
            max_attempts,
            rtol,
            boundary_factor,
            budget,
        )
        built
    else
        error("internal UPB builder is incomplete for $family")
    end

    reference_key = if family === :complete_product_basis
        if length(request.sorted_dimensions) == 1
            :elementary_product_basis
        else
            :divincenzo_mor_shor_smolin_terhal_2003
        end
    else
        _upbc_closed_form_metadata(family)
    end
    construction_kind = if family === :complete_product_basis
        :complete_product_basis
    elseif request.randomized
        :randomized_full_spark
    else
        :closed_form
    end
    return _UPBCLocalConstruction(
        family, factors, reference_key, construction_kind, attempts, request.notes
    )
end

function _upbc_tolerances(
    ::Type{T}; atol=zero(T), rtol=nothing, boundary_factor=8
) where {T<:AbstractFloat}
    checked_atol = try
        convert(T, atol)
    catch
        throw(ArgumentError("atol must be representable as $T"))
    end
    checked_rtol = if rtol === nothing
        convert(T, 256) * eps(T)
    else
        try
            convert(T, rtol)
        catch
            throw(ArgumentError("rtol must be representable as $T"))
        end
    end
    checked_boundary = try
        convert(T, boundary_factor)
    catch
        throw(ArgumentError("boundary_factor must be representable as $T"))
    end
    isfinite(checked_atol) && checked_atol >= zero(T) ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    isfinite(checked_rtol) && checked_rtol >= zero(T) ||
        throw(ArgumentError("rtol must be finite and nonnegative"))
    isfinite(checked_boundary) && checked_boundary > one(T) ||
        throw(ArgumentError("boundary_factor must be finite and greater than one"))
    return checked_atol, checked_rtol, checked_boundary
end

function _upbc_finalize(
    request::_UPBCRequest,
    built::_UPBCLocalConstruction,
    ::Type{T},
    atol::T,
    rtol::T,
    boundary_factor::T,
    max_attempts::Int,
    max_local_entries::Union{Nothing,Int},
    max_global_entries::Union{Nothing,Int},
    budget::_UPBCBudget,
) where {T<:AbstractFloat}
    sorted_factors = built.factors
    parties = length(request.sorted_dimensions)
    length(sorted_factors) == parties ||
        error("internal UPB construction returned the wrong number of local factors")
    requested_dimensions = ntuple(
        party -> request.sorted_dimensions[request.output_permutation[party]], parties
    )
    requested_dimensions == (
        if request.requested isa Tuple &&
            all(value -> value isa Integer, request.requested)
            request.requested
        else
            requested_dimensions
        end
    ) || error("internal UPB subsystem permutation failed")
    factors_in_requested_order = if built.family === :complete_product_basis
        _upbc_complete_product_basis(T, requested_dimensions)
    else
        ntuple(party -> sorted_factors[request.output_permutation[party]], parties)
    end

    promoted_type = promote_type(map(eltype, factors_in_requested_order)...)
    factors = ntuple(
        party -> Matrix{promoted_type}(factors_in_requested_order[party]), parties
    )
    for party in 1:parties
        size(factors[party]) == (requested_dimensions[party], request.cardinality) || error(
            "internal UPB factor $party has size $(size(factors[party])); " *
            "expected $((requested_dimensions[party], request.cardinality))",
        )
        all(isfinite, factors[party]) ||
            throw(ArgumentError("UPB construction produced nonfinite local entries"))
    end

    local_entries = request.cardinality * sum(requested_dimensions)
    _upbc_guard_entries(
        local_entries,
        max_local_entries,
        :local_entries,
        "local-factor storage exceeds max_local_entries",
    )
    global_dimension = _upbc_checked_product(requested_dimensions, "UPB global dimension")
    global_entries = big(global_dimension) * request.cardinality
    _upbc_guard_entries(
        global_entries,
        max_global_entries,
        :global_entries,
        "global-vector storage exceeds max_global_entries",
    )
    _upbc_consume!(budget, local_entries, "local-factor construction")
    _upbc_consume!(
        budget, Int(global_entries) * parties, "global tensor-product materialization"
    )

    global_vectors = Matrix{promoted_type}(undef, global_dimension, request.cardinality)
    product_residual = zero(T)
    for column in 1:request.cardinality
        tensor_column = _upbc_tensor_column(factors, column)
        global_vectors[:, column] = tensor_column
        product_residual = max(
            product_residual,
            convert(T, norm(@view(global_vectors[:, column]) - tensor_column)),
        )
    end

    normalization_residual = zero(T)
    for factor in factors, column in 1:request.cardinality
        normalization_residual = max(
            normalization_residual, abs(convert(T, norm(@view factor[:, column])) - one(T))
        )
    end
    for column in 1:request.cardinality
        normalization_residual = max(
            normalization_residual,
            abs(convert(T, norm(@view global_vectors[:, column])) - one(T)),
        )
    end

    pair_count = request.cardinality * (request.cardinality - 1) ÷ 2
    _upbc_consume!(
        budget,
        pair_count * sum(requested_dimensions),
        "pairwise product-state orthogonality verification",
    )
    orthogonality_residual = zero(T)
    for right in 2:request.cardinality, left in 1:(right - 1)
        overlap = one(promoted_type)
        for factor in factors
            overlap *= dot(@view(factor[:, left]), @view(factor[:, right]))
        end
        orthogonality_residual = max(orthogonality_residual, convert(T, abs(overlap)))
    end

    scale = convert(
        T, max(1, request.cardinality, global_dimension, sum(requested_dimensions))
    )
    threshold = atol + rtol * scale
    robust_threshold = boundary_factor * threshold
    normalization_residual <= robust_threshold || throw(
        ArgumentError(
            "constructed local factors are not robustly normalized: residual " *
            "$normalization_residual exceeds $robust_threshold",
        ),
    )
    orthogonality_residual <= robust_threshold || throw(
        ArgumentError(
            "constructed product states are not robustly pairwise orthogonal: " *
            "residual $orthogonality_residual exceeds $robust_threshold",
        ),
    )

    reference, reference_url = _upbc_reference(built.reference_key)
    verification_kind = if built.construction_kind === :complete_product_basis
        :exact_structure
    else
        :tolerance_robust
    end
    arithmetic = if T === BigFloat
        :arbitrary_precision_floating
    elseif T === Float32
        :float32
    else
        :float64
    end
    notes = (
        built.notes...,
        "The nested symbolic NICE helpers in the pinned file are unreachable " *
        "from the executable public UPB route; no symbolic simplification is needed.",
    )
    return UPBConstruction(
        :constructed,
        built.family,
        request.requested,
        requested_dimensions,
        request.cardinality,
        factors,
        global_vectors,
        built.reference_key,
        reference,
        reference_url,
        built.construction_kind,
        verification_kind,
        !request.randomized,
        request.randomized,
        true,
        true,
        normalization_residual,
        orthogonality_residual,
        product_residual,
        built.attempts,
        budget.random_draws,
        budget.minors_checked,
        budget.work_used,
        max_attempts,
        max_local_entries,
        max_global_entries,
        budget.max_work,
        budget.max_minors,
        arithmetic,
        :not_required_for_executable_public_route,
        _UPBC_SOURCE_REVISION,
        notes,
    )
end

function _upbc_construct(
    rng::Union{Nothing,AbstractRNG},
    request::_UPBCRequest;
    real_type=Float64,
    atol=0,
    rtol=nothing,
    boundary_factor=8,
    max_attempts=_UPBC_DEFAULT_MAX_ATTEMPTS,
    max_local_entries=_UPBC_DEFAULT_MAX_LOCAL_ENTRIES,
    max_global_entries=_UPBC_DEFAULT_MAX_GLOBAL_ENTRIES,
    max_work=_UPBC_DEFAULT_MAX_WORK,
    max_minors=_UPBC_DEFAULT_MAX_MINORS,
)
    T = _upbc_real_type(real_type)
    checked_atol, checked_rtol, checked_boundary = _upbc_tolerances(
        T; atol=atol, rtol=rtol, boundary_factor=boundary_factor
    )
    checked_attempts = _upbc_limit(max_attempts, "max_attempts"; positive=true)
    checked_attempts === nothing && throw(ArgumentError("max_attempts cannot be disabled"))
    checked_local_entries = _upbc_limit(max_local_entries, "max_local_entries")
    checked_global_entries = _upbc_limit(max_global_entries, "max_global_entries")
    checked_work = _upbc_limit(max_work, "max_work")
    checked_minors = _upbc_limit(max_minors, "max_minors")

    local_entries = big(request.cardinality) * sum(request.sorted_dimensions)
    _upbc_guard_entries(
        local_entries,
        checked_local_entries,
        :local_entries,
        "predicted local-factor storage exceeds max_local_entries",
    )
    global_dimension_big = prod(big.(request.sorted_dimensions))
    global_entries = global_dimension_big * request.cardinality
    _upbc_guard_entries(
        global_entries,
        checked_global_entries,
        :global_entries,
        "predicted global-vector storage exceeds max_global_entries",
    )

    budget = _UPBCBudget(checked_work, checked_minors, 0, 0, 0)
    built = _upbc_build_local(
        request, rng, T, checked_attempts, checked_rtol, checked_boundary, budget
    )
    return _upbc_finalize(
        request,
        built,
        T,
        checked_atol,
        checked_rtol,
        checked_boundary,
        checked_attempts,
        checked_local_entries,
        checked_global_entries,
        budget,
    )
end

"""
    upb(name, family_arguments...; kwargs...) -> UPBConstruction
    upb(dimensions; kwargs...) -> UPBConstruction
    upb(rng::AbstractRNG, name_or_dimensions, family_arguments...; kwargs...)

Construct a family from the pinned public QETLAB UPB catalog.

Keywords are `real_type=Float64`, `atol=0`, `rtol=nothing`,
`boundary_factor=8`, `max_attempts=64`, `max_local_entries=1_000_000`,
`max_global_entries=2_000_000`, `max_work=50_000_000`, and
`max_minors=250_000`. Inputs must be finite family parameters or positive
integer subsystem dimensions in the domains listed below; no subsystem
dimension, vector, or parameter is silently repaired.

Named closed-form families are:

- `:pyramid`, `:tiles`, `:minimum_4x4`, `:quad_residue`,
  `:six_parameter`, `:shifts`, and `:generalized_shifts`;
- `:generalized_tiles_1` and `:generalized_tiles_2`;
- the Feng families `:feng_2x2x3`, `:feng_2x2x5`,
  `:feng_2x2x2x2`, `:feng_4x4`, `:feng_2x2x2x4`, and
  `:feng_2x2x2x2x5`;
- `:johnston_2_power_8`, `:johnston_2_power_4k`, and
  `:chen_johnston_4x6`.

The dimension-driven Alon–Lovász, Chen–Johnston bipartite-style, and
Chen–Johnston `(2,2,4k+1)` constructions are randomized. They require the
explicit-RNG form and use bounded, recorded full-spark searches. Their expert
named forms `:alon_lovasz`, `:chen_johnston_bipartite`, and
`:chen_johnston_4k1` are also available with the corresponding dimension
argument.

`dimensions` follows the same ordered dispatch as the pinned source. The
returned factors are restored to the caller's subsystem order. If the reviewed
minimum is unknown, or is known but the pinned catalog deliberately has no
construction, [`UPBConstructionUnavailable`](@ref) is thrown with the exact
distinction and known size when available.

All vectors are materialized because the public result owns both local and
global representations. `max_local_entries` and `max_global_entries` guard
those allocations before construction; `max_work`, `max_minors`, and
`max_attempts` bound validation and randomized searches. A limit may be
disabled with `nothing` except `max_attempts`. Randomized branches currently
support `Float32` and `Float64`; deterministic branches also support
`BigFloat`.

If there are `s` product states in local dimensions `d₁,…,dₚ`, the owned
output requires `s*sum(dᵢ)` local entries and `s*prod(dᵢ)` global entries.
Final tensor reconstruction and Gram validation use at least those terms and
`O(s^2*prod(dᵢ))` scalar work. Randomized branches additionally enumerate
the recorded family-specific full-spark minors, bounded by `max_minors` and
`max_work`.

The source's nested Symbolic Math Toolbox `NICE` code is not reachable from
the executable top-level `UPB` routes. This implementation therefore performs
the complete numeric construction without a symbolic dependency and records
that disposition in every result. It also rejects `generalized_tiles_1(2)`,
which the pinned source accepts even though the resulting singleton is
extendible.

# Examples

```jldoctest
julia> using LinearAlgebra

julia> construction = upb(:tiles);

julia> (construction.dimensions, construction.cardinality)
((3, 3), 5)

julia> maximum(abs, construction.global_vectors' * construction.global_vectors - I) < 1e-14
true
```
"""
function upb(name::Union{Symbol,AbstractString}, arguments...; kwargs...)
    return _upbc_construct(nothing, _upbc_named_request(name, arguments...); kwargs...)
end

function upb(dimensions::Union{Tuple,AbstractVector}; kwargs...)
    return _upbc_construct(nothing, _upbc_dimension_request(dimensions); kwargs...)
end

function upb(rng::AbstractRNG, name::Union{Symbol,AbstractString}, arguments...; kwargs...)
    return _upbc_construct(rng, _upbc_named_request(name, arguments...); kwargs...)
end

function upb(rng::AbstractRNG, dimensions::Union{Tuple,AbstractVector}; kwargs...)
    return _upbc_construct(rng, _upbc_dimension_request(dimensions); kwargs...)
end
