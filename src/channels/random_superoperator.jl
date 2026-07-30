# Source-informed independent Julia implementation based on the specification
# and QETLAB RandomSuperoperator.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

using Random

const _RANDOM_SUPEROPERATOR_DEFAULT_MAX_ATTEMPTS = 8
const _RANDOM_SUPEROPERATOR_DEFAULT_MAX_ITERATIONS = 1_000
const _RANDOM_SUPEROPERATOR_DEFAULT_MAX_DIMENSION = 4_096
const _RANDOM_SUPEROPERATOR_DEFAULT_MAX_ENTRIES = 10_000_000
const _RANDOM_SUPEROPERATOR_DEFAULT_MAX_WORK = 1_000_000_000

"""
    RandomSuperoperatorResult

Structured output of [`random_superoperator`](@ref). `representation` has the
requested representation kind and `kraus_representation` retains the
construction certificate for complete positivity. The latter is an algebraic
guarantee; `trace_preservation_residual`, `unitality_residual`, and
`proportional_unitality_residual` are separate floating-point diagnostics.

`status` is `:success`, `:max_attempts`, or `:work_limit`. Only
`:success` sets `succeeded=true` and establishes the requested marginal
guarantees. A failed result still contains the final completely positive
random draw, but that map must not be treated as trace preserving or unital.

For an explicitly requested unequal-dimensional proportional-unital
construction, `proportional_unital_factor == input_dimension /
output_dimension` and
`Φ(I_input) ≈ proportional_unital_factor * I_output`. This is not unitality.
"""
struct RandomSuperoperatorResult{
    M<:AbstractMapRepresentation,
    K<:KrausRepresentation,
    R<:AbstractFloat,
    V<:AbstractVector{R},
    H<:AbstractVector{R},
}
    representation::M
    kraus_representation::K
    status::Symbol
    last_attempt_status::Symbol
    succeeded::Bool
    attempts::Int
    input_dimension::Int
    output_dimension::Int
    requested_kraus_rank::Int
    numerical_kraus_rank::Int
    real_output::Bool
    requested_trace_preserving::Bool
    requested_unital::Bool
    requested_proportional_unital::Bool
    complete_positivity_guaranteed::Bool
    trace_preservation_guaranteed::Bool
    unitality_guaranteed::Bool
    proportional_unitality_guaranteed::Bool
    choi_trace::R
    trace_preservation_residual::R
    unitality_residual::R
    proportional_unital_factor::R
    proportional_unitality_residual::R
    trace_preservation_tolerance::R
    unitality_tolerance::R
    factor_singular_values::V
    rank_threshold::R
    balancing_iterations::Int
    balancing_residual_history::H
    minimum_marginal_eigenvalue::R
    maximum_filter_condition::R
    failed_marginal::Union{Nothing,Symbol}
    work_used::BigInt
    message::String
end

function Base.show(io::IO, result::RandomSuperoperatorResult)
    return print(
        io,
        "RandomSuperoperatorResult(status=",
        result.status,
        ", dimensions=",
        (result.input_dimension, result.output_dimension),
        ", requested_rank=",
        result.requested_kraus_rank,
        ", numerical_rank=",
        result.numerical_kraus_rank,
        ", attempts=",
        result.attempts,
        ")",
    )
end

function _random_superoperator_flag(value, name::AbstractString)
    value isa Bool || throw(ArgumentError("$name must be Bool; got $(repr(value))"))
    return value
end

function _random_superoperator_dimensions(dim)
    if dim isa Integer
        dimension = _positive_int(dim, "dim")
        return dimension, dimension
    elseif dim isa Tuple || dim isa AbstractVector
        dim isa AbstractVector && Base.require_one_based_indexing(dim)
        length(dim) == 2 || throw(
            DimensionMismatch(
                "dim must be a scalar or contain input and output dimensions; " *
                "got $(length(dim)) entries",
            ),
        )
        return _positive_int(dim[1], "dim[1]"), _positive_int(dim[2], "dim[2]")
    end
    return throw(
        ArgumentError("dim must be a positive integer or a two-entry tuple/vector")
    )
end

function _random_superoperator_limit(value, name::AbstractString; allow_nothing::Bool=true)
    allow_nothing && value === nothing && return nothing
    value isa Bool && throw(
        ArgumentError(
            "$name must be a positive integer" *
            (allow_nothing ? " or nothing" : "") *
            ", not Bool",
        ),
    )
    value isa Integer || throw(
        ArgumentError(
            "$name must be a positive integer" * (allow_nothing ? " or nothing" : "")
        ),
    )
    value > 0 || throw(ArgumentError("$name must be positive; got $value"))
    return BigInt(value)
end

function _random_superoperator_iterations(value)
    value isa Bool &&
        throw(ArgumentError("max_iterations must be a nonnegative integer, not Bool"))
    value isa Integer ||
        throw(ArgumentError("max_iterations must be a nonnegative integer"))
    value >= 0 || throw(ArgumentError("max_iterations must be nonnegative; got $value"))
    return try
        Int(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("max_iterations=$value cannot be represented as Int"))
    end
end

function _random_superoperator_real_type(value)
    value isa Type ||
        throw(ArgumentError("T must be Float32 or Float64; got $(repr(value))"))
    value === Float32 ||
        value === Float64 ||
        throw(
            ArgumentError(
                "the bounded spectral path supports T=Float32 or T=Float64; got $value"
            ),
        )
    return value
end

function _random_superoperator_representation(value)
    value isa Symbol || throw(
        ArgumentError(
            "representation must be :kraus, :choi, or :superoperator; got $(repr(value))",
        ),
    )
    value in (:kraus, :choi, :superoperator) || throw(
        ArgumentError(
            "representation must be :kraus, :choi, or :superoperator; got $value"
        ),
    )
    return value
end

function _random_superoperator_tolerance(value, name::AbstractString, ::Type{R}) where {R}
    value isa Bool &&
        throw(ArgumentError("$name must be a finite nonnegative real number, not Bool"))
    value isa Real || throw(ArgumentError("$name must be a finite nonnegative real number"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= zero(value) || throw(ArgumentError("$name must be nonnegative; got $value"))
    converted = try
        R(value)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("$name=$value cannot be represented as $R"))
    end
    isfinite(converted) ||
        throw(ArgumentError("$name=$value is not finite when represented as $R"))
    return converted
end

function _random_superoperator_positive_real(
    value, name::AbstractString, ::Type{R}
) where {R}
    checked = _random_superoperator_tolerance(value, name, R)
    checked > zero(R) || throw(ArgumentError("$name must be strictly positive; got $value"))
    return checked
end

function _random_superoperator_reserve!(used::Base.RefValue{BigInt}, cost, limit)
    checked_cost = BigInt(cost)
    if limit !== nothing && used[] + checked_cost > limit
        return false
    end
    used[] += checked_cost
    return true
end

function _random_superoperator_peak_entries(
    product_dimension::Int,
    rank_bound::Int,
    input_dimension::Int,
    output_dimension::Int,
    representation::Symbol,
)
    factor_entries = BigInt(product_dimension) * rank_bound
    matrix_entries = BigInt(product_dimension)^2
    representation_entries = representation === :kraus ? 0 : 2 * matrix_entries
    return (
        4 * factor_entries +
        representation_entries +
        BigInt(input_dimension)^2 +
        BigInt(output_dimension)^2
    )
end

function _random_superoperator_attempt_work(
    input_dimension::Int, output_dimension::Int, rank_bound::Int
)
    product_dimension = BigInt(input_dimension) * output_dimension
    factor_entries = product_dimension * rank_bound
    input_products = BigInt(rank_bound) * output_dimension * BigInt(input_dimension)^2
    output_products = BigInt(rank_bound) * input_dimension * BigInt(output_dimension)^2
    return (
        3 * factor_entries +
        2 * input_products +
        2 * output_products +
        BigInt(input_dimension)^3 +
        BigInt(output_dimension)^3
    )
end

function _random_superoperator_final_work(
    product_dimension::Int, rank_bound::Int, representation::Symbol
)
    rank_diagnostic =
        BigInt(product_dimension) * BigInt(rank_bound)^2 + BigInt(rank_bound)^3
    representation_work =
        representation === :kraus ? 0 : BigInt(rank_bound) * BigInt(product_dimension)^2
    return rank_diagnostic + representation_work
end

function _random_superoperator_draw(
    rng::AbstractRNG,
    input_dimension::Int,
    output_dimension::Int,
    rank_bound::Int,
    real_output::Bool,
    ::Type{R},
) where {R<:AbstractFloat}
    product_dimension = _checked_product(
        (input_dimension, output_dimension), "input_dimension*output_dimension"
    )
    factor_matrix = if real_output
        randn(rng, R, product_dimension, rank_bound)
    else
        complex.(
            randn(rng, R, product_dimension, rank_bound),
            randn(rng, R, product_dimension, rank_bound),
        )
    end
    all(isfinite, factor_matrix) || return (status=:invalid_random_draw, factors=nothing)
    normalization = norm(factor_matrix)
    isfinite(normalization) && normalization > zero(R) ||
        return (status=:invalid_random_draw, factors=nothing)
    factor_matrix ./= normalization
    factors = [
        copy(reshape(view(factor_matrix, :, index), output_dimension, input_dimension)) for
        index in 1:rank_bound
    ]
    return (status=:success, factors=factors)
end

function _random_superoperator_marginals(
    factors, input_dimension::Int, output_dimension::Int
)
    input_marginal = zeros(eltype(first(factors)), input_dimension, input_dimension)
    output_marginal = zeros(eltype(first(factors)), output_dimension, output_dimension)
    for factor in factors
        input_marginal .+= adjoint(factor) * factor
        output_marginal .+= factor * adjoint(factor)
    end
    return input_marginal, output_marginal
end

function _random_superoperator_inverse_sqrt(
    marginal::AbstractMatrix, condition_limit, ::Type{R}
) where {R<:AbstractFloat}
    dimension = size(marginal, 1)
    hermiticity_residual = maximum(abs, marginal - adjoint(marginal); init=zero(R))
    scale = maximum(abs, marginal; init=zero(R))
    roundoff_tolerance = 32 * dimension * eps(R) * max(one(R), scale)
    hermiticity_residual <= roundoff_tolerance || return (
        status=:numerical_failure,
        filter=nothing,
        minimum_eigenvalue=R(NaN),
        condition=R(Inf),
    )

    # Averaging is confined to the eigensolver work matrix after a derived
    # roundoff check. The random factors and returned map are never repaired.
    work_matrix = (marginal + adjoint(marginal)) / 2
    decomposition = try
        eigen(Hermitian(work_matrix))
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        return (
            status=:numerical_failure,
            filter=nothing,
            minimum_eigenvalue=R(NaN),
            condition=R(Inf),
        )
    end
    minimum_eigenvalue = minimum(decomposition.values)
    maximum_eigenvalue = maximum(decomposition.values)
    minimum_eigenvalue > zero(R) || return (
        status=:singular_marginal,
        filter=nothing,
        minimum_eigenvalue=minimum_eigenvalue,
        condition=R(Inf),
    )
    filter_condition = sqrt(maximum_eigenvalue / minimum_eigenvalue)
    isfinite(filter_condition) && filter_condition < condition_limit || return (
        status=:ill_conditioned,
        filter=nothing,
        minimum_eigenvalue=minimum_eigenvalue,
        condition=filter_condition,
    )
    inverse_square_roots = map(value -> inv(sqrt(value)), decomposition.values)
    filter =
        decomposition.vectors *
        Diagonal(inverse_square_roots) *
        adjoint(decomposition.vectors)
    all(isfinite, filter) || return (
        status=:numerical_failure,
        filter=nothing,
        minimum_eigenvalue=minimum_eigenvalue,
        condition=filter_condition,
    )
    return (
        status=:success,
        filter=filter,
        minimum_eigenvalue=minimum_eigenvalue,
        condition=filter_condition,
    )
end

function _random_superoperator_residuals(
    factors, input_dimension::Int, output_dimension::Int, proportional_factor, ::Type{R}
) where {R<:AbstractFloat}
    input_marginal, output_marginal = _random_superoperator_marginals(
        factors, input_dimension, output_dimension
    )
    input_identity = Matrix{eltype(input_marginal)}(I, input_dimension, input_dimension)
    output_identity = Matrix{eltype(output_marginal)}(I, output_dimension, output_dimension)
    return (
        choi_trace=R(real(tr(input_marginal))),
        trace_preservation=R(norm(input_marginal - input_identity)),
        unitality=R(norm(output_marginal - output_identity)),
        proportional_unitality=R(
            norm(output_marginal - proportional_factor * output_identity)
        ),
    )
end

function _random_superoperator_single_normalization(
    factors,
    input_dimension::Int,
    output_dimension::Int,
    side::Symbol,
    condition_limit,
    ::Type{R},
) where {R<:AbstractFloat}
    input_marginal, output_marginal = _random_superoperator_marginals(
        factors, input_dimension, output_dimension
    )
    marginal = side === :input ? input_marginal : output_marginal
    data = _random_superoperator_inverse_sqrt(marginal, condition_limit, R)
    data.status === :success || return (
        factors=factors,
        status=data.status,
        iterations=0,
        residual_history=R[],
        minimum_eigenvalue=data.minimum_eigenvalue,
        maximum_condition=data.condition,
        failed_marginal=side,
    )
    normalized = if side === :input
        [factor * data.filter for factor in factors]
    else
        [data.filter * factor for factor in factors]
    end
    return (
        factors=normalized,
        status=:converged,
        iterations=1,
        residual_history=R[],
        minimum_eigenvalue=data.minimum_eigenvalue,
        maximum_condition=data.condition,
        failed_marginal=nothing,
    )
end

function _random_superoperator_balancing(
    factors,
    input_dimension::Int,
    output_dimension::Int,
    proportional_factor::R,
    trace_tolerance::R,
    output_tolerance::R,
    condition_limit::R,
    max_iterations::Int,
    used::Base.RefValue{BigInt},
    max_work,
) where {R<:AbstractFloat}
    scalar_type = eltype(first(factors))
    left_filter = Matrix{scalar_type}(I, output_dimension, output_dimension)
    right_filter = Matrix{scalar_type}(I, input_dimension, input_dimension)
    residual_history = R[]
    minimum_eigenvalue = R(Inf)
    maximum_condition = one(R)
    sweep_work = _random_superoperator_attempt_work(
        input_dimension, output_dimension, length(factors)
    )

    for iteration in 0:max_iterations
        residuals = _random_superoperator_residuals(
            factors, input_dimension, output_dimension, proportional_factor, R
        )
        push!(
            residual_history,
            residuals.trace_preservation + residuals.proportional_unitality,
        )
        if residuals.trace_preservation <= trace_tolerance &&
            residuals.proportional_unitality <= output_tolerance
            return (
                factors=factors,
                status=:converged,
                iterations=iteration,
                residual_history=residual_history,
                minimum_eigenvalue=minimum_eigenvalue,
                maximum_condition=maximum_condition,
                failed_marginal=nothing,
            )
        end
        iteration == max_iterations && break
        _random_superoperator_reserve!(used, sweep_work, max_work) || return (
            factors=factors,
            status=:work_limit,
            iterations=iteration,
            residual_history=residual_history,
            minimum_eigenvalue=minimum_eigenvalue,
            maximum_condition=maximum_condition,
            failed_marginal=nothing,
        )

        input_marginal, _ = _random_superoperator_marginals(
            factors, input_dimension, output_dimension
        )
        input_data = _random_superoperator_inverse_sqrt(input_marginal, condition_limit, R)
        minimum_eigenvalue = min(minimum_eigenvalue, input_data.minimum_eigenvalue)
        input_data.status === :success || return (
            factors=factors,
            status=input_data.status,
            iterations=iteration,
            residual_history=residual_history,
            minimum_eigenvalue=minimum_eigenvalue,
            maximum_condition=max(maximum_condition, input_data.condition),
            failed_marginal=:input,
        )
        candidate_right_filter = right_filter * input_data.filter
        candidate_right_condition = try
            R(cond(candidate_right_filter))
        catch err
            err isa LinearAlgebra.LAPACKException || rethrow()
            R(Inf)
        end
        if !isfinite(candidate_right_condition) ||
            candidate_right_condition >= condition_limit
            return (
                factors=factors,
                status=:ill_conditioned,
                iterations=iteration,
                residual_history=residual_history,
                minimum_eigenvalue=minimum_eigenvalue,
                maximum_condition=max(maximum_condition, candidate_right_condition),
                failed_marginal=:input,
            )
        end
        factors = [factor * input_data.filter for factor in factors]
        right_filter = candidate_right_filter
        maximum_condition = max(maximum_condition, candidate_right_condition)

        _, output_marginal = _random_superoperator_marginals(
            factors, input_dimension, output_dimension
        )
        output_data = _random_superoperator_inverse_sqrt(
            output_marginal, condition_limit, R
        )
        minimum_eigenvalue = min(minimum_eigenvalue, output_data.minimum_eigenvalue)
        output_data.status === :success || return (
            factors=factors,
            status=output_data.status,
            iterations=iteration,
            residual_history=residual_history,
            minimum_eigenvalue=minimum_eigenvalue,
            maximum_condition=max(maximum_condition, output_data.condition),
            failed_marginal=:output,
        )
        scaled_output_filter = sqrt(proportional_factor) * output_data.filter
        candidate_left_filter = scaled_output_filter * left_filter
        candidate_left_condition = try
            R(cond(candidate_left_filter))
        catch err
            err isa LinearAlgebra.LAPACKException || rethrow()
            R(Inf)
        end
        if !isfinite(candidate_left_condition) ||
            candidate_left_condition >= condition_limit
            return (
                factors=factors,
                status=:ill_conditioned,
                iterations=iteration,
                residual_history=residual_history,
                minimum_eigenvalue=minimum_eigenvalue,
                maximum_condition=max(maximum_condition, candidate_left_condition),
                failed_marginal=:output,
            )
        end
        factors = [scaled_output_filter * factor for factor in factors]
        left_filter = candidate_left_filter
        maximum_condition = max(maximum_condition, candidate_left_condition)
    end

    return (
        factors=factors,
        status=:max_iterations,
        iterations=max_iterations,
        residual_history=residual_history,
        minimum_eigenvalue=minimum_eigenvalue,
        maximum_condition=maximum_condition,
        failed_marginal=nothing,
    )
end

function _random_superoperator_rank_diagnostic(factors, atol, rtol, ::Type{R}) where {R}
    factor_matrix = hcat((vec(factor) for factor in factors)...)
    singular_values = R.(svdvals(factor_matrix))
    scale = max(first(singular_values), one(R))
    threshold = atol + rtol * scale
    numerical_rank = count(>(threshold), singular_values)
    return singular_values, threshold, numerical_rank
end

function _random_superoperator_convert(kraus::KrausRepresentation, kind::Symbol)
    kind === :kraus && return kraus
    kind === :choi && return choi_representation(kraus)
    return superoperator_representation(kraus)
end

function _random_superoperator_message(status::Symbol, last_attempt_status::Symbol)
    status === :success &&
        return "the requested completely positive random-map construction passed its residual checks"
    status === :work_limit &&
        return "the conservative work budget was exhausted before a valid requested normalization was obtained"
    return (
        "all bounded random draws were exhausted; the returned map is completely " *
        "positive but the requested marginal guarantees were not established " *
        "(last_attempt_status=$last_attempt_status)"
    )
end

"""
    random_superoperator(
        rng::AbstractRNG,
        dim;
        trace_preserving=true,
        unital=false,
        proportional_unital=false,
        real=false,
        kraus_rank=nothing,
        representation=:kraus,
        T=Float64,
        atol=0,
        rtol=sqrt(eps(T)),
        max_attempts=8,
        max_iterations=1000,
        max_condition_number=1/sqrt(eps(T)),
        max_dimension=4096,
        max_entries=10_000_000,
        max_work=1_000_000_000,
    ) -> RandomSuperoperatorResult

Draw a dense completely positive map from a rank-controlled Ginibre
factorization. `dim` is either one dimension or `(input_dim, output_dim)`.
The explicit `rng` is the only source of randomness.

The construction has at most `kraus_rank` Kraus factors and has exactly that
Choi rank with probability one; the returned numerical rank and singular
values remain diagnostics rather than a deterministic promise about a random
draw. `kraus_rank` defaults to `input_dim * output_dim`.

For a trace-preserving draw, the input Choi marginal is scaled to identity.
For a unital draw, the output Choi marginal is scaled to identity. Requesting
both is valid only when the input and output dimensions agree. The pinned
QETLAB routine nevertheless accepts unequal dimensions and actually produces

```math
\\Phi(I_{\\mathrm{in}})
\\approx \\frac{d_{\\mathrm{in}}}{d_{\\mathrm{out}}}I_{\\mathrm{out}},
```

which is not unitality. That executable branch is available only through the
explicit corrected keyword combination
`trace_preserving=true, unital=false, proportional_unital=true`.

Simultaneous marginal scaling uses a bounded factor-space Sinkhorn iteration.
It never projects a failed draw, and it preserves complete positivity and the
requested rank bound algebraically. A singular, ill-conditioned, nonconverged,
or work-limited attempt is retained in diagnostics and retried only up to
`max_attempts`. Exhaustion returns a failed structured result.

`representation` is `:kraus`, `:choi`, or `:superoperator`. Random draws are
dense almost surely; no sparse-output option is provided. `T` is `Float32` or
`Float64`, with `real=false` returning the corresponding complex type.
`max_dimension`, `max_entries`, and `max_work` preflight local dimensions,
conservative peak storage, and conservative arithmetic work. Set a resource
limit to `nothing` only after reviewing the cost.
"""
function random_superoperator(
    rng::AbstractRNG,
    dim;
    trace_preserving=true,
    unital=false,
    proportional_unital=false,
    real=false,
    kraus_rank=nothing,
    representation=:kraus,
    T=Float64,
    atol=0,
    rtol=nothing,
    max_attempts=_RANDOM_SUPEROPERATOR_DEFAULT_MAX_ATTEMPTS,
    max_iterations=_RANDOM_SUPEROPERATOR_DEFAULT_MAX_ITERATIONS,
    max_condition_number=nothing,
    max_dimension=_RANDOM_SUPEROPERATOR_DEFAULT_MAX_DIMENSION,
    max_entries=_RANDOM_SUPEROPERATOR_DEFAULT_MAX_ENTRIES,
    max_work=_RANDOM_SUPEROPERATOR_DEFAULT_MAX_WORK,
)
    input_dimension, output_dimension = _random_superoperator_dimensions(dim)
    requested_trace_preserving = _random_superoperator_flag(
        trace_preserving, "trace_preserving"
    )
    requested_unital = _random_superoperator_flag(unital, "unital")
    requested_proportional_unital = _random_superoperator_flag(
        proportional_unital, "proportional_unital"
    )
    real_output = _random_superoperator_flag(real, "real")
    requested_unital &&
        requested_proportional_unital &&
        throw(ArgumentError("unital and proportional_unital cannot both be true"))
    requested_proportional_unital &&
        !requested_trace_preserving &&
        throw(ArgumentError("proportional_unital=true requires trace_preserving=true"))
    requested_proportional_unital &&
        input_dimension == output_dimension &&
        throw(
            ArgumentError(
                "equal dimensions admit strict unitality; request unital=true instead of proportional_unital=true",
            ),
        )
    requested_trace_preserving &&
        requested_unital &&
        input_dimension != output_dimension &&
        throw(
            ArgumentError(
                "a map cannot be both trace preserving and unital when input and " *
                "output dimensions differ; use unital=false, " *
                "proportional_unital=true to request the corrected proportional branch",
            ),
        )

    real_type = _random_superoperator_real_type(T)
    selected_representation = _random_superoperator_representation(representation)
    product_dimension = _checked_product(
        (input_dimension, output_dimension), "input_dimension*output_dimension"
    )
    rank_bound =
        kraus_rank === nothing ? product_dimension : _positive_int(kraus_rank, "kraus_rank")
    rank_bound <= product_dimension || throw(
        ArgumentError(
            "kraus_rank=$rank_bound exceeds the maximum Choi rank $product_dimension"
        ),
    )
    requested_trace_preserving &&
        BigInt(rank_bound) * output_dimension < input_dimension &&
        throw(
            ArgumentError(
                "kraus_rank=$rank_bound is too small for a trace-preserving map " *
                "from dimension $input_dimension to $output_dimension",
            ),
        )
    (requested_unital || requested_proportional_unital) &&
        BigInt(rank_bound) * input_dimension < output_dimension &&
        throw(
            ArgumentError(
                "kraus_rank=$rank_bound is too small for the requested output " *
                "identity marginal in dimension $output_dimension",
            ),
        )

    attempt_limit = _random_superoperator_limit(
        max_attempts, "max_attempts"; allow_nothing=false
    )
    attempt_limit <= typemax(Int) ||
        throw(ArgumentError("max_attempts=$max_attempts cannot be represented as Int"))
    checked_attempts = Int(attempt_limit)
    checked_iterations = _random_superoperator_iterations(max_iterations)
    dimension_limit = _random_superoperator_limit(max_dimension, "max_dimension")
    dimension_limit !== nothing &&
        max(input_dimension, output_dimension) > dimension_limit &&
        throw(
            ArgumentError(
                "local dimension $(max(input_dimension, output_dimension)) exceeds " *
                "max_dimension=$dimension_limit",
            ),
        )
    entry_limit = _random_superoperator_limit(max_entries, "max_entries")
    required_entries = _random_superoperator_peak_entries(
        product_dimension,
        rank_bound,
        input_dimension,
        output_dimension,
        selected_representation,
    )
    entry_limit !== nothing &&
        required_entries > entry_limit &&
        throw(
            ArgumentError(
                "the conservative peak-storage estimate $required_entries exceeds " *
                "max_entries=$entry_limit",
            ),
        )
    work_limit = _random_superoperator_limit(max_work, "max_work")
    used_work = Ref(BigInt(0))
    final_work = _random_superoperator_final_work(
        product_dimension, rank_bound, selected_representation
    )
    attempt_work = _random_superoperator_attempt_work(
        input_dimension, output_dimension, rank_bound
    )
    work_limit !== nothing &&
        final_work + attempt_work > work_limit &&
        throw(
            ArgumentError(
                "one random draw plus final diagnostics require conservative work " *
                "$(final_work + attempt_work), exceeding max_work=$work_limit",
            ),
        )
    _random_superoperator_reserve!(used_work, final_work, work_limit) ||
        error("internal random-superoperator work preflight inconsistency")

    checked_atol = _random_superoperator_tolerance(atol, "atol", real_type)
    checked_rtol = if rtol === nothing
        sqrt(eps(real_type))
    else
        _random_superoperator_tolerance(rtol, "rtol", real_type)
    end
    condition_limit = if max_condition_number === nothing
        inv(sqrt(eps(real_type)))
    else
        _random_superoperator_positive_real(
            max_condition_number, "max_condition_number", real_type
        )
    end
    proportional_factor = real_type(input_dimension) / real_type(output_dimension)
    trace_tolerance =
        checked_atol + checked_rtol * max(one(real_type), sqrt(real_type(input_dimension)))
    output_target_factor =
        requested_proportional_unital ? proportional_factor : one(real_type)
    unital_tolerance =
        checked_atol +
        checked_rtol *
        max(one(real_type), abs(output_target_factor) * sqrt(real_type(output_dimension)))

    last_factors = nothing
    last_attempt_status = :not_started
    last_balancing_iterations = 0
    last_residual_history = real_type[]
    last_minimum_eigenvalue = real_type(NaN)
    last_maximum_condition = one(real_type)
    last_failed_marginal = nothing
    last_attempt = 0
    success = false
    terminal_status = :max_attempts

    for attempt in 1:checked_attempts
        if !_random_superoperator_reserve!(used_work, attempt_work, work_limit)
            terminal_status = :work_limit
            last_attempt_status = :work_limit
            break
        end
        last_attempt = attempt
        draw = _random_superoperator_draw(
            rng, input_dimension, output_dimension, rank_bound, real_output, real_type
        )
        if draw.status !== :success
            last_attempt_status = draw.status
            continue
        end
        factors = draw.factors
        last_factors = factors

        normalization =
            if requested_trace_preserving &&
                (requested_unital || requested_proportional_unital)
                _random_superoperator_balancing(
                    factors,
                    input_dimension,
                    output_dimension,
                    output_target_factor,
                    trace_tolerance,
                    unital_tolerance,
                    condition_limit,
                    checked_iterations,
                    used_work,
                    work_limit,
                )
            elseif requested_trace_preserving
                _random_superoperator_single_normalization(
                    factors,
                    input_dimension,
                    output_dimension,
                    :input,
                    condition_limit,
                    real_type,
                )
            elseif requested_unital
                _random_superoperator_single_normalization(
                    factors,
                    input_dimension,
                    output_dimension,
                    :output,
                    condition_limit,
                    real_type,
                )
            else
                (
                    factors=factors,
                    status=:converged,
                    iterations=0,
                    residual_history=real_type[],
                    minimum_eigenvalue=real_type(NaN),
                    maximum_condition=one(real_type),
                    failed_marginal=nothing,
                )
            end
        last_factors = normalization.factors
        last_attempt_status = normalization.status
        last_balancing_iterations = normalization.iterations
        last_residual_history = normalization.residual_history
        last_minimum_eigenvalue = normalization.minimum_eigenvalue
        last_maximum_condition = normalization.maximum_condition
        last_failed_marginal = normalization.failed_marginal
        normalization.status === :work_limit && begin
            terminal_status = :work_limit
            break
        end
        normalization.status === :converged || continue

        residuals = _random_superoperator_residuals(
            last_factors, input_dimension, output_dimension, proportional_factor, real_type
        )
        trace_ok =
            !requested_trace_preserving || residuals.trace_preservation <= trace_tolerance
        unital_ok = !requested_unital || residuals.unitality <= unital_tolerance
        proportional_ok =
            !requested_proportional_unital ||
            residuals.proportional_unitality <= unital_tolerance
        if trace_ok && unital_ok && proportional_ok
            success = true
            terminal_status = :success
            last_attempt_status = :success
            break
        end
        last_attempt_status = :residual_check_failed
    end

    last_factors === nothing && throw(
        ArgumentError(
            "the RNG did not produce a finite nonzero draw within max_attempts=$checked_attempts",
        ),
    )
    kraus = KrausRepresentation(last_factors)
    converted = _random_superoperator_convert(kraus, selected_representation)
    residuals = _random_superoperator_residuals(
        last_factors, input_dimension, output_dimension, proportional_factor, real_type
    )
    singular_values, rank_threshold, numerical_rank = _random_superoperator_rank_diagnostic(
        last_factors, checked_atol, checked_rtol, real_type
    )
    result = RandomSuperoperatorResult(
        converted,
        kraus,
        terminal_status,
        last_attempt_status,
        success,
        last_attempt,
        input_dimension,
        output_dimension,
        rank_bound,
        numerical_rank,
        real_output,
        requested_trace_preserving,
        requested_unital,
        requested_proportional_unital,
        true,
        success && requested_trace_preserving,
        success && requested_unital,
        success && requested_proportional_unital,
        residuals.choi_trace,
        residuals.trace_preservation,
        residuals.unitality,
        proportional_factor,
        residuals.proportional_unitality,
        trace_tolerance,
        unital_tolerance,
        singular_values,
        rank_threshold,
        last_balancing_iterations,
        last_residual_history,
        last_minimum_eigenvalue,
        last_maximum_condition,
        last_failed_marginal,
        used_work[],
        _random_superoperator_message(terminal_status, last_attempt_status),
    )
    return result
end
