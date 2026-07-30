# Source-informed independent Julia implementation based on the specification
# and numeric branches of QETLAB InducedSchattenNorm.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.

using Random: AbstractRNG, randn

"""
    InducedSchattenNormResult

Auditable result of [`induced_schatten_lower_bound`](@ref). `value` is the
Schatten `q`-norm of `map(witness)`, while `witness_input_norm` records the
Schatten `p`-norm of the witness. `bound_kind == :exact` occurs only for the
proved Frobenius-to-Frobenius branch; every alternating result remains
`:lower_bound`, including a converged one.

`status` is one of `:exact`, `:converged_lower_bound`, `:iteration_limit`,
`:work_limit`, or `:stationary_zero`. The last status means that the chosen
start reached a zero gradient; it is valid lower-bound evidence, not a proof
that the map is zero. `objective_history` contains the initial witnessed value
and every completed update.
"""
struct InducedSchattenNormResult{
    R<:AbstractFloat,P<:Real,Q<:Real,M<:AbstractMatrix,V<:AbstractVector{R}
}
    value::R
    p::P
    q::Q
    bound_kind::Symbol
    status::Symbol
    exact::Bool
    converged::Bool
    witness::M
    witness_input_norm::R
    witness_output_norm::R
    normalization_residual::R
    value_residual::R
    iteration_residual::Union{Nothing,R}
    tolerance::R
    iterations::Int
    work_used::BigInt
    max_iterations::Int
    max_work::Union{Nothing,BigInt}
    objective_history::V
    message::String
end

function Base.show(io::IO, result::InducedSchattenNormResult)
    return print(
        io,
        "InducedSchattenNormResult(value=",
        result.value,
        ", p=",
        result.p,
        ", q=",
        result.q,
        ", bound_kind=",
        result.bound_kind,
        ", status=",
        result.status,
        ", iterations=",
        result.iterations,
        ")",
    )
end

function _induced_schatten_dense_limit(value)
    value === nothing && return nothing
    value isa Bool &&
        throw(ArgumentError("max_dense_entries must be a positive integer or nothing"))
    value isa Integer ||
        throw(ArgumentError("max_dense_entries must be a positive integer or nothing"))
    value > 0 ||
        throw(ArgumentError("max_dense_entries must be positive; got $(repr(value))"))
    return BigInt(value)
end

function _induced_schatten_check_dense_entries(
    rows::Int, columns::Int, limit, context::AbstractString
)
    entries = BigInt(rows) * columns
    if limit !== nothing && entries > limit
        throw(
            ArgumentError(
                "$context requires $entries dense entries, exceeding " *
                "max_dense_entries=$limit",
            ),
        )
    end
    return entries
end

function _induced_schatten_dense_matrix(
    matrix::AbstractMatrix, allow_densify::Bool, max_dense_entries, context::AbstractString
)
    _induced_schatten_check_dense_entries(
        size(matrix, 1), size(matrix, 2), max_dense_entries, context
    )
    if issparse(matrix) && !allow_densify
        throw(
            ArgumentError(
                "$context requires a dense SVD; pass `allow_densify=true` " *
                "to permit converting this $(size(matrix)) sparse matrix",
            ),
        )
    end
    return Matrix(matrix)
end

function _induced_schatten_map_dimensions(map::AbstractMapRepresentation)
    input_rows, input_columns = input_size(map)
    output_rows, output_columns = output_size(map)
    input_rows == input_columns || throw(
        ArgumentError(
            "induced Schatten norms require a square input matrix algebra; " *
            "got input_size=$(input_size(map))",
        ),
    )
    output_rows == output_columns || throw(
        ArgumentError(
            "induced Schatten norms require a square output matrix algebra; " *
            "got output_size=$(output_size(map))",
        ),
    )
    return input_rows, output_rows
end

function _induced_schatten_random_start(
    rng::AbstractRNG, map::AbstractMapRepresentation, dimension::Int, ::Type{R}
) where {R<:AbstractFloat}
    if eltype(map) <: Real
        return randn(rng, R, dimension, dimension)
    end
    return complex.(
        randn(rng, R, dimension, dimension), randn(rng, R, dimension, dimension)
    )
end

function _induced_schatten_start(
    rng::AbstractRNG,
    map::AbstractMapRepresentation,
    dimension::Int,
    p,
    initial_matrix,
    ::Type{R},
    allow_densify::Bool,
    max_dense_entries,
) where {R<:AbstractFloat}
    matrix = if initial_matrix === nothing
        _induced_schatten_random_start(rng, map, dimension, R)
    else
        initial_matrix isa AbstractMatrix{<:Number} ||
            throw(ArgumentError("initial_matrix must be a numeric matrix or nothing"))
        Base.require_one_based_indexing(initial_matrix)
        size(initial_matrix) == (dimension, dimension) || throw(
            DimensionMismatch(
                "initial_matrix must have size ($dimension, $dimension); " *
                "got $(size(initial_matrix))",
            ),
        )
        all(isfinite, initial_matrix) ||
            throw(ArgumentError("initial_matrix must contain only finite entries"))
        initial_real_type = typeof(real(zero(eltype(initial_matrix))))
        initial_real_type <: AbstractFloat || throw(
            ArgumentError(
                "initial_matrix must use real or complex floating-point entries; " *
                "got eltype $(eltype(initial_matrix))",
            ),
        )
        initial_real_type === R || throw(
            ArgumentError(
                "initial_matrix has real precision $initial_real_type but the map " *
                "uses $R; no precision-changing conversion is performed",
            ),
        )
        if issparse(initial_matrix)
            _induced_schatten_dense_matrix(
                initial_matrix,
                allow_densify,
                max_dense_entries,
                "initial-matrix normalization",
            )
        else
            copy(initial_matrix)
        end
    end

    _induced_schatten_check_dense_entries(
        dimension, dimension, max_dense_entries, "initial-matrix normalization"
    )
    input_norm = schatten_norm(matrix, p; allow_densify=allow_densify)
    isfinite(input_norm) ||
        throw(ArgumentError("initial_matrix has a non-finite Schatten p-norm"))
    iszero(input_norm) &&
        throw(DomainError(input_norm, "initial_matrix must have positive Schatten p-norm"))
    return matrix / input_norm
end

function _induced_schatten_support(
    matrix::AbstractMatrix, order, side::Symbol; allow_densify::Bool, max_dense_entries
)
    context =
        side === :output ? "output supporting-functional update" : "input witness update"
    dense = _induced_schatten_dense_matrix(
        matrix, allow_densify, max_dense_entries, context
    )
    decomposition = svd(dense; full=false)
    singular_values = decomposition.S
    scale = maximum(singular_values; init=zero(eltype(singular_values)))
    iszero(scale) && return nothing

    weights = if side === :output
        if order == one(order)
            map(value -> iszero(value) ? zero(value) : one(value), singular_values)
        elseif isinf(order)
            result = zeros(eltype(singular_values), length(singular_values))
            result[firstindex(result)] = one(eltype(result))
            result
        else
            scaled = (singular_values ./ scale) .^ (order - one(order))
            dual_order = order / (order - one(order))
            scaled / norm(scaled, dual_order)
        end
    elseif side === :input
        if order == one(order)
            result = zeros(eltype(singular_values), length(singular_values))
            result[firstindex(result)] = one(eltype(result))
            result
        elseif isinf(order)
            map(value -> iszero(value) ? zero(value) : one(value), singular_values)
        else
            scaled = (singular_values ./ scale) .^ inv(order - one(order))
            scaled / norm(scaled, order)
        end
    else
        error("internal induced Schatten support side must be :input or :output")
    end
    return decomposition.U * Diagonal(weights) * decomposition.Vt
end

function _induced_schatten_result(
    value,
    p,
    q,
    bound_kind::Symbol,
    status::Symbol,
    witness::AbstractMatrix,
    witness_output_norm,
    iteration_residual,
    tolerance,
    iterations::Int,
    work_used::BigInt,
    max_iterations::Int,
    max_work,
    objective_history::AbstractVector,
    message::String,
)
    real_type = typeof(value)
    input_norm = convert(real_type, schatten_norm(witness, p))
    output_norm = convert(real_type, witness_output_norm)
    checked_value = convert(real_type, value)
    history = convert(Vector{real_type}, objective_history)
    return InducedSchattenNormResult(
        checked_value,
        p,
        q,
        bound_kind,
        status,
        bound_kind === :exact,
        status in (:exact, :converged_lower_bound),
        witness,
        input_norm,
        output_norm,
        abs(input_norm - one(real_type)),
        abs(output_norm - checked_value),
        iteration_residual,
        tolerance,
        iterations,
        work_used,
        max_iterations,
        max_work,
        history,
        message,
    )
end

function _induced_schatten_exact_two(
    map::AbstractMapRepresentation,
    p,
    q,
    tolerance,
    max_iterations::Int,
    max_work,
    allow_densify::Bool,
    max_dense_entries,
)
    transfer = superoperator_matrix(map)
    rows, columns = size(transfer)
    entries = _induced_schatten_check_dense_entries(
        rows, columns, max_dense_entries, "the exact 2-to-2 transfer-matrix SVD"
    )
    work = entries * max(1, min(rows, columns)) + entries + rows + columns
    _tierd_induced_require_work(
        work, max_work, "the exact 2-to-2 induced Schatten norm branch"
    )
    dense = _induced_schatten_dense_matrix(
        transfer, allow_densify, max_dense_entries, "the exact 2-to-2 transfer-matrix SVD"
    )
    decomposition = svd(dense; full=false)
    value = first(decomposition.S)
    input_dimension = input_size(map)[1]
    witness = reshape(
        Vector(adjoint(decomposition.Vt)[:, 1]), input_dimension, input_dimension
    )
    output_norm = schatten_norm(apply_channel(witness, map), q; allow_densify=allow_densify)
    return _induced_schatten_result(
        value,
        p,
        q,
        :exact,
        :exact,
        witness,
        output_norm,
        nothing,
        tolerance,
        0,
        work,
        max_iterations,
        max_work,
        [value],
        "the largest singular value of the transfer matrix is the exact " *
        "Frobenius-to-Frobenius norm",
    )
end

"""
    induced_schatten_lower_bound(
        rng, map, p;
        q=p,
        initial_matrix=nothing,
        tolerance=nothing,
        max_iterations=1_000,
        max_work=1_000_000_000,
        max_dense_entries=1_000_000,
        allow_densify=false,
    ) -> InducedSchattenNormResult

Compute the exact induced Schatten `2 -> 2` norm, or a witnessed lower bound
for every other pair of orders in `[1, Inf]`. `map` is any package
[`AbstractMapRepresentation`](@ref) between square matrix algebras. A
collection of Kraus matrices is accepted by a convenience method.

The lower-bound branch is the pinned alternating Hölder-equality iteration,
with explicit corrections to its process semantics: `rng::AbstractRNG` is
mandatory, `max_iterations` and `max_work` are deterministic limits, and a
zero-gradient start returns `:stationary_zero` instead of looping or dividing
by zero. Supplying `initial_matrix` avoids all RNG draws. No global random
stream is read or mutated.

Every update uses dense SVDs of one input- and one output-sized matrix.
Sparse matrices therefore require `allow_densify=true`.
`max_dense_entries` guards each dense matrix and the exact transfer matrix
before conversion; `nothing` explicitly disables that guard. `max_work`
counts a documented deterministic upper estimate based on matrix sizes and
map actions. Reaching a limit retains the best validated witness.

`converged_lower_bound` means only that the objective change met `tolerance`;
it never becomes an exact norm certificate. Inputs, map data, and witnesses
are never normalized except for the explicitly supplied/generated working
witness, and are never symmetrized, clipped, or projected.
"""
function induced_schatten_lower_bound(
    rng::AbstractRNG,
    map::AbstractMapRepresentation,
    p;
    q=p,
    initial_matrix=nothing,
    tolerance=nothing,
    max_iterations=_TIERD_INDUCED_DEFAULT_MAX_ITERATIONS,
    max_work=_TIERD_INDUCED_DEFAULT_MAX_WORK,
    max_dense_entries=1_000_000,
    allow_densify::Bool=false,
)
    input_dimension, output_dimension = _induced_schatten_map_dimensions(map)
    eltype(map) <: LinearAlgebra.BlasFloat || throw(
        ArgumentError(
            "induced_schatten_lower_bound requires Float32, Float64, " *
            "ComplexF32, or ComplexF64 map data; got eltype $(eltype(map)). " *
            "No implicit precision-changing conversion is performed.",
        ),
    )
    real_type = typeof(real(zero(eltype(map))))
    checked_p = _tierd_induced_order(p, "p", real_type)
    checked_q = _tierd_induced_order(q, "q", real_type)
    checked_tolerance = _tierd_induced_tolerance(tolerance, real_type)
    checked_iterations = _tierd_induced_iterations(max_iterations)
    checked_work = _tierd_induced_work_limit(max_work)
    checked_dense_entries = _induced_schatten_dense_limit(max_dense_entries)

    if checked_p == real_type(2) && checked_q == real_type(2)
        return _induced_schatten_exact_two(
            map,
            checked_p,
            checked_q,
            checked_tolerance,
            checked_iterations,
            checked_work,
            allow_densify,
            checked_dense_entries,
        )
    end

    input_entries = _induced_schatten_check_dense_entries(
        input_dimension, input_dimension, checked_dense_entries, "the input witness SVD"
    )
    output_entries = _induced_schatten_check_dense_entries(
        output_dimension,
        output_dimension,
        checked_dense_entries,
        "the output supporting-functional SVD",
    )
    map_action_work =
        BigInt(input_dimension)^2 * output_dimension^2 +
        BigInt(output_dimension)^2 * input_dimension^2
    initialization_work = input_entries + map_action_work ÷ 2 + output_entries
    _tierd_induced_require_work(
        initialization_work, checked_work, "the induced Schatten initialization"
    )

    current_witness = _induced_schatten_start(
        rng,
        map,
        input_dimension,
        checked_p,
        initial_matrix,
        real_type,
        allow_densify,
        checked_dense_entries,
    )
    dual_map = dual_channel(map)
    current_output = apply_channel(current_witness, map)
    current_value = convert(
        real_type, schatten_norm(current_output, checked_q; allow_densify=allow_densify)
    )
    best_witness = copy(current_witness)
    best_value = current_value
    history = real_type[current_value]
    work_used = initialization_work
    iteration_work =
        input_entries * max(1, input_dimension) +
        output_entries * max(1, output_dimension) +
        map_action_work +
        input_entries +
        output_entries
    last_residual = nothing
    completed_iterations = 0

    while completed_iterations < checked_iterations
        if checked_work !== nothing && work_used + iteration_work > checked_work
            return _induced_schatten_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :work_limit,
                best_witness,
                best_value,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                history,
                "the deterministic work budget was reached; the witnessed " *
                "value remains a lower bound",
            )
        end

        output_support = _induced_schatten_support(
            current_output,
            checked_q,
            :output;
            allow_densify=allow_densify,
            max_dense_entries=checked_dense_entries,
        )
        if output_support === nothing
            return _induced_schatten_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :stationary_zero,
                best_witness,
                current_value,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                history,
                "the chosen witness produced a zero output and no alternating " *
                "gradient; this is not proof that the map norm is zero",
            )
        end
        input_gradient = apply_channel(output_support, dual_map)
        new_witness = _induced_schatten_support(
            input_gradient,
            checked_p,
            :input;
            allow_densify=allow_densify,
            max_dense_entries=checked_dense_entries,
        )
        if new_witness === nothing
            return _induced_schatten_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :stationary_zero,
                best_witness,
                current_value,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                history,
                "the dual update produced a zero gradient; the current value " *
                "remains only a lower bound",
            )
        end

        new_output = apply_channel(new_witness, map)
        new_value = convert(
            real_type, schatten_norm(new_output, checked_q; allow_densify=allow_densify)
        )
        isfinite(new_value) || throw(
            DomainError(
                new_value,
                "the induced Schatten iteration produced a non-finite objective",
            ),
        )
        completed_iterations += 1
        work_used += iteration_work
        last_residual = abs(new_value - current_value)
        push!(history, new_value)
        if new_value > best_value
            best_witness = copy(new_witness)
            best_value = new_value
        end
        current_witness = new_witness
        current_output = new_output
        current_value = new_value

        if last_residual <= checked_tolerance
            return _induced_schatten_result(
                best_value,
                checked_p,
                checked_q,
                :lower_bound,
                :converged_lower_bound,
                best_witness,
                best_value,
                last_residual,
                checked_tolerance,
                completed_iterations,
                work_used,
                checked_iterations,
                checked_work,
                history,
                "the alternating objective change met tolerance; the witnessed " *
                "value is still only a lower bound",
            )
        end
    end

    return _induced_schatten_result(
        best_value,
        checked_p,
        checked_q,
        :lower_bound,
        :iteration_limit,
        best_witness,
        best_value,
        last_residual,
        checked_tolerance,
        completed_iterations,
        work_used,
        checked_iterations,
        checked_work,
        history,
        "the iteration budget was reached; the witnessed value remains a " *
        "lower bound, not an exact norm",
    )
end

function induced_schatten_lower_bound(
    rng::AbstractRNG, operators::Union{Tuple,AbstractVector}, p; kwargs...
)
    return induced_schatten_lower_bound(rng, KrausRepresentation(operators), p; kwargs...)
end
