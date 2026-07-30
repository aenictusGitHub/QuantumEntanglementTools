# Source-informed independent Julia implementation based on the specification in
# QETLAB RobkCohValue.m at d8589610f00cff106537268dee2e2a1153f3a601
# (source SHA-256 99f9eaf6c0f87ee4d72bbe50aa6dbc04d824b75fda431c4417de3e5e7fe8a7c5).
# QETLAB: Copyright 2014 Nathaniel Johnston and named coauthors,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# Mathematical specification: Theorem 1 of N. Johnston, C.-K. Li,
# S. Plosker, Y.-T. Poon, and B. Regula, Phys. Rev. A 98, 022328 (2018),
# https://doi.org/10.1103/PhysRevA.98.022328

"""
    PureKCoherenceRobustnessResult

Auditable result returned by [`pure_k_coherence_robustness`](@ref).

`value` is the common standard and generalized robustness of `k`-coherence
for the validated pure state. `branch_index` is the theorem's index `ell`.
The theorem selects the *largest* admissible index, including an equality
because its comparison is non-strict.

`branch_status` is:

- `:stable` when every adjacent branch gap exceeds `branch_tolerance`;
- `:exact_equality` when the selected comparison is exactly equal in the
  input arithmetic; or
- `:near_boundary` when an adjacent gap is nonzero but no larger than
  `branch_tolerance`.

`selected_gap` is the nonnegative margin by which branch `ell >= 2` satisfies
its comparison, and is `nothing` for branch `1`. `next_gap` is the positive
margin by which branch `ell + 1` fails, and is `nothing` for branch `k`.
These diagnostics expose branch sensitivity without changing the exact
comparison used to select `branch_index`.
"""
struct PureKCoherenceRobustnessResult{R<:AbstractFloat}
    value::R
    k::Int
    branch_index::Int
    branch_status::Symbol
    selected_gap::Union{Nothing,R}
    next_gap::Union{Nothing,R}
    branch_tolerance::R
    tail_sum::R
    tail_average::R
    norm_squared::R
    normalization_residual::R
    normalization_tolerance::R

    function PureKCoherenceRobustnessResult(
        value::R,
        k::Int,
        branch_index::Int,
        branch_status::Symbol,
        selected_gap::Union{Nothing,R},
        next_gap::Union{Nothing,R},
        branch_tolerance::R,
        tail_sum::R,
        tail_average::R,
        norm_squared::R,
        normalization_residual::R,
        normalization_tolerance::R,
    ) where {R<:AbstractFloat}
        isfinite(value) && value >= zero(R) ||
            throw(ArgumentError("robustness value must be finite and nonnegative"))
        2 <= k || throw(ArgumentError("k must be at least 2"))
        1 <= branch_index <= k || throw(ArgumentError("branch_index must lie in 1:k"))
        branch_status in (:stable, :exact_equality, :near_boundary) ||
            throw(ArgumentError("invalid branch_status $(repr(branch_status))"))
        (selected_gap === nothing) == (branch_index == 1) ||
            throw(ArgumentError("selected_gap is inconsistent with branch_index"))
        (next_gap === nothing) == (branch_index == k) ||
            throw(ArgumentError("next_gap is inconsistent with branch_index"))
        selected_gap === nothing ||
            (isfinite(selected_gap) && selected_gap >= zero(R)) ||
            throw(ArgumentError("selected_gap must be finite and nonnegative"))
        next_gap === nothing ||
            (isfinite(next_gap) && next_gap > zero(R)) ||
            throw(ArgumentError("next_gap must be finite and positive"))
        all(
            isfinite,
            (
                branch_tolerance,
                tail_sum,
                tail_average,
                norm_squared,
                normalization_residual,
                normalization_tolerance,
            ),
        ) || throw(ArgumentError("result diagnostics must be finite"))
        branch_tolerance >= zero(R) ||
            throw(ArgumentError("branch_tolerance must be nonnegative"))
        tail_sum >= zero(R) || throw(ArgumentError("tail_sum must be nonnegative"))
        tail_average >= zero(R) || throw(ArgumentError("tail_average must be nonnegative"))
        norm_squared >= zero(R) || throw(ArgumentError("norm_squared must be nonnegative"))
        normalization_residual >= zero(R) ||
            throw(ArgumentError("normalization_residual must be nonnegative"))
        normalization_tolerance >= zero(R) ||
            throw(ArgumentError("normalization_tolerance must be nonnegative"))
        return new{R}(
            value,
            k,
            branch_index,
            branch_status,
            selected_gap,
            next_gap,
            branch_tolerance,
            tail_sum,
            tail_average,
            norm_squared,
            normalization_residual,
            normalization_tolerance,
        )
    end
end

function Base.show(io::IO, result::PureKCoherenceRobustnessResult)
    return print(
        io,
        "PureKCoherenceRobustnessResult(value=",
        result.value,
        ", k=",
        result.k,
        ", branch_index=",
        result.branch_index,
        ", branch_status=",
        result.branch_status,
        ")",
    )
end

function _robk_real_type(::Type{T}) where {T<:Number}
    real_type = typeof(abs(zero(T)))
    real_type <: AbstractFloat || throw(
        ArgumentError(
            "pure_k_coherence_robustness requires a real or complex " *
            "floating-point state; got element type $T. No implicit " *
            "precision-changing conversion is performed.",
        ),
    )
    return real_type
end

function _robk_tolerance(value, default::R, name::AbstractString) where {R<:AbstractFloat}
    value === nothing && return default
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be a nonnegative finite real or `nothing`"))
    isfinite(value) || throw(ArgumentError("$name must be finite; got $(repr(value))"))
    value >= zero(value) ||
        throw(ArgumentError("$name must be nonnegative; got $(repr(value))"))
    converted = try
        convert(R, value)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("$name=$(repr(value)) cannot be represented as $R"))
    end
    isfinite(converted) ||
        throw(ArgumentError("$name=$(repr(value)) overflows the state precision $R"))
    return converted
end

function _robk_validate_state_and_tolerances(state::AbstractVector{<:Number}, atol, rtol)
    !isempty(state) || throw(ArgumentError("state must have positive length"))
    Base.require_one_based_indexing(state)
    all(isfinite, state) || throw(ArgumentError("state must contain only finite entries"))

    real_type = _robk_real_type(eltype(state))
    absolute = _robk_tolerance(atol, zero(real_type), "atol")
    relative = _robk_tolerance(rtol, sqrt(eps(real_type)), "rtol")
    norm_squared = sum(abs2, state)
    isfinite(norm_squared) ||
        throw(ArgumentError("the squared state norm overflowed the input precision"))
    normalization_tolerance = absolute + relative * max(one(real_type), abs(norm_squared))
    isfinite(normalization_tolerance) ||
        throw(ArgumentError("the requested normalization tolerance overflowed"))
    normalization_residual = abs(norm_squared - one(real_type))
    normalization_residual <= normalization_tolerance || throw(
        ArgumentError(
            "state is not normalized within atol=$absolute and rtol=$relative; " *
            "its squared norm is $norm_squared. The input is never normalized " *
            "implicitly.",
        ),
    )
    return (
        real_type=real_type,
        atol=absolute,
        rtol=relative,
        norm_squared=norm_squared,
        normalization_residual=normalization_residual,
        normalization_tolerance=normalization_tolerance,
    )
end

function _robk_validate_k(k, dimension::Int)
    k isa Integer && !(k isa Bool) ||
        throw(ArgumentError("k must be an integer in 2:length(state)"))
    2 <= k <= dimension || throw(
        DomainError(
            k,
            "k must lie in 2:length(state) = 2:$dimension, as required by " *
            "the pure-state k-coherence theorem",
        ),
    )
    return Int(k)
end

function _robk_suffix_sums(magnitudes::Vector{R}) where {R<:AbstractFloat}
    dimension = length(magnitudes)
    sums = Vector{R}(undef, dimension + 1)
    sums[dimension + 1] = zero(R)
    for index in dimension:-1:1
        sums[index] = magnitudes[index] + sums[index + 1]
        isfinite(sums[index]) ||
            throw(ArgumentError("a suffix sum overflowed the input precision"))
    end
    return sums
end

function _robk_branch(
    magnitudes::Vector{R}, suffix_sums::Vector{R}, k::Int
) where {R<:AbstractFloat}
    branch_index = 1
    for candidate in k:-1:2
        denominator = k - candidate + 1
        average = suffix_sums[candidate] / denominator
        if magnitudes[candidate - 1] >= average
            branch_index = candidate
            break
        end
    end

    selected_gap = if branch_index == 1
        nothing
    else
        denominator = k - branch_index + 1
        magnitudes[branch_index - 1] - suffix_sums[branch_index] / denominator
    end
    next_gap = if branch_index == k
        nothing
    else
        next_candidate = branch_index + 1
        denominator = k - next_candidate + 1
        suffix_sums[next_candidate] / denominator - magnitudes[next_candidate - 1]
    end
    return branch_index, selected_gap, next_gap
end

"""
    pure_k_coherence_robustness(state, k; atol=nothing, rtol=nothing)

Compute the common standard and generalized robustness of `k`-coherence for a
pure state in the computational basis, using Theorem 1 of Johnston *et al.*,
Phys. Rev. A **98**, 022328 (2018).

`state` must be a nonempty, finite, one-based real or complex floating-point
vector, and its squared Euclidean norm must equal one within the requested
tolerances. It is validated but never normalized. `k` must be an integer in
`2:length(state)`.

The theorem is evaluated on coefficient magnitudes sorted in descending order.
Consequently, complex phases and coordinate permutations do not alter the
answer. If `a[1] >= ... >= a[n] >= 0`, `s[j] = sum(a[j:n])`, and `ell` is the
largest index in `2:k` for which
`a[ell - 1] >= s[ell] / (k - ell + 1)` (or `ell = 1` if none exists), then

```math
R_k = \\frac{s_{\\ell}^2}{k - \\ell + 1}
      - \\sum_{i=\\ell}^{n} a_i^2.
```

The non-strict theorem comparison determines `result.branch_index` exactly in
the input arithmetic. `atol` and `rtol` validate normalization and diagnose
whether that branch is perturbation-stable; they never move an input between
branches. See [`PureKCoherenceRobustnessResult`](@ref) for the branch
diagnostics.

The implementation uses ascending-magnitude suffix accumulation after sorting
and a cancellation-resistant nonnegative sum for the final value. It costs
`O(n log n)` time and `O(n)` workspace. Integer and rational inputs are
rejected rather than silently converted; `Float32`, `Float64`, and
`BigFloat` precision is preserved.
"""
function pure_k_coherence_robustness(
    state::AbstractVector{<:Number}, k; atol=nothing, rtol=nothing
)
    validation = _robk_validate_state_and_tolerances(state, atol, rtol)
    dimension = length(state)
    checked_k = _robk_validate_k(k, dimension)
    real_type = validation.real_type

    # A dense magnitude vector is intentional: every coordinate, including an
    # implicit sparse zero, participates in the theorem's global ordering.
    magnitudes = Vector{real_type}(undef, dimension)
    for index in eachindex(state)
        magnitudes[index] = abs(state[index])
    end
    sort!(magnitudes; rev=true)
    suffix_sums = _robk_suffix_sums(magnitudes)
    branch_index, selected_gap, next_gap = _robk_branch(magnitudes, suffix_sums, checked_k)

    denominator = checked_k - branch_index + 1
    tail_sum = suffix_sums[branch_index]
    tail_average = tail_sum / denominator

    # This is algebraically s_ell^2/(k-ell+1) - sum(a_i^2). For the selected
    # branch every summand is nonnegative in exact arithmetic, avoiding the
    # cancellation in a direct subtraction near zero robustness.
    value = zero(real_type)
    for index in branch_index:dimension
        difference = tail_average - magnitudes[index]
        difference >= zero(real_type) || throw(
            ArithmeticError(
                "branch arithmetic became inconsistent at the input precision; " *
                "retry with higher-precision state entries",
            ),
        )
        value += magnitudes[index] * difference
    end
    isfinite(value) ||
        throw(ArithmeticError("the robustness value overflowed the input precision"))

    branch_scale = max(
        one(real_type), maximum(magnitudes; init=zero(real_type)), tail_average
    )
    branch_tolerance = validation.atol + validation.rtol * branch_scale
    isfinite(branch_tolerance) ||
        throw(ArgumentError("the requested branch tolerance overflowed"))
    exact_equality = selected_gap !== nothing && iszero(selected_gap)
    near_boundary =
        (selected_gap !== nothing && selected_gap <= branch_tolerance) ||
        (next_gap !== nothing && next_gap <= branch_tolerance)
    branch_status =
        exact_equality ? :exact_equality : (near_boundary ? :near_boundary : :stable)

    return PureKCoherenceRobustnessResult(
        value,
        checked_k,
        branch_index,
        branch_status,
        selected_gap,
        next_gap,
        branch_tolerance,
        tail_sum,
        tail_average,
        validation.norm_squared,
        validation.normalization_residual,
        validation.normalization_tolerance,
    )
end
