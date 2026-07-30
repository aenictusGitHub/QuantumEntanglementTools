# Source-informed independent Julia implementations based on the executable
# contracts of QETLAB LocalDistinguishability.m, IsSeparable.m, and
# UPBSepDistinguishable.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston and Alessandro Cosentino,
# BSD-2-Clause. Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The hierarchy and UPB formulations were independently checked against:
# A. Cosentino, arXiv:1205.1031;
# S. Bandyopadhyay et al., arXiv:1408.6981; and
# A. C. Doherty, P. A. Parrilo, and F. M. Spedalieri,
# arXiv:quant-ph/0308032.

using Random

"""
    LocalDistinguishabilityStatus

Status of a separable-measurement outer-hierarchy calculation. A successful
SDP is explicitly a relaxation result: only its dual upper bound is also an
upper bound on separable distinguishability. Its primal measurement need not
be separable.
"""
@enum LocalDistinguishabilityStatus::UInt8 begin
    LocalDistinguishabilityTrivialOptimal
    LocalDistinguishabilitySolverOptimal
    LocalDistinguishabilitySolverFeasible
    LocalDistinguishabilityNumericalBoundary
    LocalDistinguishabilityResourceLimit
    LocalDistinguishabilityBackendUnavailable
    LocalDistinguishabilityBackendFailure
    LocalDistinguishabilityInvalidCertificate
end

"""
    LocalDistinguishabilityProblem

Solver-neutral outer-hierarchy SDP for minimum-error discrimination. Each
measurement effect is the first-copy marginal of a positive, optionally PPT,
permutation-invariant or bosonic `order`-copy extension. `program` contains no
JuMP or optimizer objects.
"""
struct LocalDistinguishabilityProblem{T<:AbstractFloat,S,P,G,M,E,V,B,D}
    states::S
    priors::P
    input_kind::Symbol
    dimensions::NTuple{2,Int}
    order::Int
    ppt::Bool
    bosonic::Bool
    program::G
    measurement_views::M
    extension_views::E
    variable_views::V
    bosonic_basis::B
    tolerance::T
    preflight::D

    function LocalDistinguishabilityProblem(
        token::_ValidatedConstructorToken,
        states,
        priors,
        input_kind::Symbol,
        dimensions::NTuple{2,Int},
        order::Int,
        ppt::Bool,
        bosonic::Bool,
        program,
        measurement_views,
        extension_views,
        variable_views,
        bosonic_basis,
        tolerance::T,
        preflight,
    ) where {T<:AbstractFloat}
        _require_validated_constructor_token(token)
        input_kind in (:density_matrices, :pure_columns) ||
            throw(ArgumentError("unsupported local-discrimination input kind $input_kind"))
        all(>(0), dimensions) ||
            throw(ArgumentError("local-discrimination dimensions must be positive"))
        order > 0 || throw(ArgumentError("local-discrimination order must be positive"))
        length(states) == length(priors) ||
            throw(ArgumentError("state and prior counts must agree"))
        !isempty(priors) ||
            throw(ArgumentError("a local-discrimination problem must contain a state"))
        _sepopt_represented_sum(priors) == one(eltype(priors)) || throw(
            ArgumentError("validated local-discrimination priors must sum exactly to one"),
        )
        all(prior -> isfinite(prior) && prior >= zero(prior), priors) ||
            throw(ArgumentError("validated local-discrimination priors must be finite"))
        owned_states = _entanglement_owned_evidence(Tuple(states))
        owned_priors = _entanglement_owned_evidence(Vector(priors))
        owned_measurement_views = _entanglement_owned_evidence(Tuple(measurement_views))
        owned_extension_views = _entanglement_owned_evidence(Tuple(extension_views))
        owned_variable_views = _entanglement_owned_evidence(Tuple(variable_views))
        owned_basis = _entanglement_owned_evidence(bosonic_basis)
        owned_preflight = _entanglement_owned_evidence(preflight)
        owned_program = deepcopy(program)
        return new{
            T,
            typeof(owned_states),
            typeof(owned_priors),
            typeof(owned_program),
            typeof(owned_measurement_views),
            typeof(owned_extension_views),
            typeof(owned_variable_views),
            typeof(owned_basis),
            typeof(owned_preflight),
        }(
            owned_states,
            owned_priors,
            input_kind,
            dimensions,
            order,
            ppt,
            bosonic,
            owned_program,
            owned_measurement_views,
            owned_extension_views,
            owned_variable_views,
            owned_basis,
            tolerance,
            owned_preflight,
        )
    end
end

"""
    LocalDistinguishabilityResult

Status-rich result of [`local_distinguishability`](@ref).

`separable_lower_bound` is the success probability of the always-valid
single-guess separable POVM. `separable_upper_bound` is the universal
probability-one bound unless a validated dual bound for the configured outer
hierarchy supplies a tighter value.
`relaxation_value` is the residual-checked optimum of that hierarchy when
primal and dual bounds agree. `measurement` is the hierarchy measurement and
is therefore not advertised as a separable or LOCC measurement.
"""
struct LocalDistinguishabilityResult{T<:AbstractFloat,P,Q,M,E,D,O,R}
    status::LocalDistinguishabilityStatus
    relaxation_value::Union{Nothing,T}
    relaxation_lower_bound::Union{Nothing,T}
    relaxation_upper_bound::Union{Nothing,T}
    separable_lower_bound::T
    separable_upper_bound::Union{Nothing,T}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    dimensions::NTuple{2,Int}
    order::Int
    ppt::Bool
    bosonic::Bool
    priors::P
    problem::Q
    measurement::M
    extensions::E
    dual_solution::D
    optimization_result::O
    residuals::R
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String

    function LocalDistinguishabilityResult(
        token::_ValidatedConstructorToken,
        status::LocalDistinguishabilityStatus,
        relaxation_value::Union{Nothing,T},
        relaxation_lower_bound::Union{Nothing,T},
        relaxation_upper_bound::Union{Nothing,T},
        separable_lower_bound::T,
        separable_upper_bound::Union{Nothing,T},
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        dimensions::NTuple{2,Int},
        order::Int,
        ppt::Bool,
        bosonic::Bool,
        priors,
        problem,
        measurement,
        extensions,
        dual_solution,
        optimization_result,
        residuals,
        tolerance::T,
        warnings::Tuple{Vararg{String}},
        message::AbstractString,
    ) where {T<:AbstractFloat}
        _require_validated_constructor_token(token)
        all(>(0), dimensions) ||
            throw(ArgumentError("local-discrimination dimensions must be positive"))
        order > 0 || throw(ArgumentError("local-discrimination order must be positive"))
        isfinite(tolerance) && tolerance >= zero(T) ||
            throw(ArgumentError("local-discrimination tolerance must be finite"))
        all(prior -> isfinite(prior) && prior >= zero(prior), priors) ||
            throw(ArgumentError("local-discrimination priors must be finite"))
        _sepopt_represented_sum(priors) == one(eltype(priors)) ||
            throw(ArgumentError("local-discrimination priors must sum exactly to one"))
        for (value, name) in (
            (relaxation_value, "relaxation_value"),
            (relaxation_lower_bound, "relaxation_lower_bound"),
            (relaxation_upper_bound, "relaxation_upper_bound"),
            (separable_lower_bound, "separable_lower_bound"),
            (separable_upper_bound, "separable_upper_bound"),
        )
            value === nothing && continue
            isfinite(value) ||
                throw(ArgumentError("$name must be finite when it is present"))
        end
        if status === LocalDistinguishabilityTrivialOptimal
            certified &&
            certificate_kind === :single_guess_separable_povm &&
            relaxation_value == one(T) &&
            relaxation_lower_bound == one(T) &&
            relaxation_upper_bound == one(T) &&
            separable_lower_bound == one(T) &&
            separable_upper_bound == one(T) || throw(
                ArgumentError(
                    "the trivial-optimal status requires the exact single-guess certificate and unit bounds",
                ),
            )
            measurement !== nothing &&
            dual_solution !== nothing &&
            hasproperty(residuals, :valid) &&
            residuals.valid === true || throw(
                ArgumentError(
                    "the trivial-optimal status requires validated measurement and dual evidence",
                ),
            )
        else
            !certified && certificate_kind === nothing || throw(
                ArgumentError(
                    "nontrivial local-discrimination hierarchy results are not package certificates",
                ),
            )
        end
        owned_priors = _entanglement_owned_evidence(Vector(priors))
        owned_problem = deepcopy(problem)
        owned_measurement = _entanglement_owned_evidence(measurement)
        owned_extensions = _entanglement_owned_evidence(extensions)
        owned_dual = _entanglement_owned_evidence(dual_solution)
        owned_optimization = deepcopy(optimization_result)
        owned_residuals = _entanglement_owned_evidence(residuals)
        return new{
            T,
            typeof(owned_priors),
            typeof(owned_problem),
            typeof(owned_measurement),
            typeof(owned_extensions),
            typeof(owned_dual),
            typeof(owned_optimization),
            typeof(owned_residuals),
        }(
            status,
            relaxation_value,
            relaxation_lower_bound,
            relaxation_upper_bound,
            separable_lower_bound,
            separable_upper_bound,
            certified,
            certificate_kind,
            dimensions,
            order,
            ppt,
            bosonic,
            owned_priors,
            owned_problem,
            owned_measurement,
            owned_extensions,
            owned_dual,
            owned_optimization,
            owned_residuals,
            tolerance,
            warnings,
            String(message),
        )
    end
end

function Base.show(io::IO, result::LocalDistinguishabilityResult)
    return print(
        io,
        "LocalDistinguishabilityResult(status=",
        result.status,
        ", relaxation_value=",
        result.relaxation_value,
        ", separable_bounds=(",
        result.separable_lower_bound,
        ", ",
        result.separable_upper_bound,
        "), order=",
        result.order,
        ", ppt=",
        result.ppt,
        ", bosonic=",
        result.bosonic,
        ")",
    )
end

"""
    UPBReplacementVector

One isolated product vector orthogonal to every UPB state except
`replaced_state`. `local_factors` and `global_vector` are normalized owned
vectors. `partition` records the finite span partition that produced it.
"""
struct UPBReplacementVector{F,V,P,Q,T}
    replaced_state::Int
    local_factors::F
    global_vector::V
    projection::Q
    partition::P
    orthogonality_residual::T
    rank_margins::Tuple

    function UPBReplacementVector(
        token::_ValidatedConstructorToken,
        replaced_state::Int,
        local_factors,
        global_vector,
        projection,
        partition,
        orthogonality_residual::T,
        rank_margins::Tuple,
    ) where {T}
        _require_validated_constructor_token(token)
        replaced_state > 0 || throw(ArgumentError("replaced_state must be positive"))
        !isempty(local_factors) ||
            throw(ArgumentError("a replacement vector must contain local factors"))
        all(
            factor -> factor isa AbstractVector{<:Number} && !isempty(factor), local_factors
        ) || throw(ArgumentError("replacement local factors must be nonempty vectors"))
        global_vector isa AbstractVector{<:Number} && !isempty(global_vector) ||
            throw(ArgumentError("global_vector must be a nonempty numeric vector"))
        projection isa AbstractMatrix{<:Number} &&
        size(projection) == (length(global_vector), length(global_vector)) ||
            throw(DimensionMismatch("replacement projection has incompatible dimensions"))
        all(isfinite, global_vector) && all(isfinite, projection) ||
            throw(ArgumentError("replacement evidence must be finite"))
        isfinite(orthogonality_residual) &&
        orthogonality_residual >= zero(orthogonality_residual) ||
            throw(ArgumentError("orthogonality_residual must be finite and nonnegative"))
        reconstructed = global_vector * adjoint(global_vector)
        all(iszero, projection - reconstructed) || throw(
            ArgumentError(
                "a replacement projection must equal global_vector * global_vector' exactly",
            ),
        )
        owned_local = _entanglement_owned_evidence(Tuple(local_factors))
        owned_global = _entanglement_owned_evidence(Vector(global_vector))
        owned_projection = _entanglement_owned_evidence(Matrix(projection))
        owned_partition = deepcopy(partition)
        owned_margins = deepcopy(rank_margins)
        return new{
            typeof(owned_local),
            typeof(owned_global),
            typeof(owned_partition),
            typeof(owned_projection),
            T,
        }(
            replaced_state,
            owned_local,
            owned_global,
            owned_projection,
            owned_partition,
            orthogonality_residual,
            owned_margins,
        )
    end
end

@doc raw"""
    UPBSeparableDiscriminationProblem

Solver-neutral conic feasibility problem for the Bandyopadhyay--Cosentino--
Johnston--Russo--Watrous--Yu UPB criterion:

```math
I = \sum_r \lambda_r q_r q_r^\dagger,\qquad \lambda_r\geq 0.
```
"""
struct UPBSeparableDiscriminationProblem{U,R,P,D}
    upb_analysis::U
    replacement_vectors::R
    program::P
    dimensions::D
    state_count::Int
    candidate_count::Int

    function UPBSeparableDiscriminationProblem(
        token::_ValidatedConstructorToken,
        upb_analysis,
        replacement_vectors,
        program,
        dimensions,
        state_count::Int,
        candidate_count::Int,
    )
        _require_validated_constructor_token(token)
        state_count > 0 || throw(ArgumentError("state_count must be positive"))
        candidate_count > 0 || throw(ArgumentError("candidate_count must be positive"))
        candidate_count == length(replacement_vectors) ||
            throw(ArgumentError("candidate_count must equal the replacement-vector count"))
        all(replacement -> replacement isa UPBReplacementVector, replacement_vectors) ||
            throw(ArgumentError("replacement_vectors contains an invalid entry"))
        checked_dimensions = Tuple(dimensions)
        all(dimension -> dimension isa Int && dimension > 0, checked_dimensions) ||
            throw(ArgumentError("UPB dimensions must be positive Int values"))
        _checked_product(checked_dimensions, "UPB dimensions")
        program.variable_count == candidate_count || throw(
            ArgumentError("the conic program variable count must equal candidate_count")
        )
        owned_replacements = Tuple(replacement_vectors)
        owned_program = deepcopy(program)
        owned_analysis = deepcopy(upb_analysis)
        return new{
            typeof(owned_analysis),
            typeof(owned_replacements),
            typeof(owned_program),
            typeof(checked_dimensions),
        }(
            owned_analysis,
            owned_replacements,
            owned_program,
            checked_dimensions,
            state_count,
            candidate_count,
        )
    end
end

"""
    UPBSeparableDiscriminationResult

Auditable result returned by [`upb_sep_distinguishable`](@ref).

The current public input domain and replacement generator are floating-point
only and provide no exact proof linking the computed cone to the represented
input UPB. Consequently, `separably_distinguishable` and `certificate_kind`
are always `nothing`, and `certified` is always `false`. Numerical
reconstructions and floating-cone Farkas separators remain diagnostic evidence.
"""
struct UPBSeparableDiscriminationResult{B,U,R,P,C,M,D,O,T}
    status::Symbol
    feasibility::Symbol
    separably_distinguishable::B
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    upb_analysis::U
    replacement_vectors::R
    problem::P
    coefficients::C
    reconstruction::M
    reconstruction_residual::Union{Nothing,T}
    minimum_coefficient::Union{Nothing,T}
    dual_certificate::D
    optimization_result::O
    partitions_examined::Int
    candidates_generated::Int
    work_used::BigInt
    tolerance::T
    message::String

    function UPBSeparableDiscriminationResult(
        token::_ValidatedConstructorToken,
        status::Symbol,
        feasibility::Symbol,
        separably_distinguishable,
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        upb_analysis,
        replacement_vectors,
        problem,
        coefficients,
        reconstruction,
        reconstruction_residual::Union{Nothing,T},
        minimum_coefficient::Union{Nothing,T},
        dual_certificate,
        optimization_result,
        partitions_examined::Int,
        candidates_generated::Int,
        work_used::BigInt,
        tolerance::T,
        message::AbstractString,
    ) where {T}
        _require_validated_constructor_token(token)
        partitions_examined >= 0 ||
            throw(ArgumentError("partitions_examined must be nonnegative"))
        candidates_generated >= 0 ||
            throw(ArgumentError("candidates_generated must be nonnegative"))
        candidates_generated == length(replacement_vectors) || throw(
            ArgumentError("candidates_generated must equal the replacement-vector count"),
        )
        work_used >= 0 || throw(ArgumentError("work_used must be nonnegative"))
        isfinite(tolerance) && tolerance >= zero(tolerance) ||
            throw(ArgumentError("UPB tolerance must be finite and nonnegative"))
        certified && throw(
            ArgumentError(
                "the floating-only UPB replacement API cannot issue a linked Boolean certificate",
            ),
        )
        separably_distinguishable === nothing && certificate_kind === nothing || throw(
            ArgumentError(
                "the floating-only UPB replacement API cannot carry a Boolean or certificate kind",
            ),
        )
        for (value, name) in (
            (reconstruction_residual, "reconstruction_residual"),
            (minimum_coefficient, "minimum_coefficient"),
        )
            value === nothing && continue
            isfinite(value) || throw(ArgumentError("$name must be finite when present"))
        end
        owned_analysis = deepcopy(upb_analysis)
        owned_replacements = Tuple(replacement_vectors)
        owned_problem = deepcopy(problem)
        owned_coefficients = _entanglement_owned_evidence(coefficients)
        owned_reconstruction = _entanglement_owned_evidence(reconstruction)
        owned_dual = _entanglement_owned_evidence(dual_certificate)
        owned_optimization = deepcopy(optimization_result)
        return new{
            typeof(separably_distinguishable),
            typeof(owned_analysis),
            typeof(owned_replacements),
            typeof(owned_problem),
            typeof(owned_coefficients),
            typeof(owned_reconstruction),
            typeof(owned_dual),
            typeof(owned_optimization),
            T,
        }(
            status,
            feasibility,
            separably_distinguishable,
            certified,
            certificate_kind,
            owned_analysis,
            owned_replacements,
            owned_problem,
            owned_coefficients,
            owned_reconstruction,
            reconstruction_residual,
            minimum_coefficient,
            owned_dual,
            owned_optimization,
            partitions_examined,
            candidates_generated,
            work_used,
            tolerance,
            String(message),
        )
    end
end

function _entanglement_owned_evidence(value::OperatorSchmidtDecompositionResult)
    return OperatorSchmidtDecompositionResult(
        _entanglement_owned_evidence(value.coefficients),
        _entanglement_owned_evidence(value.left_factors),
        _entanglement_owned_evidence(value.right_factors),
        value.row_dims,
        value.column_dims,
        value.factor_convention,
        value.coordinate_imaginary_residual,
        value.coordinate_imaginary_tolerance,
    )
end

function _entanglement_owned_evidence(value::SymmetricExtensionWitness)
    return SymmetricExtensionWitness(
        _entanglement_owned_evidence(value.operator),
        value.expectation,
        value.normalization_residual,
        value.dual_stationarity_residual,
        value.tolerance,
        deepcopy(value.dual_data),
        value.separator_validated,
        value.entanglement_witness,
        value.warning,
    )
end

function _entanglement_owned_evidence(value::SymmetricExtensionResult)
    return SymmetricExtensionResult(
        value.status,
        value.verdict,
        value.certificate_kind,
        value.hierarchy,
        value.dimensions,
        value.order,
        value.ppt,
        value.bosonic,
        deepcopy(value.problem),
        _entanglement_owned_evidence(value.extension),
        _entanglement_owned_evidence(value.witness),
        deepcopy(value.optimization_result),
        _entanglement_owned_evidence(value.residuals),
        value.tolerance,
        value.warnings,
        value.message,
    )
end

function _localdisc_legacy_output_valid(result::LocalDistinguishabilityResult)
    _sepopt_represented_sum(result.priors) == one(eltype(result.priors)) || return false
    result.relaxation_value === nothing && return false
    result.measurement === nothing && return false
    hasproperty(result.residuals, :valid) && result.residuals.valid === true || return false
    dual = result.dual_solution
    dual !== nothing &&
    hasproperty(dual, :completeness_operator) &&
    dual.completeness_operator !== nothing || return false
    if result.status === LocalDistinguishabilityTrivialOptimal
        return result.certified &&
               result.certificate_kind === :single_guess_separable_povm &&
               result.relaxation_value == one(result.relaxation_value)
    end
    return result.status === LocalDistinguishabilitySolverOptimal &&
           !result.certified &&
           result.certificate_kind === nothing
end

function Base.show(io::IO, result::UPBSeparableDiscriminationResult)
    return print(
        io,
        "UPBSeparableDiscriminationResult(status=",
        result.status,
        ", feasibility=",
        result.feasibility,
        ", separably_distinguishable=",
        result.separably_distinguishable,
        ", candidates=",
        result.candidates_generated,
        ")",
    )
end

mutable struct _SeparabilityWorkBudget
    used::BigInt
    limit::Union{Nothing,BigInt}
end

function _sepopt_positive_int(value, name::AbstractString)
    value isa Integer && !(value isa Bool) && value > 0 ||
        throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _sepopt_nonnegative_int(value, name::AbstractString)
    value isa Integer && !(value isa Bool) && value >= 0 ||
        throw(ArgumentError("$name must be a nonnegative integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _sepopt_limit(value, name::AbstractString)
    value === nothing && return nothing
    value isa Integer && !(value isa Bool) && value >= 0 ||
        throw(ArgumentError("$name must be a nonnegative integer or nothing"))
    return BigInt(value)
end

function _sepopt_consume!(budget::_SeparabilityWorkBudget, amount, label::AbstractString)
    amount_big = BigInt(amount)
    amount_big >= 0 || error("internal work estimate for $label is negative")
    required = budget.used + amount_big
    if budget.limit !== nothing && required > budget.limit
        return false
    end
    budget.used = required
    return true
end

function _sepopt_guard_dense(dimension::Int, max_dense_entries, label::AbstractString)
    limit = _sepopt_limit(max_dense_entries, "max_dense_entries")
    entries = BigInt(dimension)^2
    if limit !== nothing && entries > limit
        throw(
            ArgumentError(
                "$label needs $entries dense entries, exceeding max_dense_entries=$limit"
            ),
        )
    end
    return limit
end

function _sepopt_represented_sum(values)
    total = zero(eltype(values))
    for value in values
        total += value
    end
    return total
end

function _sepopt_represented_trace(matrix::AbstractMatrix)
    total = zero(eltype(matrix))
    for index in axes(matrix, 1)
        total += matrix[index, index]
    end
    return total
end

function _localdisc_validate_dimensions(dims, dimension::Int)
    layout = _tierd_bipartite_layout(dims, dimension)
    return (layout[1], layout[2])
end

function _localdisc_prepare(
    states,
    dims;
    priors,
    atol,
    rtol,
    allow_densify::Bool,
    max_dense_entries,
    max_dimension,
    max_states,
    limits::OptimizationLimits,
)
    prepared = _discrimination_prepare(
        states;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dimension=max_dimension,
        max_states=max_states,
        # The hierarchy preflight below accounts for its actual extension
        # coordinates and converts model limits into a structured result.
        max_variables=typemax(Int),
    )
    if priors === nothing &&
        _sepopt_represented_sum(prepared.priors) != one(eltype(prepared.priors))
        # The package-generated uniform prior is not caller input. Adjust only
        # its final represented entry so the finite floating vector sums
        # exactly to one in the arithmetic used by the model.
        exact_priors = copy(prepared.priors)
        exact_priors[end] =
            one(eltype(exact_priors)) -
            _sepopt_represented_sum(@view exact_priors[1:(end - 1)])
        prepared = merge(prepared, (; priors=exact_priors))
    end
    all(prepared.states) do state
        trace_value = _sepopt_represented_trace(state)
        return iszero(imag(trace_value)) && real(trace_value) == one(real(trace_value))
    end || throw(
        ArgumentError(
            "local distinguishability requires exactly represented normalized states; " *
            "a tolerance-only state or pure-column normalization is not eligible for " *
            "an analytic certificate or probability-one bound",
        ),
    )
    _sepopt_represented_sum(prepared.priors) == one(eltype(prepared.priors)) || throw(
        ArgumentError(
            "local distinguishability requires priors whose represented sum is exactly one; " *
            "tolerance-only prior normalization is not accepted",
        ),
    )
    dimensions = _localdisc_validate_dimensions(dims, prepared.dimension)
    _sepopt_guard_dense(
        prepared.dimension, max_dense_entries, "local distinguishability input"
    )
    return merge(prepared, (; dimensions))
end

function _localdisc_preflight(
    prepared, order::Int, ppt::Bool, bosonic::Bool, limits::OptimizationLimits
)
    one_effect = _symext_preflight(prepared.dimensions, order, bosonic, ppt, limits)
    count = BigInt(prepared.state_count)
    variables = count * one_effect.variable_count
    equalities = BigInt(prepared.dimension)^2 + count * one_effect.symmetry_equalities
    psd_blocks = count * one_effect.psd_blocks
    largest_psd = max(
        one_effect.variable_matrix_dimension, ppt ? one_effect.ambient_dimension : BigInt(0)
    )
    model_entries = count * one_effect.real_block_triangle_entries
    reasons = String[]
    variables > limits.max_variables && push!(
        reasons,
        "model needs $variables variables, exceeding max_variables=$(limits.max_variables)",
    )
    equalities > limits.max_equalities && push!(
        reasons,
        "model needs $equalities equalities, exceeding max_equalities=$(limits.max_equalities)",
    )
    psd_blocks > limits.max_psd_blocks && push!(
        reasons,
        "model needs $psd_blocks PSD blocks, exceeding max_psd_blocks=$(limits.max_psd_blocks)",
    )
    largest_psd > limits.max_psd_dimension && push!(
        reasons,
        "largest PSD block has dimension $largest_psd, exceeding max_psd_dimension=$(limits.max_psd_dimension)",
    )
    model_entries > limits.max_model_entries && push!(
        reasons,
        "real-block scalarization needs $model_entries entries, exceeding max_model_entries=$(limits.max_model_entries)",
    )
    return (
        allowed=isempty(reasons),
        reason=join(reasons, "; "),
        variables=variables,
        equalities=equalities,
        psd_blocks=psd_blocks,
        largest_psd=largest_psd,
        model_entries=model_entries,
        one_effect=one_effect,
    )
end

function _localdisc_known_coordinates(
    prepared,
    order::Int,
    bosonic::Bool,
    variable_dimension::Int,
    variables_per_effect::Int,
    variable_count::Int,
)
    T = eltype(prepared.priors)
    coordinates = zeros(T, variable_count)
    state_count = prepared.state_count
    dimension_b = prepared.dimensions[2]
    diagonal_value = if bosonic
        symmetric_dimension = variable_dimension ÷ prepared.dimensions[1]
        convert(T, dimension_b) /
        (convert(T, state_count) * convert(T, symmetric_dimension))
    else
        inv(convert(T, state_count) * convert(T, BigInt(dimension_b)^(order - 1)))
    end
    for state_index in 1:state_count
        first_variable = (state_index - 1) * variables_per_effect + 1
        coordinates[first_variable:(first_variable + variable_dimension - 1)] .=
            diagonal_value
    end
    return coordinates
end

function _localdisc_build_problem(
    prepared, order::Int, ppt::Bool, bosonic::Bool, limits::OptimizationLimits, preflight
)
    T = eltype(prepared.priors)
    dimension_a, dimension_b = prepared.dimensions
    ambient = Int(preflight.one_effect.ambient_dimension)
    variable_dimension = Int(preflight.one_effect.variable_matrix_dimension)
    variables_per_effect = variable_dimension^2
    variable_count = Int(preflight.variables)
    full_dimensions = (dimension_a, ntuple(_ -> dimension_b, order)...)
    trace_out = order == 1 ? () : Tuple(3:(order + 1))
    basis = bosonic ? _symext_bosonic_lift(prepared.dimensions, order, T) : nothing

    variable_views = HermitianAffineMatrix{T}[]
    extension_views = HermitianAffineMatrix{T}[]
    measurement_views = HermitianAffineMatrix{T}[]
    psd_constraints = HermitianAffineMatrix{T}[]
    equalities = AffineEquality{T}[]

    for state_index in 1:prepared.state_count
        first_variable = (state_index - 1) * variables_per_effect + 1
        variable = hermitian_variable(
            Symbol(:local_measurement_coordinate_, state_index),
            variable_dimension;
            first_variable=first_variable,
            variable_count=variable_count,
            coefficient_type=T,
        )
        extension = if bosonic
            _symext_affine_congruence(
                variable, basis, Symbol(:local_measurement_extension_, state_index)
            )
        else
            HermitianAffineMatrix(
                Symbol(:local_measurement_extension_, state_index),
                variable.constant,
                [term.variable for term in variable.terms],
                [term.coefficient for term in variable.terms],
                variable.variable_count,
            )
        end
        measurement = if isempty(trace_out)
            HermitianAffineMatrix(
                Symbol(:local_measurement_, state_index),
                extension.constant,
                [term.variable for term in extension.terms],
                [term.coefficient for term in extension.terms],
                extension.variable_count,
            )
        else
            partial_trace_affine(
                extension,
                full_dimensions;
                trace_out=trace_out,
                name=Symbol(:local_measurement_, state_index),
            )
        end
        push!(variable_views, variable)
        push!(extension_views, extension)
        push!(measurement_views, measurement)
        push!(psd_constraints, variable)

        if !bosonic
            for first_b in 1:(order - 1)
                permuted = _symext_permuted_affine(
                    extension, prepared.dimensions, order, first_b
                )
                difference = _symext_affine_difference(
                    extension,
                    permuted,
                    Symbol(:local_measurement_symmetry_, state_index, :_, first_b),
                )
                append!(
                    equalities,
                    hermitian_equalities(
                        difference;
                        name_prefix=Symbol(
                            :local_measurement_symmetry_, state_index, :_, first_b
                        ),
                    ),
                )
            end
        end
        if ppt
            for last_system in 2:(order + 1)
                push!(
                    psd_constraints,
                    partial_transpose_affine(
                        extension,
                        full_dimensions;
                        systems=Tuple(2:last_system),
                        name=Symbol(
                            :local_measurement_ppt_, state_index, :_, last_system - 1
                        ),
                    ),
                )
            end
        end
    end

    measurement_sum = _discrimination_affine_sum(measurement_views, :local_measurement_sum)
    append!(
        equalities,
        hermitian_equalities(
            measurement_sum;
            target=Matrix{Complex{T}}(I, prepared.dimension, prepared.dimension),
            name_prefix=:local_measurement_complete,
        ),
    )
    objective = _discrimination_objective(prepared, measurement_views, variable_count)
    feasible = _localdisc_known_coordinates(
        prepared, order, bosonic, variable_dimension, variables_per_effect, variable_count
    )
    views = vcat(measurement_views, extension_views, variable_views)
    program = SemidefiniteProgram(
        :local_distinguishability_outer_hierarchy,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=psd_constraints,
        primal_views=views,
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:ppt_symmetric_extendible_measurement_outer_hierarchy,
            dimensions=prepared.dimensions,
            order=order,
            ppt=ppt,
            bosonic=bosonic,
            state_count=prepared.state_count,
            ambient_dimension=ambient,
            variable_matrix_dimension=variable_dimension,
            subsystem_order=:A_then_B_copies,
            source=:cosentino_local_discrimination_hierarchy,
        ),
    )
    return LocalDistinguishabilityProblem(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        prepared.states,
        prepared.priors,
        prepared.input_kind,
        prepared.dimensions,
        order,
        ppt,
        bosonic,
        program,
        Tuple(measurement_views),
        Tuple(extension_views),
        Tuple(variable_views),
        basis,
        prepared.tolerance,
        preflight,
    )
end

"""
    local_distinguishability_problem(states, dims; kwargs...)

Build the package-owned PPT/symmetric-extension outer relaxation for local
state discrimination. `states` follows [`state_distinguishability`](@ref):
either a tuple/vector of normalized density matrices or a matrix of normalized
pure-state columns. Priors and states are checked but never normalized. Their
represented prior sum, density traces, or pure-column squared norms must equal
one exactly; tolerance-only normalization is rejected.
"""
function local_distinguishability_problem(
    states,
    dims;
    priors=nothing,
    order=2,
    ppt::Bool=true,
    bosonic::Bool=true,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_dimension=32,
    max_states=32,
    max_order=6,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _sepopt_positive_int(order, "order")
    checked_order <= _sepopt_positive_int(max_order, "max_order") ||
        throw(ArgumentError("order=$checked_order exceeds max_order=$max_order"))
    prepared = _localdisc_prepare(
        states,
        dims;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        max_dimension=max_dimension,
        max_states=max_states,
        limits=limits,
    )
    preflight = _localdisc_preflight(prepared, checked_order, ppt, bosonic, limits)
    preflight.allowed || throw(ArgumentError(preflight.reason))
    return _localdisc_build_problem(
        prepared, checked_order, ppt, bosonic, limits, preflight
    )
end

function _localdisc_result(
    status::LocalDistinguishabilityStatus,
    prepared,
    order::Int,
    ppt::Bool,
    bosonic::Bool;
    relaxation_value=nothing,
    relaxation_lower_bound=nothing,
    relaxation_upper_bound=one(eltype(prepared.priors)),
    separable_lower_bound=maximum(prepared.priors),
    separable_upper_bound=one(eltype(prepared.priors)),
    certified::Bool=false,
    certificate_kind=nothing,
    problem=nothing,
    measurement=nothing,
    extensions=nothing,
    dual_solution=nothing,
    optimization_result=nothing,
    residuals=NamedTuple(),
    warnings=(),
    message,
)
    T = eltype(prepared.priors)
    convert_optional(value) = value === nothing ? nothing : convert(T, value)
    return LocalDistinguishabilityResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        convert_optional(relaxation_value),
        convert_optional(relaxation_lower_bound),
        convert_optional(relaxation_upper_bound),
        convert(T, separable_lower_bound),
        convert_optional(separable_upper_bound),
        certified,
        certificate_kind,
        prepared.dimensions,
        order,
        ppt,
        bosonic,
        copy(prepared.priors),
        problem,
        measurement,
        extensions,
        dual_solution,
        optimization_result,
        residuals,
        convert(T, prepared.tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

function _localdisc_trivial(prepared, order::Int, ppt::Bool, bosonic::Bool)
    T = eltype(prepared.priors)
    deterministic = findfirst(eachindex(prepared.priors)) do index
        return prepared.priors[index] == one(T) && all(
            other -> other == index || iszero(prepared.priors[other]),
            eachindex(prepared.priors),
        )
    end
    prepared.state_count == 1 || deterministic !== nothing || return nothing
    chosen = prepared.state_count == 1 ? 1 : deterministic
    measurement = Tuple(
        if index == chosen
            Matrix{Complex{T}}(I, prepared.dimension, prepared.dimension)
        else
            zeros(Complex{T}, prepared.dimension, prepared.dimension)
        end for index in 1:prepared.state_count
    )
    dual_solution = (
        completeness_operator=copy(prepared.states[chosen]),
        raw_dual=nothing,
        source=:analytic_single_guess,
        message="analytic state-discrimination dual for the deterministic single-guess branch",
    )
    return _localdisc_result(
        LocalDistinguishabilityTrivialOptimal,
        prepared,
        order,
        ppt,
        bosonic;
        relaxation_value=one(T),
        relaxation_lower_bound=one(T),
        relaxation_upper_bound=one(T),
        separable_lower_bound=one(T),
        separable_upper_bound=one(T),
        certified=true,
        certificate_kind=:single_guess_separable_povm,
        measurement=measurement,
        dual_solution,
        residuals=(
            completeness_residual=zero(T),
            positivity_violation=zero(T),
            success_probability=one(T),
            valid=true,
        ),
        message="a one-state or deterministic-prior ensemble is distinguished by the separable POVM containing the identity effect",
    )
end

function _localdisc_solution_tolerance(problem, optimization, backend)
    candidate = problem.tolerance
    if backend isa JuMPBackend
        candidate = max(candidate, convert(typeof(candidate), backend.atol + backend.rtol))
    end
    for value in (
        optimization.primal_residual, optimization.dual_residual, optimization.absolute_gap
    )
        value === nothing && continue
        isfinite(value) || continue
        candidate = max(candidate, convert(typeof(candidate), value))
    end
    return max(candidate, 64 * eps(typeof(candidate)))
end

function _localdisc_validate_solution(problem, optimization, tolerance)
    primal = optimization.primal
    primal === nothing && return (
        valid=false,
        measurement=nothing,
        extensions=nothing,
        success_probability=nothing,
        completeness_residual=nothing,
        positivity_violation=nothing,
        marginal_residual=nothing,
        permutation_residual=nothing,
        bosonic_support_residual=nothing,
        ppt_violation=nothing,
        objective_residual=nothing,
    )
    coordinates = primal.coordinates
    measurements = Tuple(
        Matrix(evaluate_affine(view, coordinates)) for view in problem.measurement_views
    )
    extensions = Tuple(
        Matrix(evaluate_affine(view, coordinates)) for view in problem.extension_views
    )
    variables = Tuple(
        Matrix(evaluate_affine(view, coordinates)) for view in problem.variable_views
    )
    T = eltype(problem.priors)
    identity_matrix = Matrix{Complex{T}}(
        I, prod(problem.dimensions), prod(problem.dimensions)
    )
    total = zeros(Complex{T}, size(identity_matrix))
    positivity_violation = zero(T)
    hermiticity_residual = zero(T)
    marginal_residual = zero(T)
    permutation_residual = zero(T)
    bosonic_support_residual = zero(T)
    ppt_violation = zero(T)
    full_dimensions = (
        problem.dimensions[1], ntuple(_ -> problem.dimensions[2], problem.order)...
    )
    trace_out = problem.order == 1 ? () : Tuple(3:(problem.order + 1))
    for (measurement, extension, variable) in zip(measurements, extensions, variables)
        total .+= measurement
        hermiticity_residual = max(
            hermiticity_residual,
            maximum(abs, measurement - adjoint(measurement); init=zero(T)),
            maximum(abs, extension - adjoint(extension); init=zero(T)),
            maximum(abs, variable - adjoint(variable); init=zero(T)),
        )
        variable_minimum = minimum(real, eigvals(Hermitian(variable)))
        positivity_violation = max(positivity_violation, max(zero(T), -variable_minimum))
        numeric_marginal = if isempty(trace_out)
            extension
        else
            partial_trace(extension, full_dimensions; trace_out=trace_out)
        end
        marginal_residual = max(
            marginal_residual, maximum(abs, numeric_marginal - measurement; init=zero(T))
        )
        if problem.bosonic
            projector = problem.bosonic_basis * adjoint(problem.bosonic_basis)
            bosonic_support_residual = max(
                bosonic_support_residual,
                maximum(abs, extension - projector * extension * projector; init=zero(T)),
            )
        else
            for first_b in 1:(problem.order - 1)
                full_count = problem.order + 1
                permutation = collect(1:full_count)
                left_position = first_b + 1
                right_position = left_position + 1
                permutation[left_position], permutation[right_position] = permutation[right_position],
                permutation[left_position]
                permuted = permute_subsystems(
                    extension, SubsystemPermutationPlan(full_dimensions, Tuple(permutation))
                )
                permutation_residual = max(
                    permutation_residual, maximum(abs, extension - permuted; init=zero(T))
                )
            end
        end
        if problem.ppt
            for last_system in 2:(problem.order + 1)
                transposed = partial_transpose(
                    extension, full_dimensions; systems=Tuple(2:last_system)
                )
                minimum_value = minimum(real, eigvals(Hermitian(Matrix(transposed))))
                ppt_violation = max(ppt_violation, max(zero(T), -minimum_value))
            end
        end
    end
    completeness_residual = maximum(abs, total - identity_matrix; init=zero(T))
    success = zero(T)
    for (prior, state, measurement) in zip(problem.priors, problem.states, measurements)
        success += prior * real(dot(state, measurement))
    end
    objective_residual = if optimization.objective_value === nothing
        nothing
    else
        abs(success - optimization.objective_value)
    end
    valid =
        hermiticity_residual <= tolerance &&
        completeness_residual <= tolerance &&
        positivity_violation <= tolerance &&
        marginal_residual <= tolerance &&
        permutation_residual <= tolerance &&
        bosonic_support_residual <= tolerance &&
        ppt_violation <= tolerance &&
        success >= -tolerance &&
        success <= one(T) + tolerance &&
        (objective_residual === nothing || objective_residual <= tolerance)
    return (
        valid,
        measurement=measurements,
        extensions,
        success_probability=success,
        completeness_residual,
        positivity_violation,
        hermiticity_residual,
        marginal_residual,
        permutation_residual,
        bosonic_support_residual,
        ppt_violation,
        objective_residual,
    )
end

function _localdisc_dual_solution(problem, optimization)
    dual = optimization.dual
    dual === nothing && return nothing
    indices = findall(
        constraint -> startswith(String(constraint.name), "local_measurement_complete"),
        problem.program.equalities,
    )
    length(indices) == prod(problem.dimensions)^2 || return (
        completeness_operator=nothing,
        raw_dual=dual,
        message="the completeness dual coordinates had an unexpected length",
    )
    values = dual.equalities[indices]
    operator = _symext_hermitian_from_duals(values, prod(problem.dimensions))
    return (
        completeness_operator=operator,
        raw_dual=dual,
        message="raw hierarchy dual retained; the completeness block is not by itself an entanglement or separability certificate",
    )
end

"""
    local_distinguishability(states, dims; backend=NoOptimizationBackend(), ...)

Bound separable-measurement discrimination with a PPT/symmetric-extension
outer hierarchy. The input is never normalized or repaired, and its represented
state and prior normalization must be exact. The returned
hierarchy POVM is not claimed to be separable. A missing backend, solver
failure, or residual boundary never becomes a mathematical negative result.
"""
function local_distinguishability(
    states,
    dims;
    priors=nothing,
    order=2,
    ppt::Bool=true,
    bosonic::Bool=true,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_dimension=32,
    max_states=32,
    max_order=6,
    limits::OptimizationLimits=OptimizationLimits(),
)
    checked_order = _sepopt_positive_int(order, "order")
    checked_max_order = _sepopt_positive_int(max_order, "max_order")
    prepared = _localdisc_prepare(
        states,
        dims;
        priors=priors,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
        max_dimension=max_dimension,
        max_states=max_states,
        limits=limits,
    )
    trivial = _localdisc_trivial(prepared, checked_order, ppt, bosonic)
    trivial === nothing || return trivial
    if checked_order > checked_max_order
        return _localdisc_result(
            LocalDistinguishabilityResourceLimit,
            prepared,
            checked_order,
            ppt,
            bosonic;
            message="order=$checked_order exceeds max_order=$checked_max_order",
        )
    end
    preflight = _localdisc_preflight(prepared, checked_order, ppt, bosonic, limits)
    if !preflight.allowed
        return _localdisc_result(
            LocalDistinguishabilityResourceLimit,
            prepared,
            checked_order,
            ppt,
            bosonic;
            residuals=(preflight=preflight,),
            message=preflight.reason,
        )
    end
    problem = _localdisc_build_problem(
        prepared, checked_order, ppt, bosonic, limits, preflight
    )
    optimization = solve_optimization(problem.program, backend)
    if optimization.status === OptimizationBackendUnavailable
        return _localdisc_result(
            LocalDistinguishabilityBackendUnavailable,
            prepared,
            checked_order,
            ppt,
            bosonic;
            problem,
            optimization_result=optimization,
            message=optimization.message,
        )
    elseif optimization.status === OptimizationLimit
        return _localdisc_result(
            LocalDistinguishabilityResourceLimit,
            prepared,
            checked_order,
            ppt,
            bosonic;
            problem,
            optimization_result=optimization,
            message=optimization.message,
        )
    elseif optimization.status in (
        OptimizationMalformedBackend,
        OptimizationUnsupported,
        OptimizationNumericalFailure,
        OptimizationInconsistent,
        OptimizationUnknown,
        OptimizationInfeasible,
        OptimizationUnbounded,
    )
        return _localdisc_result(
            LocalDistinguishabilityBackendFailure,
            prepared,
            checked_order,
            ppt,
            bosonic;
            problem,
            optimization_result=optimization,
            message="the hierarchy optimization did not produce a residual-checked primal result: $(optimization.message)",
        )
    end

    tolerance = _localdisc_solution_tolerance(problem, optimization, backend)
    residuals = _localdisc_validate_solution(problem, optimization, tolerance)
    if !residuals.valid
        return _localdisc_result(
            LocalDistinguishabilityInvalidCertificate,
            prepared,
            checked_order,
            ppt,
            bosonic;
            problem,
            optimization_result=optimization,
            residuals,
            message="the optimizer's hierarchy measurement failed package-owned residual checks",
        )
    end

    lower = residuals.success_probability
    raw_upper = if optimization.objective_bound !== nothing
        optimization.objective_bound
    else
        optimization.dual_objective_value
    end
    T = eltype(prepared.priors)
    warnings = String["the returned measurement belongs to an outer hierarchy and is not necessarily separable"]
    upper = raw_upper
    if upper !== nothing && upper > one(T)
        push!(
            warnings,
            "the solver upper bound was intersected explicitly with the analytic probability bound one",
        )
        upper = one(T)
    end
    if upper !== nothing && upper < lower - tolerance
        return _localdisc_result(
            LocalDistinguishabilityInvalidCertificate,
            prepared,
            checked_order,
            ppt,
            bosonic;
            problem,
            measurement=residuals.measurement,
            extensions=residuals.extensions,
            optimization_result=optimization,
            residuals,
            warnings,
            message="the reported hierarchy upper bound is below its checked primal value",
        )
    end
    agreed =
        upper !== nothing &&
        abs(upper - lower) <= tolerance &&
        optimization.status === OptimizationOptimal
    status = if agreed
        LocalDistinguishabilitySolverOptimal
    else
        LocalDistinguishabilitySolverFeasible
    end
    return _localdisc_result(
        status,
        prepared,
        checked_order,
        ppt,
        bosonic;
        relaxation_value=agreed ? lower : nothing,
        relaxation_lower_bound=lower,
        relaxation_upper_bound=upper,
        separable_lower_bound=maximum(prepared.priors),
        separable_upper_bound=upper === nothing ? one(T) : upper,
        problem,
        measurement=residuals.measurement,
        extensions=residuals.extensions,
        dual_solution=_localdisc_dual_solution(problem, optimization),
        optimization_result=optimization,
        residuals,
        warnings,
        message=if agreed
            "the configured outer-hierarchy primal and dual agree within residual tolerance; this is an upper relaxation of separable discrimination"
        else
            "a residual-checked outer-hierarchy primal is available, but no matching dual bound was established"
        end,
    )
end

const _SEPARABILITY_DETERMINISTIC_STRATEGIES = (
    :ppt,
    :low_rank_ppt,
    :realignment,
    :centered_realignment,
    :reduction,
    :qubit_qudit,
    :rank4_chow,
    :separable_ball,
    :rank_one_identity,
    :operator_schmidt_rank,
    :positive_maps,
    :filter_covariance,
    :symmetric_extension,
    :symmetric_inner_extension,
)

const _SEPARABILITY_ALL_STRATEGIES = (
    _SEPARABILITY_DETERMINISTIC_STRATEGIES..., :randomized_subtraction
)

function _separability_strategies(strategies)
    values = if strategies === :qetlab_deterministic
        _SEPARABILITY_DETERMINISTIC_STRATEGIES
    elseif strategies === :full
        _SEPARABILITY_ALL_STRATEGIES
    elseif strategies isa Symbol
        (strategies,)
    elseif strategies isa Tuple || strategies isa AbstractVector
        strategies isa AbstractVector && Base.require_one_based_indexing(strategies)
        Tuple(strategies)
    else
        throw(
            ArgumentError(
                "strategies must be :qetlab_deterministic, :full, a Symbol, or a collection of Symbols",
            ),
        )
    end
    isempty(values) &&
        throw(ArgumentError("at least one separability strategy is required"))
    all(value -> value isa Symbol, values) ||
        throw(ArgumentError("every separability strategy must be a Symbol"))
    length(unique(values)) == length(values) ||
        throw(ArgumentError("separability strategies must not contain duplicates"))
    unsupported = setdiff(values, _SEPARABILITY_ALL_STRATEGIES)
    isempty(unsupported) ||
        throw(ArgumentError("unsupported separability strategies: $(Tuple(unsupported))"))
    return values
end

function _separability_prepare(
    rho::AbstractMatrix{<:Number}, dims; atol, rtol, allow_densify::Bool, max_dense_entries
)
    Base.require_one_based_indexing(rho)
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("is_separable requires a square matrix"))
    ishermitian(rho) || throw(
        ArgumentError(
            "is_separable requires an exactly Hermitian state; no Hermitian part is substituted",
        ),
    )
    _sepopt_guard_dense(size(rho, 1), max_dense_entries, "separability analysis")
    analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="is_separable",
    )
    represented_trace = _sepopt_represented_trace(rho)
    iszero(imag(represented_trace)) &&
    real(represented_trace) == one(real(represented_trace)) || throw(
        ArgumentError(
            "is_separable requires an exactly represented unit-trace state; " *
            "a trace accepted only by atol/rtol is not eligible for a certificate",
        ),
    )
    layout = _tierd_bipartite_layout(dims, size(rho, 1))
    matrix = copy(rho)
    marginal_a = partial_trace(matrix, layout; trace_out=(2,))
    marginal_b = partial_trace(matrix, layout; trace_out=(1,))
    return (
        matrix,
        layout,
        dimensions=(layout[1], layout[2]),
        analysis,
        marginal_a=Matrix(marginal_a),
        marginal_b=Matrix(marginal_b),
        tolerance=analysis.tolerance,
    )
end

function _separability_attempt(
    method::Symbol,
    status::Symbol,
    certified::Bool,
    certificate_kind,
    raw_result,
    message::AbstractString;
    backend::Symbol=:native,
)
    return _validated_entanglement_attempt(
        method, backend, status, certified, certificate_kind, raw_result, message
    )
end

function _separability_report(attempts::Vector{EntanglementAttempt})
    all(_is_package_validated_entanglement_attempt, attempts) ||
        throw(ArgumentError("separability reports require package-validated attempts"))
    decisive = findlast(
        attempt -> attempt.certified && _is_package_validated_entanglement_attempt(attempt),
        attempts,
    )
    if decisive !== nothing
        attempt = attempts[decisive]
        return _validated_entanglement_report(
            attempt.status,
            true,
            attempt.certificate_kind,
            attempt.method,
            attempt.backend,
            attempt.raw_result,
            attempts,
            attempt.message,
        )
    end
    return _validated_entanglement_report(
        :unknown,
        false,
        nothing,
        :composite_separability_analysis,
        :native,
        nothing,
        attempts,
        "the configured exact, sufficient, necessary, heuristic, and hierarchy attempts produced no validated separability or entanglement certificate",
    )
end

function _separability_push!(
    attempts::Vector{EntanglementAttempt}, attempt::EntanglementAttempt
)
    push!(attempts, attempt)
    return attempt.certified
end

function _separability_rank(matrix::AbstractMatrix, tolerance)
    values = eigvals(Hermitian(Matrix(matrix)))
    scale = maximum(abs, values; init=one(tolerance))
    threshold = max(tolerance, 64 * eps(typeof(tolerance)) * max(one(scale), scale))
    rank = count(value -> value > threshold, values)
    boundary = any(value -> !iszero(value) && abs(value) <= threshold, values)
    exact_zeros = count(iszero, values)
    certain = !boundary && rank + exact_zeros == length(values)
    return (; rank, certain, threshold, eigenvalues=values)
end

function _separability_centered_realignment(prepared)
    centered = prepared.matrix - tensor_product(prepared.marginal_a, prepared.marginal_b)
    realigned = realign(centered, prepared.layout; systems=(1,))
    value = trace_norm(realigned; allow_densify=false)
    purity_a = real(tr(prepared.marginal_a * prepared.marginal_a))
    purity_b = real(tr(prepared.marginal_b * prepared.marginal_b))
    radicand = (one(purity_a) - purity_a) * (one(purity_b) - purity_b)
    tolerance = prepared.tolerance
    if radicand < -tolerance
        return (
            status=:unknown,
            value,
            bound=nothing,
            realigned,
            message="the centered-realignment bound has a negative radicand outside numerical tolerance",
        )
    end
    bound = sqrt(max(zero(radicand), radicand))
    margin = value - bound
    status = margin > tolerance ? :entangled : :unknown
    return (
        status,
        value,
        bound,
        margin,
        realigned,
        message=if status === :entangled
            "the centered realignment/covariance inequality is violated with margin"
        else
            "the centered realignment/covariance inequality gives no certificate"
        end,
    )
end

function _separability_qubit_qudit(prepared, ppt_satisfied::Bool)
    dimensions = prepared.dimensions
    minimum(dimensions) == 2 || return nothing
    ppt_satisfied || return (
        status=:unknown,
        branch=:ppt_prerequisite,
        message="qubit--qudit sufficient tests require a robust PPT prerequisite",
    )
    eigenvalues = sort(real.(prepared.analysis.eigenvalues); rev=true)
    n = maximum(dimensions)
    left = (eigenvalues[1] - eigenvalues[2n - 1])^2
    right = 4 * eigenvalues[2n - 2] * eigenvalues[2n]
    margin = right - left
    if margin > prepared.tolerance
        return (
            status=:separable,
            branch=:spectrum,
            left,
            right,
            margin,
            message="the qubit--qudit separability-from-spectrum inequality holds with margin",
        )
    end

    ordered = if dimensions[1] == 2
        prepared.matrix
    else
        swap_subsystems(prepared.matrix, prepared.layout; systems=(1, 2))
    end
    A = Matrix(@view ordered[1:n, 1:n])
    B = Matrix(@view ordered[1:n, (n + 1):(2n)])
    C = Matrix(@view ordered[(n + 1):(2n), (n + 1):(2n)])
    R = typeof(prepared.tolerance)
    skew = B - adjoint(B)
    if all(iszero, skew)
        return (
            status=:separable,
            branch=:block_hankel,
            residual=zero(prepared.tolerance),
            message="the off-diagonal block is exactly Hermitian, satisfying the block-Hankel theorem",
        )
    end
    skew_singular = svdvals(skew)
    if length(skew_singular) <= 1 || all(iszero, skew_singular[2:end])
        return (
            status=:separable,
            branch=:perturbed_block_hankel,
            singular_values=skew_singular,
            message="the off-diagonal skew-Hermitian block has exact represented rank at most one",
        )
    end

    five_sixths = R(5) / R(6)
    one_sixth = inv(R(6))
    homothetic = [
        five_sixths * A-one_sixth * C B
        adjoint(B) five_sixths * C-one_sixth * A
    ]
    homothetic_psd = is_positive_semidefinite(
        homothetic; atol=0, rtol=prepared.analysis.rtol, allow_densify=true
    )
    homothetic_ppt = ppt_criterion(
        homothetic,
        (2, n);
        systems=(2,),
        atol=0,
        rtol=prepared.analysis.rtol,
        allow_densify=true,
    )
    if homothetic_psd.status === MatrixPredicateSatisfied &&
        homothetic_ppt.status === CriterionSatisfied
        return (
            status=:separable,
            branch=:homothetic_image,
            homothetic,
            homothetic_psd,
            homothetic_ppt,
            message="the Hildebrand homothetic image is robustly PSD and PPT",
        )
    end

    min_a = minimum(real, eigvals(Hermitian(A)))
    min_c = minimum(real, eigvals(Hermitian(C)))
    block_left = opnorm(B)^2
    block_right = min_a * min_c
    block_margin = block_right - block_left
    if block_margin > prepared.tolerance
        return (
            status=:separable,
            branch=:block_norm,
            left=block_left,
            right=block_right,
            margin=block_margin,
            message="the qubit--qudit block-norm sufficient inequality holds with margin",
        )
    end
    return (
        status=:unknown,
        branch=:all_qubit_qudit_tests,
        spectrum_margin=margin,
        skew_singular_values=skew_singular,
        homothetic_psd,
        homothetic_ppt,
        block_margin,
        message="the qubit--qudit sufficient tests produced no robust certificate",
    )
end

function _separability_plucker(range_basis, indices::NTuple{4,Int})
    return det(range_basis[collect(indices), :])
end

function _separability_rank4_chow(prepared, state_rank)
    prepared.dimensions == (3, 3) || return nothing
    state_rank.certain && state_rank.rank == 4 || return (
        status=:unknown,
        determinant=nothing,
        message="the 3×3 Chow-form branch requires a numerically decisive rank-four support",
    )
    range_basis = prepared.analysis.decomposition.vectors[:, end:-1:(end - 3)]
    p(i, j, k, l) = _separability_plucker(range_basis, (i, j, k, l))
    matrix = [
        p(1, 2, 4, 5) p(1, 3, 4, 6) p(2, 3, 5, 6) p(1, 2, 4, 6)+p(1, 3, 4, 5) p(1, 2, 5, 6)+p(
            2, 3, 4, 5
        ) p(1, 3, 5, 6)+p(2, 3, 4, 6)
        p(1, 2, 7, 8) p(1, 3, 7, 9) p(2, 3, 8, 9) p(1, 2, 7, 9)+p(1, 3, 7, 8) p(1, 2, 8, 9)+p(
            2, 3, 7, 8
        ) p(1, 3, 8, 9)+p(2, 3, 7, 9)
        p(4, 5, 7, 8) p(4, 6, 7, 9) p(5, 6, 8, 9) p(4, 5, 7, 9)+p(4, 6, 7, 8) p(4, 5, 8, 9)+p(
            5, 6, 7, 8
        ) p(4, 6, 8, 9)+p(5, 6, 7, 9)
        p(1, 2, 4, 8)-p(1, 2, 5, 7) p(1, 3, 4, 9)-p(1, 3, 6, 7) p(2, 3, 5, 9)-p(2, 3, 6, 8) p(1, 2, 4, 9)-p(1, 2, 6, 7)+p(1, 3, 4, 8)-p(
            1, 3, 5, 7
        ) p(1, 2, 5, 9)-p(1, 2, 6, 8)+p(2, 3, 4, 8)-p(2, 3, 5, 7) p(1, 3, 5, 9)-p(1, 3, 6, 8)+p(2, 3, 4, 9)-p(
            2, 3, 6, 7
        )
        p(1, 4, 5, 8)-p(2, 4, 5, 7) p(1, 4, 6, 9)-p(3, 4, 6, 7) p(2, 5, 6, 9)-p(3, 5, 6, 8) p(1, 4, 5, 9)-p(2, 4, 6, 7)+p(1, 4, 6, 8)-p(
            3, 4, 5, 7
        ) p(1, 5, 6, 8)-p(2, 5, 6, 7)+p(2, 4, 5, 9)-p(3, 4, 5, 8) p(1, 5, 6, 9)-p(3, 4, 6, 8)+p(2, 4, 6, 9)-p(
            3, 5, 6, 7
        )
        p(1, 5, 7, 8)-p(2, 4, 7, 8) p(1, 6, 7, 9)-p(3, 4, 7, 9) p(2, 6, 8, 9)-p(3, 5, 8, 9) p(1, 5, 7, 9)-p(2, 4, 7, 9)+p(1, 6, 7, 8)-p(3, 4, 7, 8) p(1, 5, 8, 9)-p(2, 4, 8, 9)+p(2, 6, 7, 8)-p(3, 5, 7, 8) p(1, 6, 8, 9)-p(3, 4, 6, 8)+p(2, 6, 7, 9)-p(3, 5, 7, 9)
    ]
    determinant = det(matrix)
    scale = max(one(abs(determinant)), norm(matrix)^6)
    tolerance = prepared.analysis.atol + prepared.analysis.rtol * scale
    if abs(determinant) > tolerance
        return (
            status=:entangled,
            determinant,
            tolerance,
            chow_matrix=matrix,
            message="the 3×3 rank-four Chow determinant is nonzero with margin",
        )
    elseif iszero(determinant)
        return (
            status=:separable,
            determinant,
            tolerance,
            chow_matrix=matrix,
            message="the represented 3×3 rank-four Chow determinant is exactly zero",
        )
    end
    return (
        status=:unknown,
        determinant,
        tolerance,
        chow_matrix=matrix,
        message="the 3×3 rank-four Chow determinant lies in the numerical boundary",
    )
end

function _separability_positive_maps(prepared)
    results = Any[]
    R = typeof(real(zero(eltype(prepared.matrix))))
    if prepared.dimensions == (3, 3)
        direct_parameters = [R(index) / R(10) for index in 0:9]
        parameters = vcat(
            direct_parameters, [inv(t) for t in Iterators.drop(direct_parameters, 1)]
        )
        for t in parameters
            denominator = one(t) - t + t^2
            a = (one(t) - t)^2 / denominator
            b = t^2 / denominator
            c = one(t) / denominator
            map = choi_map(a, b, c)
            mapped = partial_map(prepared.matrix, map, 2, prepared.layout)
            predicate = is_positive_semidefinite(
                mapped;
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
            )
            evidence = (; family=:generalized_choi, t, a, b, c, mapped, predicate)
            push!(results, evidence)
            if predicate.status === MatrixPredicateViolated
                return (
                    status=:entangled,
                    evidence,
                    attempts=Tuple(results),
                    message="a generalized Choi positive map produced a robust PSD violation",
                )
            end
        end
    end

    for subsystem in 1:2
        dimension = prepared.dimensions[subsystem]
        iseven(dimension) || continue
        signs = vcat(ones(R, dimension ÷ 2), -ones(R, dimension ÷ 2))
        antisymmetric = reverse(Matrix(Diagonal(signs)); dims=2)
        maximally_entangled_vector = vec(Matrix{R}(I, dimension, dimension))
        identity_local = Matrix{R}(I, dimension, dimension)
        flip = swap_operator(
            (dimension, dimension); systems=(1, 2), T=R, sparse_output=false
        )
        lifted_antisymmetric = kron(identity_local, antisymmetric)
        choi =
            Matrix{R}(I, dimension^2, dimension^2) -
            maximally_entangled_vector * adjoint(maximally_entangled_vector) -
            lifted_antisymmetric * flip * adjoint(lifted_antisymmetric)
        map = ChoiRepresentation(choi, dimension, dimension)
        mapped = partial_map(prepared.matrix, map, subsystem, prepared.layout)
        predicate = is_positive_semidefinite(
            mapped;
            atol=prepared.analysis.atol,
            rtol=prepared.analysis.rtol,
            allow_densify=true,
        )
        evidence = (; family=:breuer_hall, subsystem, mapped, predicate)
        push!(results, evidence)
        if predicate.status === MatrixPredicateViolated
            return (
                status=:entangled,
                evidence,
                attempts=Tuple(results),
                message="a Breuer--Hall positive map produced a robust PSD violation",
            )
        end
    end
    return (
        status=:unknown,
        attempts=Tuple(results),
        message="the configured positive-map family produced no robust PSD violation",
    )
end

function _separability_filter_covariance(
    prepared; max_iterations, max_condition_number, max_entries, max_work
)
    result = filter_normal_form(
        prepared.matrix,
        prepared.layout;
        max_iterations=max_iterations,
        max_condition_number=max_condition_number,
        allow_densify=true,
        max_entries=max_entries,
        max_work=max_work,
    )
    result.status === :converged || return (
        status=:unknown,
        result,
        value=nothing,
        bound=nothing,
        message="the bounded filter-normal-form computation did not converge",
    )
    R = typeof(prepared.tolerance)
    n = minimum(prepared.dimensions)
    x = maximum(prepared.dimensions)
    n_real = R(n)
    x_real = R(x)
    value = sum(result.coefficients)
    first_bound = sqrt(n_real * x_real * (n_real - one(R)) * (x_real - one(R)))
    second_bound =
        n_real *
        x_real *
        (
            one(R) - inv(n_real) +
            (n_real^2 - one(R)) / x_real +
            min(zero(R), (x_real^2 - n_real^2) / x_real - (x_real - one(R)))
        ) / R(2)
    bound = min(first_bound, second_bound)
    margin = value - bound
    status = margin > prepared.tolerance ? :entangled : :unknown
    return (
        status,
        result,
        value,
        bound,
        margin,
        message=if status === :entangled
            "the filter covariance-matrix bound is violated with margin"
        else
            "the filter covariance-matrix test gives no certificate"
        end,
    )
end

function _separability_random_vector(
    rng::AbstractRNG, ::Type{T}, dimension::Int
) where {T<:AbstractFloat}
    for _ in 1:8
        vector = randn(rng, T, dimension) + im * randn(rng, T, dimension)
        norm_value = norm(vector)
        iszero(norm_value) || return vector / norm_value
    end
    return throw(ErrorException("the explicit RNG produced eight consecutive zero vectors"))
end

function _separability_contraction_a(rho, b, dimensions)
    dimension_a, dimension_b = dimensions
    T = eltype(rho)
    result = zeros(T, dimension_a, dimension_a)
    for column_a in 1:dimension_a, row_a in 1:dimension_a
        value = zero(T)
        for column_b in 1:dimension_b, row_b in 1:dimension_b
            row = (row_a - 1) * dimension_b + row_b
            column = (column_a - 1) * dimension_b + column_b
            value += conj(b[row_b]) * rho[row, column] * b[column_b]
        end
        result[row_a, column_a] = value
    end
    return result
end

function _separability_contraction_b(rho, a, dimensions)
    dimension_a, dimension_b = dimensions
    T = eltype(rho)
    result = zeros(T, dimension_b, dimension_b)
    for column_b in 1:dimension_b, row_b in 1:dimension_b
        value = zero(T)
        for column_a in 1:dimension_a, row_a in 1:dimension_a
            row = (row_a - 1) * dimension_b + row_b
            column = (column_a - 1) * dimension_b + column_b
            value += conj(a[row_a]) * rho[row, column] * a[column_a]
        end
        result[row_b, column_b] = value
    end
    return result
end

function _separability_best_product(
    rng::AbstractRNG,
    rho,
    dimensions;
    restarts::Int,
    iterations::Int,
    budget::_SeparabilityWorkBudget,
)
    R = typeof(real(zero(eltype(rho))))
    best_overlap = -convert(R, Inf)
    best_factors = nothing
    draws = 0
    completed = 0
    per_iteration =
        BigInt(dimensions[1])^2 * BigInt(dimensions[2])^2 +
        BigInt(dimensions[1])^3 +
        BigInt(dimensions[2])^3
    stopped = false
    for _ in 1:restarts
        if !_sepopt_consume!(budget, 2 * sum(dimensions), "random product initialization")
            stopped = true
            break
        end
        a = _separability_random_vector(rng, R, dimensions[1])
        b = _separability_random_vector(rng, R, dimensions[2])
        draws += 2 * sum(dimensions)
        for _ in 1:iterations
            if !_sepopt_consume!(budget, per_iteration, "alternating product iteration")
                stopped = true
                break
            end
            effective_a = _separability_contraction_a(rho, b, dimensions)
            decomposition_a = eigen(Hermitian(effective_a))
            a = decomposition_a.vectors[:, argmax(decomposition_a.values)]
            effective_b = _separability_contraction_b(rho, a, dimensions)
            decomposition_b = eigen(Hermitian(effective_b))
            b = decomposition_b.vectors[:, argmax(decomposition_b.values)]
            completed += 1
        end
        product = kron(a, b)
        overlap = real(dot(product, rho, product))
        if overlap > best_overlap
            best_overlap = overlap
            best_factors = (copy(a), copy(b))
        end
        stopped && break
    end
    return (
        overlap=best_overlap,
        factors=best_factors,
        draws,
        iterations=completed,
        work_used=budget.used,
        stopped,
    )
end

function _separability_subtraction_step(state, product, tolerance)
    decomposition = eigen(Hermitian(state))
    coordinates = adjoint(decomposition.vectors) * product
    support = findall(value -> value > tolerance, decomposition.values)
    null_support = setdiff(collect(eachindex(decomposition.values)), support)
    null_residual = norm(coordinates[null_support])
    null_residual <= tolerance || return (
        accepted=false, reason=:outside_range, null_residual, maximum_weight=nothing
    )
    isempty(support) &&
        return (accepted=false, reason=:zero_support, null_residual, maximum_weight=nothing)
    denominator = sum(
        abs2(coordinates[index]) / decomposition.values[index] for index in support
    )
    isfinite(denominator) && denominator > zero(denominator) || return (
        accepted=false,
        reason=:invalid_generalized_eigenvalue,
        null_residual,
        maximum_weight=nothing,
    )
    maximum_weight = inv(denominator)
    weight = maximum_weight / 2
    weight > tolerance && weight < one(weight) - tolerance ||
        return (accepted=false, reason=:weight_boundary, null_residual, maximum_weight)
    candidate = state - weight * (product * adjoint(product))
    minimum_value = minimum(real, eigvals(Hermitian(candidate)))
    minimum_value >= zero(minimum_value) || return (
        accepted=false,
        reason=:psd_boundary,
        null_residual,
        maximum_weight,
        weight,
        minimum_eigenvalue=minimum_value,
    )
    next_state = candidate / (one(weight) - weight)
    return (
        accepted=true,
        reason=:validated_subtraction,
        null_residual,
        maximum_weight,
        weight,
        minimum_eigenvalue=minimum_value,
        next_state,
    )
end

function _separability_randomized_subtraction(
    rng::AbstractRNG,
    prepared;
    max_subtractions,
    max_product_restarts,
    max_product_iterations,
    budget::_SeparabilityWorkBudget,
)
    state = copy(prepared.matrix)
    remaining_weight = one(prepared.tolerance)
    components = NamedTuple[]
    history = NamedTuple[]
    total_dimension = prod(BigInt(dimension) for dimension in prepared.dimensions)
    cubic_work = total_dimension^3
    for subtraction in 1:max_subtractions
        candidate = _separability_best_product(
            rng,
            state,
            prepared.dimensions;
            restarts=max_product_restarts,
            iterations=max_product_iterations,
            budget=budget,
        )
        candidate.factors === nothing && return (
            status=:unknown,
            final_state=state,
            components=Tuple(components),
            history=Tuple(history),
            work_used=budget.used,
            message="the explicit work budget stopped before a product candidate was produced",
        )
        product = kron(candidate.factors...)
        if !_sepopt_consume!(budget, 2cubic_work, "product-subtraction spectral validation")
            return (
                status=:unknown,
                final_state=state,
                components=Tuple(components),
                history=Tuple(history),
                work_used=budget.used,
                message="the explicit work budget stopped before validating the proposed product subtraction",
            )
        end
        step = _separability_subtraction_step(state, product, prepared.tolerance)
        push!(history, (; subtraction, candidate, step))
        step.accepted || continue
        push!(
            components,
            (
                weight=remaining_weight * step.weight,
                local_factors=candidate.factors,
                product=product,
            ),
        )
        remaining_weight *= one(step.weight) - step.weight
        state = step.next_state
        if !_sepopt_consume!(budget, cubic_work, "post-subtraction separable-ball test")
            return (
                status=:unknown,
                final_state=state,
                remaining_weight,
                components=Tuple(components),
                history=Tuple(history),
                work_used=budget.used,
                message="the explicit work budget stopped before testing the residual state against the separable ball",
            )
        end
        ball = in_separable_ball(
            state,
            prepared.layout;
            atol=prepared.analysis.atol,
            rtol=prepared.analysis.rtol,
            allow_densify=true,
        )
        if ball.status === :separable_certified
            reconstruction = remaining_weight * state
            for component in components
                reconstruction +=
                    component.weight * (component.product * adjoint(component.product))
            end
            residual = norm(reconstruction - prepared.matrix)
            tolerance =
                prepared.analysis.atol +
                prepared.analysis.rtol * max(one(residual), norm(prepared.matrix))
            if residual <= tolerance
                return (
                    status=:separable,
                    final_state=state,
                    remaining_weight,
                    ball,
                    components=Tuple(components),
                    reconstruction,
                    reconstruction_residual=residual,
                    reconstruction_tolerance=tolerance,
                    history=Tuple(history),
                    work_used=budget.used,
                    message="validated product-state subtractions reduce the residual state to the Gurvits--Barnum separable ball",
                )
            end
        end
    end
    return (
        status=:unknown,
        final_state=state,
        remaining_weight,
        components=Tuple(components),
        history=Tuple(history),
        work_used=budget.used,
        message="the bounded randomized subtraction heuristic produced no separability certificate",
    )
end

function _separability_diagonal_decomposition(prepared)
    diagonal = diag(prepared.matrix)
    off_diagonal = prepared.matrix - Diagonal(diagonal)
    all(iszero, off_diagonal) || return nothing
    all(value -> iszero(imag(value)) && real(value) >= zero(real(value)), diagonal) ||
        return nothing
    components = Tuple(
        (
            weight=real(diagonal[index]),
            basis_index=index,
            local_indices=linear_to_basis(index, prepared.layout),
        ) for index in eachindex(diagonal) if !iszero(diagonal[index])
    )
    return (
        components,
        reconstruction=Diagonal(copy(diagonal)),
        residual=zero(prepared.tolerance),
        message="the state is an exact represented convex combination of computational-basis product projectors",
    )
end

function _separability_budget_attempt(strategy::Symbol, budget)
    return _separability_attempt(
        strategy,
        :unknown,
        false,
        nothing,
        (work_used=budget.used, max_work=budget.limit),
        "the deterministic max_work budget stopped this strategy before its large allocation",
    )
end

function _is_separable(
    rng,
    rho::AbstractMatrix{<:Number},
    dims;
    strategies=:qetlab_deterministic,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    extension_orders=(2,),
    extension_ppt::Bool=true,
    extension_bosonic::Bool=true,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_dense_entries=1_000_000,
    max_work=1_000_000_000,
    max_extension_order=6,
    max_filter_iterations=10_000,
    max_filter_condition_number=nothing,
    max_filter_entries=10_000_000,
    max_filter_work=1_000_000_000,
    max_subtractions=16,
    max_product_restarts=8,
    max_product_iterations=32,
    limits::OptimizationLimits=OptimizationLimits(),
)
    selected = _separability_strategies(strategies)
    :randomized_subtraction in selected &&
        rng === nothing &&
        throw(
            ArgumentError(
                "strategy=:randomized_subtraction requires an explicit rng::AbstractRNG as the first argument",
            ),
        )
    rng === nothing ||
        rng isa AbstractRNG ||
        throw(ArgumentError("rng must be an AbstractRNG"))
    orders = if extension_orders isa Integer
        (extension_orders,)
    elseif extension_orders isa Tuple || extension_orders isa AbstractVector
        extension_orders isa AbstractVector &&
            Base.require_one_based_indexing(extension_orders)
        Tuple(extension_orders)
    else
        throw(ArgumentError("extension_orders must be an integer or collection"))
    end
    all(order -> order isa Integer && !(order isa Bool) && order >= 2, orders) ||
        throw(ArgumentError("every extension order must be an integer at least two"))
    length(unique(orders)) == length(orders) ||
        throw(ArgumentError("extension_orders must not contain duplicates"))
    checked_max_extension = _sepopt_positive_int(max_extension_order, "max_extension_order")
    checked_filter_iterations = _sepopt_nonnegative_int(
        max_filter_iterations, "max_filter_iterations"
    )
    checked_subtractions = _sepopt_nonnegative_int(max_subtractions, "max_subtractions")
    checked_restarts = _sepopt_positive_int(max_product_restarts, "max_product_restarts")
    checked_product_iterations = _sepopt_nonnegative_int(
        max_product_iterations, "max_product_iterations"
    )
    prepared = _separability_prepare(
        rho,
        dims;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
    )
    budget = _SeparabilityWorkBudget(BigInt(0), _sepopt_limit(max_work, "max_work"))
    attempts = EntanglementAttempt[]
    total_dimension = prod(prepared.dimensions)
    cubic_work = BigInt(total_dimension)^3

    if minimum(prepared.dimensions) == 1 && !prepared.analysis.structural_boundary_uncertain
        attempt = _separability_attempt(
            :trivial_local_dimension,
            :separable,
            true,
            :one_dimensional_local_factor_theorem,
            (dimensions=prepared.dimensions,),
            "every positive bipartite operator with a one-dimensional local factor is separable",
        )
        push!(attempts, attempt)
        return _separability_report(attempts)
    end

    diagonal = _separability_diagonal_decomposition(prepared)
    if diagonal !== nothing
        attempt = _separability_attempt(
            :product_decomposition,
            :separable,
            true,
            :computational_basis_product_decomposition,
            diagonal,
            diagonal.message,
        )
        push!(attempts, attempt)
        return _separability_report(attempts)
    end

    needs_ppt = any(
        strategy in (
            :ppt,
            :low_rank_ppt,
            :qubit_qudit,
            :rank4_chow,
            :rank_one_identity,
            :operator_schmidt_rank,
        ) for strategy in selected
    )
    ppt_result = nothing
    ppt_satisfied = false
    if needs_ppt
        if _sepopt_consume!(budget, cubic_work, "partial-transpose eigendecomposition")
            ppt_result = ppt_criterion(
                prepared.matrix,
                prepared.layout;
                systems=(2,),
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
            )
            if ppt_result.status === CriterionEntanglementDetected
                attempt = _separability_attempt(
                    :ppt,
                    :entangled,
                    true,
                    :negative_partial_transpose_witness,
                    ppt_result,
                    "a robust negative partial-transpose eigenpair certifies entanglement",
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            ppt_satisfied = ppt_result.status === CriterionSatisfied
            if :ppt in selected
                if ppt_satisfied && (
                    prepared.dimensions in ((2, 2), (2, 3)) ||
                    reverse(prepared.dimensions) in ((2, 2), (2, 3))
                )
                    attempt = _separability_attempt(
                        :ppt,
                        :separable,
                        true,
                        :ppt_low_dimension_theorem,
                        ppt_result,
                        "PPT is sufficient for separability in 2×2 and 2×3 systems",
                    )
                    push!(attempts, attempt)
                    return _separability_report(attempts)
                end
                push!(
                    attempts,
                    _separability_attempt(
                        :ppt,
                        :unknown,
                        false,
                        nothing,
                        ppt_result,
                        if ppt_result.status === CriterionSatisfied
                            "PPT is only a necessary condition in these dimensions"
                        else
                            "the PPT test lies on its numerical boundary"
                        end,
                    ),
                )
            end
        else
            push!(attempts, _separability_budget_attempt(:ppt, budget))
        end
    end

    state_rank = _separability_rank(prepared.matrix, prepared.tolerance)
    if :low_rank_ppt in selected
        if ppt_satisfied && _sepopt_consume!(budget, 2cubic_work, "marginal rank tests")
            rank_a = _separability_rank(prepared.marginal_a, prepared.tolerance)
            rank_b = _separability_rank(prepared.marginal_b, prepared.tolerance)
            condition =
                state_rank.certain &&
                rank_a.certain &&
                rank_b.certain &&
                (
                    state_rank.rank <= 3 ||
                    state_rank.rank <= rank_a.rank ||
                    state_rank.rank <= rank_b.rank
                )
            evidence = (; state_rank, rank_a, rank_b, ppt=ppt_result)
            if condition
                attempt = _separability_attempt(
                    :low_rank_ppt,
                    :separable,
                    true,
                    :low_rank_ppt_theorem,
                    evidence,
                    "the robust PPT and low-rank theorem certifies separability on the reduced local supports",
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :low_rank_ppt,
                    :unknown,
                    false,
                    nothing,
                    evidence,
                    "the low-rank PPT theorem did not yield a decisive certificate",
                ),
            )
        elseif !_sepopt_consume!(budget, 0, "low-rank PPT prerequisite")
            push!(attempts, _separability_budget_attempt(:low_rank_ppt, budget))
        else
            push!(
                attempts,
                _separability_attempt(
                    :low_rank_ppt,
                    :unknown,
                    false,
                    nothing,
                    (state_rank, ppt=ppt_result),
                    "the low-rank sufficient theorem lacks a robust PPT or rank prerequisite",
                ),
            )
        end
    end

    if :realignment in selected
        if _sepopt_consume!(budget, cubic_work, "realignment singular-value decomposition")
            result = realignment_criterion(
                prepared.matrix,
                prepared.layout;
                systems=(1,),
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
            )
            if result.status === CriterionEntanglementDetected
                attempt = _separability_attempt(
                    :realignment,
                    :entangled,
                    true,
                    :realignment_cross_norm_violation,
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :realignment,
                    :unknown,
                    false,
                    nothing,
                    result,
                    "the realignment necessary criterion gives no separability certificate",
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:realignment, budget))
        end
    end

    if :centered_realignment in selected
        if _sepopt_consume!(
            budget, cubic_work, "centered-realignment singular-value decomposition"
        )
            result = _separability_centered_realignment(prepared)
            if result.status === :entangled
                attempt = _separability_attempt(
                    :centered_realignment,
                    :entangled,
                    true,
                    :centered_realignment_covariance_violation,
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :centered_realignment, :unknown, false, nothing, result, result.message
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:centered_realignment, budget))
        end
    end

    if :reduction in selected
        if _sepopt_consume!(budget, 2cubic_work, "reduction-map eigendecompositions")
            result = reduction_criterion(
                prepared.matrix,
                prepared.layout;
                side=:both,
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
            )
            if result.status === CriterionEntanglementDetected
                attempt = _separability_attempt(
                    :reduction,
                    :entangled,
                    true,
                    :reduction_map_witness,
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :reduction,
                    :unknown,
                    false,
                    nothing,
                    result,
                    "the reduction necessary criterion gives no separability certificate",
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:reduction, budget))
        end
    end

    if :qubit_qudit in selected && minimum(prepared.dimensions) == 2
        if _sepopt_consume!(budget, 4cubic_work, "qubit--qudit sufficient tests")
            result = _separability_qubit_qudit(prepared, ppt_satisfied)
            if result !== nothing && result.status === :separable
                attempt = _separability_attempt(
                    :qubit_qudit,
                    :separable,
                    true,
                    Symbol(:qubit_qudit_, result.branch),
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            result === nothing || push!(
                attempts,
                _separability_attempt(
                    :qubit_qudit, :unknown, false, nothing, result, result.message
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:qubit_qudit, budget))
        end
    end

    if :rank4_chow in selected && prepared.dimensions == (3, 3)
        if !ppt_satisfied
            push!(
                attempts,
                _separability_attempt(
                    :rank4_chow,
                    :unknown,
                    false,
                    nothing,
                    (state_rank, ppt=ppt_result),
                    "the rank-four Chow-form theorem lacks a robust PPT prerequisite",
                ),
            )
        elseif _sepopt_consume!(budget, 10cubic_work, "rank-four Chow determinant")
            result = _separability_rank4_chow(prepared, state_rank)
            status = result.status
            if status in (:separable, :entangled)
                attempt = _separability_attempt(
                    :rank4_chow, status, true, :rank4_chow_form, result, result.message
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :rank4_chow, :unknown, false, nothing, result, result.message
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:rank4_chow, budget))
        end
    end

    if :separable_ball in selected
        if _sepopt_consume!(budget, cubic_work, "separable-ball eigendecomposition")
            result = in_separable_ball(
                prepared.matrix,
                prepared.layout;
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
            )
            if result.status === :separable_certified
                attempt = _separability_attempt(
                    :separable_ball,
                    :separable,
                    true,
                    :gurvits_barnum_separable_ball,
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :separable_ball, :unknown, false, nothing, result, result.message
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:separable_ball, budget))
        end
    end

    if :rank_one_identity in selected
        eigenvalues = sort(real.(prepared.analysis.eigenvalues); rev=true)
        floating_tail_equal =
            length(eigenvalues) <= 2 || all(==(eigenvalues[2]), eigenvalues[3:end])
        push!(
            attempts,
            _separability_attempt(
                :rank_one_identity,
                :unknown,
                false,
                nothing,
                (
                    eigenvalues,
                    ppt=ppt_result,
                    floating_tail_equal,
                    input_level_exact_proof=false,
                ),
                if ppt_satisfied && floating_tail_equal
                    "floating eigensolver equality is only numerical evidence; no input-level exact rank-one identity-perturbation proof is available"
                else
                    "the rank-one identity-perturbation theorem gives no certificate"
                end,
            ),
        )
    end

    if :operator_schmidt_rank in selected
        if _sepopt_consume!(budget, cubic_work, "operator-Schmidt decomposition")
            decomposition = operator_schmidt_decomposition(
                prepared.matrix, prepared.layout; allow_densify=true
            )
            exact_rank_at_most_two =
                length(decomposition.coefficients) <= 2 ||
                all(iszero, decomposition.coefficients[3:end])
            if ppt_satisfied && exact_rank_at_most_two
                attempt = _separability_attempt(
                    :operator_schmidt_rank,
                    :separable,
                    true,
                    :ppt_operator_schmidt_rank_two_theorem,
                    decomposition,
                    "the PPT state's represented operator-Schmidt rank is at most two",
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :operator_schmidt_rank,
                    :unknown,
                    false,
                    nothing,
                    decomposition,
                    "the operator-Schmidt-rank sufficient theorem gives no certificate",
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:operator_schmidt_rank, budget))
        end
    end

    if :positive_maps in selected
        map_count =
            (prepared.dimensions == (3, 3) ? 19 : 0) + count(iseven, prepared.dimensions)
        if _sepopt_consume!(
            budget, map_count * cubic_work, "positive-map eigendecompositions"
        )
            result = _separability_positive_maps(prepared)
            if result.status === :entangled
                attempt = _separability_attempt(
                    :positive_maps,
                    :entangled,
                    true,
                    :positive_map_psd_violation,
                    result,
                    result.message,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :positive_maps, :unknown, false, nothing, result, result.message
                ),
            )
        else
            push!(attempts, _separability_budget_attempt(:positive_maps, budget))
        end
    end

    if :filter_covariance in selected
        result = _separability_filter_covariance(
            prepared;
            max_iterations=checked_filter_iterations,
            max_condition_number=max_filter_condition_number,
            max_entries=max_filter_entries,
            max_work=max_filter_work,
        )
        if result.status === :entangled
            attempt = _separability_attempt(
                :filter_covariance,
                :entangled,
                true,
                :filter_covariance_matrix_violation,
                result,
                result.message,
            )
            push!(attempts, attempt)
            return _separability_report(attempts)
        end
        push!(
            attempts,
            _separability_attempt(
                :filter_covariance, :unknown, false, nothing, result, result.message
            ),
        )
    end

    if :randomized_subtraction in selected
        result = _separability_randomized_subtraction(
            rng,
            prepared;
            max_subtractions=checked_subtractions,
            max_product_restarts=checked_restarts,
            max_product_iterations=checked_product_iterations,
            budget=budget,
        )
        if result.status === :separable
            attempt = _separability_attempt(
                :randomized_subtraction,
                :separable,
                true,
                :validated_product_subtraction_and_separable_ball,
                result,
                result.message,
            )
            push!(attempts, attempt)
            return _separability_report(attempts)
        end
        push!(
            attempts,
            _separability_attempt(
                :randomized_subtraction, :unknown, false, nothing, result, result.message
            ),
        )
    end

    if :symmetric_extension in selected
        for order in orders
            if order > checked_max_extension
                push!(
                    attempts,
                    _separability_attempt(
                        :symmetric_extension,
                        :unknown,
                        false,
                        nothing,
                        (order, max_extension_order=checked_max_extension),
                        "the requested outer hierarchy order exceeds max_extension_order";
                        backend=:optimization,
                    ),
                )
                continue
            end
            result = symmetric_extension(
                prepared.matrix;
                dims=prepared.dimensions,
                order=order,
                ppt=extension_ppt,
                bosonic=extension_bosonic,
                backend=backend,
                prefer_analytic=false,
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
                max_dense_entries=max_dense_entries,
                max_order=checked_max_extension,
                limits=limits,
            )
            if result.verdict === false &&
                result.witness !== nothing &&
                result.witness.entanglement_witness
                attempt = _separability_attempt(
                    :symmetric_extension,
                    :entangled,
                    true,
                    :validated_symmetric_extension_separator,
                    result,
                    result.message;
                    backend=:optimization,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :symmetric_extension,
                    :unknown,
                    false,
                    nothing,
                    result,
                    if result.verdict === true
                        "an outer symmetric extension is only a necessary separability condition"
                    else
                        result.message
                    end;
                    backend=:optimization,
                ),
            )
        end
    end

    if :symmetric_inner_extension in selected
        for order in orders
            if order > checked_max_extension
                push!(
                    attempts,
                    _separability_attempt(
                        :symmetric_inner_extension,
                        :unknown,
                        false,
                        nothing,
                        (order, max_extension_order=checked_max_extension),
                        "the requested inner hierarchy order exceeds max_extension_order";
                        backend=:optimization,
                    ),
                )
                continue
            end
            result = symmetric_inner_extension(
                prepared.matrix;
                dims=prepared.dimensions,
                order=order,
                ppt=extension_ppt,
                backend=backend,
                atol=prepared.analysis.atol,
                rtol=prepared.analysis.rtol,
                allow_densify=true,
                max_dense_entries=max_dense_entries,
                max_order=checked_max_extension,
                limits=limits,
            )
            if result.verdict === true
                attempt = _separability_attempt(
                    :symmetric_inner_extension,
                    :separable,
                    true,
                    :validated_inner_separable_cone_membership,
                    result,
                    result.message;
                    backend=:optimization,
                )
                push!(attempts, attempt)
                return _separability_report(attempts)
            end
            push!(
                attempts,
                _separability_attempt(
                    :symmetric_inner_extension,
                    :unknown,
                    false,
                    nothing,
                    result,
                    if result.verdict === false
                        "failure of an inner-cone test does not prove entanglement"
                    else
                        result.message
                    end;
                    backend=:optimization,
                ),
            )
        end
    end
    return _separability_report(attempts)
end

"""
    is_separable(rho, dims; strategies=:qetlab_deterministic, ...)
    is_separable(rng::AbstractRNG, rho, dims; strategies=:full, ...)

Run an ordered, certificate-first reconstruction of QETLAB's composite
separability analysis. The result is an [`EntanglementReport`](@ref), never an
unsafe Boolean. Necessary tests, hierarchy membership without a separating
dual, exhausted budgets, and numerical boundaries remain `:unknown`.

The represented trace must equal one exactly. An input accepted as normalized
only through `atol`/`rtol` is rejected before any certificate branch.

The no-RNG method rejects `:randomized_subtraction`. The RNG method is the only
route that can run that bounded heuristic; it never touches Julia's default
random stream. Inputs must already be normalized density matrices and are
never symmetrized, normalized, clipped, or otherwise repaired.
"""
function is_separable(rho::AbstractMatrix{<:Number}, dims; kwargs...)
    return _is_separable(nothing, rho, dims; kwargs...)
end

function is_separable(rng::AbstractRNG, rho::AbstractMatrix{<:Number}, dims; kwargs...)
    return _is_separable(rng, rho, dims; kwargs...)
end

function _upbsep_common_type(local_factors)
    types = map(eltype, local_factors)
    promoted = promote_type(types...)
    R = typeof(real(zero(promoted)))
    R in (Float32, Float64) || throw(
        ArgumentError(
            "UPB replacement-vector optimization currently requires Float32, Float64, ComplexF32, or ComplexF64 factors; got promoted type $promoted. No implicit precision-changing conversion is performed.",
        ),
    )
    return R
end

function _upbsep_validate_factors(local_factors; max_dense_entries, max_states)
    checked_max_states = _sepopt_positive_int(max_states, "max_states")
    local_factors isa Tuple ||
        local_factors isa AbstractVector ||
        throw(ArgumentError("local_factors must be a tuple or vector of matrices"))
    local_factors isa AbstractVector && Base.require_one_based_indexing(local_factors)
    length(local_factors) >= 2 ||
        throw(ArgumentError("UPB discrimination requires at least two parties"))
    all(factor -> factor isa AbstractMatrix{<:Number}, local_factors) ||
        throw(ArgumentError("every local factor must be a numeric matrix"))
    for factor in local_factors
        Base.require_one_based_indexing(factor)
        issparse(factor) && throw(
            ArgumentError(
                "replacement-vector SVDs require dense local factors; pass owned dense matrices explicitly",
            ),
        )
        all(isfinite, factor) ||
            throw(ArgumentError("local UPB factors must contain only finite entries"))
    end
    state_count = size(first(local_factors), 2)
    state_count > 0 || throw(ArgumentError("the UPB must contain at least one state"))
    all(factor -> size(factor, 2) == state_count, local_factors) || throw(
        DimensionMismatch("all local-factor matrices must have the same column count")
    )
    dimensions = Tuple(size(factor, 1) for factor in local_factors)
    all(>(1), dimensions) ||
        throw(ArgumentError("every local dimension must be at least two"))
    global_dimension = _checked_product(
        dimensions, "the UPB global Hilbert-space dimension"
    )
    dense_limit = _sepopt_limit(max_dense_entries, "max_dense_entries")
    R = _upbsep_common_type(local_factors)
    if state_count > checked_max_states
        return (;
            factors=nothing,
            dimensions,
            global_dimension,
            state_count,
            real_type=R,
            max_states=checked_max_states,
            state_limit_exceeded=true,
        )
    end
    input_entries = sum(BigInt(length(factor)) for factor in local_factors)
    dense_limit !== nothing &&
        input_entries > dense_limit &&
        throw(
            ArgumentError(
                "UPB local-factor conversion needs $input_entries dense entries, " *
                "exceeding max_dense_entries=$dense_limit",
            ),
        )
    _sepopt_guard_dense(global_dimension, max_dense_entries, "UPB replacement projections")
    factors = Tuple(Matrix{Complex{R}}(factor) for factor in local_factors)
    return (;
        factors,
        dimensions,
        global_dimension,
        state_count,
        real_type=R,
        max_states=checked_max_states,
        state_limit_exceeded=false,
    )
end

function _upbsep_state_limit_analysis(
    data; atol, rtol, boundary_factor, max_upb_partitions, max_upb_work, max_dense_entries
)
    R = data.real_type
    options = _upb_options(
        Complex{R};
        atol=atol === nothing ? 0 : atol,
        rtol,
        boundary_factor,
        max_partitions=max_upb_partitions,
        max_work=max_upb_work,
        max_dense_entries,
        max_states=data.max_states,
        normalization=:require,
        allow_densify=false,
        materialize_witness=false,
    )
    metadata = (;
        dimensions=data.dimensions,
        state_count=data.state_count,
        input_form=:local_factors,
        options,
        input_normalized=nothing,
        analysis_rescaled=false,
        densified=false,
        product_residual=nothing,
    )
    return _upb_make_result(
        metadata;
        status=:unknown,
        reason=:state_limit,
        is_upb_value=nothing,
        message="the input has $(data.state_count) states, exceeding " *
                "max_states=$(data.max_states)",
    )
end

function _upbsep_partition_count(dimensions, state_count::Int)
    capacities = Tuple(dimension - 1 for dimension in dimensions)
    remaining = state_count - 1
    capacity_sum = sum(capacities)
    remaining == capacity_sum || return (
        supported=false,
        count=BigInt(0),
        capacities,
        message="the pinned finite replacement-vector enumeration requires a minimal UPB with state_count - 1 == sum(dimensions .- 1); got $remaining and $capacity_sum",
    )
    count = factorial(BigInt(remaining))
    for capacity in capacities
        count ÷= factorial(BigInt(capacity))
    end
    return (
        supported=true,
        count=BigInt(state_count) * count,
        capacities,
        message="finite isolated replacement-vector partitions",
    )
end

mutable struct _UPBReplacementSearch
    partitions_examined::Int
    candidates_generated::Int
    work_used::BigInt
    stopped::Union{Nothing,Symbol}
    boundary_seen::Bool
    nonisolated_seen::Bool
end

function _upbsep_reserve_work!(
    search::_UPBReplacementSearch, amount, max_work, label::AbstractString
)
    amount_big = BigInt(amount)
    amount_big >= 0 || error("internal UPB work estimate for $label is negative")
    required = search.work_used + amount_big
    if max_work !== nothing && required > max_work
        search.stopped = :work_limit
        return false
    end
    search.work_used = required
    return true
end

function _upbsep_null_factor(matrix, tolerance, boundary_factor)
    dimension, columns = size(matrix)
    columns == dimension - 1 || return (
        status=:nonisolated, vector=nothing, singular_values=nothing, margin=nothing
    )
    decomposition = svd(matrix; full=true)
    values = decomposition.S
    scale = maximum(values; init=one(tolerance))
    threshold = tolerance * max(one(tolerance), scale)
    boundary = boundary_factor * threshold
    minimum_value = minimum(values)
    if minimum_value <= threshold
        return (
            status=:nonisolated,
            vector=nothing,
            singular_values=values,
            margin=minimum_value - threshold,
        )
    elseif minimum_value <= boundary
        return (
            status=:boundary,
            vector=nothing,
            singular_values=values,
            margin=minimum_value - threshold,
        )
    end
    vector = decomposition.U[:, end]
    phase_index = findfirst(value -> !iszero(value), vector)
    if phase_index !== nothing
        phase = vector[phase_index] / abs(vector[phase_index])
        vector ./= phase
    end
    return (
        status=:isolated, vector, singular_values=values, margin=minimum_value - threshold
    )
end

function _upbsep_duplicate(candidate, replacements, tolerance)
    for replacement in replacements
        overlap = abs2(dot(replacement.global_vector, candidate))
        abs(one(overlap) - overlap) <= tolerance && return true
    end
    return false
end

function _upbsep_evaluate_partition!(
    replacements,
    search::_UPBReplacementSearch,
    factors,
    removed::Int,
    groups,
    tolerance,
    boundary_factor,
    max_candidates,
    max_work,
)
    rank_work = sum(BigInt(size(factors[party], 1))^3 for party in eachindex(factors))
    global_dimension = prod(BigInt(size(factor, 1)) for factor in factors)
    candidate_work = global_dimension^2 + BigInt(length(replacements)) * global_dimension
    _upbsep_reserve_work!(
        search,
        rank_work + candidate_work,
        max_work,
        "partition SVD, candidate construction, and duplicate comparisons",
    ) || return nothing
    search.partitions_examined += 1
    local_vectors = Vector{Vector{eltype(first(factors))}}(undef, length(factors))
    margins = Any[]
    for party in eachindex(factors)
        span = @view factors[party][:, groups[party]]
        result = _upbsep_null_factor(span, tolerance, boundary_factor)
        push!(margins, result)
        if result.status === :boundary
            search.boundary_seen = true
            return nothing
        elseif result.status === :nonisolated
            search.nonisolated_seen = true
            return nothing
        end
        local_vectors[party] = Vector(result.vector)
    end
    global_vector = foldl(kron, local_vectors)
    global_vector ./= norm(global_vector)
    residual = zero(tolerance)
    for state_index in axes(first(factors), 2)
        state_index == removed && continue
        overlap = prod(
            dot(local_vectors[party], @view(factors[party][:, state_index])) for
            party in eachindex(factors)
        )
        residual = max(residual, abs(overlap))
    end
    residual <= boundary_factor * tolerance || begin
        search.boundary_seen = true
        return nothing
    end
    _upbsep_duplicate(global_vector, replacements, boundary_factor * tolerance) &&
        return nothing
    if max_candidates !== nothing && search.candidates_generated >= max_candidates
        search.stopped = :candidate_limit
        return nothing
    end
    projection = global_vector * adjoint(global_vector)
    replacement = UPBReplacementVector(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        removed,
        Tuple(local_vectors),
        global_vector,
        projection,
        Tuple(Tuple(group) for group in groups),
        residual,
        Tuple(margin.margin for margin in margins),
    )
    push!(replacements, replacement)
    search.candidates_generated += 1
    return nothing
end

function _upbsep_assign_partitions!(
    replacements,
    search,
    factors,
    removed,
    remaining,
    capacities,
    groups,
    position,
    tolerance,
    boundary_factor,
    max_partitions,
    max_candidates,
    max_work,
    node_precharged::Bool=false,
)
    search.stopped === nothing || return nothing
    # Charge every recursion node before inspecting it. This makes deep or
    # highly branching partition searches obey max_work even before they
    # reach a leaf where the SVD/candidate estimate is reserved.
    if !node_precharged
        _upbsep_reserve_work!(search, 1, max_work, "partition-search recursion node") ||
            return nothing
    end
    if position > length(remaining)
        if max_partitions !== nothing && search.partitions_examined >= max_partitions
            search.stopped = :partition_limit
            return nothing
        end
        return _upbsep_evaluate_partition!(
            replacements,
            search,
            factors,
            removed,
            groups,
            tolerance,
            boundary_factor,
            max_candidates,
            max_work,
        )
    end
    state_index = remaining[position]
    for party in eachindex(groups)
        length(groups[party]) < capacities[party] || continue
        push!(groups[party], state_index)
        remaining_after = length(remaining) - position
        feasible = all(
            length(groups[index]) <= capacities[index] &&
                length(groups[index]) + remaining_after >= capacities[index] for
            index in eachindex(groups)
        )
        feasible && _upbsep_assign_partitions!(
            replacements,
            search,
            factors,
            removed,
            remaining,
            capacities,
            groups,
            position + 1,
            tolerance,
            boundary_factor,
            max_partitions,
            max_candidates,
            max_work,
        )
        pop!(groups[party])
        search.stopped === nothing || return nothing
    end
    return nothing
end

function _upbsep_replacement_vectors(
    data; atol, rtol, boundary_factor, max_partitions, max_candidates, max_work
)
    R = data.real_type
    absolute = atol === nothing ? zero(R) : convert(R, atol)
    relative = rtol === nothing ? 8 * eps(R) : convert(R, rtol)
    factor = convert(R, boundary_factor)
    all(isfinite, (absolute, relative, factor)) &&
    absolute >= zero(R) &&
    relative >= zero(R) &&
    factor > one(R) || throw(
        ArgumentError(
            "atol and rtol must be finite nonnegative values and boundary_factor must exceed one",
        ),
    )
    tolerance = absolute + relative
    partition_info = _upbsep_partition_count(data.dimensions, data.state_count)
    partition_info.supported || return (
        status=:unsupported_nonminimal,
        replacements=UPBReplacementVector[],
        partition_info,
        search=_UPBReplacementSearch(0, 0, BigInt(0), nothing, false, true),
        tolerance,
        message=partition_info.message,
    )
    partition_limit = _sepopt_limit(max_partitions, "max_partitions")
    candidate_limit = _sepopt_limit(max_candidates, "max_candidates")
    work_limit = _sepopt_limit(max_work, "max_work")
    replacements = UPBReplacementVector[]
    search = _UPBReplacementSearch(0, 0, BigInt(0), nothing, false, false)
    capacities = partition_info.capacities
    for removed in 1:data.state_count
        # Reserve the root node before allocating its state-index/group
        # workspaces, so max_work=0 stops without entering the recursion.
        _upbsep_reserve_work!(
            search, 1, work_limit, "partition-search root recursion node"
        ) || break
        remaining = [index for index in 1:data.state_count if index != removed]
        groups = [Int[] for _ in data.dimensions]
        _upbsep_assign_partitions!(
            replacements,
            search,
            data.factors,
            removed,
            remaining,
            capacities,
            groups,
            1,
            tolerance,
            factor,
            partition_limit,
            candidate_limit,
            work_limit,
            true,
        )
        search.stopped === nothing || break
    end
    status = if search.stopped !== nothing
        search.stopped
    elseif search.boundary_seen
        :numerical_boundary
    elseif search.nonisolated_seen
        :nonisolated_replacement_set
    elseif isempty(replacements)
        :no_replacement_vectors
    else
        :complete
    end
    return (
        status,
        replacements,
        partition_info,
        search,
        tolerance,
        message=if status === :complete
            "all finite isolated replacement-vector partitions were exhausted"
        else
            "replacement-vector enumeration ended with status $status"
        end,
    )
end

function _upbsep_build_problem(
    analysis, generation, dimensions, state_count, limits::OptimizationLimits
)
    replacements = generation.replacements
    count = length(replacements)
    count > 0 || throw(ArgumentError("at least one replacement vector is required"))
    T = typeof(generation.tolerance)
    dimension = _checked_product(dimensions, "UPB reconstruction dimensions")
    equality_count = BigInt(dimension)^2
    count <= limits.max_variables || throw(
        ArgumentError(
            "model needs $count variables, exceeding max_variables=$(limits.max_variables)",
        ),
    )
    equality_count <= limits.max_equalities || throw(
        ArgumentError(
            "model needs $equality_count equalities, exceeding max_equalities=$(limits.max_equalities)",
        ),
    )
    count <= limits.max_intervals || throw(
        ArgumentError(
            "model needs $count intervals, exceeding max_intervals=$(limits.max_intervals)",
        ),
    )
    variables = collect(1:count)
    coefficients = [replacement.projection for replacement in replacements]
    reconstruction = HermitianAffineMatrix(
        :upb_reconstruction,
        zeros(Complex{T}, dimension, dimension),
        variables,
        coefficients,
        count,
    )
    equalities = hermitian_equalities(
        reconstruction;
        target=Matrix{Complex{T}}(I, dimension, dimension),
        name_prefix=:upb_identity_reconstruction,
    )
    intervals = AffineInterval{T}[
        AffineInterval(
            AffineScalar(zero(T), [index], [one(T)], count),
            zero(T),
            nothing,
            Symbol(:upb_coefficient_nonnegative_, index),
        ) for index in 1:count
    ]
    program = SemidefiniteProgram(
        :upb_separable_discrimination,
        :feasibility,
        count,
        AffineScalar(zero(T), zeros(T, count));
        equalities,
        intervals,
        primal_views=[reconstruction],
        limits,
        metadata=(
            formulation=:upb_nonnegative_identity_reconstruction,
            dimensions,
            state_count,
            candidate_count=count,
            source=:bandyopadhyay_cosentino_johnston_russo_watrous_yu,
        ),
    )
    return UPBSeparableDiscriminationProblem(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        analysis,
        Tuple(replacements),
        program,
        dimensions,
        state_count,
        count,
    )
end

function _upbsep_prepare(
    local_factors;
    atol,
    rtol,
    boundary_factor,
    max_partitions,
    max_candidates,
    max_work,
    max_dense_entries,
    max_states,
    max_upb_partitions,
    max_upb_work,
)
    # Validate the replacement-search caps even when an earlier input-size
    # preflight returns a structured state-limit result.
    _sepopt_limit(max_partitions, "max_partitions")
    _sepopt_limit(max_candidates, "max_candidates")
    _sepopt_limit(max_work, "max_work")
    data = _upbsep_validate_factors(
        local_factors; max_dense_entries=max_dense_entries, max_states=max_states
    )
    if data.state_limit_exceeded
        analysis = _upbsep_state_limit_analysis(
            data;
            atol=atol,
            rtol=rtol,
            boundary_factor=boundary_factor,
            max_upb_partitions=max_upb_partitions,
            max_upb_work=max_upb_work,
            max_dense_entries=max_dense_entries,
        )
        return (;
            data,
            analysis,
            generation=nothing,
            message="the supplied local factors were not certified as a UPB: " *
                    analysis.message,
        )
    end
    analysis = is_upb(
        data.factors;
        atol=atol === nothing ? 0 : atol,
        rtol=rtol,
        boundary_factor=boundary_factor,
        max_partitions=max_upb_partitions,
        max_work=max_upb_work,
        max_dense_entries=max_dense_entries,
        max_states=max_states,
        normalization=:require,
        allow_densify=false,
        materialize_witness=false,
    )
    analysis.status === :upb || return (
        data,
        analysis,
        generation=nothing,
        message="the supplied local factors were not certified as a UPB: $(analysis.message)",
    )
    generation = _upbsep_replacement_vectors(
        data;
        atol=atol,
        rtol=rtol,
        boundary_factor=boundary_factor,
        max_partitions=max_partitions,
        max_candidates=max_candidates,
        max_work=max_work,
    )
    return (; data, analysis, generation, message=generation.message)
end

"""
    upb_sep_distinguishability_problem(local_factors; kwargs...)

Build the finite nonnegative identity-reconstruction problem characterizing
perfect separable discrimination of a minimal UPB. `local_factors[p][:,j]` is
party `p`'s local factor in state `j`. The input must first pass the package's
full [`is_upb`](@ref) analysis.

The pinned finite algorithm applies to minimal UPBs with
`state_count - 1 == sum(dimensions .- 1)`. Nonminimal or numerically
nonisolated replacement sets are reported by the high-level API rather than
being silently approximated by an incomplete finite list.
"""
function upb_sep_distinguishability_problem(
    local_factors;
    atol=nothing,
    rtol=nothing,
    boundary_factor=8,
    max_partitions=100_000,
    max_candidates=100_000,
    max_work=100_000_000,
    max_dense_entries=2_000_000,
    max_states=256,
    max_upb_partitions=100_000,
    max_upb_work=25_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _upbsep_prepare(
        local_factors;
        atol=atol,
        rtol=rtol,
        boundary_factor=boundary_factor,
        max_partitions=max_partitions,
        max_candidates=max_candidates,
        max_work=max_work,
        max_dense_entries=max_dense_entries,
        max_states=max_states,
        max_upb_partitions=max_upb_partitions,
        max_upb_work=max_upb_work,
    )
    prepared.analysis.status === :upb || throw(ArgumentError(prepared.message))
    prepared.generation.status === :complete || throw(ArgumentError(prepared.message))
    return _upbsep_build_problem(
        prepared.analysis,
        prepared.generation,
        prepared.data.dimensions,
        prepared.data.state_count,
        limits,
    )
end

function upb_sep_distinguishability_problem(
    first_factor::AbstractMatrix,
    second_factor::AbstractMatrix,
    remaining_factors::AbstractMatrix...;
    kwargs...,
)
    return upb_sep_distinguishability_problem(
        (first_factor, second_factor, remaining_factors...); kwargs...
    )
end

function _upbsep_result(
    status,
    feasibility,
    prepared;
    problem=nothing,
    coefficients=nothing,
    reconstruction=nothing,
    reconstruction_residual=nothing,
    minimum_coefficient=nothing,
    dual_certificate=nothing,
    optimization_result=nothing,
    message,
)
    generation = prepared.generation
    R = prepared.data.real_type
    tolerance = generation === nothing ? zero(R) : convert(R, generation.tolerance)
    search = if generation === nothing
        _UPBReplacementSearch(0, 0, BigInt(0), nothing, false, false)
    else
        generation.search
    end
    replacements = generation === nothing ? () : Tuple(generation.replacements)
    convert_optional(value) = value === nothing ? nothing : convert(R, value)
    return UPBSeparableDiscriminationResult(
        _VALIDATED_CONSTRUCTOR_TOKEN,
        status,
        feasibility,
        nothing,
        false,
        nothing,
        prepared.analysis,
        replacements,
        problem,
        coefficients,
        reconstruction,
        convert_optional(reconstruction_residual),
        convert_optional(minimum_coefficient),
        dual_certificate,
        optimization_result,
        search.partitions_examined,
        search.candidates_generated,
        search.work_used,
        tolerance,
        String(message),
    )
end

function _upbsep_validate_primal(problem, optimization, tolerance)
    optimization.primal === nothing && return (
        valid=false,
        coefficients=nothing,
        reconstruction=nothing,
        residual=nothing,
        minimum_coefficient=nothing,
        exact=false,
    )
    coefficients = copy(optimization.primal.coordinates)
    reconstruction = zeros(
        Complex{eltype(coefficients)}, prod(problem.dimensions), prod(problem.dimensions)
    )
    for (coefficient, replacement) in zip(coefficients, problem.replacement_vectors)
        reconstruction += coefficient * replacement.projection
    end
    identity_matrix = Matrix{eltype(reconstruction)}(
        I, size(reconstruction, 1), size(reconstruction, 2)
    )
    residual = norm(reconstruction - identity_matrix)
    minimum_coefficient = minimum(coefficients)
    valid =
        minimum_coefficient >= -tolerance &&
        residual <= tolerance &&
        all(isfinite, coefficients)
    exact =
        minimum_coefficient >= zero(minimum_coefficient) &&
        all(iszero, reconstruction - identity_matrix)
    return (
        valid=valid,
        coefficients=coefficients,
        reconstruction=reconstruction,
        residual=residual,
        minimum_coefficient=minimum_coefficient,
        exact=exact,
    )
end

function _upbsep_farkas_certificate(problem, optimization, tolerance)
    optimization.dual === nothing && return nothing
    y = optimization.dual.equalities
    length(y) == length(problem.program.equalities) || return nothing
    scale = norm(y)
    isfinite(scale) && scale > zero(scale) || return nothing
    normalized = y / scale
    variable_count = problem.program.variable_count
    coefficients = zeros(eltype(normalized), variable_count)
    constant = zero(eltype(normalized))
    for (weight, equality) in zip(normalized, problem.program.equalities)
        constant += weight * equality.function_data.constant
        indices, values = findnz(equality.function_data.coefficients)
        for (index, value) in zip(indices, values)
            coefficients[index] += weight * value
        end
    end
    positive_orientation = minimum(coefficients) > tolerance && constant > tolerance
    negative_orientation = maximum(coefficients) < -tolerance && constant < -tolerance
    valid = positive_orientation || negative_orientation
    operator = try
        _symext_hermitian_from_duals(normalized, prod(problem.dimensions))
    catch
        nothing
    end
    return (
        valid,
        orientation=if positive_orientation
            :positive
        else
            (negative_orientation ? :negative : :none)
        end,
        equality_dual=copy(y),
        normalized_equality_dual=normalized,
        cone_coefficients=coefficients,
        constant,
        operator,
        tolerance,
        message=if valid
            "the normalized Farkas functional strictly separates the identity from the computed floating replacement cone; no perturbation bound links it to the exact input-UPB cone"
        else
            "the optimizer dual did not pass the package-owned Farkas sign margins"
        end,
    )
end

"""
    upb_sep_distinguishable(local_factors;
                            backend=NoOptimizationBackend(), ...)

Analyze perfect separable distinguishability of a minimal UPB. Replacement
vectors are enumerated lazily with explicit partition, candidate, dense-entry,
and scalar-work limits. The work counter reserves every recursion node before
visiting it, then separately reserves the partition SVD, candidate-construction,
and duplicate-comparison estimates before those operations. Complex
orthogonality uses the Hilbert-space adjoint, correcting the pinned routine's
nonconjugating transpose.

The result never converts a small least-squares residual into a Boolean.
Numerical primal feasibility is retained as `:numerically_feasible`. Even a
bitwise reconstruction of the computed floating projectors remains numerical
unless every replacement vector is proved exactly orthogonal to the represented
input and the input UPB itself has an exact certificate. A Farkas separator for
the computed floating cone is retained as diagnostic evidence, not promoted to
a false Boolean without a rigorous perturbation margin.
"""
function upb_sep_distinguishable(
    local_factors;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    boundary_factor=8,
    max_partitions=100_000,
    max_candidates=100_000,
    max_work=100_000_000,
    max_dense_entries=2_000_000,
    max_states=256,
    max_upb_partitions=100_000,
    max_upb_work=25_000_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _upbsep_prepare(
        local_factors;
        atol=atol,
        rtol=rtol,
        boundary_factor=boundary_factor,
        max_partitions=max_partitions,
        max_candidates=max_candidates,
        max_work=max_work,
        max_dense_entries=max_dense_entries,
        max_states=max_states,
        max_upb_partitions=max_upb_partitions,
        max_upb_work=max_upb_work,
    )
    if prepared.analysis.status !== :upb
        resource_reason =
            prepared.analysis.reason in (:state_limit, :partition_limit, :work_limit)
        status = if resource_reason
            :resource_limit
        elseif prepared.analysis.status === :unknown
            :unknown_input
        else
            :invalid_input
        end
        feasibility = resource_reason ? :unknown : :not_analyzed
        return _upbsep_result(status, feasibility, prepared; message=prepared.message)
    end
    generation = prepared.generation
    if generation.status !== :complete
        status = if generation.status in (:partition_limit, :candidate_limit, :work_limit)
            :resource_limit
        else
            generation.status
        end
        return _upbsep_result(status, :unknown, prepared; message=generation.message)
    end
    problem = try
        _upbsep_build_problem(
            prepared.analysis,
            generation,
            prepared.data.dimensions,
            prepared.data.state_count,
            limits,
        )
    catch error
        if error isa ArgumentError && occursin("exceed", sprint(showerror, error))
            return _upbsep_result(
                :resource_limit, :unknown, prepared; message=sprint(showerror, error)
            )
        end
        rethrow()
    end
    optimization = solve_optimization(problem.program, backend)
    if optimization.status === OptimizationBackendUnavailable
        return _upbsep_result(
            :backend_unavailable,
            :unknown,
            prepared;
            problem,
            optimization_result=optimization,
            message=optimization.message,
        )
    elseif optimization.status === OptimizationLimit
        return _upbsep_result(
            :resource_limit,
            :unknown,
            prepared;
            problem,
            optimization_result=optimization,
            message=optimization.message,
        )
    end
    tolerance = max(
        generation.tolerance,
        if backend isa JuMPBackend
            convert(prepared.data.real_type, backend.atol + backend.rtol)
        else
            generation.tolerance
        end,
    )
    if optimization.status === OptimizationInfeasible
        certificate = _upbsep_farkas_certificate(problem, optimization, tolerance)
        return _upbsep_result(
            :infeasibility_unverified,
            :unknown,
            prepared;
            problem,
            dual_certificate=certificate,
            optimization_result=optimization,
            message=if certificate !== nothing && certificate.valid
                certificate.message
            else
                "the backend reported infeasibility without a checked separator for even the computed replacement cone"
            end,
        )
    elseif optimization.status in (
        OptimizationMalformedBackend,
        OptimizationUnsupported,
        OptimizationNumericalFailure,
        OptimizationInconsistent,
        OptimizationUnknown,
        OptimizationUnbounded,
    )
        return _upbsep_result(
            :backend_failure,
            :unknown,
            prepared;
            problem,
            optimization_result=optimization,
            message=optimization.message,
        )
    end
    primal = _upbsep_validate_primal(problem, optimization, tolerance)
    if !primal.valid
        return _upbsep_result(
            :invalid_primal,
            :unknown,
            prepared;
            problem,
            coefficients=primal.coefficients,
            reconstruction=primal.reconstruction,
            reconstruction_residual=primal.residual,
            minimum_coefficient=primal.minimum_coefficient,
            optimization_result=optimization,
            message="the backend primal failed nonnegativity or identity-reconstruction residual checks",
        )
    end
    return _upbsep_result(
        :numerically_feasible,
        :numerically_feasible,
        prepared;
        problem,
        coefficients=primal.coefficients,
        reconstruction=primal.reconstruction,
        reconstruction_residual=primal.residual,
        minimum_coefficient=primal.minimum_coefficient,
        optimization_result=optimization,
        message=if primal.exact
            "the computed floating replacement projectors reconstruct the identity bitwise, but the floating-only public input and generator domain has no exact linked proof path; the result remains numerical"
        else
            "a residual-checked numerical reconstruction is retained without converting its tolerance to a Boolean theorem"
        end,
    )
end

function upb_sep_distinguishable(
    first_factor::AbstractMatrix,
    second_factor::AbstractMatrix,
    remaining_factors::AbstractMatrix...;
    kwargs...,
)
    return upb_sep_distinguishable(
        (first_factor, second_factor, remaining_factors...); kwargs...
    )
end
