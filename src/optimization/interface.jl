# Independently designed solver-neutral optimization contracts. Optional
# modeling layers translate these immutable package-owned values without adding
# a solver dependency to the core package.

"""
    OptimizationStatus

Package-owned classification of an optimization attempt. Backend status strings
are retained separately in [`OptimizationResult`](@ref); this enum never turns
an inaccurate, limited, or failed solve into a mathematical yes/no conclusion.
"""
@enum OptimizationStatus begin
    OptimizationOptimal
    OptimizationFeasible
    OptimizationInfeasible
    OptimizationUnbounded
    OptimizationLimit
    OptimizationNumericalFailure
    OptimizationUnsupported
    OptimizationBackendUnavailable
    OptimizationMalformedBackend
    OptimizationInconsistent
    OptimizationUnknown
end

"""
    AbstractOptimizationBackend

Supertype for explicit optimization backend configurations. There is no global
or implicit default solver.
"""
abstract type AbstractOptimizationBackend end

"""
    NoOptimizationBackend()

Dependency-free backend marker. Solving with this marker returns an
[`OptimizationBackendUnavailable`](@ref OptimizationStatus) result.
"""
struct NoOptimizationBackend <: AbstractOptimizationBackend end

"""
    JuMPBackend(optimizer_factory; kwargs...)

Package-owned configuration for the optional JuMP extension.

`optimizer_factory` is passed to JuMP only after the extension is loaded.
`optimizer_options` is a `NamedTuple` or tuple of `Pair`s passed through
JuMP's public optimizer-attribute API. `time_limit_seconds` is an optional
positive finite wall-clock limit. `allow_densify=true` is required for an SDP:
the portable complex-Hermitian convention uses explicit real block PSD
matrices. `atol` and `rtol` control package-owned residual checks.

`optimizer_name` and `optimizer_version` are user-declared provenance hints.
The extension also records the solver name reported through MOI.
"""
struct JuMPBackend{F,O,T<:Real} <: AbstractOptimizationBackend
    optimizer_factory::F
    optimizer_options::O
    optimizer_name::String
    optimizer_version::Union{Nothing,VersionNumber}
    time_limit_seconds::Union{Nothing,Float64}
    silent::Bool
    allow_densify::Bool
    atol::T
    rtol::T
end

function _optimization_nonnegative_real(value, name::AbstractString)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be a finite nonnegative real number"))
    isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("$name must be a finite nonnegative real number"))
    return value
end

function _optimization_positive_int(value, name::AbstractString)
    value isa Integer && !(value isa Bool) && value > 0 ||
        throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _optimization_options(options)
    if options isa NamedTuple
        return options
    elseif options isa Tuple && all(option -> option isa Pair, options)
        names = Symbol[]
        values = Any[]
        for option in options
            key = first(option)
            key isa Union{Symbol,AbstractString} || throw(
                ArgumentError(
                    "optimizer option names must be Symbols or strings, got $(repr(key))",
                ),
            )
            symbol = Symbol(key)
            symbol in names &&
                throw(ArgumentError("duplicate optimizer option $(repr(symbol))"))
            push!(names, symbol)
            push!(values, last(option))
        end
        return NamedTuple{Tuple(names)}(Tuple(values))
    end
    return throw(
        ArgumentError("optimizer_options must be a NamedTuple or a tuple of Pair values")
    )
end

function JuMPBackend(
    optimizer_factory;
    optimizer_options=NamedTuple(),
    optimizer_name::AbstractString="unspecified",
    optimizer_version::Union{Nothing,VersionNumber}=nothing,
    time_limit_seconds=nothing,
    silent::Bool=true,
    allow_densify::Bool=false,
    atol=1.0e-7,
    rtol=1.0e-7,
)
    isempty(strip(optimizer_name)) &&
        throw(ArgumentError("optimizer_name must not be empty"))
    options = _optimization_options(optimizer_options)
    time_limit = if time_limit_seconds === nothing
        nothing
    else
        time_limit_seconds isa Real && !(time_limit_seconds isa Bool) || throw(
            ArgumentError(
                "time_limit_seconds must be nothing or a finite positive real number",
            ),
        )
        isfinite(time_limit_seconds) && time_limit_seconds > zero(time_limit_seconds) ||
            throw(
                ArgumentError(
                    "time_limit_seconds must be nothing or a finite positive real number",
                ),
            )
        converted = try
            Float64(time_limit_seconds)
        catch error
            error isa InexactError || rethrow()
            throw(ArgumentError("time_limit_seconds cannot be represented as Float64"))
        end
        isfinite(converted) && converted > 0 || throw(
            ArgumentError("time_limit_seconds must remain finite and positive as Float64"),
        )
        converted
    end
    absolute = _optimization_nonnegative_real(atol, "atol")
    relative = _optimization_nonnegative_real(rtol, "rtol")
    promoted_absolute, promoted_relative = promote(absolute, relative)
    return JuMPBackend(
        optimizer_factory,
        options,
        String(optimizer_name),
        optimizer_version,
        time_limit,
        silent,
        allow_densify,
        promoted_absolute,
        promoted_relative,
    )
end

"""
    OptimizationLimits(; kwargs...)

Deterministic pre-allocation limits for a solver-neutral conic model. Counts are
checked with `BigInt` arithmetic before the optional extension allocates a
JuMP model or a real block PSD matrix.
"""
struct OptimizationLimits
    max_variables::Int
    max_equalities::Int
    max_intervals::Int
    max_psd_blocks::Int
    max_psd_dimension::Int
    max_model_entries::Int
end

function OptimizationLimits(;
    max_variables=10_000,
    max_equalities=100_000,
    max_intervals=100_000,
    max_psd_blocks=1_000,
    max_psd_dimension=512,
    max_model_entries=10_000_000,
)
    return OptimizationLimits(
        _optimization_positive_int(max_variables, "max_variables"),
        _optimization_positive_int(max_equalities, "max_equalities"),
        _optimization_positive_int(max_intervals, "max_intervals"),
        _optimization_positive_int(max_psd_blocks, "max_psd_blocks"),
        _optimization_positive_int(max_psd_dimension, "max_psd_dimension"),
        _optimization_positive_int(max_model_entries, "max_model_entries"),
    )
end

"""
    AffineScalar(constant, coefficients)

Sparse scalar affine function `constant + dot(coefficients, x)`. Coefficients
are copied into a one-based `SparseVector`; Boolean, nonfinite, and non-real
data are rejected.
"""
struct AffineScalar{T<:Real}
    constant::T
    coefficients::SparseVector{T,Int}
end

function _optimization_validate_real(value, name::AbstractString)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must contain real non-Boolean values"))
    isfinite(value) || throw(ArgumentError("$name must contain only finite values"))
    return value
end

function AffineScalar(
    constant::Real, coefficients::AbstractVector{<:Real}; variable_count=nothing
)
    _optimization_validate_real(constant, "constant")
    firstindex(coefficients) == 1 ||
        throw(ArgumentError("coefficient vectors must use one-based indexing"))
    count = if isnothing(variable_count)
        length(coefficients)
    else
        _optimization_positive_int(variable_count, "variable_count")
    end
    length(coefficients) == count || throw(
        DimensionMismatch("expected $count coefficients, got $(length(coefficients))")
    )
    for value in coefficients
        _optimization_validate_real(value, "coefficients")
    end
    T = promote_type(typeof(constant), eltype(coefficients))
    T <: Real || throw(ArgumentError("affine coefficients must promote to a real type"))
    converted = sparsevec(
        Int[index for index in eachindex(coefficients) if !iszero(coefficients[index])],
        T[
            coefficients[index] for
            index in eachindex(coefficients) if !iszero(coefficients[index])
        ],
        count,
    )
    return AffineScalar{T}(convert(T, constant), converted)
end

function AffineScalar(
    constant::Real,
    variable_indices::AbstractVector{<:Integer},
    coefficient_values::AbstractVector{<:Real},
    variable_count::Integer,
)
    count = _optimization_positive_int(variable_count, "variable_count")
    firstindex(variable_indices) == 1 && firstindex(coefficient_values) == 1 ||
        throw(ArgumentError("affine inputs must use one-based indexing"))
    length(variable_indices) == length(coefficient_values) || throw(
        DimensionMismatch("variable_indices and coefficient_values must have equal length"),
    )
    _optimization_validate_real(constant, "constant")
    for value in coefficient_values
        _optimization_validate_real(value, "coefficient_values")
    end
    indices = Int[]
    for index in variable_indices
        index isa Bool &&
            throw(ArgumentError("variable indices must be one-based integers"))
        1 <= index <= count ||
            throw(ArgumentError("variable index $index is outside 1:$count"))
        push!(indices, Int(index))
    end
    length(unique(indices)) == length(indices) ||
        throw(ArgumentError("variable_indices must not contain duplicates"))
    T = promote_type(typeof(constant), eltype(coefficient_values))
    T <: Real || throw(ArgumentError("affine coefficients must promote to a real type"))
    keep = findall(!iszero, coefficient_values)
    return AffineScalar{T}(
        convert(T, constant),
        sparsevec(indices[keep], T[coefficient_values[index] for index in keep], count),
    )
end

"""
    AffineEquality(function, name)

Constraint `function(x) == 0`.
"""
struct AffineEquality{T<:Real}
    function_data::AffineScalar{T}
    name::Symbol
end

function AffineEquality(function_data::AffineScalar{T}, name::Symbol=:equality) where {T}
    return AffineEquality{T}(function_data, name)
end

"""
    AffineInterval(function, lower, upper, name)

Constraint `lower <= function(x) <= upper`. Either finite endpoint may be
`nothing`, but not both. Equal finite endpoints must be represented by
[`AffineEquality`](@ref).
"""
struct AffineInterval{T<:Real}
    function_data::AffineScalar{T}
    lower::Union{Nothing,T}
    upper::Union{Nothing,T}
    name::Symbol

    function AffineInterval{T}(
        function_data::AffineScalar{T},
        lower::Union{Nothing,T},
        upper::Union{Nothing,T},
        name::Symbol,
    ) where {T<:Real}
        lower === nothing &&
            upper === nothing &&
            throw(ArgumentError("an affine interval needs at least one endpoint"))
        lower !== nothing &&
            upper !== nothing &&
            lower > upper &&
            throw(ArgumentError("lower must not exceed upper"))
        lower !== nothing &&
            upper !== nothing &&
            lower == upper &&
            throw(
                ArgumentError(
                    "equal interval endpoints are an equality; use AffineEquality"
                ),
            )
        return new{T}(function_data, lower, upper, name)
    end
end

function AffineInterval(
    function_data::AffineScalar{T},
    lower::Union{Nothing,Real},
    upper::Union{Nothing,Real},
    name::Symbol=:interval,
) where {T<:Real}
    lower === nothing &&
        upper === nothing &&
        throw(ArgumentError("an affine interval needs at least one endpoint"))
    checked_lower = if lower === nothing
        nothing
    else
        _optimization_validate_real(lower, "lower")
        convert(T, lower)
    end
    checked_upper = if upper === nothing
        nothing
    else
        _optimization_validate_real(upper, "upper")
        convert(T, upper)
    end
    checked_lower !== nothing &&
        checked_upper !== nothing &&
        checked_lower > checked_upper &&
        throw(ArgumentError("lower must not exceed upper"))
    checked_lower !== nothing &&
        checked_upper !== nothing &&
        checked_lower == checked_upper &&
        throw(ArgumentError("equal interval endpoints are an equality; use AffineEquality"))
    return AffineInterval{T}(function_data, checked_lower, checked_upper, name)
end

"""
    HermitianAffineTerm(variable, coefficient)

One sparse Hermitian coefficient matrix multiplying a real scalar decision
variable.
"""
struct HermitianAffineTerm{T<:Real}
    variable::Int
    coefficient::SparseMatrixCSC{Complex{T},Int}
end

function _optimization_complex_real_type(::Type{T}) where {T<:Real}
    return T
end

function _optimization_complex_real_type(::Type{Complex{T}}) where {T<:Real}
    return T
end

function _optimization_matrix_real_type(matrix::AbstractMatrix{<:Number})
    return _optimization_complex_real_type(eltype(matrix))
end

function _optimization_check_one_based(matrix::AbstractMatrix, name::AbstractString)
    firstindex(matrix, 1) == 1 && firstindex(matrix, 2) == 1 ||
        throw(ArgumentError("$name must use one-based indexing"))
    return nothing
end

function _optimization_check_hermitian(matrix::AbstractMatrix, name::AbstractString)
    _optimization_check_one_based(matrix, name)
    size(matrix, 1) == size(matrix, 2) || throw(DimensionMismatch("$name must be square"))
    for value in matrix
        value isa Number && !(value isa Bool) ||
            throw(ArgumentError("$name must contain numeric non-Boolean values"))
        isfinite(value) || throw(ArgumentError("$name must contain only finite values"))
    end
    ishermitian(matrix) ||
        throw(ArgumentError("$name must be exactly Hermitian; no repair is applied"))
    return nothing
end

"""
    HermitianAffineMatrix(name, constant, variables, coefficients, variable_count)

Sparse solver-neutral affine Hermitian matrix
`constant + sum(coefficients[k] * x[variables[k]])`.
"""
struct HermitianAffineMatrix{T<:Real}
    name::Symbol
    dimension::Int
    variable_count::Int
    constant::SparseMatrixCSC{Complex{T},Int}
    terms::Vector{HermitianAffineTerm{T}}
end

function HermitianAffineMatrix(
    name::Symbol,
    constant::AbstractMatrix{<:Number},
    variables::AbstractVector{<:Integer},
    coefficients::AbstractVector{<:AbstractMatrix{<:Number}},
    variable_count::Integer,
)
    count = _optimization_positive_int(variable_count, "variable_count")
    _optimization_check_hermitian(constant, "constant")
    firstindex(variables) == 1 && firstindex(coefficients) == 1 ||
        throw(ArgumentError("affine matrix inputs must use one-based indexing"))
    length(variables) == length(coefficients) ||
        throw(DimensionMismatch("variables and coefficients must have equal length"))
    dimension = size(constant, 1)
    real_types = Type[_optimization_matrix_real_type(constant)]
    checked_variables = Int[]
    for (position, (variable, coefficient)) in enumerate(zip(variables, coefficients))
        variable isa Bool &&
            throw(ArgumentError("variable indices must be one-based integers"))
        1 <= variable <= count ||
            throw(ArgumentError("variable index $variable is outside 1:$count"))
        _optimization_check_hermitian(coefficient, "coefficient $position")
        size(coefficient) == (dimension, dimension) || throw(
            DimensionMismatch(
                "coefficient $position has size $(size(coefficient)); expected " *
                "($dimension, $dimension)",
            ),
        )
        push!(real_types, _optimization_matrix_real_type(coefficient))
        push!(checked_variables, Int(variable))
    end
    length(unique(checked_variables)) == length(checked_variables) ||
        throw(ArgumentError("variables must not contain duplicates"))
    T = promote_type(real_types...)
    T <: Real ||
        throw(ArgumentError("Hermitian matrix entries must have real component type"))
    converted_constant = sparse(Complex{T}.(constant))
    terms = HermitianAffineTerm{T}[]
    for (variable, coefficient) in zip(checked_variables, coefficients)
        converted = sparse(Complex{T}.(coefficient))
        nnz(converted) == 0 && continue
        push!(terms, HermitianAffineTerm{T}(variable, converted))
    end
    return HermitianAffineMatrix{T}(name, dimension, count, converted_constant, terms)
end

"""
    SemidefiniteProgram

Package-owned affine SDP over real scalar coordinates. Each PSD constraint is a
possibly complex Hermitian affine matrix. The optional JuMP extension uses the
reviewed real block embedding; no backend expressions are stored here.
"""
struct SemidefiniteProgram{T<:Real,M}
    name::Symbol
    sense::Symbol
    variable_count::Int
    objective::AffineScalar{T}
    equalities::Vector{AffineEquality{T}}
    intervals::Vector{AffineInterval{T}}
    psd_constraints::Vector{HermitianAffineMatrix{T}}
    primal_views::Vector{HermitianAffineMatrix{T}}
    initial_point::Union{Nothing,Vector{T}}
    known_feasible_point::Union{Nothing,Vector{T}}
    limits::OptimizationLimits
    metadata::M
end

function _optimization_affine_convert(function_data::AffineScalar, ::Type{T}) where {T}
    return AffineScalar{T}(
        convert(T, function_data.constant),
        SparseVector{T,Int}(
            length(function_data.coefficients),
            copy(function_data.coefficients.nzind),
            T.(function_data.coefficients.nzval),
        ),
    )
end

function _optimization_matrix_convert(
    matrix::HermitianAffineMatrix, ::Type{T}
) where {T<:Real}
    terms = HermitianAffineTerm{T}[
        HermitianAffineTerm{T}(term.variable, sparse(Complex{T}.(term.coefficient))) for
        term in matrix.terms
    ]
    return HermitianAffineMatrix{T}(
        matrix.name,
        matrix.dimension,
        matrix.variable_count,
        sparse(Complex{T}.(matrix.constant)),
        terms,
    )
end

function _optimization_checked_point(point, count::Int, ::Type{T}, name) where {T}
    point === nothing && return nothing
    point isa AbstractVector ||
        throw(ArgumentError("$name must be nothing or an AbstractVector"))
    firstindex(point) == 1 || throw(ArgumentError("$name must use one-based indexing"))
    length(point) == count || throw(DimensionMismatch("$name must have length $count"))
    converted = Vector{T}(undef, count)
    for index in 1:count
        _optimization_validate_real(point[index], name)
        converted[index] = point[index]
    end
    return converted
end

function _optimization_preflight(
    variable_count::Int,
    objective::AffineScalar,
    equalities,
    intervals,
    psd_constraints,
    views,
    limits::OptimizationLimits,
)
    variable_count <= limits.max_variables || throw(
        ArgumentError(
            "model needs $variable_count variables, exceeding max_variables=$(limits.max_variables)",
        ),
    )
    length(equalities) <= limits.max_equalities || throw(
        ArgumentError(
            "model needs $(length(equalities)) equalities, exceeding " *
            "max_equalities=$(limits.max_equalities)",
        ),
    )
    length(intervals) <= limits.max_intervals || throw(
        ArgumentError(
            "model needs $(length(intervals)) intervals, exceeding " *
            "max_intervals=$(limits.max_intervals)",
        ),
    )
    length(psd_constraints) <= limits.max_psd_blocks || throw(
        ArgumentError(
            "model needs $(length(psd_constraints)) PSD blocks, exceeding " *
            "max_psd_blocks=$(limits.max_psd_blocks)",
        ),
    )
    for matrix in psd_constraints
        matrix.dimension <= limits.max_psd_dimension || throw(
            ArgumentError(
                "PSD block $(matrix.name) has dimension $(matrix.dimension), " *
                "exceeding max_psd_dimension=$(limits.max_psd_dimension)",
            ),
        )
    end
    entries = BigInt(nnz(objective.coefficients))
    for constraint in equalities
        entries += nnz(constraint.function_data.coefficients)
    end
    for constraint in intervals
        entries += nnz(constraint.function_data.coefficients)
    end
    for matrix in Iterators.flatten((psd_constraints, views))
        entries += nnz(matrix.constant)
        for term in matrix.terms
            entries += nnz(term.coefficient)
        end
    end
    entries <= limits.max_model_entries || throw(
        ArgumentError(
            "model stores $entries coefficient entries, exceeding " *
            "max_model_entries=$(limits.max_model_entries)",
        ),
    )
    real_block_entries = sum(
        (BigInt(2) * matrix.dimension) * (BigInt(2) * matrix.dimension + 1) ÷ 2 for
        matrix in psd_constraints;
        init=BigInt(0),
    )
    real_block_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "real-block PSD scalarization needs $real_block_entries triangle " *
            "entries, exceeding max_model_entries=$(limits.max_model_entries)",
        ),
    )
    return (
        stored_entries=Int(entries), real_block_triangle_entries=Int(real_block_entries)
    )
end

function SemidefiniteProgram(
    name::Symbol,
    sense::Symbol,
    variable_count::Integer,
    objective::AffineScalar;
    equalities=AffineEquality[],
    intervals=AffineInterval[],
    psd_constraints=HermitianAffineMatrix[],
    primal_views=HermitianAffineMatrix[],
    initial_point=nothing,
    known_feasible_point=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    metadata=NamedTuple(),
)
    sense in (:minimize, :maximize, :feasibility) ||
        throw(ArgumentError("sense must be :minimize, :maximize, or :feasibility"))
    count = _optimization_positive_int(variable_count, "variable_count")
    length(objective.coefficients) == count || throw(
        DimensionMismatch(
            "objective has $(length(objective.coefficients)) coefficients; expected $count",
        ),
    )
    all(
        collection -> firstindex(collection) == 1,
        (equalities, intervals, psd_constraints, primal_views),
    ) || throw(ArgumentError("model collections must use one-based indexing"))
    functions = Any[objective]
    append!(functions, (constraint.function_data for constraint in equalities))
    append!(functions, (constraint.function_data for constraint in intervals))
    real_types = Type[typeof(function_data.constant) for function_data in functions]
    append!(
        real_types,
        (
            _optimization_complex_real_type(eltype(matrix.constant)) for
            matrix in Iterators.flatten((psd_constraints, primal_views))
        ),
    )
    T = promote_type(real_types...)
    T <: Real || throw(ArgumentError("model data must promote to a real type"))

    converted_objective = _optimization_affine_convert(objective, T)
    converted_equalities = AffineEquality{T}[]
    for constraint in equalities
        length(constraint.function_data.coefficients) == count || throw(
            DimensionMismatch("equality $(constraint.name) has the wrong variable count"),
        )
        push!(
            converted_equalities,
            AffineEquality{T}(
                _optimization_affine_convert(constraint.function_data, T), constraint.name
            ),
        )
    end
    converted_intervals = AffineInterval{T}[]
    for constraint in intervals
        length(constraint.function_data.coefficients) == count || throw(
            DimensionMismatch("interval $(constraint.name) has the wrong variable count"),
        )
        push!(
            converted_intervals,
            AffineInterval{T}(
                _optimization_affine_convert(constraint.function_data, T),
                isnothing(constraint.lower) ? nothing : convert(T, constraint.lower),
                isnothing(constraint.upper) ? nothing : convert(T, constraint.upper),
                constraint.name,
            ),
        )
    end
    converted_psd = HermitianAffineMatrix{T}[]
    for matrix in psd_constraints
        matrix.variable_count == count || throw(
            DimensionMismatch("PSD block $(matrix.name) has the wrong variable count")
        )
        push!(converted_psd, _optimization_matrix_convert(matrix, T))
    end
    converted_views = HermitianAffineMatrix{T}[]
    view_names = Symbol[]
    for matrix in primal_views
        matrix.variable_count == count || throw(
            DimensionMismatch("primal view $(matrix.name) has the wrong variable count")
        )
        matrix.name in view_names &&
            throw(ArgumentError("duplicate primal view name $(repr(matrix.name))"))
        push!(view_names, matrix.name)
        push!(converted_views, _optimization_matrix_convert(matrix, T))
    end
    preflight = _optimization_preflight(
        count,
        converted_objective,
        converted_equalities,
        converted_intervals,
        converted_psd,
        converted_views,
        limits,
    )
    checked_initial = _optimization_checked_point(initial_point, count, T, "initial_point")
    checked_feasible = _optimization_checked_point(
        known_feasible_point, count, T, "known_feasible_point"
    )
    full_metadata = if metadata isa NamedTuple
        merge(metadata, (preflight=preflight,))
    else
        throw(ArgumentError("metadata must be a NamedTuple"))
    end
    return SemidefiniteProgram{T,typeof(full_metadata)}(
        name,
        sense,
        count,
        converted_objective,
        converted_equalities,
        converted_intervals,
        converted_psd,
        converted_views,
        checked_initial,
        checked_feasible,
        limits,
        full_metadata,
    )
end

"""
    OptimizerMetadata

Stable solver/model provenance retained independently of JuMP and MOI types.
"""
struct OptimizerMetadata{O,L}
    modeling_layer::Symbol
    modeling_layer_version::Union{Nothing,VersionNumber}
    interface_version::Union{Nothing,VersionNumber}
    configured_optimizer_name::String
    reported_optimizer_name::String
    optimizer_version::Union{Nothing,VersionNumber}
    options::O
    limits::L
    coefficient_type::DataType
    complex_psd_embedding::Symbol
    raw_status::String
end

"""
    OptimizationPrimal

Owned primal coordinate vector plus named reconstructed Hermitian matrix views.
"""
struct OptimizationPrimal{T<:Real,V}
    coordinates::Vector{T}
    views::V
end

"""
    OptimizationDual

Owned equality, interval, and PSD dual values. `psd_real_blocks` retains the
solver-facing real dual matrices and `psd_hermitian_blocks` their package-owned
adjoint embedding.
"""
struct OptimizationDual{T<:Real,R,C}
    equalities::Vector{T}
    intervals::Vector{T}
    psd_real_blocks::R
    psd_hermitian_blocks::C
end

"""
    OptimizationResult

Status-rich package-owned optimization result. Missing backend data is
`nothing`, never a fabricated scalar. `certified` is reserved for a separately
validated mathematical certificate; a numerically consistent solver optimum is
not automatically marked certified.
"""
struct OptimizationResult{T<:Real,P,D,M}
    status::OptimizationStatus
    termination_status::Symbol
    primal_status::Symbol
    dual_status::Symbol
    objective_value::Union{Nothing,T}
    objective_bound::Union{Nothing,T}
    dual_objective_value::Union{Nothing,T}
    absolute_gap::Union{Nothing,T}
    relative_gap::Union{Nothing,T}
    primal_residual::Union{Nothing,T}
    dual_residual::Union{Nothing,T}
    iterations::Union{Nothing,Int}
    solve_time_seconds::Union{Nothing,Float64}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    primal::P
    dual::D
    optimizer::M
    warnings::Tuple{Vararg{String}}
    message::String

    function OptimizationResult{T,P,D,M}(
        status::OptimizationStatus,
        termination_status::Symbol,
        primal_status::Symbol,
        dual_status::Symbol,
        objective_value::Union{Nothing,T},
        objective_bound::Union{Nothing,T},
        dual_objective_value::Union{Nothing,T},
        absolute_gap::Union{Nothing,T},
        relative_gap::Union{Nothing,T},
        primal_residual::Union{Nothing,T},
        dual_residual::Union{Nothing,T},
        iterations::Union{Nothing,Int},
        solve_time_seconds::Union{Nothing,Float64},
        certified::Bool,
        certificate_kind::Union{Nothing,Symbol},
        primal::P,
        dual::D,
        optimizer::M,
        warnings::Tuple{Vararg{String}},
        message::AbstractString,
    ) where {T<:Real,P,D,M}
        certified &&
            certificate_kind === nothing &&
            throw(ArgumentError("a certified result must name its certificate_kind"))
        status in (
                OptimizationBackendUnavailable,
                OptimizationMalformedBackend,
                OptimizationUnsupported,
                OptimizationNumericalFailure,
                OptimizationLimit,
                OptimizationInconsistent,
                OptimizationUnknown,
            ) &&
            certified &&
            throw(ArgumentError("failed, limited, or unknown results cannot be certified"))
        iterations !== nothing &&
            iterations < 0 &&
            throw(ArgumentError("iterations must be nonnegative"))
        solve_time_seconds !== nothing &&
            (!isfinite(solve_time_seconds) || solve_time_seconds < 0) &&
            throw(ArgumentError("solve_time_seconds must be finite and nonnegative"))
        for (name, value) in (
            ("objective_value", objective_value),
            ("objective_bound", objective_bound),
            ("dual_objective_value", dual_objective_value),
            ("absolute_gap", absolute_gap),
            ("relative_gap", relative_gap),
            ("primal_residual", primal_residual),
            ("dual_residual", dual_residual),
        )
            value === nothing && continue
            isfinite(value) ||
                throw(ArgumentError("$name must be finite when it is available"))
            name in ("absolute_gap", "relative_gap", "primal_residual", "dual_residual") &&
                value < zero(value) &&
                throw(ArgumentError("$name must be nonnegative"))
        end
        return new{T,P,D,M}(
            status,
            termination_status,
            primal_status,
            dual_status,
            objective_value,
            objective_bound,
            dual_objective_value,
            absolute_gap,
            relative_gap,
            primal_residual,
            dual_residual,
            iterations,
            solve_time_seconds,
            certified,
            certificate_kind,
            primal,
            dual,
            optimizer,
            warnings,
            String(message),
        )
    end
end

function Base.show(io::IO, result::OptimizationResult)
    return print(
        io,
        "OptimizationResult(status=",
        result.status,
        ", termination_status=",
        result.termination_status,
        ", objective_value=",
        result.objective_value,
        ", certified=",
        result.certified,
        ")",
    )
end

function _optimization_empty_metadata(
    problem::SemidefiniteProgram; raw_status::AbstractString=""
)
    return OptimizerMetadata(
        :none,
        nothing,
        nothing,
        "none",
        "none",
        nothing,
        NamedTuple(),
        problem.limits,
        eltype(problem.objective.coefficients),
        :not_built,
        String(raw_status),
    )
end

function _optimization_failure_result(
    problem::SemidefiniteProgram{T},
    status::OptimizationStatus,
    termination_status::Symbol,
    message::AbstractString;
    primal_status::Symbol=:no_solution,
    dual_status::Symbol=:no_solution,
    optimizer=_optimization_empty_metadata(problem; raw_status=message),
    warnings=(),
) where {T<:Real}
    return OptimizationResult{T,Nothing,Nothing,typeof(optimizer)}(
        status,
        termination_status,
        primal_status,
        dual_status,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        nothing,
        nothing,
        nothing,
        optimizer,
        Tuple(String(warning) for warning in warnings),
        message,
    )
end

"""
    solve_optimization(problem, backend=NoOptimizationBackend())

Solve a package-owned problem with an explicit backend. Without the optional
JuMP extension, this returns a structured backend-unavailable result rather
than throwing or selecting a global solver.
"""
function solve_optimization(
    problem::SemidefiniteProgram, ::NoOptimizationBackend=NoOptimizationBackend()
)
    return _optimization_failure_result(
        problem,
        OptimizationBackendUnavailable,
        :backend_unavailable,
        "no optimization backend was supplied; load JuMP and pass an explicit JuMPBackend",
    )
end

function _solve_optional_optimization(
    problem::SemidefiniteProgram, ::AbstractOptimizationBackend
)
    return _optimization_failure_result(
        problem,
        OptimizationBackendUnavailable,
        :extension_unavailable,
        "the JuMP package extension is not loaded; load JuMP before solving",
    )
end

function solve_optimization(problem::SemidefiniteProgram, backend::JuMPBackend)
    return _solve_optional_optimization(problem, backend)
end
