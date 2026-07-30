# Source-informed independent Julia implementations based on the specifications
# in QETLAB IskIncoherent.m, IsAbskIncoh.m, RobustnessCoherence.m,
# TraceDistanceCoherence.m, GenRobustnesskCoherence.m, and
# helpers/has_band_k_ordering.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
#
# Pinned source SHA-256 values:
# - IskIncoherent.m:
#   d924a12b239a1ea9a784ed3252ff0f289f2281aeae6022b645b169c4d2779f6a
# - IsAbskIncoh.m:
#   6bf68ba8b87e441dd6b2ab2b718f621e2fa6a8a67a92f016b3629266496658ef
# - RobustnessCoherence.m:
#   b8af3c8f056ba49ad703e361b1b6f8554776747b0790efefed34a359a28b7986
# - TraceDistanceCoherence.m:
#   f74de030eb833da8c1bd2352961e1080c652746fa11572812000f50a0cd2048b
# - GenRobustnesskCoherence.m:
#   aa3d429fc0fb9f715391620de67d9407f638bd58307922e49bf7ab5b57e37e78
# - helpers/has_band_k_ordering.m:
#   d4b6b6be5793515b3b98f98ba585f04041e6e7875f22437a36b413f6308be8d7
#
# QETLAB is BSD-2-Clause licensed; see licenses/QETLAB-LICENSE.txt. The native
# formulations below were reconstructed from the cited mathematical sources
# and use only package-owned solver-neutral optimization contracts.

"""
    CoherenceCriterionResult

Status-rich result for [`is_k_incoherent`](@ref) and
[`is_absolutely_k_incoherent`](@ref).

`verdict` is `true` or `false` only when an analytic theorem or a completed
feasibility solve supports that conclusion. It is `nothing` at numerical
boundaries, after resource or solver limits, and when no implemented theorem
decides the property. `certificate_kind` and `exact` distinguish mathematical
theorem certificates from tolerance-dependent solver evidence.
"""
struct CoherenceCriterionResult{T<:AbstractFloat,S,D,O,P,R,Q}
    property::Symbol
    status::Symbol
    verdict::Union{Nothing,Bool}
    k::Int
    dimension::Int
    method::Symbol
    certificate_kind::Union{Nothing,Symbol}
    exact::Bool
    margin::Union{Nothing,T}
    purity::T
    spectrum::S
    decomposition::D
    band_ordering::O
    problem::P
    optimization::R
    diagnostics::Q
    message::String
end

function Base.show(io::IO, result::CoherenceCriterionResult)
    return print(
        io,
        "CoherenceCriterionResult(property=",
        result.property,
        ", status=",
        result.status,
        ", verdict=",
        result.verdict,
        ", k=",
        result.k,
        ", method=",
        result.method,
        ")",
    )
end

"""
    CoherenceOptimizationResult

Package-owned result for coherence optimization measures.

`value` is absent unless an analytic formula or a usable optimizer primal is
available. `free_state`, `unnormalized_noise`, `noise_state`, and
`decomposition` retain the optimizer/certificate data when the formulation
defines them. For zero robustness, `noise_state === nothing`; the zero matrix
is retained in `unnormalized_noise`, so no division by zero is attempted.
"""
struct CoherenceOptimizationResult{T<:AbstractFloat,F,N,U,D,A,B,P,R,Q}
    quantity::Symbol
    status::Symbol
    value::Union{Nothing,T}
    k::Union{Nothing,Int}
    dimension::Int
    method::Symbol
    exact::Bool
    free_state::F
    noise_state::N
    unnormalized_noise::U
    decomposition::D
    positive_part::A
    negative_part::B
    problem::P
    optimization::R
    diagnostics::Q
    message::String
end

function Base.show(io::IO, result::CoherenceOptimizationResult)
    return print(
        io,
        "CoherenceOptimizationResult(quantity=",
        result.quantity,
        ", status=",
        result.status,
        ", value=",
        result.value,
        ", method=",
        result.method,
        ")",
    )
end

function _cohopt_real_type(::Type{T}) where {T<:Number}
    R = typeof(abs(zero(T)))
    R <: AbstractFloat || throw(
        ArgumentError(
            "coherence optimization requires real or complex floating-point data; " *
            "got element type $T. No precision-changing conversion is applied.",
        ),
    )
    return R
end

function _cohopt_tolerance(value, default::R, name::AbstractString) where {R<:AbstractFloat}
    value === nothing && return default
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$name must be a finite nonnegative real or `nothing`"))
    isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("$name must be a finite nonnegative real or `nothing`"))
    converted = try
        convert(R, value)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("$name=$(repr(value)) cannot be represented as $R"))
    end
    isfinite(converted) ||
        throw(ArgumentError("$name=$(repr(value)) overflows the input precision $R"))
    return converted
end

function _cohopt_state(
    state::AbstractVector{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    Base.require_one_based_indexing(state)
    isempty(state) && throw(ArgumentError("state must have positive length"))
    issparse(state) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "validating a sparse pure state and constructing its density matrix " *
                "would densify; pass allow_densify=true explicitly",
            ),
        )
    all(isfinite, state) || throw(ArgumentError("state must contain only finite entries"))
    R = _cohopt_real_type(eltype(state))
    absolute = _cohopt_tolerance(atol, zero(R), "atol")
    relative = _cohopt_tolerance(rtol, sqrt(eps(R)), "rtol")
    vector = collect(state)
    norm_squared = sum(abs2, vector)
    normalization_tolerance = absolute + relative * max(one(R), abs(norm_squared))
    abs(norm_squared - one(R)) <= normalization_tolerance || throw(
        ArgumentError(
            "state is not normalized within atol=$absolute and rtol=$relative; " *
            "its squared norm is $norm_squared. The input is never normalized.",
        ),
    )
    rho = vector * adjoint(vector)
    spectrum = vcat(fill(zero(R), length(vector) - 1), R(norm_squared))
    return (
        rho=rho,
        pure_vector=vector,
        dimension=length(vector),
        real_type=R,
        atol=absolute,
        rtol=relative,
        trace_value=R(norm_squared),
        trace_residual=abs(R(norm_squared) - one(R)),
        trace_tolerance=normalization_tolerance,
        spectrum=spectrum,
        minimum_eigenvalue=zero(R),
        psd_tolerance=absolute + relative,
        rank_lower=1,
        rank_upper=1,
        rank_tolerance=absolute + relative,
        input_kind=:pure_vector,
        densified=issparse(state),
    )
end

function _cohopt_state(
    state::AbstractMatrix{<:Number}; atol=nothing, rtol=nothing, allow_densify::Bool=false
)
    Base.require_one_based_indexing(state)
    size(state, 1) == size(state, 2) ||
        throw(DimensionMismatch("density matrix must be square"))
    size(state, 1) > 0 ||
        throw(ArgumentError("density matrix must have positive dimension"))
    issparse(state) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "density validation requires a dense Hermitian eigensolve; " *
                "pass allow_densify=true explicitly for sparse input",
            ),
        )
    all(isfinite, state) ||
        throw(ArgumentError("density matrix must contain only finite entries"))
    ishermitian(state) || throw(
        ArgumentError(
            "density matrix must be exactly Hermitian; no symmetrization is applied"
        ),
    )
    R = _cohopt_real_type(eltype(state))
    absolute = _cohopt_tolerance(atol, zero(R), "atol")
    relative = _cohopt_tolerance(rtol, sqrt(eps(R)), "rtol")
    rho = Matrix(state)
    trace_value = real(tr(rho))
    trace_tolerance = absolute + relative * max(one(R), abs(trace_value))
    abs(trace_value - one(R)) <= trace_tolerance || throw(
        ArgumentError(
            "density matrix trace is $trace_value, outside atol=$absolute and " *
            "rtol=$relative of one. The input is never normalized.",
        ),
    )
    eig = eigen(Hermitian(rho))
    spectrum = R.(real.(eig.values))
    spectral_scale = max(one(R), maximum(abs, spectrum; init=zero(R)))
    psd_tolerance = absolute + relative * spectral_scale
    minimum_eigenvalue = minimum(spectrum)
    minimum_eigenvalue >= -psd_tolerance || throw(
        ArgumentError(
            "density matrix is not positive semidefinite within the requested " *
            "tolerance; minimum eigenvalue is $minimum_eigenvalue. No clipping " *
            "or projection is applied.",
        ),
    )
    rank_tolerance = absolute + relative * spectral_scale
    boundary_rank_tolerance = 8rank_tolerance
    rank_lower = count(value -> value > boundary_rank_tolerance, spectrum)
    rank_upper = count(value -> value > rank_tolerance, spectrum)
    pure_vector = if rank_lower == rank_upper == 1
        index = argmax(spectrum)
        sqrt(spectrum[index]) .* eig.vectors[:, index]
    else
        nothing
    end
    return (
        rho=rho,
        pure_vector=pure_vector,
        dimension=size(rho, 1),
        real_type=R,
        atol=absolute,
        rtol=relative,
        trace_value=R(trace_value),
        trace_residual=abs(R(trace_value) - one(R)),
        trace_tolerance=R(trace_tolerance),
        spectrum=spectrum,
        minimum_eigenvalue=R(minimum_eigenvalue),
        psd_tolerance=R(psd_tolerance),
        rank_lower,
        rank_upper,
        rank_tolerance=R(rank_tolerance),
        input_kind=:density_matrix,
        densified=issparse(state),
    )
end

function _cohopt_k(k, dimension::Int)
    k isa Integer && !(k isa Bool) ||
        throw(ArgumentError("k must be an integer in 1:$dimension"))
    1 <= k <= dimension || throw(DomainError(k, "k must lie in 1:$dimension"))
    return Int(k)
end

function _cohopt_positive_limit(value, name::AbstractString)
    value isa Integer && !(value isa Bool) && value > 0 ||
        throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _cohopt_strategy(strategy)
    strategy isa Symbol && strategy in (:auto, :sdp) ||
        throw(ArgumentError("strategy must be :auto or :sdp"))
    return strategy
end

function _cohopt_state_diagnostics(data)
    return (
        input_kind=data.input_kind,
        trace_value=data.trace_value,
        trace_residual=data.trace_residual,
        trace_tolerance=data.trace_tolerance,
        minimum_eigenvalue=data.minimum_eigenvalue,
        psd_tolerance=data.psd_tolerance,
        rank_lower=data.rank_lower,
        rank_upper=data.rank_upper,
        rank_tolerance=data.rank_tolerance,
        densified=data.densified,
    )
end

function _cohopt_criterion(
    data,
    property::Symbol,
    status::Symbol,
    verdict,
    k::Int,
    method::Symbol,
    certificate_kind,
    exact::Bool,
    message::AbstractString;
    margin=nothing,
    decomposition=nothing,
    band_ordering=nothing,
    problem=nothing,
    optimization=nothing,
    diagnostics=NamedTuple(),
)
    T = data.real_type
    converted_margin = margin === nothing ? nothing : T(margin)
    full_diagnostics = merge(_cohopt_state_diagnostics(data), diagnostics)
    return CoherenceCriterionResult(
        property,
        status,
        verdict,
        k,
        data.dimension,
        method,
        certificate_kind,
        exact,
        converted_margin,
        T(real(tr(data.rho * data.rho))),
        copy(data.spectrum),
        decomposition,
        band_ordering,
        problem,
        optimization,
        full_diagnostics,
        String(message),
    )
end

function _cohopt_measure(
    data,
    quantity::Symbol,
    status::Symbol,
    value,
    method::Symbol,
    exact::Bool,
    message::AbstractString;
    k=nothing,
    free_state=nothing,
    noise_state=nothing,
    unnormalized_noise=nothing,
    decomposition=nothing,
    positive_part=nothing,
    negative_part=nothing,
    problem=nothing,
    optimization=nothing,
    diagnostics=NamedTuple(),
)
    T = data.real_type
    converted_value = value === nothing ? nothing : T(value)
    converted_value === nothing ||
        (isfinite(converted_value) && converted_value >= zero(T)) ||
        throw(ArgumentError("coherence measure values must be finite and nonnegative"))
    full_diagnostics = merge(_cohopt_state_diagnostics(data), diagnostics)
    return CoherenceOptimizationResult{
        T,
        typeof(free_state),
        typeof(noise_state),
        typeof(unnormalized_noise),
        typeof(decomposition),
        typeof(positive_part),
        typeof(negative_part),
        typeof(problem),
        typeof(optimization),
        typeof(full_diagnostics),
    }(
        quantity,
        status,
        converted_value,
        k,
        data.dimension,
        method,
        exact,
        free_state,
        noise_state,
        unnormalized_noise,
        decomposition,
        positive_part,
        negative_part,
        problem,
        optimization,
        full_diagnostics,
        String(message),
    )
end

function _cohopt_rebuild(
    name::Symbol,
    matrix::HermitianAffineMatrix;
    constant=matrix.constant,
    coefficient_scale=one(real(eltype(matrix.constant))),
)
    return HermitianAffineMatrix(
        name,
        constant,
        [term.variable for term in matrix.terms],
        [coefficient_scale .* term.coefficient for term in matrix.terms],
        matrix.variable_count,
    )
end

function _cohopt_affine_sum(
    name::Symbol, matrices::AbstractVector{<:HermitianAffineMatrix}; constant=nothing
)
    isempty(matrices) && throw(ArgumentError("at least one affine matrix is required"))
    dimension = first(matrices).dimension
    variable_count = first(matrices).variable_count
    all(
        matrix -> matrix.dimension == dimension && matrix.variable_count == variable_count,
        matrices,
    ) || throw(DimensionMismatch("affine matrices must have matching dimensions"))
    T = promote_type(
        (typeof(real(zero(eltype(matrix.constant)))) for matrix in matrices)...
    )
    base = if constant === nothing
        spzeros(Complex{T}, dimension, dimension)
    else
        sparse(Complex{T}.(constant))
    end
    for matrix in matrices
        base += matrix.constant
    end
    coefficient_map = Dict{Int,SparseMatrixCSC{Complex{T},Int}}()
    for matrix in matrices, term in matrix.terms
        coefficient = sparse(Complex{T}.(term.coefficient))
        coefficient_map[term.variable] =
            get(coefficient_map, term.variable, spzeros(Complex{T}, dimension, dimension)) +
            coefficient
    end
    variables = sort!(collect(keys(coefficient_map)))
    coefficients = [coefficient_map[index] for index in variables]
    return HermitianAffineMatrix(name, base, variables, coefficients, variable_count)
end

function _cohopt_embed(
    matrix::HermitianAffineMatrix,
    support::AbstractVector{<:Integer},
    dimension::Int,
    name::Symbol,
)
    length(support) == matrix.dimension ||
        throw(DimensionMismatch("support length must match block dimension"))
    support_int = Int.(support)
    length(unique(support_int)) == length(support_int) ||
        throw(ArgumentError("support indices must be distinct"))
    all(index -> 1 <= index <= dimension, support_int) ||
        throw(ArgumentError("support index lies outside 1:$dimension"))
    T = typeof(real(zero(eltype(matrix.constant))))
    embed_one = function (local_matrix)
        rows, columns, values = findnz(sparse(local_matrix))
        return sparse(
            support_int[rows],
            support_int[columns],
            Complex{T}.(values),
            dimension,
            dimension,
        )
    end
    return HermitianAffineMatrix(
        name,
        embed_one(matrix.constant),
        [term.variable for term in matrix.terms],
        [embed_one(term.coefficient) for term in matrix.terms],
        matrix.variable_count,
    )
end

function _cohopt_diagonal_affine(
    name::Symbol, dimension::Int, first_variable::Int, variable_count::Int, ::Type{T}
) where {T<:Real}
    variables = collect(first_variable:(first_variable + dimension - 1))
    coefficients = AbstractMatrix{<:Number}[
        sparse([index], [index], Complex{T}[one(T)], dimension, dimension) for
        index in 1:dimension
    ]
    return HermitianAffineMatrix(
        name,
        spzeros(Complex{T}, dimension, dimension),
        variables,
        coefficients,
        variable_count,
    )
end

function _cohopt_scalar_sum(functions::AbstractVector{<:AffineScalar}, count::Int)
    isempty(functions) && return AffineScalar(0.0, zeros(count))
    T = promote_type((typeof(function_data.constant) for function_data in functions)...)
    constant = zero(T)
    coefficients = zeros(T, count)
    for function_data in functions
        constant += function_data.constant
        indices, values = findnz(function_data.coefficients)
        coefficients[indices] .+= values
    end
    return AffineScalar(constant, coefficients)
end

function _cohopt_nonnegative_intervals(
    first_variable::Int, count::Int, variable_count::Int, ::Type{T}, prefix::Symbol
) where {T<:Real}
    constraints = AffineInterval{T}[]
    for offset in 0:(count - 1)
        index = first_variable + offset
        push!(
            constraints,
            AffineInterval(
                AffineScalar(zero(T), [index], T[one(T)], variable_count),
                zero(T),
                nothing,
                Symbol(prefix, :_, offset + 1),
            ),
        )
    end
    return constraints
end

function _cohopt_combinations(dimension::Int, k::Int, max_subsets::Int)
    total = binomial(BigInt(dimension), BigInt(k))
    total <= max_subsets || return nothing, total
    supports = Vector{Vector{Int}}()
    current = collect(1:k)
    while true
        push!(supports, copy(current))
        position = k
        while position >= 1 && current[position] == dimension - k + position
            position -= 1
        end
        position == 0 && break
        current[position] += 1
        for index in (position + 1):k
            current[index] = current[index - 1] + 1
        end
    end
    return supports, total
end

function _cohopt_try_model(builder)
    try
        return builder(), nothing
    catch error
        if error isa ArgumentError && occursin("exceeding max_", sprint(showerror, error))
            return nothing, sprint(showerror, error)
        end
        rethrow()
    end
end

function _cohopt_structural_limit_message(
    variable_count::BigInt,
    equality_count::BigInt,
    interval_count::BigInt,
    psd_families,
    limits::OptimizationLimits,
)
    variable_count <= limits.max_variables ||
        return "model needs $variable_count variables, exceeding " *
               "max_variables=$(limits.max_variables)"
    equality_count <= limits.max_equalities ||
        return "model needs $equality_count equalities, exceeding " *
               "max_equalities=$(limits.max_equalities)"
    interval_count <= limits.max_intervals ||
        return "model needs $interval_count intervals, exceeding " *
               "max_intervals=$(limits.max_intervals)"
    psd_block_count = sum(BigInt(count) for (_, count) in psd_families; init=BigInt(0))
    psd_block_count <= limits.max_psd_blocks ||
        return "model needs $psd_block_count PSD blocks, exceeding " *
               "max_psd_blocks=$(limits.max_psd_blocks)"
    for (dimension, count) in psd_families
        iszero(count) && continue
        dimension <= limits.max_psd_dimension ||
            return "PSD blocks have dimension $dimension, exceeding " *
                   "max_psd_dimension=$(limits.max_psd_dimension)"
    end
    real_block_entries = sum(
        BigInt(count) * (2BigInt(dimension)) * (2BigInt(dimension) + 1) ÷ 2 for
        (dimension, count) in psd_families;
        init=BigInt(0),
    )
    real_block_entries <= limits.max_model_entries ||
        return "real-block PSD scalarization needs $real_block_entries triangle " *
               "entries, exceeding max_model_entries=$(limits.max_model_entries)"
    return nothing
end

function _cohopt_k_incoherence_model(
    rho::AbstractMatrix{<:Number},
    k::Int,
    ::Type{T},
    limits::OptimizationLimits,
    max_subsets::Int,
) where {T<:Real}
    dimension = size(rho, 1)
    supports, total = _cohopt_combinations(dimension, k, max_subsets)
    supports === nothing &&
        return nothing, "model needs $total k-subsets, exceeding max_subsets=$max_subsets"
    variable_count_big = total * BigInt(k)^2
    variable_count_big <= typemax(Int) ||
        return nothing, "model variable count $variable_count_big exceeds Int"
    structural_limit = _cohopt_structural_limit_message(
        variable_count_big, BigInt(dimension)^2, BigInt(0), ((k, total),), limits
    )
    structural_limit === nothing || return nothing, structural_limit
    variable_count = Int(variable_count_big)
    local_blocks = HermitianAffineMatrix[]
    embedded_blocks = HermitianAffineMatrix[]
    next_variable = 1
    for (index, support) in enumerate(supports)
        local_block = hermitian_variable(
            Symbol(:k_block_, index),
            k;
            first_variable=next_variable,
            variable_count,
            coefficient_type=T,
        )
        push!(local_blocks, local_block)
        push!(
            embedded_blocks,
            _cohopt_embed(
                local_block, support, dimension, Symbol(:embedded_k_block_, index)
            ),
        )
        next_variable += k^2
    end
    reconstructed = _cohopt_affine_sum(:k_incoherent_reconstruction, embedded_blocks)
    equalities = hermitian_equalities(
        reconstructed; target=rho, name_prefix=:k_incoherent_balance
    )
    objective = AffineScalar(zero(T), zeros(T, variable_count))
    views = vcat(HermitianAffineMatrix[reconstructed], local_blocks)
    problem, limit_message = _cohopt_try_model() do
        return SemidefiniteProgram(
            :k_incoherence_feasibility,
            :feasibility,
            variable_count,
            objective;
            equalities,
            psd_constraints=local_blocks,
            primal_views=views,
            limits,
            metadata=(
                formulation=:factor_width_block_decomposition,
                dimension,
                k,
                supports=Tuple(Tuple(support) for support in supports),
            ),
        )
    end
    return problem, limit_message
end

function _cohopt_generalized_model(
    rho::AbstractMatrix{<:Number},
    k::Int,
    ::Type{T},
    limits::OptimizationLimits,
    max_subsets::Int,
) where {T<:Real}
    dimension = size(rho, 1)
    supports, total = _cohopt_combinations(dimension, k, max_subsets)
    supports === nothing &&
        return nothing, "model needs $total k-subsets, exceeding max_subsets=$max_subsets"
    variable_count_big = BigInt(dimension)^2 + total * BigInt(k)^2
    variable_count_big <= typemax(Int) ||
        return nothing, "model variable count $variable_count_big exceeds Int"
    structural_limit = _cohopt_structural_limit_message(
        variable_count_big,
        BigInt(dimension)^2,
        BigInt(0),
        ((dimension, BigInt(1)), (k, total)),
        limits,
    )
    structural_limit === nothing || return nothing, structural_limit
    variable_count = Int(variable_count_big)
    noise = hermitian_variable(
        :generalized_coherence_noise,
        dimension;
        first_variable=1,
        variable_count,
        coefficient_type=T,
    )
    local_blocks = HermitianAffineMatrix[]
    embedded_blocks = HermitianAffineMatrix[]
    next_variable = dimension^2 + 1
    for (index, support) in enumerate(supports)
        local_block = hermitian_variable(
            Symbol(:generalized_k_block_, index),
            k;
            first_variable=next_variable,
            variable_count,
            coefficient_type=T,
        )
        push!(local_blocks, local_block)
        push!(
            embedded_blocks,
            _cohopt_embed(
                local_block,
                support,
                dimension,
                Symbol(:embedded_generalized_k_block_, index),
            ),
        )
        next_variable += k^2
    end
    free_unnormalized = _cohopt_affine_sum(:generalized_free_unnormalized, embedded_blocks)
    negative_noise = _cohopt_rebuild(
        :negative_generalized_noise, noise; coefficient_scale=(-one(T))
    )
    balance = _cohopt_affine_sum(
        :generalized_coherence_balance,
        HermitianAffineMatrix[free_unnormalized, negative_noise],
    )
    equalities = hermitian_equalities(
        balance; target=rho, name_prefix=:generalized_coherence_balance
    )
    objective = trace_affine(noise)
    views = vcat(HermitianAffineMatrix[noise, free_unnormalized], local_blocks)
    problem, limit_message = _cohopt_try_model() do
        return SemidefiniteProgram(
            :generalized_k_coherence_robustness,
            :minimize,
            variable_count,
            objective;
            equalities,
            psd_constraints=vcat(HermitianAffineMatrix[noise], local_blocks),
            primal_views=views,
            limits,
            metadata=(
                formulation=:generalized_factor_width_robustness,
                dimension,
                k,
                supports=Tuple(Tuple(support) for support in supports),
            ),
        )
    end
    return problem, limit_message
end

function _cohopt_robustness_model(
    rho::AbstractMatrix{<:Number}, ::Type{T}, limits::OptimizationLimits
) where {T<:Real}
    dimension = size(rho, 1)
    variable_count_big = BigInt(dimension)^2 + dimension
    structural_limit = _cohopt_structural_limit_message(
        variable_count_big,
        BigInt(dimension)^2,
        BigInt(dimension),
        ((dimension, BigInt(1)),),
        limits,
    )
    structural_limit === nothing || return nothing, structural_limit
    variable_count = Int(variable_count_big)
    noise = hermitian_variable(
        :coherence_noise, dimension; first_variable=1, variable_count, coefficient_type=T
    )
    free_unnormalized = _cohopt_diagonal_affine(
        :incoherent_free_unnormalized, dimension, dimension^2 + 1, variable_count, T
    )
    negative_free = _cohopt_rebuild(
        :negative_incoherent_free, free_unnormalized; coefficient_scale=(-one(T))
    )
    balance = _cohopt_affine_sum(
        :coherence_robustness_balance,
        HermitianAffineMatrix[noise, negative_free];
        constant=rho,
    )
    equalities = hermitian_equalities(balance; name_prefix=:coherence_robustness_balance)
    intervals = _cohopt_nonnegative_intervals(
        dimension^2 + 1, dimension, variable_count, T, :free_diagonal
    )
    problem, limit_message = _cohopt_try_model() do
        return SemidefiniteProgram(
            :robustness_of_coherence,
            :minimize,
            variable_count,
            trace_affine(noise);
            equalities,
            intervals,
            psd_constraints=[noise],
            primal_views=[noise, free_unnormalized],
            limits,
            metadata=(formulation=:incoherent_cone_robustness, dimension),
        )
    end
    return problem, limit_message
end

function _cohopt_trace_model(
    rho::AbstractMatrix{<:Number}, ::Type{T}, limits::OptimizationLimits
) where {T<:Real}
    dimension = size(rho, 1)
    variable_count_big = 2BigInt(dimension)^2 + dimension
    structural_limit = _cohopt_structural_limit_message(
        variable_count_big,
        BigInt(dimension)^2 + 1,
        BigInt(dimension),
        ((dimension, BigInt(2)),),
        limits,
    )
    structural_limit === nothing || return nothing, structural_limit
    variable_count = Int(variable_count_big)
    positive = hermitian_variable(
        :trace_positive_part,
        dimension;
        first_variable=1,
        variable_count,
        coefficient_type=T,
    )
    negative = hermitian_variable(
        :trace_negative_part,
        dimension;
        first_variable=dimension^2 + 1,
        variable_count,
        coefficient_type=T,
    )
    free = _cohopt_diagonal_affine(
        :trace_closest_incoherent, dimension, 2dimension^2 + 1, variable_count, T
    )
    negative_negative = _cohopt_rebuild(
        :negative_trace_negative_part, negative; coefficient_scale=(-one(T))
    )
    balance = _cohopt_affine_sum(
        :trace_distance_balance,
        HermitianAffineMatrix[positive, negative_negative, free];
        constant=(-rho),
    )
    equalities = hermitian_equalities(balance; name_prefix=:trace_distance_balance)
    intervals = _cohopt_nonnegative_intervals(
        2dimension^2 + 1, dimension, variable_count, T, :free_probability
    )
    probability_coefficients = zeros(T, variable_count)
    probability_coefficients[(2dimension ^ 2 + 1):end] .= one(T)
    push!(
        equalities,
        AffineEquality(
            AffineScalar(-one(T), probability_coefficients), :free_probability_sum
        ),
    )
    objective = _cohopt_scalar_sum(
        [trace_affine(positive), trace_affine(negative)], variable_count
    )
    problem, limit_message = _cohopt_try_model() do
        return SemidefiniteProgram(
            :trace_distance_of_coherence,
            :minimize,
            variable_count,
            objective;
            equalities,
            intervals,
            psd_constraints=[positive, negative],
            primal_views=[positive, negative, free],
            limits,
            metadata=(formulation=:hermitian_trace_norm_epigraph, dimension),
        )
    end
    return problem, limit_message
end

function _cohopt_coordinate_equality(
    matrix::HermitianAffineMatrix{T},
    row::Int,
    column::Int,
    component::Symbol,
    target::T,
    name::Symbol,
) where {T<:Real}
    constant_entry = matrix.constant[row, column]
    constant = (component === :real ? real(constant_entry) : imag(constant_entry)) - target
    variables = Int[]
    coefficients = T[]
    for term in matrix.terms
        entry = term.coefficient[row, column]
        value = component === :real ? real(entry) : imag(entry)
        iszero(value) && continue
        push!(variables, term.variable)
        push!(coefficients, value)
    end
    return AffineEquality(
        AffineScalar(constant, variables, coefficients, matrix.variable_count), name
    )
end

function _cohopt_absolute_nminus1_model(
    spectrum::AbstractVector{T}, limits::OptimizationLimits
) where {T<:AbstractFloat}
    dimension = length(spectrum)
    variable_count_big = BigInt(dimension)^2
    equality_count = BigInt(dimension) * (dimension - 1) ÷ 2 + dimension
    structural_limit = _cohopt_structural_limit_message(
        variable_count_big, equality_count, BigInt(0), ((dimension, BigInt(1)),), limits
    )
    structural_limit === nothing || return nothing, structural_limit
    variable_count = Int(variable_count_big)
    matrix = hermitian_variable(
        :absolute_nminus1_matrix, dimension; variable_count, coefficient_type=T
    )
    equalities = AffineEquality{T}[]
    for column in 2:dimension, row in 1:(column - 1)
        push!(
            equalities,
            _cohopt_coordinate_equality(
                matrix,
                row,
                column,
                :imag,
                zero(T),
                Symbol(:absolute_real_, row, :_, column),
            ),
        )
    end
    for index in 2:dimension
        push!(
            equalities,
            _cohopt_coordinate_equality(
                matrix,
                index,
                index,
                :real,
                spectrum[index],
                Symbol(:absolute_diagonal_, index),
            ),
        )
    end
    variables = Int[]
    coefficients = T[]
    constant = spectrum[1]
    for term in matrix.terms
        coefficient =
            real(term.coefficient[1, 1]) +
            2sum(real(term.coefficient[1, column]) for column in 2:dimension)
        iszero(coefficient) && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    push!(
        equalities,
        AffineEquality(
            AffineScalar(constant, variables, coefficients, variable_count),
            :absolute_first_row_balance,
        ),
    )
    objective = AffineScalar(zero(T), zeros(T, variable_count))
    problem, limit_message = _cohopt_try_model() do
        return SemidefiniteProgram(
            :absolute_nminus1_incoherence,
            :feasibility,
            variable_count,
            objective;
            equalities,
            psd_constraints=[matrix],
            primal_views=[matrix],
            limits,
            metadata=(
                formulation=:johnston_moein_pereira_plosker_theorem_8,
                dimension,
                sorted_spectrum=Tuple(spectrum),
            ),
        )
    end
    return problem, limit_message
end

function _cohopt_optimization_status(status::OptimizationStatus)
    status === OptimizationOptimal && return :optimal
    status === OptimizationFeasible && return :feasible
    status === OptimizationInfeasible && return :infeasible
    status === OptimizationUnbounded && return :unbounded
    status === OptimizationLimit && return :solver_limit
    status === OptimizationNumericalFailure && return :numerical_failure
    status === OptimizationUnsupported && return :unsupported
    status === OptimizationBackendUnavailable && return :backend_unavailable
    status === OptimizationMalformedBackend && return :malformed_backend
    status === OptimizationInconsistent && return :inconsistent
    return :unknown
end

function _cohopt_solver_criterion(
    data,
    property::Symbol,
    k::Int,
    problem::SemidefiniteProgram,
    optimization::OptimizationResult;
    decomposition_builder=nothing,
)
    solver_status = _cohopt_optimization_status(optimization.status)
    if optimization.status in (OptimizationOptimal, OptimizationFeasible)
        optimization.primal === nothing && return _cohopt_criterion(
            data,
            property,
            :malformed_backend,
            nothing,
            k,
            :sdp,
            nothing,
            false,
            "optimizer reported a usable status without a primal solution";
            problem,
            optimization,
        )
        decomposition = if decomposition_builder === nothing
            nothing
        else
            decomposition_builder(optimization.primal.views)
        end
        return _cohopt_criterion(
            data,
            property,
            :solver_feasible,
            true,
            k,
            :sdp,
            :numerical_primal_feasibility,
            false,
            "the solver returned a residual-checked feasible primal";
            decomposition,
            problem,
            optimization,
            diagnostics=(
                solver_status=solver_status, primal_residual=optimization.primal_residual
            ),
        )
    elseif optimization.status === OptimizationInfeasible
        return _cohopt_criterion(
            data,
            property,
            :solver_infeasible,
            false,
            k,
            :sdp,
            :numerical_infeasibility_status,
            false,
            "the solver reported the exact feasibility model infeasible";
            problem,
            optimization,
            diagnostics=(solver_status=solver_status,),
        )
    end
    return _cohopt_criterion(
        data,
        property,
        solver_status,
        nothing,
        k,
        :sdp,
        nothing,
        false,
        "the solver did not return a conclusive feasibility status";
        problem,
        optimization,
        diagnostics=(solver_status=solver_status,),
    )
end

function _cohopt_bandwidth(adjacency::AbstractMatrix{Bool}, ordering)
    inverse = zeros(Int, length(ordering))
    for (position, vertex) in enumerate(ordering)
        inverse[vertex] = position
    end
    maximum_distance = -1
    for column in axes(adjacency, 2), row in 1:(column - 1)
        adjacency[row, column] || continue
        maximum_distance = max(maximum_distance, abs(inverse[row] - inverse[column]))
    end
    return maximum_distance < 0 ? 0 : maximum_distance + 1
end

# The pinned MATLAB helper calls
# order_is_reversed(unselected_new, num_placed, candidate), although its nested
# function is declared as (unselected, candidate, num_placed). Rather than
# reproduce that argument-order defect or the resulting invalid pruning, this
# bounded recognizer independently enumerates graph layouts and checks exactly
# the package/QETLAB one-based bandwidth convention.
function _cohopt_band_ordering(
    matrix::AbstractMatrix{<:Number}, k::Int; max_search_nodes::Int=100_000
)
    max_search_nodes > 0 || throw(ArgumentError("max_band_search_nodes must be positive"))
    dimension = size(matrix, 1)
    adjacency = falses(dimension, dimension)
    for column in 2:dimension, row in 1:(column - 1)
        supported = !iszero(matrix[row, column]) || !iszero(matrix[column, row])
        adjacency[row, column] = supported
        adjacency[column, row] = supported
    end
    identity_order = collect(1:dimension)
    _cohopt_bandwidth(adjacency, identity_order) <= k && return (
        status=:found,
        ordering=identity_order,
        nodes=0,
        exact=true,
        convention=:one_based_bandwidth,
    )
    k == 0 && return (
        status=:not_found,
        ordering=Int[],
        nodes=0,
        exact=true,
        convention=:one_based_bandwidth,
    )
    ordering = zeros(Int, dimension)
    selected = falses(dimension)
    degrees = vec(sum(adjacency; dims=1))
    nodes = Ref(0)
    limited = Ref(false)
    function search(position::Int)
        position > dimension && return true
        candidates = sort(
            [vertex for vertex in 1:dimension if !selected[vertex]];
            by=vertex -> (-degrees[vertex], vertex),
        )
        for candidate in candidates
            nodes[] += 1
            if nodes[] > max_search_nodes
                limited[] = true
                return false
            end
            compatible = true
            last_far_position = position - k
            if last_far_position >= 1
                for prior_position in 1:last_far_position
                    if adjacency[ordering[prior_position], candidate]
                        compatible = false
                        break
                    end
                end
            end
            compatible || continue
            ordering[position] = candidate
            selected[candidate] = true
            search(position + 1) && return true
            selected[candidate] = false
            limited[] && return false
        end
        return false
    end
    found = search(1)
    return (
        status=if found
            :found
        elseif limited[]
            :limit
        else
            :not_found
        end,
        ordering=found ? copy(ordering) : Int[],
        nodes=nodes[],
        exact=(!limited[]),
        convention=:one_based_bandwidth,
    )
end

function _cohopt_block_decomposition(views, supports, prefix::Symbol)
    blocks = Tuple(
        Matrix(getproperty(views, Symbol(prefix, index))) for index in eachindex(supports)
    )
    return (supports=Tuple(Tuple(support) for support in supports), blocks)
end

function _cohopt_comparison_matrix(rho, ::Type{T}) where {T<:Real}
    dimension = size(rho, 1)
    comparison = zeros(T, dimension, dimension)
    for column in 1:dimension, row in 1:dimension
        comparison[row, column] =
            row == column ? abs(rho[row, column]) : -abs(rho[row, column])
    end
    return comparison
end

@doc raw"""
    is_k_incoherent(state, k; backend=NoOptimizationBackend(), strategy=:auto,
                    atol=nothing, rtol=nothing, limits=OptimizationLimits(),
                    max_subsets=10_000, max_band_search_nodes=100_000,
                    allow_densify=false)

Determine whether a density matrix has coherence number at most `k`.

The automatic path applies exact or one-sided theorem certificates first:
diagonality/trivial endpoints, the comparison-matrix criterion, the purity
ball, a corrected bounded bandwidth recognizer, and the dephasing criterion.
If these do not decide the question, it solves the exact factor-width
feasibility model

```math
\rho = \sum_{|S|=k} E_S A_S E_S^*, \qquad A_S \succeq 0.
```

The missing `IskCoherent` dependency referenced by pinned QETLAB is therefore
reconstructed directly from the primary definition: `k`-incoherent positive
matrices are precisely the factor-width-`k` cone. No CVX expression is
injected at runtime.
"""
function is_k_incoherent(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    strategy::Symbol=:auto,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    max_subsets::Integer=10_000,
    max_band_search_nodes::Integer=100_000,
    allow_densify::Bool=false,
)
    strategy = _cohopt_strategy(strategy)
    subset_limit = _cohopt_positive_limit(max_subsets, "max_subsets")
    band_node_limit = _cohopt_positive_limit(max_band_search_nodes, "max_band_search_nodes")
    data = _cohopt_state(state; atol, rtol, allow_densify)
    coherence_level = _cohopt_k(k, data.dimension)
    rho = data.rho
    T = data.real_type
    purity_value = T(real(tr(rho * rho)))
    tolerance = data.atol + data.rtol * max(one(T), maximum(abs, data.spectrum))

    if strategy === :auto
        if coherence_level == data.dimension
            decomposition = (supports=(Tuple(1:data.dimension),), blocks=(copy(rho),))
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :trivial_full_dimension,
                :full_dimension_factor_width,
                true,
                "every density matrix has factor width at most its dimension";
                decomposition,
            )
        end
        offdiagonal_max = maximum(
            (
                abs(rho[row, column]) for column in 2:data.dimension for
                row in 1:(column - 1)
            );
            init=zero(T),
        )
        if iszero(offdiagonal_max)
            supports = Tuple((index,) for index in 1:data.dimension)
            blocks = Tuple(
                reshape(T[real(rho[index, index])], 1, 1) for index in 1:data.dimension
            )
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :diagonal,
                :explicit_incoherent_decomposition,
                true,
                "the density matrix is exactly diagonal";
                decomposition=(supports=supports, blocks=blocks),
                margin=zero(T),
            )
        elseif coherence_level == 1
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_false,
                false,
                coherence_level,
                :offdiagonal_entry,
                :exact_one_incoherence_characterization,
                true,
                "a 1-incoherent density matrix must be diagonal";
                margin=offdiagonal_max,
            )
        end

        comparison = _cohopt_comparison_matrix(rho, T)
        comparison_minimum = minimum(eigvals(Symmetric(comparison)))
        if comparison_minimum > tolerance || iszero(comparison_minimum)
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :comparison_matrix,
                :comparison_matrix_psd,
                true,
                "the comparison matrix is positive semidefinite";
                margin=comparison_minimum,
                diagnostics=(comparison_minimum, comparison_tolerance=tolerance),
            )
        elseif coherence_level == 2 && comparison_minimum < -tolerance
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_false,
                false,
                coherence_level,
                :comparison_matrix,
                :comparison_matrix_2_incoherence_characterization,
                true,
                "the comparison matrix violates the necessary-and-sufficient " *
                "2-incoherence criterion";
                margin=(-comparison_minimum),
                diagnostics=(comparison_minimum, comparison_tolerance=tolerance),
            )
        end

        if coherence_level > 2 && data.dimension > 1
            threshold = inv(T(data.dimension - 1))
            if purity_value < threshold - tolerance || purity_value == threshold
                return _cohopt_criterion(
                    data,
                    :k_incoherent,
                    :certified_true,
                    true,
                    coherence_level,
                    :purity_ball,
                    :purity_sufficient_condition,
                    true,
                    "the state lies in the theorem-backed purity ball";
                    margin=threshold - purity_value,
                    diagnostics=(purity_threshold=threshold,),
                )
            end
        end

        band = _cohopt_band_ordering(rho, coherence_level; max_search_nodes=band_node_limit)
        if band.status === :found
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :band_ordering,
                :factor_width_bandwidth_upper_bound,
                true,
                "a symmetric permutation with bandwidth at most k was found";
                band_ordering=band.ordering,
                diagnostics=(band_search_status=band.status, band_search_nodes=band.nodes),
            )
        end

        alpha = T(data.dimension - coherence_level) / T(data.dimension - 1)
        dephased = Diagonal(diag(rho))
        dephasing_matrix = rho - alpha * dephased
        dephasing_minimum = minimum(real, eigvals(Hermitian(dephasing_matrix)))
        if dephasing_minimum > tolerance || iszero(dephasing_minimum)
            return _cohopt_criterion(
                data,
                :k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :dephasing_criterion,
                :dephasing_psd_sufficient_condition,
                true,
                "the dephasing sufficient condition is positive semidefinite";
                margin=dephasing_minimum,
                diagnostics=(
                    dephasing_coefficient=alpha,
                    dephasing_minimum,
                    band_search_status=band.status,
                    band_search_nodes=band.nodes,
                ),
            )
        end
    end

    problem, limit_message = _cohopt_k_incoherence_model(
        rho, coherence_level, T, limits, subset_limit
    )
    if problem === nothing
        return _cohopt_criterion(
            data,
            :k_incoherent,
            :resource_limit,
            nothing,
            coherence_level,
            :sdp,
            nothing,
            false,
            limit_message;
            diagnostics=(max_subsets=subset_limit,),
        )
    end
    optimization = solve_optimization(problem, backend)
    supports = problem.metadata.supports
    return _cohopt_solver_criterion(
        data,
        :k_incoherent,
        coherence_level,
        problem,
        optimization;
        decomposition_builder=views ->
            _cohopt_block_decomposition(views, supports, :k_block_),
    )
end

function _cohopt_absolute_threshold(value, threshold, tolerance)
    value == threshold && return :equal
    value < threshold - tolerance && return :below
    value > threshold + tolerance && return :above
    return :boundary
end

"""
    is_absolutely_k_incoherent(state, k; backend=NoOptimizationBackend(), ...)

Classify absolute `k`-incoherence from the spectrum using the rank, maximal
eigenvalue, purity, and exact `k = d-1` semidefinite criteria of Johnston,
Moein, Pereira, and Plosker. The function returns `unknown` where the cited
results are only one-sided. Invalid density matrices are rejected before any
theorem branch; unlike the pinned MATLAB control flow, an invalid input cannot
be overwritten by a later positive branch.
"""
function is_absolutely_k_incoherent(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    strategy::Symbol=:auto,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    allow_densify::Bool=false,
)
    strategy = _cohopt_strategy(strategy)
    data = _cohopt_state(state; atol, rtol, allow_densify)
    coherence_level = _cohopt_k(k, data.dimension)
    T = data.real_type
    dimension = data.dimension
    spectrum = sort(copy(data.spectrum); rev=true)
    purity_value = sum(abs2, spectrum)
    tolerance = data.atol + data.rtol * max(one(T), maximum(abs, spectrum))

    if strategy === :auto
        coherence_level == dimension && return _cohopt_criterion(
            data,
            :absolutely_k_incoherent,
            :certified_true,
            true,
            coherence_level,
            :trivial_full_dimension,
            :full_dimension_factor_width,
            true,
            "every density matrix is absolutely d-incoherent",
        )
        if coherence_level == 1
            deviations = abs.(spectrum .- inv(T(dimension)))
            maximum_deviation = maximum(deviations)
            if iszero(maximum_deviation)
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :certified_true,
                    true,
                    coherence_level,
                    :maximally_mixed,
                    :absolute_one_incoherence_characterization,
                    true,
                    "the state is exactly maximally mixed";
                    margin=zero(T),
                )
            elseif maximum_deviation > tolerance
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :certified_false,
                    false,
                    coherence_level,
                    :maximally_mixed,
                    :absolute_one_incoherence_characterization,
                    true,
                    "only the maximally mixed state is absolutely 1-incoherent";
                    margin=maximum_deviation,
                )
            end
            return _cohopt_criterion(
                data,
                :absolutely_k_incoherent,
                :boundary,
                nothing,
                coherence_level,
                :maximally_mixed,
                nothing,
                false,
                "the spectrum is within the numerical boundary of the maximally " *
                "mixed state";
                diagnostics=(maximum_deviation, boundary_tolerance=tolerance),
            )
        end

        if data.rank_lower == data.rank_upper
            rank_value = data.rank_lower
            if rank_value <= dimension - coherence_level
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :certified_false,
                    false,
                    coherence_level,
                    :rank,
                    :absolute_k_rank_necessary_condition,
                    true,
                    "the rank is too small for absolute k-incoherence";
                    margin=T(dimension - coherence_level + 1 - rank_value),
                    diagnostics=(rank=rank_value,),
                )
            elseif rank_value == dimension - coherence_level + 1
                nonzero = spectrum[spectrum .> data.rank_tolerance]
                spread = maximum(nonzero) - minimum(nonzero)
                if iszero(spread)
                    return _cohopt_criterion(
                        data,
                        :absolutely_k_incoherent,
                        :certified_true,
                        true,
                        coherence_level,
                        :equal_nonzero_spectrum,
                        :tight_rank_absolute_k_state,
                        true,
                        "the state has the tight rank and equal nonzero eigenvalues";
                        margin=zero(T),
                    )
                elseif spread < tolerance
                    return _cohopt_criterion(
                        data,
                        :absolutely_k_incoherent,
                        :boundary,
                        nothing,
                        coherence_level,
                        :equal_nonzero_spectrum,
                        nothing,
                        false,
                        "the nonzero eigenvalues are numerically near the tight-rank " *
                        "equality branch";
                        diagnostics=(nonzero_spread=spread, boundary_tolerance=tolerance),
                    )
                end
            end
        elseif data.rank_upper <= dimension - coherence_level
            return _cohopt_criterion(
                data,
                :absolutely_k_incoherent,
                :certified_false,
                false,
                coherence_level,
                :rank,
                :absolute_k_rank_necessary_condition,
                true,
                "even the tolerance upper rank is too small";
                diagnostics=(rank_lower=data.rank_lower, rank_upper=data.rank_upper),
            )
        end

        largest = first(spectrum)
        sufficient_threshold = inv(T(dimension - coherence_level + 1))
        largest_class = _cohopt_absolute_threshold(largest, sufficient_threshold, tolerance)
        if largest_class in (:below, :equal)
            return _cohopt_criterion(
                data,
                :absolutely_k_incoherent,
                :certified_true,
                true,
                coherence_level,
                :maximum_eigenvalue,
                :maximum_eigenvalue_sufficient_condition,
                true,
                "the maximal eigenvalue satisfies the absolute k-incoherence " *
                "sufficient condition";
                margin=sufficient_threshold - largest,
                diagnostics=(maximum_eigenvalue_threshold=sufficient_threshold,),
            )
        end

        if coherence_level == 2
            purity_threshold = inv(T(dimension - 1))
            purity_class = _cohopt_absolute_threshold(
                purity_value, purity_threshold, tolerance
            )
            if purity_class in (:below, :equal)
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :certified_true,
                    true,
                    coherence_level,
                    :purity,
                    if dimension <= 3
                        :absolute_two_incoherence_characterization
                    else
                        :absolute_two_incoherence_sufficient_condition
                    end,
                    true,
                    "the spectrum satisfies the absolute 2-incoherence purity bound";
                    margin=purity_threshold - purity_value,
                    diagnostics=(purity_threshold=purity_threshold,),
                )
            elseif dimension <= 3 && purity_class === :above
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :certified_false,
                    false,
                    coherence_level,
                    :purity,
                    :absolute_two_incoherence_characterization,
                    true,
                    "the low-dimensional necessary-and-sufficient purity bound is " *
                    "violated";
                    margin=purity_value - purity_threshold,
                    diagnostics=(purity_threshold=purity_threshold,),
                )
            elseif dimension <= 3
                return _cohopt_criterion(
                    data,
                    :absolutely_k_incoherent,
                    :boundary,
                    nothing,
                    coherence_level,
                    :purity,
                    nothing,
                    false,
                    "the purity lies within the numerical boundary of the exact " *
                    "low-dimensional criterion";
                    diagnostics=(purity_threshold, boundary_tolerance=tolerance),
                )
            end
        end
    end

    if coherence_level != dimension - 1
        return _cohopt_criterion(
            data,
            :absolutely_k_incoherent,
            :unknown,
            nothing,
            coherence_level,
            :known_theorems,
            nothing,
            false,
            "the implemented primary-source criteria are one-sided for this " *
            "(dimension, k) pair",
        )
    end

    necessary_threshold = one(T) - inv(T(dimension))
    largest = first(spectrum)
    largest_class = _cohopt_absolute_threshold(largest, necessary_threshold, tolerance)
    if strategy === :auto && largest_class === :above
        return _cohopt_criterion(
            data,
            :absolutely_k_incoherent,
            :certified_false,
            false,
            coherence_level,
            :maximum_eigenvalue,
            :absolute_nminus1_necessary_condition,
            true,
            "the maximal eigenvalue violates the absolute (d-1)-incoherence " *
            "necessary condition";
            margin=largest - necessary_threshold,
            diagnostics=(maximum_eigenvalue_threshold=necessary_threshold,),
        )
    end
    problem, limit_message = _cohopt_absolute_nminus1_model(spectrum, limits)
    if problem === nothing
        return _cohopt_criterion(
            data,
            :absolutely_k_incoherent,
            :resource_limit,
            nothing,
            coherence_level,
            :sdp,
            nothing,
            false,
            limit_message,
        )
    end
    optimization = solve_optimization(problem, backend)
    return _cohopt_solver_criterion(
        data,
        :absolutely_k_incoherent,
        coherence_level,
        problem,
        optimization;
        decomposition_builder=views ->
            (theorem_matrix=Matrix(getproperty(views, :absolute_nminus1_matrix)),),
    )
end

function _cohopt_fixed_cardinality_distribution(
    marginals::AbstractVector{T}, cardinality::Int
) where {T<:AbstractFloat}
    isempty(marginals) && return if cardinality == 0
        (weights=T[one(T)], subsets=[Int[]])
    else
        throw(ArgumentError("nonzero cardinality needs nonempty marginals"))
    end
    cumulative = cumsum(marginals)
    tolerance = 64eps(T) * max(one(T), T(cardinality))
    abs(last(cumulative) - T(cardinality)) <= tolerance || throw(
        ArgumentError(
            "internal k-support marginals sum to $(last(cumulative)), not $cardinality"
        ),
    )
    all(value -> -tolerance <= value <= one(T) + tolerance, marginals) ||
        throw(ArgumentError("internal k-support marginals must lie in [0, 1]"))
    boundaries = T[zero(T), one(T)]
    for value in @view cumulative[1:(end - 1)]
        fractional = value - floor(value)
        (iszero(fractional) || fractional == one(T)) && continue
        push!(boundaries, fractional)
    end
    sort!(unique!(boundaries))
    weights = T[]
    subsets = Vector{Vector{Int}}()
    for interval in 1:(length(boundaries) - 1)
        left = boundaries[interval]
        right = boundaries[interval + 1]
        right > left || continue
        sample = (left + right) / T(2)
        subset = Int[]
        for offset in 0:(cardinality - 1)
            target = sample + T(offset)
            index = searchsortedfirst(cumulative, target)
            index <= length(marginals) ||
                throw(ArgumentError("internal systematic sampling exceeded its support"))
            push!(subset, index)
        end
        length(unique(subset)) == cardinality || throw(
            ArgumentError("internal systematic sampling produced a repeated support index"),
        )
        existing = findfirst(==(subset), subsets)
        if existing === nothing
            push!(subsets, subset)
            push!(weights, right - left)
        else
            weights[existing] += right - left
        end
    end
    return (weights=weights, subsets=subsets)
end

function _cohopt_pure_factor_width_certificate(
    vector::AbstractVector{<:Number}, k::Int, tolerance
)
    dimension = length(vector)
    T = typeof(abs(zero(eltype(vector))))
    magnitudes = T.(abs.(vector))
    permutation = sortperm(magnitudes; rev=true)
    sorted_magnitudes = magnitudes[permutation]
    norm_squared = sum(abs2, sorted_magnitudes)
    support_count = count(!iszero, sorted_magnitudes)
    if support_count <= k
        support = sort(permutation[1:support_count])
        local_vector = vector[support]
        block = local_vector * adjoint(local_vector)
        zero_noise = zeros(eltype(vector), dimension, dimension)
        return (
            value=zero(T),
            free_state=vector * adjoint(vector),
            noise_state=nothing,
            unnormalized_noise=zero_noise,
            decomposition=(
                supports=(Tuple(support),),
                blocks=(block,),
                weights=(one(T),),
                atoms=(copy(vector),),
            ),
            branch_index=max(1, min(k, support_count)),
            branch_status=:free,
            certificate_minimum_eigenvalue=zero(T),
            k_support_norm_squared=norm_squared,
        )
    end

    branch_index = 1
    suffix_sums = reverse(cumsum(reverse(sorted_magnitudes)))
    for candidate in k:-1:2
        average = suffix_sums[candidate] / T(k - candidate + 1)
        if sorted_magnitudes[candidate - 1] >= average
            branch_index = candidate
            break
        end
    end
    tail_cardinality = k - branch_index + 1
    tail_sum = suffix_sums[branch_index]
    beta = tail_sum / T(tail_cardinality)
    beta > zero(T) ||
        throw(ArgumentError("internal pure-state k-support tail average vanished"))
    k_support_norm_squared =
        sum(abs2, @view sorted_magnitudes[1:(branch_index - 1)]; init=zero(T)) +
        tail_sum^2 / T(tail_cardinality)
    value = k_support_norm_squared - norm_squared
    if value < zero(T)
        value >= -tolerance || throw(
            ArgumentError(
                "internal pure-state robustness formula produced negative value $value"
            ),
        )
        value = zero(T)
    end

    tail_positions = collect(branch_index:dimension)
    marginals = sorted_magnitudes[tail_positions] ./ beta
    distribution = _cohopt_fixed_cardinality_distribution(marginals, tail_cardinality)
    large_original = permutation[1:(branch_index - 1)]
    phases = similar(vector)
    for index in eachindex(vector)
        phases[index] = if iszero(magnitudes[index])
            one(eltype(vector))
        else
            vector[index] / magnitudes[index]
        end
    end
    free_unnormalized = zeros(eltype(vector), dimension, dimension)
    supports = Tuple[]
    blocks = Matrix{eltype(vector)}[]
    atoms = Vector{eltype(vector)}[]
    for (weight, tail_subset) in zip(distribution.weights, distribution.subsets)
        selected_sorted_positions = tail_positions[tail_subset]
        selected_original = permutation[selected_sorted_positions]
        support = sort!(vcat(collect(large_original), collect(selected_original)))
        atom = zeros(eltype(vector), dimension)
        for original in large_original
            atom[original] = vector[original]
        end
        for original in selected_original
            atom[original] = beta * phases[original]
        end
        local_block = weight .* (atom[support] * adjoint(atom[support]))
        free_unnormalized[support, support] .+= local_block
        push!(supports, Tuple(support))
        push!(blocks, local_block)
        push!(atoms, atom)
    end
    rho = vector * adjoint(vector)
    noise = free_unnormalized - rho
    certificate_minimum = try
        minimum(real, eigvals(Hermitian(noise)))
    catch error
        error isa MethodError || rethrow()
        nothing
    end
    certificate_minimum === nothing ||
        certificate_minimum >= -tolerance ||
        throw(
            ArgumentError(
                "constructed pure-state robustness certificate is not PSD within " *
                "tolerance; minimum eigenvalue is $certificate_minimum",
            ),
        )
    free_trace = real(tr(free_unnormalized))
    free_state = free_unnormalized / free_trace
    noise_state = value > tolerance ? noise / value : nothing
    branch_gap = if branch_index == 1
        nothing
    else
        sorted_magnitudes[branch_index - 1] - beta
    end
    branch_status = if branch_gap === nothing
        :stable
    elseif iszero(branch_gap)
        :exact_equality
    elseif branch_gap <= tolerance
        :near_boundary
    else
        :stable
    end
    return (
        value,
        free_state,
        noise_state,
        unnormalized_noise=noise,
        decomposition=(
            supports=Tuple(supports),
            blocks=Tuple(blocks),
            weights=Tuple(distribution.weights),
            atoms=Tuple(atoms),
        ),
        branch_index,
        branch_status,
        certificate_minimum_eigenvalue=(
            certificate_minimum === nothing ? nothing : T(certificate_minimum)
        ),
        k_support_norm_squared=T(k_support_norm_squared),
    )
end

function _cohopt_is_exactly_diagonal(rho)
    dimension = size(rho, 1)
    return all(iszero(rho[row, column]) for column in 2:dimension for row in 1:(column - 1))
end

function _cohopt_solver_measure_failure(
    data, quantity::Symbol, method::Symbol, problem, optimization; k=nothing
)
    status = _cohopt_optimization_status(optimization.status)
    return _cohopt_measure(
        data,
        quantity,
        status,
        nothing,
        method,
        false,
        "the solver did not return a usable primal optimum";
        k,
        problem,
        optimization,
        diagnostics=(
            solver_status=status,
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
        ),
    )
end

function _cohopt_solver_value(data, optimization)
    optimization.status in (OptimizationOptimal, OptimizationFeasible) ||
        return nothing, nothing, false
    optimization.primal === nothing && return nothing, :missing_primal, false
    optimization.objective_value === nothing && return nothing, :missing_objective, false
    value = data.real_type(optimization.objective_value)
    tolerance = max(
        data.atol + data.rtol,
        if optimization.primal_residual === nothing
            zero(data.real_type)
        else
            data.real_type(optimization.primal_residual)
        end,
        64eps(data.real_type),
    )
    value < -tolerance && return nothing, :negative_objective, false
    corrected = value < zero(value)
    if value < zero(value)
        value = zero(value)
    end
    return value, tolerance, corrected
end

@doc raw"""
    robustness_coherence(state; backend=NoOptimizationBackend(), strategy=:auto, ...)

Compute the robustness of coherence. Exact pure-state and qubit certificates
retain an optimal incoherent state and positive noise. General mixed states use
the solver-neutral SDP

```math
\min \mathrm{tr}(X)
\quad\mathrm{such\ that}\quad
X \succeq 0,\quad \rho + X\ \mathrm{is\ diagonal}.
```
"""
function robustness_coherence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}};
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    strategy::Symbol=:auto,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    allow_densify::Bool=false,
)
    strategy = _cohopt_strategy(strategy)
    data = _cohopt_state(state; atol, rtol, allow_densify)
    rho = data.rho
    T = data.real_type
    tolerance = data.atol + data.rtol * max(one(T), maximum(abs, data.spectrum))
    if strategy === :auto
        if _cohopt_is_exactly_diagonal(rho)
            zero_noise = zeros(eltype(rho), data.dimension, data.dimension)
            return _cohopt_measure(
                data,
                :robustness_coherence,
                :analytic_exact,
                zero(T),
                :diagonal,
                true,
                "an incoherent state has zero robustness";
                free_state=copy(rho),
                noise_state=nothing,
                unnormalized_noise=zero_noise,
                diagnostics=(zero_case=true,),
            )
        elseif data.pure_vector !== nothing
            certificate = _cohopt_pure_factor_width_certificate(
                data.pure_vector, 1, tolerance
            )
            return _cohopt_measure(
                data,
                :robustness_coherence,
                :analytic_exact,
                certificate.value,
                :pure_state_k_support_norm,
                true,
                "the pure-state robustness equals the l1 coherence and an " *
                "explicit optimal decomposition was constructed";
                free_state=certificate.free_state,
                noise_state=certificate.noise_state,
                unnormalized_noise=certificate.unnormalized_noise,
                decomposition=certificate.decomposition,
                diagnostics=(
                    zero_case=iszero(certificate.value),
                    branch_index=certificate.branch_index,
                    branch_status=certificate.branch_status,
                    certificate_minimum_eigenvalue=certificate.certificate_minimum_eigenvalue,
                ),
            )
        elseif data.dimension == 2
            coherence = abs(rho[1, 2])
            noise = Matrix{eltype(rho)}([coherence -rho[1, 2]; -rho[2, 1] coherence])
            free_unnormalized = rho + noise
            value = T(2coherence)
            free_state = free_unnormalized / real(tr(free_unnormalized))
            noise_state = iszero(value) ? nothing : noise / value
            return _cohopt_measure(
                data,
                :robustness_coherence,
                :analytic_exact,
                value,
                :qubit,
                true,
                "the qubit robustness equals twice the off-diagonal magnitude";
                free_state,
                noise_state,
                unnormalized_noise=noise,
                diagnostics=(
                    zero_case=iszero(value),
                    certificate_minimum_eigenvalue=minimum(real, eigvals(Hermitian(noise))),
                ),
            )
        end
    end
    problem, limit_message = _cohopt_robustness_model(rho, T, limits)
    problem === nothing && return _cohopt_measure(
        data,
        :robustness_coherence,
        :resource_limit,
        nothing,
        :sdp,
        false,
        limit_message,
    )
    optimization = solve_optimization(problem, backend)
    value, zero_tolerance, objective_roundoff_correction = _cohopt_solver_value(
        data, optimization
    )
    value === nothing && return _cohopt_solver_measure_failure(
        data, :robustness_coherence, :sdp, problem, optimization
    )
    views = optimization.primal.views
    noise = Matrix(getproperty(views, :coherence_noise))
    free_unnormalized = Matrix(getproperty(views, :incoherent_free_unnormalized))
    free_trace = real(tr(free_unnormalized))
    free_trace > zero(T) || return _cohopt_measure(
        data,
        :robustness_coherence,
        :inconsistent,
        nothing,
        :sdp,
        false,
        "the optimizer returned a nonpositive free-state trace";
        problem,
        optimization,
    )
    free_state = free_unnormalized / free_trace
    zero_case = value <= zero_tolerance
    noise_state = zero_case ? nothing : noise / value
    return _cohopt_measure(
        data,
        :robustness_coherence,
        _cohopt_optimization_status(optimization.status),
        value,
        :sdp,
        false,
        "the solver returned a residual-checked robustness primal";
        free_state,
        noise_state,
        unnormalized_noise=noise,
        problem,
        optimization,
        diagnostics=(
            zero_case,
            zero_tolerance,
            objective_roundoff_correction,
            free_trace,
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
        ),
    )
end

function _cohopt_trace_pure(vector, tolerance)
    T = typeof(abs(zero(eltype(vector))))
    magnitudes = T.(abs.(vector))
    permutation = sortperm(magnitudes; rev=true)
    sorted = magnitudes[permutation]
    dimension = length(vector)
    support_count = count(!iszero, sorted)
    support_count <= 1 && return (
        value=zero(T),
        probabilities=real.(diag(vector * adjoint(vector))),
        branch_index=1,
        branch_status=:free,
        q=zero(T),
    )
    selected = nothing
    near_boundary = false
    for candidate in 1:dimension
        s = sum(@view sorted[1:candidate])
        m = if candidate == dimension
            zero(T)
        else
            sum(abs2, @view sorted[(candidate + 1):dimension])
        end
        r = s^2 - one(T) - T(candidate) * m
        discriminant = r^2 + 4T(candidate) * m * s^2
        discriminant >= -tolerance || throw(
            ArgumentError(
                "pure trace-distance formula produced negative discriminant " *
                "$discriminant",
            ),
        )
        discriminant < zero(T) && (discriminant = zero(T))
        q = (r + sqrt(discriminant)) / (2T(candidate) * s)
        gap = sorted[candidate] - q
        abs(gap) <= tolerance && (near_boundary = true)
        if gap > zero(T)
            denominator = s - T(candidate) * q
            denominator > zero(T) || continue
            probabilities_sorted = zeros(T, dimension)
            for index in 1:candidate
                probabilities_sorted[index] = (sorted[index] - q) / denominator
            end
            value = 2q / denominator
            selected = (
                value,
                probabilities_sorted,
                branch_index=candidate,
                branch_status=near_boundary ? :near_boundary : :stable,
                q,
            )
        end
    end
    selected === nothing &&
        throw(ArgumentError("pure trace-distance theorem found no admissible branch"))
    probabilities = zeros(T, dimension)
    probabilities[permutation] .= selected.probabilities_sorted
    return (
        value=selected.value,
        probabilities,
        branch_index=selected.branch_index,
        branch_status=selected.branch_status,
        q=selected.q,
    )
end

"""
    trace_distance_coherence(state; backend=NoOptimizationBackend(),
                             strategy=:auto, ...)

Compute the minimum trace norm distance to a diagonal density matrix. Pure
states use the exact Chen--Grogan--Johnston--Li--Plosker formula, qubit mixed
states use their closed form, and general mixed states use a Hermitian
positive/negative-part SDP. The closest diagonal state is always retained.
"""
function trace_distance_coherence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}};
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    strategy::Symbol=:auto,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    allow_densify::Bool=false,
)
    strategy = _cohopt_strategy(strategy)
    data = _cohopt_state(state; atol, rtol, allow_densify)
    rho = data.rho
    T = data.real_type
    tolerance = data.atol + data.rtol * max(one(T), maximum(abs, data.spectrum))
    if strategy === :auto
        if data.dimension == 1
            return _cohopt_measure(
                data,
                :trace_distance_coherence,
                :analytic_exact,
                zero(T),
                :dimension_one,
                true,
                "the only one-dimensional state is incoherent";
                free_state=copy(rho),
                diagnostics=(branch_status=:free,),
            )
        elseif data.pure_vector !== nothing
            pure = _cohopt_trace_pure(data.pure_vector, tolerance)
            free_state = Matrix(Diagonal(pure.probabilities))
            return _cohopt_measure(
                data,
                :trace_distance_coherence,
                :analytic_exact,
                pure.value,
                :pure_state,
                true,
                "the exact pure-state trace-distance formula was used";
                free_state,
                diagnostics=(
                    branch_index=pure.branch_index,
                    branch_status=pure.branch_status,
                    q=pure.q,
                    free_trace=sum(pure.probabilities),
                ),
            )
        elseif data.dimension == 2
            value = T(2abs(rho[1, 2]))
            free_state = Matrix(Diagonal(real.(diag(rho))))
            return _cohopt_measure(
                data,
                :trace_distance_coherence,
                :analytic_exact,
                value,
                :qubit,
                true,
                "the qubit trace distance equals twice the off-diagonal magnitude";
                free_state,
                diagnostics=(free_trace=real(tr(free_state)),),
            )
        elseif _cohopt_is_exactly_diagonal(rho)
            return _cohopt_measure(
                data,
                :trace_distance_coherence,
                :analytic_exact,
                zero(T),
                :diagonal,
                true,
                "the state is already incoherent";
                free_state=copy(rho),
                diagnostics=(branch_status=:free,),
            )
        end
    end
    problem, limit_message = _cohopt_trace_model(rho, T, limits)
    problem === nothing && return _cohopt_measure(
        data,
        :trace_distance_coherence,
        :resource_limit,
        nothing,
        :sdp,
        false,
        limit_message,
    )
    optimization = solve_optimization(problem, backend)
    value, objective_tolerance, objective_roundoff_correction = _cohopt_solver_value(
        data, optimization
    )
    value === nothing && return _cohopt_solver_measure_failure(
        data, :trace_distance_coherence, :sdp, problem, optimization
    )
    views = optimization.primal.views
    positive = Matrix(getproperty(views, :trace_positive_part))
    negative = Matrix(getproperty(views, :trace_negative_part))
    free_state = Matrix(getproperty(views, :trace_closest_incoherent))
    return _cohopt_measure(
        data,
        :trace_distance_coherence,
        _cohopt_optimization_status(optimization.status),
        value,
        :sdp,
        false,
        "the solver returned a residual-checked trace-norm primal";
        free_state,
        positive_part=positive,
        negative_part=negative,
        problem,
        optimization,
        diagnostics=(
            objective_tolerance,
            objective_roundoff_correction,
            free_trace=real(tr(free_state)),
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
        ),
    )
end

@doc raw"""
    generalized_robustness_k_coherence(state, k; backend=NoOptimizationBackend(),
                                       strategy=:auto, ...)

Compute generalized robustness relative to the free set of states with
coherence number at most `k`. This supplies the capability missing from the
pinned QETLAB checkout: its `IskCoherent` constraint is replaced by the
primary-definition factor-width cone

```math
\rho + X = \sum_{|S|=k} E_S A_S E_S^*,
\qquad X \succeq 0,\quad A_S \succeq 0.
```

Pure states use an exact `k`-support-norm formula together with an explicit
fixed-cardinality marginal decomposition, so the optimal free state, noise,
and `k`-sparse blocks are all retained.
"""
function generalized_robustness_k_coherence(
    state::Union{AbstractVector{<:Number},AbstractMatrix{<:Number}},
    k;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    strategy::Symbol=:auto,
    atol=nothing,
    rtol=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
    max_subsets::Integer=10_000,
    max_band_search_nodes::Integer=100_000,
    allow_densify::Bool=false,
)
    strategy = _cohopt_strategy(strategy)
    subset_limit = _cohopt_positive_limit(max_subsets, "max_subsets")
    band_node_limit = _cohopt_positive_limit(max_band_search_nodes, "max_band_search_nodes")
    data = _cohopt_state(state; atol, rtol, allow_densify)
    coherence_level = _cohopt_k(k, data.dimension)
    rho = data.rho
    T = data.real_type
    tolerance = data.atol + data.rtol * max(one(T), maximum(abs, data.spectrum))
    if strategy === :auto
        if coherence_level == data.dimension
            zero_noise = zeros(eltype(rho), data.dimension, data.dimension)
            return _cohopt_measure(
                data,
                :generalized_robustness_k_coherence,
                :analytic_exact,
                zero(T),
                :full_dimension,
                true,
                "every state is free at k equal to the dimension";
                k=coherence_level,
                free_state=copy(rho),
                noise_state=nothing,
                unnormalized_noise=zero_noise,
                decomposition=(supports=(Tuple(1:data.dimension),), blocks=(copy(rho),)),
                diagnostics=(zero_case=true,),
            )
        elseif data.pure_vector !== nothing
            certificate = _cohopt_pure_factor_width_certificate(
                data.pure_vector, coherence_level, tolerance
            )
            return _cohopt_measure(
                data,
                :generalized_robustness_k_coherence,
                :analytic_exact,
                certificate.value,
                :pure_state_k_support_norm,
                true,
                "the exact pure-state k-support-norm certificate was constructed";
                k=coherence_level,
                free_state=certificate.free_state,
                noise_state=certificate.noise_state,
                unnormalized_noise=certificate.unnormalized_noise,
                decomposition=certificate.decomposition,
                diagnostics=(
                    zero_case=iszero(certificate.value),
                    branch_index=certificate.branch_index,
                    branch_status=certificate.branch_status,
                    certificate_minimum_eigenvalue=certificate.certificate_minimum_eigenvalue,
                    k_support_norm_squared=certificate.k_support_norm_squared,
                ),
            )
        else
            free_test = is_k_incoherent(
                rho,
                coherence_level;
                strategy=:auto,
                backend=NoOptimizationBackend(),
                atol=data.atol,
                rtol=data.rtol,
                limits,
                max_subsets=subset_limit,
                max_band_search_nodes=band_node_limit,
                allow_densify=true,
            )
            if free_test.verdict === true
                zero_noise = zeros(eltype(rho), data.dimension, data.dimension)
                return _cohopt_measure(
                    data,
                    :generalized_robustness_k_coherence,
                    :analytic_exact,
                    zero(T),
                    :free_state_certificate,
                    free_test.exact,
                    "a k-incoherence certificate proves zero generalized robustness";
                    k=coherence_level,
                    free_state=copy(rho),
                    noise_state=nothing,
                    unnormalized_noise=zero_noise,
                    decomposition=free_test.decomposition,
                    diagnostics=(
                        zero_case=true,
                        free_test_method=free_test.method,
                        free_test_certificate_kind=free_test.certificate_kind,
                    ),
                )
            end
        end
    end
    problem, limit_message = _cohopt_generalized_model(
        rho, coherence_level, T, limits, subset_limit
    )
    problem === nothing && return _cohopt_measure(
        data,
        :generalized_robustness_k_coherence,
        :resource_limit,
        nothing,
        :sdp,
        false,
        limit_message;
        k=coherence_level,
    )
    optimization = solve_optimization(problem, backend)
    value, zero_tolerance, objective_roundoff_correction = _cohopt_solver_value(
        data, optimization
    )
    value === nothing && return _cohopt_solver_measure_failure(
        data,
        :generalized_robustness_k_coherence,
        :sdp,
        problem,
        optimization;
        k=coherence_level,
    )
    views = optimization.primal.views
    noise = Matrix(getproperty(views, :generalized_coherence_noise))
    free_unnormalized = Matrix(getproperty(views, :generalized_free_unnormalized))
    free_trace = real(tr(free_unnormalized))
    free_trace > zero(T) || return _cohopt_measure(
        data,
        :generalized_robustness_k_coherence,
        :inconsistent,
        nothing,
        :sdp,
        false,
        "the optimizer returned a nonpositive free-state trace";
        k=coherence_level,
        problem,
        optimization,
    )
    free_state = free_unnormalized / free_trace
    zero_case = value <= zero_tolerance
    noise_state = zero_case ? nothing : noise / value
    supports = problem.metadata.supports
    raw_blocks = Tuple(
        Matrix(getproperty(views, Symbol(:generalized_k_block_, index))) for
        index in eachindex(supports)
    )
    decomposition = (
        supports,
        unnormalized_blocks=raw_blocks,
        normalized_blocks=Tuple(block / free_trace for block in raw_blocks),
    )
    return _cohopt_measure(
        data,
        :generalized_robustness_k_coherence,
        _cohopt_optimization_status(optimization.status),
        value,
        :sdp,
        false,
        "the solver returned a residual-checked generalized robustness primal";
        k=coherence_level,
        free_state,
        noise_state,
        unnormalized_noise=noise,
        decomposition,
        problem,
        optimization,
        diagnostics=(
            zero_case,
            zero_tolerance,
            objective_roundoff_correction,
            free_trace,
            primal_residual=optimization.primal_residual,
            dual_residual=optimization.dual_residual,
        ),
    )
end
