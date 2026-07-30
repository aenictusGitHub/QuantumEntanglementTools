# Source-informed independent Julia implementation based on the executable
# contracts of QETLAB DiamondNorm.m, CBNorm.m,
# ChannelDistinguishability.m, and MaximumOutputFidelity.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# The solver-neutral formulations follow John Watrous, The Theory of Quantum
# Information, Sections 3.3.3--3.3.4. This file does not depend on JuMP or on
# a concrete optimizer.

"""
    ChannelOptimizationStatus

Package-owned status for channel-norm, channel-discrimination, and
maximum-output-fidelity computations. Analytic certificates, residual-checked
solver results, resource limits, unavailable backends, and numerical failures
remain distinct.
"""
@enum ChannelOptimizationStatus::UInt8 begin
    ChannelOptimizationAnalyticOptimal
    ChannelOptimizationSolverOptimal
    ChannelOptimizationSolverFeasible
    ChannelOptimizationBackendUnavailable
    ChannelOptimizationResourceLimit
    ChannelOptimizationNumericalBoundary
    ChannelOptimizationInvalidCertificate
    ChannelOptimizationBackendFailure
end

"""
    DiamondNormProblem

Solver-neutral Watrous primal SDP for the completely bounded trace norm
(diamond norm). `program` is a package-owned [`SemidefiniteProgram`](@ref);
`density_views` and `block_view` identify the two input density operators and
the positive block matrix used by the formulation.
"""
struct DiamondNormProblem{M,C,P,D,B,T<:Real}
    map::M
    choi_matrix::C
    input_dimension::Int
    output_dimension::Int
    program::P
    density_views::D
    block_view::B
    tolerance::T
end

"""
    ChannelNormResult

Status-rich result returned by [`diamond_norm`](@ref) and [`cb_norm`](@ref).
`value` is present only for an analytic branch or when residual-checked lower
and upper bounds agree. A numerical optimizer result is never relabeled as an
exact certificate.
"""
struct ChannelNormResult{T<:Real,P,D,W,O,R,C}
    quantity::Symbol
    status::ChannelOptimizationStatus
    value::Union{Nothing,T}
    lower_bound::Union{Nothing,T}
    upper_bound::Union{Nothing,T}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    input_dimension::Int
    output_dimension::Int
    problem::P
    density_operators::D
    witness_operator::W
    optimization_result::O
    residuals::R
    certificate::C
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

"""
    ChannelDistinguishabilityResult

Maximum-success-probability result for discriminating two validated quantum
channels with the recorded priors. `diamond_norm_result` retains the
prior-weighted channel-difference problem and all solver evidence.
"""
struct ChannelDistinguishabilityResult{T<:Real,N,R,C}
    status::ChannelOptimizationStatus
    success_probability::Union{Nothing,T}
    lower_bound::Union{Nothing,T}
    upper_bound::Union{Nothing,T}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    priors::Vector{T}
    diamond_norm_result::N
    residuals::R
    certificate::C
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

"""
    MaximumOutputFidelityProblem

Solver-neutral SDP for the maximum root fidelity between the outputs of two
completely positive maps. For the default channel contract, both maps are
also required to be trace preserving.
"""
struct MaximumOutputFidelityProblem{M0,M1,P,D,B,T<:Real}
    first_map::M0
    second_map::M1
    input_dimension::Int
    output_dimension::Int
    program::P
    density_views::D
    block_view::B
    tolerance::T
    trace_preserving_required::Bool
end

"""
    MaximumOutputFidelityResult

Status-rich maximum-output-root-fidelity result. Residual-checked input and
output states and the SDP coupling are retained when a backend supplies them.
"""
struct MaximumOutputFidelityResult{T<:Real,P,I,O,W,S,R,C}
    status::ChannelOptimizationStatus
    value::Union{Nothing,T}
    lower_bound::Union{Nothing,T}
    upper_bound::Union{Nothing,T}
    certified::Bool
    certificate_kind::Union{Nothing,Symbol}
    input_dimension::Int
    output_dimension::Int
    problem::P
    input_states::I
    output_states::O
    coupling::W
    optimization_result::S
    residuals::R
    certificate::C
    tolerance::T
    warnings::Tuple{Vararg{String}}
    message::String
end

function Base.show(io::IO, result::ChannelNormResult)
    return print(
        io,
        "ChannelNormResult(quantity=",
        result.quantity,
        ", status=",
        result.status,
        ", value=",
        result.value,
        ", certified=",
        result.certified,
        ")",
    )
end

function Base.show(io::IO, result::ChannelDistinguishabilityResult)
    return print(
        io,
        "ChannelDistinguishabilityResult(status=",
        result.status,
        ", success_probability=",
        result.success_probability,
        ", certified=",
        result.certified,
        ")",
    )
end

function Base.show(io::IO, result::MaximumOutputFidelityResult)
    return print(
        io,
        "MaximumOutputFidelityResult(status=",
        result.status,
        ", value=",
        result.value,
        ", certified=",
        result.certified,
        ")",
    )
end

function _channel_opt_positive_integer(value, name::AbstractString)
    value isa Integer && !(value isa Bool) && value > 0 ||
        throw(ArgumentError("$name must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$name is too large for Int"))
    return Int(value)
end

function _channel_opt_real_float_type(::Type{T}) where {T<:Number}
    R = typeof(float(real(zero(T))))
    R <: AbstractFloat ||
        throw(ArgumentError("channel optimization needs a floating real component type"))
    isconcretetype(R) ||
        throw(ArgumentError("channel optimization needs a concrete coefficient type"))
    return R
end

function _channel_opt_tolerances(::Type{T}, atol, rtol) where {T<:AbstractFloat}
    absolute = isnothing(atol) ? zero(T) : atol
    relative = isnothing(rtol) ? sqrt(eps(T)) : rtol
    for (value, name) in ((absolute, "atol"), (relative, "rtol"))
        value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
            throw(ArgumentError("$name must be a finite nonnegative real number"))
    end
    converted_absolute = try
        convert(T, absolute)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("atol cannot be represented in the model coefficient type"))
    end
    converted_relative = try
        convert(T, relative)
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("rtol cannot be represented in the model coefficient type"))
    end
    all(isfinite, (converted_absolute, converted_relative)) ||
        throw(ArgumentError("converted tolerances must remain finite"))
    return converted_absolute, converted_relative
end

function _channel_opt_require_square_algebras(
    map::AbstractMapRepresentation, operation::AbstractString
)
    space = operator_space(map)
    input_rows, input_columns = input_size(space)
    output_rows, output_columns = output_size(space)
    input_rows == input_columns || throw(
        ArgumentError(
            "$operation requires a square input matrix algebra; got " *
            "input_size=$(input_size(space))",
        ),
    )
    output_rows == output_columns || throw(
        ArgumentError(
            "$operation requires a square output matrix algebra; got " *
            "output_size=$(output_size(space))",
        ),
    )
    return input_rows, output_rows
end

function _channel_opt_prepare_map(
    map::AbstractMapRepresentation;
    operation::AbstractString,
    atol,
    rtol,
    max_input_dimension,
    max_output_dimension,
    max_choi_dimension,
    max_choi_entries,
)
    input_limit = _channel_opt_positive_integer(max_input_dimension, "max_input_dimension")
    output_limit = _channel_opt_positive_integer(
        max_output_dimension, "max_output_dimension"
    )
    choi_limit = _channel_opt_positive_integer(max_choi_dimension, "max_choi_dimension")
    entry_limit = _channel_opt_positive_integer(max_choi_entries, "max_choi_entries")
    input_dimension, output_dimension = _channel_opt_require_square_algebras(map, operation)
    input_dimension <= input_limit || throw(
        ArgumentError(
            "$operation input dimension $input_dimension exceeds " *
            "max_input_dimension=$input_limit",
        ),
    )
    output_dimension <= output_limit || throw(
        ArgumentError(
            "$operation output dimension $output_dimension exceeds " *
            "max_output_dimension=$output_limit",
        ),
    )
    choi_dimension_big = BigInt(input_dimension) * BigInt(output_dimension)
    choi_dimension_big <= choi_limit || throw(
        ArgumentError(
            "$operation Choi dimension $choi_dimension_big exceeds " *
            "max_choi_dimension=$choi_limit",
        ),
    )
    choi_entries = choi_dimension_big^2
    choi_entries <= entry_limit || throw(
        ArgumentError(
            "$operation Choi representation needs $choi_entries scalar entries, " *
            "exceeding max_choi_entries=$entry_limit",
        ),
    )
    raw_choi = choi_matrix(map)
    choi_dimension = Int(choi_dimension_big)
    size(raw_choi) == (choi_dimension, choi_dimension) || throw(
        DimensionMismatch(
            "$operation Choi matrix has size $(size(raw_choi)); expected " *
            "($choi_dimension, $choi_dimension)",
        ),
    )
    T = _channel_opt_real_float_type(eltype(raw_choi))
    absolute, relative = _channel_opt_tolerances(T, atol, rtol)
    converted_choi = if issparse(raw_choi)
        sparse(Complex{T}.(raw_choi))
    else
        Matrix{Complex{T}}(raw_choi)
    end
    scale = maximum(abs, converted_choi; init=one(T))
    tolerance = absolute + relative * max(one(T), scale)
    isfinite(tolerance) ||
        throw(ArgumentError("the combined channel-optimization tolerance is not finite"))
    return (
        map=map,
        choi=converted_choi,
        input_dimension=input_dimension,
        output_dimension=output_dimension,
        choi_dimension=choi_dimension,
        coefficient_type=T,
        atol=absolute,
        rtol=relative,
        requested_atol=atol,
        requested_rtol=rtol,
        tolerance=tolerance,
        max_choi_entries=entry_limit,
    )
end

function _channel_opt_affine_combination(
    name::Symbol, weighted_matrices::Tuple{Vararg{Tuple}}
)
    isempty(weighted_matrices) &&
        throw(ArgumentError("an affine combination needs at least one matrix"))
    first_matrix = weighted_matrices[1][2]
    first_matrix isa HermitianAffineMatrix || throw(
        ArgumentError("affine combination terms must be HermitianAffineMatrix values")
    )
    dimension = first_matrix.dimension
    variable_count = first_matrix.variable_count
    real_types = Type[]
    for (weight, matrix) in weighted_matrices
        weight isa Real && !(weight isa Bool) && isfinite(weight) ||
            throw(ArgumentError("affine combination weights must be finite real values"))
        matrix isa HermitianAffineMatrix ||
            throw(ArgumentError("affine combination terms must be Hermitian affine"))
        matrix.dimension == dimension ||
            throw(DimensionMismatch("affine combination matrix dimensions do not agree"))
        matrix.variable_count == variable_count ||
            throw(DimensionMismatch("affine combination variable counts do not agree"))
        push!(real_types, typeof(weight))
        push!(real_types, _optimization_complex_real_type(eltype(matrix.constant)))
    end
    T = promote_type(real_types...)
    T <: Real || throw(ArgumentError("affine combination data must promote to real"))
    constant = spzeros(Complex{T}, dimension, dimension)
    coefficients = Dict{Int,SparseMatrixCSC{Complex{T},Int}}()
    for (weight, matrix) in weighted_matrices
        converted_weight = convert(T, weight)
        constant += converted_weight .* sparse(Complex{T}.(matrix.constant))
        for term in matrix.terms
            contribution = converted_weight .* sparse(Complex{T}.(term.coefficient))
            if haskey(coefficients, term.variable)
                coefficients[term.variable] += contribution
            else
                coefficients[term.variable] = contribution
            end
        end
    end
    variables = sort!(collect(keys(coefficients)))
    kept_variables = Int[]
    kept_coefficients = SparseMatrixCSC{Complex{T},Int}[]
    for variable in variables
        coefficient = coefficients[variable]
        dropzeros!(coefficient)
        nnz(coefficient) == 0 && continue
        push!(kept_variables, variable)
        push!(kept_coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name, constant, kept_variables, kept_coefficients, variable_count
    )
end

function _channel_opt_principal_block(
    matrix::HermitianAffineMatrix, indices::UnitRange{Int}; name::Symbol
)
    first(indices) >= 1 && last(indices) <= matrix.dimension ||
        throw(BoundsError(matrix, indices))
    constant = matrix.constant[indices, indices]
    variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    for term in matrix.terms
        coefficient = term.coefficient[indices, indices]
        nnz(coefficient) == 0 && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name, constant, variables, coefficients, matrix.variable_count
    )
end

function _channel_opt_tensor_output_identity(
    density::HermitianAffineMatrix{T}, output_dimension::Int; name::Symbol
) where {T<:Real}
    identity_output = spdiagm(0 => fill(one(Complex{T}), output_dimension))
    return _optimization_apply_matrix_map(
        coefficient -> kron(coefficient, identity_output), density, name
    )
end

function _channel_opt_trace_one_equality(
    matrix::HermitianAffineMatrix{T}, name::Symbol
) where {T<:Real}
    trace_data = trace_affine(matrix)
    shifted = AffineScalar(
        trace_data.constant - one(T),
        trace_data.coefficients;
        variable_count=matrix.variable_count,
    )
    return AffineEquality(shifted, name)
end

function _channel_opt_write_hermitian_coordinates!(
    coordinates::AbstractVector{T}, first_variable::Int, matrix::AbstractMatrix{<:Number}
) where {T<:Real}
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("coordinate source must be square"))
    ishermitian(matrix) ||
        throw(ArgumentError("coordinate source must be exactly Hermitian"))
    dimension = size(matrix, 1)
    last_variable = BigInt(first_variable) + BigInt(dimension)^2 - 1
    last_variable <= length(coordinates) ||
        throw(BoundsError(coordinates, first_variable:Int(last_variable)))
    position = first_variable
    for index in 1:dimension
        coordinates[position] = convert(T, real(matrix[index, index]))
        position += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        coordinates[position] = convert(T, real(matrix[row, column]))
        position += 1
    end
    for column in 2:dimension, row in 1:(column - 1)
        coordinates[position] = convert(T, imag(matrix[row, column]))
        position += 1
    end
    return coordinates
end

function _channel_opt_diamond_objective(
    block::HermitianAffineMatrix{T}, choi::AbstractMatrix{<:Number}, choi_dimension::Int
) where {T<:Real}
    size(choi) == (choi_dimension, choi_dimension) ||
        throw(DimensionMismatch("Choi matrix has the wrong objective dimension"))
    variables = Int[]
    coefficients = T[]
    for term in block.terms
        off_diagonal = @view term.coefficient[
            1:choi_dimension, (choi_dimension + 1):(2choi_dimension)
        ]
        value = convert(T, real(dot(choi, off_diagonal)))
        iszero(value) && continue
        push!(variables, term.variable)
        push!(coefficients, value)
    end
    return AffineScalar(zero(T), variables, coefficients, block.variable_count)
end

function _channel_opt_build_diamond_problem(
    prepared; max_variables, limits::OptimizationLimits
)
    variable_limit = _channel_opt_positive_integer(max_variables, "max_variables")
    d_in = prepared.input_dimension
    d_out = prepared.output_dimension
    choi_dimension = prepared.choi_dimension
    variable_count_big = 2BigInt(d_in)^2 + (2BigInt(choi_dimension))^2
    variable_count_big <= variable_limit || throw(
        ArgumentError(
            "diamond-norm SDP needs $variable_count_big variables, exceeding " *
            "max_variables=$variable_limit",
        ),
    )
    variable_count = Int(variable_count_big)
    T = prepared.coefficient_type
    first_left = 1
    first_right = first_left + d_in^2
    first_block = first_right + d_in^2
    left_density = hermitian_variable(
        :diamond_rho_left,
        d_in;
        first_variable=first_left,
        variable_count=variable_count,
        coefficient_type=T,
    )
    right_density = hermitian_variable(
        :diamond_rho_right,
        d_in;
        first_variable=first_right,
        variable_count=variable_count,
        coefficient_type=T,
    )
    block = hermitian_variable(
        :diamond_block,
        2choi_dimension;
        first_variable=first_block,
        variable_count=variable_count,
        coefficient_type=T,
    )
    top = _channel_opt_principal_block(block, 1:choi_dimension; name=:diamond_block_top)
    bottom = _channel_opt_principal_block(
        block, (choi_dimension + 1):(2choi_dimension); name=:diamond_block_bottom
    )
    top_target = _channel_opt_tensor_output_identity(
        left_density, d_out; name=:diamond_left_tensor_identity
    )
    bottom_target = _channel_opt_tensor_output_identity(
        right_density, d_out; name=:diamond_right_tensor_identity
    )
    top_relation = _channel_opt_affine_combination(
        :diamond_top_relation, ((one(T), top), (-one(T), top_target))
    )
    bottom_relation = _channel_opt_affine_combination(
        :diamond_bottom_relation, ((one(T), bottom), (-one(T), bottom_target))
    )
    equalities = AffineEquality{T}[
        _channel_opt_trace_one_equality(left_density, :diamond_left_trace),
        _channel_opt_trace_one_equality(right_density, :diamond_right_trace),
    ]
    append!(
        equalities, hermitian_equalities(top_relation; name_prefix=:diamond_top_structure)
    )
    append!(
        equalities,
        hermitian_equalities(bottom_relation; name_prefix=:diamond_bottom_structure),
    )
    objective = _channel_opt_diamond_objective(block, prepared.choi, choi_dimension)
    feasible = zeros(T, variable_count)
    maximally_mixed = Matrix{Complex{T}}(I, d_in, d_in) / convert(T, d_in)
    diagonal_block = kron(maximally_mixed, Matrix{Complex{T}}(I, d_out, d_out))
    feasible_block = [
        diagonal_block zeros(Complex{T}, choi_dimension, choi_dimension)
        zeros(Complex{T}, choi_dimension, choi_dimension) diagonal_block
    ]
    _channel_opt_write_hermitian_coordinates!(feasible, first_left, maximally_mixed)
    _channel_opt_write_hermitian_coordinates!(feasible, first_right, maximally_mixed)
    _channel_opt_write_hermitian_coordinates!(feasible, first_block, feasible_block)
    program = SemidefiniteProgram(
        :diamond_norm_watrous_primal,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=[left_density, right_density, block],
        primal_views=[left_density, right_density, block],
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:watrous_completely_bounded_trace_norm_primal,
            input_dimension=d_in,
            output_dimension=d_out,
            choi_order=:input_tensor_output,
            source=:watrous_tqi_sections_3_3_3_and_3_3_4,
        ),
    )
    return DiamondNormProblem(
        prepared.map,
        prepared.choi,
        d_in,
        d_out,
        program,
        (left_density, right_density),
        block,
        prepared.tolerance,
    )
end

"""
    diamond_norm_problem(map; kwargs...)

Build the Watrous primal SDP

```math
\\begin{aligned}
\\text{maximize}\\quad & \\operatorname{Re}\\langle J(\\Phi),X\\rangle,\\\\
\\text{subject to}\\quad &
\\begin{bmatrix}
\\rho_0\\otimes I_{\\mathrm{out}} & X\\\\
X^\\dagger & \\rho_1\\otimes I_{\\mathrm{out}}
\\end{bmatrix}\\succeq0,\\\\
&\\rho_0,\\rho_1\\succeq0,\\quad
\\operatorname{tr}(\\rho_0)=\\operatorname{tr}(\\rho_1)=1.
\\end{aligned}
```

The package Choi order is input tensor output, so the diagonal blocks are
`rho ⊗ I`. The map may be a general two-sided representation and its input
and output dimensions may differ, but each operator space must be square.
No map is projected onto a Hermiticity-preserving or completely positive
cone.
"""
function diamond_norm_problem(
    map::AbstractMapRepresentation;
    atol=nothing,
    rtol=nothing,
    max_input_dimension=32,
    max_output_dimension=32,
    max_choi_dimension=64,
    max_choi_entries=4_096,
    max_variables=50_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _channel_opt_prepare_map(
        map;
        operation="diamond_norm",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    return _channel_opt_build_diamond_problem(
        prepared; max_variables=max_variables, limits=limits
    )
end

function _channel_opt_dense_matrix(
    matrix::AbstractMatrix,
    allow_densify::Bool,
    max_dense_entries,
    operation::AbstractString,
)
    entry_limit = _channel_opt_positive_integer(max_dense_entries, "max_dense_entries")
    BigInt(length(matrix)) <= entry_limit || throw(
        ArgumentError(
            "$operation needs $(length(matrix)) dense entries, exceeding " *
            "max_dense_entries=$entry_limit",
        ),
    )
    issparse(matrix) &&
        !allow_densify &&
        throw(
            ArgumentError(
                "$operation would densify sparse data; pass allow_densify=true explicitly"
            ),
        )
    return Matrix(matrix)
end

function _channel_opt_trace_output(
    choi::AbstractMatrix{T}, input_dimension::Int, output_dimension::Int
) where {T}
    result = zeros(T, input_dimension, input_dimension)
    @inbounds for input_column in 1:input_dimension,
        input_row in 1:input_dimension,
        output_index in 1:output_dimension

        result[input_row, input_column] += choi[
            output_index + (input_row - 1) * output_dimension,
            output_index + (input_column - 1) * output_dimension,
        ]
    end
    return result
end

function _channel_opt_trace_output(
    choi::SparseMatrixCSC{T,Int}, input_dimension::Int, output_dimension::Int
) where {T}
    rows, columns, values = findnz(choi)
    reduced_rows = Int[]
    reduced_columns = Int[]
    reduced_values = T[]
    for index in eachindex(values)
        output_row = mod(rows[index] - 1, output_dimension) + 1
        output_column = mod(columns[index] - 1, output_dimension) + 1
        output_row == output_column || continue
        push!(reduced_rows, div(rows[index] - 1, output_dimension) + 1)
        push!(reduced_columns, div(columns[index] - 1, output_dimension) + 1)
        push!(reduced_values, values[index])
    end
    return sparse(
        reduced_rows, reduced_columns, reduced_values, input_dimension, input_dimension
    )
end

function _channel_opt_cp_diagnostic(prepared; allow_densify::Bool)
    if issparse(prepared.choi) && !(prepared.map isa KrausRepresentation) && !allow_densify
        return nothing
    end
    return is_completely_positive(
        prepared.map;
        atol=prepared.requested_atol,
        rtol=prepared.requested_rtol,
        allow_densify=allow_densify,
    )
end

function _channel_opt_analytic_norm(prepared; allow_densify::Bool, max_dense_entries)
    T = prepared.coefficient_type
    choi = prepared.choi
    d_in = prepared.input_dimension
    d_out = prepared.output_dimension
    if all(iszero, choi)
        return (
            value=zero(T),
            certificate_kind=:zero_map,
            certificate=(kind=:exact_zero_choi,),
            residuals=(exact_zero=true,),
            message="the Choi matrix is exactly zero",
        )
    end
    issparse(choi) && !allow_densify && return nothing
    if d_in == 1
        dense = _channel_opt_dense_matrix(
            choi, allow_densify, max_dense_entries, "the scalar-input diamond-norm shortcut"
        )
        singular_values = svdvals(dense)
        value = sum(singular_values)
        return (
            value=convert(T, value),
            certificate_kind=:scalar_input_trace_norm,
            certificate=(kind=:scalar_input_trace_norm, singular_values=singular_values),
            residuals=(trace_norm=value,),
            message="a one-dimensional input reduces the diamond norm to a trace norm",
        )
    elseif d_out == 1
        dense = _channel_opt_dense_matrix(
            choi,
            allow_densify,
            max_dense_entries,
            "the scalar-output diamond-norm shortcut",
        )
        singular_values = svdvals(dense)
        value = maximum(singular_values)
        return (
            value=convert(T, value),
            certificate_kind=:scalar_output_operator_norm,
            certificate=(
                kind=:scalar_output_operator_norm, singular_values=singular_values
            ),
            residuals=(operator_norm=value,),
            message="a one-dimensional output reduces the diamond norm to an operator norm",
        )
    end
    complete_positivity = _channel_opt_cp_diagnostic(prepared; allow_densify=allow_densify)
    complete_positivity === nothing && return nothing
    complete_positivity.status === MatrixPredicateSatisfied || return nothing
    reduced = _channel_opt_trace_output(choi, d_in, d_out)
    dense = _channel_opt_dense_matrix(
        reduced,
        allow_densify,
        max_dense_entries,
        "the completely-positive diamond-norm shortcut",
    )
    ishermitian(dense) || return nothing
    eigenvalues = eigvals(Hermitian(dense))
    minimum_eigenvalue = minimum(real, eigenvalues)
    minimum_eigenvalue >= -prepared.tolerance || return nothing
    value = maximum(real, eigenvalues)
    return (
        value=convert(T, value),
        certificate_kind=:completely_positive_adjoint_identity,
        certificate=(
            kind=:completely_positive_adjoint_identity,
            complete_positivity=complete_positivity,
            adjoint_identity_spectrum=copy(eigenvalues),
        ),
        residuals=(
            minimum_adjoint_identity_eigenvalue=minimum_eigenvalue,
            complete_positivity_tolerance=complete_positivity.tolerance,
        ),
        message="complete positivity reduces the diamond norm to the operator norm of the adjoint applied to the identity",
    )
end

function _channel_opt_norm_result(
    quantity::Symbol,
    status::ChannelOptimizationStatus,
    T::Type{<:Real},
    input_dimension::Int,
    output_dimension::Int,
    tolerance;
    value=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    certified=false,
    certificate_kind=nothing,
    problem=nothing,
    density_operators=nothing,
    witness_operator=nothing,
    optimization_result=nothing,
    residuals=NamedTuple(),
    certificate=nothing,
    warnings=(),
    message,
)
    return ChannelNormResult(
        quantity,
        status,
        value === nothing ? nothing : convert(T, value),
        lower_bound === nothing ? nothing : convert(T, lower_bound),
        upper_bound === nothing ? nothing : convert(T, upper_bound),
        certified,
        certificate_kind,
        input_dimension,
        output_dimension,
        problem,
        density_operators,
        witness_operator,
        optimization_result,
        residuals,
        certificate,
        convert(T, tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

function _channel_opt_spectral_residual(
    matrix::AbstractMatrix,
    allow_densify::Bool,
    max_dense_entries,
    operation::AbstractString,
)
    dense = _channel_opt_dense_matrix(matrix, allow_densify, max_dense_entries, operation)
    R = typeof(real(zero(eltype(dense))))
    hermiticity = maximum(abs, dense - adjoint(dense); init=zero(R))
    iszero(hermiticity) || return (
        minimum_eigenvalue=(-convert(R, Inf)),
        positivity_violation=convert(R, Inf),
        hermiticity_residual=hermiticity,
    )
    eigenvalues = eigvals(Hermitian(dense))
    minimum_eigenvalue = minimum(real, eigenvalues)
    return (
        minimum_eigenvalue=minimum_eigenvalue,
        positivity_violation=max(zero(minimum_eigenvalue), -minimum_eigenvalue),
        hermiticity_residual=hermiticity,
    )
end

function _channel_opt_diamond_primal(
    problem::DiamondNormProblem,
    optimization::OptimizationResult,
    tolerance,
    allow_densify::Bool,
    max_dense_entries,
)
    optimization.primal === nothing && return nothing
    views = optimization.primal.views
    all(
        name -> hasproperty(views, name),
        (:diamond_rho_left, :diamond_rho_right, :diamond_block),
    ) || return nothing
    left = Matrix(getproperty(views, :diamond_rho_left))
    right = Matrix(getproperty(views, :diamond_rho_right))
    block = Matrix(getproperty(views, :diamond_block))
    choi_dimension = size(problem.choi_matrix, 1)
    size(block) == (2choi_dimension, 2choi_dimension) || return nothing
    witness = Matrix(@view block[1:choi_dimension, (choi_dimension + 1):(2choi_dimension)])
    left_psd = _channel_opt_spectral_residual(
        left, allow_densify, max_dense_entries, "diamond-norm left-density validation"
    )
    right_psd = _channel_opt_spectral_residual(
        right, allow_densify, max_dense_entries, "diamond-norm right-density validation"
    )
    block_psd = _channel_opt_spectral_residual(
        block, allow_densify, max_dense_entries, "diamond-norm block validation"
    )
    identity_output = Matrix{eltype(block)}(
        I, problem.output_dimension, problem.output_dimension
    )
    top_target = kron(left, identity_output)
    bottom_target = kron(right, identity_output)
    top = @view block[1:choi_dimension, 1:choi_dimension]
    bottom = @view block[
        (choi_dimension + 1):(2choi_dimension), (choi_dimension + 1):(2choi_dimension)
    ]
    trace_residual = max(abs(real(tr(left)) - 1), abs(real(tr(right)) - 1))
    structure_residual = max(
        maximum(abs, top - top_target; init=zero(tolerance)),
        maximum(abs, bottom - bottom_target; init=zero(tolerance)),
    )
    value = real(dot(problem.choi_matrix, witness))
    objective_residual = if optimization.objective_value === nothing
        nothing
    else
        abs(value - optimization.objective_value)
    end
    valid =
        left_psd.hermiticity_residual <= tolerance &&
        right_psd.hermiticity_residual <= tolerance &&
        block_psd.hermiticity_residual <= tolerance &&
        left_psd.positivity_violation <= tolerance &&
        right_psd.positivity_violation <= tolerance &&
        block_psd.positivity_violation <= tolerance &&
        trace_residual <= tolerance &&
        structure_residual <= tolerance &&
        (objective_residual === nothing || objective_residual <= 8tolerance)
    return (
        valid=valid,
        value=value,
        density_operators=(left, right),
        witness_operator=witness,
        residuals=(
            left_density=left_psd,
            right_density=right_psd,
            block=block_psd,
            trace_residual=trace_residual,
            structure_residual=structure_residual,
            solver_objective_residual=objective_residual,
        ),
    )
end

function _channel_opt_backend_tolerance(problem, optimization, backend)
    T = typeof(problem.tolerance)
    if backend isa JuMPBackend
        scale = if optimization.objective_value === nothing
            one(T)
        else
            max(one(T), abs(convert(T, optimization.objective_value)))
        end
        return max(
            problem.tolerance, convert(T, backend.atol) + convert(T, backend.rtol) * scale
        )
    end
    return problem.tolerance
end

function _channel_opt_solver_norm_result(
    problem::DiamondNormProblem,
    optimization::OptimizationResult,
    backend;
    quantity::Symbol,
    source_input_dimension::Int=problem.input_dimension,
    source_output_dimension::Int=problem.output_dimension,
    allow_densify::Bool,
    max_dense_entries,
)
    T = typeof(problem.tolerance)
    tolerance = _channel_opt_backend_tolerance(problem, optimization, backend)
    primal = try
        _channel_opt_diamond_primal(
            problem,
            optimization,
            8tolerance,
            allow_densify || (backend isa JuMPBackend && backend.allow_densify),
            max_dense_entries,
        )
    catch error
        error isa ArgumentError || rethrow()
        nothing
    end
    lower = if primal !== nothing && primal.valid
        max(zero(T), convert(T, primal.value))
    else
        zero(T)
    end
    upper = if optimization.objective_bound === nothing
        nothing
    else
        convert(T, optimization.objective_bound)
    end
    bound_order_valid = upper === nothing || upper + 8tolerance >= lower
    matched =
        upper !== nothing &&
        bound_order_valid &&
        abs(lower - upper) <= 16tolerance * max(one(T), abs(lower), abs(upper))
    value = matched ? (lower + upper) / 2 : nothing
    residuals = if primal === nothing
        (primal_available=false,)
    else
        merge((primal_available=true,), primal.residuals)
    end
    density_operators =
        primal === nothing || !primal.valid ? nothing : primal.density_operators
    witness_operator =
        primal === nothing || !primal.valid ? nothing : primal.witness_operator

    if optimization.status === OptimizationBackendUnavailable
        return _channel_opt_norm_result(
            quantity,
            ChannelOptimizationBackendUnavailable,
            T,
            source_input_dimension,
            source_output_dimension,
            tolerance;
            lower_bound=lower,
            problem=problem,
            optimization_result=optimization,
            residuals=residuals,
            message=optimization.message,
        )
    elseif optimization.status === OptimizationLimit
        return _channel_opt_norm_result(
            quantity,
            ChannelOptimizationResourceLimit,
            T,
            source_input_dimension,
            source_output_dimension,
            tolerance;
            lower_bound=lower,
            upper_bound=bound_order_valid ? upper : nothing,
            problem=problem,
            density_operators=density_operators,
            witness_operator=witness_operator,
            optimization_result=optimization,
            residuals=residuals,
            message="the optimizer stopped at a resource limit; validated bounds are retained",
        )
    elseif optimization.status === OptimizationOptimal &&
        primal !== nothing &&
        primal.valid &&
        matched
        return _channel_opt_norm_result(
            quantity,
            ChannelOptimizationSolverOptimal,
            T,
            source_input_dimension,
            source_output_dimension,
            tolerance;
            value=value,
            lower_bound=lower,
            upper_bound=upper,
            certified=false,
            problem=problem,
            density_operators=density_operators,
            witness_operator=witness_operator,
            optimization_result=optimization,
            residuals=residuals,
            certificate=(
                kind=:numerical_primal_dual_agreement,
                primal_residual=optimization.primal_residual,
                dual_residual=optimization.dual_residual,
                absolute_gap=optimization.absolute_gap,
            ),
            message="the optimizer returned residual-checked matching primal and dual bounds",
        )
    elseif optimization.status in (OptimizationOptimal, OptimizationFeasible) &&
        primal !== nothing &&
        primal.valid
        status = if optimization.status === OptimizationFeasible
            ChannelOptimizationSolverFeasible
        else
            ChannelOptimizationInvalidCertificate
        end
        return _channel_opt_norm_result(
            quantity,
            status,
            T,
            source_input_dimension,
            source_output_dimension,
            tolerance;
            lower_bound=lower,
            upper_bound=bound_order_valid ? upper : nothing,
            problem=problem,
            density_operators=density_operators,
            witness_operator=witness_operator,
            optimization_result=optimization,
            residuals=residuals,
            warnings=(
                "the available primal and dual bounds do not establish a common value",
            ),
            message="a residual-checked primal point is available, but optimality is not independently established",
        )
    end
    failure = optimization.status in (
        OptimizationMalformedBackend,
        OptimizationNumericalFailure,
        OptimizationUnsupported,
        OptimizationUnknown,
        OptimizationInfeasible,
        OptimizationUnbounded,
        OptimizationInconsistent,
    )
    return _channel_opt_norm_result(
        quantity,
        failure ? ChannelOptimizationBackendFailure : ChannelOptimizationInvalidCertificate,
        T,
        source_input_dimension,
        source_output_dimension,
        tolerance;
        lower_bound=lower,
        upper_bound=bound_order_valid ? upper : nothing,
        problem=problem,
        density_operators=density_operators,
        witness_operator=witness_operator,
        optimization_result=optimization,
        residuals=residuals,
        message="the optimization outcome does not establish a valid channel norm",
    )
end

"""
    diamond_norm(map; backend=NoOptimizationBackend(), kwargs...)

Compute or bound the completely bounded trace norm of a finite-dimensional
map. Exact zero maps, one-dimensional input/output spaces, and maps whose
complete positivity is established use solver-free theorems. Every other map
uses [`diamond_norm_problem`](@ref) and an explicit optional backend.

The returned [`ChannelNormResult`](@ref) retains lower and upper bounds,
primal density operators and block witness, generic solver statuses, residuals,
and optimizer metadata. Sparse spectral shortcuts require
`allow_densify=true`; the optional complex-to-real PSD translation separately
requires `JuMPBackend(...; allow_densify=true)`.
"""
function diamond_norm(
    map::AbstractMapRepresentation;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_input_dimension=32,
    max_output_dimension=32,
    max_choi_dimension=64,
    max_choi_entries=4_096,
    max_dense_entries=1_000_000,
    max_variables=50_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    prepared = _channel_opt_prepare_map(
        map;
        operation="diamond_norm",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    analytic = _channel_opt_analytic_norm(
        prepared; allow_densify=allow_densify, max_dense_entries=max_dense_entries
    )
    if analytic !== nothing
        return _channel_opt_norm_result(
            :diamond,
            ChannelOptimizationAnalyticOptimal,
            prepared.coefficient_type,
            prepared.input_dimension,
            prepared.output_dimension,
            prepared.tolerance;
            value=analytic.value,
            lower_bound=analytic.value,
            upper_bound=analytic.value,
            certified=true,
            certificate_kind=analytic.certificate_kind,
            residuals=analytic.residuals,
            certificate=analytic.certificate,
            message=analytic.message,
        )
    end
    problem = _channel_opt_build_diamond_problem(
        prepared; max_variables=max_variables, limits=limits
    )
    optimization = solve_optimization(problem.program, backend)
    return _channel_opt_solver_norm_result(
        problem,
        optimization,
        backend;
        quantity=:diamond,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
    )
end

"""
    cb_norm(map; kwargs...)

Compute the finite-dimensional completely bounded operator norm through the
duality identity
`cb_norm(map) == diamond_norm(dual_channel(map))`.
The returned result records `quantity == :completely_bounded`; any SDP in
`result.problem` is the explicit diamond-norm problem for the adjoint map.
"""
function cb_norm(
    map::AbstractMapRepresentation;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    kwargs...,
)
    source_input_dimension, source_output_dimension = _channel_opt_require_square_algebras(
        map, "cb_norm"
    )
    adjoint_map = dual_channel(map)
    inner = diamond_norm(adjoint_map; backend=backend, kwargs...)
    certificate = (
        transformation=:hilbert_schmidt_adjoint, diamond_certificate=inner.certificate
    )
    return _channel_opt_norm_result(
        :completely_bounded,
        inner.status,
        typeof(inner.tolerance),
        source_input_dimension,
        source_output_dimension,
        inner.tolerance;
        value=inner.value,
        lower_bound=inner.lower_bound,
        upper_bound=inner.upper_bound,
        certified=inner.certified,
        certificate_kind=inner.certificate_kind,
        problem=inner.problem,
        density_operators=inner.density_operators,
        witness_operator=inner.witness_operator,
        optimization_result=inner.optimization_result,
        residuals=inner.residuals,
        certificate=certificate,
        warnings=inner.warnings,
        message="completely bounded norm via the diamond norm of the Hilbert--Schmidt adjoint: " *
                inner.message,
    )
end

function _channel_opt_validate_channel(
    prepared; allow_densify::Bool, operation::AbstractString
)
    complete_positivity = _channel_opt_cp_diagnostic(prepared; allow_densify=allow_densify)
    complete_positivity === nothing && throw(
        ArgumentError(
            "$operation must establish complete positivity, which would densify " *
            "sparse Choi data; pass allow_densify=true",
        ),
    )
    complete_positivity.status === MatrixPredicateSatisfied || throw(
        DomainError(
            complete_positivity,
            "$operation requires completely positive maps; the diagnostic is " *
            "$(complete_positivity.status)",
        ),
    )
    trace_matrix = _channel_opt_trace_output(
        prepared.choi, prepared.input_dimension, prepared.output_dimension
    )
    identity_input = if issparse(trace_matrix)
        spdiagm(0 => fill(one(eltype(trace_matrix)), prepared.input_dimension))
    else
        Matrix{eltype(trace_matrix)}(I, prepared.input_dimension, prepared.input_dimension)
    end
    trace_residual = maximum(
        abs, trace_matrix - identity_input; init=zero(prepared.tolerance)
    )
    if !iszero(trace_residual)
        if trace_residual <= prepared.tolerance
            throw(
                DomainError(
                    trace_residual,
                    "$operation trace preservation lies inside a numerical boundary; " *
                    "the map is not repaired or silently accepted",
                ),
            )
        end
        throw(
            DomainError(
                trace_residual,
                "$operation requires trace-preserving channels; residual " *
                "$trace_residual exceeds tolerance $(prepared.tolerance)",
            ),
        )
    end
    return (
        complete_positivity=complete_positivity, trace_preservation_residual=trace_residual
    )
end

function _channel_opt_priors(priors, ::Type{T}, tolerance) where {T<:Real}
    values = if priors === nothing
        T[one(T) / 2, one(T) / 2]
    else
        priors isa Tuple ||
            priors isa AbstractVector ||
            throw(ArgumentError("priors must be a tuple or vector with two entries"))
        priors isa AbstractVector &&
            firstindex(priors) != 1 &&
            throw(ArgumentError("priors must use one-based indexing"))
        length(priors) == 2 ||
            throw(DimensionMismatch("priors must contain exactly two entries"))
        converted = Vector{T}(undef, 2)
        for index in 1:2
            value = priors[index]
            value isa Real && !(value isa Bool) && isfinite(value) ||
                throw(ArgumentError("priors must contain finite real values"))
            value >= zero(value) ||
                throw(DomainError(value, "priors must be nonnegative"))
            converted[index] = try
                convert(T, value)
            catch error
                error isa InexactError || rethrow()
                throw(ArgumentError("prior $index cannot be represented by $T"))
            end
        end
        converted
    end
    prior_sum = sum(values)
    residual = abs(prior_sum - one(T))
    residual <= tolerance || throw(
        ArgumentError(
            "priors must sum to one within tolerance $tolerance; residual is $residual"
        ),
    )
    return values, prior_sum, residual
end

function _channel_opt_distinguishability_result(
    status,
    ::Type{T},
    priors,
    tolerance;
    success_probability=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    certified=false,
    certificate_kind=nothing,
    diamond_norm_result=nothing,
    residuals=NamedTuple(),
    certificate=nothing,
    warnings=(),
    message,
) where {T<:Real}
    return ChannelDistinguishabilityResult(
        status,
        success_probability === nothing ? nothing : convert(T, success_probability),
        lower_bound === nothing ? nothing : convert(T, lower_bound),
        upper_bound === nothing ? nothing : convert(T, upper_bound),
        certified,
        certificate_kind,
        Vector{T}(priors),
        diamond_norm_result,
        residuals,
        certificate,
        convert(T, tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

"""
    channel_distinguishability(first_map, second_map;
                               priors=(1/2, 1/2), kwargs...)

Return the optimal success probability for one-use discrimination of two
validated channels:

```math
p_{\\mathrm{success}}
= \\frac{p_1+p_2+\\lVert p_1\\Phi_1-p_2\\Phi_2\\rVert_\\diamond}{2}.
```

Inputs must be completely positive, exactly trace preserving outside numerical
boundaries, and have equal input and output dimensions. Priors are validated
but never normalized. Missing backends retain the always-guess lower bound and
the probability upper bound.

The pinned QETLAB routine omits the affine `(1 + norm)/2` conversion for
non-deterministic priors. This implementation records and corrects that defect.
"""
function channel_distinguishability(
    first_map::AbstractMapRepresentation,
    second_map::AbstractMapRepresentation;
    priors=nothing,
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    max_input_dimension=32,
    max_output_dimension=32,
    max_choi_dimension=64,
    max_choi_entries=4_096,
    max_dense_entries=1_000_000,
    max_variables=50_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    first = _channel_opt_prepare_map(
        first_map;
        operation="channel_distinguishability",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    second = _channel_opt_prepare_map(
        second_map;
        operation="channel_distinguishability",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    (first.input_dimension, first.output_dimension) ==
    (second.input_dimension, second.output_dimension) || throw(
        DimensionMismatch(
            "channel input/output dimensions must agree; got " *
            "$(first.input_dimension)→$(first.output_dimension) and " *
            "$(second.input_dimension)→$(second.output_dimension)",
        ),
    )
    first_validation = _channel_opt_validate_channel(
        first; allow_densify=allow_densify, operation="channel_distinguishability"
    )
    second_validation = _channel_opt_validate_channel(
        second; allow_densify=allow_densify, operation="channel_distinguishability"
    )
    T = promote_type(first.coefficient_type, second.coefficient_type)
    tolerance = convert(T, max(first.tolerance, second.tolerance))
    prior_values, prior_sum, prior_residual = _channel_opt_priors(priors, T, tolerance)
    baseline = maximum(prior_values)
    validation = (
        first_channel=first_validation,
        second_channel=second_validation,
        prior_sum=prior_sum,
        prior_sum_residual=prior_residual,
    )
    if any(iszero, prior_values)
        value = prior_sum
        return _channel_opt_distinguishability_result(
            ChannelOptimizationAnalyticOptimal,
            T,
            prior_values,
            tolerance;
            success_probability=value,
            lower_bound=value,
            upper_bound=value,
            certified=true,
            certificate_kind=:deterministic_prior,
            residuals=validation,
            certificate=(
                kind=:deterministic_prior, ignored_index=findfirst(iszero, prior_values)
            ),
            message="one channel has zero prior probability, so the other is guessed",
        )
    end
    if first.choi == second.choi
        value = baseline
        return _channel_opt_distinguishability_result(
            ChannelOptimizationAnalyticOptimal,
            T,
            prior_values,
            tolerance;
            success_probability=value,
            lower_bound=value,
            upper_bound=value,
            certified=true,
            certificate_kind=:identical_channels,
            residuals=validation,
            certificate=(kind=:exact_identical_choi,),
            message="the channels have exactly equal Choi matrices, so the more likely prior is optimal",
        )
    end
    difference_matrix = prior_values[1] .* first.choi - prior_values[2] .* second.choi
    difference = ChoiRepresentation(
        difference_matrix,
        OperatorSpace(
            first.input_dimension,
            first.input_dimension,
            first.output_dimension,
            first.output_dimension,
        ),
    )
    norm_result = diamond_norm(
        difference;
        backend=backend,
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
        max_dense_entries=max_dense_entries,
        max_variables=max_variables,
        limits=limits,
    )
    lower = if norm_result.lower_bound === nothing
        baseline
    else
        max(baseline, (prior_sum + norm_result.lower_bound) / 2)
    end
    upper = if norm_result.upper_bound === nothing
        prior_sum
    else
        (prior_sum + norm_result.upper_bound) / 2
    end
    probability =
        norm_result.value === nothing ? nothing : (prior_sum + norm_result.value) / 2
    return _channel_opt_distinguishability_result(
        norm_result.status,
        T,
        prior_values,
        tolerance;
        success_probability=probability,
        lower_bound=lower,
        upper_bound=upper,
        certified=norm_result.certified,
        certificate_kind=norm_result.certified ? :channel_holevo_helstrom : nothing,
        diamond_norm_result=norm_result,
        residuals=validation,
        certificate=(
            kind=:channel_holevo_helstrom,
            weighted_difference_norm_certificate=norm_result.certificate,
        ),
        warnings=norm_result.warnings,
        message="channel Holevo--Helstrom reduction: " * norm_result.message,
    )
end

function _channel_opt_apply_affine(
    density::HermitianAffineMatrix, map::AbstractMapRepresentation; name::Symbol
)
    return _optimization_apply_matrix_map(
        function (coefficient)
            output = apply_channel(coefficient, map)
            ishermitian(output) || throw(
                DomainError(
                    maximum(abs, output - adjoint(output)),
                    "the validated completely positive map did not produce an exactly " *
                    "Hermitian affine coefficient; no symmetrization is applied",
                ),
            )
            return output
        end,
        density,
        name,
    )
end

function _channel_opt_fidelity_objective(
    block::HermitianAffineMatrix{T}, output_dimension::Int
) where {T<:Real}
    variables = Int[]
    coefficients = T[]
    for term in block.terms
        coupling = @view term.coefficient[
            1:output_dimension, (output_dimension + 1):(2output_dimension)
        ]
        value = convert(T, real(tr(coupling)))
        iszero(value) && continue
        push!(variables, term.variable)
        push!(coefficients, value)
    end
    return AffineScalar(zero(T), variables, coefficients, block.variable_count)
end

function _channel_opt_prepare_fidelity_maps(
    first_map,
    second_map;
    atol,
    rtol,
    allow_densify,
    max_input_dimension,
    max_output_dimension,
    max_choi_dimension,
    max_choi_entries,
    require_trace_preserving,
)
    first = _channel_opt_prepare_map(
        first_map;
        operation="maximum_output_fidelity",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    second = _channel_opt_prepare_map(
        second_map;
        operation="maximum_output_fidelity",
        atol=atol,
        rtol=rtol,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
    )
    (first.input_dimension, first.output_dimension) ==
    (second.input_dimension, second.output_dimension) || throw(
        DimensionMismatch(
            "maximum_output_fidelity requires equal input and output dimensions; " *
            "got $(first.input_dimension)→$(first.output_dimension) and " *
            "$(second.input_dimension)→$(second.output_dimension)",
        ),
    )
    if require_trace_preserving
        first_validation = _channel_opt_validate_channel(
            first; allow_densify=allow_densify, operation="maximum_output_fidelity"
        )
        second_validation = _channel_opt_validate_channel(
            second; allow_densify=allow_densify, operation="maximum_output_fidelity"
        )
    else
        first_cp = _channel_opt_cp_diagnostic(first; allow_densify=allow_densify)
        second_cp = _channel_opt_cp_diagnostic(second; allow_densify=allow_densify)
        first_cp !== nothing && first_cp.status === MatrixPredicateSatisfied ||
            throw(DomainError(first_cp, "the first map must be completely positive"))
        second_cp !== nothing && second_cp.status === MatrixPredicateSatisfied ||
            throw(DomainError(second_cp, "the second map must be completely positive"))
        first_validation = (complete_positivity=first_cp,)
        second_validation = (complete_positivity=second_cp,)
    end
    return first, second, first_validation, second_validation
end

function _channel_opt_build_fidelity_problem(
    first,
    second,
    first_validation,
    second_validation;
    max_variables,
    limits::OptimizationLimits,
    require_trace_preserving::Bool,
)
    variable_limit = _channel_opt_positive_integer(max_variables, "max_variables")
    d_in = first.input_dimension
    d_out = first.output_dimension
    variable_count_big = 2BigInt(d_in)^2 + (2BigInt(d_out))^2
    variable_count_big <= variable_limit || throw(
        ArgumentError(
            "maximum-output-fidelity SDP needs $variable_count_big variables, " *
            "exceeding max_variables=$variable_limit",
        ),
    )
    variable_count = Int(variable_count_big)
    T = promote_type(first.coefficient_type, second.coefficient_type)
    first_density = hermitian_variable(
        :fidelity_rho_first,
        d_in;
        first_variable=1,
        variable_count=variable_count,
        coefficient_type=T,
    )
    second_density = hermitian_variable(
        :fidelity_rho_second,
        d_in;
        first_variable=1 + d_in^2,
        variable_count=variable_count,
        coefficient_type=T,
    )
    block = hermitian_variable(
        :fidelity_block,
        2d_out;
        first_variable=1 + 2d_in^2,
        variable_count=variable_count,
        coefficient_type=T,
    )
    first_output = _channel_opt_apply_affine(
        first_density, first.map; name=:fidelity_first_output
    )
    second_output = _channel_opt_apply_affine(
        second_density, second.map; name=:fidelity_second_output
    )
    top = _channel_opt_principal_block(block, 1:d_out; name=:fidelity_block_top)
    bottom = _channel_opt_principal_block(
        block, (d_out + 1):(2d_out); name=:fidelity_block_bottom
    )
    top_relation = _channel_opt_affine_combination(
        :fidelity_top_relation, ((one(T), top), (-one(T), first_output))
    )
    bottom_relation = _channel_opt_affine_combination(
        :fidelity_bottom_relation, ((one(T), bottom), (-one(T), second_output))
    )
    equalities = AffineEquality{T}[
        _channel_opt_trace_one_equality(first_density, :fidelity_first_trace),
        _channel_opt_trace_one_equality(second_density, :fidelity_second_trace),
    ]
    append!(
        equalities, hermitian_equalities(top_relation; name_prefix=:fidelity_top_structure)
    )
    append!(
        equalities,
        hermitian_equalities(bottom_relation; name_prefix=:fidelity_bottom_structure),
    )
    objective = _channel_opt_fidelity_objective(block, d_out)
    feasible = zeros(T, variable_count)
    maximally_mixed = Matrix{Complex{T}}(I, d_in, d_in) / convert(T, d_in)
    first_output_value = Matrix{Complex{T}}(apply_channel(maximally_mixed, first.map))
    second_output_value = Matrix{Complex{T}}(apply_channel(maximally_mixed, second.map))
    ishermitian(first_output_value) && ishermitian(second_output_value) || throw(
        DomainError(
            nothing,
            "validated completely positive maps did not produce exactly Hermitian " *
            "known-feasible outputs; no symmetrization is applied",
        ),
    )
    feasible_block = [
        first_output_value zeros(Complex{T}, d_out, d_out)
        zeros(Complex{T}, d_out, d_out) second_output_value
    ]
    _channel_opt_write_hermitian_coordinates!(feasible, 1, maximally_mixed)
    _channel_opt_write_hermitian_coordinates!(feasible, 1 + d_in^2, maximally_mixed)
    _channel_opt_write_hermitian_coordinates!(feasible, 1 + 2d_in^2, feasible_block)
    program = SemidefiniteProgram(
        :maximum_output_fidelity_primal,
        :maximize,
        variable_count,
        objective;
        equalities=equalities,
        psd_constraints=[first_density, second_density, block],
        primal_views=[first_density, second_density, first_output, second_output, block],
        initial_point=feasible,
        known_feasible_point=feasible,
        limits=limits,
        metadata=(
            formulation=:watrous_maximum_output_fidelity_primal,
            input_dimension=d_in,
            output_dimension=d_out,
            trace_preserving_required=require_trace_preserving,
            first_validation=first_validation,
            second_validation=second_validation,
            source=:watrous_tqi_definition_3_57_and_fidelity_sdp,
        ),
    )
    return MaximumOutputFidelityProblem(
        first.map,
        second.map,
        d_in,
        d_out,
        program,
        (first_density, second_density),
        block,
        max(first.tolerance, second.tolerance),
        require_trace_preserving,
    )
end

"""
    maximum_output_fidelity_problem(first_map, second_map; kwargs...)

Build the solver-neutral SDP

```math
\\begin{aligned}
\\text{maximize}\\quad & \\operatorname{Re}\\operatorname{tr}(X),\\\\
\\text{subject to}\\quad &
\\begin{bmatrix}
\\Phi(\\rho) & X\\\\
X^\\dagger & \\Psi(\\sigma)
\\end{bmatrix}\\succeq0,\\\\
&\\rho,\\sigma\\succeq0,\\quad
\\operatorname{tr}(\\rho)=\\operatorname{tr}(\\sigma)=1.
\\end{aligned}
```

By default, both inputs must be validated quantum channels. Setting
`require_trace_preserving=false` explicitly enables the more general
completely-positive-map definition from Watrous; values then need not lie in
`[0,1]`.
"""
function maximum_output_fidelity_problem(
    first_map::AbstractMapRepresentation,
    second_map::AbstractMapRepresentation;
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    require_trace_preserving::Bool=true,
    max_input_dimension=32,
    max_output_dimension=32,
    max_choi_dimension=64,
    max_choi_entries=4_096,
    max_variables=50_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    first, second, first_validation, second_validation = _channel_opt_prepare_fidelity_maps(
        first_map,
        second_map;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
        require_trace_preserving=require_trace_preserving,
    )
    return _channel_opt_build_fidelity_problem(
        first,
        second,
        first_validation,
        second_validation;
        max_variables=max_variables,
        limits=limits,
        require_trace_preserving=require_trace_preserving,
    )
end

function _channel_opt_replacer_output(prepared)
    d_in = prepared.input_dimension
    d_out = prepared.output_dimension
    choi = prepared.choi
    reference = copy(@view choi[1:d_out, 1:d_out])
    for input_column in 1:d_in, input_row in 1:d_in
        rows = ((input_row - 1) * d_out + 1):(input_row * d_out)
        columns = ((input_column - 1) * d_out + 1):(input_column * d_out)
        block = @view choi[rows, columns]
        expected =
            input_row == input_column ? reference : spzeros(eltype(choi), d_out, d_out)
        block == expected || return nothing
    end
    return reference
end

function _channel_opt_root_fidelity(
    first::AbstractMatrix,
    second::AbstractMatrix,
    tolerance,
    allow_densify::Bool,
    max_dense_entries,
)
    dense_first = _channel_opt_dense_matrix(
        first,
        allow_densify,
        max_dense_entries,
        "maximum-output-fidelity analytic first output",
    )
    dense_second = _channel_opt_dense_matrix(
        second,
        allow_densify,
        max_dense_entries,
        "maximum-output-fidelity analytic second output",
    )
    ishermitian(dense_first) && ishermitian(dense_second) || return nothing
    first_eigen = eigen(Hermitian(dense_first))
    second_eigen = eigen(Hermitian(dense_second))
    minimum_first = minimum(real, first_eigen.values)
    minimum_second = minimum(real, second_eigen.values)
    (minimum_first >= -tolerance && minimum_second >= -tolerance) || return nothing
    any(value -> value < zero(value), first_eigen.values) && return nothing
    any(value -> value < zero(value), second_eigen.values) && return nothing
    root_first =
        first_eigen.vectors *
        Diagonal(sqrt.(first_eigen.values)) *
        adjoint(first_eigen.vectors)
    root_second =
        second_eigen.vectors *
        Diagonal(sqrt.(second_eigen.values)) *
        adjoint(second_eigen.vectors)
    singular_values = svdvals(root_first * root_second)
    return (
        value=sum(singular_values),
        singular_values=singular_values,
        first_spectrum=copy(first_eigen.values),
        second_spectrum=copy(second_eigen.values),
    )
end

function _channel_opt_fidelity_result(
    status,
    ::Type{T},
    input_dimension,
    output_dimension,
    tolerance;
    value=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    certified=false,
    certificate_kind=nothing,
    problem=nothing,
    input_states=nothing,
    output_states=nothing,
    coupling=nothing,
    optimization_result=nothing,
    residuals=NamedTuple(),
    certificate=nothing,
    warnings=(),
    message,
) where {T<:Real}
    return MaximumOutputFidelityResult(
        status,
        value === nothing ? nothing : convert(T, value),
        lower_bound === nothing ? nothing : convert(T, lower_bound),
        upper_bound === nothing ? nothing : convert(T, upper_bound),
        certified,
        certificate_kind,
        input_dimension,
        output_dimension,
        problem,
        input_states,
        output_states,
        coupling,
        optimization_result,
        residuals,
        certificate,
        convert(T, tolerance),
        Tuple(String(warning) for warning in warnings),
        String(message),
    )
end

function _channel_opt_analytic_fidelity(
    first, second; require_trace_preserving::Bool, allow_densify::Bool, max_dense_entries
)
    T = promote_type(first.coefficient_type, second.coefficient_type)
    tolerance = convert(T, max(first.tolerance, second.tolerance))
    d_in = first.input_dimension
    d_out = first.output_dimension
    if require_trace_preserving && first.choi == second.choi
        density = Matrix{Complex{T}}(I, d_in, d_in) / convert(T, d_in)
        output = Matrix{Complex{T}}(apply_channel(density, first.map))
        return (
            value=one(T),
            certificate_kind=:identical_channels,
            input_states=(density, copy(density)),
            output_states=(output, copy(output)),
            certificate=(kind=:exact_identical_choi,),
            residuals=(exact_identical_choi=true,),
            message="identical trace-preserving channels attain root fidelity one",
        )
    end
    if require_trace_preserving
        basis_states = Tuple(
            sparse([index], [index], Complex{T}[one(T)], d_in, d_in) for index in 1:d_in
        )
        first_outputs = Tuple(
            Matrix{Complex{T}}(apply_channel(state, first.map)) for state in basis_states
        )
        second_outputs = Tuple(
            Matrix{Complex{T}}(apply_channel(state, second.map)) for state in basis_states
        )
        for first_index in 1:d_in, second_index in 1:d_in
            first_outputs[first_index] == second_outputs[second_index] || continue
            return (
                value=one(T),
                certificate_kind=:common_basis_output,
                input_states=(
                    Matrix{Complex{T}}(basis_states[first_index]),
                    Matrix{Complex{T}}(basis_states[second_index]),
                ),
                output_states=(first_outputs[first_index], second_outputs[second_index]),
                certificate=(
                    kind=:exact_common_basis_output,
                    first_basis_index=first_index,
                    second_basis_index=second_index,
                ),
                residuals=(exact_common_output=true,),
                message="computational-basis inputs produce exactly equal channel outputs, attaining root fidelity one",
            )
        end
    end
    first_replacer = _channel_opt_replacer_output(first)
    second_replacer = _channel_opt_replacer_output(second)
    if first_replacer !== nothing && second_replacer !== nothing
        fidelity_data = try
            _channel_opt_root_fidelity(
                first_replacer,
                second_replacer,
                tolerance,
                allow_densify,
                max_dense_entries,
            )
        catch error
            error isa ArgumentError || rethrow()
            nothing
        end
        fidelity_data === nothing && return nothing
        density = Matrix{Complex{T}}(I, d_in, d_in) / convert(T, d_in)
        return (
            value=convert(T, fidelity_data.value),
            certificate_kind=:replacer_channels,
            input_states=(density, copy(density)),
            output_states=(
                Matrix{Complex{T}}(first_replacer), Matrix{Complex{T}}(second_replacer)
            ),
            certificate=(
                kind=:exact_replacer_choi_structure,
                output_fidelity_singular_values=fidelity_data.singular_values,
            ),
            residuals=(
                first_output_spectrum=fidelity_data.first_spectrum,
                second_output_spectrum=fidelity_data.second_spectrum,
            ),
            message="both channels are exact replacer channels, so their fixed-output fidelity is optimal",
        )
    end
    return nothing
end

function _channel_opt_fidelity_primal(
    problem::MaximumOutputFidelityProblem,
    optimization::OptimizationResult,
    tolerance,
    allow_densify::Bool,
    max_dense_entries,
)
    optimization.primal === nothing && return nothing
    views = optimization.primal.views
    names = (
        :fidelity_rho_first,
        :fidelity_rho_second,
        :fidelity_first_output,
        :fidelity_second_output,
        :fidelity_block,
    )
    all(name -> hasproperty(views, name), names) || return nothing
    first_density = Matrix(getproperty(views, :fidelity_rho_first))
    second_density = Matrix(getproperty(views, :fidelity_rho_second))
    first_output = Matrix(getproperty(views, :fidelity_first_output))
    second_output = Matrix(getproperty(views, :fidelity_second_output))
    block = Matrix(getproperty(views, :fidelity_block))
    d_out = problem.output_dimension
    coupling = Matrix(@view block[1:d_out, (d_out + 1):(2d_out)])
    first_psd = _channel_opt_spectral_residual(
        first_density,
        allow_densify,
        max_dense_entries,
        "maximum-output-fidelity first-density validation",
    )
    second_psd = _channel_opt_spectral_residual(
        second_density,
        allow_densify,
        max_dense_entries,
        "maximum-output-fidelity second-density validation",
    )
    block_psd = _channel_opt_spectral_residual(
        block, allow_densify, max_dense_entries, "maximum-output-fidelity block validation"
    )
    trace_residual = max(
        abs(real(tr(first_density)) - 1), abs(real(tr(second_density)) - 1)
    )
    top = @view block[1:d_out, 1:d_out]
    bottom = @view block[(d_out + 1):(2d_out), (d_out + 1):(2d_out)]
    structure_residual = max(
        maximum(abs, top - first_output; init=zero(tolerance)),
        maximum(abs, bottom - second_output; init=zero(tolerance)),
    )
    value = real(tr(coupling))
    objective_residual = if optimization.objective_value === nothing
        nothing
    else
        abs(value - optimization.objective_value)
    end
    valid =
        first_psd.hermiticity_residual <= tolerance &&
        second_psd.hermiticity_residual <= tolerance &&
        block_psd.hermiticity_residual <= tolerance &&
        first_psd.positivity_violation <= tolerance &&
        second_psd.positivity_violation <= tolerance &&
        block_psd.positivity_violation <= tolerance &&
        trace_residual <= tolerance &&
        structure_residual <= tolerance &&
        (objective_residual === nothing || objective_residual <= 8tolerance)
    return (
        valid=valid,
        value=value,
        input_states=(first_density, second_density),
        output_states=(first_output, second_output),
        coupling=coupling,
        residuals=(
            first_density=first_psd,
            second_density=second_psd,
            block=block_psd,
            trace_residual=trace_residual,
            structure_residual=structure_residual,
            solver_objective_residual=objective_residual,
        ),
    )
end

function _channel_opt_solver_fidelity_result(
    problem, optimization, backend; allow_densify, max_dense_entries
)
    T = typeof(problem.tolerance)
    tolerance = _channel_opt_backend_tolerance(problem, optimization, backend)
    primal = try
        _channel_opt_fidelity_primal(
            problem,
            optimization,
            8tolerance,
            allow_densify || (backend isa JuMPBackend && backend.allow_densify),
            max_dense_entries,
        )
    catch error
        error isa ArgumentError || rethrow()
        nothing
    end
    lower = if primal !== nothing && primal.valid
        max(zero(T), convert(T, primal.value))
    else
        zero(T)
    end
    solver_upper = if optimization.objective_bound === nothing
        nothing
    else
        convert(T, optimization.objective_bound)
    end
    upper = if problem.trace_preserving_required
        if solver_upper === nothing || solver_upper > one(T)
            one(T)
        else
            solver_upper
        end
    else
        solver_upper
    end
    bound_order_valid = upper === nothing || upper + 8tolerance >= lower
    matched =
        upper !== nothing &&
        bound_order_valid &&
        abs(lower - upper) <= 16tolerance * max(one(T), abs(lower), abs(upper))
    value = matched ? (lower + upper) / 2 : nothing
    residuals = if primal === nothing
        (primal_available=false,)
    else
        merge((primal_available=true,), primal.residuals)
    end
    input_states = primal === nothing || !primal.valid ? nothing : primal.input_states
    output_states = primal === nothing || !primal.valid ? nothing : primal.output_states
    coupling = primal === nothing || !primal.valid ? nothing : primal.coupling
    if optimization.status === OptimizationBackendUnavailable
        return _channel_opt_fidelity_result(
            ChannelOptimizationBackendUnavailable,
            T,
            problem.input_dimension,
            problem.output_dimension,
            tolerance;
            lower_bound=lower,
            upper_bound=upper,
            problem=problem,
            optimization_result=optimization,
            residuals=residuals,
            message=optimization.message,
        )
    elseif optimization.status === OptimizationLimit
        return _channel_opt_fidelity_result(
            ChannelOptimizationResourceLimit,
            T,
            problem.input_dimension,
            problem.output_dimension,
            tolerance;
            lower_bound=lower,
            upper_bound=bound_order_valid ? upper : nothing,
            problem=problem,
            input_states=input_states,
            output_states=output_states,
            coupling=coupling,
            optimization_result=optimization,
            residuals=residuals,
            message="the optimizer stopped at a resource limit; validated bounds are retained",
        )
    elseif optimization.status === OptimizationOptimal &&
        primal !== nothing &&
        primal.valid &&
        matched
        return _channel_opt_fidelity_result(
            ChannelOptimizationSolverOptimal,
            T,
            problem.input_dimension,
            problem.output_dimension,
            tolerance;
            value=value,
            lower_bound=lower,
            upper_bound=upper,
            certified=false,
            problem=problem,
            input_states=input_states,
            output_states=output_states,
            coupling=coupling,
            optimization_result=optimization,
            residuals=residuals,
            certificate=(
                kind=:numerical_primal_dual_agreement,
                primal_residual=optimization.primal_residual,
                dual_residual=optimization.dual_residual,
                absolute_gap=optimization.absolute_gap,
            ),
            message="the optimizer returned residual-checked matching fidelity bounds",
        )
    elseif optimization.status in (OptimizationOptimal, OptimizationFeasible) &&
        primal !== nothing &&
        primal.valid
        status = if optimization.status === OptimizationFeasible
            ChannelOptimizationSolverFeasible
        else
            ChannelOptimizationInvalidCertificate
        end
        return _channel_opt_fidelity_result(
            status,
            T,
            problem.input_dimension,
            problem.output_dimension,
            tolerance;
            lower_bound=lower,
            upper_bound=bound_order_valid ? upper : nothing,
            problem=problem,
            input_states=input_states,
            output_states=output_states,
            coupling=coupling,
            optimization_result=optimization,
            residuals=residuals,
            warnings=(
                "the available primal and upper bounds do not establish a common value",
            ),
            message="a residual-checked fidelity witness is available, but optimality is not established",
        )
    end
    failure = optimization.status in (
        OptimizationMalformedBackend,
        OptimizationNumericalFailure,
        OptimizationUnsupported,
        OptimizationUnknown,
        OptimizationInfeasible,
        OptimizationUnbounded,
        OptimizationInconsistent,
    )
    return _channel_opt_fidelity_result(
        failure ? ChannelOptimizationBackendFailure : ChannelOptimizationInvalidCertificate,
        T,
        problem.input_dimension,
        problem.output_dimension,
        tolerance;
        lower_bound=lower,
        upper_bound=bound_order_valid ? upper : nothing,
        problem=problem,
        input_states=input_states,
        output_states=output_states,
        coupling=coupling,
        optimization_result=optimization,
        residuals=residuals,
        message="the optimization outcome does not establish a valid maximum output fidelity",
    )
end

"""
    maximum_output_fidelity(first_map, second_map;
                            backend=NoOptimizationBackend(), kwargs...)

Compute or bound
`max fidelity(first_map(rho), second_map(sigma))` over independent input
density operators. Identical channels and exact replacer channels use
solver-free certificates; other inputs use
[`maximum_output_fidelity_problem`](@ref).

Unlike the pinned complementary-map implementation, this direct fidelity SDP
does not truncate unequal Kraus ranks and does not depend on a chosen Kraus
representation. Inputs are never normalized, padded, symmetrized, or clipped.
"""
function maximum_output_fidelity(
    first_map::AbstractMapRepresentation,
    second_map::AbstractMapRepresentation;
    backend::AbstractOptimizationBackend=NoOptimizationBackend(),
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    require_trace_preserving::Bool=true,
    max_input_dimension=32,
    max_output_dimension=32,
    max_choi_dimension=64,
    max_choi_entries=4_096,
    max_dense_entries=1_000_000,
    max_variables=50_000,
    limits::OptimizationLimits=OptimizationLimits(),
)
    first, second, first_validation, second_validation = _channel_opt_prepare_fidelity_maps(
        first_map,
        second_map;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        max_input_dimension=max_input_dimension,
        max_output_dimension=max_output_dimension,
        max_choi_dimension=max_choi_dimension,
        max_choi_entries=max_choi_entries,
        require_trace_preserving=require_trace_preserving,
    )
    analytic = _channel_opt_analytic_fidelity(
        first,
        second;
        require_trace_preserving=require_trace_preserving,
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
    )
    T = promote_type(first.coefficient_type, second.coefficient_type)
    tolerance = convert(T, max(first.tolerance, second.tolerance))
    if analytic !== nothing
        return _channel_opt_fidelity_result(
            ChannelOptimizationAnalyticOptimal,
            T,
            first.input_dimension,
            first.output_dimension,
            tolerance;
            value=analytic.value,
            lower_bound=analytic.value,
            upper_bound=analytic.value,
            certified=true,
            certificate_kind=analytic.certificate_kind,
            input_states=analytic.input_states,
            output_states=analytic.output_states,
            residuals=analytic.residuals,
            certificate=analytic.certificate,
            message=analytic.message,
        )
    end
    problem = _channel_opt_build_fidelity_problem(
        first,
        second,
        first_validation,
        second_validation;
        max_variables=max_variables,
        limits=limits,
        require_trace_preserving=require_trace_preserving,
    )
    optimization = solve_optimization(problem.program, backend)
    return _channel_opt_solver_fidelity_result(
        problem,
        optimization,
        backend;
        allow_densify=allow_densify,
        max_dense_entries=max_dense_entries,
    )
end
