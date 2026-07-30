# Central solver-neutral affine maps and residual calculations for the optional
# optimization layer.

"""
    hermitian_variable(name, dimension; first_variable=1, variable_count=...)

Construct a sparse affine complex-Hermitian matrix using `dimension^2` real
coordinates. Coordinates are ordered as diagonal entries, real upper-triangle
entries, then imaginary upper-triangle entries. The return value is a
[`HermitianAffineMatrix`](@ref).
"""
function hermitian_variable(
    name::Symbol,
    dimension::Integer;
    first_variable::Integer=1,
    variable_count::Integer=BigInt(first_variable) + BigInt(dimension)^2 - 1,
    coefficient_type::Type{T}=Float64,
) where {T<:Real}
    d = _optimization_positive_int(dimension, "dimension")
    first = _optimization_positive_int(first_variable, "first_variable")
    count = _optimization_positive_int(variable_count, "variable_count")
    needed = BigInt(first) + BigInt(d)^2 - 1
    needed <= count || throw(
        DimensionMismatch(
            "Hermitian variable needs coordinates $first:$needed but variable_count=$count",
        ),
    )
    isconcretetype(T) ||
        throw(ArgumentError("coefficient_type must be a concrete real type"))
    zero_matrix = spzeros(Complex{T}, d, d)
    variables = Int[]
    coefficients = SparseMatrixCSC{Complex{T},Int}[]
    next_variable = first
    for index in 1:d
        push!(variables, next_variable)
        push!(coefficients, sparse([index], [index], Complex{T}[one(T)], d, d))
        next_variable += 1
    end
    for column in 2:d, row in 1:(column - 1)
        push!(variables, next_variable)
        push!(
            coefficients,
            sparse([row, column], [column, row], Complex{T}[one(T), one(T)], d, d),
        )
        next_variable += 1
    end
    for column in 2:d, row in 1:(column - 1)
        push!(variables, next_variable)
        push!(
            coefficients,
            sparse(
                [row, column],
                [column, row],
                Complex{T}[complex(zero(T), one(T)), complex(zero(T), -one(T))],
                d,
                d,
            ),
        )
        next_variable += 1
    end
    next_variable == first + d^2 ||
        error("internal Hermitian coordinate construction is inconsistent")
    return HermitianAffineMatrix(name, zero_matrix, variables, coefficients, count)
end

"""
    evaluate_affine(function_or_matrix, coordinates)

Evaluate solver-neutral scalar or Hermitian affine data. The matrix method
preserves sparse storage.
"""
function evaluate_affine(
    function_data::AffineScalar{T}, coordinates::AbstractVector{<:Real}
) where {T<:Real}
    firstindex(coordinates) == 1 ||
        throw(ArgumentError("coordinates must use one-based indexing"))
    length(coordinates) == length(function_data.coefficients) || throw(
        DimensionMismatch(
            "expected $(length(function_data.coefficients)) coordinates, got " *
            "$(length(coordinates))",
        ),
    )
    value = function_data.constant
    indices, values = findnz(function_data.coefficients)
    for (index, coefficient) in zip(indices, values)
        _optimization_validate_real(coordinates[index], "coordinates")
        value += coefficient * coordinates[index]
    end
    return value
end

function evaluate_affine(
    matrix::HermitianAffineMatrix{T}, coordinates::AbstractVector{<:Real}
) where {T<:Real}
    firstindex(coordinates) == 1 ||
        throw(ArgumentError("coordinates must use one-based indexing"))
    length(coordinates) == matrix.variable_count || throw(
        DimensionMismatch(
            "expected $(matrix.variable_count) coordinates, got $(length(coordinates))"
        ),
    )
    R = promote_type(T, eltype(coordinates))
    result = sparse(Complex{R}.(matrix.constant))
    for term in matrix.terms
        value = coordinates[term.variable]
        _optimization_validate_real(value, "coordinates")
        iszero(value) && continue
        result += value * term.coefficient
    end
    return result
end

"""
    trace_affine(matrix)

Return the real scalar affine trace of a Hermitian affine matrix.
"""
function trace_affine(matrix::HermitianAffineMatrix{T}) where {T<:Real}
    indices = Int[]
    values = T[]
    for term in matrix.terms
        value = real(tr(term.coefficient))
        iszero(value) && continue
        push!(indices, term.variable)
        push!(values, value)
    end
    return AffineScalar(real(tr(matrix.constant)), indices, values, matrix.variable_count)
end

function _optimization_apply_matrix_map(
    map_function, matrix::HermitianAffineMatrix{T}, name::Symbol
) where {T<:Real}
    constant = map_function(matrix.constant)
    variables = Int[]
    coefficients = AbstractMatrix{<:Number}[]
    for term in matrix.terms
        coefficient = map_function(term.coefficient)
        nnz(sparse(coefficient)) == 0 && continue
        push!(variables, term.variable)
        push!(coefficients, coefficient)
    end
    return HermitianAffineMatrix(
        name, constant, variables, coefficients, matrix.variable_count
    )
end

"""
    partial_trace_affine(matrix, dims; trace_out, name=...)

Apply the package subsystem convention coefficientwise to an affine Hermitian
matrix. Sparse coefficients remain sparse.
"""
function partial_trace_affine(
    matrix::HermitianAffineMatrix,
    dims;
    trace_out,
    name::Symbol=Symbol(matrix.name, :_partial_trace),
)
    plan = PartialTracePlan(dims, trace_out)
    plan.layout.total_dimension == matrix.dimension || throw(
        DimensionMismatch(
            "subsystem dimension $(plan.layout.total_dimension) does not match " *
            "affine matrix dimension $(matrix.dimension)",
        ),
    )
    return _optimization_apply_matrix_map(
        coefficient -> partial_trace(coefficient, plan), matrix, name
    )
end

"""
    partial_transpose_affine(matrix, dims; systems, name=...)

Apply partial transposition coefficientwise using the package subsystem order.
"""
function partial_transpose_affine(
    matrix::HermitianAffineMatrix,
    dims;
    systems,
    name::Symbol=Symbol(matrix.name, :_partial_transpose),
)
    plan = PartialTransposePlan(dims, systems)
    plan.row_layout.total_dimension == matrix.dimension || throw(
        DimensionMismatch(
            "subsystem dimension $(plan.row_layout.total_dimension) does not match " *
            "affine matrix dimension $(matrix.dimension)",
        ),
    )
    return _optimization_apply_matrix_map(
        coefficient -> partial_transpose(coefficient, plan), matrix, name
    )
end

"""
    tensor_affine(matrix; left=nothing, right=nothing, name=...)

Tensor fixed numeric factors to the left and/or right of every coefficient.
At least one factor is required. Sparse factors and coefficients remain sparse.
"""
function tensor_affine(
    matrix::HermitianAffineMatrix;
    left=nothing,
    right=nothing,
    name::Symbol=Symbol(matrix.name, :_tensor),
)
    left === nothing &&
        right === nothing &&
        throw(ArgumentError("tensor_affine needs a left or right factor"))
    for (factor_name, factor) in (("left", left), ("right", right))
        factor === nothing && continue
        factor isa AbstractMatrix{<:Number} ||
            throw(ArgumentError("$factor_name factor must be a numeric matrix"))
        _optimization_check_hermitian(factor, "$factor_name factor")
    end
    map_function = function (coefficient)
        value = coefficient
        left === nothing || (value = tensor_product(left, value))
        right === nothing || (value = tensor_product(value, right))
        return value
    end
    return _optimization_apply_matrix_map(map_function, matrix, name)
end

"""
    choi_trace_preserving_affine(choi, input_dimension, output_dimension)

Return `Tr_output(choi) - I_input` as an affine Hermitian matrix. The Choi
ordering is output ⊗ input, matching the package's channel convention.
Setting this map equal to zero imposes trace preservation.
"""
function choi_trace_preserving_affine(
    choi::HermitianAffineMatrix,
    input_dimension::Integer,
    output_dimension::Integer;
    name::Symbol=Symbol(choi.name, :_trace_preserving),
)
    d_in = _optimization_positive_int(input_dimension, "input_dimension")
    d_out = _optimization_positive_int(output_dimension, "output_dimension")
    BigInt(d_in) * BigInt(d_out) == choi.dimension || throw(
        DimensionMismatch(
            "Choi dimension $(choi.dimension) does not equal " *
            "output_dimension * input_dimension = $(BigInt(d_out) * d_in)",
        ),
    )
    reduced = partial_trace_affine(choi, (d_out, d_in); trace_out=(1,), name=name)
    shifted_constant =
        reduced.constant -
        sparse(1:d_in, 1:d_in, fill(one(eltype(reduced.constant)), d_in), d_in, d_in)
    return HermitianAffineMatrix(
        name,
        shifted_constant,
        [term.variable for term in reduced.terms],
        [term.coefficient for term in reduced.terms],
        reduced.variable_count,
    )
end

"""
    hermitian_equalities(matrix; target=zeros(...), name_prefix=...)

Convert a Hermitian affine matrix equality into the independent real scalar
equalities: diagonal, real upper triangle, and imaginary upper triangle.
"""
function hermitian_equalities(
    matrix::HermitianAffineMatrix{T};
    target=spzeros(Complex{T}, matrix.dimension, matrix.dimension),
    name_prefix::Symbol=matrix.name,
) where {T<:Real}
    _optimization_check_hermitian(target, "target")
    size(target) == (matrix.dimension, matrix.dimension) ||
        throw(DimensionMismatch("target must match the affine matrix dimension"))
    d = matrix.dimension
    target_converted = sparse(Complex{T}.(target))
    coordinate_specs = Tuple{Int,Int,Symbol}[]
    append!(coordinate_specs, ((index, index, :real) for index in 1:d))
    append!(
        coordinate_specs, ((row, column, :real) for column in 2:d for row in 1:(column - 1))
    )
    append!(
        coordinate_specs, ((row, column, :imag) for column in 2:d for row in 1:(column - 1))
    )
    constraints = AffineEquality{T}[]
    for (row, column, component) in coordinate_specs
        constant_entry = matrix.constant[row, column] - target_converted[row, column]
        constant = component === :real ? real(constant_entry) : imag(constant_entry)
        variables = Int[]
        coefficients = T[]
        for term in matrix.terms
            entry = term.coefficient[row, column]
            value = component === :real ? real(entry) : imag(entry)
            iszero(value) && continue
            push!(variables, term.variable)
            push!(coefficients, value)
        end
        function_data = AffineScalar(
            constant, variables, coefficients, matrix.variable_count
        )
        constraint_name = Symbol(name_prefix, :_, component, :_, row, :_, column)
        push!(constraints, AffineEquality(function_data, constraint_name))
    end
    return constraints
end

"""
    real_block_embedding(matrix)

For Hermitian `H = A + im*B`, return the real symmetric matrix
`[A -B; B A]`. Positive semidefiniteness is equivalent in both
representations. No Hermitian repair is applied.
"""
function real_block_embedding(matrix::AbstractMatrix{<:Number})
    _optimization_check_hermitian(matrix, "matrix")
    A = real.(matrix)
    B = imag.(matrix)
    return [A -B; B A]
end

function real_block_embedding(
    matrix::HermitianAffineMatrix{T}; name::Symbol=Symbol(matrix.name, :_real_block)
) where {T<:Real}
    constant = sparse(real_block_embedding(matrix.constant))
    return HermitianAffineMatrix(
        name,
        constant,
        [term.variable for term in matrix.terms],
        [sparse(real_block_embedding(term.coefficient)) for term in matrix.terms],
        matrix.variable_count,
    )
end

"""
    hermitian_dual_from_real_block(matrix)

Apply the adjoint of [`real_block_embedding`](@ref) to a real symmetric dual
matrix. If `Y` is PSD, the returned complex Hermitian matrix is PSD.
"""
function hermitian_dual_from_real_block(matrix::AbstractMatrix{<:Real})
    _optimization_check_hermitian(matrix, "real block dual")
    iseven(size(matrix, 1)) ||
        throw(DimensionMismatch("real block dual dimension must be even"))
    d = size(matrix, 1) ÷ 2
    Y11 = view(matrix, 1:d, 1:d)
    Y12 = view(matrix, 1:d, (d + 1):(2d))
    Y21 = view(matrix, (d + 1):(2d), 1:d)
    Y22 = view(matrix, (d + 1):(2d), (d + 1):(2d))
    return Hermitian((Y11 + Y22) + im * (Y21 - Y12))
end

function _optimization_spectral_violation(
    matrix::AbstractMatrix{<:Number}; allow_densify::Bool
)
    isempty(matrix) && return zero(real(eltype(matrix)))
    if issparse(matrix) && !allow_densify
        throw(
            ArgumentError(
                "PSD residual evaluation would densify sparse data; pass allow_densify=true"
            ),
        )
    end
    values = eigvals(Hermitian(Matrix(matrix)))
    minimum_value = minimum(real, values)
    return max(zero(minimum_value), -minimum_value)
end

"""
    primal_residual(problem, coordinates; allow_densify=false)

Maximum equality, interval, and PSD violation. PSD spectral work on sparse
model data requires `allow_densify=true`.
"""
function primal_residual(
    problem::SemidefiniteProgram{T},
    coordinates::AbstractVector{<:Real};
    allow_densify::Bool=false,
) where {T<:Real}
    firstindex(coordinates) == 1 ||
        throw(ArgumentError("coordinates must use one-based indexing"))
    length(coordinates) == problem.variable_count ||
        throw(DimensionMismatch("coordinate vector has the wrong length"))
    for value in coordinates
        _optimization_validate_real(value, "coordinates")
    end
    R = promote_type(T, eltype(coordinates))
    residual = zero(R)
    for constraint in problem.equalities
        residual = max(
            residual, abs(evaluate_affine(constraint.function_data, coordinates))
        )
    end
    for constraint in problem.intervals
        value = evaluate_affine(constraint.function_data, coordinates)
        constraint.lower === nothing ||
            (residual = max(residual, max(zero(R), constraint.lower - value)))
        constraint.upper === nothing ||
            (residual = max(residual, max(zero(R), value - constraint.upper)))
    end
    for matrix in problem.psd_constraints
        residual = max(
            residual,
            _optimization_spectral_violation(
                evaluate_affine(matrix, coordinates); allow_densify
            ),
        )
    end
    return residual
end

"""
    dual_stationarity_residual(problem, equality_duals, interval_duals,
                               psd_real_block_duals)

Compute the MOI/JuMP conic stationarity residual. PSD duals are the real-block
matrices returned by JuMP. Feasibility-sense models have a zero objective.
"""
function dual_stationarity_residual(
    problem::SemidefiniteProgram{T},
    equality_duals::AbstractVector{<:Real},
    interval_duals::AbstractVector{<:Real},
    psd_real_block_duals::AbstractVector{<:AbstractMatrix{<:Real}},
) where {T<:Real}
    length(equality_duals) == length(problem.equalities) ||
        throw(DimensionMismatch("wrong number of equality duals"))
    length(interval_duals) == length(problem.intervals) ||
        throw(DimensionMismatch("wrong number of interval duals"))
    length(psd_real_block_duals) == length(problem.psd_constraints) ||
        throw(DimensionMismatch("wrong number of PSD duals"))
    R = promote_type(
        T,
        eltype(equality_duals),
        eltype(interval_duals),
        isempty(psd_real_block_duals) ? T : eltype(first(psd_real_block_duals)),
    )
    gradient = zeros(R, problem.variable_count)
    if problem.sense !== :feasibility
        indices, values = findnz(problem.objective.coefficients)
        gradient[indices] .= values
    end
    contribution = zeros(R, problem.variable_count)
    cone_violation = zero(R)
    for (constraint, dual) in zip(problem.equalities, equality_duals)
        _optimization_validate_real(dual, "equality_duals")
        indices, values = findnz(constraint.function_data.coefficients)
        contribution[indices] .+= dual .* values
    end
    for (constraint, dual) in zip(problem.intervals, interval_duals)
        _optimization_validate_real(dual, "interval_duals")
        indices, values = findnz(constraint.function_data.coefficients)
        contribution[indices] .+= dual .* values
        if constraint.lower !== nothing && constraint.upper === nothing
            cone_violation = max(cone_violation, max(zero(R), -dual))
        elseif constraint.lower === nothing && constraint.upper !== nothing
            cone_violation = max(cone_violation, max(zero(R), dual))
        end
    end
    for (matrix, dual) in zip(problem.psd_constraints, psd_real_block_duals)
        size(dual) == (2matrix.dimension, 2matrix.dimension) ||
            throw(DimensionMismatch("PSD dual real block has the wrong size"))
        for value in dual
            _optimization_validate_real(value, "psd_real_block_duals")
        end
        cone_violation = max(
            cone_violation,
            maximum(abs, dual - transpose(dual); init=zero(R)),
            _optimization_spectral_violation(dual; allow_densify=true),
        )
        for term in matrix.terms
            contribution[term.variable] += real(
                dot(dual, real_block_embedding(term.coefficient))
            )
        end
    end
    residual_vector = if problem.sense === :maximize
        gradient + contribution
    else
        gradient - contribution
    end
    return max(cone_violation, maximum(abs, residual_vector; init=zero(R)))
end

function _optimization_residual_tolerance(
    problem::SemidefiniteProgram, backend::JuMPBackend, objective_value
)
    scale = one(
        promote_type(
            eltype(problem.objective.coefficients),
            typeof(backend.atol),
            typeof(backend.rtol),
        ),
    )
    objective_value === nothing || (scale = max(scale, abs(objective_value)))
    return backend.atol + backend.rtol * scale
end

function _optimization_consistent_status(
    proposed::OptimizationStatus, primal_residual_value, dual_residual_value, tolerance
)
    proposed in (OptimizationOptimal, OptimizationFeasible) || return proposed
    primal_residual_value !== nothing &&
        primal_residual_value > tolerance &&
        return OptimizationInconsistent
    dual_residual_value !== nothing &&
        dual_residual_value > tolerance &&
        return OptimizationInconsistent
    return proposed
end
