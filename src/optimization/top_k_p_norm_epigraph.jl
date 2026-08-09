# Source-informed independent Julia implementation based on the CVX-expression
# contract of QETLAB kpNorm.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# This independently specified solver-neutral replacement leaves numeric
# vectors and matrices on top_k_p_norm. This file defines a composable affine
# epigraph atom; its JuMP/MOI materialization lives in the optional extension.

"""
    ComplexAffineTerm(variable, coefficient)

One arbitrary complex sparse coefficient matrix multiplying a real model
coordinate.
"""
struct ComplexAffineTerm{T<:Real}
    variable::Int
    coefficient::SparseMatrixCSC{Complex{T},Int}
end

"""
    ComplexAffineMatrix(name, constant, variables, coefficients, variable_count)

Package-owned rectangular complex affine matrix
`constant + sum(coefficients[j] * x[variables[j]])`.

Coordinates are real, while constants and coefficients may be complex. Inputs
are copied into sparse storage and must be finite and one-based. Unlike
[`HermitianAffineMatrix`](@ref), no square or Hermiticity condition is imposed.
"""
struct ComplexAffineMatrix{T<:Real}
    name::Symbol
    row_dimension::Int
    column_dimension::Int
    variable_count::Int
    constant::SparseMatrixCSC{Complex{T},Int}
    terms::Vector{ComplexAffineTerm{T}}
end

function _complex_affine_check_matrix(matrix, name::AbstractString)
    matrix isa AbstractMatrix{<:Number} ||
        throw(ArgumentError("$name must be a numeric matrix"))
    firstindex(matrix, 1) == 1 && firstindex(matrix, 2) == 1 ||
        throw(ArgumentError("$name must use one-based indexing"))
    isempty(matrix) && throw(ArgumentError("$name must not be empty"))
    for value in matrix
        value isa Number && !(value isa Bool) ||
            throw(ArgumentError("$name must contain numeric non-Boolean values"))
        isfinite(value) || throw(ArgumentError("$name must contain only finite values"))
    end
    return nothing
end

function ComplexAffineMatrix(
    name::Symbol,
    constant::AbstractMatrix{<:Number},
    variables::AbstractVector{<:Integer},
    coefficients::AbstractVector{<:AbstractMatrix{<:Number}},
    variable_count::Integer,
)
    count = _optimization_positive_int(variable_count, "variable_count")
    _complex_affine_check_matrix(constant, "constant")
    firstindex(variables) == 1 && firstindex(coefficients) == 1 ||
        throw(ArgumentError("affine matrix inputs must use one-based indexing"))
    length(variables) == length(coefficients) ||
        throw(DimensionMismatch("variables and coefficients must have equal length"))
    shape = size(constant)
    real_types = Type[_optimization_matrix_real_type(constant)]
    checked_variables = Int[]
    for (position, (variable, coefficient)) in enumerate(zip(variables, coefficients))
        variable isa Bool &&
            throw(ArgumentError("variable indices must be one-based integers"))
        1 <= variable <= count ||
            throw(ArgumentError("variable index $variable is outside 1:$count"))
        _complex_affine_check_matrix(coefficient, "coefficient $position")
        size(coefficient) == shape || throw(
            DimensionMismatch(
                "coefficient $position has size $(size(coefficient)); expected $shape"
            ),
        )
        push!(checked_variables, Int(variable))
        push!(real_types, _optimization_matrix_real_type(coefficient))
    end
    length(unique(checked_variables)) == length(checked_variables) ||
        throw(ArgumentError("variables must not contain duplicates"))
    T = promote_type(real_types...)
    isconcretetype(T) && T <: Real ||
        throw(ArgumentError("affine matrix entries need a concrete real component type"))
    terms = ComplexAffineTerm{T}[]
    for (variable, coefficient) in zip(checked_variables, coefficients)
        converted = sparse(Complex{T}.(coefficient))
        nnz(converted) == 0 && continue
        push!(terms, ComplexAffineTerm{T}(variable, converted))
    end
    return ComplexAffineMatrix{T}(
        name, shape[1], shape[2], count, sparse(Complex{T}.(constant)), terms
    )
end

"""
    complex_affine_variable(name, rows, columns; kwargs...)

Construct a rectangular complex affine variable. Real entry coordinates come
first in column-major order, followed by imaginary entry coordinates. A
`rows`-by-`columns` variable therefore uses `2 * rows * columns` real
coordinates.
"""
function complex_affine_variable(
    name::Symbol,
    rows::Integer,
    columns::Integer;
    first_variable::Integer=1,
    variable_count::Integer=BigInt(first_variable) +
                            BigInt(2) * BigInt(rows) * BigInt(columns) - 1,
    coefficient_type::Type{T}=Float64,
) where {T<:Real}
    m = _optimization_positive_int(rows, "rows")
    n = _optimization_positive_int(columns, "columns")
    first = _optimization_positive_int(first_variable, "first_variable")
    count = _optimization_positive_int(variable_count, "variable_count")
    needed = BigInt(first) + BigInt(2) * m * n - 1
    needed <= count || throw(
        DimensionMismatch(
            "complex affine variable needs coordinates $first:$needed but " *
            "variable_count=$count",
        ),
    )
    isconcretetype(T) ||
        throw(ArgumentError("coefficient_type must be a concrete real type"))
    variables = Int[]
    coefficients = SparseMatrixCSC{Complex{T},Int}[]
    next_variable = first
    for column in 1:n, row in 1:m
        push!(variables, next_variable)
        push!(coefficients, sparse([row], [column], Complex{T}[one(T)], m, n))
        next_variable += 1
    end
    for column in 1:n, row in 1:m
        push!(variables, next_variable)
        push!(
            coefficients,
            sparse([row], [column], Complex{T}[complex(zero(T), one(T))], m, n),
        )
        next_variable += 1
    end
    return ComplexAffineMatrix(
        name, spzeros(Complex{T}, m, n), variables, coefficients, count
    )
end

"""
    evaluate_affine(matrix::ComplexAffineMatrix, coordinates)

Evaluate a rectangular complex affine matrix without densifying its sparse
coefficient storage.
"""
function evaluate_affine(
    matrix::ComplexAffineMatrix{T}, coordinates::AbstractVector{<:Real}
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

function _top_k_p_vector_diagonal(matrix::ComplexAffineMatrix{T}) where {T<:Real}
    length_vector = max(matrix.row_dimension, matrix.column_dimension)
    constant_values = vec(matrix.constant)
    constant = spdiagm(0 => constant_values)
    variables = Int[]
    coefficients = SparseMatrixCSC{Complex{T},Int}[]
    for term in matrix.terms
        push!(variables, term.variable)
        push!(coefficients, spdiagm(0 => vec(term.coefficient)))
    end
    return ComplexAffineMatrix(
        Symbol(matrix.name, :_diagonal),
        constant,
        variables,
        coefficients,
        matrix.variable_count,
    )
end

function _top_k_p_hermitian_dilation(
    matrix::ComplexAffineMatrix{T}, name::Symbol
) where {T<:Real}
    zero_rows = spzeros(Complex{T}, matrix.row_dimension, matrix.row_dimension)
    zero_columns = spzeros(Complex{T}, matrix.column_dimension, matrix.column_dimension)
    constant = sparse([zero_rows matrix.constant; adjoint(matrix.constant) zero_columns])
    coefficients = SparseMatrixCSC{Complex{T},Int}[]
    for term in matrix.terms
        push!(
            coefficients,
            sparse([
                zero_rows term.coefficient
                adjoint(term.coefficient) zero_columns
            ]),
        )
    end
    return HermitianAffineMatrix(
        name,
        constant,
        [term.variable for term in matrix.terms],
        coefficients,
        matrix.variable_count,
    )
end

"""
    TopKPNormEpigraph

Solver-neutral epigraph atom for the top-`k`, Schatten/entrywise `p`-norm of a
rectangular complex affine matrix.

For a matrix, the atom acts on singular values. A one-row or one-column input
is converted to a diagonal affine matrix so it acts on entry magnitudes,
matching the documented vector contract. The optional JuMP extension
materializes exact Ky Fan semidefinite constraints and an MOI `NormCone`.
"""
struct TopKPNormEpigraph{T<:Real,P<:Real}
    matrix::ComplexAffineMatrix{T}
    dilation::HermitianAffineMatrix{T}
    k::Int
    p::P
    singular_value_count::Int
    auxiliary_variable_count::Int
    psd_block_count::Int
    limits::OptimizationLimits
end

function _top_k_p_epigraph_order(value)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("p must be a real number in [1, Inf]"))
    isnan(value) && throw(ArgumentError("p must not be NaN"))
    value >= one(value) || throw(ArgumentError("p must lie in [1, Inf]; got $value"))
    return value
end

function _top_k_p_epigraph_preflight(
    matrix::ComplexAffineMatrix,
    dilation::HermitianAffineMatrix,
    singular_value_count::Int,
    limits::OptimizationLimits,
)
    n = BigInt(singular_value_count)
    d = BigInt(dilation.dimension)
    auxiliary_variables = n * (d^2 + 1) + n + 1
    total_variables = BigInt(matrix.variable_count) + auxiliary_variables
    total_variables <= limits.max_variables || throw(
        ArgumentError(
            "top-k p-norm epigraph needs $total_variables total variables, " *
            "exceeding max_variables=$(limits.max_variables)",
        ),
    )
    psd_blocks = BigInt(2) * n
    psd_blocks <= limits.max_psd_blocks || throw(
        ArgumentError(
            "top-k p-norm epigraph needs $psd_blocks PSD blocks, exceeding " *
            "max_psd_blocks=$(limits.max_psd_blocks)",
        ),
    )
    d <= limits.max_psd_dimension || throw(
        ArgumentError(
            "top-k p-norm auxiliary blocks have dimension $d, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    real_dimension = BigInt(2) * d
    real_triangle_entries = psd_blocks * real_dimension * (real_dimension + 1) ÷ 2
    real_triangle_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "top-k p-norm real PSD blocks need $real_triangle_entries triangle " *
            "entries, exceeding max_model_entries=$(limits.max_model_entries)",
        ),
    )
    stored_input_entries = BigInt(nnz(dilation.constant))
    for term in dilation.terms
        stored_input_entries += nnz(term.coefficient)
    end
    stored_input_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "top-k p-norm affine dilation stores $stored_input_entries entries, " *
            "exceeding max_model_entries=$(limits.max_model_entries)",
        ),
    )
    return Int(auxiliary_variables), Int(psd_blocks)
end

"""
    top_k_p_norm_epigraph(matrix::ComplexAffineMatrix, k, p; kwargs...)

Construct a composable, solver-neutral epigraph atom for the CVX/model
expression branch of QETLAB `kpNorm`.

`k` is clipped to the number of entries for vectors or to the available
singular-value count for matrices, matching [`top_k_p_norm`](@ref). The atom
does not create a solver or mutate a model. Call
[`add_top_k_p_norm_epigraph!`](@ref) after loading JuMP and a suitable conic
solver.
"""
function top_k_p_norm_epigraph(
    matrix::ComplexAffineMatrix,
    k::Integer,
    p::Real;
    name::Symbol=Symbol(matrix.name, :_top_k_p_norm),
    limits::OptimizationLimits=OptimizationLimits(),
)
    k isa Bool && throw(ArgumentError("k must be a positive integer"))
    k > 0 || throw(ArgumentError("k must be a positive integer; got $k"))
    order = _top_k_p_epigraph_order(p)
    vector_input = matrix.row_dimension == 1 || matrix.column_dimension == 1
    spectral_matrix = vector_input ? _top_k_p_vector_diagonal(matrix) : matrix
    singular_value_count = min(
        spectral_matrix.row_dimension, spectral_matrix.column_dimension
    )
    checked_k = min(Int(k), singular_value_count)
    dilation = _top_k_p_hermitian_dilation(spectral_matrix, Symbol(name, :_dilation))
    auxiliary_variables, psd_blocks = _top_k_p_epigraph_preflight(
        matrix, dilation, singular_value_count, limits
    )
    return TopKPNormEpigraph(
        matrix,
        dilation,
        checked_k,
        order,
        singular_value_count,
        auxiliary_variables,
        psd_blocks,
        limits,
    )
end

"""
    add_top_k_p_norm_epigraph!(model, atom, coordinates; allow_densify=false)

Materialize a [`TopKPNormEpigraph`](@ref) in an optional modeling backend and
return a backend handle whose `epigraph` field is the norm upper-bound
expression. JuMP supplies the concrete method through the package extension.
"""
function add_top_k_p_norm_epigraph!(
    model, atom::TopKPNormEpigraph, coordinates::AbstractVector; allow_densify::Bool=false
)
    return throw(
        ArgumentError(
            "top-k p-norm epigraph materialization requires the optional JuMP extension"
        ),
    )
end
