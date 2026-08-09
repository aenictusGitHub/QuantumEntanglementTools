# Source-informed independent Julia implementation based on the executable
# contract of QETLAB MatsumotoFidelity.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# This independently specified solver-neutral replacement covers the pinned
# CVX-expression branch. The numeric branch remains in
# measures/scalar_measures.jl. This file exposes the exact semidefinite
# hypograph lift instead of mutating an ambient CVX model or returning a
# backend-owned expression.

"""
    MatsumotoFidelityModel

Solver-neutral semidefinite lift of the Matsumoto fidelity.

For Hermitian affine matrices `rho` and `sigma`, the model introduces a
Hermitian coupling `X` and represents

```math
F_{\\mathrm{M}}(\\rho,\\sigma)
=
\\max_X\\left\\{
\\mathop{\\mathrm{tr}}X:
\\begin{bmatrix}
\\rho & X\\\\
X & \\sigma
\\end{bmatrix}\\succeq 0
\\right\\}.
```

`objective` is `tr(X)`, `psd_constraint` is the displayed block, and
`coupling` is a named primal view. All fields contain package-owned numeric
affine data; no JuMP, MathOptInterface, or optimizer object is retained.
"""
struct MatsumotoFidelityModel{T<:Real}
    dimension::Int
    input_variable_count::Int
    variable_count::Int
    coupling_variables::UnitRange{Int}
    objective::AffineScalar{T}
    coupling::HermitianAffineMatrix{T}
    psd_constraint::HermitianAffineMatrix{T}
end

function _matsumoto_model_block(
    top_left::SparseMatrixCSC{Complex{T},Int},
    top_right::SparseMatrixCSC{Complex{T},Int},
    bottom_left::SparseMatrixCSC{Complex{T},Int},
    bottom_right::SparseMatrixCSC{Complex{T},Int},
) where {T<:Real}
    return sparse([top_left top_right; bottom_left bottom_right])
end

function _matsumoto_model_add_term!(
    terms::Dict{Int,SparseMatrixCSC{Complex{T},Int}},
    variable::Int,
    coefficient::SparseMatrixCSC{Complex{T},Int},
) where {T<:Real}
    if haskey(terms, variable)
        terms[variable] = terms[variable] + coefficient
    else
        terms[variable] = coefficient
    end
    return nothing
end

function _matsumoto_model_preflight(
    dimension::Int,
    input_variable_count::Int,
    rho_constant,
    rho_terms,
    sigma_constant,
    sigma_terms,
    limits::OptimizationLimits,
)
    total_variables = BigInt(input_variable_count) + BigInt(dimension)^2
    total_variables <= limits.max_variables || throw(
        ArgumentError(
            "Matsumoto-fidelity lift needs $total_variables variables, exceeding " *
            "max_variables=$(limits.max_variables)",
        ),
    )
    block_dimension = BigInt(2) * dimension
    block_dimension <= limits.max_psd_dimension || throw(
        ArgumentError(
            "Matsumoto-fidelity PSD block has dimension $block_dimension, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )
    limits.max_psd_blocks >= 1 ||
        throw(ArgumentError("Matsumoto-fidelity lift needs one PSD block"))

    stored_entries = BigInt(nnz(rho_constant)) + BigInt(nnz(sigma_constant))
    for term in rho_terms
        stored_entries += nnz(term.coefficient)
    end
    for term in sigma_terms
        stored_entries += nnz(term.coefficient)
    end
    # Each independent Hermitian coupling coefficient appears in both
    # off-diagonal blocks.
    stored_entries += BigInt(2) * dimension^2
    stored_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "Matsumoto-fidelity lift stores at least $stored_entries entries, " *
            "exceeding max_model_entries=$(limits.max_model_entries)",
        ),
    )

    # The optional extension realifies a 2d-by-2d complex block into a
    # 4d-by-4d real symmetric block.
    real_dimension = BigInt(4) * dimension
    real_triangle_entries = real_dimension * (real_dimension + 1) ÷ 2
    real_triangle_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "Matsumoto-fidelity real block needs $real_triangle_entries triangle " *
            "entries, exceeding max_model_entries=$(limits.max_model_entries)",
        ),
    )
    return Int(total_variables)
end

function _matsumoto_fidelity_model(
    rho_constant::AbstractMatrix{<:Number},
    rho_terms,
    sigma_constant::AbstractMatrix{<:Number},
    sigma_terms,
    input_variable_count::Int;
    name::Symbol,
    limits::OptimizationLimits,
)
    dimension = size(rho_constant, 1)
    size(rho_constant) == (dimension, dimension) ||
        throw(DimensionMismatch("rho must be square"))
    size(sigma_constant) == (dimension, dimension) || throw(
        DimensionMismatch(
            "rho and sigma must have the same square size; got " *
            "$(size(rho_constant)) and $(size(sigma_constant))",
        ),
    )
    T = promote_type(
        _optimization_matrix_real_type(rho_constant),
        _optimization_matrix_real_type(sigma_constant),
        ((
            _optimization_matrix_real_type(term.coefficient) for
            term in Iterators.flatten((rho_terms, sigma_terms))
        ))...,
    )
    isconcretetype(T) && T <: Real ||
        throw(ArgumentError("model coefficient type must be a concrete real type"))
    rho0 = sparse(Complex{T}.(rho_constant))
    sigma0 = sparse(Complex{T}.(sigma_constant))
    total_variables = _matsumoto_model_preflight(
        dimension, input_variable_count, rho0, rho_terms, sigma0, sigma_terms, limits
    )
    first_coupling_variable = input_variable_count + 1
    coupling = hermitian_variable(
        Symbol(name, :_coupling),
        dimension;
        first_variable=first_coupling_variable,
        variable_count=total_variables,
        coefficient_type=T,
    )

    zero_block = spzeros(Complex{T}, dimension, dimension)
    constant = _matsumoto_model_block(rho0, zero_block, zero_block, sigma0)
    block_terms = Dict{Int,SparseMatrixCSC{Complex{T},Int}}()
    for term in rho_terms
        coefficient = sparse(Complex{T}.(term.coefficient))
        _matsumoto_model_add_term!(
            block_terms,
            term.variable,
            _matsumoto_model_block(coefficient, zero_block, zero_block, zero_block),
        )
    end
    for term in sigma_terms
        coefficient = sparse(Complex{T}.(term.coefficient))
        _matsumoto_model_add_term!(
            block_terms,
            term.variable,
            _matsumoto_model_block(zero_block, zero_block, zero_block, coefficient),
        )
    end
    for term in coupling.terms
        coefficient = sparse(Complex{T}.(term.coefficient))
        _matsumoto_model_add_term!(
            block_terms,
            term.variable,
            _matsumoto_model_block(zero_block, coefficient, coefficient, zero_block),
        )
    end
    variables = sort!(collect(keys(block_terms)))
    coefficients = [block_terms[variable] for variable in variables]
    psd_constraint = HermitianAffineMatrix(
        Symbol(name, :_block), constant, variables, coefficients, total_variables
    )
    objective = trace_affine(coupling)
    return MatsumotoFidelityModel(
        dimension,
        input_variable_count,
        total_variables,
        first_coupling_variable:total_variables,
        objective,
        coupling,
        psd_constraint,
    )
end

"""
    matsumoto_fidelity_model(rho, sigma; kwargs...)

Construct the package-owned semidefinite lift for two numeric density
matrices. Inputs receive the same strict validation as
[`matsumoto_fidelity`](@ref): they are never normalized, symmetrized, clipped,
or repaired. Sparse spectral validation requires `allow_densify=true`.

This function only constructs model data. Use
[`matsumoto_fidelity_problem`](@ref) for a directly solvable fixed-input SDP.
"""
function matsumoto_fidelity_model(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    name::Symbol=:matsumoto_fidelity,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    support_boundary_policy::Symbol=:reject,
    limits::OptimizationLimits=OptimizationLimits(),
)
    size(rho) == size(sigma) || throw(
        DimensionMismatch(
            "rho and sigma must have the same size; got $(size(rho)) and $(size(sigma))"
        ),
    )
    size(rho, 1) == size(rho, 2) ||
        throw(DimensionMismatch("rho and sigma must be square; got size $(size(rho))"))
    policy = _tierd_matsumoto_support_policy(support_boundary_policy)
    ishermitian(rho) ||
        throw(ArgumentError("rho must be exactly Hermitian; no repair is applied"))
    ishermitian(sigma) ||
        throw(ArgumentError("sigma must be exactly Hermitian; no repair is applied"))
    rho_analysis = _tierd_density_analysis(
        rho;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="matsumoto_fidelity_model",
    )
    sigma_analysis = _tierd_density_analysis(
        sigma;
        atol=atol,
        rtol=rtol,
        allow_densify=allow_densify,
        boundary_policy=:record,
        operation="matsumoto_fidelity_model",
    )
    _tierd_matsumoto_support(rho_analysis, policy, "rho")
    _tierd_matsumoto_support(sigma_analysis, policy, "sigma")
    return _matsumoto_fidelity_model(
        rho_analysis.matrix, (), sigma_analysis.matrix, (), 0; name, limits
    )
end

"""
    matsumoto_fidelity_model(rho::HermitianAffineMatrix,
                              sigma::HermitianAffineMatrix; kwargs...)

Construct a composable semidefinite lift for affine Hermitian model data.
`rho` and `sigma` must have the same dimension and coordinate count. The
builder appends `dimension^2` real coordinates for the Hermitian coupling.

No density-matrix claim can be checked for symbolic affine inputs. The caller
must add any required positivity, trace, and outer-model constraints
explicitly.
"""
function matsumoto_fidelity_model(
    rho::HermitianAffineMatrix,
    sigma::HermitianAffineMatrix;
    name::Symbol=:matsumoto_fidelity,
    limits::OptimizationLimits=OptimizationLimits(),
)
    rho.dimension == sigma.dimension || throw(
        DimensionMismatch(
            "rho and sigma dimensions differ: $(rho.dimension) and $(sigma.dimension)"
        ),
    )
    rho.variable_count == sigma.variable_count || throw(
        DimensionMismatch(
            "rho and sigma must use the same affine coordinate count; got " *
            "$(rho.variable_count) and $(sigma.variable_count)",
        ),
    )
    return _matsumoto_fidelity_model(
        rho.constant,
        rho.terms,
        sigma.constant,
        sigma.terms,
        rho.variable_count;
        name,
        limits,
    )
end

function _matsumoto_lift_affine(function_data::AffineScalar, total_variables::Int)
    indices, values = findnz(function_data.coefficients)
    return AffineScalar(function_data.constant, indices, values, total_variables)
end

function _matsumoto_lift_constraint(constraint::AffineEquality, total_variables::Int)
    return AffineEquality(
        _matsumoto_lift_affine(constraint.function_data, total_variables), constraint.name
    )
end

function _matsumoto_lift_constraint(constraint::AffineInterval, total_variables::Int)
    return AffineInterval(
        _matsumoto_lift_affine(constraint.function_data, total_variables),
        constraint.lower,
        constraint.upper,
        constraint.name,
    )
end

function _matsumoto_lift_matrix(matrix::HermitianAffineMatrix, total_variables::Int)
    return HermitianAffineMatrix(
        matrix.name,
        matrix.constant,
        [term.variable for term in matrix.terms],
        [term.coefficient for term in matrix.terms],
        total_variables,
    )
end

function _matsumoto_lift_point(
    point, input_variables::Int, total_variables::Int, name::AbstractString
)
    point === nothing && return nothing
    point isa AbstractVector ||
        throw(ArgumentError("$name must be nothing or an AbstractVector"))
    firstindex(point) == 1 || throw(ArgumentError("$name must use one-based indexing"))
    if length(point) == total_variables
        return collect(point)
    end
    length(point) == input_variables || throw(
        DimensionMismatch(
            "$name must contain either $input_variables input coordinates or all " *
            "$total_variables lifted coordinates",
        ),
    )
    return vcat(collect(point), zeros(eltype(point), total_variables - input_variables))
end

"""
    matsumoto_fidelity_problem(rho, sigma; kwargs...)

Return a directly solvable [`SemidefiniteProgram`](@ref) for fixed numeric
density matrices. The zero coupling is retained as a known feasible point.
With no explicit backend, [`solve_optimization`](@ref) returns
`OptimizationBackendUnavailable`.
"""
function matsumoto_fidelity_problem(
    rho::AbstractMatrix{<:Number},
    sigma::AbstractMatrix{<:Number};
    name::Symbol=:matsumoto_fidelity,
    atol=nothing,
    rtol=nothing,
    allow_densify::Bool=false,
    support_boundary_policy::Symbol=:reject,
    limits::OptimizationLimits=OptimizationLimits(),
)
    model = matsumoto_fidelity_model(
        rho, sigma; name, atol, rtol, allow_densify, support_boundary_policy, limits
    )
    T = typeof(model.objective.constant)
    zero_point = zeros(T, model.variable_count)
    return SemidefiniteProgram(
        Symbol(name, :_problem),
        :maximize,
        model.variable_count,
        model.objective;
        psd_constraints=[model.psd_constraint],
        primal_views=[model.coupling],
        initial_point=zero_point,
        known_feasible_point=zero_point,
        limits,
        metadata=(
            upstream_function="MatsumotoFidelity",
            formulation=:matrix_geometric_mean_hypograph,
            value_semantics=:solver_optimum_not_theorem_certificate,
            input_dimension=model.dimension,
            input_variable_count=0,
        ),
    )
end

"""
    matsumoto_fidelity_problem(rho::HermitianAffineMatrix,
                                sigma::HermitianAffineMatrix; kwargs...)

Build an SDP that maximizes Matsumoto fidelity jointly with caller-supplied
constraints on affine `rho` and `sigma`. Base constraints and views using the
inputs' coordinate count are lifted automatically when the coupling
coordinates are appended.
"""
function matsumoto_fidelity_problem(
    rho::HermitianAffineMatrix,
    sigma::HermitianAffineMatrix;
    name::Symbol=:matsumoto_fidelity,
    equalities=AffineEquality[],
    intervals=AffineInterval[],
    psd_constraints=HermitianAffineMatrix[],
    primal_views=HermitianAffineMatrix[],
    initial_point=nothing,
    known_feasible_point=nothing,
    limits::OptimizationLimits=OptimizationLimits(),
)
    model = matsumoto_fidelity_model(rho, sigma; name, limits)
    base_count = model.input_variable_count
    for constraint in Iterators.flatten((equalities, intervals))
        length(constraint.function_data.coefficients) == base_count || throw(
            DimensionMismatch(
                "base constraint $(constraint.name) has the wrong coordinate count"
            ),
        )
    end
    for matrix in Iterators.flatten((psd_constraints, primal_views))
        matrix.variable_count == base_count || throw(
            DimensionMismatch("base matrix $(matrix.name) has the wrong coordinate count"),
        )
    end
    lifted_equalities = [
        _matsumoto_lift_constraint(constraint, model.variable_count) for
        constraint in equalities
    ]
    lifted_intervals = [
        _matsumoto_lift_constraint(constraint, model.variable_count) for
        constraint in intervals
    ]
    lifted_psd = [
        _matsumoto_lift_matrix(matrix, model.variable_count) for matrix in psd_constraints
    ]
    lifted_views = [
        _matsumoto_lift_matrix(matrix, model.variable_count) for matrix in primal_views
    ]
    push!(lifted_psd, model.psd_constraint)
    push!(lifted_views, model.coupling)
    lifted_initial = _matsumoto_lift_point(
        initial_point, base_count, model.variable_count, "initial_point"
    )
    lifted_feasible = _matsumoto_lift_point(
        known_feasible_point, base_count, model.variable_count, "known_feasible_point"
    )
    return SemidefiniteProgram(
        Symbol(name, :_problem),
        :maximize,
        model.variable_count,
        model.objective;
        equalities=lifted_equalities,
        intervals=lifted_intervals,
        psd_constraints=lifted_psd,
        primal_views=lifted_views,
        initial_point=lifted_initial,
        known_feasible_point=lifted_feasible,
        limits,
        metadata=(
            upstream_function="MatsumotoFidelity",
            formulation=:matrix_geometric_mean_hypograph,
            value_semantics=:solver_optimum_not_theorem_certificate,
            input_dimension=model.dimension,
            input_variable_count=base_count,
        ),
    )
end
