# Source-informed independent Julia implementation based on the CVX-expression
# contract of QETLAB kpNormDual.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# This independently specified solver-neutral replacement leaves the numeric
# array method on `top_k_p_norm_dual`. This file supplies an exact composable
# epigraph for affine vector and matrix expressions. Its JuMP/MOI realization
# is isolated in the optional package extension.

@doc raw"""
    TopKPNormDualEpigraph

Solver-neutral epigraph atom for the dual of the top-`k`,
Schatten/entrywise-`p` norm of a rectangular [`ComplexAffineMatrix`](@ref).

The atom uses singular-value majorization followed by the exact perspective
form of the vector `k`-support norm. For `1 < p < Inf`, with conjugate order
`q`, the lifted variables satisfy

```math
0 \le z_i \le t,\qquad \sum_i z_i \le kt,\qquad
u_i^{1/q}z_i^{1-1/q} \ge s_i,\qquad \sum_i u_i \le t.
```

The endpoint `p == 1` is
`max(maximum(s), sum(s)/k)`, and `p == Inf` is `sum(s)`. A one-row or
one-column input acts on entry magnitudes, matching the numeric vector API.
"""
struct TopKPNormDualEpigraph{T<:Real,P<:Real,Q<:Real}
    matrix::ComplexAffineMatrix{T}
    dilation::HermitianAffineMatrix{T}
    k::Int
    p::P
    conjugate_order::Q
    singular_value_count::Int
    auxiliary_variable_count::Int
    psd_block_count::Int
    scalar_constraint_count::Int
    limits::OptimizationLimits
end

function _top_k_p_dual_conjugate_order(p::Real)
    p == one(p) && return oftype(float(p), Inf)
    isinf(p) && return one(float(p))
    return p / (p - one(p))
end

function _top_k_p_dual_checked_k(k::Integer, count::Int)
    k isa Bool && throw(ArgumentError("k must be a positive integer"))
    k > 0 || throw(ArgumentError("k must be a positive integer; got $k"))
    return k >= count ? count : Int(k)
end

function _top_k_p_dual_epigraph_preflight(
    matrix::ComplexAffineMatrix,
    dilation::HermitianAffineMatrix,
    singular_value_count::Int,
    p::Real,
    limits::OptimizationLimits,
)
    n = BigInt(singular_value_count)
    d = BigInt(dilation.dimension)
    perspective_variables = one(p) < p < Inf ? 2n : BigInt(0)
    auxiliary_variables = n * (d^2 + 2) + 1 + perspective_variables
    total_variables = BigInt(matrix.variable_count) + auxiliary_variables
    total_variables <= limits.max_variables || throw(
        ArgumentError(
            "top-k p-norm dual epigraph needs $total_variables total variables, " *
            "exceeding max_variables=$(limits.max_variables)",
        ),
    )

    psd_blocks = 2n
    psd_blocks <= limits.max_psd_blocks || throw(
        ArgumentError(
            "top-k p-norm dual epigraph needs $psd_blocks PSD blocks, exceeding " *
            "max_psd_blocks=$(limits.max_psd_blocks)",
        ),
    )
    d <= limits.max_psd_dimension || throw(
        ArgumentError(
            "top-k p-norm dual auxiliary blocks have dimension $d, exceeding " *
            "max_psd_dimension=$(limits.max_psd_dimension)",
        ),
    )

    endpoint_constraints = if p == one(p)
        n + 1
    elseif isinf(p)
        one(n)
    else
        2n + 2
    end
    scalar_constraints = max(n - 1, 0) + n + endpoint_constraints
    scalar_constraints <= limits.max_intervals || throw(
        ArgumentError(
            "top-k p-norm dual epigraph needs $scalar_constraints scalar/conic " *
            "constraints, exceeding max_intervals=$(limits.max_intervals)",
        ),
    )

    real_dimension = 2d
    real_triangle_entries = psd_blocks * real_dimension * (real_dimension + 1) ÷ 2
    stored_input_entries = BigInt(nnz(dilation.constant))
    for term in dilation.terms
        stored_input_entries += nnz(term.coefficient)
    end
    total_entries = real_triangle_entries + stored_input_entries
    total_entries <= limits.max_model_entries || throw(
        ArgumentError(
            "top-k p-norm dual real-block formulation needs $total_entries stored " *
            "or scalarized entries, exceeding max_model_entries=" *
            "$(limits.max_model_entries)",
        ),
    )
    return (
        auxiliary_variables=Int(auxiliary_variables),
        psd_blocks=Int(psd_blocks),
        scalar_constraints=Int(scalar_constraints),
    )
end

"""
    top_k_p_norm_dual_epigraph(matrix::ComplexAffineMatrix, k, p; kwargs...)

Construct an exact solver-neutral epigraph atom for the model-expression
branch of QETLAB `kpNormDual`.

`k` is clipped to the number of vector entries or available matrix singular
values. The constructor allocates no solver object and checks all model-size
budgets before the optional extension creates dense real-block expressions.
Use [`add_top_k_p_norm_dual_epigraph!`](@ref) to materialize the atom.
"""
function top_k_p_norm_dual_epigraph(
    matrix::ComplexAffineMatrix,
    k::Integer,
    p::Real;
    name::Symbol=Symbol(matrix.name, :_top_k_p_norm_dual),
    limits::OptimizationLimits=OptimizationLimits(),
)
    order = _top_k_p_epigraph_order(p)
    vector_input = matrix.row_dimension == 1 || matrix.column_dimension == 1
    spectral_matrix = vector_input ? _top_k_p_vector_diagonal(matrix) : matrix
    singular_value_count = min(
        spectral_matrix.row_dimension, spectral_matrix.column_dimension
    )
    checked_k = _top_k_p_dual_checked_k(k, singular_value_count)
    dilation = _top_k_p_hermitian_dilation(spectral_matrix, Symbol(name, :_dilation))
    counts = _top_k_p_dual_epigraph_preflight(
        matrix, dilation, singular_value_count, order, limits
    )
    return TopKPNormDualEpigraph(
        matrix,
        dilation,
        checked_k,
        order,
        _top_k_p_dual_conjugate_order(order),
        singular_value_count,
        counts.auxiliary_variables,
        counts.psd_blocks,
        counts.scalar_constraints,
        limits,
    )
end

"""
    add_top_k_p_norm_dual_epigraph!(
        model, atom, coordinates; allow_densify=false
    )

Materialize a [`TopKPNormDualEpigraph`](@ref) in an optional modeling backend.
The returned backend handle has an `epigraph` field. JuMP supplies the concrete
method through the package extension; no optimizer is selected implicitly.
"""
function add_top_k_p_norm_dual_epigraph!(
    model,
    atom::TopKPNormDualEpigraph,
    coordinates::AbstractVector;
    allow_densify::Bool=false,
)
    return throw(
        ArgumentError(
            "top-k p-norm dual epigraph materialization requires the optional " *
            "JuMP extension",
        ),
    )
end
