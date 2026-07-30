# Included by QuantumEntanglementToolsJuMPExt. This file uses only public JuMP
# and MathOptInterface APIs to materialize the package-owned top-k p-norm
# epigraph atom.

function _top_k_p_zero_affine(::Type{T}, ::Type{V}) where {T<:Real,V}
    return JuMP.GenericAffExpr(zero(T), Pair{V,T}[])
end

function _top_k_p_hermitian_real_block(
    model::JuMP.GenericModel{T}, dimension::Int
) where {T<:Real}
    diagonal = JuMP.@variable(model, [1:dimension])
    triangle_count = dimension * (dimension - 1) ÷ 2
    real_upper = JuMP.@variable(model, [1:triangle_count])
    imaginary_upper = JuMP.@variable(model, [1:triangle_count])
    V = eltype(diagonal)
    real_part = [_top_k_p_zero_affine(T, V) for _ in 1:dimension, _ in 1:dimension]
    imaginary_part = [_top_k_p_zero_affine(T, V) for _ in 1:dimension, _ in 1:dimension]
    for index in 1:dimension
        JuMP.add_to_expression!(real_part[index, index], one(T), diagonal[index])
    end
    position = 1
    for column in 2:dimension, row in 1:(column - 1)
        JuMP.add_to_expression!(real_part[row, column], one(T), real_upper[position])
        JuMP.add_to_expression!(real_part[column, row], one(T), real_upper[position])
        JuMP.add_to_expression!(
            imaginary_part[row, column], one(T), imaginary_upper[position]
        )
        JuMP.add_to_expression!(
            imaginary_part[column, row], -one(T), imaginary_upper[position]
        )
        position += 1
    end

    real_dimension = 2dimension
    block = Matrix{eltype(real_part)}(undef, real_dimension, real_dimension)
    for column in 1:dimension, row in 1:dimension
        block[row, column] = real_part[row, column]
        block[row, column + dimension] = -imaginary_part[row, column]
        block[row + dimension, column] = imaginary_part[row, column]
        block[row + dimension, column + dimension] = real_part[row, column]
    end
    return (
        block=LinearAlgebra.Symmetric(block),
        trace=sum(diagonal),
        coordinates=(; diagonal, real_upper, imaginary_upper),
    )
end

function _top_k_p_slack_block(
    auxiliary_block, dilation_block, threshold, ::Type{T}
) where {T<:Real}
    dimension = size(auxiliary_block, 1)
    size(dilation_block) == (dimension, dimension) ||
        throw(DimensionMismatch("real-block dimensions are inconsistent"))
    block = Matrix{eltype(auxiliary_block)}(undef, dimension, dimension)
    for column in 1:dimension, row in 1:dimension
        entry = auxiliary_block[row, column] - dilation_block[row, column]
        row == column && (entry += threshold)
        block[row, column] = entry
    end
    return LinearAlgebra.Symmetric(block)
end

function QET.add_top_k_p_norm_epigraph!(
    model::JuMP.GenericModel{T},
    atom::QET.TopKPNormEpigraph{T},
    coordinates::AbstractVector{<:JuMP.AbstractJuMPScalar};
    allow_densify::Bool=false,
) where {T<:Real}
    allow_densify || throw(
        ArgumentError(
            "top-k p-norm materialization allocates dense real-block expression " *
            "matrices; pass allow_densify=true",
        ),
    )
    firstindex(coordinates) == 1 ||
        throw(ArgumentError("coordinates must use one-based indexing"))
    length(coordinates) == atom.matrix.variable_count || throw(
        DimensionMismatch(
            "expected $(atom.matrix.variable_count) model coordinates, got " *
            "$(length(coordinates))",
        ),
    )
    n = atom.singular_value_count
    dimension = atom.dilation.dimension
    singular_values = JuMP.@variable(model, [1:n], lower_bound = 0)
    thresholds = JuMP.@variable(model, [1:n])
    epigraph = JuMP.@variable(model, lower_bound = 0)

    ordering_constraints = JuMP.ConstraintRef[]
    for index in 1:(n - 1)
        push!(
            ordering_constraints,
            JuMP.@constraint(model, singular_values[index] >= singular_values[index + 1]),
        )
    end

    dilation_block = _jump_real_block(atom.dilation, coordinates)
    psd_constraints = JuMP.ConstraintRef[]
    majorization_constraints = JuMP.ConstraintRef[]
    auxiliaries = NamedTuple[]
    for index in 1:n
        auxiliary = _top_k_p_hermitian_real_block(model, dimension)
        push!(auxiliaries, auxiliary.coordinates)
        push!(psd_constraints, JuMP.@constraint(model, auxiliary.block in JuMP.PSDCone()))
        slack = _top_k_p_slack_block(auxiliary.block, dilation_block, thresholds[index], T)
        push!(psd_constraints, JuMP.@constraint(model, slack in JuMP.PSDCone()))
        push!(
            majorization_constraints,
            JuMP.@constraint(
                model,
                index * thresholds[index] + auxiliary.trace <=
                    sum(singular_values[1:index]),
            ),
        )
    end

    norm_constraint = if isinf(atom.p)
        JuMP.@constraint(model, singular_values[1] <= epigraph)
    elseif atom.p == one(atom.p)
        JuMP.@constraint(model, sum(singular_values[1:(atom.k)]) <= epigraph)
    else
        JuMP.@constraint(
            model,
            [epigraph; singular_values[1:(atom.k)]] in
                MOI.NormCone(Float64(atom.p), atom.k + 1),
        )
    end
    return (
        epigraph,
        singular_value_majorant=singular_values,
        thresholds,
        hermitian_auxiliaries=auxiliaries,
        constraints=(
            ordering=ordering_constraints,
            psd=psd_constraints,
            majorization=majorization_constraints,
            norm=norm_constraint,
        ),
        formulation=:ky_fan_majorization_and_norm_cone,
        certified_epigraph=true,
    )
end
