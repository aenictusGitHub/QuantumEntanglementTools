# Included by QuantumEntanglementToolsJuMPExt after top_k_p_norm_epigraph.jl.
# Only public JuMP and MathOptInterface APIs are used.

function _top_k_p_dual_singular_majorant!(
    model::JuMP.GenericModel{T},
    atom::QET.TopKPNormDualEpigraph{T},
    coordinates::AbstractVector{<:JuMP.AbstractJuMPScalar},
) where {T<:Real}
    n = atom.singular_value_count
    dimension = atom.dilation.dimension
    singular_values = JuMP.@variable(model, [1:n], lower_bound = 0)
    thresholds = JuMP.@variable(model, [1:n])
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
    return (
        singular_values,
        thresholds,
        hermitian_auxiliaries=auxiliaries,
        constraints=(
            ordering=ordering_constraints,
            psd=psd_constraints,
            majorization=majorization_constraints,
        ),
    )
end

function QET.add_top_k_p_norm_dual_epigraph!(
    model::JuMP.GenericModel{T},
    atom::QET.TopKPNormDualEpigraph{T},
    coordinates::AbstractVector{<:JuMP.AbstractJuMPScalar};
    allow_densify::Bool=false,
) where {T<:Real}
    allow_densify || throw(
        ArgumentError(
            "top-k p-norm dual materialization allocates dense real-block " *
            "expression matrices; pass allow_densify=true",
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

    majorant = _top_k_p_dual_singular_majorant!(model, atom, coordinates)
    singular_values = majorant.singular_values
    n = atom.singular_value_count
    epigraph = JuMP.@variable(model, lower_bound = 0)
    perspective_variables = nothing
    gauge_constraints = JuMP.ConstraintRef[]

    formulation = if atom.p == one(atom.p)
        for index in 1:n
            push!(
                gauge_constraints,
                JuMP.@constraint(model, singular_values[index] <= epigraph),
            )
        end
        push!(
            gauge_constraints,
            JuMP.@constraint(model, sum(singular_values) <= atom.k * epigraph),
        )
        :linf_and_scaled_l1
    elseif isinf(atom.p)
        push!(gauge_constraints, JuMP.@constraint(model, sum(singular_values) <= epigraph))
        :l1
    else
        allocation = JuMP.@variable(model, [1:n], lower_bound = 0)
        perspective = JuMP.@variable(model, [1:n], lower_bound = 0)
        power = Float64(inv(atom.conjugate_order))
        zero(power) < power < one(power) || throw(
            ArgumentError(
                "the conjugate-order power-cone exponent is not representable " *
                "strictly inside (0, 1) as Float64",
            ),
        )
        for index in 1:n
            push!(gauge_constraints, JuMP.@constraint(model, allocation[index] <= epigraph))
            push!(
                gauge_constraints,
                JuMP.@constraint(
                    model,
                    [perspective[index], allocation[index], singular_values[index]] in
                        MOI.PowerCone(power),
                ),
            )
        end
        push!(
            gauge_constraints,
            JuMP.@constraint(model, sum(allocation) <= atom.k * epigraph),
        )
        push!(gauge_constraints, JuMP.@constraint(model, sum(perspective) <= epigraph))
        perspective_variables = (; allocation, perspective)
        :k_support_perspective_power_cones
    end

    return (
        epigraph,
        singular_value_majorant=singular_values,
        thresholds=majorant.thresholds,
        hermitian_auxiliaries=majorant.hermitian_auxiliaries,
        perspective_variables,
        constraints=merge(majorant.constraints, (gauge=gauge_constraints,)),
        formulation,
        certified_epigraph=true,
    )
end
