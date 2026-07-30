# Explicit JuMP materialization for the package-owned affine replacement of
# QETLAB IsPSD's CVX-expression branch.

function QET.positive_semidefinite_constraint(
    model::JuMP.AbstractModel,
    matrix::QET.HermitianAffineMatrix,
    variables::AbstractVector{<:JuMP.AbstractJuMPScalar};
    allow_densify::Bool=false,
)
    allow_densify || throw(
        ArgumentError(
            "JuMP materialization of a complex PSD block requires allow_densify=true"
        ),
    )
    length(variables) == matrix.variable_count || throw(
        DimensionMismatch(
            "the affine PSD block expects $(matrix.variable_count) coordinates; " *
            "got $(length(variables))",
        ),
    )
    block = _jump_real_block(matrix, variables)
    return JuMP.@constraint(model, block in JuMP.PSDCone())
end
