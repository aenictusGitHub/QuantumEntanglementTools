# Source-informed independent Julia materialization based on the CVX-expression
# contract of QETLAB IsPSD.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
# QETLAB: Copyright 2014 Nathaniel Johnston, BSD-2-Clause.
# Full upstream terms: licenses/QETLAB-LICENSE.txt.
#
# This is the explicit JuMP materialization for the package-owned affine
# replacement of the pinned CVX-expression branch.

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
