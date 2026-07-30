# Independently designed solver-neutral replacement for the CVX-expression
# branch of QETLAB IsPSD.m at
# d8589610f00cff106537268dee2e2a1153f3a601.
#
# Numeric matrices continue to use is_positive_semidefinite. Affine model data
# uses this explicit constraint builder so a predicate is never confused with a
# mutation of an ambient optimization model.

"""
    positive_semidefinite_constraint(matrix::HermitianAffineMatrix)

Return an owned solver-neutral positive-semidefinite constraint block.

The result is a [`HermitianAffineMatrix`](@ref) suitable for the
`psd_constraints` field of [`SemidefiniteProgram`](@ref). It contains no JuMP,
MathOptInterface, or optimizer object. The input constant and coefficient
matrices are copied, so later mutation of the caller's sparse storage cannot
change the reviewed model.

This function is the Julia-native no-loss replacement for calling QETLAB
`IsPSD` on a CVX expression. Numeric arrays instead use
[`is_positive_semidefinite`](@ref), which returns a tri-state
[`MatrixPredicateResult`](@ref).
"""
function positive_semidefinite_constraint(matrix::HermitianAffineMatrix{T}) where {T<:Real}
    return _optimization_matrix_convert(matrix, T)
end
