using LinearAlgebra
using QuantumEntanglementTools
using SparseArrays
using Test

const KpEpigraphQET = QuantumEntanglementTools

if !isdefined(KpEpigraphQET, :TopKPNormEpigraph)
    Base.include(
        KpEpigraphQET,
        joinpath(@__DIR__, "..", "src", "optimization", "top_k_p_norm_epigraph.jl"),
    )
end

@testset "WP3 top-k p-norm solver-neutral epigraph" begin
    variable = KpEpigraphQET.complex_affine_variable(:X, 2, 3; coefficient_type=Float32)
    @test variable isa KpEpigraphQET.ComplexAffineMatrix{Float32}
    @test variable.row_dimension == 2
    @test variable.column_dimension == 3
    @test variable.variable_count == 12
    @test length(variable.terms) == 12
    @test eltype(variable.constant) === ComplexF32
    @test all(issparse(term.coefficient) for term in variable.terms)

    coordinates = Float32[1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6]
    value = Matrix(KpEpigraphQET.evaluate_affine(variable, coordinates))
    @test value == ComplexF32[
        1-im 3-3im 5-5im
        2-2im 4-4im 6-6im
    ]

    atom = KpEpigraphQET.top_k_p_norm_epigraph(variable, 20, 3)
    @test atom isa KpEpigraphQET.TopKPNormEpigraph{Float32,Int}
    @test atom.matrix === variable
    @test atom.k == 2
    @test atom.p == 3
    @test atom.singular_value_count == 2
    @test atom.dilation.dimension == 5
    @test atom.dilation.variable_count == variable.variable_count
    @test atom.psd_block_count == 4
    @test atom.auxiliary_variable_count == 55

    dilation_value = Matrix(KpEpigraphQET.evaluate_affine(atom.dilation, coordinates))
    @test dilation_value == [
        zeros(ComplexF32, 2, 2) value
        adjoint(value) zeros(ComplexF32, 3, 3)
    ]
    @test ishermitian(dilation_value)
    @test sort(real.(eigvals(Hermitian(dilation_value)))) ≈
        sort(vcat(-svdvals(value), zeros(Float32, 1), svdvals(value)))

    vector_variable = KpEpigraphQET.complex_affine_variable(
        :v, 1, 3; coefficient_type=Float64
    )
    vector_atom = KpEpigraphQET.top_k_p_norm_epigraph(vector_variable, 2, Inf)
    @test vector_atom.k == 2
    @test vector_atom.singular_value_count == 3
    @test vector_atom.dilation.dimension == 6
    vector_coordinates = [1.0, -4.0, 2.0, 0.5, 0.0, -1.0]
    vector_value = vec(
        Matrix(KpEpigraphQET.evaluate_affine(vector_variable, vector_coordinates))
    )
    vector_dilation = Matrix(
        KpEpigraphQET.evaluate_affine(vector_atom.dilation, vector_coordinates)
    )
    @test sort(abs.(eigvals(Hermitian(vector_dilation)))) ≈
        sort(repeat(abs.(vector_value), inner=2))
    @test KpEpigraphQET.top_k_p_norm(vector_value, 2, Inf) == maximum(abs, vector_value)

    constant = sparse(ComplexF64[1 0; 0 2])
    coefficient = sparse(ComplexF64[0 1im; 0 0])
    affine = KpEpigraphQET.ComplexAffineMatrix(:affine, constant, [1], [coefficient], 1)
    @test affine.constant !== constant
    @test affine.terms[1].coefficient !== coefficient
    constant[1, 1] = 9
    coefficient[1, 2] = 4
    @test affine.constant[1, 1] == 1
    @test affine.terms[1].coefficient[1, 2] == im

    @test_throws ArgumentError KpEpigraphQET.complex_affine_variable(:bad, 0, 2)
    @test_throws DimensionMismatch KpEpigraphQET.complex_affine_variable(
        :bad, 2, 2; variable_count=7
    )
    @test_throws ArgumentError KpEpigraphQET.ComplexAffineMatrix(
        :duplicate, zeros(2, 2), [1, 1], [ones(2, 2), ones(2, 2)], 2
    )
    @test_throws ArgumentError KpEpigraphQET.ComplexAffineMatrix(
        :nonfinite, fill(Inf, 2, 2), Int[], Matrix{Float64}[], 1
    )
    @test_throws DimensionMismatch KpEpigraphQET.evaluate_affine(affine, [1.0, 2.0])
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(affine, 0, 2)
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(affine, 1, 0.5)
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(affine, 1, NaN)
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(
        variable, 1, 2; limits=KpEpigraphQET.OptimizationLimits(max_variables=60)
    )
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(
        variable, 1, 2; limits=KpEpigraphQET.OptimizationLimits(max_psd_blocks=3)
    )
    @test_throws ArgumentError KpEpigraphQET.top_k_p_norm_epigraph(
        variable, 1, 2; limits=KpEpigraphQET.OptimizationLimits(max_psd_dimension=4)
    )
    @test_throws ArgumentError KpEpigraphQET.add_top_k_p_norm_epigraph!(
        nothing, atom, coordinates; allow_densify=true
    )
end
