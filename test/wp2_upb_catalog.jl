using LinearAlgebra
using Random
using SparseArrays
using Test
using QuantumEntanglementTools: QuantumEntanglementTools

module IsolatedUPBCatalog
using LinearAlgebra
using Random
using SparseArrays

include("../src/dimensions.jl")
include("../src/tensor_products.jl")
include("../src/entanglement/upb_size.jl")
include("../src/entanglement/upb_catalog.jl")
include("../src/entanglement/upb_analysis.jl")
end

const CatalogUPB = IsolatedUPBCatalog
const PublicUPB = QuantumEntanglementTools

struct _UPBCZeroBasedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    storage::V
end
Base.size(vector::_UPBCZeroBasedVector) = size(vector.storage)
Base.axes(vector::_UPBCZeroBasedVector) = (0:(length(vector.storage) - 1),)
Base.IndexStyle(::Type{<:_UPBCZeroBasedVector}) = IndexLinear()
Base.getindex(vector::_UPBCZeroBasedVector, index::Int) = vector.storage[index + 1]

function _wp2_upbc_global_column(factors, column)
    return foldl(kron, (factor[:, column] for factor in factors))
end

function _wp2_upbc_check_structure(construction; atol=2e-11)
    @test construction isa CatalogUPB.UPBConstruction
    @test construction.status === :constructed
    @test construction.cardinality == size(construction.global_vectors, 2)
    @test size(construction.global_vectors, 1) == prod(construction.dimensions)
    @test length(construction.local_factors) == length(construction.dimensions)
    @test all(
        size(construction.local_factors[party]) ==
        (construction.dimensions[party], construction.cardinality) for
        party in eachindex(construction.dimensions)
    )
    @test all(
        isapprox(
            construction.global_vectors[:, column],
            _wp2_upbc_global_column(construction.local_factors, column);
            atol=atol,
            rtol=atol,
        ) for column in 1:construction.cardinality
    )
    gram = adjoint(construction.global_vectors) * construction.global_vectors
    @test gram ≈ I atol=atol rtol=atol
    @test construction.normalized
    @test construction.pairwise_orthogonal
    @test construction.normalization_residual <= atol
    @test construction.orthogonality_residual <= atol
    @test construction.product_residual <= atol
    @test construction.source_revision == "d8589610f00cff106537268dee2e2a1153f3a601"
    @test construction.symbolic_simplification === :not_required_for_executable_public_route
    @test !isempty(construction.reference)
    @test occursin("UPBConstruction", sprint(show, construction))
    return construction
end

function _wp2_upbc_check_unextendible(construction; max_partitions=2_000_000)
    analysis = CatalogUPB.is_upb(
        construction.local_factors; max_partitions=max_partitions, max_work=1_000_000_000
    )
    @test analysis.status === :upb
    @test analysis.reason === :unextendible
    @test analysis.is_upb === true
    return analysis
end

function _wp2_upbc_ray_supports(factor; atol=1e-10)
    supports = Vector{UInt64}()
    used = falses(size(factor, 2))
    for column in axes(factor, 2)
        used[column] && continue
        reference = @view factor[:, column]
        support = UInt64(0)
        for candidate in column:size(factor, 2)
            used[candidate] && continue
            overlap = abs(dot(reference, @view factor[:, candidate]))
            if isapprox(overlap, one(overlap); atol=atol, rtol=atol)
                used[candidate] = true
                support |= UInt64(1) << (candidate - 1)
            end
        end
        push!(supports, support)
    end
    return supports
end

function _wp2_upbc_qubit_cover_bound(construction)
    construction.cardinality <= 64 || error("test bit mask supports at most 64 states")
    reachable = Set([UInt64(0)])
    for factor in construction.local_factors
        next_reachable = Set{UInt64}()
        for prior in reachable, support in _wp2_upbc_ray_supports(factor)
            push!(next_reachable, prior | support)
        end
        reachable = next_reachable
    end
    full = (UInt64(1) << construction.cardinality) - UInt64(1)
    maximum_covered = maximum(count_ones, reachable)
    return full in reachable, maximum_covered
end

@testset "WP2 UPB full construction catalog" begin
    @testset "fixed and parameterized named families" begin
        cases = [
            ("Pyramid", (), (3, 3), 5, :pyramid),
            ("Tiles", (), (3, 3), 5, :tiles),
            ("GenTiles1", (4,), (4, 4), 9, :generalized_tiles_1),
            ("GenTiles2", (3, 4), (3, 4), 7, :generalized_tiles_2),
            ("Min4x4", (), (4, 4), 8, :minimum_4x4),
            ("QuadRes", (3,), (3, 3), 5, :quad_residue),
            ("SixParam", ([0.31, 0.47, -0.2, 0.53, 0.61, 0.4],), (3, 3), 5, :six_parameter),
            ("Shifts", (), (2, 2, 2), 4, :generalized_shifts),
            ("GenShifts", (5,), (2, 2, 2, 2, 2), 6, :generalized_shifts),
            ("Feng2x2x2x2", (), (2, 2, 2, 2), 6, :feng_2x2x2x2),
            ("John2^8", (), ntuple(_ -> 2, 8), 11, :johnston_2_power_8),
            ("John2^4k", (8,), ntuple(_ -> 2, 8), 12, :johnston_2_power_4k),
            ("CJBip46", (), (4, 6), 10, :chen_johnston_4x6),
            ("Feng2x2x3", (), (2, 2, 3), 6, :feng_2x2x3),
            ("Feng2x2x5", (), (2, 2, 5), 8, :feng_2x2x5),
            ("Feng4x4", (), (4, 4), 8, :feng_4x4),
            ("Feng2x2x2x4", (), (2, 2, 2, 4), 8, :feng_2x2x2x4),
            ("Feng2x2x2x2x5", (), (2, 2, 2, 2, 5), 10, :feng_2x2x2x2x5),
        ]

        constructions = Dict{Symbol,Any}()
        for (name, arguments, dimensions, cardinality, family) in cases
            construction = _wp2_upbc_check_structure(CatalogUPB.upb(name, arguments...))
            @test construction.family === family
            @test construction.dimensions == dimensions
            @test construction.cardinality == cardinality
            @test construction.construction_kind in (:closed_form, :complete_product_basis)
            @test construction.deterministic
            @test !construction.rng_used
            @test construction.attempts == 0
            @test construction.random_draws == 0
            constructions[family] = construction
        end

        # Exhaustive independent partition/rank analysis is practical for
        # every named family except the two large qubit graph catalogs.
        for family in (
            :pyramid,
            :tiles,
            :generalized_tiles_1,
            :generalized_tiles_2,
            :minimum_4x4,
            :quad_residue,
            :six_parameter,
            :generalized_shifts,
            :feng_2x2x2x2,
            :chen_johnston_4x6,
            :feng_2x2x3,
            :feng_2x2x5,
            :feng_4x4,
            :feng_2x2x2x4,
            :feng_2x2x2x2x5,
        )
            max_partitions = family === :feng_2x2x2x2x5 ? 500_000 : 2_000_000
            _wp2_upbc_check_unextendible(
                constructions[family]; max_partitions=max_partitions
            )
        end

        # The independent qubit ray-cover dynamic program implements the
        # paper's orthogonality-graph argument without enumerating 30M+
        # labelled partitions.
        for family in (:johnston_2_power_8, :johnston_2_power_4k)
            extendible, maximum_covered = _wp2_upbc_qubit_cover_bound(constructions[family])
            @test !extendible
            @test maximum_covered < constructions[family].cardinality
        end
        @test _wp2_upbc_qubit_cover_bound(constructions[:johnston_2_power_4k])[2] == 11
    end

    @testset "randomized theorem constructions and RNG semantics" begin
        random_cases = [
            ("AlonLovasz", ((3, 4),), (3, 4), 6, :alon_lovasz, 101),
            ("CJBip", ((6, 6),), (6, 6), 12, :chen_johnston_bipartite, 202),
            ("CJ4k1", (5,), (2, 2, 5), 8, :chen_johnston_4k1, 303),
        ]
        for (name, arguments, dimensions, cardinality, family, seed) in random_cases
            first = _wp2_upbc_check_structure(
                CatalogUPB.upb(MersenneTwister(seed), name, arguments...)
            )
            second = CatalogUPB.upb(MersenneTwister(seed), name, arguments...)
            @test first.family === family
            @test first.dimensions == dimensions
            @test first.cardinality == cardinality
            @test first.construction_kind === :randomized_full_spark
            @test !first.deterministic
            @test first.rng_used
            @test first.attempts >= 1
            @test first.random_draws > 0
            @test first.minors_checked > 0
            @test first.global_vectors == second.global_vectors
            @test first.local_factors == second.local_factors
            _wp2_upbc_check_unextendible(first)
        end

        float32 = CatalogUPB.upb(
            MersenneTwister(404), :alon_lovasz, (3, 4); real_type=Float32
        )
        @test eltype(float32.global_vectors) == Float32
        @test float32.normalization_residual isa Float32
        @test float32.orthogonality_residual isa Float32

        @test_throws ArgumentError CatalogUPB.upb(:alon_lovasz, (3, 4))
        @test_throws ArgumentError CatalogUPB.upb(
            MersenneTwister(1), :alon_lovasz, (3, 4); real_type=BigFloat
        )

        # Supplying a private RNG never touches Julia's global random stream.
        Random.seed!(0x51a7)
        expected = rand(UInt64)
        Random.seed!(0x51a7)
        CatalogUPB.upb(MersenneTwister(55), :alon_lovasz, (3, 4))
        @test rand(UInt64) == expected

        # A deterministic family does not consume an explicitly supplied RNG.
        deterministic_rng = MersenneTwister(91)
        control_rng = copy(deterministic_rng)
        deterministic = CatalogUPB.upb(deterministic_rng, :tiles)
        @test !deterministic.rng_used
        @test rand(deterministic_rng, UInt64) == rand(control_rng, UInt64)
    end

    @testset "dimension dispatch and original subsystem order" begin
        deterministic_cases = [
            ((3,), :complete_product_basis, 3),
            ((3, 2), :complete_product_basis, 6),
            ((3, 3), :tiles, 5),
            ((4, 4), :feng_4x4, 8),
            ((2, 2, 2), :generalized_shifts, 4),
            ((2, 2, 3), :feng_2x2x3, 6),
            ((2, 2, 5), :feng_2x2x5, 8),
            ((2, 2, 2, 2), :feng_2x2x2x2, 6),
            ((4, 2, 2, 2), :feng_2x2x2x4, 8),
            ((5, 2, 2, 2, 2), :feng_2x2x2x2x5, 10),
            (ntuple(_ -> 2, 8), :johnston_2_power_8, 11),
            (ntuple(_ -> 2, 12), :johnston_2_power_4k, 16),
            (ntuple(_ -> 2, 5), :generalized_shifts, 6),
            ((4, 6), :chen_johnston_4x6, 10),
        ]
        for (dimensions, family, cardinality) in deterministic_cases
            construction = _wp2_upbc_check_structure(CatalogUPB.upb(dimensions); atol=5e-11)
            @test construction.dimensions == dimensions
            @test construction.family === family
            @test construction.cardinality == cardinality
        end

        random_cases = [
            ((4, 3), :alon_lovasz, 6),
            ((6, 6), :chen_johnston_bipartite, 12),
            ((9, 2, 2), :chen_johnston_4k1, 12),
        ]
        for (index, (dimensions, family, cardinality)) in enumerate(random_cases)
            construction = _wp2_upbc_check_structure(
                CatalogUPB.upb(MersenneTwister(index), dimensions)
            )
            @test construction.dimensions == dimensions
            @test construction.family === family
            @test construction.cardinality == cardinality
            _wp2_upbc_check_unextendible(construction)
        end

        hard = try
            CatalogUPB.upb((2, 2, 2, 2, 2, 2))
            nothing
        catch err
            err
        end
        @test hard isa CatalogUPB.UPBConstructionUnavailable
        @test hard.reason === :known_but_not_in_catalog
        @test hard.known_minimum == 8
        @test occursin("known_but_not_in_catalog", sprint(showerror, hard))

        unknown = try
            CatalogUPB.upb((2, 4, 7))
            nothing
        catch err
            err
        end
        @test unknown isa CatalogUPB.UPBConstructionUnavailable
        @test unknown.reason === :minimum_unknown
        @test unknown.known_minimum === nothing
    end

    @testset "generic precision, validation, and resource limits" begin
        float32 = _wp2_upbc_check_structure(
            CatalogUPB.upb(:tiles; real_type=Float32); atol=2e-5
        )
        @test eltype(float32.global_vectors) == Float32
        @test float32.arithmetic === :float32

        setprecision(BigFloat, 192) do
            big_result = _wp2_upbc_check_structure(
                CatalogUPB.upb(:quad_residue, 3; real_type=BigFloat); atol=big"1e-50"
            )
            @test eltype(big_result.global_vectors) == Complex{BigFloat}
            @test big_result.arithmetic === :arbitrary_precision_floating
            @test big_result.normalization_residual isa BigFloat
        end

        @test_throws ArgumentError CatalogUPB.upb(:generalized_tiles_1, 2)
        @test_throws ArgumentError CatalogUPB.upb(:generalized_tiles_1, 5)
        @test_throws ArgumentError CatalogUPB.upb(:generalized_tiles_2, 4, 3)
        @test_throws ArgumentError CatalogUPB.upb(:quad_residue, 5)
        @test_throws ArgumentError CatalogUPB.upb(:generalized_shifts, 4)
        @test_throws ArgumentError CatalogUPB.upb(:johnston_2_power_4k, 10)
        @test_throws ArgumentError CatalogUPB.upb(:six_parameter, zeros(5))
        @test_throws ArgumentError CatalogUPB.upb(
            :six_parameter, [0, 0.2, 0.3, 0.4, 0.5, 0.6]
        )
        @test_throws ArgumentError CatalogUPB.upb(:not_a_family)
        @test_throws ArgumentError CatalogUPB.upb(())
        @test_throws ArgumentError CatalogUPB.upb((2, 0))
        @test_throws ArgumentError CatalogUPB.upb((2.0, 3))
        @test_throws ArgumentError CatalogUPB.upb(_UPBCZeroBasedVector([3, 3]))
        @test_throws ArgumentError CatalogUPB.upb(:tiles; real_type=Int)
        @test_throws ArgumentError CatalogUPB.upb(:tiles; atol=-1)
        @test_throws ArgumentError CatalogUPB.upb(:tiles; rtol=-1)
        @test_throws ArgumentError CatalogUPB.upb(:tiles; boundary_factor=1)
        @test_throws ArgumentError CatalogUPB.upb(:tiles; max_attempts=0)

        @test_throws CatalogUPB.UPBResourceLimitError CatalogUPB.upb(
            :tiles; max_local_entries=29
        )
        @test_throws CatalogUPB.UPBResourceLimitError CatalogUPB.upb(
            :tiles; max_global_entries=44
        )
        @test_throws CatalogUPB.UPBResourceLimitError CatalogUPB.upb(:tiles; max_work=0)
        @test_throws CatalogUPB.UPBResourceLimitError CatalogUPB.upb(
            MersenneTwister(1), :alon_lovasz, (3, 4); max_minors=0
        )
        limit_error = try
            CatalogUPB.upb(:tiles; max_global_entries=44)
            nothing
        catch err
            err
        end
        @test limit_error.resource === :global_entries
        @test limit_error.required == 45
        @test occursin("global_entries", sprint(showerror, limit_error))
    end

    @testset "one-factorization helper supersession and corrected 4k construction" begin
        for count in (2, 4, 6, 8)
            factorization = CatalogUPB._upbc_one_factorization(count)
            @test size(factorization) == (count - 1, count)
            pair_counts = Dict{Tuple{Int,Int},Int}()
            for round in axes(factorization, 1), pair in 1:(count ÷ 2)
                left = factorization[round, 2pair - 1]
                right = factorization[round, 2pair]
                edge = minmax(left, right)
                pair_counts[edge] = get(pair_counts, edge, 0) + 1
            end
            @test length(pair_counts) == binomial(count, 2)
            @test all(==(1), values(pair_counts))
        end
        labels = [:a, :b, :c, :d]
        labelled = CatalogUPB._upbc_one_factorization(labels)
        @test Set(vec(labelled)) == Set(labels)
        @test_throws ArgumentError CatalogUPB._upbc_one_factorization(3)
        @test_throws ArgumentError CatalogUPB._upbc_one_factorization([1, 1, 2, 3])

        corrected = CatalogUPB.upb(:johnston_2_power_4k, 8)
        @test corrected.orthogonality_residual < 1e-12
        @test _wp2_upbc_qubit_cover_bound(corrected) == (false, 11)
    end

    @testset "public exports and MATLAB compatibility" begin
        native = PublicUPB.upb(:tiles)
        @test native isa PublicUPB.UPBConstruction
        @test native.cardinality == 5
        @test PublicUPB.UPBConstructionUnavailable <: Exception
        @test PublicUPB.UPBResourceLimitError <: Exception

        global_output = PublicUPB.MATLABCompat.UPB("Tiles"; output=:global)
        local_output = PublicUPB.MATLABCompat.UPB(:tiles; output=:local)
        structured = PublicUPB.MATLABCompat.UPB(:tiles; output=:structured)
        @test global_output == native.global_vectors
        @test local_output == native.local_factors
        @test structured isa PublicUPB.UPBConstruction

        single_party = PublicUPB.MATLABCompat.UPB(3, 0; output=:structured)
        @test single_party.dimensions == (3,)
        @test single_party.cardinality == 3

        dimension_result, reference_text = mktemp() do _, reference_output
            result = redirect_stdout(reference_output) do
                PublicUPB.MATLABCompat.UPB((3, 2), 1; output=:structured)
            end
            flush(reference_output)
            seekstart(reference_output)
            return result, read(reference_output, String)
        end
        @test dimension_result.dimensions == (3, 2)
        @test !isempty(reference_text)

        seeded = PublicUPB.MATLABCompat.UPB(
            MersenneTwister(606), (3, 4), 0; output=:structured
        )
        @test seeded.family === :alon_lovasz
        @test seeded.rng_used
        @test_throws ArgumentError PublicUPB.MATLABCompat.UPB((3, 4), 0; output=:structured)
        @test_throws ArgumentError PublicUPB.MATLABCompat.UPB(:tiles; output=:invalid)
        @test_throws ArgumentError PublicUPB.MATLABCompat.UPB(:tiles; output="global")
        @test_throws ArgumentError PublicUPB.MATLABCompat.UPB((3, 2), 0, 1)
        @test_throws ArgumentError PublicUPB.MATLABCompat.UPB((3, 2), 2)
    end
end
