using LinearAlgebra
using SparseArrays

const QETPipeline = QuantumEntanglementTools

@testset "Tier D native entanglement certificate pipeline" begin
    @testset "backend metadata and method configuration" begin
        backends = QETPipeline.available_entanglement_backends()
        @test length(backends) == 1
        absent_backend_error = try
            QETPipeline.detect_entanglement(
                Matrix{Float64}(I, 4, 4) / 4,
                (2, 2),
                QETPipeline.EntanglementDetectionSearch(; max_iteration=1),
            )
            nothing
        catch error
            error
        end
        @test only(backends) isa QETPipeline.NativeEntanglementBackend &&
            absent_backend_error isa ArgumentError &&
            occursin("is not loaded", sprint(showerror, absent_backend_error))

        capabilities = QETPipeline.backend_capabilities(only(backends))
        @test capabilities.name === :native
        @test capabilities.side_effect_free
        @test !capabilities.optional_dependency
        @test :ppt in capabilities.methods
        @test capabilities.conclusions == (:entangled, :separable, :unknown)

        method = QETPipeline.NativePPT(
            systems=(1,), atol=1.0f-7, rtol=2.0f-6, allow_densify=true
        )
        @test method.systems == (1,)
        @test method.atol === 1.0f-7
        @test method.rtol === 2.0f-6
        @test method.allow_densify
        @test_throws ArgumentError QETPipeline.NativePPT(atol=-1)
        @test_throws ArgumentError QETPipeline.NativePPT(rtol=Inf)
    end

    @testset "explicit native PPT dispatch" begin
        bell = ComplexF64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * bell'
        detected = QETPipeline.detect_entanglement(
            bell_density, (2, 2), QETPipeline.NativePPT(atol=1e-12, rtol=0)
        )
        @test detected.status === :entangled
        @test detected.certified
        @test detected.certificate_kind === :negative_partial_transpose_witness
        @test detected.method === :ppt
        @test detected.backend === :native
        @test length(detected.attempts) == 1
        @test detected.evidence.status === QETPipeline.CriterionEntanglementDetected

        maximally_mixed = Matrix{Float64}(I, 4, 4) / 4
        low_dimensional = QETPipeline.detect_entanglement(
            maximally_mixed, (2, 2), QETPipeline.NativePPT(atol=1e-12, rtol=0)
        )
        @test low_dimensional.status === :separable
        @test low_dimensional.certified
        @test low_dimensional.certificate_kind === :ppt_low_dimension_theorem

        mixed_2x3 = Matrix{Float64}(I, 6, 6) / 6
        theorem_2x3 = QETPipeline.detect_entanglement(
            mixed_2x3, (2, 3), QETPipeline.NativePPT(atol=1e-12, rtol=0)
        )
        @test theorem_2x3.status === :separable
        @test theorem_2x3.certified

        mixed_3x3 = Matrix{Float64}(I, 9, 9) / 9
        higher_dimensional = QETPipeline.detect_entanglement(
            mixed_3x3, (3, 3), QETPipeline.NativePPT(atol=1e-12, rtol=0)
        )
        @test higher_dimensional.status === :unknown
        @test !higher_dimensional.certified
        @test higher_dimensional.certificate_kind === nothing

        product_density = Diagonal([1.0, 0.0, 0.0, 0.0])
        boundary = QETPipeline.detect_entanglement(
            Matrix(product_density), (2, 2), QETPipeline.NativePPT(atol=1e-12, rtol=0)
        )
        @test boundary.status === :unknown
        @test occursin("boundary", boundary.message)

        @test_throws ArgumentError QETPipeline.detect_entanglement(
            bell_density, (2, 2), QETPipeline.NativePPT(systems=())
        )
    end

    @testset "certificate-first density analysis" begin
        bell = Float64[1, 0, 0, 1] / sqrt(2)
        bell_density = bell * bell'
        report = QETPipeline.analyze_entanglement(bell_density, (2, 2); atol=1e-12, rtol=0)
        @test report.status === :entangled
        @test report.certified
        @test length(report.attempts) == 1

        maximally_mixed_2 = Matrix{Float64}(I, 4, 4) / 4
        separable = QETPipeline.analyze_entanglement(
            maximally_mixed_2, (2, 2); atol=1e-12, rtol=0
        )
        @test separable.status === :separable
        @test separable.certified
        @test length(separable.attempts) == 1

        maximally_mixed_3 = Matrix{Float64}(I, 9, 9) / 9
        inconclusive = QETPipeline.analyze_entanglement(
            maximally_mixed_3, (3, 3); atol=1e-12, rtol=0
        )
        @test inconclusive.status === :unknown
        @test !inconclusive.certified
        @test [attempt.method for attempt in inconclusive.attempts] == [:ppt, :realignment, :reduction]
        @test all(!attempt.certified for attempt in inconclusive.attempts)

        @test_throws ArgumentError QETPipeline.analyze_entanglement(
            maximally_mixed_2, (2, 2); strategy=:guess
        )

        sparse_state = sparse(maximally_mixed_2)
        @test_throws ArgumentError QETPipeline.analyze_entanglement(sparse_state, (2, 2))
        sparse_report = QETPipeline.analyze_entanglement(
            sparse_state, (2, 2); atol=1e-12, rtol=0, allow_densify=true
        )
        @test sparse_report.status === :separable
    end

    @testset "pure-state Schmidt certificates" begin
        bell = ComplexF64[1, 0, 0, im] / sqrt(2)
        entangled = QETPipeline.analyze_entanglement(bell, (2, 2); atol=0, rtol=1e-12)
        @test entangled.status === :entangled
        @test entangled.certified
        @test entangled.certificate_kind === :pure_state_schmidt_rank
        @test entangled.evidence.rank == 2

        product = ComplexF64[1, 0, 0, 0]
        separable = QETPipeline.analyze_entanglement(product, (2, 2))
        @test separable.status === :separable
        @test separable.certified
        @test separable.certificate_kind === :pure_product_decomposition
        @test separable.evidence.rank == 1

        weakly_entangled = ComplexF64[sqrt(1 - 1e-20), 0, 0, 1e-10]
        boundary = QETPipeline.analyze_entanglement(
            weakly_entangled, (2, 2); atol=0, rtol=sqrt(eps(Float64))
        )
        @test boundary.status === :unknown
        @test !boundary.certified
        @test boundary.certificate_kind === nothing
        @test boundary.evidence.rank == 1

        oversized_tolerance = QETPipeline.analyze_entanglement(
            product, (2, 2); atol=0, rtol=2
        )
        @test oversized_tolerance.status === :unknown
        @test !oversized_tolerance.certified
        @test oversized_tolerance.evidence.rank == 0

        sparse_product = sparsevec([1], ComplexF64[1], 6)
        @test_throws ArgumentError QETPipeline.analyze_entanglement(sparse_product, (2, 3))
        sparse_separable = QETPipeline.analyze_entanglement(
            sparse_product, (2, 3); allow_densify=true
        )
        @test sparse_separable.status === :separable

        @test_throws ArgumentError QETPipeline.analyze_entanglement(0.9product, (2, 2))
        @test_throws DimensionMismatch QETPipeline.analyze_entanglement(product, (2, 3))
    end

    @testset "result invariants and display" begin
        @test_throws ArgumentError QETPipeline.EntanglementAttempt(
            :test, :native, :invalid, false, nothing, nothing, "invalid"
        )
        @test_throws ArgumentError QETPipeline.EntanglementAttempt(
            :test, :native, :entangled, true, nothing, nothing, "missing certificate"
        )
        @test_throws ArgumentError QETPipeline.EntanglementReport(
            :unknown,
            true,
            :impossible,
            :test,
            :native,
            nothing,
            QETPipeline.EntanglementAttempt[],
            "invalid",
        )

        report = QETPipeline.analyze_entanglement(Matrix{Float64}(I, 4, 4) / 4, (2, 2))
        @test occursin("EntanglementReport", sprint(show, report))
        @test occursin("EntanglementAttempt", sprint(show, only(report.attempts)))
    end
end
