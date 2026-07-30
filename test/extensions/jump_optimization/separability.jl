using Hypatia
using JuMP
using LinearAlgebra
using QuantumEntanglementTools
using SCS
using Test

const SeparabilityOptQET = QuantumEntanglementTools

if !isdefined(SeparabilityOptQET, :StateDiscriminationResult)
    Base.include(
        SeparabilityOptQET,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "optimization", "state_discrimination.jl"
        ),
    )
end
if !isdefined(SeparabilityOptQET, :LocalDistinguishabilityResult)
    Base.include(
        SeparabilityOptQET,
        joinpath(
            @__DIR__,
            "..",
            "..",
            "..",
            "src",
            "entanglement",
            "separability_optimization.jl",
        ),
    )
end

const SeparabilityOptExt = Base.get_extension(
    SeparabilityOptQET, :QuantumEntanglementToolsJuMPExt
)
isnothing(SeparabilityOptExt) && error("JuMP extension is not loaded")
if !isdefined(SeparabilityOptExt, :_SEPARABILITY_OPTIMIZATION_EXTENSION_LOADED)
    Base.include(
        SeparabilityOptExt,
        joinpath(@__DIR__, "..", "..", "..", "ext", "separability_optimization.jl"),
    )
end
const SeparabilityOptCompat = SeparabilityOptQET.MATLABCompat
if !isdefined(SeparabilityOptCompat, :_SEPARABILITY_OPTIMIZATION_COMPAT_LOADED)
    Core.eval(
        SeparabilityOptCompat,
        :(using ..QuantumEntanglementTools:
            local_distinguishability, is_separable, upb_sep_distinguishable),
    )
    Base.include(
        SeparabilityOptCompat,
        joinpath(
            @__DIR__, "..", "..", "..", "src", "compat", "separability_optimization.jl"
        ),
    )
end

function separability_hypatia_backend(; options=NamedTuple(), kwargs...)
    return SeparabilityOptQET.JuMPBackend(
        Hypatia.Optimizer;
        optimizer_options=options,
        optimizer_name="Hypatia",
        optimizer_version=Base.pkgversion(Hypatia),
        allow_densify=true,
        atol=5.0e-7,
        rtol=5.0e-7,
        kwargs...,
    )
end

function separability_scs_backend(; kwargs...)
    return SeparabilityOptQET.JuMPBackend(
        SCS.Optimizer;
        optimizer_options=(eps_abs=1.0e-7, eps_rel=1.0e-7),
        optimizer_name="SCS",
        optimizer_version=Base.pkgversion(SCS),
        allow_densify=true,
        atol=1.0e-5,
        rtol=1.0e-5,
        kwargs...,
    )
end

function separability_opt_bell()
    return [
        0.5 0.0 0.0 0.5
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
        0.5 0.0 0.0 0.5
    ]
end

@testset "WP6 separability optional optimization" begin
    product_columns = hcat([1.0, 0, 0, 0], [0.0, 0, 0, 1])

    @testset "local-discrimination outer hierarchy" begin
        hypatia = SeparabilityOptQET.local_distinguishability(
            product_columns, (2, 2); backend=separability_hypatia_backend()
        )
        positional = SeparabilityOptQET.local_distinguishability(
            product_columns, (2, 2), separability_hypatia_backend()
        )
        @test hypatia.status === SeparabilityOptQET.LocalDistinguishabilitySolverOptimal
        @test positional.status === hypatia.status
        @test positional.relaxation_value ≈ hypatia.relaxation_value atol = 1.0e-7
        @test hypatia.optimization_result.status === SeparabilityOptQET.OptimizationOptimal
        @test hypatia.relaxation_value ≈ 1 atol = 5.0e-7
        @test hypatia.relaxation_lower_bound ≈ 1 atol = 5.0e-7
        @test hypatia.relaxation_upper_bound ≈ 1 atol = 5.0e-7
        @test hypatia.separable_lower_bound == 1 / 2
        @test hypatia.separable_upper_bound ≈ 1 atol = 5.0e-7
        @test hypatia.residuals.valid
        @test hypatia.measurement !== nothing
        @test length(hypatia.measurement) == 2
        @test !hypatia.certified
        @test occursin("outer hierarchy", first(hypatia.warnings))

        legacy = SeparabilityOptCompat.LocalDistinguishability(
            product_columns,
            nothing,
            (2, 2),
            2,
            1,
            1,
            eps(Float64)^(1 / 4);
            structured=false,
            backend=separability_hypatia_backend(),
        )
        @test legacy.dist ≈ 1 atol = 5.0e-7
        @test length(legacy.meas) == 2
        @test size(legacy.dual_sol) == (4, 4)
        @test ishermitian(legacy.dual_sol)

        identical_columns = hcat(product_columns[:, 1], product_columns[:, 1])
        identical = SeparabilityOptQET.local_distinguishability(
            identical_columns, (2, 2); backend=separability_hypatia_backend()
        )
        @test identical.status === SeparabilityOptQET.LocalDistinguishabilitySolverOptimal
        @test identical.relaxation_value ≈ 1 / 2 atol = 5.0e-7
        @test identical.separable_lower_bound == 1 / 2
        @test identical.separable_upper_bound ≈ 1 / 2 atol = 5.0e-7

        plus_plus = fill(inv(2.0), 4)
        nonorthogonal = hcat(product_columns[:, 1], plus_plus)
        global_relaxation = SeparabilityOptQET.local_distinguishability(
            nonorthogonal,
            (2, 2);
            order=1,
            ppt=false,
            bosonic=false,
            backend=separability_hypatia_backend(),
        )
        helstrom_value = (1 + sqrt(1 - abs2(dot(nonorthogonal[:, 1], plus_plus)))) / 2
        @test global_relaxation.status ===
            SeparabilityOptQET.LocalDistinguishabilitySolverOptimal
        @test global_relaxation.relaxation_value ≈ helstrom_value atol = 5.0e-7
        @test global_relaxation.residuals.valid
        @test global_relaxation.problem.order == 1

        scs = SeparabilityOptQET.local_distinguishability(
            product_columns, (2, 2); backend=separability_scs_backend()
        )
        @test scs.status in (
            SeparabilityOptQET.LocalDistinguishabilitySolverOptimal,
            SeparabilityOptQET.LocalDistinguishabilitySolverFeasible,
        )
        @test scs.relaxation_lower_bound ≈ 1 atol = 5.0e-5
        @test scs.residuals.valid
        @test scs.optimization_result.optimizer.reported_optimizer_name == "SCS"
    end

    @testset "outer and inner separability certificates" begin
        bell = separability_opt_bell()
        outer = SeparabilityOptQET.is_separable(
            bell,
            (2, 2),
            separability_hypatia_backend();
            strategies=(:symmetric_extension,),
            extension_ppt=false,
            extension_bosonic=false,
        )
        @test outer.status === :entangled
        @test outer.certified
        @test outer.certificate_kind === :validated_symmetric_extension_separator
        @test outer.evidence.verdict === false
        @test outer.evidence.witness.entanglement_witness
        @test outer.evidence.witness.separator_validated

        pauli_x = [0.0 1.0; 1.0 0.0]
        interior = Matrix{Float64}(I, 4, 4) / 4 + 0.01 * kron(pauli_x, pauli_x) / 4
        inner = SeparabilityOptQET.is_separable(
            interior,
            (2, 2),
            separability_hypatia_backend();
            strategies=(:symmetric_inner_extension,),
            extension_ppt=true,
        )
        @test inner.status === :separable
        @test inner.certified
        @test inner.certificate_kind === :validated_inner_separable_cone_membership
        @test inner.evidence.verdict === true
        @test inner.evidence.residuals.valid
    end

    @testset "UPB reconstruction preserves numerical status" begin
        tiles = SeparabilityOptQET.upb(:tiles)
        result = SeparabilityOptQET.upb_sep_distinguishable(
            tiles.local_factors, separability_hypatia_backend()
        )
        @test result.status === :numerically_feasible
        @test result.feasibility === :numerically_feasible
        @test result.separably_distinguishable === nothing
        @test !result.certified
        @test result.certificate_kind === nothing
        @test result.optimization_result.status === SeparabilityOptQET.OptimizationOptimal
        @test result.coefficients !== nothing
        @test result.minimum_coefficient >= -result.tolerance
        @test result.reconstruction_residual <= result.tolerance
        @test result.candidates_generated == 30
        @test_throws Exception (result.coefficients[1] = 0)
        @test_throws Exception (result.reconstruction[1, 1] = 0)
        @test_throws DomainError SeparabilityOptCompat.UPBSepDistinguishable(
            tiles.local_factors...; structured=false, backend=separability_hypatia_backend()
        )

        impossible_program = SeparabilityOptQET.SemidefiniteProgram(
            :upb_farkas_sign_probe,
            :feasibility,
            1,
            SeparabilityOptQET.AffineScalar(0.0, [0.0]);
            equalities=[
                SeparabilityOptQET.AffineEquality(
                    SeparabilityOptQET.AffineScalar(1.0, [1.0]), :impossible
                ),
            ],
            intervals=[
                SeparabilityOptQET.AffineInterval(
                    SeparabilityOptQET.AffineScalar(0.0, [1.0]), 0.0, nothing, :nonnegative
                ),
            ],
        )
        infeasible = SeparabilityOptQET.solve_optimization(
            impossible_program, separability_hypatia_backend()
        )
        @test infeasible.status === SeparabilityOptQET.OptimizationInfeasible
        certificate = SeparabilityOptQET._upbsep_farkas_certificate(
            (program=impossible_program, dimensions=(1, 1)), infeasible, 1.0e-6
        )
        @test certificate.valid
        @test certificate.orientation in (:positive, :negative)
        @test minimum(abs, certificate.cone_coefficients) > certificate.tolerance
        @test abs(certificate.constant) > certificate.tolerance
        @test occursin("computed floating replacement cone", certificate.message)
    end

    @testset "backend limits and failures stay structured" begin
        limited = SeparabilityOptQET.local_distinguishability(
            product_columns,
            (2, 2);
            backend=separability_hypatia_backend(options=(iter_limit=0,)),
        )
        @test limited.status === SeparabilityOptQET.LocalDistinguishabilityResourceLimit
        @test limited.optimization_result.status === SeparabilityOptQET.OptimizationLimit
        @test limited.relaxation_value === nothing
        @test limited.measurement === nothing

        malformed_backend = SeparabilityOptQET.JuMPBackend(
            () -> error("intentional separability optimizer factory failure");
            optimizer_name="intentional malformed factory",
            allow_densify=true,
        )
        malformed = SeparabilityOptQET.local_distinguishability(
            product_columns, (2, 2); backend=malformed_backend
        )
        @test malformed.status === SeparabilityOptQET.LocalDistinguishabilityBackendFailure
        @test malformed.optimization_result.status ===
            SeparabilityOptQET.OptimizationMalformedBackend
        @test malformed.relaxation_value === nothing
        @test malformed.measurement === nothing

        malformed_upb = SeparabilityOptQET.upb_sep_distinguishable(
            SeparabilityOptQET.upb(:tiles).local_factors, malformed_backend
        )
        @test malformed_upb.status === :backend_failure
        @test malformed_upb.feasibility === :unknown
        @test malformed_upb.separably_distinguishable === nothing
        @test malformed_upb.optimization_result.status ===
            SeparabilityOptQET.OptimizationMalformedBackend
    end
end
