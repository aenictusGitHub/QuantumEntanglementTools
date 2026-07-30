using Test

const UPBCompat = QuantumEntanglementTools.MATLABCompat

@testset "WP2 minimum UPB size" begin
    @testset "theorem branches" begin
        two_by_three = minimum_upb_size((2, 3))
        @test two_by_three isa MinimumUPBSizeResult
        @test two_by_three.status === :known
        @test two_by_three.size == 6
        @test two_by_three.lower_bound == 4
        @test two_by_three.reference_key === :divincenzo_mor_shor_smolin_terhal_2003
        @test startswith(two_by_three.reference_url, "https://")

        all_odd = minimum_upb_size((3, 3, 3))
        @test all_odd.size == 7
        @test all_odd.reference_key === :alon_lovasz_2001
        @test minimum_upb_size((2, 3, 3)).size == 6
        @test minimum_upb_size((4, 6)).size == 10
        @test minimum_upb_size((6, 4)).dimensions == (4, 6)
        @test minimum_upb_size((4, 6)).reference_key === :chen_johnston_2015

        four_qubits = minimum_upb_size(fill(2, 4))
        @test four_qubits.size == 6
        @test four_qubits.reference_key === :feng_2006
        @test minimum_upb_size(fill(2, 6)).size == 8
        @test minimum_upb_size(fill(2, 8)).size == 11
        @test minimum_upb_size(fill(2, 12)).size == 16
        @test minimum_upb_size(fill(2, 12)).reference_key === :johnston_2013

        @test minimum_upb_size((2, 2, 3)).size == 6
        @test minimum_upb_size((5, 2, 2)).size == 8
        @test minimum_upb_size((2, 2, 9)).size == 12
        @test minimum_upb_size((2, 2, 9)).reference_key === :chen_johnston_2015
    end

    @testset "explicit unknown" begin
        unknown = minimum_upb_size((2, 3, 4))
        @test unknown.status === :unknown
        @test unknown.size === nothing
        @test unknown.lower_bound == 7
        @test unknown.reference_key === :unresolved_by_reviewed_table
        @test unknown.reference === nothing
        @test unknown.reference_url === nothing
        @test occursin("does not determine", unknown.message)
        @test occursin("status=unknown", sprint(show, unknown))
    end

    @testset "compatibility behavior" begin
        @test UPBCompat.MinUPBSize([3, 3], 0) == 5
        value, printed_reference = mktemp() do _, output
            value = redirect_stdout(output) do
                UPBCompat.MinUPBSize([2, 3])
            end
            flush(output)
            seekstart(output)
            return value, read(output, String)
        end
        @test value == 6
        @test occursin("DiVincenzo", printed_reference)
        error = try
            UPBCompat.MinUPBSize([2, 3, 4], 0)
            nothing
        catch err
            err
        end
        @test error isa DomainError
        @test error.val isa MinimumUPBSizeResult
        @test error.val.status === :unknown
        @test_throws ArgumentError UPBCompat.MinUPBSize([3, 3], 2)
    end

    @testset "validation and checked arithmetic" begin
        @test_throws ArgumentError minimum_upb_size(3)
        @test_throws ArgumentError minimum_upb_size((3,))
        @test_throws ArgumentError minimum_upb_size((1, 3))
        @test_throws ArgumentError minimum_upb_size((true, 3))
        @test_throws ArgumentError minimum_upb_size((2.0, 3))
        @test_throws ArgumentError minimum_upb_size((typemax(Int), 3, 3))
        @test_throws ArgumentError minimum_upb_size([2, 0, 3])
    end
end
