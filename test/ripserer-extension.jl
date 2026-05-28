using CyclingSignatures
using Distances: euclidean, normalize
using Ripserer
using Test

import CyclingSignatures: colspace_normal_form

include("test-util.jl")

@testset "Ripserer extension" begin
    boxsize = 0.2
    tsn, circle_data = circle_trajectory_space(40, 1; boxsize)

    @testset "barcode dispatch" begin
        points = CyclingSignatures.evaluate_interval(
            get_trajectory(tsn),
            1,
            size(circle_data, 2),
        )
        barcode = CyclingSignatures.trajectory_barcode(
            Val(:Ripserer),
            points,
            euclidean,
            boxsize,
            FF{2},
        )

        @test !isempty(barcode)
        @test all(bar -> bar.death == Inf, barcode)
        @test all(bar -> length(bar.simplex_list) == length(bar.coeff_list), barcode)
        @test all(bar -> all(simplex -> length(simplex) == 2, bar.simplex_list), barcode)
        @test all(bar -> eltype(bar.coeff_list) == FF{2}, barcode)
    end

    @testset "backend invariants" begin
        for alg in (Val(:DistanceMatrix), Val(:Ripserer))
            @testset "$(alg)" begin
                sig = cycling_signature(alg, tsn, (1, size(circle_data, 2)), boxsize)

                @test dimension(sig) == 1
                @test size(cycling_matrix(sig), 1) == betti_1(tsn)
                @test size(cycling_matrix(sig), 2) == length(birth_vector(sig))
                @test !all(iszero, cycling_matrix(sig))

                empty_sig = cycling_signature(alg, tsn, (1, 3), boxsize)
                @test dimension(empty_sig) == 0
                @test size(cycling_matrix(empty_sig)) == (betti_1(tsn), 0)
                @test birth_vector(empty_sig) == Float64[]
            end
        end
    end

    @testset "same cycling space as distance matrix" begin
        dm_sig = cycling_signature(Val(:DistanceMatrix), tsn, (1, size(circle_data, 2)), boxsize)
        ripserer_sig = cycling_signature(Val(:Ripserer), tsn, (1, size(circle_data, 2)), boxsize)

        @test dimension(ripserer_sig) == dimension(dm_sig)
        @test colspace_normal_form(cycling_matrix(ripserer_sig)) ==
              colspace_normal_form(cycling_matrix(dm_sig))
    end

    @testset "coefficient fields" begin
        for field in (FF{3}, FF{5})
            sig = cycling_signature(
                Val(:Ripserer),
                tsn,
                (1, size(circle_data, 2)),
                boxsize;
                field,
            )

            @test dimension(sig) == 1
            @test eltype(cycling_matrix(sig)) == field
            @test !all(iszero, cycling_matrix(sig))
        end
    end
end
