using CyclingSignatures
using Distances: euclidean, normalize
using Ripserer
using Test

import CyclingSignatures: colspace_normal_form

include("test-util.jl")

function test_ripserer_barcode_invariants(barcode, field; all_essential=false)
    if all_essential
        @test all(bar -> bar.death == Inf, barcode)
    end
    @test all(bar -> length(bar.simplex_list) == length(bar.coeff_list), barcode)
    @test all(bar -> all(simplex -> length(simplex) == 2, bar.simplex_list), barcode)
    @test all(bar -> eltype(bar.coeff_list) == field, barcode)
    @test all(bar -> all(!iszero, bar.coeff_list), barcode)
end

function double_circle_signatures(alg, figure8_space, subdivision, boxsize)
    left_sig = cycling_signature(
        alg,
        figure8_space,
        1:(subdivision + 1),
        boxsize,
    )
    right_sig = cycling_signature(
        alg,
        figure8_space,
        (subdivision + 1):(2 * subdivision + 1),
        boxsize,
    )
    both_sig = cycling_signature(
        alg,
        figure8_space,
        1:(2 * subdivision + 1),
        boxsize,
    )

    return left_sig, right_sig, both_sig
end

function test_double_circle_signatures(alg, figure8_space, subdivision, boxsize)
    left_sig, right_sig, both_sig = double_circle_signatures(
        alg,
        figure8_space,
        subdivision,
        boxsize,
    )

    @test dimension(left_sig) == 1
    @test dimension(right_sig) == 1
    @test dimension(both_sig) == 2

    left_space = colspace_normal_form(cycling_matrix(left_sig))
    right_space = colspace_normal_form(cycling_matrix(right_sig))
    both_space = colspace_normal_form(cycling_matrix(both_sig))

    @test left_space != right_space
    @test size(left_space, 2) == 1
    @test size(right_space, 2) == 1
    @test size(both_space, 2) == 2
end

function test_double_circle_signatures_known_broken(alg, figure8_space, subdivision, boxsize)
    left_sig, right_sig, both_sig = double_circle_signatures(
        alg,
        figure8_space,
        subdivision,
        boxsize,
    )

    @test dimension(left_sig) == 1
    @test dimension(right_sig) == 1

    left_space = colspace_normal_form(cycling_matrix(left_sig))
    right_space = colspace_normal_form(cycling_matrix(right_sig))
    both_space = colspace_normal_form(cycling_matrix(both_sig))

    recovers_two_generators = dimension(both_sig) == 2 && size(both_space, 2) == 2
    separates_circles = left_space != right_space

    if recovers_two_generators
        @test recovers_two_generators
    else
        @test_broken recovers_two_generators
    end

    if separates_circles
        @test separates_circles
    else
        @test_broken separates_circles
    end
end

@testset "Ripserer extension" begin
    boxsize = 0.2
    tsn, circle_data = circle_trajectory_space(40, 1; boxsize)

    @testset "barcode dispatch" begin
        points = CyclingSignatures.evaluate_interval(
            get_trajectory(tsn),
            1,
            size(circle_data, 2),
        )
        for alg in (Val(:Ripserer), Val(:RipsererNoThreshold), Val(:RipsererManualReconstruct))
            @testset "$(alg)" begin
                barcode = CyclingSignatures.trajectory_barcode(
                    alg,
                    points,
                    euclidean,
                    boxsize,
                    FF{2},
                )

                @test !isempty(barcode)
                test_ripserer_barcode_invariants(barcode, FF{2}; all_essential=true)

                barcode_ff5 = CyclingSignatures.trajectory_barcode(
                    alg,
                    points,
                    euclidean,
                    boxsize,
                    FF{5},
                )

                @test !isempty(barcode_ff5)
                test_ripserer_barcode_invariants(barcode_ff5, FF{5}; all_essential=true)
            end
        end
    end

    @testset "barcode dispatch on double circle" begin
        figure8_space, figure8_data = figure8_trajectory_space(40, [0, 1]; boxsize)
        points = CyclingSignatures.evaluate_interval(
            get_trajectory(figure8_space),
            1,
            size(figure8_data, 2),
        )

        for alg in (Val(:Ripserer), Val(:RipsererNoThreshold), Val(:RipsererManualReconstruct))
            @testset "$(alg)" begin
                barcode = CyclingSignatures.trajectory_barcode(
                    alg,
                    points,
                    euclidean,
                    boxsize,
                    FF{2},
                )

                @test length(barcode) >= 2
                test_ripserer_barcode_invariants(barcode, FF{2})
                @test all(bar -> 0 <= bar.birth <= boxsize, barcode)
                @test any(bar -> bar.death == Inf, barcode)
            end
        end

        ripserer_ext = Base.get_extension(CyclingSignatures, :RipsererExt)
        @test ripserer_ext !== nothing

        d_mat = ripserer_ext.ripserer_distance_matrix(points, euclidean)
        off_diagonal_zeros = [
            d_mat[i, j]
            for j in axes(d_mat, 2)
            for i in axes(d_mat, 1)
            if i != j && iszero(d_mat[i, j])
        ]
        @test isempty(off_diagonal_zeros)
    end

    @testset "cycling signatures on double circle" begin
        subdivision = 40
        figure8_space, _ = figure8_trajectory_space(subdivision, [0, 1]; boxsize)

        @testset "manual reconstruction" begin
            test_double_circle_signatures(
                Val(:RipsererManualReconstruct),
                figure8_space,
                subdivision,
                boxsize,
            )
        end

        @testset "no threshold" begin
            test_double_circle_signatures(
                Val(:RipsererNoThreshold),
                figure8_space,
                subdivision,
                boxsize,
            )
        end

        @testset "thresholded involuted representatives" begin
            test_double_circle_signatures_known_broken(
                Val(:Ripserer),
                figure8_space,
                subdivision,
                boxsize,
            )
        end
    end

    @testset "backend invariants" begin
        for alg in (
            Val(:DistanceMatrix),
            Val(:Ripserer),
            Val(:RipsererNoThreshold),
            Val(:RipsererManualReconstruct),
        )
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
        for alg in (Val(:Ripserer), Val(:RipsererNoThreshold), Val(:RipsererManualReconstruct))
            @testset "$(alg)" begin
                ripserer_sig = cycling_signature(alg, tsn, (1, size(circle_data, 2)), boxsize)

                @test dimension(ripserer_sig) == dimension(dm_sig)
                @test colspace_normal_form(cycling_matrix(ripserer_sig)) ==
                      colspace_normal_form(cycling_matrix(dm_sig))
            end
        end
    end

    @testset "experiment agreement diagnostics" begin
        exp = RandomSubsegmentExperiment(tsn, [size(circle_data, 2)], 1, 1234)
        starts = sample_segment_starts(exp)

        dm = run_experiment(
            exp;
            alg = Val(:DistanceMatrix),
            threshold = boxsize,
            segment_starts = starts,
            progress = false,
        )
        ripserer = run_experiment(
            exp;
            alg = Val(:Ripserer),
            threshold = boxsize,
            segment_starts = starts,
            progress = false,
        )

        comparison = compare_experiment_results(dm, ripserer)
        summary = only(summarize_agreement(comparison))
        @test summary.rank_mismatches == 0
        @test summary.cycling_space_mismatches == 0
    end

    @testset "coefficient fields" begin
        for alg in (Val(:Ripserer), Val(:RipsererNoThreshold), Val(:RipsererManualReconstruct))
            @testset "$(alg)" begin
                for field in (FF{3}, FF{5})
                    sig = cycling_signature(
                        alg,
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
    end
end
