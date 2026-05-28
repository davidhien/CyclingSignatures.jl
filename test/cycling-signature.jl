using CyclingSignatures
using Test
using Distances: euclidean, normalize
include("test-util.jl")

@testset "cycling signature w/o utb" begin
    @testset "simple circle data" begin
        boxsize = 0.2
        tsn, circle_data = circle_trajectory_space(40, 2; boxsize)

        test_dimension_profile(
            Val(:DistanceMatrix),
            tsn,
            boxsize,
            [1:39 => 0, 40:size(circle_data, 2) => 1],
        )
    end

    @testset "figure 8" begin
        boxsize = 0.2
        tsn, _ = figure8_trajectory_space(40, [0, 1, 0, 0, 1, 1]; boxsize)

        test_dimension_profile(
            Val(:DistanceMatrix),
            tsn,
            boxsize,
            [1:39 => 0, 40:78 => 1, 79:100 => 2],
        )

        # check generators for first and second turn
        sig1 = cycling_signature(Val(:DistanceMatrix), tsn, (1,45), boxsize)
        sig2 = cycling_signature(Val(:DistanceMatrix), tsn, (37,85), boxsize)
        sig3 = cycling_signature(Val(:DistanceMatrix), tsn, (75,125), boxsize)

        @test cycling_matrix(sig1) != cycling_matrix(sig2)
        @test cycling_matrix(sig1) == cycling_matrix(sig3)
    end
end


@testset "cycling signature w utb" begin
    @testset "simple circle utb" begin
        boxsize = 0.2
        tsn, utb_circle_data = utb_circle_trajectory_space(41, 2; boxsize)

        test_dimension_profile(
            Val(:DistanceMatrix),
            tsn,
            boxsize,
            [1:40 => 0, 41:size(utb_circle_data, 2) => 1],
        )
    end

    @testset "figure 8" begin
        boxsize = 0.2
        tsn, _ = utb_figure8_trajectory_space(80, [0, 1, 0, 0, 1, 1]; boxsize)

        test_dimension_profile(
            Val(:DistanceMatrix),
            tsn,
            boxsize,
            [1:78 => 0, 79:154 => 1, 155:200 => 2],
        )

        sig1 = cycling_signature(Val(:DistanceMatrix), tsn, (  1, 85), boxsize)
        sig2 = cycling_signature(Val(:DistanceMatrix), tsn, ( 75,165), boxsize)
        sig3 = cycling_signature(Val(:DistanceMatrix), tsn, (155,245), boxsize)

        @test cycling_matrix(sig1) != cycling_matrix(sig2)
        @test cycling_matrix(sig1) == cycling_matrix(sig3)
    end
end

@testset "unit tests" begin
    @testset "CyclingSignature" begin
        F = FF{2}
        birth_vector = [1.0, 2.0, 3.0]
        cycling_matrix = [1 0 1; 0 1 2]
        cs = CyclingSignature(cycling_matrix, birth_vector)
        @test cs.cycling_matrix == cycling_matrix
        @test cs.birth_vector == birth_vector

        # Test with incompatible sizes
        @test_throws ArgumentError CyclingSignature(cycling_matrix, [1.0, 2.0])
        @test dimension(cs) == 3

        dim_fct = CyclingSignatures.dimension_function(cs)
        @test dim_fct(0.5) == 0
        @test dim_fct(1.5) == 1
        @test dim_fct(2.5) == 2
        @test dim_fct(3.5) == 3
    end

    @testset "TrajectorySpace" begin
        # generate traj
        traj = Matrix(collect(0:0.1:1)'[:, :])
        rt = RefinedEquidistantTrajectory(traj, [1, 6, 12])

        # generate comparison space
        comp_space = cubical_vr_comparison_space_via_cover(traj, 0.2)

        # generate distance
        d = euclidean

        tsn = TrajectorySpace(rt, comp_space, d)

        @test get_trajectory(tsn) == rt
        @test get_comparison_space(tsn) == comp_space
        @test get_metric(tsn) == d
    end

    @testset "cycling_signature defaults and empty signatures" begin
        boxsize = 0.2
        tsn, circle_data = circle_trajectory_space(40, 1; boxsize, flt_max_heuristic=boxsize)
        sig_default = cycling_signature(tsn, (1, size(circle_data, 2)))
        sig_explicit = cycling_signature(Val(:DistanceMatrix), tsn, (1, size(circle_data, 2)), boxsize)

        @test cycling_matrix(sig_default) == cycling_matrix(sig_explicit)
        @test birth_vector(sig_default) == birth_vector(sig_explicit)

        tsn_without_heuristic, _ = circle_trajectory_space(40, 1; boxsize)
        @test_throws ArgumentError cycling_signature(tsn_without_heuristic, (1, size(circle_data, 2)))

        empty_sig = cycling_signature(Val(:DistanceMatrix), tsn_without_heuristic, (1, 1), boxsize)
        @test dimension(empty_sig) == 0
        @test size(cycling_matrix(empty_sig)) == (betti_1(tsn_without_heuristic), 0)
        @test birth_vector(empty_sig) == Float64[]
    end

    @testset "coefficient fields" begin
        boxsize = 0.2
        tsn, circle_data = circle_trajectory_space(40, 1; boxsize)

        for field in (FF{3}, FF{5})
            sig = cycling_signature(
                Val(:DistanceMatrix),
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

    @testset "barcode to cycling vectors" begin
        boxsize = 0.2
        tsn, circle_data = circle_trajectory_space(40, 1; boxsize)
        comp_space = get_comparison_space(tsn)

        edges, coeffs = CyclingSignatures.curve_cycle(1, size(circle_data, 2); F=FF{2})
        infinite_bar = CyclingSignatures.TrajectoryBar(0.2, Inf, edges, coeffs)
        finite_bar = CyclingSignatures.TrajectoryBar(0.1, 0.3, edges, coeffs)

        births, vectors = CyclingSignatures.births_cycling_vectors_from_trajectory_barcode(
            circle_data,
            [finite_bar, infinite_bar],
            comp_space,
        )

        @test births == [0.2]
        @test length(vectors) == 1
        @test !all(iszero, vectors[1])

        dependent_births, dependent_vectors =
            CyclingSignatures.births_cycling_vectors_from_trajectory_barcode(
                circle_data,
                [infinite_bar, CyclingSignatures.TrajectoryBar(0.25, Inf, edges, coeffs)],
                comp_space,
            )

        @test dependent_births == [0.2]
        @test length(dependent_vectors) == 1
    end
end
