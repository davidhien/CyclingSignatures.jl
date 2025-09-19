using CyclingSignatures
using Test
using Distances: euclidean, normalize
include("test-util.jl")

@testset "cycling_signature w/o utb" begin
    @testset "simple circle data" begin
        # simple circle data testset
        circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(40, 2), dims=2)

        # generate comparison space
        boxsize = 0.2
        comp_space = cubical_vr_comparison_space_via_cover(circle_data, boxsize)
        traj = RefinedEquidistantTrajectory(circle_data)

        tsn = TrajectorySpaceNew(traj, comp_space, euclidean)

        for i = 1:size(circle_data, 2)
            sig = cycling_signature(Val(:DistanceMatrix), tsn, (1,i), boxsize)
            if i <= 39
                @test dimension(sig) == 0
            else
                @test dimension(sig) == 1
            end
        end
    end
    @testset "figure 8" begin
        circle_data = double_circle_time_series(40, [0,1,0,0,1,1])
        boxsize = 0.2
        comp_space = cubical_vr_comparison_space_via_cover(circle_data, boxsize)
        traj = RefinedEquidistantTrajectory(circle_data)

        tsn = TrajectorySpaceNew(traj, comp_space, euclidean)

        # check for dimension 0, 1, 2 at the start
        for i = 1:100
            sig = cycling_signature(Val(:DistanceMatrix), tsn, (1,i), boxsize)
            if i <= 39
                @test dimension(sig) == 0
            elseif i <= 78
                @test dimension(sig) == 1
            else
                @test dimension(sig) == 2
            end
        end

        # check generators for first and second turn
        sig1 = cycling_signature(Val(:DistanceMatrix), tsn, (1,45), boxsize)
        sig2 = cycling_signature(Val(:DistanceMatrix), tsn, (37,85), boxsize)
        sig3 = cycling_signature(Val(:DistanceMatrix), tsn, (75,125), boxsize)

        @test cycling_matrix(sig1) != cycling_matrix(sig2)
        @test cycling_matrix(sig1) == cycling_matrix(sig3)
    end
end


@testset "cycling-signature w utb" begin
        @testset "simple circle utb" begin
        # l_infty circle data
        circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(41, 2), dims=2)
        # discrete derivative
        circle_data_dd = circle_data[:, 2:end] - circle_data[:, 1:end-1]

        # data in utb
        utb_circle_data = [circle_data[:, 1:end-1]; circle_data_dd]

        # generate comparison space
        boxsize = 0.2
        sb_radius = 1
        comp_space = sb_cubical_vr_comparison_space_via_cover(utb_circle_data, boxsize, sb_radius)
        traj = RefinedEquidistantTrajectory(utb_circle_data)

        tsn = TrajectorySpaceNew(traj, comp_space, DynamicDistance(2, sb_radius))

        for i = 1:size(utb_circle_data, 2)
            sig = cycling_signature(Val(:DistanceMatrix), tsn, (1,i), boxsize)
            if i <= 40
                @test dimension(sig) == 0
            else
                @test dimension(sig) == 1
            end
        end

    end
    @testset "figure 8" begin
        circle_data = double_circle_time_series(80, [0,1,0,0,1,1])

        # discrete derivative
        circle_data_dd = circle_data[:, 2:end] - circle_data[:, 1:end-1]

        # data in utb
        utb_circle_data = [circle_data[:, 1:end-1]; circle_data_dd]

        boxsize = 0.2
        sb_radius = 1
        comp_space = sb_cubical_vr_comparison_space_via_cover(utb_circle_data, boxsize, sb_radius)
        traj = RefinedEquidistantTrajectory(utb_circle_data)

        tsn = TrajectorySpaceNew(traj, comp_space, DynamicDistance(2, sb_radius))

        # check for dimension 0, 1, 2 at the start
        for i = 1:200
            sig = cycling_signature(Val(:DistanceMatrix), tsn, (1,i), boxsize)
            if i <= 78
                @test dimension(sig) == 0
            elseif i <= 154
                @test dimension(sig) == 1
            else
                @test dimension(sig) == 2
            end
        end
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
        cycling_matrix = [1 0; 0 1]
        cs = CyclingSignature(cycling_matrix, birth_vector)
        @test cs.cycling_matrix == cycling_matrix
        @test cs.birth_vector == birth_vector

        # Test with incompatible sizes
        @test_throws ArgumentError CyclingSignature(cycling_matrix, [1.0, 2.0])
        @test dimension(cs) == 2
    end

    @testset "TrajectorySpaceNew" begin
        # generate traj
        traj = Matrix(collect(0:0.1:1)'[:, :])
        rt = RefinedEquidistantTrajectory(traj, [1, 6, 12])

        # generate comparison space
        comp_space = cubical_vr_comparison_space_via_cover(traj, 0.2)

        # generate distance
        d = euclidean

        tsn = TrajectorySpaceNew(rt, comp_space, d)

        @test trajectory(tsn) == rt
        @test comparison_space(tsn) == comp_space
        @test metric(tsn) == d
    end
end
