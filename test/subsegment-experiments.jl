using CyclingSignatures
using StepFunctions
using Test
using LinearAlgebra
using Random
using CyclingSignatures: sample_segment_starts


include("test-util.jl")

@testset "test helper methods" begin
    @testset "sample_segment_starts" begin
    rng1 = MersenneTwister(1234)
    rng2 = MersenneTwister(1234)

    # test
    s1 = sample_segment_starts(100, 10, 10, rng1)
    @test length(s1) == 10
    @test all(s1 .>= 1) && all(s1 .<= 91)

    # test that same exp gives same result
    s2 = sample_segment_starts(100, 10, 10, rng2)
    @test s1 == s2
    end

    @testset "cycspace_inclusion_matrix" begin
        F = FF{47}

        # Zero subspace (0-dimensional): 3×0 matrix
        U0 = F.(reshape(Int[], 3, 0))

        # Vectors in R^3 over FF{47}:
        e1 = F.([1, 0, 0])
        e2 = F.([0, 1, 0])
        e12 = F.([1, 1, 0])

        # Subspaces with differing column counts
        U1 = hcat(e1)                # 3×1 : span{e1}
        U2 = hcat(e1, e2)            # 3×2 : span{e1, e2}
        U3 = hcat(e12, e1)           # 3×2 : span{e1 + e2, e1} = U2
        U4 = hcat(e12)               # 3×1 : span{e1 + e2}

        # V and W purposely built with mixed shapes
        V = [U0, U1, U4]
        W = [U0, U1, U2, U3, U4]

        # Inclusion analysis:
        #
        # U0 ⊆ everything
        #
        # U1 ⊆ U1, U2, U3 but not U4
        #
        # U4 = span{e1 + e2}
        #   ⊆ U2 and U3 (since they span {e1,e2})
        #   not included in U1
        #   included in U4
        #
        expected = F.([
            1 1 1 1 1;   # U0 included everywhere
            0 1 1 1 0;   # U1
            0 0 1 1 1    # U4
        ])

        incl_mat = cycspace_inclusion_matrix(V, W)
        @test incl_mat == expected
    end
end

@testset "experiment pipeline" begin
    # ten turns around the circle
    circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(20, 10), dims=2)
    boxsize = 0.5
    # test cubical_vr_comparison_space_via_cover
    traj_space = trajectory_space_from_trajectory(circle_data, boxsize)

    # test both constructors
    segment_lengths = [10, 20, 30]
    n_runs = 5
    seed = 1234

    # constructor without seed
    test_exp = RandomSubsegmentExperiment(traj_space, segment_lengths, n_runs)
    @test get_trajectory_space(test_exp) === traj_space
    @test get_segment_lengths(test_exp) === segment_lengths
    @test get_n_experiments(test_exp) == n_runs

    # constructor with seed
    exp = RandomSubsegmentExperiment(traj_space, segment_lengths, n_runs, seed)
    @test get_trajectory_space(exp) === traj_space
    @test get_segment_lengths(exp) === segment_lengths
    @test get_n_experiments(exp) == n_runs
    @test exp.seed == seed

    # test run experiment
    # need to test: resample_segment_start, threshold and that seed is applied correctly
    result = run_experiment(exp::RandomSubsegmentExperiment;
        field=FF{2},
        alg=Val(:DistanceMatrix),
        threshold=nothing,
        resample_segment_start=true,
        progress=false,
        parallel_inner=false)

    @test result.trajectory_space === traj_space
    @test result.segment_lengths == segment_lengths
    @test result.n_runs == n_runs
    @test length(result.segment_starts) == length(segment_lengths)
    @test all(length.(result.segment_starts) .== n_runs)
    @test result.flt_threshold === boxsize # this is the default heuristic

    # test that running again with same seed gives same result
    result2 = run_experiment(exp::RandomSubsegmentExperiment;
        field=FF{2},
        alg=Val(:DistanceMatrix),
        threshold=nothing,
        resample_segment_start=true,
        progress=false,
        parallel_inner=false)

    @test result.segment_starts == result2.segment_starts
end


@testset "test analysis methods" begin
    # ten turns around the circle
    circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(20, 10), dims=2)
    boxsize = 0.5
    # test cubical_vr_comparison_space_via_cover
    traj_space = trajectory_space_from_trajectory(circle_data, boxsize)

    # test both constructors
    segment_lengths = [10, 20, 30] # TODO: include length 19 here
    n_runs = 5
    seed = 1234

    # constructor without seed
    test_exp = RandomSubsegmentExperiment(traj_space, segment_lengths, n_runs)
    result = run_experiment(test_exp; progress=false, parallel_inner=false)

    r_dist_0 = rank_distribution(result, 0)
    @test r_dist_0[1] == StepFunction([0], 0, [5])
    @test r_dist_0[2](0) == 5 && r_dist_0[2](0.7) == 0
    @test r_dist_0[3](0) == 5 && r_dist_0[2](0.4) == 0

    # test rank distribution
    r_dist_1 = rank_distribution(result, 1)
    @test r_dist_1[1] == StepFunction(Int[], 0, [])
    @test r_dist_1[2](0) == 0 && r_dist_1[2](0.7) == 5
    @test r_dist_1[3](0) == 0 && r_dist_1[2](0.4) == 5

    # test cycspace intervals
    cs_intervals = cycspace_intervals(result, [1][:,:])
    @test isempty(cs_intervals[1])  # no rank 1 cycles for length 10 segments
    @test length(cs_intervals[2]) == 5  # all length 20
    @test all(==(Inf), last.(cs_intervals[2]))
    @test length(cs_intervals[3]) == 5  # all length 30
    @test all(==(Inf), last.(cs_intervals[3]))

    # cycspace_distribution
    cs_dist = cycspace_distribution(result, [1][:,:])
    @test cs_dist[1] == StepFunction(Int[], 0, [])
    @test cs_dist[2](0) == 0 && cs_dist[2](0.7) == 5
    @test cs_dist[3](0) == 0 && cs_dist[3](0.4) == 5

    # cycspace_radius_distribution and cycspace_timespan_distribution

    # cycspace_segments
    segs = cycspace_segments(result, [1][:,:])
    @test isempty(segs[1])  # no rank 1 cycles for length 10 segments
    @test length(segs[2]) == 5
    @test all(length.(segs[2]) .== 20)
    @test length(segs[3]) == 5  # all length 20
    @test all(length.(segs[3]) .== 30)

    # cycspace_segments_at_r
    segs_r = cycspace_segments_at_r(result, [1][:,:], 0.5)
    @test isempty(segs_r[1])  # no rank 1 cycles for length 10 segments
    @test length(segs[2]) == 5
    @test all(length.(segs[2]) .== 20)
end
