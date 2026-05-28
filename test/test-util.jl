function circle_time_series(subdivision, turns)
    a = LinRange(0, 2*π*turns, subdivision*turns + 1)
    return [cos.(a)'; sin.(a)']
end

function double_circle_time_series(subdivision, bits)
    a = zeros(2,0)
    template = circle_time_series(subdivision, 1)[:,1:end-1]
    l = template .+ [-1,0]
    r =  [-1,1] .*template .+ [1,0]
    for b in bits
        if b == 0
            a = hcat(a, l)
        elseif b == 1
            a = hcat(a, r)
        end
    end
    #note: l[:,1] == r[:,1]
    return [a l[:,1]]
end

function normalized_circle_time_series(subdivision, turns)
    return mapslices(v -> normalize(v, Inf), circle_time_series(subdivision, turns), dims=2)
end

function trajectory_space_for_points(points, boxsize, metric; flt_max_heuristic=nothing)
    comp_space = cubical_vr_comparison_space_via_cover(points, boxsize)
    traj = RefinedEquidistantTrajectory(points)
    if flt_max_heuristic === nothing
        return TrajectorySpace(traj, comp_space, metric)
    end
    return TrajectorySpace(traj, comp_space, metric, flt_max_heuristic)
end

function circle_trajectory_space(subdivision, turns; boxsize=0.2, metric=euclidean, flt_max_heuristic=nothing)
    points = normalized_circle_time_series(subdivision, turns)
    return trajectory_space_for_points(points, boxsize, metric; flt_max_heuristic), points
end

function figure8_trajectory_space(subdivision, bits; boxsize=0.2, metric=euclidean, flt_max_heuristic=nothing)
    points = double_circle_time_series(subdivision, bits)
    return trajectory_space_for_points(points, boxsize, metric; flt_max_heuristic), points
end

function utb_circle_trajectory_space(subdivision, turns; boxsize=0.2, sb_radius=1)
    circle_points = normalized_circle_time_series(subdivision, turns)
    circle_derivatives = circle_points[:, 2:end] - circle_points[:, 1:end-1]
    points = [circle_points[:, 1:end-1]; circle_derivatives]
    comp_space = sb_cubical_vr_comparison_space_via_cover(points, boxsize, sb_radius)
    traj = RefinedEquidistantTrajectory(points)
    return TrajectorySpace(traj, comp_space, DynamicDistance(2, sb_radius)), points
end

function utb_figure8_trajectory_space(subdivision, bits; boxsize=0.2, sb_radius=1)
    circle_points = double_circle_time_series(subdivision, bits)
    circle_derivatives = circle_points[:, 2:end] - circle_points[:, 1:end-1]
    points = [circle_points[:, 1:end-1]; circle_derivatives]
    comp_space = sb_cubical_vr_comparison_space_via_cover(points, boxsize, sb_radius)
    traj = RefinedEquidistantTrajectory(points)
    return TrajectorySpace(traj, comp_space, DynamicDistance(2, sb_radius)), points
end

function test_dimension_profile(alg, trajectory_space, r_max, ranges)
    for (range, expected_dimension) in ranges
        for i in range
            sig = cycling_signature(alg, trajectory_space, (1, i), r_max)
            @test dimension(sig) == expected_dimension
        end
    end
end
