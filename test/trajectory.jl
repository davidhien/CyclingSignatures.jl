using CyclingSignatures
using Test

@testset "RefinedEquidistantTrajectory" begin
    traj = Matrix(collect(0:0.1:1)'[:, :])
    rt = RefinedEquidistantTrajectory(traj, [1, 6, 12])
    rt2 = RefinedEquidistantTrajectory(traj)

    @test rt2.t_vec == RefinedEquidistantTrajectory(traj, collect(1:12)).t_vec

    @test_throws ArgumentError RefinedEquidistantTrajectory(traj, [1, 5, 11])
    @test_throws ArgumentError RefinedEquidistantTrajectory(traj, [5, 1, 12])

    @test time_domain(rt) == 1:2

    @test evaluate_interval(rt, 1, 1) == traj[:, 1:5]
    @test evaluate_interval(rt, 2, 2) == traj[:, 6:11]
    @test evaluate_interval(rt, 1, 2) == traj[:, 1:11]

    @test isapprox(max_consecutive_distance(rt, (x, y) -> abs(x[1] - y[1])), 0.1)
    @test curve_hypothesis_violations(rt, (x, y) -> abs(x[1] - y[1]), 0.15) == 0
    @test curve_hypothesis_violations(rt, (x, y) -> abs(x[1] - y[1]), 0.05) == 10
end
