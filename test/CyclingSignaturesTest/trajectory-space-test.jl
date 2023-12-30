include("../../src/CyclingSignatures-include-file.jl")
include("circle-test-util.jl")
using Test

function testRemoveDependentIntervals()
    M = [
        1 1 1 1 1;
        0 0 1 0 1;
        0 0 0 1 0 ]
    int_list = map(enumerate(eachcol(M))) do t
        i,col = t
        return PersistenceInterval(i,Inf,inclusion_representative=col)
    end
    int_list_mod = removeDependentIntervals(int_list)
    birth_list = map(x->x.birth, int_list_mod)
    return birth_list == [1;3;4]
end

@testset "Test removeDependentIntervals" begin
    @test testRemoveDependentIntervals()
end

@testset "Test ResampledTrajectory" begin
    A = reshape(collect(1:10),(2,5))
    @test_throws ArgumentError ResampledTrajectory(A, [1;3;7])
    @test_throws ArgumentError ResampledTrajectory(A, [1;3;2;6])
    rt = ResampledTrajectory(A, [1;3;6])
    @test get(rt,1:3) == A
end

@testset "Test BoxSpace" begin
    # TODO
end

@testset "Test SBBoxSpace" begin
    # TODO
end

function getCircleBoxSpace()
    
end

###
### Test evaluateCycling
###

function longCircleTimeSeries()
    a = 0:.05:20*pi
    return [cos.(a)'; sin.(a)']
end


@testset "evaluateCycling" begin
    @testset "circle time series" begin
        Y = circleTimeSeries(10)
        boxsize = 0.1
        ts = trajectoryToTrajectorySpace(Y, boxsize)
        t_ranges = map(v->1:v,50:50:250)
        dm_pds = map(t_ranges) do range
            evaluateCycling(Val(:DistanceMatrix), ts, range, .1)
        end
        expected_ranks = [0;0;1;1;1]
        dm_ranks = map(length, dm_pds)
        @test expected_ranks == dm_ranks
    end
end

