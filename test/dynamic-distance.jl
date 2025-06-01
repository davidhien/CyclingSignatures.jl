using Test
using Distances
using CyclingSignatures: DynamicDistance, get_dim, get_c

@testset "dynamic distance" begin
    x0 = [0,0]
    x1 = [4,0]
    x2 = [0,3]
    x3 = [0,1]

    dist = DynamicDistance(1,2)
    @test dist(x1,x0) == 4
    @test dist(x2,x0) == 6
    @test dist(x2,x1) == 6
    @test dist(x3,x1) == 4

    @test result_type(dist, Int, Float64) == Float64
end
