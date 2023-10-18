include("../../src/CyclingSignatures/trajectory.jl")
using Test

@testset "test countDynamicInconsistencies" begin
    # box iterator input
    @test countDynamicInconsistencies(eachrow([1;4;2;3;4;10]),2) == 2

    # vector of boxes input
    @test countDynamicInconsistencies([[1];[4];[2];[3];[4];[10]],2) == 2

    # edge cases
    @test countDynamicInconsistencies([],2) == 0
    @test countDynamicInconsistencies([1],2) == 0
end


