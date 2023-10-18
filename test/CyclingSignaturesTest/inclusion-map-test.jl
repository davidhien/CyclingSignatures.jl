include("../../src/CyclingSignatures-include-file.jl")
include("../../examples/DadrasUtil.jl")
include("circle-test-util.jl")

using Test
using DynamicalSystems
using ForwardDiff
using Distances
using LinearAlgebra

function testEdgeBoxes1()
    p1 = [0;0]
    p2 = [1;0]
    
    eb = edgeBoxes(p1,p2,.1)
    res = [[i;0] for i = 0:10]
    return eb == res
end

@testset "edgeBoxes" begin
    @test testEdgeBoxes1()
    # TODO: add more
end

@testset "testNextMaxSwitchLambda" begin
    # line a over line b
    p1 = [2;1]
    p2 = [3;2]
    @test nextMaxSwitchLambda(p1,p2) == 1
    @test nextMaxSwitchLambda(-p1,-p2) == 1

    # line a intersects line b
    p1 = [2;1]
    p2 = [1;2]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    # line a over line b, b has irrelevant zero intersection 
    p1 = [2;1]
    p2 = [1;-.5]
    @test nextMaxSwitchLambda(p1,p2) == 1
    @test nextMaxSwitchLambda(-p1,-p2) == 1

    # line a positive, line b negative and takes max
    p1 = [2;0]
    p2 = [4;-6]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    # like previous but flipped sign at second component
    p1 = [2;0]
    p2 = [4;6]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    p1 = [2;0]
    p2 = [-1;1]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    p1 = [1;1]
    p2 = [3;-5]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    p1 = [+1;0]
    p2 = [-1;0]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.5)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.5)

    p1 = [0;1;0.8]
    p2 = [1;0;0.8]
    @test isapprox(nextMaxSwitchLambda(p1,p2),.2)
    @test isapprox(nextMaxSwitchLambda(-p1,-p2),.2)

    # TODO: more test cases with multiple entries and multiple phenomena
    # TODO: document phenomena, make sure they are tested, also in combination
end


@testset "allMaxSwitchLambda" begin
    p1 = [0;1;0.8]
    p2 = [1;0;0.8]
    bpl = allMaxSwitchLambda(p1,p2)
    
    @test length(bpl) == 2 && isapprox(bpl[1],.2) && isapprox(bpl[2],.8)
end

@testset "edgeBoxesSphereBundle" begin
    # TODO!
end

@testset "getInclusionEdgeBoxIndices" begin
    # TODO!
end

###
### Test for inclusion helpers
###

"""
Generates points for the circle dataset with angles 0:.05:2*pi and an inclusion helper for a cover with radius .1. 
"""
function generateCircleInclusionHelper()
    pts = circleTimeSeries()
    boxsize = .1
    X = quantize(pts, boxsize)
    X = unique(X,dims=2)
    dm_X = pairwise(Distances.chebyshev, X)
    cplx, h1_gen, _ = boxSpaceH1Info(dm_X)
    return pts,InclusionHelper(cplx, X, h1_gen, boxsize; filter_missing=true)
end

function testIncludeCycleCircleData1()
    # cycle wraps around circle once
    pts, inc_helper = generateCircleInclusionHelper()
    edgeList = map(i-> (i,i+1), 1:size(pts,2)-1)
    push!(edgeList, (1,size(pts,2)))
    coeffs = FF{2}.(ones(Int,length(edgeList)))
    inc_gen = includeCycle(inc_helper, pts, edgeList, coeffs)
    res = inc_helper.h1_gen_t * inc_gen
    return res[1,1] == 1
end

function testIncludeCycleCircleData2()
    # cycle is trivial
    pts, inc_helper = generateCircleInclusionHelper()
    edgeList = map(i-> (i,i+1), 1:10)
    push!(edgeList, (1,10))
    coeffs = FF{2}.(ones(Int,length(edgeList)))
    inc_gen = includeCycle(inc_helper, pts, edgeList, coeffs)
    res = inc_helper.h1_gen_t * inc_gen
    return res[1,1] == 0
end

@testset "InclusionHelper: Circle Dataset" begin
    @testset "Constructor Test"  begin
        pts, inc_helper = generateCircleInclusionHelper()
        @test issorted(inc_helper.cplxPointsSorted)
        # TODO 
    end
    @testset "Nontrivial Generator" begin
        @test testIncludeCycleCircleData1()
    end
    @testset "Trivial Generator" begin
        @test testIncludeCycleCircleData2()
    end
end

# Tests for SBInclusionHelper
function generateCircleSBInclusionHelper()
    pts = circleTimeSeries()
    pts_VF = mapslices(normalize, [-pts[1,:]';pts[2,:]'], dims=1)
    pts_all = [pts;pts_VF]
    boxsize = .1
    sb_radius = 3
    X = quantize(pts, boxsize)
    VF = quantizeSB(pts_VF, sb_radius)
    SB_pts = unique([X;VF], dims=2)
    dm_X = sphereBundleDistanceMatrix(SB_pts)
    cplx, h1_gen, _ = boxSpaceH1Info(dm_X)
    inc_helper = SBInclusionHelper(cplx, SB_pts, h1_gen, boxsize, sb_radius; filter_missing=true)
    return pts_all,inc_helper
end

function testGetInclusionEdgeBoxIndices()
    _,inc_helper = generateCircleSBInclusionHelper()
    incBoxInd = getInclusionEdgeBoxIndices([[10, 0, -3, 0], [10, 1, -3, 0]], inc_helper.cplxPointsSorted; incKwargs(inc_helper)...)
    # TODO: improve this check:
    return length(incBoxInd) == 2
end

function testIncludeCycleCircleDataSB1()
    # cycle wraps around circle once
    SB_pts, inc_helper = generateCircleSBInclusionHelper()
    edgeList = map(i-> (i,i+1), 1:size(SB_pts,2)-1)
    push!(edgeList, (1,size(SB_pts,2)))
    coeffs = FF{2}.(ones(Int,length(edgeList)))
    inc_gen = includeCycle(inc_helper, SB_pts, edgeList, coeffs)
    res = inc_helper.h1_gen_t * inc_gen
    return res[1,1] == 1
end

function testIncludeCycleCircleDataSB2()
    # cycle is trivial
    SB_pts, inc_helper = generateCircleSBInclusionHelper()
    edgeList = map(i-> (i,i+1), 1:10)
    push!(edgeList, (1,10))
    coeffs = FF{2}.(ones(Int,length(edgeList)))
    inc_gen = includeCycle(inc_helper, SB_pts, edgeList, coeffs)
    res = inc_helper.h1_gen_t * inc_gen
    return res[1,1] == 0
end

@testset "SBInclusionHelper: Circle Dataset" begin
    @testset "Constructor Test"  begin
        SB_pts, inc_helper = generateCircleSBInclusionHelper()
        @test issorted(inc_helper.cplxPointsSorted)
        # TODO 
    end
    @testset "getInclusionEdgeBoxIndices" begin
        @test testGetInclusionEdgeBoxIndices()
        # TODO
    end
    @testset "Nontrivial Generator" begin
        @test testIncludeCycleCircleDataSB1()
    end
    @testset "Trivial Generator" begin
        @test testIncludeCycleCircleDataSB2()
    end
end

# TODO:
# - InclusionHelpers
#   - test constructors, check that everything is sorted. Maybe ensure this in inner constructors
# - Test shortest path heuristic