include("../../src/CyclingSignatures-include-file.jl")
include("circle-test-util.jl")
using Test
using JLD2

@testset "small unit tests" begin
    field = FF{3}
    edges, coeffs = curveGeneratorData(1,5;F=field)

    @test edges == [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
    @test coeffs == field.([1;1;1;1;-1])
end

@testset "circle datasets" begin
    @testset "simple circle 4 points" begin
        trajPoints = [0 1 2 1; 0 1 0 -1]
        metric = chebyshev 
        fltThreshold = 1.1

        cc, d_mat,cc_labels = twoPassDMAnnotationFirstPass(trajPoints, metric, fltThreshold)

        @test size(d_mat) == size(cc_labels)
        @test num_groups(cc) == 2

        smallestNodesDict = twoPassDMAnnotationSecondPass(cc,d_mat,cc_labels)

        @test length(collect(values(smallestNodesDict))) == 2

        smallestNodes = connectedComponentRepresentatives(smallestNodesDict, cc_labels)

        pd = persistenceDiagramFromNodes(smallestNodes, trajPoints, metric)
        @test length(pd) == 1
        bar = pd[1]
        @test birth(bar) == 1 && death(bar) == Inf
    end

    @testset "circle time series 10 turns" begin
        # TODO: make the following into an actual test
        Y = circleTimeSeries(10)
        trajPoints = Y[:,1:150]
        metric = chebyshev 
        threshold = .1
        pd = trajectoryBarcode(Val(:DistanceMatrix), trajPoints, metric, threshold)
        @test length(pd) == 1
    end
end

@testset "generator test" begin
    @testset "4 point circle test" begin
        trajPoints = [0 1 2 1; 0 1 0 -1]
        metric = chebyshev 
        fltThreshold = 1.1

        k = 23
        field = FF{k}        
        traj_barc = trajectoryBarcode(Val(:DistanceMatrix), trajPoints, metric, fltThreshold, field)
        @test traj_barc[1].simplex_list == [ (1, 2), (2, 3), (3, 4), (1, 4)]
        @test traj_barc[1].coeff_list == field.([1,1,1,-1])
    end
end

function loadLorenzTestData()
    #Y_lorenz, Z_lorenz, t_lorenz = lorenzTrajForCycling()

    #jldsave("test/CyclingSignaturesTest/lorenz-data.jld2", Y_lorenz=Y_lorenz, Z_lorenz=Z_lorenz, t_lorenz=t_lorenz)
    Y_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Y_lorenz")
    Z_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Z_lorenz")
    t_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "t_lorenz")
    return Y_lorenz, Z_lorenz, t_lorenz
end

@testset "Comparison DistanceMatrix and DistanceMatrixOld passes" begin
    Y_lorenz, Z_lorenz, t_lorenz = loadLorenzTestData()
    boxsize_lorenz = 8
    sb_radius_lorenz=1
    ts_lorenz = trajectoryToTrajectorySpaceSB(Y_lorenz, Z_lorenz, boxsize_lorenz, sb_radius_lorenz, t_vec = t_lorenz; filter_missing=true, shortest_path_fix=false)
    rt_lorenz = getTrajectory(ts_lorenz)
    d_lorenz = getMetric(ts_lorenz)
    rngs = [282670:282670+80,100:260,34130:34130+320,1200:3200]
    @testset "first pass" for rng in rngs
        trajPoints = get(rt_lorenz, rng)

        cc, _, cc_labels = twoPassDMAnnotationFirstPass(trajPoints, d_lorenz, boxsize_lorenz)
        cc2, cc_labels2 = twoPassDMAnnotationFirstPass2(trajPoints, d_lorenz, boxsize_lorenz)

        @test cc.parents == cc2.parents
        @test cc.ranks == cc2.ranks
        @test cc.ngroups == cc2.ngroups

        @test cc_labels == cc_labels2
    end
    @testset "second pass" for rng in rngs
        trajPoints = get(rt_lorenz, rng)

        cc, dmat, cc_labels = twoPassDMAnnotationFirstPass(trajPoints, d_lorenz, boxsize_lorenz)
        cc2, cc_labels2 = twoPassDMAnnotationFirstPass2(trajPoints, d_lorenz, boxsize_lorenz)

        dict_1 = twoPassDMAnnotationSecondPass(cc, dmat, cc_labels)
        dict_2 = twoPassDMAnnotationSecondPass2(cc, cc_labels, trajPoints, d_lorenz)
        @test dict_1 == dict_2
    end
end

@testset "Comparison Test DistanceMatrix and DistanceMatrixOld" begin
    Y_lorenz, Z_lorenz, t_lorenz = loadLorenzTestData()
    
    boxsize_lorenz = 8
    sb_radius_lorenz=1
    ts_lorenz = trajectoryToTrajectorySpaceSB(Y_lorenz, Z_lorenz, boxsize_lorenz, sb_radius_lorenz, t_vec = t_lorenz; filter_missing=true, shortest_path_fix=false)
    cyc1 = evaluateCycling(Val(:DistanceMatrix), ts_lorenz, 282670:282670+80, 8)[:]
    cyc2 = evaluateCycling(Val(:DistanceMatrixOld), ts_lorenz, 282670:282670+80, 8)[:]
    @test length(cyc1) == length(cyc2) == 1
    @test isapprox(cyc1[1].birth, cyc2[1].birth) 
end