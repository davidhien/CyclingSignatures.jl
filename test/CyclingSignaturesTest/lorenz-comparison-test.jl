include("../../src/CyclingSignatures-include-file-Ripserer.jl")
using DynamicalSystems
using Test
using JLD2

function lorenzTrajForCycling()
    ds = PredefinedDynamicalSystems.lorenz()
    dt = 0.1
    sol,_ = trajectory(ds, 10000, Δt=0.01)

    Y_init = Matrix(sol[:,:])'[:,100:end]
    resample_boxsize = 1.0
    resample_sb_radius = 3.0

    #resample_r = 1
    lorenz_f(v) = Vector(ds.integ.f(v, ds.integ.p, 0))
    
    Y_res, t_vec = resampleToConsistent(ds, Y_init, resample_boxsize, dt, sb_radius = resample_sb_radius, sb_fct=lorenz_f)
    #Y_res, t_vec = resampleToDistance(ds, dt, Y_init, resample_r; sb_r = .25, sb_fct=lorenz_f)
    Z_res = mapslices(v->normalize(lorenz_f(v),2 ), Y_res, dims=[1])
    return Y_res, Z_res, t_vec
end

function loadLorenzTestData()
    #Y_lorenz, Z_lorenz, t_lorenz = lorenzTrajForCycling()

    #jldsave("test/CyclingSignaturesTest/lorenz-data.jld2", Y_lorenz=Y_lorenz, Z_lorenz=Z_lorenz, t_lorenz=t_lorenz)
    Y_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Y_lorenz")
    Z_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Z_lorenz")
    t_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "t_lorenz")
    return Y_lorenz, Z_lorenz, t_lorenz
end

@testset "Comparison Test DistanceMatrix and Ripserer" begin
    Y_lorenz, Z_lorenz, t_lorenz = loadLorenzTestData()

    boxsize_lorenz = 8
    sb_radius_lorenz=3
    ts_lorenz = trajectoryToTrajectorySpaceSB(Y_lorenz, Z_lorenz, boxsize_lorenz, sb_radius_lorenz, t_vec = t_lorenz; filter_missing=true, shortest_path_fix=false)
    n = 20
    lorenzExperimentParams = map(i->SubsegmentSampleParameter(i,n), 10:10:500)

    lCompExp = map(lorenzExperimentParams) do p
        RandomSubsegmentExperiment(ts_lorenz, p, boxsize_lorenz)
    end

    lCompRes = runExperiments(lCompExp, Val(:DistanceMatrix);verbose=false);
    lCompResR = runExperiments(lCompExp, Val(:DistanceMatrix);verbose=false,field=Mod{2});

    lCompDiag = getDiagrams.(lCompRes)
    lCompDiagR = getDiagrams.(lCompResR)
    comp_vec = map(zip(lCompDiag,lCompDiagR)) do t
        d_v1, d_v2 = t
        any(length.(d_v1) .!= length.(d_v2))
    end
    @test sum(comp_vec) == 0

    # for r > ch_val the curve hypothesis guarantees the same subsapce normal form (and birth values)
    # TODO: test birth values for r > ch_val
    @show ch_val = curveHypothesis(getTrajectory(ts_lorenz), getMetric(ts_lorenz))
    eval_r = 2/3ch_val + 1/3*boxsize_lorenz

    ssNFs = map(vcat(lCompDiag...)) do diag
        diag_f = filter(int -> birth(int) <eval_r, collect(diag))
        cycMat = diagToCyclingMat(diag_f)
        return Int.(subspaceNormalForm(cycMat))
    end

    ssNFsR = map(vcat(lCompDiagR...)) do diag
        diag_f = filter(int -> birth(int) <eval_r, collect(diag))
        cycMat = diagToCyclingMat(diag_f)
        return Int.(subspaceNormalForm(cycMat))
    end

    @test all(ssNFs .== ssNFsR)
end
