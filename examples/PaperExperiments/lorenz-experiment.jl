include("../../src/CyclingSignatures-include-file.jl")

using DynamicalSystems
using JLD2

function lorenzTrajForCycling(u0=[0.0, 10.0,0.0])
    ds = PredefinedDynamicalSystems.lorenz(u0)
    dt = 0.1
    sol,_ = trajectory(ds, 10000, Δt=0.01)

    Y_init = Matrix(sol[:,:])'[:,100:end]
    resample_boxsize = 1.0
    resample_sb_radius = 3.0
    lorenz_f(v) = Vector(ds.integ.f(v, ds.integ.p, 0))
    
    Y_res, t_vec = resampleToConsistent(ds, Y_init, resample_boxsize, dt, sb_radius = resample_sb_radius, sb_fct=lorenz_f)
    Z_res = mapslices(v->normalize(lorenz_f(v),2 ), Y_res, dims=[1])
    return Y_res, Z_res, t_vec
end

Y_lorenz, Z_lorenz, t_lorenz = lorenzTrajForCycling()

boxsize_lorenz=8
sb_radius_lorenz = 3

ts_lorenz = trajectoryToTrajectorySpaceSB(Y_lorenz, Z_lorenz, boxsize_lorenz, sb_radius_lorenz, t_vec = t_lorenz; filter_missing=true, shortest_path_fix=true)

n = 1000
lorenzExperimentParams = map(i->SubsegmentSampleParameter(i,n), 10:10:1000)

lorenzExperiments = map(lorenzExperimentParams) do p
    RandomSubsegmentExperiment(ts_lorenz, p, boxsize_lorenz)
end

times, lorenzResults = runExperimentsWithTimer(lorenzExperiments)

dir_name = "examples/PaperExperiments/"
jldsave(dir_name*"lorenz-results.jld2", lorenzResults=lorenzResults, times=times)