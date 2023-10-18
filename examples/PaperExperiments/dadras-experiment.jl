include("../DadrasUtil.jl")
include("../../src/CyclingSignatures-include-file.jl")

using DynamicalSystems
using DifferentialEquations
using ForwardDiff
using JLD2

# data generation
function dadrasTimeSeries()
    p = [8;40;14.9]
    ic = [10.;1.;10.;1.]
    dt = .01
    ds = ContinuousDynamicalSystem(eom_dadras!, ic, p)
    sol,_ = trajectory(ds, 100000, Δt=dt)
    dat = Matrix(sol)'[:,:]
    m = 1000000
    resample_boxsize = .8
    resample_sb_r = .2
    Y_res, t_vec = resampleToDistance(ds, dt, dat[:,1:m], resample_boxsize; pp=rescale, sb_r=resample_sb_r, sb_fct=x->rescaled_eom(rescale(x), p, 0), max_depth=512, verbose=false)
    Y_res_rescaled = mapslices(rescale, Y_res, dims=[1])
    Z_res_rescaled = mapslices(x -> normalize(rescaled_eom(x, p, 0),2), Y_res_rescaled, dims=[1])

    return Y_res_rescaled, Z_res_rescaled, t_vec
end

Y_dadras, Z_dadras, t_vec = dadrasTimeSeries()
boxsize_dadras = 4
sb_radius_dadras = 3

# preprocessing
ts_dadras = trajectoryToTrajectorySpaceSB(Y_dadras, Z_dadras, boxsize_dadras, sb_radius_dadras, t_vec = t_vec; filter_missing=true, shortest_path_fix=true)

# experiment parameters
n = 1000
dadrasExperimentParams = map(i->SubsegmentSampleParameter(i,n), 10:10:1000)

dadrasExperiments = map(dadrasExperimentParams) do p
    RandomSubsegmentExperiment(ts_dadras, p, boxsize_dadras)
end

# run experiments
times, dadrasResults = runExperimentsWithTimer(dadrasExperiments, Val(:DistanceMatrix))

# save data
dir_name = "examples/PaperExperiments/"
jldsave(dir_name*"dadras-results.jld2", dadrasResults=dadrasResults, times=times)