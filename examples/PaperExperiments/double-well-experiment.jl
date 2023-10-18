include("../../src/CyclingSignatures-include-file.jl")

using JLD2

ts_dir_name = "examples/PaperExperiments/"
    
data_doublewell = load(ts_dir_name*"dw-timeseries.jld2", "data_doublewell")

boxsize_dw, sb_radius_dw = .2, 3
Y_dw = data_doublewell[:,1:end-1]
Z_dw = mapslices(normalize, data_doublewell[:,2:end]-data_doublewell[:,1:end-1], dims=1)

ts_dw = trajectoryToTrajectorySpaceSB(Y_dw,Z_dw, boxsize_dw, sb_radius_dw, filter_missing=true,shortest_path_fix=true)

n = 1000
dwExperimentParams = map(i->SubsegmentSampleParameter(i,n), 10:10:1000)

dwExperiments = map(dwExperimentParams) do p
    RandomSubsegmentExperiment(ts_dw, p, boxsize_dw)
end

times, dwResults = runExperimentsWithTimer(dwExperiments, Val(:DistanceMatrix));

dir_name = "examples/PaperExperiments/"
jldsave(dir_name*"dw-results.jld2", dwResults=dwResults, times=times)