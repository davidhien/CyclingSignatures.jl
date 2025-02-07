include("../../src/CyclingSignatures-include-file.jl")
using JLD2
using Colors
using ColorSchemes

dir_name = "examples/PaperExperiments/"
lorenzResults = SubsegmentResultReduced.(load(dir_name*"lorenz-results.jld2", "lorenzResults"));

lorenz_r = 6
rankFig = plotRanksWithLegend(lorenzResults, lorenz_r)
sig1Fig = plotSignaturesWithLegend(lorenzResults, 1, lorenz_r; cutoff=50)
sig2Fig = plotSignaturesWithLegend(lorenzResults, 2, lorenz_r; cutoff=10, legend_kwargs=(;nbanks=10))
inclFig = plotSubspaceInclusion(lorenzResults, 1, 2, lorenz_r)

figCombined = Figure(res=(1000,600),fontsize=18)
g11 = figCombined[1,1] = GridLayout()
g12 = figCombined[1,2] = GridLayout()
g21 = figCombined[2,1] = GridLayout()
g22 = figCombined[2,2]
plotRanksWithLegend(lorenzResults, lorenz_r; gl=g11)
plotSignaturesWithLegend(lorenzResults, 1, lorenz_r; gl=g12)
plotSignaturesWithLegend(lorenzResults, 2, lorenz_r; gl=g21, legend_kwargs=(;nbanks=5))
plotSubspaceInclusion(lorenzResults, 1, 2, lorenz_r; gp=g22, axis_kwargs=(;limits=(.5,3.5,.8,2.2)))


# plot a segment
trajLorenz = getTrajectory(getTrajectorySpace(lorenzResults[10].experiment))
sig1_l,_ = subspaceFrequencyMatrix(lorenzResults, 1, lorenz_r)
sig_1_range = getSignatureRanges(lorenzResults[10], sig1_l[1], lorenz_r)

background_pts  = get(trajLorenz, 1:10000)[1:3,:] 
segment_pts = get(trajLorenz, sig_1_range[1])[1:3,:]
scatter(background_pts, color=:black, markersize=1)
scatter!(segment_pts, color=:red)