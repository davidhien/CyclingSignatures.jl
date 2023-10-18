include("../../src/CyclingSignatures-include-file.jl")
using JLD2
using Colors
using ColorSchemes

dir_name = "examples/PaperExperiments/"
dwResults = SubsegmentResultReduced.(load(dir_name*"dw-results.jld2", "dwResults"));

dw_r = 0.18
rankFig = plotRanksWithLegend(dwResults, dw_r)
sig1Fig = plotSignaturesWithLegend(dwResults, 1, dw_r; cutoff=50)
sig2Fig = plotSignaturesWithLegend(dwResults, 2, dw_r; cutoff=100, legend_kwargs=(;nbanks=10))
inclFig = plotSubspaceInclusion(dwResults, 1, 2, dw_r; cutoff1=200, cutoff2=200)

figCombined = Figure(res=(1000,600),fontsize=18)
g11 = figCombined[1,1] = GridLayout()
g12 = figCombined[1,2] = GridLayout()
g21 = figCombined[2,1] = GridLayout()
g22 = figCombined[2,2]
plotRanksWithLegend(dwResults, dw_r; gl=g11)
plotSignaturesWithLegend(dwResults, 1, dw_r; gl=g12, cutoff=200)
plotSignaturesWithLegend(dwResults, 2, dw_r; gl=g21, cutoff=500, legend_kwargs=(;nbanks=5))
plotSubspaceInclusion(dwResults, 1, 2, dw_r; gp=g22, cutoff1=200, cutoff2=500, axis_kwargs=(;limits=(.5,3.5,.8,2.2)))
