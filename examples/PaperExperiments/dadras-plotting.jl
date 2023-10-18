include("../../src/CyclingSignatures-include-file.jl")
using JLD2
using Colors
using ColorSchemes

dir_name = "examples/PaperExperiments/"
dadrasResults = SubsegmentResultReduced.(load(dir_name*"dadras-results.jld2", "dadrasResults"));

rankFig = plotRanksWithLegend(dadrasResults, 2.5)
sig1Fig = plotSignaturesWithLegend(dadrasResults, 1, 2.5; cutoff=200)
sig2Fig = plotSignaturesWithLegend(dadrasResults, 2, 2.5;cutoff=100, legend_kwargs=(;nbanks=10))
inclFig = plotSubspaceInclusion(dadrasResults, 1, 2, 2.5; cutoff1=200, cutoff2=500)

figCombined = Figure(res=(1000,600),fontsize=18)
g11 = figCombined[1,1] = GridLayout()
g12 = figCombined[1,2] = GridLayout()
g21 = figCombined[2,1] = GridLayout()
g22 = figCombined[2,2]
plotRanksWithLegend(dadrasResults, 2.5; gl=g11)
plotSignaturesWithLegend(dadrasResults, 1, 2.5; gl=g12,cutoff=200)
plotSignaturesWithLegend(dadrasResults, 2, 2.5; gl=g21,cutoff=500, legend_kwargs=(;nbanks=5))
plotSubspaceInclusion(dadrasResults, 1, 2, 2.5; gp=g22,cutoff1=200, cutoff2=500, axis_kwargs=(;limits=(.5,6.5,.8,2.2)))
