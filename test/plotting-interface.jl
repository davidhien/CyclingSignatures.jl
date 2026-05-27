using Test
using CyclingSignatures

@testset "plot_cycspace_distribution API" begin
    exported_names = names(CyclingSignatures)
    @test :plot_cycspace_distribution in exported_names
    @test :plot_cycspace_distribution! in exported_names
    @test :plot_cycspace_level_contours in exported_names
    @test :plot_cycspace_level_contours! in exported_names

    makie_src = read(joinpath(dirname(@__DIR__), "ext", "MakiePltExt.jl"), String)
    plots_src = read(joinpath(dirname(@__DIR__), "ext", "PlotsPltExt.jl"), String)

    @test occursin("function CyclingSignatures.plot_cycspace_distribution!", makie_src)
    @test occursin("function CyclingSignatures.plot_cycspace_distribution(", makie_src)
    @test occursin("CyclingSignatures._cycspace_distribution_heatmap_data(", makie_src)

    @test occursin("function CyclingSignatures.plot_cycspace_distribution!", plots_src)
    @test occursin("function CyclingSignatures.plot_cycspace_distribution(", plots_src)
    @test occursin("CyclingSignatures._cycspace_distribution_heatmap_data(", plots_src)

    @test occursin("function CyclingSignatures.plot_cycspace_level_contours!", makie_src)
    @test occursin("function CyclingSignatures.plot_cycspace_level_contours(", makie_src)
    @test occursin("CyclingSignatures.cycspace_level_intervals(", makie_src)

    @test occursin("function CyclingSignatures.plot_cycspace_level_contours!", plots_src)
    @test occursin("function CyclingSignatures.plot_cycspace_level_contours(", plots_src)
    @test occursin("CyclingSignatures.cycspace_level_intervals(", plots_src)
end
