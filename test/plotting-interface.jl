using Test
using CyclingSignatures

struct PlotTestTrajectorySpace
    flt_max_heuristic::Float64
end

CyclingSignatures.betti_1(::PlotTestTrajectorySpace) = 1

function synthetic_plot_result()
    F = FF{2}
    e1 = reshape(F.([1, 0]), 2, 1)
    e2 = reshape(F.([0, 1]), 2, 1)
    zero_sig = CyclingSignature(zeros(F, 2, 0), Float64[])

    signatures = [
        CyclingSignature[CyclingSignature(e1, [0.2]), zero_sig],
        CyclingSignature[CyclingSignature(e1, [0.1]), CyclingSignature(e2, [0.3])],
        CyclingSignature[CyclingSignature(e1, [0.4]), CyclingSignature(e1, [0.6])],
    ]

    return RandomSubsegmentResult(
        PlotTestTrajectorySpace(1.0),
        [10, 20, 20],
        2,
        [[1, 2], [1, 2], [1, 2]],
        1.0,
        signatures,
    )
end

function synthetic_cycling_spaces()
    F = FF{2}
    e1 = reshape(F.([1, 0]), 2, 1)
    e2 = reshape(F.([0, 1]), 2, 1)
    return e1, [zeros(F, 2, 0), e1, e2], [zeros(F, 2, 0), hcat(e1, e2)]
end

@testset "plotting interface fallback" begin
    @test :plot_cycspace_distribution in names(CyclingSignatures)
    @test :plot_cycspace_distribution! in names(CyclingSignatures)
    @test :plot_cycspace_level_contours in names(CyclingSignatures)
    @test :plot_cycspace_level_contours! in names(CyclingSignatures)

    @test_throws ErrorException CyclingSignatures.plot_rank_heatmap(synthetic_plot_result(), 1)
end

@testset "Plots extension" begin
    using Plots
    using StatsPlots

    @test Base.get_extension(CyclingSignatures, :PlotsPltExt) !== nothing

    result = synthetic_plot_result()
    cycling_space, V0, V1 = synthetic_cycling_spaces()

    @test plot_rank_distribution(result, 1) isa Plots.Plot
    @test plot_rank_heatmap(result, 1) isa Plots.Plot
    @test plot_all_rank_heatmaps(result) isa Plots.Plot
    @test plot_rank_distribution_at_r(result, 0.5) isa Plots.Plot

    sig_at_r, plt_at_r = plot_subspace_frequency_at_r(result, 1, 0.5; n_subspaces = 2)
    @test length(sig_at_r) == 2
    @test plt_at_r isa Plots.Plot

    @test plot_cycspace_inclusion(V0, V1) isa Plots.Plot

    sig_radius, plt_radius = plot_cycspace_radius_frequency(result, 1; n_subspaces = 2)
    @test length(sig_radius) == 2
    @test plt_radius isa Plots.Plot

    @test plot_cycspace_distribution(result, cycling_space; radius_bins = 3) isa Plots.Plot

    sig_length, plt_length = plot_cycspace_length_frequency(result, 1; n_subspaces = 2)
    @test length(sig_length) == 2
    @test plt_length isa Plots.Plot

    @test plot_cycspace_level_contours(result, cycling_space, [1, 2]) isa Plots.Plot
    @test plot_cycspace_level_contours(result, cycling_space, 1; mode = :bands) isa Plots.Plot
    @test_throws ArgumentError plot_cycspace_level_contours(result, cycling_space, 1; mode = :bad)
    @test_throws ArgumentError plot_cycspace_level_contours(result, cycling_space, 1; x_mode = :scaled)
end

@testset "Makie extension" begin
    using Makie

    @test Base.get_extension(CyclingSignatures, :MakiePltExt) !== nothing

    result = synthetic_plot_result()
    cycling_space, V0, V1 = synthetic_cycling_spaces()

    @test plot_rank_distribution(result, 1) isa Makie.Figure
    @test plot_rank_heatmap(result, 1) isa Makie.Figure
    @test plot_all_rank_heatmaps(result) isa Makie.Figure
    @test plot_rank_distribution_at_r(result, 0.5) isa Makie.Figure

    sig_at_r, fig_at_r = plot_subspace_frequency_at_r(result, 1, 0.5; n_subspaces = 2)
    @test length(sig_at_r) == 2
    @test fig_at_r isa Makie.Figure

    @test plot_cycspace_inclusion(V0, V1) isa Makie.Figure

    sig_radius, fig_radius = plot_cycspace_radius_frequency(result, 1; n_subspaces = 2)
    @test length(sig_radius) == 2
    @test fig_radius isa Makie.Figure

    @test plot_cycspace_distribution(result, cycling_space; radius_bins = 3) isa Makie.Figure

    sig_length, fig_length = plot_cycspace_length_frequency(result, 1; n_subspaces = 2)
    @test length(sig_length) == 2
    @test fig_length isa Makie.Figure

    @test plot_cycspace_level_contours(result, cycling_space, [1, 2]) isa Makie.Figure
    @test plot_cycspace_level_contours(result, cycling_space, 1; mode = :bands) isa Makie.Figure
    @test_throws ArgumentError plot_cycspace_level_contours(result, cycling_space, 1; mode = :bad)
    @test_throws ArgumentError plot_cycspace_level_contours(result, cycling_space, 1; x_mode = :scaled)
end
