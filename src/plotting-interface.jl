const _PLOTTING_EXT_ERROR =
    "Load Plots/StatsPlots or Makie (with a backend) before calling plotting functions."

plot_rank_distribution(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_rank_distribution!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_rank_heatmap(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_rank_heatmap!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_all_rank_heatmaps(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_all_rank_heatmaps!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_rank_distribution_at_r(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_rank_distribution_at_r!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_subspace_frequency_at_r(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_subspace_frequency_at_r!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_cycspace_inclusion(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_inclusion!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_cycspace_radius_frequency(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_radius_frequency!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_cycspace_length_frequency(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_length_frequency!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
