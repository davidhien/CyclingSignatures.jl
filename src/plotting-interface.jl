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

"""
    plot_cycspace_distribution(results, cycling_space; r_max=nothing, radius_bins=nothing, kwargs...)

Plot a heatmap showing how often `cycling_space` appears across segments. The x-axis uses the
distinct segment time spans in `results`. The y-axis divides `[0, r_max]` into `radius_bins`
equal bins, defaulting to the number of distinct time spans. Each cell stores the number of
segments in that time-span column whose cycling-space distribution is active at the center of
that radius bin.
"""
plot_cycspace_distribution(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_distribution!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

plot_cycspace_length_frequency(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_length_frequency!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)

"""
    plot_cycspace_level_contours(results, cycling_space, levels;
        relation=:geq, r_min=0.0, r_max=nothing, mode=:boundaries, kwargs...)

Plot level intervals from `cycspace_level_intervals` over segment length and radius. `levels`
may be a scalar or an iterable. `mode=:boundaries` draws lower and upper interval boundaries;
`mode=:bands` draws each interval as a vertical segment.
"""
plot_cycspace_level_contours(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
plot_cycspace_level_contours!(args...; kwargs...) = error(_PLOTTING_EXT_ERROR)
