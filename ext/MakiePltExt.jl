module MakiePltExt

using CyclingSignatures
using Makie

function _categorical_colorscheme(colorscheme, n::Integer)
    if n <= 0
        return Any[]
    end
    return Makie.cgrad(colorscheme, n; categorical = true)
end

function _rank_heatmap_data(results, rank; r_max, r_subdivisions)
    segment_lengths = results.segment_lengths
    radii = range(0, r_max; length = r_subdivisions)

    r_dists = rank_distribution(results, rank)
    Z = [r_dists[j](r) for r in radii, j in eachindex(segment_lengths)]

    return segment_lengths, radii, Z
end

function CyclingSignatures.plot_rank_distribution!(
    gp,
    results::RandomSubsegmentResult,
    rank;
    r_max = nothing,
    r_subdivisions = nothing,
)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end

    segment_lengths, radii, Z =
        _rank_heatmap_data(results, rank; r_max = r_max, r_subdivisions = r_subdivisions)

    ax = Axis(
        gp;
        xlabel = "Segment length",
        ylabel = "Radius",
        title = "Rank distribution heatmap",
    )
    hm = heatmap!(ax, segment_lengths, radii, Z'; colormap = :viridis)
    return ax, hm
end

function CyclingSignatures.plot_rank_distribution(
    results::RandomSubsegmentResult,
    rank;
    r_max = nothing,
    r_subdivisions = nothing,
)
    fig = Figure()
    _, hm = plot_rank_distribution!(
        fig[1, 1],
        results,
        rank;
        r_max = r_max,
        r_subdivisions = r_subdivisions,
    )
    Colorbar(fig[1, 2], hm; label = "Frequency")
    return fig
end

function CyclingSignatures.plot_rank_heatmap!(
    gp,
    results,
    rank;
    r_max = nothing,
    r_subdivisions = nothing,
)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end

    segment_lengths, radii, Z =
        _rank_heatmap_data(results, rank; r_max = r_max, r_subdivisions = r_subdivisions)

    ax = Axis(
        gp;
        xlabel = "Segment length",
        ylabel = "Radius",
        title = "Rank $rank",
    )
    hm = heatmap!(ax, segment_lengths, radii, Z'; colormap = :viridis)
    return ax, hm
end

function CyclingSignatures.plot_rank_heatmap(results, rank; r_max = nothing, r_subdivisions = nothing)
    fig = Figure()
    _, hm = plot_rank_heatmap!(
        fig[1, 1],
        results,
        rank;
        r_max = r_max,
        r_subdivisions = r_subdivisions,
    )
    Colorbar(fig[1, 2], hm; label = "Value")
    return fig
end

function CyclingSignatures.plot_all_rank_heatmaps!(gp, results; r_max = nothing, r_subdivisions = nothing)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end

    max_rank = betti_1(results.trajectory_space)
    ranks = 0:max_rank

    heatmap_data = map(ranks) do rank
        _rank_heatmap_data(results, rank; r_max = r_max, r_subdivisions = r_subdivisions)
    end

    zmax = maximum(last.(heatmap_data)) do Z
        maximum(Z)
    end
    colorrange = (0.0, Float64(zmax))

    first_hm = nothing
    for (i, rank) in enumerate(ranks)
        segment_lengths, radii, Z = heatmap_data[i]
        ax = Axis(
            gp[1, i];
            xlabel = "Segment length",
            ylabel = (i == 1 ? "Radius" : ""),
            title = "Rank $rank",
        )
        hm = heatmap!(ax, segment_lengths, radii, Z'; colormap = :viridis, colorrange = colorrange)
        first_hm === nothing && (first_hm = hm)
    end
    Colorbar(gp[1, length(ranks) + 1], first_hm; label = "Value")
    return gp
end

function CyclingSignatures.plot_all_rank_heatmaps(results; r_max = nothing, r_subdivisions = nothing)
    fig = Figure()
    plot_all_rank_heatmaps!(fig[1, 1], results; r_max = r_max, r_subdivisions = r_subdivisions)
    return fig
end

function CyclingSignatures.plot_rank_distribution_at_r!(
    gp,
    results,
    evaluation_radius;
    legend = true,
    legend_position = :rt,
    colorscheme = :magma,
    markersize = 6,
)
    segment_lengths = results.segment_lengths
    max_rank = betti_1(results.trajectory_space)

    n_seg = length(segment_lengths)
    n_ranks = max_rank + 1

    M = zeros(n_seg, n_ranks)
    for (rank_idx, rank) in enumerate(0:max_rank)
        dists = rank_distribution(results, rank)
        @assert length(dists) == n_seg
        for (seg_idx, sf) in enumerate(dists)
            M[seg_idx, rank_idx] = sf(evaluation_radius)
        end
    end

    rank_colors = _categorical_colorscheme(colorscheme, n_ranks)

    ax = Axis(
        gp;
        xlabel = "Segment length",
        ylabel = "Value at radius = $evaluation_radius",
        title = "Rank distribution at r = $evaluation_radius",
    )

    positions = repeat(collect(segment_lengths), inner = n_ranks)
    heights = vec(M')
    stack = repeat(collect(1:n_ranks), outer = n_seg)
    colors = [rank_colors[i] for i in stack]

    bp = barplot!(ax, positions, heights; stack = stack, color = colors)

    legend_elements = [PolyElement(color = rank_colors[i]) for i in 1:n_ranks]
    legend_labels = ["Rank $(i - 1)" for i in 1:n_ranks]
    leg = legend ? axislegend(ax, legend_elements, legend_labels; position = legend_position) : nothing
    return ax, bp, leg
end

function CyclingSignatures.plot_rank_distribution_at_r(
    results,
    evaluation_radius;
    legend = true,
    legend_position = :rt,
    colorscheme = :magma,
    markersize = 6,
)
    fig = Figure()
    plot_rank_distribution_at_r!(
        fig[1, 1],
        results,
        evaluation_radius;
        legend = legend,
        legend_position = legend_position,
        colorscheme = colorscheme,
        markersize = markersize,
    )
    return fig
end

function CyclingSignatures.plot_subspace_frequency_at_r!(
    gp,
    results::RandomSubsegmentResult,
    rank,
    evaluation_radius;
    n_subspaces = nothing,
    legend = true,
    legend_position = :rt,
    colorscheme = :viridis,
    linewidth = 2,
)
    sig, M = cycspace_length_countmatrix_at_r(
        results,
        rank,
        evaluation_radius;
        n_subspaces = n_subspaces,
    )
    segment_lengths = results.segment_lengths

    n_subspaces_found = size(M, 1)
    colors = _categorical_colorscheme(colorscheme, n_subspaces_found)

    ax = Axis(
        gp;
        xlabel = "Segment length",
        ylabel = "Frequency",
        title = "$(rank)d cyc. spaces at r = $evaluation_radius",
    )

    plots = Any[]
    for (i, row) in enumerate(eachrow(M))
        p = lines!(
            ax,
            segment_lengths,
            collect(row);
            color = colors[i],
            linewidth = linewidth,
            label = "Subspace $i",
        )
        push!(plots, p)
    end
    leg = legend ? axislegend(ax; position = legend_position) : nothing
    return sig, ax, plots, leg
end

function CyclingSignatures.plot_subspace_frequency_at_r(
    results::RandomSubsegmentResult,
    rank,
    evaluation_radius;
    n_subspaces = nothing,
    legend = true,
    legend_position = :rt,
    colorscheme = :viridis,
    linewidth = 2,
)
    fig = Figure()
    sig, _, _, _ = plot_subspace_frequency_at_r!(
        fig[1, 1],
        results,
        rank,
        evaluation_radius;
        n_subspaces = n_subspaces,
        legend = legend,
        legend_position = legend_position,
        colorscheme = colorscheme,
        linewidth = linewidth,
    )
    return sig, fig
end

function subspace_inclusion_points(ms...)
    max_size = maximum(ms)
    return map(ms) do m
        xs = collect(1:m)
        if m < max_size
            xs .+= (max_size - m) / 2
        end
        return xs
    end
end

function CyclingSignatures.plot_cycspace_inclusion!(
    gp,
    V0,
    V1;
    legend = true,
    legend_position = :rt,
    markersize = 25,
    text_color = :black,
    text_fontsize = 10,
    colorscheme = :viridis,
)
    M = cycspace_inclusion_matrix(V0, V1)
    m, n = size(M)

    xs_inc = subspace_inclusion_points(m, n)
    x_1d = xs_inc[1]
    x_2d = xs_inc[2]

    y_1d = fill(1.0, m)
    y_2d = fill(2.0, n)

    colors = _categorical_colorscheme(colorscheme, 2)

    ax = Axis(
        gp;
        xlabel = "Cyc. space index",
        ylabel = "Dimension",
        yticks = ([1, 2], ["1D", "2D"]),
        xlabelvisible = false,
        xticklabelsvisible = false,
        xticksvisible = false,
        xgridvisible = false,
        xminorgridvisible = false,
        xminorticksvisible = false,
        aspect = DataAspect(),
    )

    lx = Float64[]
    ly = Float64[]
    for i in 1:m
        for j in 1:n
            if M[i, j]
                push!(lx, x_1d[i])
                push!(lx, x_2d[j])
                push!(ly, 1.0)
                push!(ly, 2.0)
            end
        end
    end

    inc = linesegments!(ax, lx, ly; color = (:gray, 0.5), linewidth = 1, label = "Inclusion")
    sc1 = scatter!(ax, x_1d, y_1d; markersize = markersize, color = colors[1], label = "1D cyc. spaces")
    sc2 = scatter!(ax, x_2d, y_2d; markersize = markersize, color = colors[2], label = "2D cyc. spaces")

    for i in 1:m
        text!(
            ax,
            "$i";
            position = (x_1d[i], 1.0),
            align = (:center, :center),
            color = text_color,
            fontsize = text_fontsize,
        )
    end
    for j in 1:n
        text!(
            ax,
            "$j";
            position = (x_2d[j], 2.0),
            align = (:center, :center),
            color = text_color,
            fontsize = text_fontsize,
        )
    end

    leg = legend ? axislegend(ax; position = legend_position) : nothing
    return ax, inc, sc1, sc2, leg
end

function CyclingSignatures.plot_cycspace_inclusion(
    V0,
    V1;
    legend = true,
    legend_position = :rt,
    markersize = 25,
    text_color = :black,
    text_fontsize = 10,
    colorscheme = :viridis,
)
    fig = Figure()
    plot_cycspace_inclusion!(
        fig[1, 1],
        V0,
        V1;
        legend = legend,
        legend_position = legend_position,
        markersize = markersize,
        text_color = text_color,
        text_fontsize = text_fontsize,
        colorscheme = colorscheme,
    )
    return fig
end

function CyclingSignatures.plot_cycspace_radius_frequency!(
    gp,
    results::RandomSubsegmentResult,
    k;
    r_max = nothing,
    r_max_for_sorting = nothing,
    n_subspaces = 10,
    n_highlight = nothing,
    filter_shorter_as = 0,
    colorscheme = :viridis,
    default_color = :gray60,
    linewidth = 1,
    legend = true,
    legend_position = :rt,
)
    if r_max === nothing
        r_max = results.flt_threshold
    end

    sig, fs = CyclingSignatures.cycspace_radius_frequency_functions(
        results,
        k;
        r_max_for_sorting = r_max_for_sorting,
        filter_shorter_as = filter_shorter_as,
        max_n_sig = n_subspaces,
    )

    ax = Axis(gp; xlabel = "thickening radius", ylabel = "#segments")
    if isempty(sig)
        return sig, ax, Any[], nothing
    end

    if n_highlight === nothing
        n_highlight = min(length(sig), 6)
    end

    highlight_colors = _categorical_colorscheme(colorscheme, n_highlight)

    plots = Any[]
    for i in eachindex(sig)
        xs, ys = CyclingSignatures._step_post_xy(fs[i], r_max; x_start = 0.0)
        if i <= n_highlight
            p = stairs!(
                ax,
                xs,
                ys;
                step = :post,
                color = highlight_colors[i],
                linewidth = linewidth,
                label = "V_$i",
            )
        else
            p = stairs!(
                ax,
                xs,
                ys;
                step = :post,
                color = default_color,
                linewidth = linewidth,
                label = "",
            )
        end
        push!(plots, p)
    end

    leg = (legend && n_highlight > 0) ? axislegend(ax; position = legend_position) : nothing
    return sig, ax, plots, leg
end

function CyclingSignatures.plot_cycspace_radius_frequency(
    results::RandomSubsegmentResult,
    k;
    r_max = nothing,
    r_max_for_sorting = nothing,
    n_subspaces = 10,
    n_highlight = nothing,
    filter_shorter_as = 0,
    colorscheme = :viridis,
    default_color = :gray60,
    linewidth = 1,
    legend = true,
    legend_position = :rt,
)
    fig = Figure()
    sig, _, _, _ = plot_cycspace_radius_frequency!(
        fig[1, 1],
        results,
        k;
        r_max = r_max,
        r_max_for_sorting = r_max_for_sorting,
        n_subspaces = n_subspaces,
        n_highlight = n_highlight,
        filter_shorter_as = filter_shorter_as,
        colorscheme = colorscheme,
        default_color = default_color,
        linewidth = linewidth,
        legend = legend,
        legend_position = legend_position,
    )
    return sig, fig
end

function CyclingSignatures.plot_cycspace_length_frequency!(
    gp,
    results::RandomSubsegmentResult,
    k;
    n_subspaces = 10,
    n_highlight = nothing,
    sort_by_tp_with_rmax = nothing,
    filter_shorter_than = 0,
    filter_shorter_r_max = Inf,
    colorscheme = :viridis,
    default_color = :gray60,
    markersize = 3,
    legend = true,
    legend_position = :rt,
)
    sig, M = CyclingSignatures.cycspace_length_interval_countmatrix(
        results,
        k;
        n_subspaces = n_subspaces,
        sort_by_tp_with_rmax = sort_by_tp_with_rmax,
        filter_shorter_than = filter_shorter_than,
        filter_shorter_r_max = filter_shorter_r_max,
    )
    segment_lengths = results.segment_lengths

    ax = Axis(gp; xlabel = "time span", ylabel = "#segments")
    if isempty(sig)
        return sig, ax, Any[], nothing
    end

    if n_highlight === nothing
        n_highlight = min(length(sig), 6)
    end

    highlight_colors = _categorical_colorscheme(colorscheme, n_highlight)

    plots = Any[]
    for i in 1:size(M, 1)
        if i <= n_highlight
            p = scatter!(
                ax,
                segment_lengths,
                vec(M[i, :]);
                color = highlight_colors[i],
                markersize = markersize,
                label = "V_$i",
            )
        else
            p = scatter!(
                ax,
                segment_lengths,
                vec(M[i, :]);
                color = default_color,
                markersize = markersize,
                label = "",
            )
        end
        push!(plots, p)
    end

    leg = (legend && n_highlight > 0) ? axislegend(ax; position = legend_position) : nothing
    return sig, ax, plots, leg
end

function CyclingSignatures.plot_cycspace_length_frequency(
    results::RandomSubsegmentResult,
    k;
    n_subspaces = 10,
    n_highlight = nothing,
    sort_by_tp_with_rmax = nothing,
    filter_shorter_than = 0,
    filter_shorter_r_max = Inf,
    colorscheme = :viridis,
    default_color = :gray60,
    markersize = 3,
    legend = true,
    legend_position = :rt,
)
    fig = Figure()
    sig, _, _, _ = plot_cycspace_length_frequency!(
        fig[1, 1],
        results,
        k;
        n_subspaces = n_subspaces,
        n_highlight = n_highlight,
        sort_by_tp_with_rmax = sort_by_tp_with_rmax,
        filter_shorter_than = filter_shorter_than,
        filter_shorter_r_max = filter_shorter_r_max,
        colorscheme = colorscheme,
        default_color = default_color,
        markersize = markersize,
        legend = legend,
        legend_position = legend_position,
    )
    return sig, fig
end

end
