module PlotsPltExt

using CyclingSignatures
using Plots
using StatsPlots

function CyclingSignatures.plot_rank_distribution(results::RandomSubsegmentResult, rank; r_max=nothing, r_subdivisions=nothing)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end
    r_dists = rank_distribution(results, rank)
    radii = range(0, r_max; length = r_subdivisions)
    segment_lengths = results.segment_lengths

    Z = [r_dists[j](r) for r in radii, j in eachindex(segment_lengths)]

    return heatmap(
        segment_lengths,
        radii,
        Z;
        xlabel = "Segment length",
        ylabel = "Radius",
        colorbar_title = "Frequency",   # or whatever the StepFunction represents
        title = "Rank distribution heatmap",
    )
end

function _rank_heatmap_data(results, rank; r_max, r_subdivisions)
    segment_lengths = results.segment_lengths
    radii = range(0, r_max; length = r_subdivisions)

    r_dists = rank_distribution(results, rank)
    Z = [r_dists[j](r) for r in radii, j in eachindex(segment_lengths)]

    return segment_lengths, radii, Z
end

function CyclingSignatures.plot_rank_heatmap(results, rank; r_max=nothing, r_subdivisions=nothing)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end
    segment_lengths, radii, Z =
        _rank_heatmap_data(results, rank; r_max=r_max, r_subdivisions=r_subdivisions)

    heatmap(segment_lengths, radii, Z;
        xlabel="Segment length",
        ylabel="Radius",
        title="Rank $rank",
        colorbar_title="Value",
    )
end

function CyclingSignatures.plot_all_rank_heatmaps(results; r_max=nothing, r_subdivisions=nothing)
    if r_max === nothing
        r_max = results.trajectory_space.flt_max_heuristic
    end
    if r_subdivisions === nothing
        r_subdivisions = length(results.segment_lengths)
    end

    max_rank = betti_1(results.trajectory_space)

    plots = map(0:max_rank) do rank
        seg_lengths, radii, Z =
            _rank_heatmap_data(results, rank; r_max=r_max, r_subdivisions=r_subdivisions)

        heatmap(seg_lengths, radii, Z;
            xlabel="Segment length",
            ylabel="Radius",
            title="Rank $rank",
            colorbar_title="Value",
        )
    end

    plot(plots...; layout=(1, max_rank + 1))
end

function CyclingSignatures.plot_rank_distribution_at_r(results, evaluation_radius)
    segment_lengths = results.segment_lengths
    max_rank        = betti_1(results.trajectory_space)

    n_seg   = length(segment_lengths)
    n_ranks = max_rank + 1

    # M[i, j] = value at segment i, rank j
    M = zeros(n_seg, n_ranks)

    for (rank_idx, rank) in enumerate(0:max_rank)
        dists = rank_distribution(results, rank)
        @assert length(dists) == n_seg
        for (seg_idx, sf) in enumerate(dists)
            M[seg_idx, rank_idx] = sf(evaluation_radius)
        end
    end

    labels = permutedims(["Rank $r" for r in 0:max_rank])  # 1×n_ranks

    groupedbar(
        segment_lengths,
        M;
        bar_position = :stack,           # grouped bars
        xlabel       = "Segment length",
        ylabel       = "Value at radius = $evaluation_radius",
        label        = labels,
        legend       = :topright,
        title        = "Rank distribution at r = $evaluation_radius"
    )
end

function CyclingSignatures.plot_subspace_frequency_at_r(results::RandomSubsegmentResult, rank, evaluation_radius; n_subspaces=nothing)
    sig, M = cycspace_length_countmatrix_at_r(results, rank, evaluation_radius, n_subspaces=n_subspaces)
    segment_lengths = results.segment_lengths
    subspace_labels = ["Subspace $i" for i in 1:size(M, 2)]

    plt = plot(
        xlabel = "Segment length",
        ylabel = "Frequency",
        title  = "$(rank)d cyc. spaces at r = $evaluation_radius",
        legend = :topright,
    )

    # each row of M is one subspace → one line
    for (i, row) in enumerate(eachrow(M))
        plot!(plt, segment_lengths, collect(row); label = subspace_labels[i])
    end

    return sig, plt

    return sig, plt
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

function CyclingSignatures.plot_cycspace_inclusion(V0, V1)
    M = cycspace_inclusion_matrix(V0, V1)
    m, n = size(M)

    xs_inc = subspace_inclusion_points(m, n)
    x_1d = xs_inc[1]
    x_2d = xs_inc[2]

    y_1d = fill(1.0, m)
    y_2d = fill(2.0, n)

    p = plot(
        legend  = :topright,
        xlabel  = "Cyc. space index",
        ylabel  = "Dimension",
        yticks  = ([1, 2], ["1D", "2D"]),
        xaxis   = false,
        aspect_ratio = 1,
    )

    # -------------------------------------------------------------------------
    # Build line segments via NaN separators
    lx = Float64[]
    ly = Float64[]

    for i in 1:m
        for j in 1:n
            if M[i,j]
                append!(lx, [x_1d[i], x_2d[j], NaN])
                append!(ly, [1.0,      2.0,      NaN])
            end
        end
    end

    plot!(p, lx, ly;
        seriestype = :path,
        color      = :gray,
        alpha      = 0.5,
        lw         = 1,
        label      = "Inclusion",
    )

    # -------------------------------------------------------------------------
    # Draw nodes: all markers as circles
    scatter!(p, x_1d, y_1d;
        markershape = :circle,
        markersize  = 15,
        label       = "1D cyc. spaces",
    )
    scatter!(p, x_2d, y_2d;
        markershape = :circle,
        markersize  = 15,
        label       = "2D cyc. spaces",
    )

    # -------------------------------------------------------------------------
    # Add index numbers next to each marker
    for i in 1:m
        annotate!(p, x_1d[i] + 0*0.15, 1.0, text("$i", :black, 8))
    end
    for j in 1:n
        annotate!(p, x_2d[j] + 0*0.15, 2.0, text("$j", :black, 8))
    end

    return p
end

end
