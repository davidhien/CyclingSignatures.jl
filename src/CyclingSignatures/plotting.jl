### rank plots

function rankMatToStackedBarplotData(xs, rankMat)
    m,n = size(rankMat)
    xs_barplot = repeat(xs, inner=m)
    ys_barplot = vcat(eachcol(rankMat)...)
    color_barplot = repeat(1:m, outer=n)
    stack_barplot = repeat(1:m, inner=n)
    return xs_barplot, ys_barplot, color_barplot, stack_barplot
end

function getRankColors(k=6)
    slightBend(x,p) = x -p*x*(x-1)
    rpColors = ColorSchemes.magma[slightBend.(LinRange(0,1,k),.75)]
    return RGBA.(rpColors,1.0)
end

function plotRanks(results::Vector, r; gp = nothing, colors=nothing, axis_kwargs=(;), barplot_kwargs=(;gap=0, strokewidth=.2))
    segLengths = getSegmentLength.(results)
    rankMat = rankDistributionMatrix(results, r)
    xs, ys, c, st = rankMatToStackedBarplotData(segLengths, rankMat)

    if gp === nothing
        fig = Figure()
        gp = fig[1,1]
    end

    nRanks = size(rankMat,2)
    if colors === nothing
        colors = getRankColors(nRanks)
    end
    Axis(gp; axis_kwargs...)
    barplot!(gp, xs,ys; color=colors[c], stack=st, barplot_kwargs...)

    return gp.layout.parent
end

function plotRanksWithLegend(results::Vector, r;gl = nothing, colors = nothing, legend_kwargs=(;), kwargs...)
    if gl === nothing
        fig = Figure()
        gl = fig.layout
    end
    cms = rankCountmap.(results, r)
    nRanks = maximum(cm -> maximum(keys(cm)), cms)+1
    if colors === nothing
        colors = getRankColors(nRanks)
    end

    plotRanks(results, r; gp=gl[1,1], colors=colors, kwargs...)

    def_legend_kwa = (; orientation = :vertical, framevisible = true, titleposition=:left, nbanks=3,patchsize=(10,10), padding=(8f0,8f0,6f0,6f0),tellwidth=false,tellheight=true)
    

    elements = [PolyElement(polycolor = c) for c in colors[1:nRanks]]
    labels = string.(0:nRanks-1)
    Makie.Legend(gl[2,1], elements, labels, "rank"; merge(def_legend_kwa, legend_kwargs)...)    

    return gl.parent
end

### signature plots

function getSig1Colors(k=6)
    slightBend(x,p) = x -p*x*(x-1)
    return ColorSchemes.viridis[slightBend.(LinRange(0,1,max(k,2)),.75)]
end

function getSig2Colors(k=6)
    slightBend(x,p) = x -p*x*(x-1)
    return ColorSchemes.hawaii[slightBend.(LinRange(0,1,max(k,2)),.75)]
end

function plotSignatures(results::Vector, rk, r; cutoff=0, gp = nothing, colors = nothing,axis_kwargs=(;),stairs_kwargs=(;) )
    xs = getSegmentLength.(results)
    sig, sig_mat = subspaceFrequencyMatrix(results,rk, r; cutoff=cutoff)

    if gp === nothing
        fig = Figure()
        gp = fig[1,1]
    end

    nSig = length(sig)
    if colors === nothing
        colors = getSig1Colors(nSig)
    end
    Axis(gp; axis_kwargs...)
    foreach(enumerate(zip(eachrow(sig_mat),colors))) do t
        i=t[1]
        v,c=t[2]
        stairs!(gp, xs, v, color=c, step=:post, label=L"v_{%$i}"; stairs_kwargs...)
    end

    return gp.layout.parent
end

function plotSignaturesWithLegend(results::Vector, rk, r; gl = nothing, colors = nothing, legend_kwargs=(;), kwargs...)
    if gl=== nothing
        fig = Figure()
        gl = fig.layout
    end

    plotSignatures(results, rk, r; gp=gl[1,1], colors=colors, kwargs...)

    def_legend_kwa = (; orientation = :vertical, framevisible = true, titleposition=:left, nbanks=3, patchsize=(10,10), padding=(8f0,8f0,6f0,6f0), tellwidth=false,tellheight=true)
    Makie.Legend(gl[2,1], content(gl[1,1]), "rank $(rk)\n signature "; merge(def_legend_kwa, legend_kwargs)...)

    return gl.parent
end

### inclusion graph plots

function subspaceInclusionPoints(m,n)
    xs1 = collect(Float64,1:n)
    xs2 = collect(Float64,1:m)
    if m>n
        xs1 .+= (m-n)/2
    elseif n>m
        xs2 .+= (n-m)/2
    end
    return xs1, xs2
end

function plotSubspaceInclusion(ax, inc_mat, colors1, colors2)
    linesegments!(ax, lines_mat[1,:], lines_mat[2,:], color=:black, linewidth=1)

    scatter!(ax, xs1, ones(length(xs1)), color=colors1[1:length(xs1)], markersize=17)
    scatter!(ax, xs1, ones(length(xs1)), color=:white, markersize=12)

    scatter!(ax, xs2, 2*ones(length(xs2)), color=colors2[1:length(xs2)], markersize=17)
    scatter!(ax, xs2, 2*ones(length(xs2)), color=:white, markersize=12)

    markers1 = [L"v_%$i" for i = 1:length(xs1)]
    markers2 = [L"w_%$i" for i = 1:length(xs1)]
    foreach(zip(xs1, ones(length(xs1)),markers1)) do t
        x,y,m = t
        text!(ax, x, 0.018+y, text=m, align=(:center,:center))
    end
    foreach(zip(xs2, 2*ones(length(xs2)),markers2)) do t
        x,y,m = t
        text!(ax, x, 0.018+y, text=m, align=(:center,:center))
    end

    return ax
end

function plotSubspaceInclusion(results::Vector, rk1, rk2, r,field=FF{2}; cutoff1=0, cutoff2=0, gp=nothing,colors1 = nothing, colors2 = nothing, axis_kwargs=(;),linesegments_kwargs=(;),kwargs_scatter_outer=(;),kwargs_scatter_inner=(;),kwargs_text=(;))
    sig1, _ = subspaceFrequencyMatrix(results,rk1, r; cutoff=cutoff1)
    sig2, _ = subspaceFrequencyMatrix(results,rk2, r; cutoff=cutoff2)

    inc_mat = subspaceInclusionMatrix(sig1, sig2, field)

    if colors1 === nothing
        colors1 = getSig1Colors(length(sig1))
    end
    if colors2 === nothing
        colors2 = getSig1Colors(length(sig2))
    end
    if gp === nothing
        fig = Figure()
        gp = fig[1,1]
    end
    def_axis_kwargs = (;xticks=1:max(length(sig1),length(sig2)) ,yticks=1:2)
    ax = Axis(gp;merge(def_axis_kwargs,axis_kwargs)...)

    # actual plotting
    m,n = size(inc_mat)
    xs1, xs2 = subspaceInclusionPoints(m,n)

    line_ind = findall(==(1), inc_mat)
    lines_mat = hcat(map(line_ind) do i
        return [xs1[i[2]] xs2[i[1]]; 1 2]
    end...)

    def_arg_linesegments = (; color=:black, linewidth=1)
    def_arg_scatter_outer = (; markersize=40)
    def_arg_scatter_inner = (; color=:white, markersize=30)

    linesegments!(ax, lines_mat[1,:], lines_mat[2,:]; merge(def_arg_linesegments,linesegments_kwargs)...)

    scatter!(ax, xs1, ones(length(xs1)), color=colors1[1:length(xs1)]; merge(def_arg_scatter_outer,kwargs_scatter_outer)...)
    scatter!(ax, xs1, ones(length(xs1)); merge(def_arg_scatter_inner,kwargs_scatter_inner)...)

    scatter!(ax, xs2, 2*ones(length(xs2)), color=colors2[1:length(xs2)]; merge(def_arg_scatter_outer,kwargs_scatter_outer)...)
    scatter!(ax, xs2, 2*ones(length(xs2)); merge(def_arg_scatter_inner,kwargs_scatter_inner)...)

    markers1 = [L"v_%$i" for i = 1:length(xs1)]
    markers2 = [L"w_%$i" for i = 1:length(xs1)]
    foreach(zip(xs1, ones(length(xs1)),markers1)) do t
        x,y,m = t
        text!(ax, x, y, text=m, align=(:center,:center); kwargs_text...)
    end
    foreach(zip(xs2, 2*ones(length(xs2)),markers2)) do t
        x,y,m = t
        text!(ax, x, y, text=m, align=(:center,:center);kwargs_text...)
    end
    return gp.layout.parent
end
