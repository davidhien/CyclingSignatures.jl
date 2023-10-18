function plotDadrasSlice(X, slice=0, slice_dim=4, scene=nothing; slice_color=:nothing, show_axis=true, kwargs...)
    fig = nothing
    if scene == nothing
        fig = Figure(resolution=(800,600))
        scene = LScene(fig[1, 1], scenekw = (camera = cam3d!, raw = false), show_axis=true)
    end
    plot_dims = [collect(1:slice_dim-1); collect(slice_dim+1:4)]
    plot_timesteps = filter(x-> X[slice_dim,x] == slice, 1:size(X,2))
    if slice_color != :nothing
        GLMakie.scatter!(scene, X[plot_dims,plot_timesteps]', color = slice_color[plot_timesteps], shading = false; kwargs...)
    else
        GLMakie.scatter!(scene, X[plot_dims,plot_timesteps]', shading = false; kwargs...)
    end
    #update_cam!(scene.scene, Vec3f(13,46,7), Vec3f(0, 0, 0))
    # TODO: fix weir return behavior
    if fig === nothing
        cameracontrols(content(scene).scene).attributes.fov = 45
        return scene
    else
        cameracontrols(content(fig[1,1]).scene).attributes.fov = 45
        return fig
    end
end

function plotDadrasCoordinate(X, c,slice=0, slice_dim=4)
    fig = Figure(resolution=(800,600))
    plot_dims = [collect(1:slice_dim-1); collect(slice_dim+1:4)]
    plot_timesteps = filter(x-> X[slice_dim,x] == slice, 1:size(X,2))
    lscene = LScene(fig[1, 1], scenekw = (camera = cam3d!, raw = false), show_axis=false)
    GLMakie.scatter!(lscene, X[plot_dims,plot_timesteps]', colormap=:balance, colorrange=(0,1), color = c[plot_timesteps],
            markersize = 15, shading = false)
#    update_cam!(fig.content[1].scene, Vec3f(13,46,7), Vec3f(0, 0, 0))
    update_cam!(content(fig[1,1]).scene, Vec3f(8.75,25,-1), Vec3f(0,0.5,0))
    cameracontrols(content(fig[1,1]).scene).attributes.fov = 35
    return fig
end

function plotDadrasTSWithGenerator(X, boxsize, time_series, gen)
    figtsgen = plotDadrasSlice(X, markersize=5, show_axis=false)
    scatter!(figtsgen[1,1], 1/boxsize*time_series[:,1:end], markersize=5, color=(:orange,.5))
    cameracontrols(content(figtsgen[1,1]).scene).attributes.fov[] = 35
    update_cam!(content(figtsgen[1,1]).scene, Vec3f(8.75,25,-1), Vec3f(0,0.5,0))
    plotEdges(figtsgen[1,1], 1/boxsize*time_series, ripsererBarToEdgeList(gen), color=:red, linewidth=4, plot_dims=1:3)
    return figtsgen
end