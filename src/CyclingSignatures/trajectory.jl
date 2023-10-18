"""
quantizes each column of x to a grid of size r.

If r is an integer, boxes have size r in each dimension,
if r is a vector, it has to contain a length for every dimension.
"""
function quantize(x, r)
    return round.(Int, x./r)
end

"""
Quantizes with s giving the box size per dimension, then maps
points outside [q_min,q-max] to nearest box
"""
function quantize_subset(data,s,q_min,q_max)
    dataq = quantize(data,s)
    dataq = mapslices(y -> max.(q_min,y), dataq, dims=1)
    dataq = mapslices(y -> min.(q_max,y), dataq, dims=1)

    return dataq
end

"""
    function dynamicIndices(dataq::Matrix)

Returns a vector v and a matrix M s.t. the columns of M are pairwise distinct
and M[:,v[t]] = dataq[:,t] for all t.
In other words, M contains the unique boxes and v[t] the index of the box at
time step t.
"""
function dynamicIndices(dataq::Matrix)
    # TODO: this is probably not optimal
    # each column is a data-point
    n = size(dataq,2)
    indices = zeros(Int,n)
    datamap = Dict{Vector{Int}, Int}()

    m = 1
    for i = 1:n
        v = dataq[:,i]
        if haskey(datamap, v)
            indices[i] = datamap[v]
        else
            datamap[v] = m
            indices[i] = m
            m = m+1
        end
    end

    return indices, unique(dataq, dims = 2)
end

function sortedDynamicIndices(dataq::Matrix)
    t, X = dynamicIndices(dataq)
    perm = sortperm(collect(eachcol(X)))
    return invperm(perm)[t], X[:,perm]
end

function boxTrajectory(t)
    t_new = [t[1]]
    t_map = [1]

    for i = 2:length(t)
        if t[i] != t[i-1]
            push!(t_new, t[i])
            push!(t_map, i)
        end
    end
    return t_new, t_map
end

function quantizeTimeSeries(x, r)
    Xt = quantize(x, r)
    t, X = dynamicIndices(Xt)

    return t, X
end

"""
    function resampleToDistance(ds, Y, r, dt; kwargs...)

Adds additional sample points to a time series Y from a dynamical system ds at time step dt such that the resulting time series satisfies the specified distance requirements

Instead of the actual trajectory points, one can also require pp(Y) to satisfy the requirement, where pp is a postprocessing applied to the data (for example a rescaling).

Furthermore, it is possible to specify a sb_radius and a sb_fct which ensures consistency in the sphere bundle.
Note: this works as follows: the output of the integrator is applied to sb_fct, this will be scaled to sbRadius (in l_Inf norm) and then checked for consistency

"""
function resampleToDistance(ds, dt::Number, Y, r, d=euclidean; pp=identity, norm_sb=norm, sb_r=nothing, sb_fct=nothing, verbose = false, max_depth=512)
    v = Vector{Matrix{Float64}}(undef, size(Y,2))
    #dsint = integrator(ds, Y[:,end-1], diffeq = (dt = dt,))
    for i = 1:size(Y,2)-1
        v[i] = resampleToDistance(ds, dt, Y[:,i], Y[:,i+1], r, d; pp=pp, norm_sb=norm_sb, sb_r=sb_r, sb_fct=sb_fct, verbose = verbose, max_depth=max_depth)
    end
    v[end] = Y[:,end][:,:]
    new_lengths = map(m->size(m,2),v)
    t_vec = zeros(Int,length(new_lengths)+1)
    t_vec[1] = 1
    for i = 2:length(t_vec)
        t_vec[i] = t_vec[i-1] + new_lengths[i-1]
    end
    Y_new = reduce(hcat,v)
    return Y_new, t_vec
end

function resampleToDistance(ds, dt::Number, ic, ec, r, d; pp=identity, norm_sb=norm, sb_r=nothing, sb_fct=nothing, verbose=false, max_depth=512)
    # ec is end condition
    function isSpaceConsistent(x,y)
        return d(pp(x), pp(y))<=r
    end
    function isSBConsistent(x,y)
        if isnothing(sb_r)
            return true
        end
        x_sb = sb_fct(x)
        y_sb = sb_fct(y)
        return norm_sb(x_sb/norm_sb(x_sb)-y_sb/norm_sb(y_sb)) <= sb_r
    end
    X = ic[:,:]
    ct = 1

    while !(isSpaceConsistent(X[:,end],ec) && isSBConsistent(X[:,end], ec)) && (ct += 1) < 100*max_depth
        X = [X ec]
        v = 2
        while !(isSpaceConsistent(X[:,end-1], X[:,end]) && isSBConsistent(X[:,end-1], X[:,end]))
            reinit!(ds, X[:,end-1])
            step!(ds, dt/v, true)
            X[:,end] = get_state(ds)
            v *= 2
            if v > max_depth
                @warn "Went over max depth"
                break
            end
            verbose && print(v)
        end
    end
    if ct >= 100*max_depth
        @warn "Counter exhausted!"
    end
    verbose && println()
    return X
end

"""
    function resampleToConsistent(ds, Y, r, dt)

Adds additional sample points to a time series Y from a dynamical system ds at time step dt such that
the resulting time series quantized at r is consistent.

Instead of the actual trajectory points, one can also make pp(Y) dynamically consistent, 
where pp is a postprocessing applied to the data (for example a rescaling).

Furthermore, it is possible to specify a sb_radius and a sb_fct which ensures consistency in the sphere bundle.
Note: this works as follows: the output of the integrator is applied to sb_fct, this will be scaled to sbRadius (in l_Inf norm) and then checked for consistency
"""
function resampleToConsistent(ds, Y, r, dt; pp=identity, sb_radius=nothing, sb_fct=nothing, verbose = false, max_depth=512)
    # TODO: reorder arguments ds,dt,Y,r
    v = Vector{Matrix{Float64}}(undef, size(Y,2))
    for i = 1:size(Y,2)-1
        v[i] = resampleToConsistent(ds, Y[:,i], Y[:,i+1], r, dt; pp=pp, sb_radius=sb_radius, sb_fct=sb_fct, verbose = verbose, max_depth=max_depth)
    end
    v[end] = Y[:,end][:,:]

    new_lengths = map(m->size(m,2),v)
    t_vec = zeros(Int,length(new_lengths)+1)
    t_vec[1] = 1
    for i = 2:length(t_vec)
        t_vec[i] = t_vec[i-1] + new_lengths[i-1]
    end
    Y_new = reduce(hcat,v)
    return Y_new, t_vec
end

"""
    function resampleToConsistent(ds, ic, ec, r, dt; pp=identity, verbose=false, max_depth=16)

Returns a matrix M such that [M ec] is dynamically consistent and M[:,1] = ic. Note that M[:,end] != ec.
"""
function resampleToConsistent(ds, ic, ec, r, dt; pp=identity, sb_radius=nothing, sb_fct=nothing, verbose=false, max_depth=512)
    # ec is end condition
    X = ic[:,:]
    ct = 1
    function sb_consistency(x,y)
        if sb_radius == nothing
            return true
        end
        x_sb = sb_fct(x)
        y_sb = sb_fct(y)
        x_scaled = sb_radius*x_sb/norm(x_sb,Inf)
        y_scaled = sb_radius*y_sb/norm(y_sb,Inf)
        return isDynConsistent(x_scaled, y_scaled, 1)
    end

    while !(isDynConsistent(pp(X[:,end]), pp(ec), r) && sb_consistency(X[:,end], ec)) && (ct += 1) < 100*max_depth
        X = [X ec]
        v = 2
        while !(isDynConsistent(pp(X[:,end-1]), pp(X[:,end]), r) && sb_consistency(X[:,end-1], X[:,end]))
            reinit!(ds, X[:,end-1]);
            step!(ds, dt/v, true)
            X[:,end] = get_state(ds)
            v *= 2
            if v > max_depth
                @warn "Went over max depth"
                break
            end
            verbose && print(v)
        end
    end
    if ct >= 100*max_depth
        @warn "Counter exhausted!"
    end
    verbose && println()
    return X
end

function isDynConsistent(p1,p2,r)
    p1q = round.(Int, p1/r)
    p2q = round.(Int, p2/r)

    return norm(p2q-p1q, Inf) <= 1
end

function countDynamicInconsistencies(boxIt, d=1)
    boxIt = Iterators.Stateful(boxIt)
    a = peek(boxIt)
    if a === nothing
        return 0
    end
    a = first(boxIt)
    return count(boxIt) do b
        isInconsistent = chebyshev(a,b)>d
        a = b
        return isInconsistent
    end
end

function findInconsistentBoxes(X, t, d=1)
    Xt = X[:,t]
    diff = Xt[:,2:end] - Xt[:,1:end-1]
    isConsistent = mapslices(x -> norm(x,Inf)>d, diff, dims=[1])[:]
    @show size([isConsistent;false])
    return unique( [t[[isConsistent;false]] t[[false;isConsistent]]], dims=2)
end

function quantizeSB(VF_t, r)
    VF_t_rescaled = mapslices(v -> r*v/norm(v,Inf), VF_t, dims=1)
    return quantize(VF_t_rescaled, 1)
end

function removeSuccessiveDuplicates(t)
    t_new = [t[1]]
    for i = 2:length(t)
        if t[i] != t[i-1]
            push!(t_new, t[i])
        end
    end
    return t_new
end