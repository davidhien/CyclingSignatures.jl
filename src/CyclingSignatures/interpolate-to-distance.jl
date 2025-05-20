function getInterpolationFunction(Y)
    curves = map(eachrow(Y)) do r
        return CubicSpline(r, 1:size(Y,2))
    end
    return t -> map(c-> c(t), curves), t -> map(c->DataInterpolations.derivative(c,t), curves)
end

function interpolateToDistance(Y, r, d=euclidean; norm_sb=norm, sb_r=nothing, verbose = false, max_depth=50000)
    v = Vector{Matrix{Float64}}(undef, size(Y,2))
    w = Vector{Matrix{Float64}}(undef, size(Y,2))
    c, d_c = getInterpolationFunction(Y)
    for i = 1:size(Y,2)-1
        v[i],w[i] = interpolateToDistance(c, d_c, i, i+1, r, d; norm_sb=norm_sb, sb_r=sb_r, verbose = verbose, max_depth=max_depth)
    end

    v[end] = Y[:,end][:,:]
    w[end] = d_c(size(Y,2))[:,:]
    new_lengths = map(m->size(m,2),v)
    t_vec = zeros(Int,length(new_lengths)+1)
    t_vec[1] = 1
    for i = 2:length(t_vec)
        t_vec[i] = t_vec[i-1] + new_lengths[i-1]
    end
    Y_new = reduce(hcat,v)
    Z_new = reduce(hcat, w)
    return Y_new, Z_new, t_vec
end

function interpolateToDistance(c, d_c, a, b, r, d; norm_sb=norm, sb_r=nothing, verbose=false, max_depth=50000)
    # ec is end condition
    function isSpaceConsistent(x,y)
        return d(x, y)<=r
    end
    function isSBConsistent(x,y)
        if isnothing(sb_r)
            return true
        end
        return norm_sb(x/norm_sb(x)-y/norm_sb(y)) <= sb_r
    end

    counter = 1 # count to not end in inf loop

    t_vec = Float64[a]

    ec = c(b)
    d_ec = d_c(b)

    while !(isSpaceConsistent(c(t_vec[end]),ec) && isSBConsistent(d_c(t_vec[end]),d_ec)) && (counter += 1) < max_depth
        t_vec = [t_vec b]
        v = 2
        while !(isSpaceConsistent(c(t_vec[end-1]), c(t_vec[end])) && isSBConsistent(d_c(t_vec[end-1]), d_c(t_vec[end])))
            t_vec[end] = t_vec[end-1] + 1/v*(b-t_vec[end-1])

            v *= 2
            if v > max_depth
                @warn "Went over max depth"
                break
            end
            verbose && print(v)
        end
    end
    if counter >= max_depth
        @warn "Counter exhausted!"
    end
    verbose && println()
    X = reduce(hcat, map(t->c(t), t_vec))
    X_bar = reduce(hcat,map(t->d_c(t), t_vec))
    return X[:,:],X_bar[:,:]
end

