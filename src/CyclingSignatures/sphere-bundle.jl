"""
    function blowUpSphereBundle(pts)
Pts is a matrix with sorted columns
"""
function blowUpSphereBundle(pts, sb_radius)
    m = size(pts,1)
    d = div(m,2)
    nDispacements = m^3
    displacements = map(0:nDispacements) do i
        digits(i, base=3, pad=m) .- 1
    end
    col_iterator = eachcol(pts)
    displaced_iterators = map(displacements) do dp
        # here, a transducer for filtering after map could help
        Iterators.map(p->dp+p, col_iterator)
    end
    all_pts = foldl(union, displaced_iterators, init=collect(first(displaced_iterators)))
    all_pts = Vector{Vector{Int}}(all_pts) # NOTE: this is since eltype inference appears to not work for nested iterators.
    filter!(v-> maximum(abs.(v[d+1:2*d])) == sb_radius, all_pts)

    perm = sortperm(all_pts)
    return hcat(all_pts[perm]...), perm
end

function sphereBundleDistanceMatrix(SB_pts)
    d = div(size(SB_pts,1),2)
    SB_X = SB_pts[1:d,:]
    SB_VF= SB_pts[d+1:end,:]
    println(size(SB_X))
    @time dm = max.(pairwise(chebyshev, SB_X), pairwise(chebyshev, SB_VF))
    return dm
end

struct SBDistance <: Metric
    dim::Int
    c::Float64
end

getdim(dist::SBDistance) = getfield(dist, :dim)
getc(dist::SBDistance) = getfield(dist, :c)

@inline function (dist::SBDistance)(x, y)
    d = getdim(dist) 
    c = getc(dist)

    return max(euclidean(@view(x[1:d]), @view(y[1:d])), c*euclidean(@view(x[d+1:end]),@view(y[d+1:end])))
end

function result_type(f::SBDistance, a::Type, b::Type)
    d = f.dim
    x = fill(oneunit(a),2*d)
    y = fill(oneunit(b),2*d)
    return typeof(f(x,y))
end
# TODO: implement specialized pairwise method for this thing here.
