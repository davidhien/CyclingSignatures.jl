"""
    DynamicDistance{C} <: Metric

Calculates the distance between points in the tangent bundle according to the formula

``d_{dyn}((p,v),(q,w)) = \\max(\\|p - q\\|, C \\|v - w\\|)``

where `C` is a constant.

# Fields
- `dim::Int`: The dimension of the space (i.e. half the tangent space dimension).
- `c::C`: The constant used in the distance calculation.
"""
struct DynamicDistance{C} <: Metric
    dim::Int
    c::C
end

get_dim(dist::DynamicDistance) = getfield(dist, :dim)
get_c(dist::DynamicDistance) = getfield(dist, :c)

@inline function (dist::DynamicDistance)(x, y)
    d = get_dim(dist)
    c = get_c(dist)

    return max(euclidean(@view(x[1:d]), @view(y[1:d])), c*euclidean(@view(x[d+1:end]),@view(y[d+1:end])))
end

function result_type(f::DynamicDistance, a::Type, b::Type)
    d = f.dim
    x = fill(oneunit(a),2*d)
    y = fill(oneunit(b),2*d)
    return typeof(f(x,y))
end
