struct PointIdentificationDistanceMatrix{S<:Integer, D} <: AbstractMatrix{Float64}
    pts::Matrix{S}
    identifications::Vector{Vector{S}}
    point_distance::D

    function PointIdentificationDistanceMatrix(pts::Matrix{S}, identifications::Vector{Vector{S}}, point_distance::D) where {S<:Integer, D}
        npts = size(pts, 2)

        @inbounds for (bi, v) in pairs(identifications)
            isempty(v) && throw(ArgumentError("identifications[$bi] is empty"))
            for idx in v
                (1 <= idx <= npts) || throw(ArgumentError(
                    "identifications[$bi] contains invalid index $idx; must be in 1:$npts"
                ))
            end
        end

        return new{S, D}(pts, identifications, point_distance)
    end
end

IndexStyle(::Type{<:PointIdentificationDistanceMatrix}) = IndexCartesian()
size(M::PointIdentificationDistanceMatrix) = (length(M.identifications), length(M.identifications))
axes(M::PointIdentificationDistanceMatrix) = (Base.OneTo(size(M,1)), Base.OneTo(size(M,2)))

@inline function getindex(M::PointIdentificationDistanceMatrix, i::Int, j::Int)::Float64
    pts = M.pts
    i1  = M.identifications[i]
    j1  = M.identifications[j]
    d   = M.point_distance

    v = Inf
    @inbounds @views for a in i1
        xa = view(pts, :, a)   # avoids allocation vs pts[:,a]
        for b in j1
            xb = view(pts, :, b)
            v = min(v, Float64(d(xa, xb)))
        end
    end
    return v
end

# Optional: linear indexing M[k]
@inline function getindex(M::PointIdentificationDistanceMatrix, k::Int)::Float64
    i, j = Base._ind2sub(axes(M), k)
    return getindex(M, i, j)
end
