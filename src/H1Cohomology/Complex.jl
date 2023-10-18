"""
    function oneSkeleton(D::AbstractArray{T,2}, r) where T <: Number

Returns the 1-skeleton of the VR-complex with distance matrix D. D only need the lower left entries. 
"""
function oneSkeleton(D::AbstractArray{T,2}, r) where T <: Number
    n = size(D,1)
    cells = Dict{Int,Vector{Simplex}}()
    cells[0] = map(x ->Simplex([x]), collect(1:n))[:]

    # fix 1-cells
    cells[1] = Simplex[]
    sizehint!(cells[1], div(count(<=(r), D),2))
    for j = 1:n
        for i = (j+1):n
            if D[i,j]<=r # i > n
                push!(cells[1], Simplex([j;i]))
            end
        end
    end

    return cells
end

"""
    function vrIncremental(X::AbstractArray{T,2}, d, r; maxdim=size(X,1)) where T <: Number

Compute VR complex for a given point cloud X, distance function d and radius r.
"""
function vrIncremental(X::AbstractArray{T,2}, d, r; maxdim=2) where T <: Number
    D = pairwise(d, X, dims=2) # can maybe be improved by only computing lower left part. 
    return vrIncremental(D, r, maxdim=maxdim)
end

"""
    function vrIncremental(D::AbstractArray{T,2}, r; maxdim=size(D,1)) where T <: Number

Constructs the VR-complex using the incremental construction. 
"""
function vrIncremental(D::AbstractArray{T,2}, r; maxdim=2) where T <: Number
    n = size(D,1)
    cells = oneSkeleton(D,r)
    
    # TODO: potential speed up: instead of creating a vector of simplices, first just have an integer matrix and then convert
    for i = 2:maxdim
        cells[i] = Vector{Simplex}()

        for cell in cells[i-1]
            lastVertex = cell.vertices[end]
            for j = lastVertex+1:n
                if all(<=(r), (@view D[cell.vertices, j]) )
                    push!(cells[i],Simplex([cell.vertices;j]))
                end
            end
        end
    end
    return CellComplex(cells)
end
