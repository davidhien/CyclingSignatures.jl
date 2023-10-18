struct Cube <: AbstractCell
    base::Vector{Int}
    extent::Vector{Bool}

    function Cube(base::Vector{Int}, extent::Vector{Bool})
        if length(base) != length(extent)
            error("Length of base and length of extent must match")
        end
        return new(base, extent)
    end
end

mutable struct CubicalSet
    emb::Int
    cubes::Vector{Cube}

    CubicalSet(emb) = new(emb, Cube[])
end

# -----------------------------------------------------------------------------
# Start of methods for cubes

"""
    function dim(x::Cube)
Compute (true) dimension of an cube, aka the number of non-trivial intervals in its interval representation.
"""
function dim(x::Cube)
    return count(x.extent)
end

"""
    function emb(x::Cube)
Compute embedding dimension of an cube, aka the dim of the ambient space
"""
function emb(x::Cube)
    return length(x.base)
end

"""
    function stringRep(Q::Cube)

String representation of an cube as product of its intervals
"""
function stringRep(Q::Cube)
    if emb(Q) == 0
        return ""
    end
    s = ""
    for i = 1:length(Q.base)
        if i != 1
            s *= "x"
        end
        coord = Q.base[i]
        if Q.extent[i]
            s = string(s, "[", coord, ",", coord+1, "]")
        else
            s = string(s, "[", coord, "]")
        end
    end
    return s
end

function Base.show(io::IO, Q::Cube)
    return print(io, stringRep(Q))
end

"""
    function ==(Q1::Cube, Q2::Cube)

Two cubes are identical if
    - their ambient spaces have identical dimension,
    - they have the same elemntary  intervals.
"""
function Base.:(==)(Q1::Cube, Q2::Cube)
    return Q1.base == Q2.base && Q1.extent == Q2.extent
end

function Base.hash(Q::Cube, h::UInt)
    h = hash(Q.base, h)
    return hash(Q.extent, h)
end

function Base.isless(Q1::Cube, Q2::Cube)
    return (Q1.base,Q1.extent) < (Q2.base,Q2.extent)
end

"""
    function makeCube1(x::Int ...)

Convenience function for creating an cube.
A call with arguments a1, b1, a2, b2, ... generates the cube
[a1,b1]x[a2,b2]x...

TODO: example
"""
function makeCube1(x::Int ...)
    makeCube1(collect(x))
end

"""
    function makeCube1(x::Vector{Int})

Convenience function for creating an cube.
A call with arguments a1, b1, a2, b2, ... generates the cube
[a1,b1]x[a2,b2]x...

"""
function makeCube1(x::Vector{Int})
    if length(x)%2 == 1
        error("Interval needs start and end point even if they coincide")
    end
    emb = div(length(x), 2)

    base = zeros(Int, emb)
    extent = zeros(Bool, emb)

    for i = 1:emb
        base[i] = x[2*i-1]
        extent[i] = x[2*i-1] != x[2*i]
    end
    return Cube(base, extent)
end

"""
    function makeCube2(x::Int ...)

Convenience function for creating an cube.
A call with arguments a1, a2, ...,  b1, b2, ... generates the cube
[a1,b1]x[a2,b2]x...

TODO: example
"""
function makeCube2(x::Int ...)
    makeCube(collect(x))
end

"""
    function makeCube2(x::Vector{Int})

Convenience function for creating an cube.
A call with arguments a1, b1, a2, b2, ... generates the cube
[a1,b1]x[a2,b2]x...

TODO: example
"""
function makeCube2(x::Vector{Int})
    if length(x)%2 == 1
        error("Interval needs start and end point even if they coincide")
    end

    emb = div(length(x),2)
    base = Vector{Int}()
    extent = Vector{Bool}()
    for i = 1:emb
        base[i] = x[i]
        extent[i] = x[i] == x[i+emb]
    end

    return Cube(base, extent)
end

"""
    function *(Q1::Cube, Q2::Cube)

Given two elemetary cubes Q1 and Q2, returns the product Q1 x Q2.
"""
function Base.:*(Q1::Cube, Q2::Cube)
    return Cube([copy(Q1.base);copy(Q2.base)], [copy(Q1.extent);copy(Q2.extent)])
end

"""
    function primaryFaces(Q::Cube)

Return primary faces (aka with dimension one less) of a given cube.
"""
function primaryFaces(Q::Cube)
    faces = Vector{Cube}()
    for i in 1:length(Q.base)
        if Q.extent[i] != 0
            Q1,Q2 = primaryFaces(Q, i)
            push!(faces, Q1)
            push!(faces, Q2)
        end
    end

    return CubicalSet(faces)
end

function primaryFaces(Q::Cube, i)
    if Q.extent[i] == 0
        return Q
    end
    b1 = copy(Q.base)
    b2 = copy(Q.base)
    e1 = copy(Q.extent)
    e2 = copy(Q.extent)
    b2[i] += 1
    e1[i] = 0
    e2[i] = 0

    return Cube(b1,e1), Cube(b2,e2)
end

"""
    function vertices(Q::Cube)

Return vertices of a given cube.
"""
function vertices(Q::Cube)
    m = dim(Q)
    d = emb(Q)
    v = zeros(Int,d,2^m)

    for i = 0:2^m-1
        s = bitstring(i)[end-m+1:end]
        s = split(s,"")
        s = parse.(Int,s)

        dimcount = 1
        for j = 1:emb(Q)
            if Q.extent[j] == 0
                v[j,i+1] = Q.intervals[j].lower
            else
                if s[dimcount] == 0
                    v[j,i+1] = Q.base[j]
                else
                    v[j,i+1] = Q.base[j]+1
                end
                dimcount += 1
            end
        end
    end

    return v
end

"""
    function boundaryOperator(Q::Cube)

Apply the boundary operator to an cube
"""
function boundaryOperator(Q::Cube)
    # computest the boundary operator of an elementry cube
    sgn = 1
    chain = Dict{Cube,Int}()
    for i = 1:length(Q.base)
        if Q.extent[i] == 1
            Q1,Q2 = primaryFaces(Q, i)
            chain[Q1] = -sgn
            chain[Q2] = sgn
            sgn = -sgn
        end
    end
    return chain
end

# -----------------------------------------------------------------------------
# Start of methods for cubical sets

function CubicalSet(cubes::Vector{Cube})
    if isempty(cubes)
        error("For creating an empty cubical set use CubicalSet(emb::Int) instead")
    end
    embnew = emb(cubes[1])

    for c in cubes
        if emb(c) != embnew
            error("All cubes must be embedded in same dimension")
        end
    end

    cubicalSet = CubicalSet(embnew)
    cubicalSet.cubes = cubes

    return cubicalSet
end

"""
    function emb(x::CubicalSet)
Compute embedding dimension of a cubical set, aka the dim of the ambient space
"""
function emb(X::CubicalSet)
    return X.emb
end

"""
    function addCube(X::CubicalSet, x::Cube)

Adds an cube to a cubical set
"""
function Base.push!(X::CubicalSet, x::Cube)
    if emb(x) == emb(X)
        push!(X.cubes, x)
    else
        error("Dimension mismatch: emb(X) = $(emb(X)), emb(x) = $(emb(x)) " )
    end
    return X
end

function Base.push!(X::CubicalSet, x::Vector{Cube})
    for c in x
        addCube(X,c)
    end
    return X
end

"""
    function product(X::CubicalSet, Y::CubicalSet)

Given two cubical sets returns the product cubical set, i.e. the cubical set consisting of all products Q1xQ2
for Q1 in X, Q2 in Y
"""
function Base.:*(X::CubicalSet, Y::CubicalSet)
    cubes = permutedims(X.cubes) .* Y.cubes

    return CubicalSet(cubes[:])
end

"""
    function sum(X::CubicalSet, Y::CubicalSet)

Given two cubical sets returns the union of the cubical sets. TODO: improve
"""
function Base.:+(X::CubicalSet, Y::CubicalSet)
    cubes = [copy(X.cubes);copy(Y.cubes)]

    return CubicalSet(cubes[:])
end

function createCellComplex(X::CubicalSet)
    # TODO: this can be made much more efficient
    d = emb(X)
    cells=Dict{Int,Set{Cube}}()
    for i = 0:d
        cells[i] = Set{Cube}()
    end

    for Q in X.cubes
        push!(cells[dim(Q)],Q)
    end

    for i = d:-1:1
        for Q in cells[i]
            union!(cells[i-1],primaryFaces(Q).cubes)
        end
    end

    return CellComplex(Dict(map( x -> x[1]=>collect(x[2]),collect(cells))))
end
