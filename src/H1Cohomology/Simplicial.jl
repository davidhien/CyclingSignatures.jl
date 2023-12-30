struct Simplex <: AbstractCell
    vertices::Vector{Int}
    Simplex(vertices::Vector{Int}) = new(sort(vertices))
end

#-----------------
# hash and equals
#-----------------

Base.hash(s::Simplex, h::UInt) = hash(s.vertices, h)

function Base.:(==)(s1::Simplex,s2::Simplex)
    return s1.vertices == s2.vertices
end

Base.isless(s1::Simplex, s2::Simplex) = isless(s1.vertices, s2.vertices)

function dim(Q::Simplex)
    return length(Q.vertices)-1
end

#-----------------
# other methods
#-----------------

function boundaryOperator(Q::Simplex)
    chain = Dict{Simplex,Int}()
    for i = 1:length(Q.vertices)
        i1 = Q.vertices[1:i-1]
        i2 = Q.vertices[i+1:end]
        chain[Simplex([i1;i2])] = (-1)^(i+1)
    end
    return chain
end

function Base.show(io::IO, S::Simplex)
    s = "["
    for i = 1:length(S.vertices)
        if i < length(S.vertices)
            s = string(s, S.vertices[i],",")
        else
            s = string(s, S.vertices[i])
        end

    end
    s = string(s,"]")

    return print(io, s)
end
