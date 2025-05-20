using Random

function lmsTestdata(n,k;bound=100000)
    # generate boundary matrix of 1-skeleton on n vertices
    n_edges = div(n*(n-1),2)
    V = zeros(Int,2,n_edges)
    counter = 1
    for i = 1:n
        for j = 1:(i-1)
            V[:,counter] = [i;j]
            counter += 1
        end
    end
    D0 = spzeros(Int, n_edges,n)
    for i = 1:n_edges
        D0[i,V[1,i]] = 1
        D0[i,V[2,i]] = -1
    end
    # generate boundary matrix for triangles
    T = zeros(Int,3,k)
    d = Dict{Vector{Int},Bool}()

    counter = 0
    counter2 = 0
    f = x -> div((x-1)*(x-2),2)
    while counter < k && counter2<bound
        v = sort(rand(1:n,3))
        if !haskey(d,v) && v[1] != v[2] && v[2] != v[3]
            # v contains edge inidces
            w = zeros(Int,3)
            w[1] = f(v[2])+v[1]
            w[2] = f(v[3])+v[1]
            w[3] = f(v[3])+v[2]
            counter += 1
            T[:,counter] = w
            d[v] = 1
        end
        counter2 += 1
    end

    D1 = spzeros(Int,k,n_edges)
    for i = 1:k
        v = T[:,i]
        D1[i,v[1]] = 1
        if maximum(abs.(D0[v[1],:] + D0[v[2],:])) == 2
            D1[i,v[2]] = -1
        else
            D1[i,v[2]] = 1
        end
        if maximum(abs.(D0[v[1],:] + D1[i,v[2]]*D0[v[2],:] + D0[v[3],:])) ==2
            D1[i,v[3]] = -1
        else
            D1[i,v[3]] = 1
        end
    end
    return D0, D1
end

function lmsTestdataWithRandomTriangles(n,k;bound=100000)
    # generate boundary matrix of 1-skeleton on n vertices
    n_edges = div(n*(n-1),2)
    V = zeros(Int,2,n_edges)
    counter = 1
    for i = 1:n
        for j = 1:(i-1)
            V[:,counter] = [i;j]
            counter += 1
        end
    end
    D0 = spzeros(Int, n_edges,n)
    for i = 1:n_edges
        D0[i,V[1,i]] = 1
        D0[i,V[2,i]] = -1
    end
    # generate boundary matrix for triangles
    T = zeros(Int,3,k)
    vertexT = zeros(Int,3,k)
    d = Dict{Vector{Int},Int}()

    counter = 0
    counter2 = 0
    f = x -> div((x-1)*(x-2),2)
    while counter < k && counter2<bound
        v = sort(rand(1:n,3))
        if !haskey(d,v) && v[1] != v[2] && v[2] != v[3]
            # v contains edge inidces
            w = zeros(Int,3)
            w[1] = f(v[2])+v[1]
            w[2] = f(v[3])+v[1]
            w[3] = f(v[3])+v[2]
            counter += 1
            T[:,counter] = w
            vertexT[:,counter] = v
            d[v] = counter
        end
        counter2 += 1
    end

    D1 = spzeros(Int,k,n_edges)
    for i = 1:k
        v = T[:,i]
        D1[i,v[1]] = 1
        if maximum(abs.(D0[v[1],:] + D0[v[2],:])) == 2
            D1[i,v[2]] = -1
        else
            D1[i,v[2]] = 1
        end
        if maximum(abs.(D0[v[1],:] + D1[i,v[2]]*D0[v[2],:] + D0[v[3],:])) ==2
            D1[i,v[3]] = -1
        else
            D1[i,v[3]] = 1
        end
    end

    # compute list of all possible tetrahedrons
    tetList = Vector{Int}[]
    ind = push!([[2;3;4]], [1;3;4], [1;2;4], [1;2;3])
    @show ind
    for i = 1:k # triangle index
        t = vertexT[:,i]
        for j = 1:n # vertex index
            tet = sort([t;j])

            if all(map(j -> haskey(d, tet[j]) , ind))
                push!(tetList, tet)
            end
        end
    end
    @show length(tetList)
    D2 = spzeros(Int, length(tetList), k)
    for i in 1:length(tetList)
        tet = tetList[i]
        for j = 1:length(ind)
            tri = tet[ind[j]]
            D2[i,d[tri]] = (-1)^j
        end
    end

    return D0, D1, D2
end
#=
D0, D1, D2 = lmsTestdataWithRandomTriangles(40,1650)
firstCohomology(D0,D1)

d = Dict(0 => sparse(D0'), 1=> sparse(D1'), 2=> sparse(D2'), -1 => spzeros(Int, 0,size(D0,2)), 3=> spzeros(Int, size(D2,1),0))

cplx = ChainComplex(d)
redCplx = reduceComplex(cplx)

bm2 = boundaryMatrix(redCplx, 2)
bm1 = boundaryMatrix(redCplx, 1)
smithForm(ZMatrix(Matrix(bm1)))
=#
