"""
A chain complex which is made by reducing a different chain complex.

TODO: keep in mind that a reduced complex is not required to have
all boundary matrices of the original complex. In particular, for the
homology implementation, it is important that everything behaves nicely
"""
mutable struct ReducedComplex <: AbstractChainComplex
    complex::ChainComplex
    unredcomplex::AbstractChainComplex
    up_maps
    down_maps
end

function boundaryMatrix(cplx::ReducedComplex, i)
    return boundaryMatrix(cplx.complex, i)
end

function boundaryMatrices(cplx::ReducedComplex)
    return boundaryMatrices(cplx.complex)
end
# TODO: write function boundary_matrix for reduced complex which also handles
# empty boundaries

function cohomologyGeneratorsWithReduction(D0::SparseMatrixCSC{T0,Int},
                                        D1::SparseMatrixCSC{T1,Int}) where {T0 <: Integer, T1 <: Integer}
    # TODO: refactor this and the methods in computeH1_
    bm = Dict(1 =>D1, 2=>D0) #not a typo!
    complex = ChainComplex(bm)
    red_complex = reduceComplex(complex)
    bm_new = boundaryMatrices(red_complex.complex)
    if !haskey(bm_new,2)
        return spzeros(Int,size(D1,2) ,0)
    end
    gen_down = cohomologyGenerators(bm_new[2], bm_new[1])
    if size(gen_down, 2) == 0
        return zeros(size(D1,2), 0)
    end
    gen_up = mapslices(red_complex.up_maps[1],gen_down, dims=[1] )

    return gen_up
end

"""
Given a chain complex, this generates a reduced complex which is chain equivalent
to the original complex as well as the respective chain equivalence.
"""
function reduceComplex(complex::AbstractChainComplex)
    # do the matching
    oldBoundaryMatrices = boundaryMatrices(complex)
    A, w, orders = coreductionMatching(oldBoundaryMatrices)

    # get the gammas
    gammas = getgammas(w,orders,oldBoundaryMatrices)

    # get the up and down maps
    up_maps, down_maps = getupdownmaps(A, oldBoundaryMatrices, gammas)

    # finally compute the boundary matrices
    bm_New = reducedBoundaryMatrices(oldBoundaryMatrices,up_maps,down_maps,A)

    return ReducedComplex(ChainComplex(bm_New), complex, up_maps, down_maps)
end


"""
computes reduced boundary matrices
"""
function reducedBoundaryMatrices(oldBoundaryMatrices,up_maps,down_maps, A)
    boundaryMatrices = Dict{Int,SparseMatrixCSC{Int,Int}}()
    for (k,D) in oldBoundaryMatrices
        if haskey(up_maps,k) && haskey(down_maps,k-1)
            boundaryMatrices[k] = reduceBoundaryMatrix(D,up_maps[k], down_maps[k-1], length(A[k-1]), length(A[k]))
        end
    end

    return boundaryMatrices
end

"""
reduces a single boundary matrix
"""
function reduceBoundaryMatrix(oldBoundaryMatrix,up_map, down_map, m, n)
    boundaryMatrix = spzeros(Int,m,n)

    for i = 1:n
        c = sp_ei(i,n)
        boundaryMatrix[:,i] = down_map(oldBoundaryMatrix*up_map(c))
    end

    return boundaryMatrix
end

"""
Gets all up maps which are needed to define the new boundary matrices.

Given:
- the lists A
- boundary matrices
- gammas

Computes:
- the up_maps from reduced complex to old complex
- the down_maps from old complex to reduced complex
"""
function getupdownmaps(A, boundaryMatrices, gammas)
    up_maps = Dict{Int,Function}()
    down_maps = Dict{Int,Function}()

    i_A = getembedmatrices(A, boundaryMatrices)

    for (k,D) in boundaryMatrices
        if haskey(A,k) && !isempty(A[k])
            up_maps[k] = getupmap(i_A[k],D,gammas[k-1])
        end
        if haskey(A, k-1) && !isempty(A[k-1])
            down_maps[k-1] = getdownmap(i_A[k-1],D,gammas[k-1])
        end
    end

    return up_maps, down_maps
end

"""
given:
- the k-th embedding map
- the k-th boundary map
- the (k-1)-th gamma map

returns:
the map phi_k : C_k -> C^red_k
"""
function getupmap(i_A, D, gamma)
    return c -> i_A*c+gamma(D*i_A*c)
#    return c -> @show size(i_A), length(c), size(D)
end

"""
given:
- the k-th embedding map
- the k+1-th boundary map
- the k-th gamma map

returns:
the map psi_k : C^red_k -> C_k
"""
function getdownmap(i_A, D, gamma)
    return c -> i_A' *(c + D*gamma(c))
end

"""
given A, this generates the matrices
"""
function getembedmatrices(A, boundaryMatrices)
    # TODO: need to make sure that every chain group of A has a boundary map,
    # even if the image of this map is trivial. This is important for
    # correct boundary map computation
    i_A = Dict{Int,SparseMatrixCSC{Int,Int}}()
    for (k,D) in boundaryMatrices
        m,n = size(D)
        i_A[k] = getembedmatrix(A[k], n)
        i_A[k-1] = getembedmatrix(A[k-1],m)
    end

    return i_A
end

"""
given an array A of length n, returns a m-by-n matrix which has
the A[i] th unit vector in the i-th column (1≦i≦n)
"""
function getembedmatrix(A, m)
    n = length(A)
    i_A = spzeros(Int,m,n)

    for i = 1:n
        i_A[A[i],i] = 1
    end
    return i_A
end

"""
generates a vector of all gammas, a function is generated for every
boundary matrix even though it may not be needed
"""
function getgammas(w,orders,boundaryMatrices)
    gammas = Dict{Int,Function}()
    for (k,D) in boundaryMatrices
        gammas[k-1] = getgamma(w[k-1],orders[k-1],D)
    end

    return gammas
end

function getgamma(w,o,D)
    return c -> gammaAlgorithm(c,w,o,D)
end

"""
Given boundary matrices, it returns a dictionary w and a vector o

Requires the topand bottom matrix of the complex to zero/trivial with suitable dimension
"""
function coreductionMatching(boundaryMatrices)
    # TODO: improve the implementation of this method!!!
    # TODO: make sure that empty/top bdy matrix gets dealt with w/o error
    status = Dict{Int,Vector{Bool}}() # true iff cell still in complex
    n_faces = Dict{Int,Vector{Int}}() # number of faces
    queues = Dict{Int,Vector{Int}}() # the queue for each dimension
    n_cells = Dict{Int,Int}()
    total_n_cells = 0

    size_q = 0 # counts

    # declare return  variables
    A = Dict{Int,Vector{Int}}() # contains the A-cells
    w = Dict{Int,Dict{Int,Int}}() # contains the map w
    order = Dict{Int,Vector{Int}}() # contains the order vectors

    # initialze
    min_cell_dim = minimum(keys(boundaryMatrices))-1
    max_cell_dim = maximum(keys(boundaryMatrices))

    for i = min_cell_dim:max_cell_dim
        if haskey(boundaryMatrices, i)
            m,n = size(boundaryMatrices[i])
            n_cells[i] = n
            n_cells[i-1] = m
            n_faces[i] = zeros(Int,n)
            for j = 1:n
                n_faces[i][j] = nnz(boundaryMatrices[i][:,j])
            end
            if !haskey(n_faces,i-1) || isempty(n_faces[i-1])
                n_faces[i-1] = zeros(m)
            end
        else
            n_cells[i] = 0
            n_faces[i] = Int[]
        end
    end

    for i = min_cell_dim:max_cell_dim
        m = n_cells[i]
        status[i] = ones(Bool,m)
        queues[i] = Int[]
        A[i] = Int[]
        w[i] = Dict{Int,Int}()
        order[i] = zeros(m)
        total_n_cells += m
    end

    # main loop
    while total_n_cells > 0
        for n = min_cell_dim:max_cell_dim
            while !isempty(queues[n])
                k = pop!(queues[n])
                if n_faces[n][k] == 1 && status[n][k]
                    q = faces(k, status[n-1], boundaryMatrices[n])[1]

                    removeFromReductionComplex!(n,k,status,n_faces,queues,boundaryMatrices)

                    removeFromReductionComplex!(n-1,q,status,n_faces,queues,boundaryMatrices)

                    # update counter variables
                    size_q += 1
                    n_cells[n] -= 1
                    n_cells[n-1] -= 1
                    total_n_cells -= 2

                    # update return variables
                    w[n-1][q] = k
                    order[n-1][q] = size_q
                    order[n][k] = 0
                end
            end
        end

        # at this point, no coreduction pair should exist anymore
        if total_n_cells > 0
            n = min_cell_dim
            for i = min_cell_dim:max_cell_dim
                if n_cells[i] > 0
                    break
                end
                n += 1
            end
            for i = 1:length(status[n])
                if status[n][i]
                    removeFromReductionComplex!(n, i, status, n_faces,
                        queues, boundaryMatrices)

                    # update counter variables
                    total_n_cells -= 1
                    n_cells[n] -= 1

                    # update return variables
                    order[n][i] = -1
                    push!(A[n], i)
                    break
                end
            end
        end
    end

    return A, w, order
end

"""
given an index k of a cell, computes the indices of the cells which  k
has as face and have their status as true
"""
function faces(k, status, D)
    fcs = Int[]
    ind, val = findnz(D[:,k])
    for i = 1:length(ind)
        if val[i] != 0 && status[ind[i]]
            push!(fcs, ind[i])
        end
    end
    return fcs
end


"""
given an index k of a cell, computes the indices of the cells which have k
as a face and have their status as true
"""
function cofaces(k, status, D)
    cofcs = Int[]
    ind, val = findnz(D[k,:])
    for i = 1:length(ind)
        if val[i] != 0 && status[ind[i]]
            push!(cofcs, ind[i])
        end
    end
    return cofcs
end

"""
n is the index of the chain group
k is the index of the chain
status is the status variable
n_faces is the face count variable
D is the relevant boundary matrix
"""
function removeFromReductionComplex!(n, k, status, n_faces, queues, D)
    status[n][k] = false

    if haskey(D,n+1)
        for f in cofaces(k, status[n+1], D[n+1])
            n_faces[n+1][f] -= 1
            if n_faces[n+1][f] == 1
                push!(queues[n+1], f)
            end
        end
    end
end

"""
Implements the gamma algorithm, where c is a k-chain, w is the
matching function for k to k+1 cells and o is a realization of
the partial order on the cells, D is the k+1 boundatry matrix

the order vector o is supposed to be:
    -1   iff i-th basis element is in A
    x>0  indicating the order else
"""
function gammaAlgorithm(c::AbstractSparseArray, w, o, D)
    cc = copy(c)
    # TODO: Check if chains should be sparse
    m,n = size(D)
    x = spzeros(Int, n)

    # find inclusion maximal
    q, weight = findMaxQinChain(cc,o)

    while weight > 0
        k = w[q]
        omega = -div(cc[q], D[q,k])

        q, weight = findMaxQinChain(cc,o)
        cc = cc + omega*D[:,k]
        x[k]+= omega
    end

    return x
end

"""
given a chain c, finds index which maximizes o and has nnz c
"""
function findMaxQinChain(c::AbstractSparseArray, o)
    index = 0
    weight = -1
    for i in findnz(c)[1]
        if o[i]>weight
            weight = o[i]
            index = i
        end
    end

    return index, weight
end

function sp_ei(k, n)
    return sparse([k], [1], [1], n, 1)
end
