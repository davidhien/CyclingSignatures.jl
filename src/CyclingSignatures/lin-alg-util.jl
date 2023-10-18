function basicMatrixReduction!(M::AbstractMatrix{T}) where T
    # NOTE: vanilla standard reduction without any speedups
    # TODO: more efficient implementation
    m,n = size(M)
    V = Matrix(one(T)*I,n,n)
    piv = Int[]

    for i = m:-1:1
        lowest = findall(!=(0), M[i,:])
        setdiff!(lowest, piv)
        if length(lowest) == 0
            push!(piv, 0)
        else
            k = lowest[1]
            c = M[i,k]

            push!(piv,k)
            for j in lowest[2:end]
                V[:,j] -= M[i,j]/c*V[:,k] # do V first since M[i,j] gets changed!
                M[:,j] -= M[i,j]/c*M[:,k]
            end
        end
    end
    return V
end

function subspaceNormalForm(M)
    M = copy(M)
    basicMatrixReduction!(M)
    pivots = collect(Iterators.map(v -> findlast(v .!= 0), eachcol(M)))
    
    # filter columns without pivots
    hasPivot = findall(i->!isnothing(i), pivots)
    pivots = pivots[hasPivot]
    M = M[:,hasPivot]

    # sort by pivot
    sp = sortperm(pivots)
    M .= M[:,sp]
    pivots = pivots[sp]

    # fully reduce pivot rows
    _,n = size(M)
    for (i,p) in enumerate(pivots)
        piv_val = M[p,i]
        for j = i+1:n
            if M[p,j] != 0
                M[:,j] -= M[p,j]/piv_val*M[:,i]
            end
        end
    end
    return M
end