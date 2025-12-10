"""
    basic_reduction!(M::AbstractMatrix{T}) where T

Performs the basic persistence matrix reduction of `M` in place.
"""
function basic_reduction!(M::AbstractMatrix{T}) where T
    # NOTE: standard reduction without any speedups
    m,n = size(M)
    V = Matrix(one(T)*I,n,n)
    pivots = Int[]

    for i = m:-1:1
        lowest = findall(!=(0), M[i,:])
        setdiff!(lowest, pivots)
        if length(lowest) == 0
            push!(pivots, 0)
        else
            k = lowest[1]
            c = M[i,k]

            push!(pivots,k)
            for j in lowest[2:end]
                V[:,j] -= M[i,j]/c*V[:,k] # do V first since M[i,j] gets changed!
                M[:,j] -= M[i,j]/c*M[:,k]
            end
        end
    end
    return V
end


function reduced_column_echelon_form!(M)
    basic_reduction!(M)
    pivots = compute_pivots(M)

    # sort by pivot
    sp = sortperm(pivots)
    M .= M[:,sp]
    pivots = pivots[sp]

    # rescale such that pivot is 1
    for (i,t) in enumerate(pivots)
        if t != 0
            M[:,i] ./= M[t,i]
        end
    end

    # fully reduce pivot rows
    _,n = size(M)
    for (i,p) in enumerate(pivots)
        if p != 0
            exhaust_pivot_row!(M, i, p, n)
        end
    end

    return M
end

function exhaust_pivot_row!(M, i, p, n)
    piv_val = M[p,i]
    for j = i+1:n
        if M[p,j] != 0
            M[:,j] -= M[p,j]/piv_val*M[:,i]
        end
    end
end

function compute_pivots(M)
    return map(v -> something(findlast(!iszero, v), 0), eachcol(M))
end

"""
    colspace_normal_form!(M)

Computes and returns the (nonzero part of the) reduced column echelon form of the matrix `M`.

Elementary column operations are performed to transform `M` into a form where:
- Each pivot value is 1 and is the only non-zero entry in its column.
- The leading entry of each non-zero column is below the leading entry of the previous row.
- columns without leading entries are removed from the matrix.
"""
function colspace_normal_form(M)
    N = copy(M)
    reduced_column_echelon_form!(N)
    # remove columns without pivots
    k = findfirst(!iszero, eachcol(N))
    if k === nothing
        return N[:,1:0]  # return empty matrix with correct number of rows
    end
    return N[:,k:end]
end

function is_subspace(V, W)
   M = copy(hcat(W, V))
   basic_reduction!(M)
   return all(iszero, M[:, size(W, 2)+1:end])
end
