mutable struct ZMatrix{T <: Integer}
    # the matrices should satisfy B = Q^-1AR
    B::Array{T,2}
    Q::Array{T,2}
    Qi::Array{T,2}
    R::Array{T,2}
    Ri::Array{T,2}
end

function ZMatrix(B::Array{T,2}) where T <: Integer
    m,n = size(B)
    ZMatrix(copy(B), eye(T, m), eye(T, m), eye(T,n), eye(T, n))
end

function eye(T, k)
    return Matrix{T}(I, k, k)
end

# Helper methods for Z-matrices

function rowExchange!(A::Array{T,2} , i, j) where T <: Integer
    v = A[i, :]
    A[i, :] = A[j, :]
    A[j,:] = v
    return A
end

function rowMultiply!(A::Array{T,2}, i) where T <: Integer
    A[i, : ] = -A[i, :]
    return A
end

function rowAdd!(A::Array{T,2}, i, j, q) where T <: Integer
    # Adds q times row j to row i
    A[i, :] += q*A[j, :]
    return A
end

function columnExchange!(A::Array{T,2}, i, j) where T <: Integer
    v = A[:, i]
    A[:, i] = A[:, j]
    A[:, j] = v
    return A
end

function columnMultiply!(A::Array{T,2}, i) where T <: Integer
    A[:, i] = -A[:, i]
    return A
end

function columnAdd!(A::Array{T,2}, i, j, q) where T <: Integer
    A[:, i] += q*A[:, j]
    return A
end

function rowExchangeOperation(A::ZMatrix, i, j)
    rowExchange!(A.B, i, j)
    rowExchange!(A.Qi, i, j)
    columnExchange!(A.Q, i, j)
    return A
end

function rowMultiplyOperation(A::ZMatrix, i)
    rowMultiply!(A.B, i)
    rowMultiply!(A.Qi, i)
    columnMultiply!(A.Q, i)
    return A
end

function rowAddOperation(A::ZMatrix, i, j, q)
    rowAdd!(A.B, i, j, q)
    rowAdd!(A.Qi, i, j, q)
    columnAdd!(A.Q, j, i, -q)
    return A
end

function columnExchangeOperation(A::ZMatrix, i, j)
    columnExchange!(A.B, i, j)
    columnExchange!(A.R, i, j)
    rowExchange!(A.Ri, i, j)
    return A
end

function columnMultiplyOperation(A::ZMatrix, i)
    columnMultiply!(A.B, i)
    columnMultiply!(A.R, i)
    rowMultiply!(A.Ri, i)
    return A
end

function columnAddOperation(A::ZMatrix, i, j, q)
    columnAdd!(A.B, i, j, q)
    columnAdd!(A.R, i, j, q)
    rowAdd!(A.Ri, j, i, -q)
    return A
end

#reduces the l-th column starting from the k-th row
function partRowReduce(A::ZMatrix{T}, k, l) where T <: Integer
    for i=k+1:size(A.B)[1]
        q = floor(T, A.B[i,l]/A.B[k,l])
        rowAddOperation(A, i, k, -q)
    end
    return A
end

# reduces the k-th row starting from the l-th column
function partColumnReduce(A::ZMatrix{T}, k, l) where T <: Integer
    for i=l+1:size(A.B)[2]
        q = floor(T, A.B[k,i]/A.B[k,l])
        columnAddOperation(A, i, l, -q)
    end
    return A
end

function rowPrepare(A::ZMatrix{T}, k, l) where T <: Integer
    index = k
    value = maximum(abs.(A.B[k:end, l]))
    for i = k:size(A.B)[1]
        if abs(A.B[i,l])<=abs(value) && A.B[i,l]!=0
            index = i
            value = A.B[i,l]
        end
    end
    rowExchangeOperation(A, k, index)
end

function rowReduce(A::ZMatrix{T}, k, l) where T <: Integer
    m = size(A.B)[1]
    while !all(i -> i==0, A.B[k+1:m, l])
        rowPrepare(A, k, l)
        partRowReduce(A, k, l)
    end
    return A
end

function rowEchelon(A::ZMatrix{T}) where T <: Integer
    m = size(A.B)[1]
    n = size(A.B)[2]
    k = 0
    l = 1
    while k < m
        while l <= n && all( i -> i==0, A.B[k+1:m,l])
            l = l+1
        end
        if l == n+1
            break
        end
        k = k+1
        rowReduce(A, k, l)
    end
    return A, k
end
#=
function kernelImage(A::Array{Int,2})
    # returns a basis of the kernel and a basis of the image
    BT = A'
    C, k = rowEchelon(ZMatrix(BT))
    BT = C.B'
    RT = C.Qi'
    return RT[1:end,k+1:end], BT[1:end, 1:k]
end=#

function minNonzero(A::Array{T,2}, r) where T <: Integer
    k, l = 1, 1
    value = maximum(abs.(A[r:end, r:end]))
    for i = r:size(A,1)
        for j=r:size(A,2)
            if abs(A[i,j])<= abs(value) && A[i,j] != 0
                value = A[i,j]
                k,l = i,j
            end
        end
    end
    return value, k, l
end

function moveMinNonzero(A::ZMatrix{T}, k) where T <: Integer
    value, i, j = minNonzero(A.B, k)
    rowExchangeOperation(A, k, i)
    columnExchangeOperation(A, k, j)
    return A
end

function checkForDivisibility(A::Array{T,2}, k) where T <: Integer
    for i=k+1:size(A)[1]
        for j=k+1:size(A)[2]
            q = floor(T, A[i,j]/A[k,k])
            if q*A[k,k] != A[i,j]
                return false, i, j, q
            end
        end
    end
    return true, 0, 0, 0
end

function partSmithForm(A::ZMatrix{T}, k) where T <: Integer
    m = size(A.B)[1]
    n = size(A.B)[2]
    divisible = false
    while !divisible
        moveMinNonzero(A, k)
        partRowReduce(A, k, k)
        if any(i -> i != 0, A.B[k+1:m,k])
            continue
        end
        partColumnReduce(A, k, k)
        if any(i -> i != 0, A.B[k,k+1:n])
            continue
        end
        divisible, i, j, q = checkForDivisibility(A.B, k)
        if !divisible
            rowAddOperation(A, i, k, 1)
            columnAddOperation(A, j, k, -q)
        end
    end
    return A
end

function smithForm(A::ZMatrix{T}) where T <: Integer
    # computest the Smith-Form of A
    # returns a tupel containing
    # 1) A in Smith-Form
    # 2) number of diagonal entries which are 1
    # 3) number of nonzero columns of A.B
    s, t = 0,0
    while any(i-> i!=0, A.B[t+1:end,t+1:end])
        t+=1
        partSmithForm(A,t)
        if A.B[t,t] <0
            rowMultiplyOperation(A, t)
        end
        if A.B[t,t] == 1
            s+=1
        end
    end
    return A, s, t
end

#=
function solve(A::Array{Int, 2}, b::Vector{Int})
    B, s, t = smithForm(ZMatrix(copy(A)))
    c = B.Qi*b
    u = zeros(Int, size(A)[2])
    for i = 1:t
        if c[i] % B.B[t,t] == 0
            u[i] = c[i]/B.B[i,i]
        else
            throw(DomainError())
        end
    end
    if !all(i -> i==0 ,c[t+1:end])
        throw(DomainError())
    end
    return B.R * u
end

function quotientGroup(W::Array{Int,2}, V::Array{Int,2})
    # Computes the quotient W/V
    # Returns a tuple containing
    # 1) A basis of Z^p such that the first s columns a base of V
    # 2) A matrix which on the i,i-th entry conains the order of the group spanned by U[1:end,i]
    # 3) Columns 1, ..., s of U are in 0 + V
    # 4) Columns s+1, ..., s+t have order infinity
    if size(V)[2] == 0
        return W
    end
    A = zeros(Int, size(W)[2], size(V)[2])
    for i = 1:size(V)[2]
        A[:,i] = solve(W, V[:,i])
    end
    B, s, t = smithForm(ZMatrix(A))
    U = W*B.Q
    return U, B.B, s, t
end
=#
