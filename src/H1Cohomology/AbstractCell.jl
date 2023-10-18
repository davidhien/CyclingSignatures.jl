"""
The following methods must be implemented:
- boundaryOperator
-
"""
abstract type AbstractCell end


"""
    function boundaryMatrix(C1, C2)

Given a pair of cubical chain groups with bdy(C2) ⊂ c1 this computes the boundary matrix
"""
function boundaryMatrix(C1::Vector{T}, C2::Vector{S}) where {S<:AbstractCell, T<:AbstractCell}
    m = length(C1)
    n = length(C2)

    if n == 0
        return spzeros(Int, m, 0)
    end
    if m == 0
        return spzeros(Int, 0, n)
    end

    I = Int[]
    J = Int[]
    V = Int[]

    for j = 1:n
        chain = boundaryOperator(C2[j])
        for cell in sort(collect(keys(chain)))
            i = collect(searchsorted(C1, cell))
            if length(i) != 1
                error("C1 should contain the cell $cell once (is $(length(i)) many times)")
            end
            i = i[1]
            push!(I, i)
            push!(J, j)
            push!(V, chain[cell])
        end
    end
    return sparse(I,J,V,m,n)
end

function coboundaryMatrix(C1::Vector{T}, C2::Vector{S}) where {S<:AbstractCell, T<:AbstractCell}
    m = length(C1)
    n = length(C2)
    if m == 0
        return spzeros(Int, n, 0)
    end
    if n == 0
        return spzeros(Int, 0, m)
    end
    return sparse(boundaryMatrix(C1,C2)')
end

struct CellComplex{T <: AbstractCell} <: AbstractChainComplex
    cells::Dict{Int,Vector{T}}
    boundaryMatrices::Dict{Int,SparseMatrixCSC{Int64,Int64}}

    function CellComplex(cells::Dict{Int,Vector{T}},
        boundaryMatrices::Dict{Int,SparseMatrixCSC{Int64,Int64}}) where T <: AbstractCell
        for k in keys(cells)
            if !issorted(cells[k])
                sort!(cells[k])
            end
        end
        return new{T}(cells, boundaryMatrices)
    end
end

# getter methods
cells(cplx::CellComplex) = getfield(cplx, :cells)
boundaryMatrices(cplx::CellComplex) = getfield(cplx, :boundaryMatrices)

function CellComplex(cells::Dict{Int,Vector{T}}) where T <: AbstractCell
    return CellComplex(cells, Dict{Int,SparseMatrixCSC{Int64,Int64}}())
end

"""
    function allBoundaryMatrices(C::CellComplex)
First, computes all boundary matrices. Then, returns them.
"""
function allBoundaryMatrices(cplx::CellComplex)
    for i in collect(keys(cells(cplx)))
        boundaryMatrix(cplx,i)
    end
    return cplx.boundaryMatrices
end

"""
    function boundaryMatrix(C::CellComplex, i)

Returns the i-th boundary matrix of the complex C
"""
function boundaryMatrix(cplx::CellComplex, i::Integer)
    if haskey(cplx.boundaryMatrices, i)
        return cplx.boundaryMatrices[i]
    elseif haskey(cplx.cells, i) && haskey(cplx.cells, i-1)
        cplx.boundaryMatrices[i] = boundaryMatrix(cplx.cells[i-1], cplx.cells[i])
        return cplx.boundaryMatrices[i]
    elseif haskey(cplx.cells, i)
        return spzeros(Int, 0, length(cplx.cells[i]))
    end
    return spzeros(Int,0,0)
end

"""
    function coboundaryMatrix(C::CellComplex, i)
Returns the i-th coboundary matrix of the complex C
"""
function coboundaryMatrix(cplx::CellComplex, i)
    if haskey(cplx.boundaryMatrices, i+1)
        return sparse(cplx.boundaryMatrices[i+1]')
    elseif haskey(cplx.cells, i+1) && haskey(cplx.cells, i)
        cplx.boundaryMatrices[i+1] = boundaryMatrix(cplx.cells[i], cplx.cells[i+1])
        return sparse(cplx.boundaryMatrices[i+1]')
    elseif haskey(cplx.cells, i)
        return spzeros(Int,0,length(cplx.cells[i]))
    end
    return spzeros(Int,0,0)
end

function deleteAllBoundaryMatrices(cplx::CellComplex)
    bm = boundaryMatrices(cplx)
    for k in keys(bm)
        delete!(bm, k)
    end
end
