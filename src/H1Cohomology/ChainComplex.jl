#
#
# Implementations of AbstractChainComplex can support a lazy computation of
# boundary matrices,
#
abstract type AbstractChainComplex end

boundaryMatrices(cplx::AbstractChainComplex) = error()

mutable struct ChainComplex <: AbstractChainComplex
    boundaryMatrices::Dict{Int,SparseMatrixCSC{Int64,Int64}}
end

#getter methodss
boundaryMatrices(cplx::ChainComplex) = getfield(cplx, :boundaryMatrices)

function boundaryMatrix(cplx::ChainComplex, i)
    return boundaryMatrices(cplx)[i]
end
