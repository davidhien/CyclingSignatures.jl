module ATTools

export boundaryMatrix, coboundaryMatrix, CellComplex, cells, boundaryMatrices, allBoundaryMatrices, deleteAllBoundaryMatrices
export ChainComplex

export Simplex, boundaryOperator
export firstCohomology
export H1, betti_1, circularCoordinates

export vr_incremental

import Base: ==, hash
import SparseArrays: SparseMatrixCSC

using DataStructures
using Distances: chebyshev, pairwise
using Graphs
using LinearAlgebra
using MetaGraphs
using SparseArrays

include("ChainComplex.jl")
include("AbstractCell.jl")
include("Simplicial.jl")
include("H1.jl")
include("Complex.jl")
include("Coordinate.jl")

end
