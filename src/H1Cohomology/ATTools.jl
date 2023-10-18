module ATTools

# TODO:
# - export for cubical
# - remove the old weird H1 type, more generally, look over Coordinate.jl
# - create proper chain complex etc API

export boundaryMatrix, coboundaryMatrix, CellComplex, cells, boundaryMatrices, allBoundaryMatrices, deleteAllBoundaryMatrices
export ChainComplex
# export for cubical here
export Simplex, boundaryOperator
export cohomologyGenerators
export reduceComplex
export firstCohomology
export H1, betti_1, circularCoordinates

import Base: hash, ==
using SmithNormalForm
using LinearAlgebra
using SparseArrays
using DataStructures
using LightGraphs, MetaGraphs

include("ChainComplex.jl")
include("AbstractCell.jl")
include("Cubical.jl")
include("Simplicial.jl")
include("Homology.jl")
include("Reduction.jl")
include("H1.jl")
include("Complex.jl")
include("Coordinate.jl")

end
