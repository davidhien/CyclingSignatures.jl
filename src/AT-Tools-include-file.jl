import Base: hash, ==
using LinearAlgebra
using SparseArrays
using DataStructures
using Graphs, MetaGraphs
using IterativeSolvers: lsqr, lsqr! # for circular coordinates

include("H1Cohomology/ChainComplex.jl")
include("H1Cohomology/AbstractCell.jl")
include("H1Cohomology/Simplicial.jl")
include("H1Cohomology/H1.jl")
include("H1Cohomology/Complex.jl")
include("H1Cohomology/Coordinate.jl")
include("H1Cohomology/Cubical.jl")
