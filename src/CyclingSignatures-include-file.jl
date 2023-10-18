import Base.get
import Distances.result_type
using Distances: PreMetric, Metric, pairwise, chebyshev, euclidean
using PersistenceDiagrams
using LinearAlgebra
using Graphs, SimpleWeightedGraphs
using ProgressBars
using StatsBase: countmap

include("AT-Tools-include-file.jl")
include("CyclingSignatures/ff.jl")
include("CyclingSignatures/inclusion-map.jl")
include("CyclingSignatures/trajectory-space.jl")
include("CyclingSignatures/lin-alg-util.jl")
include("CyclingSignatures/sphere-bundle.jl")
include("CyclingSignatures/trajectory.jl")
include("CyclingSignatures/distance-matrix-persistence.jl")
include("CyclingSignatures/subsegment-experiments.jl")

using GLMakie
include("CyclingSignatures/plotting.jl")

const DEFAULT_FIELD = FF{2}