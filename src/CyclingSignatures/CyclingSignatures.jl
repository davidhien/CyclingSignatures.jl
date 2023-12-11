module CyclingSignatures

# InclusionMap.jl
export InclusionHelper, SBInclusionHelper

# QuantizedTrajectory.jl
export ResampledTrajectory, tRange, nTimeSteps, getGridPoints
export getInclusionHelper, BoxSpace, SBBoxSpace, getPlotPoints
export TrajectorySpace, trajectoryToTrajectorySpace, trajectoryToTrajectorySpaceSB, getTrajectory, getBoxSpace, getMetric
export maxInclusionThreshold, evaluateCycling, trajectoryBarcode

# LinAlg.jl
export basicMatrixReduction, subspaceNormalForm

# SphereBundle.jl
export SBDistance

# Trajectory.jl
export quantize, resampleToConsistent, resampleToDistance

# SubsegmentExperiments.jl
export SubsegmentSampleParameter,sampleSegments,RandomSubsegmentExperiment,getTrajectorySpace,getSegmentRanges,RandomSubsegmentResult
export runExperiment,runExperiments,getDiagrams
export SubsegmentResultReduced
export rankCountmap, getRankDistribution, rankDistributionMatrix, subspaceDistribution, subspaceFrequencyMatrix
export inclusionVectors, subspaceInclusionDistribution, pairInclusionData, subspaceInclusionMatrix
export getSignatureRanges, isSubspace
export SubsegmentResultSummary


import Base: get, show
import Distances.result_type
using Distances: PreMetric, Metric, pairwise, chebyshev
using PersistenceDiagrams
using LinearAlgebra
using StatsBase: countmap
using ATTools
using Graphs, SimpleWeightedGraphs
using ProgressBars

include("ff.jl")
include("inclusion-map.jl")
include("trajectory-space.jl")
include("lin-alg-util.jl")
include("sphere-bundle.jl")
include("trajectory.jl")
include("distance-matrix-persistence.jl")
include("subsegment-experiments.jl")

DEFAULT_FIELD = FF{2}
 
end
