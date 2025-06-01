module CyclingSignatures

# ff.jl
export is_prime, FF

# InclusionMap.jl
export InclusionHelper, SBInclusionHelper

# QuantizedTrajectory.jl
export ResampledTrajectory, tRange, nTimeSteps, getGridPoints
export getInclusionHelper, BoxSpace, SBBoxSpace, getPlotPoints
export TrajectorySpace, trajectoryToTrajectorySpace, trajectoryToTrajectorySpaceSB, getTrajectory, getBoxSpace, getMetric
export maxInclusionThreshold, evaluateCycling, trajectoryBarcode

# lin-alg-util.jl
export basic_reduction!, colspace_normal_form

# SphereBundle.jl
export DynamicDistance

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
using Distances: PreMetric, Metric, pairwise, chebyshev, euclidean
using PersistenceDiagrams
using LinearAlgebra
using StatsBase: countmap
using Graphs, SimpleWeightedGraphs
using ProgressBars
using DataInterpolations
using DataStructures: IntDisjointSets, find_root!

include("H1Cohomology/ATTools.jl")
using .ATTools

include("ff.jl")
include("CyclingSignatures/inclusion-map.jl")
include("CyclingSignatures/trajectory-space.jl")
include("lin-alg-util.jl")
include("dynamic-distance.jl")
include("CyclingSignatures/trajectory.jl")
include("distance-matrix-persistence.jl")
include("CyclingSignatures/subsegment-experiments.jl")
include("CyclingSignatures/interpolate-to-distance.jl")

DEFAULT_FIELD = FF{2}
end
