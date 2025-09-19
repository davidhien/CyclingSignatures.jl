module CyclingSignatures

# ff.jl
export is_prime, FF

# trajectory.jl
export AbstractSampleableTrajectory, RefinedEquidistantTrajectory, time_domain,
    evaluate_interval, max_consecutive_distance, curve_hypothesis_violations

# comparison-space.jl
export AbstractCubicalAcyclicCarrier, induced_one_chain, annotate_chain, betti_1
export AbstractComparisonSpace, map_cycle, CubicalComparisonSpace, cubical_vr_comparison_space_via_cover, sb_cubical_vr_comparison_space_via_cover

# cycling-signature.jl
export CyclingSignature, dimension, cycling_matrix, birth_vector
export TrajectorySpaceNew, trajectory, comparison_space, metric
export cycling_signature

# trajectory-space.jl
export InclusionHelper, SBInclusionHelper

# QuantizedTrajectory.jl
export ResampledTrajectory, tRange, nTimeSteps, getGridPoints
export getInclusionHelper, BoxSpace, SBBoxSpace, getPlotPoints
export TrajectorySpace, trajectoryToTrajectorySpace, trajectoryToTrajectorySpaceSB, getTrajectory, getBoxSpace, getMetric
export maxInclusionThreshold, evaluateCycling, trajectory_barcode

# lin-alg-util.jl
export basic_reduction!, colspace_normal_form

# dynamic-distance.jl
export DynamicDistance

# ???
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
using SparseArrays: spzeros
using StatsBase: countmap
using Graphs, SimpleWeightedGraphs
using ProgressBars
using DataInterpolations
using DataStructures: IntDisjointSets, find_root!

include("H1Cohomology/ATTools.jl")
using .ATTools

include("ff.jl")
include("comparison-space.jl")
include("CyclingSignatures/inclusion-map-old.jl")
include("CyclingSignatures/trajectory-space.jl")
include("lin-alg-util.jl")
include("dynamic-distance.jl")
include("trajectory.jl")
include("cycling-signature.jl")
include("CyclingSignatures/sample-tools.jl")
include("distance-matrix-persistence.jl")
include("CyclingSignatures/subsegment-experiments.jl")
include("CyclingSignatures/interpolate-to-distance.jl")

DEFAULT_FIELD = FF{2}
end
