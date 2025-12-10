module CyclingSignatures

# ff.jl
export is_prime, FF

# trajectory.jl
export AbstractSampleableTrajectory, RefinedEquidistantTrajectory, time_domain,
    evaluate_interval, t_vec_segment, max_consecutive_distance, curve_hypothesis_violations

# comparison-space.jl
export AbstractCubicalAcyclicCarrier, induced_one_chain, annotate_chain, betti_1
export AbstractComparisonSpace, map_cycle, CubicalComparisonSpace, cubical_vr_comparison_space_via_cover, sb_cubical_vr_comparison_space_via_cover

# cycling-signature.jl
export CyclingSignature, dimension, cycling_matrix, birth_vector
export TrajectorySpace, get_trajectory, get_comparison_space, get_metric
export trajectory_space_from_trajectory, utb_trajectory_space_from_trajectory
export cycling_signature

# lin-alg-util.jl
export basic_reduction!, colspace_normal_form

# dynamic-distance.jl
export DynamicDistance

# subsegment-experiments.jl
export RandomSubsegmentExperiment, get_trajectory_space, get_segment_lengths, get_n_experiments
export RandomSubsegmentResult, run_experiment, sample_segment_starts
export rank_distribution, cycspace_intervals, cycspace_distribution
export cycspace_length_count, cycspace_length_countmatrix, cycspace_length_count_at_r, cycspace_length_countmatrix_at_r
export cycspace_segments, cycspace_segments_at_r
export cycspace_inclusion_matrix

import Base: get, show
import Distances.result_type
using Distances: PreMetric, Metric, pairwise, chebyshev, euclidean
using DataStructures
using LinearAlgebra
using SparseArrays: spzeros
using StatsBase: countmap
using Graphs, SimpleWeightedGraphs
using ProgressBars
using DataInterpolations
using StepFunctions
using Random: MersenneTwister
using Base.Threads

include("H1Cohomology/ATTools.jl")
using .ATTools

include("ff.jl")
include("comparison-space.jl")
include("lin-alg-util.jl")
include("dynamic-distance.jl")
include("trajectory.jl")
include("cycling-signature.jl")
include("CyclingSignatures/sample-tools.jl")
include("distance-matrix-persistence.jl")
include("CyclingSignatures/interpolate-to-distance.jl")
include("subsegment-experiments.jl")

include("plotting-interface.jl")
export plot_rank_distribution,
       plot_rank_heatmap,
       plot_all_rank_heatmaps,
       plot_rank_distribution_at_r,
       plot_subspace_frequency_at_r,
       plot_cycspace_inclusion

include("dm-persistence-birth-curves.jl")
export birth_curves, CyclingBirthCurve

DEFAULT_FIELD = FF{2}
end
