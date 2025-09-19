"""
    struct ResampledTrajectory

ResampledTrajectory models a refinement of a trajectory on an equidistant grid.

More precisely, the points
    [trajectory[:,i] for i in t_vec[1:end-1]]
are assumed to be the sample points and the other points are the refined points.

It holds for ALL i=1:length(t_vec)-1 that the points in the i-th interval are given by the indiced in
    t_vec[i]:t_vec[i+1]-1
"""
struct ResampledTrajectory
    trajectory::Matrix{Float64}
    t_vec::Vector{Int}

    function ResampledTrajectory(trajectory, t_vec::Vector{Int})
        _, n = size(trajectory)
        if !issorted(t_vec) || t_vec[end] != n+1
            throw(ArgumentError("Invalid t_vec"))
        end
        return new(trajectory, t_vec)
    end
end

function ResampledTrajectory(trajectory)
    _, n = size(trajectory)
    return ResampledTrajectory(trajectory, collect(1:n+1))
end

function Base.show(io::IO, t::ResampledTrajectory)
    print(io, "ResampledTrajectory $(size(t.trajectory,2)) points in ")
    print(io, "$(length(t.t_vec)-1) steps.")
end

# TODO: get for single index, range and interval
function Base.get(rt::ResampledTrajectory, t_range::AbstractUnitRange{Int})
    a,b=first(t_range),last(t_range)
    traj_range = rt.t_vec[a]:rt.t_vec[b+1]-1
    return rt.trajectory[:,traj_range]
end

function tRange(rt::ResampledTrajectory)
    return rt.t_vec[[1,end]]
end

function nTimeSteps(rt::ResampledTrajectory)
    return length(rt.t_vec)-1
end

function getGridPoints(rt::ResampledTrajectory, t_range::AbstractUnitRange{Int})
    grid_t = rt.t_vec[t_range]
    return rt.trajectory[:,grid_t]
end

function curveHypothesis(rt::ResampledTrajectory, d)
    i1 = eachcol(rt.trajectory)
    i2 = Iterators.drop(eachcol(rt.trajectory),1)
    return maximum(zip(i1,i2)) do t
        d(t[1],t[2])
    end
end

function countCurveHypothesisViolations(rt::ResampledTrajectory, d, r)
    i1 = eachcol(rt.trajectory)
    i2 = Iterators.drop(eachcol(rt.trajectory),1)
    return count(zip(i1,i2)) do t
        d(t[1],t[2])>= r
    end
end

abstract type AbstractBoxSpace end

getInclusionHelper(bs::T) where {T<:AbstractBoxSpace} = getfield(bs, :inc_helper)

struct BoxSpace{S<:AbstractVector{Int},T<:Real} <: AbstractBoxSpace
    X::Matrix{Int}
    boxsize::T
    h1::H1
    inc_helper::InclusionHelper{S,T}
end

function BoxSpace(X, boxsize; kwargs...)
    X = unique(X,dims=2)
    dm_X = pairwise(chebyshev, X)
    cplx, h1_gen, h1 = boxSpaceH1Info(dm_X)
    inc_helper = InclusionHelper(cplx, X, h1_gen, boxsize; kwargs...)
    return BoxSpace(dm_X, boxsize, h1, inc_helper)
end

function Base.show(io::IO, bs::BoxSpace)
    print(io, "BoxSpace: ")
    print(io, "$(size(bs.X,2)) boxes of size $(bs.boxsize), ")
    print(io, "with β1 = $(betti_1(bs.h1))")
end

struct SBBoxSpace{S<:AbstractVector{Int},T<:Real} <: AbstractBoxSpace
    X::Matrix{Int}          # quantization in space
    VF::Matrix{Int}         # quantization in bundle
    SB_ind::Matrix{Int}     # points in sphere bundle are [ [X[:,i];VF[:,j]] for (i,j) in SB_ind ]
    boxsize::T               # boxsize for space
    sphereBundleRadius::Int # radius of sphere for sphere bundle
    h1::H1                  # h1 object for coordinates
    inc_helper::SBInclusionHelper{S,T}
end

function SBBoxSpace(SB_pts, boxsize, sphereBundleRadius; kwargs...)
    SB_pts = unique(SB_pts, dims=2)
    dm_pts = pairwise(chebyshev, SB_pts, dims=2)
    cplx, h1_gen, h1 = boxSpaceH1Info(dm_pts)

    X,VF,SB_ind = sphereBundleToComponents(SB_pts)

    inc_helper = SBInclusionHelper(cplx, SB_pts, h1_gen, boxsize, sphereBundleRadius; kwargs...)

    return SBBoxSpace(X,VF,SB_ind,boxsize, sphereBundleRadius,h1, inc_helper)
end

function Base.show(io::IO, bs::SBBoxSpace)
    print(io, "SBBoxSpace: ")
    print(io, "$(size(bs.SB_ind,2)) boxes of size ($(bs.boxsize),1/$(bs.sphereBundleRadius)) ")
    print(io, "with β1 = $(betti_1(bs.h1))")
end

function getPlotPoints(bs::SBBoxSpace,offset=.1)
    X = bs.X
    VF = bs.VF
    return reduce(hcat,map(eachcol(bs.SB_ind)) do v
        return X[:,v[1]] + offset*VF[:,v[2]]
    end)
end

function boxSpaceH1Info(dm_X)
    cplx = vr_incremental(dm_X, 1.5)
    @info string("#edges=", length(cplx.cells[1]), ", #triangles=", length(cplx.cells[2]))
    D0 = coboundaryMatrix(cplx, 0)
    D1 = coboundaryMatrix(cplx, 1)

    @time h1_gen = firstCohomology(D0,D1)
    h1 = H1(h1_gen,D0,D1)
    return cplx, h1_gen, h1
end

"""
    function sphereBundleToComponents(SB_pts)

Returns X, VF, ind such that
SB_pts[:,i] = [X[:,ind[1,i]]; VF[:,ind[2,i]]]
"""
function sphereBundleToComponents(SB_pts)
    m = size(SB_pts,1)
    d = div(m,2)
    Xt = SB_pts[1:d,:]
    VFt = SB_pts[d+1:2*d,:]
    tX, X = sortedDynamicIndices(Xt)
    tVF, VF = sortedDynamicIndices(VFt)

    return X, VF, [tX';tVF']
end

struct TrajectorySpace{T<: AbstractBoxSpace,S<:PreMetric}
    trajectory::ResampledTrajectory
    boxSpace::T
    metric::S
end

function trajectoryToTrajectorySpace(Y, boxsize; t_vec=nothing, metric=euclidean, kwargs...)
    rt = nothing
    if isnothing(t_vec)
        rt = ResampledTrajectory(Y)
    else
        rt = ResampledTrajectory(Y, t_vec)
    end

    X = quantize(Y,boxsize)
    boxSpace = BoxSpace(X, boxsize; kwargs...)

    return TrajectorySpace(rt, boxSpace, metric)
end

function trajectoryToTrajectorySpaceSB(Y,Z, boxsize, sb_radius; t_vec=nothing, metric=nothing, kwargs...)
    for c in eachcol(Z)
        if !isapprox(norm(c),1)
            error("Sphere Bundle component is not an unit vectors.")
        end
    end
    rt = nothing
    ts = [Y;Z]
    if isnothing(t_vec)
        rt = ResampledTrajectory(ts)
    else
        rt = ResampledTrajectory(ts, t_vec)
    end

    SB_pts = unique([quantize(Y, boxsize);quantizeSB(Z,sb_radius)],dims=2)
    boxSpace = SBBoxSpace(SB_pts, boxsize, sb_radius; kwargs...)

    if isnothing(metric)
        d = div(size(ts,1),2)
        # nerve theorem allows a blow up to boxsize in space and 1/sb_radius in the sphere
        metric = DynamicDistance(d, boxsize*sb_radius)
    end

    return TrajectorySpace(rt, boxSpace, metric)
end

getTrajectory(ts::TrajectorySpace) = ts.trajectory
getBoxSpace(ts::TrajectorySpace) = ts.boxSpace
getMetric(ts::TrajectorySpace) = ts.metric

function Base.show(io::IO, ts::TrajectorySpace)
    print(io, "TrajectorySpace with\n ∘ ")
    show(io, getTrajectory(ts))
    print(io, "\n ∘ ")
    show(io, getBoxSpace(ts))
    print(io, "\n ∘ ")
    show(io, getMetric(ts))
end

function curveHypothesis(ts::TrajectorySpace)
    return curveHypothesis(getTrajectory(ts), getMetric(ts))
end

function countCurveHypothesisViolations(ts::TrajectorySpace, r)
    return countCurveHypothesisViolations(getTrajectory(ts), getMetric(ts), r)
end

"""
function maxInclusionThreshold(boxsize, sb_radius, C)

Computes the maximal threshold such that a general inclusion is well-defined.
It is assumed that the metric is of the form
d((p,v),(q,w)) = d(p,q) + C*d(v,w)
"""
function maxInclusionThreshold(boxsize, sb_radius, C)
    return min(boxsize, C/sb_radius)
end

"""
    function evaluateCycling(alg::Val, trajSpace::TrajectorySpace, range, fltThreshold)

Computes the cycling barcode for a segment of trajSpace.

The range indicates the segment of trajSpace to be used, fltThreshold the maximal threshold for the persistence calculation.
"""
function evaluateCycling(alg::Val, trajSpace::TrajectorySpace, range, fltThreshold, field=DEFAULT_FIELD)
    traj = getTrajectory(trajSpace)
    trajPoints = get(traj, range)
    trajBars = trajectory_barcode(alg, trajPoints, getMetric(trajSpace), fltThreshold, field)
    return inclusionDiagram(trajPoints, trajBars, getBoxSpace(trajSpace))
end

"""
    function trajectory_barcode(alg, trajPoints, metric, fltThreshold)

Computes a barcode for the trajectory specified by `trajPoints` using the given `metric`.
The parameter `fltThreshold` determines how far the filtrations are evaluated.

# Arguments
- `alg`: Specifies the algorithm to use. Currently supports `:Ripserer` or `:DistanceMatrix`.
- `trajPoints`: A matrix where each column represents a point on the trajectory.
- `metric`: A function of two arguments that satisfies the usual metric axioms.
- `fltThreshold`: A number indicating the filtration threshold.

# Note
The returned persistence diagram is not guaranteed to be a true persistence diagram for the given curve.
"""
trajectory_barcode(alg::Any, trajPoints, metric, fltThreshold, field = DEFAULT_FIELD) = error("Not implemented for specified alg.")

"""
    function inclusionDiagram(trajPoints, trajDiag, boxSpace)

trajDiag is expected to be a persistence diagram where each interval has meta with: simplex_list, coeff_list
"""
function inclusionDiagram(trajPoints, trajDiag, boxSpace)
    infinite_intervals = filter(x -> x.death==Inf, trajDiag)
    # TODO: potential speed-up: dont use PersistenceInterval as carrier for representatives
    included_intervals = map(infinite_intervals) do int
        inc_rep = includeCycle(getInclusionHelper(boxSpace), trajPoints, int.simplex_list, int.coeff_list)
#        inc_rep = includeCycle(getInclusionHelper(boxSpace), trajPoints, int.simplex_list::Vector{Tuple{Int,Int}}, int.coeff_list::Vector{<:Integer})
        sig_rep = inclusionMatrix(boxSpace.inc_helper) * inc_rep
        return PersistenceInterval(int.birth,int.death, inclusion_representative=sig_rep) #, simplex_list=int.simplex_list, coeff_list=int.coeff_list)
    end
    sort!(included_intervals, by=birth)
    independent_intervals = removeDependentIntervals(included_intervals)
    return PersistenceDiagram(independent_intervals)
end

"""
    function removeDependentIntervals(int)

int is required to have the keyword inclusion_generator
"""
function removeDependentIntervals(int)
    int = filter(v -> any(!=(0), v.inclusion_representative) , int)
    if isempty(int)
        return int
    end
    int = filter(v -> any(!=(0), v.inclusion_representative) , int)
    M = reduce(hcat, map(v-> v.inclusion_representative, int))
    basic_reduction!(M)
    keep_indices = findall(v-> any(!=(0), v), eachcol(M))::Vector{Int}
    return int[keep_indices]
end

function getTrajDiag(alg::Val, trajSpace::TrajectorySpace, range, fltThreshold)
    traj = getTrajectory(trajSpace)
    trajPoints = get(traj, range)
    return trajectory_barcode(alg, trajPoints, getMetric(trajSpace), fltThreshold)
end

# TODOs:
# - allow to specify metric in BoxSpace and SBBoxSpace
# - include a metric for the filtration
# - construct VR complex of BoxSpace directly from metric or even better, compute h_1 directly from
# - improve API for TrajectorySpace
# - remove PersistenceDiagrams dependency, or at least dont use keyword-argument fields to store representative data
