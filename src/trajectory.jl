"""
Abstract base type for trajectories that can be used in cycling computations.

Subtypes must implement:
- `evaluate_interval(trajectory, a, b)`: Returns points that the trajectory visits between time `a` and `b`.
- `time_domain(trajectory)`: Returns time domain of the trajectory
"""
abstract type AbstractSampleableTrajectory end

"""
    struct RefinedEquidistantTrajectory

Models a refinement of a trajectory on an equidistant grid.

More precisely, the points `[trajectory[:, i] for i in t_vec[1:end-1]]` are assumed to be the sample points, and the other points are the refined points.

It holds for all `i in 1:length(t_vec)-1` that the points in the `i`-th interval are given by the indices in `t_vec[i]:t_vec[i+1]-1`.
"""
struct RefinedEquidistantTrajectory{T,M<:AbstractMatrix{T}} <: AbstractSampleableTrajectory
    trajectory::M                  # columns are points
    t_vec::Vector{Int}             # t_vec[i] = time of i-th point

    function RefinedEquidistantTrajectory(traj::M, t::Vector{Int}) where {T,M<:AbstractMatrix{T}}
        _, n = size(traj)
        if !issorted(t) || last(t) != n + 1
            throw(ArgumentError("Invalid t_vec"))
        end
        return new{T,M}(traj, t)
    end
end

function RefinedEquidistantTrajectory(trajectory)
    _, n = size(trajectory)
    return RefinedEquidistantTrajectory(trajectory, collect(1:n+1))
end

# Interface methods for AbstractSampleableTrajectory
function evaluate_interval(trajectory::RefinedEquidistantTrajectory, a::Int, b::Int)
    traj_range = trajectory.t_vec[a]:trajectory.t_vec[b+1]-1
    return trajectory.trajectory[:, traj_range]
end

function time_domain(trajectory::RefinedEquidistantTrajectory)
    return 1:length(trajectory.t_vec)-1
end

function Base.show(io::IO, t::RefinedEquidistantTrajectory)
    print(io, "RefinedEquidistantTrajectory $(size(t.trajectory,2)) points in ")
    print(io, "$(length(t.t_vec)-1) steps.")
end

function max_consecutive_distance(trajectory::RefinedEquidistantTrajectory, d)
    i1 = eachcol(trajectory.trajectory)
    i2 = Iterators.drop(eachcol(trajectory.trajectory), 1)
    return maximum(zip(i1, i2)) do t
        d(t[1], t[2])
    end
end

function curve_hypothesis_violations(tr::RefinedEquidistantTrajectory, d, r)
    i1 = eachcol(tr.trajectory)
    i2 = Iterators.drop(eachcol(tr.trajectory), 1)
    return count(zip(i1, i2)) do t
        d(t[1], t[2]) >= r
    end
end
