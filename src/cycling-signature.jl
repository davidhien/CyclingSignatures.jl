struct CyclingSignature{F,M<:AbstractMatrix{F},T}
    cycling_matrix::M
    birth_vector::Vector{T}

    function CyclingSignature(cycling_matrix::M, birth_vector::Vector{T}) where {F,M<:AbstractMatrix{F},T}
        if size(cycling_matrix, 2) == length(birth_vector)
            new{F,M,T}(cycling_matrix, birth_vector)
        else
            throw(ArgumentError("#cols cycling matrix and length of birth_vector must be equal. Got: $(size(cycling_matrix, 2)) and $(length(birth_vector))"))
        end
    end
end

function dimension(cs::CyclingSignature, r=Inf)
    return count(<=(r), cs.birth_vector)
end

function cycling_matrix(cs::CyclingSignature; r=Inf)
    idx = findall(<=(r), cs.birth_vector)
    return cs.cycling_matrix[:, idx]
end

function birth_vector(cs::CyclingSignature; r=Inf)
    return filter(<=(r), cs.birth_vector)
end

struct TrajectorySpaceNew{T<:AbstractSampleableTrajectory,S<:AbstractComparisonSpace,D<:PreMetric}
    trajectory::T
    comparison_space::S
    metric::D
end

trajectory(trajectory_space::TrajectorySpaceNew) = trajectory_space.trajectory
comparison_space(trajectory_space::TrajectorySpaceNew) = trajectory_space.comparison_space
metric(trajectory_space::TrajectorySpaceNew) = trajectory_space.metric

function cycling_signature(trajectory_space::TrajectorySpaceNew, range, r_max, field=DEFAULT_FIELD)
    return cycling_signature(Val(:DistanceMatrix), trajectory_space, range, r_max, field)
end

"""
    cycling_signature([alg::Val, ]trajectory_space::TrajectorySpaceNew, range, r_max, field=DEFAULT_FIELD)

Computes the cycling signature of the segment specified by `range` inside of `trajectory_space` for filtration values `` [0,r_{max}] ``.
Optionally, a field and an algorithm can be specified, see REF:TODO.

# Returns
TODO
"""
function cycling_signature(alg::Val, trajectory_space::TrajectorySpaceNew, range, r_max, field=DEFAULT_FIELD)
    # compute trajectory points and barcode
    traj = trajectory(trajectory_space)
    traj_pts = evaluate_interval(traj, range[1], range[2])
    traj_barc = trajectory_barcode(alg, traj_pts, metric(trajectory_space), r_max, field)

    # compute cycling signature in the homological comparison space
    births, cyc_vectors = births_cycling_vectors_from_trajectory_barcode(traj_pts, traj_barc, comparison_space(trajectory_space))
    length(births) != length(cyc_vectors) && @warn "Inconsistent lengths of births and cycling vectors"
    if isempty(births)
        return CyclingSignature(zeros(field,0,0), Float64[])
    end
    cycling_matrix = hcat(cyc_vectors...)

    return CyclingSignature(cycling_matrix, births)
end

function births_cycling_vectors_from_trajectory_barcode(points, barcode, comparison_space)
    infinite_bars = filter(x -> x.death == Inf, barcode)
    included_bars = map(infinite_bars) do bar
        sig_rep = map_cycle(comparison_space, points, bar.simplex_list, bar.coeff_list)
        return (bar.birth, sig_rep)
    end
    sort!(included_bars, by=x -> x[1])

    birth_vector = map(t -> t[1], included_bars)
    cycling_vectors = map(t -> t[2], included_bars)

    birth_vector, cycling_vectors = remove_dependent_bars(birth_vector, cycling_vectors)

    return birth_vector, cycling_vectors
end

function remove_dependent_bars(birth_vector, cycling_vectors)
    # bars is assumed to be a vector of pairs (birth, signature)
    nz_vec = findall(v -> any(!=(0), v), cycling_vectors)
    if isempty(nz_vec)
        return birth_vector[1:0], birth_vector[1:0]
    end

    cycling_vectors = cycling_vectors[nz_vec]
    birth_vector = birth_vector[nz_vec]

    M = reduce(hcat, cycling_vectors)
    basic_reduction!(M)
    keep_indices = findall(v -> any(!=(0), v), eachcol(M))::Vector{Int}
    return birth_vector[keep_indices], cycling_vectors[keep_indices]
end
