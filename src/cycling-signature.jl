struct CyclingSignature{F,M<:AbstractMatrix{F},T}
    cycling_matrix::M
    birth_vector::Vector{T}

    function CyclingSignature(cycling_matrix::M, birth_vector::Vector{T}) where {F,M<:AbstractMatrix{F},T}
        if size(cycling_matrix, 2) == length(birth_vector)
            return new{F,M,T}(cycling_matrix, birth_vector)
        else
            throw(ArgumentError("#cols cycling matrix and length of birth_vector must be equal. Got: $(size(cycling_matrix, 2)) and $(length(birth_vector))"))
        end
    end
end

function dimension(cs::CyclingSignature, r=Inf)
    return count(<=(r), cs.birth_vector)
end

function dimension_interval(cs::CyclingSignature, k)
    if k < 0 || k > length(cs.birth_vector)
        return (0, -1)
    end
    if k == 0 && length(cs.birth_vector) == 0
        return (0, Inf)
    elseif k == 0
        return (0, cs.birth_vector[1])
    elseif k == length(cs.birth_vector)
        return (cs.birth_vector[end], Inf)
    else
        return (cs.birth_vector[k], cs.birth_vector[k+1])
    end
end

function dimension_function(cs::CyclingSignature)
    return StepFunction(cs.birth_vector, 0:length(birth_vector))
end

function cycling_matrix(cs::CyclingSignature; r=Inf)
    idx = findall(<=(r), cs.birth_vector)
    return cs.cycling_matrix[:, idx]
end

function birth_vector(cs::CyclingSignature; r=Inf)
    return filter(<=(r), cs.birth_vector)
end

struct TrajectorySpace{T<:AbstractSampleableTrajectory,S<:AbstractComparisonSpace,D<:PreMetric,F}
    trajectory::T
    comparison_space::S
    metric::D
    flt_max_heuristic::F
end

function TrajectorySpace(trajectory, comparison_space, metric)
    return TrajectorySpace(trajectory, comparison_space, metric, nothing)
end

function trajectory_space_from_trajectory(x, boxsize; t_vec=collect(1:size(x, 2)+1), metric = euclidean, flt_max_heuristic=nothing)
    traj = RefinedEquidistantTrajectory(x, t_vec)
    comparison_space = cubical_vr_comparison_space_via_cover(x, boxsize)
    if flt_max_heuristic === nothing
        flt_max_heuristic = boxsize
    end
    return TrajectorySpace(traj, comparison_space, metric, flt_max_heuristic)
end

function utb_trajectory_space_from_trajectory(x, tx, boxsize, sb_radius; t_vec=collect(1:size(x, 2)+1), metric=nothing, flt_max_heuristic=nothing)
    for (i, t) in enumerate(eachcol(tx))
        if !isapprox(norm(t), 1)
            error("Sphere Bundle component $i is not a unit vector.")
        end
    end
    M = [x; tx]
    traj = RefinedEquidistantTrajectory(M, t_vec)
    comparison_space = sb_cubical_vr_comparison_space_via_cover(M, boxsize, sb_radius)

    if metric === nothing
        d = size(x, 1)
        # nerve theorem allows a blow up to boxsize in space and 1/sb_radius in the sphere
        metric = DynamicDistance(d, boxsize * sb_radius)
    end
    if flt_max_heuristic === nothing
        flt_max_heuristic = boxsize
    end
    return TrajectorySpace(traj, comparison_space, metric, flt_max_heuristic)
end

get_trajectory(trajectory_space::TrajectorySpace) = trajectory_space.trajectory
get_comparison_space(trajectory_space::TrajectorySpace) = trajectory_space.comparison_space
get_metric(trajectory_space::TrajectorySpace) = trajectory_space.metric
betti_1(trajectory_space::TrajectorySpace) = betti_1(get_comparison_space(trajectory_space))

function cycling_signature(trajectory_space::TrajectorySpace, range, r_max=nothing; field=DEFAULT_FIELD)
    return cycling_signature(Val(:DistanceMatrix), trajectory_space, range, r_max, field=field)
end

"""
    cycling_signature([alg::Val, ]trajectory_space::TrajectorySpace, range; r_max, field=DEFAULT_FIELD)

Computes the cycling signature of the segment specified by `range` inside of `trajectory_space` for filtration values `` [0,r_{max}] ``.
Optionally, a field and an algorithm can be specified.

# Arguments
- `alg`: currently only :DistanceMatrix works
- `trajectory_space`: the trajectory space
- `range`: currently, this evaluates the cycling signature on the interval `[first(range):last(range)]`
- `r_max`: if unspecified, the default in `trajectory_space` is used
- `field`: coefficient field for homology computation

# Returns
TODO
"""
function cycling_signature(alg::Val, trajectory_space::TrajectorySpace, range, r_max=nothing; field=DEFAULT_FIELD)
    if r_max === nothing
        if trajectory_space.flt_max_heuristic === nothing
            throw(ArgumentError("r_max must be specified, if trajectory_space.flt_max_heuristic is nothing."))
        end
        r_max = trajectory_space.flt_max_heuristic
    end
    # compute trajectory points and barcode
    traj = get_trajectory(trajectory_space)
    traj_pts = evaluate_interval(traj, first(range), last(range))
    traj_barc = trajectory_barcode(alg, traj_pts, get_metric(trajectory_space), r_max, field)

    # compute cycling signature in the homological comparison space
    births, cyc_vectors = births_cycling_vectors_from_trajectory_barcode(traj_pts, traj_barc, get_comparison_space(trajectory_space))
    length(births) != length(cyc_vectors) && @warn "Inconsistent lengths of births and cycling vectors"
    if isempty(births)
        return CyclingSignature(zeros(field, 0, 0), Float64[])
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
