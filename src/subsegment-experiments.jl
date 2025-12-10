struct RandomSubsegmentExperiment{T,S}
    trajectory_space::T
    segment_lengths::S
    n_runs::Int
    seed::Int
end

function RandomSubsegmentExperiment(trajectory_space, segment_lengths, n_runs)
    return RandomSubsegmentExperiment(trajectory_space, segment_lengths, n_runs, rand(Int))
end

get_trajectory_space(exp::RandomSubsegmentExperiment) = exp.trajectory_space
get_segment_lengths(exp::RandomSubsegmentExperiment) = exp.segment_lengths
get_n_experiments(exp::RandomSubsegmentExperiment) = exp.n_runs

struct RandomSubsegmentResult{T,S,R<:Real}
    trajectory_space::T
    segment_lengths::S
    n_runs::Int
    segment_starts::Vector{Vector{Int}}
    flt_threshold::R
    signatures::Vector{Vector{CyclingSignature}}
end

function Base.show(io::IO, result::RandomSubsegmentResult)
    println(io, "Result($(result.n_runs) runs for $(length(result.segment_lengths)) segment lengths)")
end

function run_experiment(exp::RandomSubsegmentExperiment;
    field=DEFAULT_FIELD,
    alg=Val(:DistanceMatrix),
    threshold=nothing,
    resample_segment_start=true,
    progress=true,
    parallel_inner=false)

    traj_space = get_trajectory_space(exp)
    traj = get_trajectory(traj_space)
    t_indices = time_domain(traj)
    n_time_steps = length(t_indices)
    rng = MersenneTwister(exp.seed)

    seg_lengths = get_segment_lengths(exp)
    n_runs = exp.n_runs

    segment_starts = map(seg_lengths) do len
        sample_segment_starts(n_time_steps, len, n_runs, rng)
    end

    if !resample_segment_start
        _, k = findmax(seg_lengths)
        for i in eachindex(segment_starts)
            if i != k
                segment_starts[i] = segment_starts[k]
            end
        end
    end

    signatures = [Vector{CyclingSignature}(undef, n_runs) for _ in seg_lengths]

    outer_iter = progress ? ProgressBar(enumerate(exp.segment_lengths)) : enumerate(exp.segment_lengths)

    for (i, len) in outer_iter
        if parallel_inner
            @threads for j in 1:n_runs
                cur_range = segment_starts[i][j]:(segment_starts[i][j]+len-1)
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, threshold, field=field)
            end
        else
            inner_iter = 1:n_runs

            for j in inner_iter
                cur_range = segment_starts[i][j]:(segment_starts[i][j]+len-1)
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, threshold, field=field)
            end
        end
    end

    if threshold === nothing
        threshold = traj_space.flt_max_heuristic
    end

    return RandomSubsegmentResult(traj_space, seg_lengths, n_runs, segment_starts, threshold, signatures)
end

function sample_segment_starts(n_time_steps::Int, segment_length::Int, n_experiments::Int, rng=MersenneTwister())
    max_start = n_time_steps - segment_length + 1
    return rand(rng, 1:max_start, n_experiments)
end


##
## Analysis Code
##

"""
    rank_distribution(result::RandomSubsegmentResult, k)

Compute the rank distribution for dimension `k` from the given `RandomSubsegmentResult`.
"""
function rank_distribution(result::RandomSubsegmentResult, k)
    signatures = result.signatures
    return map(signatures) do sigs
        return sum(sigs, init=StepFunction(Int[], [0])) do sig
            p1, p2 = dimension_interval(sig, k)
            indicator_function(p1, p2, 1)
        end
    end
end

function cycspace_length_count(result::RandomSubsegmentResult, k)
    ct_dict = Dict{Matrix{Int},Vector{Int}}()
    for (i,sigs) in enumerate(result.signatures)
        for sig in sigs
            if length(sig.birth_vector) >= k
                V = Int.(colspace_normal_form(sig.cycling_matrix[:, 1:k]))
                if !haskey(ct_dict, V)
                    ct_dict[V] = zeros(Int, length(result.segment_lengths))
                end
                ct_dict[V][i] += 1
            end
        end
    end
    return ct_dict
end

function cycspace_length_countmatrix(result::RandomSubsegmentResult, k; n_subspaces=nothing)
    ct_dict = cycspace_length_count(result::RandomSubsegmentResult, k)
    return lenght_countdict_to_countmatrix(ct_dict; n_relevant=n_subspaces)
end

function cycspace_length_count_at_r(result::RandomSubsegmentResult, k, r)
    ct_dict = Dict{Matrix{Int},Vector{Int}}()
    for (i,sigs) in enumerate(result.signatures)
        for sig in sigs
            V = cycling_matrix(sig, r=r)
            if size(V, 2) == k
                V = Int.(colspace_normal_form(V))
                if !haskey(ct_dict, V)
                    ct_dict[V] = zeros(Int, length(result.segment_lengths))
                end
                ct_dict[V][i] += 1
            end
        end
    end
    return ct_dict
end

function cycspace_length_countmatrix_at_r(result::RandomSubsegmentResult, k, r; n_subspaces=nothing)
    ct_dict = cycspace_length_count_at_r(result, k, r)
    return lenght_countdict_to_countmatrix(ct_dict; n_relevant=n_subspaces)
end

"""
    lenght_countdict_to_countmatrix(ct_dict; n_relevant = nothing)

Convert a length count dictionary (i.e. key maps to a vector of counts per segment length) to a count matrix for the n_relevant most appearing keys.
"""
function lenght_countdict_to_countmatrix(ct_dict; n_relevant = nothing)
    if n_relevant === nothing
        n_relevant = length(collect(keys(ct_dict)))
    end

    ks = collect(keys(ct_dict))
    p = sortperm(map(k-> sum(ct_dict[k]), ks), rev=true)[1:n_relevant]
    relevant_keys = ks[p]

    M = zeros(Int, n_relevant, length(ct_dict[relevant_keys[1]]))
    for (i,k) in enumerate(relevant_keys)
        M[i, :] = ct_dict[k]
    end
    return relevant_keys, M
end

function cycspace_intervals(result::RandomSubsegmentResult, V)
    sig_vecs = result.signatures
    k = size(V, 2)
    V = colspace_normal_form(V)
    int = Vector{Vector{Tuple{Float64,Float64}}}(undef, length(sig_vecs))

    for i in eachindex(sig_vecs)
        sigs = sig_vecs[i]
        int[i] = Tuple{Int,Int}[]
        for sig in sigs
            if length(sig.birth_vector) >= k && colspace_normal_form(sig.cycling_matrix[:, 1:k]) == V
                p1, p2 = dimension_interval(sig, k)
                if p1 <= p2
                    push!(int[i], (p1, p2))
                end
            end
        end
    end
    return int
end

function cycspace_distribution(result::RandomSubsegmentResult, V)
    int_vecs = cycspace_intervals(result, V)

    return map(int_vecs) do vec
        sum(vec, init=indicator_function(1, 0)) do (p1, p2)
            indicator_function(p1, p2, 1.0)
        end
    end
end

function cycspace_radius_distribution(result::RandomSubsegmentResult, V)
    return sum.(cycspace_distribution(result::RandomSubsegmentResult, V))
end

function cycspace_timespan_distribution(result::RandomSubsegmentResult, V)
    return map(cycspace_distribution(result::RandomSubsegmentResult, V)) do f
        return integrate(f)
    end
end

function cycspace_timespan_count(result::RandomSubsegmentResult, V)
    int_vec = cycspace_intervals(result, V)
    return length.(int_vec)
end

"""
    cycspace_segments(result::RandomSubsegmentResult, V)

Get the segments from the result that generate the cycling space V for some radius.
"""
function cycspace_segments(result::RandomSubsegmentResult, V)
    sig_vecs = result.signatures
    seg_starts = result.segment_starts
    k = size(V, 2)
    V = colspace_normal_form(V)

    bool_vecs = map(sig_vecs) do sigs
        map(sigs) do sig
            return length(sig.birth_vector) >= k && colspace_normal_form(sig.cycling_matrix[:, 1:k]) == V
        end
    end

    seg_start_vecs = map(zip(seg_starts, bool_vecs)) do (starts, bools)
        return [starts[i] for i in eachindex(bools) if bools[i]]
    end

    return map(zip(seg_start_vecs, result.segment_lengths)) do (starts, l)
        return [starts[i]:starts[i]+l-1 for i in eachindex(starts)]
    end
end

"""
    cycspace_segments_at_r(result::RandomSubsegmentResult, V, r)

Get the segments from the result have cycling space V at radius r.
"""
function cycspace_segments_at_r(result::RandomSubsegmentResult, V, r)
    sig_vecs = result.signatures
    seg_starts = result.segment_starts
    k = size(V, 2)
    V = colspace_normal_form(V)

    bool_vecs = map(sig_vecs) do sigs
        map(sigs) do sig
            return length(sig.birth_vector) >= k && colspace_normal_form(sig.cycling_matrix[:, 1:k]) == V
        end
    end

    seg_start_vecs = map(zip(seg_starts, bool_vecs)) do (starts, bools)
        return [starts[i] for i in eachindex(bools) if bools[i]]
    end

    return map(zip(seg_start_vecs, result.segment_lengths)) do (starts, l)
        return [starts[i]:starts[i]+l-1 for i in eachindex(starts)]
    end
end

"""
    cycspace_inclusion_matrix(V, W)

For a list of matrices V and W, compute the inclusion matrix where entry (i,j) is 1 if the image of V[i] is included in the image of W[j], and 0 otherwise.
"""
function cycspace_inclusion_matrix(V0, V1)
    n0 = length(V0)
    n1 = length(V1)
    M = zeros(Bool, n0, n1)

    for i in 1:n0
        V_i = colspace_normal_form(V0[i])
        for j in 1:n1
            W_j = colspace_normal_form(V1[j])
            if is_subspace(V_i, W_j)
                M[i, j] = true
            end
        end
    end
    return M
end
