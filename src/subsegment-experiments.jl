struct RandomSubsegmentExperiment{T,S}
    trajectory_space::T
    segment_lengths::S
    n_runs::Int
    seed::Int
end

function RandomSubsegmentExperiment(trajectory_space, segment_lengths, n_runs)
    return RandomSubsegmentExperiment(trajectory_space, segment_lengths, n_runs, rand(1:typemax(Int)))
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
    parallel_inner=false,
    segment_starts=nothing)

    traj_space = get_trajectory_space(exp)
    seg_lengths = get_segment_lengths(exp)
    n_runs = exp.n_runs
    resolved_threshold = _resolve_r_max(traj_space, threshold)

    if segment_starts === nothing
        segment_starts = sample_segment_starts(exp; resample_segment_start)
    else
        segment_starts = _validate_segment_starts(traj_space, seg_lengths, n_runs, segment_starts)
    end

    signatures = [Vector{CyclingSignature}(undef, n_runs) for _ in seg_lengths]

    outer_iter = progress ? ProgressBar(enumerate(exp.segment_lengths)) : enumerate(exp.segment_lengths)

    for (i, len) in outer_iter
        if parallel_inner
            @threads for j in 1:n_runs
                cur_range = segment_starts[i][j]:(segment_starts[i][j]+len-1)
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, resolved_threshold, field=field)
            end
        else
            inner_iter = 1:n_runs

            for j in inner_iter
                cur_range = segment_starts[i][j]:(segment_starts[i][j]+len-1)
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, resolved_threshold, field=field)
            end
        end
    end

    return RandomSubsegmentResult(traj_space, seg_lengths, n_runs, segment_starts, resolved_threshold, signatures)
end

function sample_segment_starts(n_time_steps::Int, segment_length::Int, n_experiments::Int, rng=MersenneTwister())
    max_start = n_time_steps - segment_length + 1
    max_start >= 1 || throw(ArgumentError("segment_length must not exceed the number of time steps. Got $segment_length and $n_time_steps."))
    return rand(rng, 1:max_start, n_experiments)
end

function sample_segment_starts(
    trajectory_space::TrajectorySpace,
    segment_lengths,
    n_runs::Integer,
    rng::AbstractRNG=MersenneTwister();
    resample_segment_start::Bool=true,
)
    n_runs >= 0 || throw(ArgumentError("n_runs must be non-negative. Got $n_runs."))
    n_time_steps = length(time_domain(get_trajectory(trajectory_space)))
    starts = map(segment_lengths) do len
        sample_segment_starts(n_time_steps, len, n_runs, rng)
    end

    if !resample_segment_start && !isempty(starts)
        _, k = findmax(segment_lengths)
        for i in eachindex(starts)
            if i != k
                starts[i] = copy(starts[k])
            end
        end
    end

    return _validate_segment_starts(trajectory_space, segment_lengths, n_runs, starts)
end

function sample_segment_starts(exp::RandomSubsegmentExperiment; resample_segment_start::Bool=true)
    return sample_segment_starts(
        get_trajectory_space(exp),
        get_segment_lengths(exp),
        exp.n_runs,
        MersenneTwister(exp.seed);
        resample_segment_start,
    )
end

function _validate_segment_starts(trajectory_space::TrajectorySpace, segment_lengths, n_runs, segment_starts)
    length(segment_starts) == length(segment_lengths) ||
        throw(ArgumentError("segment_starts must have one entry per segment length. Got $(length(segment_starts)) and $(length(segment_lengths))."))

    n_time_steps = length(time_domain(get_trajectory(trajectory_space)))
    starts = Vector{Vector{Int}}(undef, length(segment_lengths))

    for i in eachindex(segment_lengths)
        len = segment_lengths[i]
        len >= 1 || throw(ArgumentError("segment lengths must be positive. Got $len."))
        max_start = n_time_steps - len + 1
        max_start >= 1 || throw(ArgumentError("segment_length must not exceed the number of time steps. Got $len and $n_time_steps."))

        cur_starts = collect(Int, segment_starts[i])
        length(cur_starts) == n_runs ||
            throw(ArgumentError("segment_starts[$i] must have length $n_runs. Got $(length(cur_starts))."))

        for start in cur_starts
            if start < 1 || start > max_start
                throw(ArgumentError("segment start $start for segment length $len is outside 1:$max_start."))
            end
        end

        starts[i] = cur_starts
    end

    return starts
end

struct TimedRandomSubsegmentResult{T<:RandomSubsegmentResult}
    result::T
    backend::Symbol
    elapsed_ns::Vector{Vector{UInt64}}
    nthreads::Int
end

get_experiment_result(timed::TimedRandomSubsegmentResult) = timed.result
get_elapsed_ns(timed::TimedRandomSubsegmentResult) = timed.elapsed_ns

function elapsed_seconds(timed::TimedRandomSubsegmentResult)
    return map(timed.elapsed_ns) do ns
        Float64.(ns) ./ 1.0e9
    end
end

function total_time_seconds(timed::TimedRandomSubsegmentResult)
    return sum(sum, timed.elapsed_ns; init = UInt64(0)) / 1.0e9
end

function _all_elapsed_seconds(timed::TimedRandomSubsegmentResult)
    seconds = elapsed_seconds(timed)
    return reduce(vcat, seconds; init = Float64[])
end

function mean_time_seconds(timed::TimedRandomSubsegmentResult)
    seconds = _all_elapsed_seconds(timed)
    return isempty(seconds) ? NaN : mean(seconds)
end

function median_time_seconds(timed::TimedRandomSubsegmentResult)
    seconds = _all_elapsed_seconds(timed)
    return isempty(seconds) ? NaN : median(seconds)
end

function _backend_symbol(::Val{B}) where {B}
    return B
end

function _backend_symbol(alg)
    return Symbol(string(alg))
end

function _warmup_signature(alg, traj_space, segment_lengths, segment_starts, threshold, field)
    isempty(segment_lengths) && return nothing
    isempty(segment_starts[1]) && return nothing

    len = segment_lengths[1]
    start = segment_starts[1][1]
    cur_range = start:(start + len - 1)
    return cycling_signature(alg, traj_space, cur_range, threshold; field)
end

function run_timed_experiment(exp::RandomSubsegmentExperiment;
    field=DEFAULT_FIELD,
    alg=Val(:DistanceMatrix),
    threshold=nothing,
    segment_starts=nothing,
    resample_segment_start=true,
    progress=true,
    parallel_inner=false,
    warmup=true)

    traj_space = get_trajectory_space(exp)
    seg_lengths = get_segment_lengths(exp)
    n_runs = exp.n_runs
    resolved_threshold = _resolve_r_max(traj_space, threshold)

    if segment_starts === nothing
        segment_starts = sample_segment_starts(exp; resample_segment_start)
    else
        segment_starts = _validate_segment_starts(traj_space, seg_lengths, n_runs, segment_starts)
    end

    if warmup
        _warmup_signature(alg, traj_space, seg_lengths, segment_starts, resolved_threshold, field)
    end

    signatures = [Vector{CyclingSignature}(undef, n_runs) for _ in seg_lengths]
    elapsed_ns = [Vector{UInt64}(undef, n_runs) for _ in seg_lengths]

    outer_iter = progress ? ProgressBar(enumerate(exp.segment_lengths)) : enumerate(exp.segment_lengths)

    for (i, len) in outer_iter
        if parallel_inner
            @threads for j in 1:n_runs
                cur_range = segment_starts[i][j]:(segment_starts[i][j] + len - 1)
                start_ns = time_ns()
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, resolved_threshold; field)
                elapsed_ns[i][j] = time_ns() - start_ns
            end
        else
            for j in 1:n_runs
                cur_range = segment_starts[i][j]:(segment_starts[i][j] + len - 1)
                start_ns = time_ns()
                signatures[i][j] = cycling_signature(alg, traj_space, cur_range, resolved_threshold; field)
                elapsed_ns[i][j] = time_ns() - start_ns
            end
        end
    end

    result = RandomSubsegmentResult(
        traj_space,
        seg_lengths,
        n_runs,
        segment_starts,
        resolved_threshold,
        signatures,
    )

    return TimedRandomSubsegmentResult(result, _backend_symbol(alg), elapsed_ns, nthreads())
end

function run_paired_timed_experiments(exp::RandomSubsegmentExperiment, backends;
    field=DEFAULT_FIELD,
    threshold=nothing,
    resample_segment_start=true,
    progress=true,
    parallel_inner=false,
    warmup=true)

    starts = sample_segment_starts(exp; resample_segment_start)
    results = Dict{Symbol,TimedRandomSubsegmentResult}()

    for alg in backends
        timed = run_timed_experiment(
            exp;
            alg,
            field,
            threshold,
            segment_starts = starts,
            progress,
            parallel_inner,
            warmup,
        )
        results[timed.backend] = timed
    end

    return results
end

struct SignatureAgreement
    dimension_match::Bool
    birth_vector_match::Bool
    cycling_space_match::Bool
    left_dimension::Int
    right_dimension::Int
    max_birth_discrepancy::Float64
    reason::String
end

struct ExperimentAgreement{S}
    segment_lengths::S
    n_runs::Int
    segment_starts::Vector{Vector{Int}}
    comparisons::Vector{NamedTuple}
end

function _birth_discrepancy(left, right)
    length(left) == length(right) || return Inf
    isempty(left) && return 0.0

    discrepancies = map(zip(left, right)) do (l, r)
        if l == r
            return 0.0
        end
        return abs(Float64(l) - Float64(r))
    end
    return maximum(discrepancies)
end

function _birth_vectors_match(left, right, atol, rtol)
    length(left) == length(right) || return false
    return all(isapprox(l, r; atol, rtol) for (l, r) in zip(left, right))
end

function compare_signatures(left::CyclingSignature, right::CyclingSignature;
    r=nothing,
    atol=1e-8,
    rtol=1e-8)

    left_dimension = r === nothing ? dimension(left) : dimension(left, r)
    right_dimension = r === nothing ? dimension(right) : dimension(right, r)
    dimension_match = left_dimension == right_dimension

    left_births = r === nothing ? birth_vector(left) : birth_vector(left; r)
    right_births = r === nothing ? birth_vector(right) : birth_vector(right; r)
    birth_vector_match = _birth_vectors_match(left_births, right_births, atol, rtol)
    max_birth_discrepancy = _birth_discrepancy(left_births, right_births)

    left_matrix = r === nothing ? cycling_matrix(left) : cycling_matrix(left; r)
    right_matrix = r === nothing ? cycling_matrix(right) : cycling_matrix(right; r)
    cycling_space_match = colspace_normal_form(left_matrix) == colspace_normal_form(right_matrix)

    reasons = String[]
    dimension_match || push!(reasons, "dimension mismatch: $left_dimension != $right_dimension")
    birth_vector_match || push!(reasons, "birth-vector mismatch")
    cycling_space_match || push!(reasons, "cycling-space mismatch")
    reason = isempty(reasons) ? "match" : join(reasons, "; ")

    return SignatureAgreement(
        dimension_match,
        birth_vector_match,
        cycling_space_match,
        left_dimension,
        right_dimension,
        max_birth_discrepancy,
        reason,
    )
end

function _validate_comparable_results(left::RandomSubsegmentResult, right::RandomSubsegmentResult)
    left.segment_lengths == right.segment_lengths ||
        throw(ArgumentError("results must have identical segment_lengths."))
    left.n_runs == right.n_runs ||
        throw(ArgumentError("results must have identical n_runs."))
    left.segment_starts == right.segment_starts ||
        throw(ArgumentError("results must have identical segment_starts."))
    length(left.signatures) == length(right.signatures) ||
        throw(ArgumentError("results must have the same number of signature groups."))

    for i in eachindex(left.signatures)
        length(left.signatures[i]) == length(right.signatures[i]) ||
            throw(ArgumentError("signature groups must have matching lengths at index $i."))
    end
end

function compare_experiment_results(left::RandomSubsegmentResult, right::RandomSubsegmentResult;
    r=nothing,
    atol=1e-8,
    rtol=1e-8)

    _validate_comparable_results(left, right)
    comparisons = NamedTuple[]

    for (i, len) in enumerate(left.segment_lengths)
        for j in 1:left.n_runs
            agreement = compare_signatures(
                left.signatures[i][j],
                right.signatures[i][j];
                r,
                atol,
                rtol,
            )
            push!(
                comparisons,
                (
                    segment_length_index = i,
                    segment_length = len,
                    run_index = j,
                    start_index = left.segment_starts[i][j],
                    agreement = agreement,
                ),
            )
        end
    end

    return ExperimentAgreement(left.segment_lengths, left.n_runs, left.segment_starts, comparisons)
end

function compare_experiment_results(left::TimedRandomSubsegmentResult, right::TimedRandomSubsegmentResult; kwargs...)
    return compare_experiment_results(left.result, right.result; kwargs...)
end

function _summary_stat(f, xs)
    return isempty(xs) ? NaN : f(xs)
end

function summarize_timings(timed::TimedRandomSubsegmentResult)
    seconds = elapsed_seconds(timed)

    return map(eachindex(timed.result.segment_lengths)) do i
        cur = seconds[i]
        (
            backend = timed.backend,
            segment_length = timed.result.segment_lengths[i],
            n_runs = length(cur),
            total_time_seconds = sum(cur; init = 0.0),
            mean_time_seconds = _summary_stat(mean, cur),
            median_time_seconds = _summary_stat(median, cur),
            std_time_seconds = length(cur) <= 1 ? 0.0 : std(cur),
            min_time_seconds = isempty(cur) ? NaN : minimum(cur),
            max_time_seconds = isempty(cur) ? NaN : maximum(cur),
            nthreads = timed.nthreads,
        )
    end
end

function summarize_agreement(comparison::ExperimentAgreement)
    return map(eachindex(comparison.segment_lengths)) do i
        rows = filter(row -> row.segment_length_index == i, comparison.comparisons)
        discrepancies = map(row -> row.agreement.max_birth_discrepancy, rows)
        finite_discrepancies = filter(isfinite, discrepancies)
        max_birth_discrepancy = if isempty(discrepancies)
            NaN
        elseif isempty(finite_discrepancies)
            Inf
        else
            maximum(finite_discrepancies)
        end

        (
            segment_length = comparison.segment_lengths[i],
            n_runs = length(rows),
            rank_mismatches = count(row -> !row.agreement.dimension_match, rows),
            birth_vector_mismatches = count(row -> !row.agreement.birth_vector_match, rows),
            cycling_space_mismatches = count(row -> !row.agreement.cycling_space_match, rows),
            max_birth_discrepancy = max_birth_discrepancy,
        )
    end
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
    if isempty(ct_dict)
        return Matrix{Int}[], zeros(Int, 0, length(result.segment_lengths))
    end
    return lenght_countdict_to_countmatrix(ct_dict; n_relevant = n_subspaces)
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
    if isempty(ct_dict)
        return Matrix{Int}[], zeros(Int, 0, length(result.segment_lengths))
    end
    return lenght_countdict_to_countmatrix(ct_dict; n_relevant = n_subspaces)
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
        int[i] = Tuple{Float64,Float64}[]
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

function _level_predicate(relation::Symbol)
    if relation == :geq
        return >=
    elseif relation == :leq
        return <=
    end
    throw(ArgumentError("relation must be :geq or :leq. Got $relation."))
end

function _level_intervals(
    f::StepFunction,
    level,
    predicate,
    r_min,
    r_max,
)::Vector{Tuple{Float64,Float64}}
    lo = Float64(r_min)
    hi = Float64(r_max)
    lo <= hi || throw(ArgumentError("r_min must be <= r_max. Got $lo and $hi."))

    if lo == hi
        return predicate(f(lo), level) ? [(lo, hi)] : Tuple{Float64,Float64}[]
    end

    breakpoints = Float64[lo]
    for x in f.xs
        xx = Float64(x)
        if isfinite(xx) && xx > lo && xx < hi
            push!(breakpoints, xx)
        end
    end
    push!(breakpoints, hi)

    intervals = Tuple{Float64,Float64}[]
    open_start = nothing
    for i in 1:(length(breakpoints) - 1)
        a = breakpoints[i]
        b = breakpoints[i + 1]
        if predicate(f(a), level)
            if open_start === nothing
                open_start = a
            end
        elseif open_start !== nothing
            push!(intervals, (open_start, a))
            open_start = nothing
        end

        if i == length(breakpoints) - 1 && open_start !== nothing
            push!(intervals, (open_start, b))
        end
    end

    return intervals
end

"""
    cycspace_level_intervals(fs::AbstractVector{<:StepFunction}, level;
        relation=:geq, r_min=0.0, r_max=Inf)

Return one vector of radius intervals for each step function in `fs`. With `relation=:geq`,
intervals are the radii where `f(r) >= level`; with `relation=:leq`, intervals are the radii
where `f(r) <= level`. Intervals are clipped to `[r_min, r_max]`.
"""
function cycspace_level_intervals(
    fs::AbstractVector{<:StepFunction},
    level;
    relation::Symbol = :geq,
    r_min = 0.0,
    r_max = Inf,
)
    predicate = _level_predicate(relation)
    return Vector{Vector{Tuple{Float64,Float64}}}(
        [_level_intervals(f, level, predicate, r_min, r_max) for f in fs],
    )
end

"""
    cycspace_level_intervals(result::RandomSubsegmentResult, V, level;
        relation=:geq, r_min=0.0, r_max=result.flt_threshold)

Compute `cycspace_distribution(result, V)` and return the radius intervals where each segment
length entry satisfies the requested level relation.
"""
function cycspace_level_intervals(
    result::RandomSubsegmentResult,
    V,
    level;
    relation::Symbol = :geq,
    r_min = 0.0,
    r_max = result.flt_threshold,
)
    return cycspace_level_intervals(
        cycspace_distribution(result, V),
        level;
        relation = relation,
        r_min = r_min,
        r_max = r_max,
    )
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

function _heatmap_bin_centers(x_min::Real, x_max::Real, n_bins::Integer)
    n_bins > 0 || throw(ArgumentError("n_bins must be positive. Got $n_bins."))

    x0 = Float64(x_min)
    x1 = Float64(x_max)
    x0 <= x1 || throw(ArgumentError("x_min must be <= x_max. Got $x0 and $x1."))

    if isapprox(x0, x1)
        return fill(x0, n_bins)
    end

    edges = range(x0, x1; length = n_bins + 1)
    return collect(@. (edges[1:end-1] + edges[2:end]) / 2)
end

"""
    _cycspace_distribution_heatmap_data(result::RandomSubsegmentResult, V;
        r_max=nothing, radius_bins=nothing)

Return `(segment_spans, radii, Z)` for the heatmap used by `plot_cycspace_distribution`.
`segment_spans` are the distinct segment lengths in their first-occurrence order. `radii`
contains the centers of `radius_bins` equally spaced bins over `[0, r_max]`. `Z[i, j]` is the
number of segments with time span `segment_spans[j]` whose cycling-space step function for `V`
is active at `radii[i]`. Repeated entries in `result.segment_lengths` are aggregated into the
same column.
"""
function _cycspace_distribution_heatmap_data(
    result::RandomSubsegmentResult,
    V;
    r_max = nothing,
    radius_bins = nothing,
)
    if r_max === nothing
        r_max = result.flt_threshold
    end
    if radius_bins === nothing
        radius_bins = length(unique(result.segment_lengths))
    end

    segment_spans = unique(result.segment_lengths)
    radii = _heatmap_bin_centers(0.0, r_max, radius_bins)
    dists = cycspace_distribution(result, V)

    span_to_idx = Dict(span => i for (i, span) in enumerate(segment_spans))
    Z = zeros(Int, length(radii), length(segment_spans))

    for (seg_idx, span) in enumerate(result.segment_lengths)
        col_idx = span_to_idx[span]
        f = dists[seg_idx]
        for (row_idx, r) in enumerate(radii)
            Z[row_idx, col_idx] += Int(f(r))
        end
    end

    return segment_spans, radii, Z
end

function _cycspace_interval_dict(result::RandomSubsegmentResult, k)
    int_dict = Dict{Matrix{Int},Vector{Tuple{Float64,Float64}}}()
    for sigs in result.signatures
        for sig in sigs
            if length(sig.birth_vector) >= k
                V = Int.(colspace_normal_form(sig.cycling_matrix[:, 1:k]))
                p1, p2 = dimension_interval(sig, k)
                if p1 <= p2
                    push!(get!(int_dict, V, Tuple{Float64,Float64}[]), (Float64(p1), Float64(p2)))
                end
            end
        end
    end
    return int_dict
end

function _cycspace_total_persistence(int_vec::Vector{Tuple{Float64,Float64}}, y_max::Real)
    y = Float64(y_max)
    return sum(int_vec, init = 0.0) do (p1, p2)
        return max(0.0, min(p2, y) - p1)
    end
end

"""
    cycspace_radius_frequency_functions(result::RandomSubsegmentResult, k;
        r_max_for_sorting=nothing, filter_shorter_as=0, max_n_sig=nothing)

Return `(sig, fs)` where `sig` are `k`-dimensional cycling subspaces ordered as in the old
`allSignatureRadiusFunctions` pipeline and `fs[i]` is the total radius-frequency step function
for `sig[i]` (summed over all segment lengths).
"""
function cycspace_radius_frequency_functions(
    result::RandomSubsegmentResult,
    k;
    r_max_for_sorting = nothing,
    filter_shorter_as = 0,
    max_n_sig = nothing,
)
    int_dict = _cycspace_interval_dict(result, k)
    keys_vec = collect(keys(int_dict))

    if r_max_for_sorting !== nothing
        scores = map(keys_vec) do V
            _cycspace_total_persistence(int_dict[V], r_max_for_sorting)
        end
        p = sortperm(scores, rev = true)
        keys_vec = keys_vec[p]
    else
        counts = map(V -> length(int_dict[V]), keys_vec)
        p = sortperm(counts, rev = true)
        keys_vec = keys_vec[p]
    end

    if max_n_sig !== nothing
        keys_vec = keys_vec[1:min(length(keys_vec), max_n_sig)]
    end

    fs = map(keys_vec) do V
        int_vec = int_dict[V]
        if filter_shorter_as > 0
            int_vec = filter(t -> t[2] - t[1] >= filter_shorter_as, int_vec)
        end
        return sum(int_vec, init = indicator_function(1, 0)) do (p1, p2)
            indicator_function(p1, p2, 1.0)
        end
    end

    return keys_vec, fs
end

"""
    cycspace_length_interval_countmatrix(result::RandomSubsegmentResult, k;
        n_subspaces=nothing, sort_by_tp_with_rmax=nothing,
        filter_shorter_than=0, filter_shorter_r_max=Inf)

Return `(sig, M)` where `M[i,j]` is the number of segments of length `segment_lengths[j]` whose
`k`-dimensional cycling subspace equals `sig[i]`, counting only intervals with length at least
`filter_shorter_than` (with interval endpoints clipped at `filter_shorter_r_max`).
Ordering matches the old `signaturesLengthToIntervalCount` behavior.
"""
function cycspace_length_interval_countmatrix(
    result::RandomSubsegmentResult,
    k;
    n_subspaces = nothing,
    sort_by_tp_with_rmax = nothing,
    filter_shorter_than = 0,
    filter_shorter_r_max = Inf,
)
    int_dict = _cycspace_interval_dict(result, k)
    keys_vec = collect(keys(int_dict))
    n_lengths = length(result.segment_lengths)

    ct_dict = Dict{Matrix{Int},Vector{Int}}(V => zeros(Int, n_lengths) for V in keys_vec)
    for (i, sigs) in enumerate(result.signatures)
        for sig in sigs
            if length(sig.birth_vector) >= k
                p1, p2 = dimension_interval(sig, k)
                if p1 <= p2
                    d = max(0.0, min(Float64(p2), Float64(filter_shorter_r_max)) - Float64(p1))
                    if d >= filter_shorter_than
                        V = Int.(colspace_normal_form(sig.cycling_matrix[:, 1:k]))
                        ct_dict[V][i] += 1
                    end
                end
            end
        end
    end

    if sort_by_tp_with_rmax !== nothing
        scores = map(keys_vec) do V
            _cycspace_total_persistence(int_dict[V], sort_by_tp_with_rmax)
        end
        p = sortperm(scores, rev = true)
        keys_vec = keys_vec[p]
    else
        scores = map(V -> sum(ct_dict[V]), keys_vec)
        p = sortperm(scores, rev = true)
        keys_vec = keys_vec[p]
    end

    if n_subspaces !== nothing
        keys_vec = keys_vec[1:min(length(keys_vec), n_subspaces)]
    end

    M = zeros(Int, length(keys_vec), n_lengths)
    for (i, V) in enumerate(keys_vec)
        M[i, :] = ct_dict[V]
    end
    return keys_vec, M
end

"""
    cycspace_total_distribution(result::RandomSubsegmentResult, V)

Return a step function `f(r)` giving the total number of segments (summed over all segment lengths)
whose `k = size(V,2)`-dimensional cycling space equals `V` at radius `r`.
"""
function cycspace_total_distribution(result::RandomSubsegmentResult, V)
    int_vecs = cycspace_intervals(result, V)
    all_ints = reduce(vcat, int_vecs; init = Tuple{Float64,Float64}[])
    return sum(all_ints, init = indicator_function(1, 0)) do (p1, p2)
        indicator_function(p1, p2, 1.0)
    end
end

"""
    cycspace_total_distributions(result::RandomSubsegmentResult, k; n_subspaces=nothing)

Compute total radius-distributions for the most frequent `k`-dimensional cycling spaces.
Returns `(sig, fs)` where `sig` is a vector of subspace matrices and `fs[i]` is the corresponding
total distribution step function.
"""
function cycspace_total_distributions(result::RandomSubsegmentResult, k; n_subspaces = nothing)
    sig, _ = cycspace_length_countmatrix(result, k; n_subspaces = n_subspaces)
    fs = map(sig) do V
        cycspace_total_distribution(result, V)
    end
    return sig, fs
end

function _step_post_xy(f::StepFunction, x_stop::Real; x_start::Real = 0.0)
    x_start > x_stop && throw(ArgumentError("x_start must be <= x_stop. Got $x_start and $x_stop."))

    xs = Float64[x_start]
    ys = Float64[float(f(x_start))]

    for x in f.xs
        xx = float(x)
        if isfinite(xx) && xx > x_start && xx < x_stop
            push!(xs, xx)
            push!(ys, float(f(xx)))
        end
    end

    push!(xs, float(x_stop))
    push!(ys, ys[end])
    return xs, ys
end

function _signature_start_indices(result::RandomSubsegmentResult, V)
    sig_vecs = result.signatures
    seg_starts = result.segment_starts
    k = size(V, 2)
    V = colspace_normal_form(V)

    bool_vecs = map(sig_vecs) do sigs
        map(sigs) do sig
            return length(sig.birth_vector) >= k && colspace_normal_form(sig.cycling_matrix[:, 1:k]) == V
        end
    end

    return map(zip(seg_starts, bool_vecs)) do (starts, bools)
        [starts[i] for i in eachindex(bools) if bools[i]]
    end
end

function _segment_span(start::Integer, len::Integer)
    return start:start + len - 1
end

"""
    signature_time_spans(result::RandomSubsegmentResult, V)

Return the time spans of all sampled segments whose cycling signature has the same column space as
`V`. The outer vector matches `result.segment_lengths`; each entry contains the matching sampled
time spans for that segment length.
"""
function signature_time_spans(result::RandomSubsegmentResult, V)
    seg_start_vecs = _signature_start_indices(result, V)
    return map(zip(seg_start_vecs, result.segment_lengths)) do (starts, len)
        [_segment_span(start, len) for start in starts]
    end
end

"""
    signature_segments(result::RandomSubsegmentResult, V)

Return the sampled trajectory segments whose cycling signature has the same column space as `V`.
The outer vector matches `result.segment_lengths`; each entry contains the extracted segments for
that segment length, as returned by `evaluate_interval`.
"""
function signature_segments(result::RandomSubsegmentResult, V)
    traj = get_trajectory(result.trajectory_space)
    spans = signature_time_spans(result, V)

    return map(spans) do cur_spans
        [evaluate_interval(traj, first(span), last(span)) for span in cur_spans]
    end
end

"""
    cycspace_segments(result::RandomSubsegmentResult, V)

Get the segments from the result that generate the cycling space V for some radius.
"""
function cycspace_segments(result::RandomSubsegmentResult, V)
    return signature_time_spans(result, V)
end

"""
    cycspace_segments_at_r(result::RandomSubsegmentResult, V, r)

Get the segments from the result have cycling space V at radius r.
"""
function cycspace_segments_at_r(result::RandomSubsegmentResult, V, r)
    # TODO: use `r` to filter by active radii once the semantics of this API are clarified.
    return signature_time_spans(result, V)
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
