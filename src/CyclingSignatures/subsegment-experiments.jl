struct SubsegmentSampleParameter
    segLength::Int
    nSegments::Int
end

Base.show(io::IO, p::SubsegmentSampleParameter) = print(io, "Parameter: segLength=$(p.segLength), nSegments=$(p.nSegments)")

function sampleSegments(p::SubsegmentSampleParameter, range)
    a,b = first(range), last(range)-p.segLength+1
    
    starts = rand(a:b, p.nSegments)
    return map(starts) do i
        return i:i+p.segLength-1
    end
end

struct RandomSubsegmentExperiment
    trajectorySpace::TrajectorySpace
    sampleParameter::SubsegmentSampleParameter
    fltThreshold::Real

    segmentRanges::Vector{UnitRange}
end

function getTrajectorySpace(ex::RandomSubsegmentExperiment)
    return ex.trajectorySpace
end

function getSegmentRanges(ex::RandomSubsegmentExperiment)
    return ex.segmentRanges
end

function RandomSubsegmentExperiment(trajectorySpace, segmentsParam, fltThreshold=nothing)
    if fltThreshold === nothing
        fltThreshold = getBoxSpace(trajectorySpace).boxsize
    end
    segments = sampleSegments(segmentsParam, 1:nTimeSteps(getTrajectory(trajectorySpace)))
    return RandomSubsegmentExperiment(trajectorySpace, segmentsParam, fltThreshold, segments)
end

function sameSegmentExperiment(trajectorySpace, ex::RandomSubsegmentExperiment, fltThreshold=nothing)
    if fltThreshold === nothing
        fltThreshold = getBoxSpace(trajectorySpace).boxsize
    end
    return RandomSubsegmentExperiment(trajectorySpace, ex.sampleParameter, fltThreshold, ex.segmentRanges)
end

function Base.show(io::IO, ex::RandomSubsegmentExperiment)
    print(io, "RandomSubsegmentExperiment: ")
    print(io, "fltThreshold=$(ex.fltThreshold), ")
    show(io, ex.sampleParameter)
end

function getParameter(ex::RandomSubsegmentExperiment)
    return ex.sampleParameter
end

struct RandomSubsegmentResult
    experiment
    algorithm
    cyclingDiagrams
end

function Base.show(io::IO, res::RandomSubsegmentResult)
    print(io, "Result for ")
    show(io, res.experiment)
end

function runExperiment(ex::RandomSubsegmentExperiment, algorithm=Val(:DistanceMatrix); verbose=true, kwargs...)
    ts = getTrajectorySpace(ex)
    segments = getSegmentRanges(ex)
    iter = segments

    verbose && (iter = ProgressBar(segments))
    cyclingDiagrams = map(iter) do seg
        verbose && (set_description(iter, string(algorithm)*", "*string(getParameter(ex))))
        if haskey(kwargs, :field)
            return evaluateCycling(algorithm, ts, seg, ex.fltThreshold, kwargs[:field])
        else
            return evaluateCycling(algorithm, ts, seg, ex.fltThreshold)
        end
    end

    return RandomSubsegmentResult(ex, algorithm, cyclingDiagrams)
end

function runExperiments(exs::Vector{RandomSubsegmentExperiment}, algorithm=Val(:DistanceMatrix); kwargs...)
    return map(exs) do ex
        runExperiment(ex, algorithm; kwargs...)
    end
end

function runExperimentsWithTimer(exs::Vector{RandomSubsegmentExperiment}, algorithm=Val(:DistanceMatrix); kwargs...)
    times = zeros(length(exs))
    results = map(enumerate(exs)) do t
        i,ex = t
        times[i] = @elapsed res = runExperiment(ex, algorithm; kwargs...)
        return res
    end
    return times, results
end

function getDiagrams(res::RandomSubsegmentResult)
    return res.cyclingDiagrams
end

struct SubsegmentResultReduced
    experiment::RandomSubsegmentExperiment
    algorithm
    cyclingMatrices::Vector{Matrix{<:Integer}}
    birthVectors::Vector{Vector{Float64}}
    SubsegmentResultReduced(experiment, algorithm, cyclingMatrices, birthVectors) = all(issorted.(birthVectors)) ? new(experiment, algorithm, cyclingMatrices, birthVectors) : error()
end

function SubsegmentResultReduced(res::RandomSubsegmentResult)
    exp = res.experiment
    alg = res.algorithm
    
    cyclingMatrices = map(res.cyclingDiagrams) do diag
        return diagToCyclingMat(diag)
    end
    birthVectors = map(res.cyclingDiagrams) do diag
        return map(birth, diag)
    end
    return SubsegmentResultReduced(exp, alg, cyclingMatrices, birthVectors)
end

function diagToCyclingMat(diag)
    if length(diag) == 0
        return Array{Int}(undef,0,0)
    end
    reps = map(diag) do int
        return int.inclusion_representative
    end
    return Matrix(hcat(reps...))
end

function getSampleParameter(res::RandomSubsegmentResult)
    return res.experiment.sampleParameter
end

function getSampleParameter(res::SubsegmentResultReduced)
    return res.experiment.sampleParameter
end

function getSegmentLength(res::RandomSubsegmentResult)
    return getSampleParameter(res).segLength
end

function getSegmentLength(res::SubsegmentResultReduced)
    return getSampleParameter(res).segLength
end

function getnSegments(res::RandomSubsegmentResult)
    return getSampleParameter(res).nSegments
end

function getnSegments(res::SubsegmentResultReduced)
    return getSampleParameter(res).nSegments
end

function rankCountmap(res, r=Inf)
    rankVec = map(v -> count(<=(r), v), res.birthVectors)
    return countmap(rankVec)
end

function getRankDistribution(res, r=Inf)
    cm = rankCountmap(res, r)
    maxRank = maximum(keys(cm))
    return map(i -> get(cm, i, 0), 0:maxRank)
end

function rankDistributionMatrix(results::Vector, r=Inf)
    cms = rankCountmap.(results, r)
    maxRank = maximum(cm -> maximum(keys(cm)), cms)

    rankMat = zeros(Int, maxRank+1, length(results))

    for (j,cm) in enumerate(cms)
        for i in 0:maxRank
            rankMat[i+1,j] = get(cm, i, 0)
        end
    end
    return rankMat
end

function subspaceDistribution(res, d, r=Inf)
    rank_d_ind = findall(res.birthVectors) do v
        return count(<=(r), v) == d
    end

    rank_d_subspaces = map(res.cyclingMatrices[rank_d_ind]) do mat
        return Int.(colspace_normal_form(mat[:,1:d]))
    end
    
    return isempty(rank_d_ind) ? Dict{Vector{Int},Int}() : countmap(rank_d_subspaces,alg=:dict)
end

function subspaceFrequencyMatrix(results::Vector,d,r=Inf;cutoff=0)
    cms = subspaceDistribution.(results, d, r)
    allkeys = union(keys.(cms)...)
    keys_vec = collect(allkeys)
    filter!(keys_vec) do k
        sum(cm-> get(cm,k,0), cms) >= cutoff
    end

    subspaceMat = zeros(Int, length(keys_vec), length(results))
    for i = 1:size(subspaceMat,1)
        for j = 1:size(subspaceMat,2)
            subspaceMat[i,j] = get(cms[j], keys_vec[i], 0)
        end
    end
    p = sortperm(map(sum,eachrow(subspaceMat)), rev=true)
    return keys_vec[p], subspaceMat[p,:]
end

function inclusionVectors(res::SubsegmentResultReduced, signatures)
    map(res.cyclingMatrices) do mat
        map(signatures) do v
            A = [mat v]
            basic_reduction!(A)
            return all(==(0), A[:,end])
        end
    end
end

function subspaceInclusionDistribution(res::SubsegmentResultReduced, signatures, d, r=Inf)
    ssD = subspaceDistribution(res, d, r)
    subspaces = collect(keys(ssD))
    subspaceFrequencies = map(k -> ssD[k], subspaces)
    inclusionVectors = map(subspaces) do mat
        map(signatures) do sig
            A = [mat sig]
            basic_reduction!(A)
            return all(==(0), A[:,d+1:end])
        end
    end
    return subspaces, subspaceFrequencies, inclusionVectors
end

function pairInclusionData(res::SubsegmentResultReduced, signatures, d, r=Inf)
    # return a matrix A such that A[i,j] contains the number of
    # d dimensional subspaces containing signatures i and j
    _, subspaceFrequencies, inclusionVectors = subspaceInclusionDistribution(res, signatures, d, r)
    n = length(signatures)

    A = zeros(Int,n,n)
    for i = 1:n
        for j = 1:n
            inc_ind = findall(inclusionVectors) do v
                v[i] == 1 && v[j] == 1
            end
            A[i,j] = sum(subspaceFrequencies[inc_ind])
        end
    end
    return A
end

function subspaceInclusionMatrix(sig1, sig2, F)
    m,n = length(sig2), length(sig1)
    A = zeros(Bool, m,n)

    for i = 1:m
        for j = 1:n
            A[i,j] = isSubspace(F.(sig2[i]), F.(sig1[j]))
        end
    end
    return A
end

function isSubspace(V, W)
    A = [V W]
    basic_reduction!(A)
    return all(==(0), A[:,size(V,2)+1:end])
end

function getSignatureRanges(res::SubsegmentResultReduced, sig, r)
    _,n = size(sig)
    ind1 = findall(v -> count(v .<= r) == n, res.birthVectors)
    ind2 = findall(i -> res.cyclingMatrices[i][:,1:n] == sig, ind1)
    return res.experiment.segmentRanges[ind1[ind2]]
end

function getSubspaceRanges(res::SubsegmentResultReduced, sig, r)
    _,n = size(sig)
    ind1 = findall(v -> count(v .<= r) == n, res.birthVectors)
    ind2 = findall(i -> colspace_normal_form(res.cyclingMatrices[i][:,1:n]) == sig, ind1)
    return res.experiment.segmentRanges[ind1[ind2]]
end

function getSubspaceRangesAndBars(res::SubsegmentResultReduced, sig)
    _,n = size(sig)
    ind = findall(M -> size(M,2)>=n && colspace_normal_form(M[:,1:n]) == sig, res.cyclingMatrices)
    ranges = res.experiment.segmentRanges[ind]
    bars = map(ind) do i
        b = res.birthVectors[i][n]
        d = length(res.birthVectors[i])==n ? Inf : res.birthVectors[i][n+1]
        return (b,d)
    end
    return ranges,bars
end

struct SubsegmentResultSummary
    sampleParameter::SubsegmentSampleParameter
    fltThreshold::Real
    cyclingMatrices::Vector{Matrix{<:Integer}}
    birthVectors::Vector{Vector{Float64}}

    SubsegmentResultSummary(sampleParameter, fltThreshold, cyclingMatrices, birthVectors) = all(issorted.(birthVectors)) ? new(sampleParameter, fltThreshold, cyclingMatrices, birthVectors) : error()
end

function SubsegmentResultSummary(res::SubsegmentResultReduced)
    param = getSampleParameter(res)
    thrsh = res.experiment.fltThreshold
    cycMat = res.cyclingMatrices
    bV = res.birthVectors

    return SubsegmentResultSummary(param, thrsh, cycMat, bV)
end

function SubsegmentResultSummary(res::RandomSubsegmentResult)
    return SubsegmentResultSummary(SubsegmentResultReduced(res))
end

# TODO: do not use ranges for indexing QuantizedTrajectories but rather some sort of interval description
# TODO: for future make the following changes
# an experiment should be able to deal with a variety of parameters.

###
### Code for radius dependent evaluation
###

function intervalsToFrequencyFunction(ints)
    births = map(t->t[1], ints)
    deaths = filter(isfinite, map(t->t[2],ints))
    xs_unsorted = [births;deaths]
    p = sortperm(xs_unsorted)
    sign_vec = [ones(Int,length(births));-ones(Int,length(deaths)) ]
    xs_sorted = xs_unsorted[p]
    ys_sorted = cumsum(sign_vec[p])

    # take care that xs is unique, for this, need last occurence for each x
    i_unique = unique(i->xs_sorted[i] , Iterators.reverse(eachindex(xs_sorted)))
    sort!(i_unique)
    return xs_sorted[i_unique], ys_sorted[i_unique]
end

"""
    function evaluateStepFct(x, xs, ys)

Evaluates the expression ∑ ys[i] * 1_{ [xs[i],xs[i+1]) }(x)
"""
function evaluateStepFct(x, xs, ys)
    if length(ys) == 0 || x <= xs[1]
        return 0
    end
    i = last(searchsorted(xs, x))
    return ys[i]
end

function rankIntervals(res, d)
    ind = findall(v->length(v) >= d, res.birthVectors)
    if d == 0
        return map(bv -> length(bv) == 0 ? (0,Inf) : (0,bv[1]),res.birthVectors[ind])
    else
        return map(bv -> length(bv) == d ? (bv[d],Inf) : (bv[d],bv[d+1]),res.birthVectors[ind])
    end
end

function rankIntervalStepFct(res, d)
    ints = rankIntervals(res, d)
    return intervalsToFrequencyFunction(ints)
end

function rankHeatmapData(results,d, y_rng)
    rk_fcts = rankIntervalStepFct.(results,d)

    #rk_fcts_all_max = maximum(rk_fcts) do fct
    #    maximum(fct[2],init=0)
    #end

    M = map(CartesianIndices((1:length(rk_fcts),1:length(y_rng)))) do t
        return evaluateStepFct(y_rng[t[2]], rk_fcts[t[1]]...)
    end
end

# for all d-dimensional subspaces, find all intervals where the subspace is taken
function subspaceIntervalDict(res,d)
    ind = findall(v->length(v) >= d, res.birthVectors)
    intervals = map(bv -> length(bv) == d ? (bv[d],Inf) : (bv[d],bv[d+1]),res.birthVectors[ind])
    sigMats = map(M-> M[:,1:d],res.cyclingMatrices[ind])
    d = Dict{Matrix{Int},Vector{Tuple{Float64,Float64}}}()
    for (int,M) in zip(intervals,sigMats)
        M_nf = Int.(colspace_normal_form(M))
        if !haskey(d,M_nf)
            d[M_nf] = [int]
        else
            push!(d[M_nf],int)
        end
    end
    return d
end

function subspaceFrequencyRadiusFct(res,d)
    d = subspaceIntervalDict(res,d)
    k = collect(keys(d))
    fs = map(ints -> intervalsToFrequencyFunction(ints), values(d))
    return k, fs
end

function stepFunctionIntegal(xs, ys)
    if length(xs) != length(ys)+1
        error("Invalid step function")
    end
    return sum(i-> ys[i]*(xs[i+1]-xs[i]), 1:length(ys))
end


function signatureHeatmapMatrix(results, sig, y_rng)
    si_dicts = subspaceIntervalDict.(results,size(sig,2))
    sig_intervals = map(d->haskey(d, sig) ? getindex(d, sig) : Float64[],si_dicts)
    cs_f_fcts = intervalsToFrequencyFunction.(sig_intervals)
    M = map(CartesianIndices((1:length(cs_f_fcts),1:length(y_rng)))) do t
        return evaluateStepFct(y_rng[t[2]], cs_f_fcts[t[1]]...)
    end
    return M
end

function allSignatureHeatmapMatrices(results, k, y_rng)
    si_dicts = subspaceIntervalDict.(results,k)
    keys_vec = collect(union(keys.(si_dicts)...))
    mats = map(keys_vec) do sig
        sig_intervals = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
        cs_f_fcts = intervalsToFrequencyFunction.(sig_intervals)
        M = map(CartesianIndices((1:length(cs_f_fcts),1:length(y_rng)))) do t
            return evaluateStepFct(y_rng[t[2]], cs_f_fcts[t[1]]...)
        end
        return M
    end
    p = sortperm(mats, by=sum, rev=true)
    return keys_vec[p], mats[p]
end

function signaturesTotalPersistence(results, k, y_max; p=1)
    si_dicts = subspaceIntervalDict.(results,k)
    keys_vec = collect(union(keys.(si_dicts)...))

    tp = map(keys_vec) do sig
        int_vecs = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
        all_intervals = vcat(int_vecs...)
        all_lengths = map(t->min(t[2],y_max)-t[1],all_intervals)
        return norm(all_lengths, p)
    end
    p = sortperm(tp, rev=true)
    return keys_vec[p], tp[p]
end

function allSignatureRadiusFunctions(results, k; r_max_for_sorting=nothing, filter_shorter_as=0, max_n_sig=nothing)
    si_dicts = subspaceIntervalDict.(results,k)
    keys_vec = collect(union(keys.(si_dicts)...))
    if r_max_for_sorting !== nothing
        keys_vec, _ = signaturesTotalPersistence(results, k, r_max_for_sorting)
    else
        keys_vec, _ = signaturesTotalIntervalCount(results, k)
    end
    if max_n_sig !== nothing
        keys_vec = keys_vec[1:min(length(keys_vec),max_n_sig)]
    end

    r_fcts = map(keys_vec) do sig
        int_vecs = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
        all_intervals = vcat(int_vecs...)
        filter!(i-> i[2]-i[1]>= filter_shorter_as, all_intervals)
        return intervalsToFrequencyFunction(all_intervals)
    end
    return keys_vec, r_fcts
end



function signaturesLengthToTotalPersistence(results, k, r_max)
    si_dicts = subspaceIntervalDict.(results,k)
    keys_vec = collect(union(keys.(si_dicts)...))
    l_fcts = map(keys_vec) do sig
        int_vecs = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
        # for every Vector{Tuple{Float64,Float64}} that represents a collection of intervals, calculate total persistence
        return map(int_vecs) do int_vec
            return sum(int_vec, init=0) do int
                return min(int[2], r_max) - int[1]
            end
        end
    end
    p = sortperm(l_fcts, rev=true, by=sum)
    return keys_vec[p], l_fcts[p]
end

function signaturesLengthToIntervalCount(results, k; sort_by_tp_with_rmax=nothing, filter_shorter_than=0, filter_shorter_r_max=Inf)
    si_dicts = subspaceIntervalDict.(results,k)
    keys_vec = collect(union(keys.(si_dicts)...))
    l_fcts = map(keys_vec) do sig
        int_vecs = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
        # for every Vector{Tuple{Float64,Float64}} that represents a collection of intervals, calculate total persistence
        int_vecs = map(int_vecs) do v
            filter(t-> min(t[2],filter_shorter_r_max)-t[1]>= filter_shorter_than, v)
        end
        return map(length,int_vecs)
    end
    p = sortperm(l_fcts, rev=true, by=sum)

    if !isnothing(sort_by_tp_with_rmax)
        l_fcts_sort_vec = map(keys_vec) do sig
            int_vecs = map(d->haskey(d, sig) ? getindex(d, sig) : Tuple{Float64,Float64}[],si_dicts)
            # for every Vector{Tuple{Float64,Float64}} that represents a collection of intervals, calculate total persistence
            return sum(int_vecs, init=0) do int_vec
                return sum(int_vec, init=0) do int
                    return min(int[2], sort_by_tp_with_rmax) - int[1]
                end
            end
        end    
        p = sortperm(l_fcts_sort_vec, rev=true, by=sum)
    end
    return keys_vec[p], l_fcts[p]
end

function signaturesTotalIntervalCount(results, k; filter_shorter_than=0, filter_shorter_r_max=Inf)
    keys_vec, l_fcts = signaturesLengthToIntervalCount(results, k; filter_shorter_than=filter_shorter_than, filter_shorter_r_max=filter_shorter_r_max)
    total_number = sum.(l_fcts, init=0)
    return keys_vec, total_number
end

function cyclingSpaceSegmentIndices(res, sig)
    d = size(sig,2)
    ind = findall(v->length(v) >= d, res.birthVectors)

    filter!(i-> colspace_normal_form(res.cyclingMatrices[i][:,1:d] ) == sig, ind)
    return ind
end