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
        return Int.(subspaceNormalForm(mat[:,1:d]))
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
            basicMatrixReduction!(A)
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
            basicMatrixReduction!(A)
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
    basicMatrixReduction!(A)
    return all(==(0), A[:,size(V,2)+1:end])
end

function getSignatureRanges(res::SubsegmentResultReduced, sig, r)
    _,n = size(sig)
    ind1 = findall(v -> count(v .<= r) == n, res.birthVectors)
    ind2 = findall(i -> res.cyclingMatrices[i][:,1:n] == sig, ind1)
    return res.experiment.segmentRanges[ind1[ind2]]
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
