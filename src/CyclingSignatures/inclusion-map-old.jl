############################################################################################
############################################################################################
#
# OLD CODE:
#
############################################################################################
############################################################################################




function edge_boxes_sphere_bundle(p1, p2, v1, v2, sb_radius)
    lambdas = projection_max_switches(sb_radius*v1, sb_radius*v2)

    # in each interval [lambdas[i],lambdas[i+1]], the edge is projected to the same face
    #   => normal edge_boxes can be used
    boxes = Vector{Int}[]
    x1 = 1.0 * [p1;v1]
    x2 = 1.0 * [p1; sb_radius*normalize(v1,Inf)]
    for i in 2:length(lambdas)
        l = lambdas[i]
        x1 .= x2
        x2 .= [ (1-l)*p1 + l*p2; sb_radius*normalize((1-l)*v1 + l*v2,Inf) ]
        eb = edge_boxes(x1, x2)
        if !isempty(boxes) && boxes[end] == eb[1]
            append!(boxes, eb[2:end])
        else
            append!(boxes, eb)
        end
    end
    return boxes
end

@doc raw"""
    projection_max_switches(v1, v2)

Computes a vector of times `T` which partition the interval `[0,1]` such that
`λ \mapsto (1-λ)v_1 + λ v_2` has the same maximal component in each interval `[T[0];T[1]]`.
"""
function projection_max_switches(v1, v2)
    Ts = zeros(1)
    max_i = projection_initial_max_component(v1, v2)
    cur_t = 0
    t_swap_min = 1.0
    while Ts[end] < 1
        for i in 1:length(v1)
            t1 = intersect_lines(v1[max_i], v1[i], v2[max_i], v2[i])
            if cur_t < t1 && t1 < t_swap_min
                t_swap_min = t1
            end
            t2 = intersect_lines(v1[max_i], -v1[i], v2[max_i], -v2[i])
            if cur_t < t2 && t2 < t_swap_min
                t_swap_min = t2
            end
            # if either intersection point is in [t_swap_min,1] make it the new intersection point
        end
        # update
        v = (1-t_swap_min)*v1 + t_swap_min*v2
        max_i = projection_initial_max_component(v,v2)
        # reset
        push!(Ts, t_swap_min)
        cur_t = t_swap_min
        t_swap_min = 1.0
    end

    return Ts
end

function projection_initial_max_component(v1, v2)
    max_v1 = maximum(abs, v1)
    is = findall(x-> isapprox(abs(x), max_v1), v1)
    return argmax(i-> (sign(v1[i])>= 0 ? 1 : -1)*v2[i], is)
end

"""
    intersect_lines(a,b)

Computes the intersection of the lines `(1-l)*a[i]+l*b[i]`
"""
function intersect_lines(a,b)
    return intersect_lines(a[1,], a[2], b[1], b[2]) # note condition not great, maybe use geometrybasics implementation
end

function intersect_lines(a1, a2, b1, b2)
    return (a1-a2)/(b2-b1-a2+a1) # note condition not great, maybe use geometrybasics implementation
end











############################################################################################
############################################################################################
#
# OLD CODE:
#
############################################################################################
############################################################################################




"""
    edgeBoxesSphereBundle(p1, p2, v1, v2, radius, sphereBundleRadius)

Returns a vector of all boxes between (p1,v1) and (p2,v2) in the sphere bundle.

The sphere bundle consists of
- boxes in space have edge length `boxsize`,
- boxes in sphere bundle are taken in the sphere with radius `sphereBundleRadius` and have edge length one.

In the computation, boxsize and sphereBundleRadius are accounted for as follows:
- boxsize is the edge length of the box, it is assumed that
"""
function edgeBoxesSphereBundle(p1, p2, v1, v2, sphereBundleRadius=1)
    allBoxes = Vector{Vector{Int}}[]
    # the projection of a line in the sphere bundle is a piecwise linear path
    while (λ = nextMaxSwitchLambda(v1, v2)) < 1
        pNew = (1 - λ) * p1 + λ * p2
        vNew = normalize((1 - λ) * v1 + λ * v2, Inf)
        push!(allBoxes, edge_boxes([p1; sphereBundleRadius * v1], [pNew; sphereBundleRadius * vNew])[1:end-1])
        p1 = pNew
        v1 = vNew
    end
    push!(allBoxes, edge_boxes([p1; sphereBundleRadius * v1], [p2; sphereBundleRadius * v2]))
    return reduce(vcat, allBoxes)
end

"""
    nextMaxSwitchLambda(p1, p2)

Given p1 and p2, returns the first value lambda in (0,1) such that two components of abs.(p1 + lambda*(p2-p1))
are maximal, or 0 if such a value does not exist.


"""
function nextMaxSwitchLambda(p1, p2)
    sign_v = replace(sign.(p1), 0 => 1)
    q1 = sign_v .* p1
    q2 = sign_v .* p2

    max_q1 = maximum(q1)
    max_q1_ind = findall(v -> isapprox(v, max_q1), q1)
    v_max, i_max_pre = findmax(k -> q2[k], max_q1_ind) # this line is max for [0,λ), need to find λ now

    i_max = max_q1_ind[i_max_pre]

    otherInd = findall(x -> abs(x) > v_max, q2)

    v1 = [q1[i_max]; 0]
    v2 = [q2[i_max]; 0]
    lambdas = map(otherInd) do j
        s = 1 - 2 * signbit(q2[j])
        v1[2] = s * q1[j]
        v2[2] = s * q2[j]
        return intersectionLambda(v1, v2)
    end
    filter!(l -> l > 0 && l < 1, lambdas)
    return minimum(lambdas, init=1)
end

function intersectionLambda(a, b)
    # computes the intersection of the lines (1-l)*a[i]+l*b[i]
    return (a[1] - a[2]) / (b[2] - b[1] - a[2] + a[1]) # note condition not great, maybe use geometrybasics implementation
end

"""
    allMaxSwitchLambda(p1,p2)

Given p1 and p2, returns all values lambda in (0,1) such that two components of abs.(p1 + lambda*(p2-p1))
are maximal. Returns a vector of Float64 which may be empty.
"""
function allMaxSwitchLambda(p1, p2)
    # TODO: refactor this method!!!
    lambdas = Float64[]
    lambdaNew = nextMaxSwitchLambda(p1, p2)
    while lambdaNew < 1
        if isempty(lambdas)
            push!(lambdas, lambdaNew)
        else
            lambdaRescaled = lambdas[end] + lambdaNew * (1 - lambdas[end])
            push!(lambdas, lambdaRescaled)
        end
        p1 += lambdaNew * (p2 - p1)
        lambdaNew = nextMaxSwitchLambda(p1, p2)
    end
    return lambdas
end

function countMaxSwitchLambdas(p1, p2)
    lCount = 0
    lambdaNew = nextMaxSwitchLambda(p1, p2)
    while lambdaNew < 1
        lCount += 1
        p1 += lambdaNew * (p2 - p1)
        lambdaNew = nextMaxSwitchLambda(p1, p2)
    end
    return lCount
end

###
### Code for inclusion helpers
###

abstract type AbstractInclusionHelper end

function incKwargs(h::AbstractInclusionHelper)
    return (; :filter_missing => h.filter_missing, :shortest_path_fix => h.shortest_path_fix, :safe_output => h.safe_output)
end

function inclusionMatrix(inc_helper::AbstractInclusionHelper)
    return inc_helper.h1_gen_t
end

function getComplex(inc_helper::AbstractInclusionHelper)
    return inc_helper.cplx
end

function getCplxPointsSorted(inc_helper::AbstractInclusionHelper)
    return inc_helper.cplxPointsSorted
end

function getCplxSortPerm(inc_helper::AbstractInclusionHelper)
    return inc_helper.cplxPointsSortPerm
end


mutable struct InclusionHelper{S<:AbstractVector{Int},T<:Real} <: AbstractInclusionHelper
    cplx::CellComplex
    cplxPointsSorted::Vector{S}
    cplxPointsSortPerm::Vector{Int}
    h1_gen_t

    # box radii
    boxsize::T

    filter_missing::Bool
    shortest_path_fix::Bool
    safe_output::Bool
end

function InclusionHelper(cplx, cplxPoints, h1_gen, boxsize; filter_missing=false, shortest_path_fix=false, safe_output=true)
    cplxPointsVec = collect(eachcol(cplxPoints))
    cpSortPerm = sortperm(cplxPointsVec)
    cplxPointsSorted = cplxPointsVec[cpSortPerm]

    return InclusionHelper(cplx, cplxPointsSorted, cpSortPerm, sparse(h1_gen'), boxsize, filter_missing, shortest_path_fix, safe_output)
end

"""
    struct SphereBundleInclusionHelper

Encapsules all complex information in order to find inclusion indices.

More precisely, an object of this type contains the information of a cubical complex in a sphere bundle such that we can solve:
given an edge in the union of the boxes, find an inclusion into the Nerve of the cubical cover.

Constructor is
    function SphereBundleInclusionHelper(cplx, cplxPoints, radius, sphereBundleRadius;
                                        filter_missing=false, shortest_path_fix=false, safe_output=true)
where
- cplx is a simplicial complex,
- cplxPoints is a not necessarily sorted vector of points in the complex,
- h1_gen is a matrix whose columns are generators for first cohomology of cplx,
- radius is the box radius in space direction,
- sphereBundleRadius is the radius of the sphere.
"""
mutable struct SBInclusionHelper{S<:AbstractVector{Int},T<:Real} <: AbstractInclusionHelper
    # elements for handling cplx and its points
    cplx::CellComplex
    cplxPointsSorted::Vector{S}
    cplxPointsSortPerm::Vector{Int}
    h1_gen_t

    # box radii
    boxsize::T
    sphereBundleRadius::Int

    filter_missing::Bool
    shortest_path_fix::Bool
    safe_output::Bool
end

function SBInclusionHelper(cplx, cplxPoints, h1_gen, boxsize, sphereBundleRadius;
    filter_missing=false, shortest_path_fix=false, safe_output=true)
    cplxPointsVec = collect(eachcol(cplxPoints))
    cpSortPerm = sortperm(cplxPointsVec)
    cplxPointsSorted = cplxPointsVec[cpSortPerm]

    return SBInclusionHelper(cplx, cplxPointsSorted, cpSortPerm, sparse(h1_gen'), boxsize, sphereBundleRadius,
        filter_missing, shortest_path_fix, safe_output)
end

function Base.show(io::IO, bs::SBInclusionHelper)
    print(io, "SBInclusionHelper")
end

function Base.show(io::IO, bs::SBInclusionHelper)
    print(io, "SBInclusionHelper")
end

function sbPointToComponents(p)
    d = div(length(p), 2)
    return p[1:d], p[d+1:end]
end

###
### Inclusion Map Code
###

function includeCycle(h, cyclePoints, cycleEdges::Vector{Tuple{Int,Int}}, coeffs::Vector{T}) where {T<:Integer}
    chainVec = spzeros(T, length(cells(getComplex(h))[1]))
    foreach(zip(cycleEdges, coeffs)) do t
        edge, coeff = t
        a, b = edge
        boxIndices = inclusionIndices(h, cyclePoints[:, a], cyclePoints[:, b])::Vector{Int}
        addEdgeBoxesToChain!(chainVec, boxIndices, coeff, getComplex(h))
    end
    return chainVec
end

"""
    function addEdgeBoxesToChain!(chain, boxIndices, coeff, cplx)

Chain is a dictionary indexed by cplx.cells[1]. This method adds the chain
sum_i coeff[i]*[boxIndices[i-1];boxIndices[i]]
"""
function addEdgeBoxesToChain!(chain, boxIndices, coeff, cplx)
    for i = 1:length(boxIndices)-1
        v1 = boxIndices[i]
        v2 = boxIndices[i+1]
        edge = Simplex([v1, v2])
        orientation = sign(v2 - v1)
        edgeIndex = collect(searchsorted(cplx.cells[1], edge))
        if length(edgeIndex) != 1
            @warn length(edgeIndex)
            @warn boxIndices[[i, i + 1]]
            error("A required edge is missing from the complex!")
        end
        chain[edgeIndex[1]] += orientation * coeff
    end
    return chain
end

###
### Inclusion indices
###

function inclusionIndices(h::AbstractInclusionHelper, p1, p2)
    boxes = edge_boxes(h, p1, p2)
    return getInclusionEdgeBoxIndices(boxes, getCplxPointsSorted(h), getCplxSortPerm(h); incKwargs(h)...)
end

function edge_boxes(h::AbstractInclusionHelper, p1, p2)
    @info typeof(h)
    error("Check type specialization in edgeBoxes")
end

function edge_boxes(h::InclusionHelper, p1, p2)
    return edge_boxes(p1, p2, h.boxsize)
end

function edge_boxes(h::SBInclusionHelper, sb_p1, sb_p2)
    p1, v1 = sbPointToComponents(sb_p1)
    p2, v2 = sbPointToComponents(sb_p2)
    p1 ./= h.boxsize
    p2 ./= h.boxsize
    normalize!(v1, Inf)
    normalize!(v2, Inf)

    # check if they are in same box and return the box if true

    p_dist = maximum(i -> abs(round(Int, p1[i]) - round(Int, p2[i])), 1:length(p1))
    v_dist = maximum(i -> abs(round(Int, h.sphereBundleRadius * v1[i]) - round(Int, h.sphereBundleRadius * v2[i])), 1:length(v1))

    if p_dist == 0 && v_dist == 0
        return [round.(Int, [p1; h.sphereBundleRadius * v1])]
    elseif p_dist <= 1 && v_dist <= 1
        return [round.(Int, [p1; h.sphereBundleRadius * v1]), round.(Int, [p2; h.sphereBundleRadius * v2])]
    end

    m_p, m_v = (p1 + p2) / 2, normalize!(v1 + v2, Inf)
    boxes1 = edgeBoxesSphereBundle(p1, m_p, v1, m_v, h.sphereBundleRadius)
    boxes2 = edgeBoxesSphereBundle(m_p, p2, m_v, v2, h.sphereBundleRadius)
    return vcat(boxes1, boxes2)
end

"""
    function getEdgeBoxIndices(boxes, cplxPointsVecSorted; allow_outside_detour=false, verbose=false)

Computes a vector of boxes indices in cplxPointsVecSorted which corresponds to the boxes in boxes.
By default, only the indices of the boxes are computed, and an error is returned if any box is not found.
The keyword arguments allow different levels of fixing this error:
- if filter_missing is set to true, boxes which can't be included will be ignored
  Example: boxes = [ [0;1], [1;1], [1,0] ], cplxPointsVecSorted = [ [0;1], [1,0] ] will return [1;2] even though [1;1] not in cplx
- if in addition, shortest_path_fix is set to true, attempt to fix dynamic inconsistencies using shortest path heuristic
- if safe_output, it is guaranteed that the output is dynamically consistent, otherwise an error is thrown

Boxes is assumed to be a vector of boxes along an edge.

Throws ErrorException in error cases. TODO: maybe fix this
"""
function getInclusionEdgeBoxIndices(boxes::Vector{Vector{Int}}, cplxPointsVecSorted, cplxPointsSortPerm=nothing; filter_missing=false, shortest_path_fix=false, safe_output=true)
    boxesInd = getBoxIndexRanges(boxes, cplxPointsVecSorted)
    isMissingFromCplx = isempty.(boxesInd)
    boxesIndFiltered = first.(boxesInd[.!isMissingFromCplx])

    if any(isMissingFromCplx) && !filter_missing
        a_missing_box = boxes[findfirst(isMissingFromCplx)]
        throw(ErrorException("Box $a_missing_box missing."))
    end

    boxesFiltered = boxes[.!isMissingFromCplx]

    if countDynamicInconsistencies(boxesFiltered) > 0
        if shortest_path_fix
            fixedBoxes = fixInclusionWithShortestPaths(boxesFiltered, cplxPointsVecSorted)
            boxesRngFixed = getBoxIndexRanges(fixedBoxes, cplxPointsVecSorted)
            boxesIndFiltered = first.(boxesRngFixed)

            if safe_output && countDynamicInconsistencies(fixedBoxes) > 0
                throw(ErrorException("Unskippable boxes in the edge."))
            end
        elseif safe_output
            # This will pop up until I use the Nerve
            throw(ErrorException("Unskippable boxes in the edge."))
        end
    end
    # due to the geometry of the projection of the l_2 onto the l_infty sphere,
    # a path may hit the same box multiple times.
    # removing missing boxes may lead to the same box appearing several times in succession.
    boxesIndFiltered = removeSuccessiveDuplicates(boxesIndFiltered)

    if isnothing(cplxPointsSortPerm)
        return boxesIndFiltered
    else
        return cplxPointsSortPerm[boxesIndFiltered]
    end
end

function getBoxIndexRanges(boxes, cplxPointsVecSorted)
    return map(boxes) do b
        return searchsorted(cplxPointsVecSorted, b)
    end
end

###
### Shortest Path Fix
### This Heuristic is needed since the cubes only cover the data points while a thickening of the cubes
### (which still has the same homotopy type) also covers a thickening of the data.
### In order to construct the inclusion map we then sometimes have to fix the map.
###

function fixInclusionWithShortestPaths(boxes, cplxPointsVecSorted)
    boxVecs = Vector{Vector{Int}}()
    for i = 1:length(boxes)-1
        curBoxVecs = shortestPathInclusion(boxes[i], boxes[i+1], cplxPointsVecSorted)
        if i == 1
            push!(boxVecs, curBoxVecs...)
        else
            push!(boxVecs, curBoxVecs[2:end]...)
        end
    end
    return boxVecs
end

"""
    function shortestPathInclusion(box1, box2, cplxPointsVecSorted)

    Given two boxes, returns a shortest path between them using the boxes in cplxPointsVecSorted.
"""
function shortestPathInclusion(box1, box2, cplxPointsVecSorted)
    boxDist = chebyshev(box1, box2)
    if boxDist == 0
        return [box1]
    elseif boxDist == 1
        return push!([box1], box2)
    else
        boxesFiltered = filter(v -> chebyshev(v, box1) <= boxDist + 2 && chebyshev(v, box2) <= boxDist + 2, cplxPointsVecSorted)
        return bfsShortestPath(box1, box2, boxesFiltered)
    end
end

"""
    function bfsShortestPath(box1, box2, boxesSorted)

Finds a shortets path between box1 and box2 in boxesSorted. Requires that box1 and box2 are in boxesSorted.
Returns a vector of box indices.
"""
function bfsShortestPath(box1, box2, boxesSorted)
    # TODO: check if it is necessary to choose a special shortest path (one with small l1 variance)
    if box1 == box2
        return [box1]
    end
    predecessor = zeros(Int, length(boxesSorted))
    ind1 = searchsorted(boxesSorted, box1)
    if length(ind1) != 1
        error("Box has to be unique and in the points vec. It is at position: $ind1")
    end
    ind1 = ind1[1]
    queue = [ind1]
    predecessor[ind1] = -1

    while !isempty(queue)
        i = popfirst!(queue)
        cur_box = boxesSorted[i]
        for (j, box) in enumerate(boxesSorted)
            if chebyshev(cur_box, box) <= 1 && predecessor[j] == 0
                predecessor[j] = i
                push!(queue, j)

                if box == box2
                    # if we are here, we have reached the box we want
                    # backtrack and return vector of indices
                    retVal = [box]
                    i = predecessor[j]
                    while i != -1
                        cur_box = boxesSorted[i]
                        push!(retVal, cur_box)
                        i = predecessor[i]
                    end
                    return reverse(retVal)
                end
            end
        end
    end
    error("Could not reach")
end
