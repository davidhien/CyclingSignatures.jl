function trajectoryBarcode(::Val{:DistanceMatrixOld}, trajPoints, metric, fltThreshold, field=DEFAULT_FIELD)
    # only compute representatives for infinite bars
    # TODO: check if the following significantly improve performance
    # - only one single allocation of d_mat and cc_labels
    # - implicit use of the distance matrix (by ad hoc computing the relevant distances)
    # - adding inlines and other performance enhancing drugs
    smallestNodes = twoPassDMAnnotation(trajPoints, metric, fltThreshold)    

    # generatePersistence diagram
    return trajectoryBarsFromNodes(smallestNodes, trajPoints, metric, field=field)
end

function trajectoryBarcode(::Val{:DistanceMatrix}, trajPoints, metric, fltThreshold, field=DEFAULT_FIELD)
    # only compute representatives for infinite bars
    # TODO: check if the following significantly improve performance
    # - only one single allocation of d_mat and cc_labels
    # - implicit use of the distance matrix (by ad hoc computing the relevant distances)
    # - adding inlines and other performance enhancing drugs
    smallestNodes = twoPassDMAnnotation2(trajPoints, metric, fltThreshold)    

    # generatePersistence diagram
    return trajectoryBarsFromNodes(smallestNodes, trajPoints, metric, field=field)
    #return persistenceDiagramFromNodes(smallestNodes, trajPoints, metric, field=field)
end

function trajectoryBarsFromNodes(smallestNodes, trajPoints, metric; field=DEFAULT_FIELD)
    return map(smallestNodes) do v
        i,j = v[1],v[2]
        b = metric(trajPoints[:,i],trajPoints[:,j])
        d = Inf
        edges, coeffs = curveGeneratorData(i,j; F=field)
        return TrajectoryBar(b,d,edges,coeffs)
    end
end

function persistenceDiagramFromNodes(smallestNodes, trajPoints, metric; field=DEFAULT_FIELD)
    intervals = map(smallestNodes) do v
        i,j = v[1],v[2]
        b = metric(trajPoints[:,i],trajPoints[:,j])
        d = Inf
        edges, coeffs = curveGeneratorData(i,j; F=field)
        return PersistenceInterval(b,d,simplex_list=edges, coeff_list=coeffs)
    end
    return PersistenceDiagram(intervals)
end

function curveGeneratorData(i,j;F=DEFAULT_FIELD)
    a, b = min(i,j),max(i,j)

    edges = map(a:(b-1)) do i
        return (i,i+1)
    end
    push!(edges, (a,b))
    coeffs = F.([ones(Int,length(a:b)-1);-1])
    return edges, coeffs
end

function twoPassDMAnnotation(trajPoints, metric, fltThreshold)
    cc, d_mat,cc_labels = twoPassDMAnnotationFirstPass(trajPoints, metric, fltThreshold)
    smallestNodesDict = twoPassDMAnnotationSecondPass(cc, d_mat,cc_labels)
    smallestNodes = connectedComponentRepresentatives(smallestNodesDict, cc_labels)
    return collect(values(smallestNodes))
end

function twoPassDMAnnotation2(trajPoints, metric, fltThreshold)
    cc, cc_labels = twoPassDMAnnotationFirstPass2(trajPoints, metric, fltThreshold)
    smallestNodesDict = twoPassDMAnnotationSecondPass2(cc, cc_labels, trajPoints, metric)
    smallestNodes = connectedComponentRepresentatives(smallestNodesDict, cc_labels)
    return collect(values(smallestNodes))
end

function connectedComponentRepresentatives(smallestNodesDict, cc_labels)
    diag_label = cc_labels[1,1]
    keys_filtered = filter!(i-> i != diag_label,collect(keys(smallestNodesDict)))

    return map(k-> smallestNodesDict[k], keys_filtered)
end

function twoPassDMAnnotationSecondPass(cc::IntDisjointSets{Int}, d_mat,cc_labels::Matrix{Int})
    smallestNodes = Dict{Int,Tuple{Int,Int}}()
    _,n = size(d_mat)
    @inbounds for j = 1:n
        for i = 1:j
            if cc_labels[i,j] > 0
                tmp_label = cc_labels[i,j]
                label = find_root(cc,tmp_label)
                # now fix representant of component label
                if !haskey(smallestNodes, label)
                    smallestNodes[label] = (i,j)
                elseif d_mat[i,j] <d_mat[smallestNodes[label]...]
                    smallestNodes[label] = (i,j)
                end
            end
        end
    end
    return smallestNodes
end

function twoPassDMAnnotationSecondPass2(cc::IntDisjointSets{Int}, cc_labels::Matrix{Int}, trajPoints, metric)
    smallestNodes = Dict{Int,Tuple{Int,Int}}()
    smallestNodeVal = Dict{Int,Float64}()
    n = size(trajPoints,2)
    @inbounds for j = 1:n
        for i = 1:j
            if cc_labels[i,j] > 0
                tmp_label = cc_labels[i,j]
                label = find_root(cc,tmp_label)
                cur_metric = metric((@view trajPoints[:,i]),(@view trajPoints[:,j]))
                # now fix representant of component label
                if !haskey(smallestNodes, label)
                    smallestNodes[label] = (i,j)
                    smallestNodeVal[label] = cur_metric
                elseif cur_metric < smallestNodeVal[label]
                    smallestNodes[label] = (i,j)
                    smallestNodeVal[label] = cur_metric
                end
            end
        end
    end
    return smallestNodes
end


"""
    function twoPassDMAnnotationFirstPass(trajPoints, metric, fltThreshold)

Given a matrix with trajectory points, a metric and a filtrationThreshold, computes the connected components of the
graph defined by the matrix pairwise(metric,trajPoints) .<= fltThreshold where
- vertices are matrix entries,
- edges are nearest 4.

Returns a disjoint set with labels, the distance matrix and the label of each entry.
Note: entries with different labels may be in the same component (iff their labels are in the same set in cc) 
"""
function twoPassDMAnnotationFirstPass(trajPoints, metric, fltThreshold)
    cc = IntDisjointSets(0)
    d_mat = pairwise(metric, trajPoints)
    cc_labels = zeros(Int, size(d_mat))

    n,_ = size(d_mat)
    @inbounds for j = 1:n
        for i = 1:j
            # iterate over the upper right part of d_mat
            if d_mat[i,j] <= fltThreshold
                west = j>i && cc_labels[i,j-1] != 0
                north = i>1 && cc_labels[i-1,j] != 0 
                if west && north
                    w_label = cc_labels[i,j-1]
                    n_label = cc_labels[i-1,j]
                    labels_sorted = w_label<n_label ? [w_label;n_label] : [n_label;w_label]
                    cc_labels[i,j] = labels_sorted[1]
                    union!(cc, labels_sorted...)
                elseif north
                    cc_labels[i,j] = cc_labels[i-1,j]
                elseif west
                    cc_labels[i,j] = cc_labels[i,j-1]
                else
                    # new cc
                    new_label = push!(cc)
                    cc_labels[i,j] = new_label
                end
            end
        end
    end
    return cc, d_mat,cc_labels
end

function twoPassDMAnnotationFirstPass2(trajPoints, metric, fltThreshold)
    cc = IntDisjointSets(0)
    n = size(trajPoints,2)
    cc_labels = zeros(Int, n,n)
    @inbounds for j = 1:n
        # Note: using views gives significant speed-up
        for i = 1:j
            # iterate over the upper right part of d_mat
            if metric((@view trajPoints[:,j]),(@view trajPoints[:,i])) <= fltThreshold
                west = j>i && cc_labels[i,j-1] != 0
                north = i>1 && cc_labels[i-1,j] != 0 
                if west && north
                    w_label = cc_labels[i,j-1]
                    n_label = cc_labels[i-1,j]
                    labels_sorted = w_label<n_label ? [w_label;n_label] : [n_label;w_label]
                    cc_labels[i,j] = labels_sorted[1]
                    union!(cc, labels_sorted...)
                elseif north
                    cc_labels[i,j] = cc_labels[i-1,j]
                elseif west
                    cc_labels[i,j] = cc_labels[i,j-1]
                else
                    # new cc
                    new_label = push!(cc)
                    cc_labels[i,j] = new_label
                end
            end
        end
    end
    return cc,cc_labels
end
