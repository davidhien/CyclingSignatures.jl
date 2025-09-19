struct TrajectoryBar{T<:Integer,N}
    # helper struct with detailed type information for persistence bar
    birth::Float64
    death::Float64
    simplex_list::Vector{NTuple{N,Int}}
    coeff_list::Vector{T}
end

function TrajectoryBar(birth, death, simplex_list, coeff_list)
    return TrajectoryBar(Float64(birth), Float64(death), simplex_list, coeff_list)
end

"""
    trajectory_barcode(::Val{:DistanceMatrix}, points, metric, fltThreshold, field=DEFAULT_FIELD)

Computes a trajectory barcode using the distance matrix method, for arguments see [`trajectory_barcode`](@ref).
This assumes (and formall makes use of) a curve hypothesis.

The distance complex of a the points in `points` with respect to the `metric` is the filtered complex where
- vertices are pairs `(i, j)` with `i < j`, filtered by `metric(points[:, i], points[:, j])` and
- edges are nearest 4, with filtration value of an edge being the max of its adjacent vertices.

# Note:
In general, the returned collection of bars is not a persistence diagram for the point cloud.
"""
function trajectory_barcode(::Val{:DistanceMatrix}, points, metric, flt_threshold, field=DEFAULT_FIELD)
    min_vertices = dm_components_implicit(points, metric, flt_threshold)

    # generatePersistence diagram
    return map(node->vertex_to_essential_bar(node, points, metric, field=field), min_vertices)
    #return persistenceDiagramFromNodes(smallestNodes, points, metric, field=field)
end

function trajectory_barcode(::Val{:DistanceMatrixOld}, points, metric, flt_threshold, field=DEFAULT_FIELD)
    min_vertices = dm_components_explicit(points, metric, flt_threshold)

    # generatePersistence diagram
    return map(node->vertex_to_essential_bar(node, points, metric, field=field), min_vertices)
end

"""
    vertex_to_essential_bar(node, points, metric; field=DEFAULT_FIELD)

Given a node `(i,j)` in the distance matrix, returns an essential bar corresponding to the induced cycle in the Vietoris--Rips complex.
"""
function vertex_to_essential_bar(node, points, metric; field=DEFAULT_FIELD)
    i,j = node[1],node[2]
    b = metric(points[:,i],points[:,j])
    d = Inf
    edges, coeffs = curve_cycle(i,j; F=field)
    return TrajectoryBar(b,d,edges,coeffs)
end

"""
    curve_cycle(i,j;F=DEFAULT_FIELD)

Returns the curve cycle corresponding to `i` and `j` where `i<j`.

# Returns
- `edges`: vector of edges `[i,i+1], [i+1,i+2], ..., [j-1,j], [i,j]`.
- `coeffs`: coefficients `[1,...,1,-1]`.
"""
function curve_cycle(i,j;F=DEFAULT_FIELD)
    a, b = min(i,j),max(i,j)

    edges = map(a:(b-1)) do i
        return (i,i+1)
    end
    push!(edges, (a,b))
    coeffs = F.([ones(Int,length(a:b)-1);-1])
    return edges, coeffs
end

"""
    dm_components_explicit(points, metric, flt_threshold)

Computes a filtration-minimal vertex of each connected component of the distance complex up to the filtration threshold `flt_threshold`.
The distance complex is specified by `points` and `metric`.

This is implmented using a two pass annotation algorithm for images, the image being the upper right part of the distance matrix.
In this implementation, the distance matrix is computed explicitly.
"""
function dm_components_explicit(points, metric, flt_threshold)
    # implmented using a two pass distance matrix algorithm
    # First pass: provisional labels are assigned to each entry, a union-structure for the labels is maintained
    # Second pass: union-find structure is used to find the smallest node in each connected component
    cc, cc_labels, d_mat = dm_components_first_pass_explicit(points, metric, flt_threshold)
    min_vertex_dict = dm_components_second_pass_explicit(cc, d_mat,cc_labels)

    min_vertices = component_representatives(min_vertex_dict, cc_labels)

    return collect(values(min_vertices))
end

"""
    filtration_min_vertices_implicit_dm(points, metric, flt_threshold)

Computes a filtration-minimal vertex of each connected component of the distance complex specified by `points` and `metric` up to the filtration threshold `flt_threshold`.
This is implmented using a two pass annotation algorithm for images, the image being the upper right part of the distance matrix.
In this implementation, the distance matrix is *not* computed explicitly.
"""
function dm_components_implicit(points, metric, flt_threshold)
    # implmented using a two pass distance matrix algorithm
    # First pass: provisional labels are assigned to each entry, a union-structure for the labels is maintained
    # Second pass: union-find structure is used to find the smallest node in each connected component
    cc, cc_labels = dm_components_first_pass_implicit(points, metric, flt_threshold)
    min_vertex_dict = dm_components_second_pass_implicit(cc, cc_labels, points, metric)

    min_vertices = component_representatives(min_vertex_dict, cc_labels)

    return collect(values(min_vertices))
end

function component_representatives(min_vertex_dict, cc_labels)
    diag_label = cc_labels[1,1]
    keys_filtered = filter!(!=(diag_label), collect(keys(min_vertex_dict)))

    return map(k-> min_vertex_dict[k], keys_filtered)
end

function dm_components_second_pass_explicit(cc::IntDisjointSets{Int}, d_mat,cc_labels::Matrix{Int})
    smallestNodes = Dict{Int,Tuple{Int,Int}}()
    _,n = size(d_mat)
    @inbounds for j = 1:n
        for i = 1:j
            if cc_labels[i,j] > 0
                tmp_label = cc_labels[i,j]
                label = find_root!(cc,tmp_label)
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

function dm_components_second_pass_implicit(cc::IntDisjointSets{Int}, cc_labels::Matrix{Int}, points, metric)
    smallestNodes = Dict{Int,Tuple{Int,Int}}()
    smallestNodeVal = Dict{Int,Float64}()
    n = size(points,2)
    @inbounds for j = 1:n
        for i = 1:j
            if cc_labels[i,j] > 0
                tmp_label = cc_labels[i,j]
                label = find_root!(cc,tmp_label)
                cur_metric = metric((@view points[:,i]),(@view points[:,j]))
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
    dm_components_first_pass_explicit(points, metric, threshold)

Given a matrix with trajectory points, a metric and a filtrationThreshold, computes the connected components of the
graph defined by the matrix `pairwise(metric,points) .<= threshold` where
- vertices are matrix entries,
- edges are nearest 4.

Returns a disjoint set with labels, the label of each entry and the distance matrix.
Note: entries with different labels may be in the same component (iff their labels are in the same set in cc)
"""
function dm_components_first_pass_explicit(points, metric, threshold)
    cc = IntDisjointSets(0)
    d_mat = pairwise(metric, points, dims=2)
    cc_labels = zeros(Int, size(d_mat))

    n,_ = size(d_mat)
    @inbounds for j = 1:n
        for i = 1:j
            # iterate over the upper right part of d_mat
            if d_mat[i,j] <= threshold
                assign_initial_label!(cc_labels, cc, i, j)
            end
        end
    end
    return cc, cc_labels, d_mat
end

"""
    dm_components_first_pass_explicit(points, metric, threshold)

Given a matrix with trajectory points, a metric and a filtration threshold, computes the connected components of the
graph defined by the matrix `pairwise(metric,points) .<= threshold` where
- vertices are matrix entries,
- edges are nearest 4.

Returns a disjoint set with label and the label of each entry.
Note: entries with different labels may be in the same component (iff their labels are in the same set in cc)
"""
function dm_components_first_pass_implicit(points, metric, threshold)
    cc = IntDisjointSets(0)
    n = size(points,2)
    cc_labels = zeros(Int, n,n)
    @inbounds for j = 1:n
        # Note: using views gives significant speed-up
        for i = 1:j
            # iterate over the upper right part of d_mat
            if metric((@view points[:,j]),(@view points[:,i])) <= threshold
                assign_initial_label!(cc_labels, cc, i, j)
            end
        end
    end
    return cc, cc_labels
end

function assign_initial_label!(cc_labels, cc, i, j)
    west = j>i && cc_labels[i,j-1] != 0
    north = i>1 && cc_labels[i-1,j] != 0
    if west && north
        w_label = cc_labels[i,j-1]
        n_label = cc_labels[i-1,j]
        if w_label < n_label
            cc_labels[i,j] = w_label
            union!(cc, w_label, n_label)
        else
            cc_labels[i,j] = n_label
            union!(cc, n_label, w_label)
        end
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
