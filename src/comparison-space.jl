# ==============================================================================
# CUBICAL ACYCLIC CARRIER INTERFACE
# ==============================================================================

"""
Abstract base type for cubical acyclic carriers.
Subtypes represent an implementation of a cell complex representing a union of cubes which
allows the construction of chains from acyclic subsets and the annotation of cycles using a cohomology basis.

Subtypes must implement:
- `induced_one_chain(carrier, edge_boxes)`: Given a vector of boxes covering an edge, computes the induced one-chain.
- `annotate_chain(carrier, chain)`: Annotates a chain using a cohomology basis. For cycles, the result is independent of the representative.
- `betti_1(carrier)`: Returns the first betti number of the carrier.
"""
abstract type AbstractCubicalAcyclicCarrier end

"""
    induced_one_chain(carrier::AbstractCubicalAcyclicCarrier, edge_boxes)

Implements the conversion of an acyclic cubical cover to an induced chain.

# Arguments
- `carrier::AbstractCubicalAcyclicCarrier`: The cubical acyclcic carrier
- `edge_boxes`: Boxes that cover an one chain

# Returns
A representation of an induced chain
"""
function induced_one_chain(carrier::AbstractCubicalAcyclicCarrier, edge_boxes)
    error("induced_one_cycle not implemented for $(typeof(carrier))")
end

function annotate_chain(carrier::AbstractCubicalAcyclicCarrier, chain)
    error("annotate_chain not implemented for $(typeof(carrier))")
end

"""
    betti_1(carrier)

Returns the first betti number of the carrier.
"""
function betti_1(carrier::AbstractCubicalAcyclicCarrier)
    error("betti_1 not implemented for $(typeof(carrier))")
end

# ==============================================================================
# CUBICAL ACYCLIC CARRIER IMPLEMENTATIONS
# ==============================================================================

"""

# Fields
- `pts::Matrix{S}`: Points in the carrier, where `S` is an integer type. The columns are assumed to be sorted.
- `cplx::CellComplex{Simplex}`: The Nerve of the cubical cover, it is the ``l_\\infty``-distance VR-complex of the box centers.
- `h1::V`: The transpose of h1 is a basis for the first cohomology of the complex.
"""
struct CubicalVRCarrier{S<:Integer,V} <: AbstractCubicalAcyclicCarrier
    pts::Matrix{S}
    cplx::CellComplex{Simplex}
    h1::V
end

function Base.show(io::IO, carrier::CubicalVRCarrier)
    println(io, "CubicalVRCarrier")
    println(io, "  Number of points: ", size(carrier.pts, 2))
    println(io, "  Number of edges:  ", length(carrier.cplx.cells[1]))
    println(io, "  β_1:              ", size(carrier.h1, 1))
end


points(C::CubicalVRCarrier) = C.pts
complex(C::CubicalVRCarrier) = C.cplx
h1(C::CubicalVRCarrier) = C.h1

"""
    CubicalVRCarrier(pts::Matrix{S}) where S

TODO: write docstring
"""
function CubicalVRCarrier(pts::AbstractMatrix{S}) where S<:Integer
    if !issorted(eachcol(pts))
        pts = pts[:, sortperm(eachcol(pts))]
    end

    cplx = vr_incremental(pts, chebyshev, 1)
    @info "Complex built with $(length(get(cplx.cells, 1, []))) edges and $(length(get(cplx.cells, 2, []))) triangles"
    @info "Generating boundary matrices..."
    D0 = coboundaryMatrix(cplx, 0)
    D1 = coboundaryMatrix(cplx, 1)
    @info "Computing cohomology..."
    h1 = copy(transpose(firstCohomology(D0, D1)))
    @info "Cohomology computation complete."
    deleteAllBoundaryMatrices(cplx) # can be recomputed if necessary
    @info "CubicalVRCarrier created with $(size(pts, 2)) points and $(size(h1,1)) cohomology basis vectors."
    size(h1,1) == 0 && @warn "The carrier has trivial first cohomology."
    return CubicalVRCarrier{S,typeof(h1)}(pts, cplx, h1)
end

function induced_one_chain(carrier::CubicalVRCarrier, cover_boxes)
    indices = Int[]
    filter!(cover_boxes) do box
        ind = searchsorted(eachcol(carrier.pts), box)
        if length(ind) == 1
            push!(indices, ind[1])
        end
        return !isempty(ind)
    end

    # remove successive duplicates
    no_duplicate_indices = filter(i -> i == 1 || indices[i] != indices[i-1], 1:length(indices))
    cover_boxes = cover_boxes[no_duplicate_indices]
    indices = indices[no_duplicate_indices]

    # generate chain
    v = spzeros(Int, length(carrier.cplx.cells[1]))
    for i in 1:length(cover_boxes)-1
        if chebyshev(cover_boxes[i], cover_boxes[i+1]) == 1
            # find edge in cplx
            e = Simplex([indices[i], indices[i+1]])
            e_ind = searchsortedfirst(carrier.cplx.cells[1], e)
            v[e_ind] += sign(indices[i+1] - indices[i])
        else
            fix_indices = fix_edge(carrier, cover_boxes[i], cover_boxes[i+1])

            for j in 1:length(fix_indices)-1
                e = Simplex([fix_indices[j], fix_indices[j+1]])
                e_ind = searchsortedfirst(carrier.cplx.cells[1], e)
                v[e_ind] += sign(fix_indices[j+1] - fix_indices[j])
            end
        end
    end
    return v
end

"""
    fix_edge(carrier::CubicalVRCarrier, box1, box2)

Searches for a (short) shortest path between the boxes.
Note that current implementation is not well-defined since it may be random which way a missing box is bypassed.
"""
function fix_edge(carrier::CubicalVRCarrier, box1, box2)
    dist = chebyshev(box1, box2)
    if dist == 0
        return [box1]
    elseif dist == 1
        return push!([box1], box2)
    else
        boxes_filtered = filter(v -> max(chebyshev(v, box1), chebyshev(v, box2)) <= dist + 2, eachcol(carrier.pts))
        edge_boxes = boxes_bfs_shortest_path(boxes_filtered, box1, box2)
        return map(b -> searchsortedfirst(eachcol(carrier.pts), b), edge_boxes)
    end
end

"""
    grid_bfs_shortest_path(points, p1, p2)

Given a sorted matrix of points `points` on an integer grid and two points `p1` and `p2`, finds a shortest path between them using BFS.
"""
function boxes_bfs_shortest_path(points, p1, p2)
    if p1 == p2
        return [p1]
    end
    predecessor = zeros(Int, length(points))
    ind1 = searchsorted(points, p1)
    if length(ind1) != 1
        error("Box has to be unique and in points. It is (not) at position: $ind1")
    end
    ind1 = ind1[1]
    queue = [ind1]
    predecessor[ind1] = -1

    while !isempty(queue)
        i = popfirst!(queue)
        cur_pt = points[i]
        for (j, pt) in enumerate(points)
            if chebyshev(cur_pt, pt) <= 1 && predecessor[j] == 0
                predecessor[j] = i
                push!(queue, j)

                if pt == p2
                    # if we are here, we have reached the box we want
                    # backtrack and return vector of indices
                    pt_path = [pt]
                    i = predecessor[j]
                    while i != -1
                        cur_box = points[i]
                        push!(pt_path, cur_box)
                        i = predecessor[i]
                    end
                    return reverse(pt_path)
                end
            end
        end
    end
    error("Could not reach")
end

function annotate_chain(carrier::CubicalVRCarrier, chain)
    return carrier.h1 * chain
end

function betti_1(carrier::CubicalVRCarrier)
    return size(carrier.h1, 1)
end

struct CubicalCarrier{S<:Integer,C,V} <: AbstractCubicalAcyclicCarrier
    pts::Matrix{S}
    cplx::CellComplex{C}
    h1::V
end

# Interface methods for AbstractCubicalAcyclicCarrier
# TODO: Add required interface methods here:
# get_complex(carrier::AbstractCubicalAcyclicCarrier) = error("Not implemented")
# get_h1(carrier::AbstractCubicalAcyclicCarrier) = error("Not implemented")

# ==============================================================================
# COMPARISON SPACE INTERFACE
# ==============================================================================

"""
Abstract base type for comparison spaces.

Subtypes must implement:
- `map_cycle(cs, points, simplices, coeffs)`: Maps a cycle into the comparison space
- `betti_1(cs)`
"""
abstract type AbstractComparisonSpace end

"""
    map_cycle(cs::AbstractComparisonSpace, points, simplices, coeffs)

Annotates the cycle specified by points, simplices and coeffs into a comparison
space using the given comparison space.

# Arguments
- `cs::AbstractComparisonSpace`: The comparison space
- `points`: Points defining the cycle
- `simplices`: Simplices in the cycle
- `coeffs`: Coefficients for the cycle

# Returns
Mapped cycle representation (implementation-dependent)
"""
function map_cycle(cs::AbstractComparisonSpace, points, simplices, coeffs)
    error("map_cycle not implemented for $(typeof(cs))") # TODO: change this and make it return the proper exception type
end

"""
    betti_1(cs::AbstractComparisonSpace)

Returns the first betti number of the comparison space.
"""
betti_1(cs::AbstractComparisonSpace) = error("betti_1 not implemente for $(typeof(cs))")

# Specialized abstract type for cubical comparison spaces
"""
Abstract type for cubical-based comparison spaces.
Inherits from AbstractComparisonSpace.

Subtypes must implement:
- `edge_boxes(comparison_space, p1, p2)`: Computes the boxes covering the edge between points `p1` and `p2`.
- `carrier(comparison_space)`: Returns an associated cubical acyclic carrier.
- `betti_1(comparison_space)`: Returns the first betti number of the comparison space.
"""
abstract type AbstractCubicalComparisonSpace <: AbstractComparisonSpace end

function map_cycle(cs::AbstractCubicalComparisonSpace, points::AbstractMatrix, simplices, coeffs)
    # TODO: save allocations by adding to this chain only
    v = spzeros(eltype(coeffs), length(carrier(cs).cplx.cells[1]))

    for (a, simplex) in zip(coeffs, simplices)
        v += a * map_simplex_into_cubes(cs, points, simplex)
    end

    return annotate_chain(carrier(cs), v)
end

function map_simplex_into_cubes(cover::AbstractCubicalComparisonSpace, points, simplex)
    if length(simplex) == 2
        return map_edge_to_cubes(cover, @view(points[:, simplex[1]]), @view(points[:, simplex[2]]))
    end
    error("Not supported yet.")
end

function map_edge_to_cubes(cover::AbstractCubicalComparisonSpace, p1, p2)
    eb = edge_boxes(cover, p1, p2)
    return induced_one_chain(carrier(cover), eb)
end

# ==============================================================================
# CONCRETE COMPARISON SPACE IMPLEMENTATIONS
# ==============================================================================

"""
Standard cubical comparison space.

# Fields
- `boxsize::T`: Size of cubical boxes
- `carrier::C`: Associated cubical acyclic carrier
"""
struct CubicalComparisonSpace{C<:AbstractCubicalAcyclicCarrier,T} <: AbstractCubicalComparisonSpace
    carrier::C
    boxsize::T
end

function Base.show(io::IO, space::CubicalComparisonSpace)
    println(io, "CubicalComparisonSpace (boxsize $(space.boxsize), β_1: $(betti_1(space)))")
    print(io, "  Carrier: ")
    show(IOContext(io, :compact => true), space.carrier)
end


function CubicalComparisonSpace(box_centers::AbstractMatrix{P}, boxsize) where P<:Integer
    if boxsize <= 0
        error("Box size must be positive.")
    end

    carrier = CubicalVRCarrier(box_centers)
    return CubicalComparisonSpace{typeof(carrier),typeof(boxsize)}(carrier, boxsize)
end

function cubical_vr_comparison_space_via_cover(pts::AbstractMatrix, boxsize)
    if boxsize <= 0
        error("Box size must be positive.")
    end

    quantized_pts = round.(Int, pts ./ boxsize)
    comp_space_pts = sortslices(unique(quantized_pts, dims=2), dims=2)

    return CubicalComparisonSpace(comp_space_pts, boxsize)
end


"""
Sphere bundle cubical comparison space.

# Fields
- `boxsize::T`: Size of cubical boxes
- `sb_radius::TT`: Sub-block radius parameter
- `carrier::C`: Associated cubical acyclic carrier
"""
struct SBCubicalComparisonSpace{C<:AbstractCubicalAcyclicCarrier,T,TT} <: AbstractCubicalComparisonSpace
    carrier::C
    boxsize::T
    sb_radius::TT
end

function Base.show(io::IO, space::SBCubicalComparisonSpace)
    println(io, "SBCubicalComparisonSpace (boxsize: $(space.boxsize), sb_radius: $(space.sb_radius), β_1: $(betti_1(space)))")
    print(io, "  Carrier: ")
    show(IOContext(io, :compact => true), space.carrier)
end


function SBCubicalComparisonSpace(box_centers::AbstractMatrix{P}, boxsize, sb_radius) where P<:Integer
    if boxsize <= 0
        error("Box size must be positive.")
    elseif sb_radius <= 0
        error("Sphere bundle radius must be positive.")
    end

    carrier = CubicalVRCarrier(box_centers)
    return SBCubicalComparisonSpace{typeof(carrier),typeof(boxsize),typeof(sb_radius)}(carrier, boxsize, sb_radius)
end

function sb_cubical_vr_comparison_space_via_cover(pts::AbstractMatrix, boxsize, sb_radius)
    if boxsize <= 0
        error("Box size must be positive.")
    elseif sb_radius <= 0
        error("Sphere bundle radius must be positive.")
    elseif isodd(size(pts, 1))
        error("sphere bundle points have to be even-dimensional.")
    end

    d = div(size(pts, 1), 2)

    quantized_pts = Matrix{Int}(undef, size(pts))
    quantized_pts[1:d, :] = round.(Int, pts[1:d, :] ./ boxsize)
    quantized_pts[d+1:2*d, :] = mapslices(v -> round.(Int, normalize(v, Inf) * sb_radius), pts[d+1:2*d, :], dims=1)
    comp_space_pts = sortslices(unique(quantized_pts, dims=2), dims=2)

    return SBCubicalComparisonSpace(comp_space_pts, boxsize, sb_radius)
end


# ==============================================================================
# INTERFACE IMPLEMENTATIONS
# ==============================================================================

carrier(cover::CubicalComparisonSpace) = cover.carrier
carrier(cover::SBCubicalComparisonSpace) = cover.carrier

betti_1(cover::CubicalComparisonSpace) = betti_1(cover.carrier)
betti_1(cover::SBCubicalComparisonSpace) = betti_1(cover.carrier)

function edge_boxes(cover::CubicalComparisonSpace, p1, p2)
    q1 = p1 / cover.boxsize
    q2 = p2 / cover.boxsize
    box1 = round.(Int, q1)
    box2 = round.(Int, q2)

    if box1 == box2
        return [box1]
    elseif chebyshev(box1, box2) == 1
        return [box1, box2]
    else
        return edge_boxes(q1, q2)
    end
end

function edge_boxes(cover::SBCubicalComparisonSpace, x1, x2)
    d = div(length(x1), 2)
    b1 = Vector{Int}(undef, 2d)
    b2 = Vector{Int}(undef, 2d)
    utb_normalizer_1 = cover.sb_radius / norm(x1[d+1:end],Inf)
    utb_normalizer_2 = cover.sb_radius / norm(x2[d+1:end],Inf)
    @inbounds for i in 1:d
        b1[i] = round(Int, x1[i] / cover.boxsize)
        b2[i] = round(Int, x2[i] / cover.boxsize)
        b1[d+i] = round(Int, utb_normalizer_1 * x1[d+i])
        b2[d+i] = round(Int, utb_normalizer_2 * x2[d+i])
    end

    # Compute distances
    dist = maximum(i -> abs(b1[i] - b2[i]), 1:2*d)

    if dist == 0
        return [b1]
    elseif dist == 1
        return [b1, b2]
    end

    return edge_boxes_sphere_bundle(x1[1:d], x2[1:d], normalize(x1[d+1:end], Inf), normalize(x2[d+1:end], Inf), cover.sb_radius)
end

# ==============================================================================
# BOX COVER METHODS
# ==============================================================================

"""
    edge_boxes(p1, p2)

Computes the centers of all integer grid boxes which intersect the edge `[p1;p2]`.
"""
function edge_boxes(p1, p2)
    n = length(p1)
    if n != length(p2)
        error("p1 and p2 have to be of the same dimension.")
    end
    b1 = round.(Int, p1)
    b2 = round.(Int, p2)

    if b1 == b2
        return [b1]
    end

    # setup variables
    t_next_crossing = zeros(length(p1)) # t of the next crossing in the component
    t_between_crossings = zeros(length(p1))
    dir = p2 - p1
    dir_sign = Int.(sign.(dir))

    for i in 1:n
        if b1[i] != b2[i]
            t_next_crossing[i] = (2 * (b1[i] - p1[i]) + dir_sign[i]) / (2 * (dir[i]))
            t_between_crossings[i] = 1 / abs(dir[i])
        else
            t_next_crossing[i] = Inf
            t_between_crossings[i] = Inf
        end
    end

    # main loop
    boxes = Vector{Int}[copy(b1)]
    bnext = b1
    while bnext != b2
        i = argmin(t_next_crossing)
        bnext[i] += dir_sign[i]
        push!(boxes, copy(bnext))
        # update t_next_crossing
        t_next_crossing[i] += t_between_crossings[i]
    end
    return boxes
end

function edge_boxes(p1, p2, boxsize)
    error("This is deprecated due to inconsistencies")
    if boxsize != 1
        p1 /= boxsize
        p2 /= boxsize
    end
    return edge_boxes(p1, p2)
end

function edge_boxes_sphere_bundle(p1, p2, v1, v2, sb_radius)
    lambdas = projection_max_switches(sb_radius * v1, sb_radius * v2)

    # in each interval [lambdas[i],lambdas[i+1]], the edge is projected to the same face
    #   => normal edge_boxes can be used
    boxes = Vector{Int}[]
    x1 = 1.0 * [p1; v1]
    x2 = 1.0 * [p1; sb_radius * normalize(v1, Inf)]
    for i in 2:length(lambdas)
        l = lambdas[i]
        x1 .= x2
        x2 .= [(1 - l) * p1 + l * p2; sb_radius * normalize((1 - l) * v1 + l * v2, Inf)]
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
        v = (1 - t_swap_min) * v1 + t_swap_min * v2
        max_i = projection_initial_max_component(v, v2)
        # reset
        push!(Ts, t_swap_min)
        cur_t = t_swap_min
        t_swap_min = 1.0
    end

    return Ts
end

function projection_initial_max_component(v1, v2)
    max_v1 = maximum(abs, v1)
    is = findall(x -> isapprox(abs(x), max_v1), v1)
    return argmax(i -> (sign(v1[i]) >= 0 ? 1 : -1) * v2[i], is)
end

"""
    intersect_lines(a,b)

Computes the intersection of the lines `(1-l)*a[i]+l*b[i]`
"""
function intersect_lines(a, b)
    return intersect_lines(a[1,], a[2], b[1], b[2]) # note condition not great, maybe use geometrybasics implementation
end

function intersect_lines(a1, a2, b1, b2)
    return (a1 - a2) / (b2 - b1 - a2 + a1) # note condition not great, maybe use geometrybasics implementation
end
