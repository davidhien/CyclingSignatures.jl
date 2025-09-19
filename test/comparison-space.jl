using Test
using LinearAlgebra
using CyclingSignatures: edge_boxes, projection_initial_max_component, projection_max_switches
using CyclingSignatures: edge_boxes_sphere_bundle
using CyclingSignatures: boxes_bfs_shortest_path

using CyclingSignatures: CubicalVRCarrier, induced_one_chain, annotate_chain, betti_1
using CyclingSignatures: fix_edge
using CyclingSignatures.ATTools: boundaryMatrix

using CyclingSignatures: cubical_vr_comparison_space_via_cover, map_cycle
using CyclingSignatures: sb_cubical_vr_comparison_space_via_cover
using CyclingSignatures: FF
include("test-util.jl")

@testset "edge_boxes" begin
    p1 = [0; 0]
    p2 = [10; 0]

    res = [[i; 0] for i = 0:10]
    @test res == edge_boxes(p1, p2)

    p1 = [0; 0.25]
    p2 = [5; 5.25]
    res = [[div(i, 2); div(i + 1, 2)] for i = 0:10]
    @test res == edge_boxes(p1, p2)

    p1 = [4.2; 4.2]
    p2 = [4.2; 4.2]
    @test length(edge_boxes(p1, p2)) == 1
end

@testset "edge_boxes_utb" begin
    p1 = [0, 0]
    p2 = [0, 0]
    v1 = [1, 0]
    v2 = [1, 1]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 1)
    @test eb == [[0, 0, 1, 0], [0, 0, 1, 1]]

    for k = 1:10
        p1 = [0, 0]
        p2 = [0, 0]
        v1 = [1, 0]
        v2 = [0, 1]
        eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, k)
        res = Vector{Int}[]
        for i = 0:k
            push!(res, [0; 0; k; i])
        end
        for i = k-1:-1:0
            push!(res, [0; 0; i; k])
        end
        @test res == eb
    end

    p1 = [0, 0]
    p2 = [0, 0]
    v1 = [1, 0]
    v2 = [0, 1]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 2)
    @test eb == [[0, 0, 2, 0], [0, 0, 2, 1], [0, 0, 2, 2], [0, 0, 1, 2], [0, 0, 0, 2]]


    p1 = [0, 0]
    p2 = [10, 0]
    v1 = [1, 0]
    v2 = [1, 0]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 1)
    res = [[i; 0; 1; 0] for i = 0:10]
    @test res == eb
end

@testset "projection_initial_max_component" begin
    v1 = [0.8, 0.9, 0.9]
    v2 = [1.0, 0.6, 0.5]
    @test projection_initial_max_component(v1, v2) == 2

    v1 = [0.8, 0.9, 0.9]
    v2 = [1.0, -0.6, -0.5]
    @test projection_initial_max_component(v1, v2) == 3

    v1 = [0.8, 0.9, -0.9]
    v2 = [1.0, 0.6, -0.5]
    @test projection_initial_max_component(v1, v2) == 2

    v1 = [0.8, -0.9, 0.9]
    v2 = [1.0, 0.6, -0.5]
    @test projection_initial_max_component(v1, v2) == 3
end

@testset "projection_max_switches" begin
    @testset "explicit cases" begin
        # line a over line b
        v1 = [2; 1]
        v2 = [3; 2]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 1.0]; atol=1e-8)
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 1.0]; atol=1e-8)
        ls = projection_max_switches(v1, -v2)
        @test isapprox(ls, [0.0, 0.375, 0.5, 1.0]; atol=1e-8)
        ls = projection_max_switches(-v1, v2)
        @test isapprox(ls, [0.0, 0.375, 0.5, 1.0]; atol=1e-8)

        # line a intersects line b
        v1 = [2; 1]
        v2 = [1; 2]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])
        ls = projection_max_switches(v1, -v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])

        # line a over line b, b has irrelevant zero intersection
        v1 = [2; 1]
        v2 = [1; -0.5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 1.0])

        # line a positive, line b negative and takes max
        v1 = [2; 0]
        v2 = [4; -6]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])

        # like previous but flipped sign at second component
        v1 = [2; 0]
        v2 = [4; 6]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])

        # degenerate intersection on right endpoint
        v1 = [2; 0]
        v2 = [-1; 1]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])

        # start at same point, one in positive one in negative direction
        v1 = [1; 1]
        v2 = [3; -5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 0.5, 1.0])

        # two intersections
        v1 = [1; 3]
        v2 = [1; -5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 0.25, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 0.25, 0.5, 1.0])

        # intersection with constant zero function
        v1 = [1; 0]
        v2 = [1; 0]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 1.0])

        # multiple intersections
        # TODO!!!
    end
    @testset "random tests" begin
        for _ in 1:10
            v1 = rand(10)
            v2 = rand(10)
            ls = projection_max_switches(-v1, -v2)
            for l in ls[2:end-1]
                v = (1 - l) * v1 + l * v2
                m = maximum(abs, v)
                m_i = findall(isapprox(m), v)
                @test length(m_i) >= 2
            end
        end
    end
end

@testset "boxes_bfs_shortest_path" begin
    pts = [[0, 0], [1, 1], [2, 0]]
    boxes = boxes_bfs_shortest_path(pts, [0, 0], [2, 0])
    @test pts == boxes

    pts = [[0, 0], [1, 1], [2, 0]]
    boxes = boxes_bfs_shortest_path(pts, [2, 0], [0, 0])
    @test reverse(pts) == boxes
end

# ==============================================================================
# CUBICAL ACYCLIC CARRIER TESTS
# ==============================================================================

@testset "CubicalVRCarrier" begin
    @testset "basic functionality" begin
        # generate circle with holes
        pts_square = [collect(p) for p in Iterators.product(fill(-5:5, 2)...)]
        pts_unsorted = filter(!=([0, 0]), pts_square)
        pts_m_unsorted = hcat(pts_unsorted...)
        pts_sorted = pts_m_unsorted[:, sortperm(eachcol(pts_m_unsorted))]

        # test fields
        carrier = CubicalVRCarrier(pts_m_unsorted)
        @test size(carrier.h1, 1) == 1
        @test pts_sorted == carrier.pts

        ### setup for interface methods

        # generate cycle boxes
        boxes = sort(filter(p -> norm(p, Inf) == 2, pts_square), by=p -> atan(p[1], p[2])) # counter-clockwise order
        push!(boxes, boxes[1]) # close cycle

        c = induced_one_chain(carrier, boxes)

        # test counter-clockwise order
        edge_index = findfirst(!=(0), c)
        e_verts = carrier.cplx.cells[1][edge_index].vertices
        e_boxes = map(i -> pts_sorted[:, i], e_verts)
        e_boxes_ind = map(b -> findfirst(==(b), boxes), e_boxes)
        @test issorted(e_boxes_ind) == (sign(c[edge_index]) == 1)

        # test that its a chain
        D = boundaryMatrix(carrier.cplx, 1)
        @test all(==(0), D * c)

        # test annotation
        @test annotate_chain(carrier, c) == [1] || annotate_chain(carrier, c) == [-1]
        c2 = induced_one_chain(carrier, reverse(boxes))
        @test annotate_chain(carrier, c) == -annotate_chain(carrier, c2)

        # test betti_1
        @test betti_1(carrier) == 1
    end

    @testset "shortest path fix" begin
        pts_square = [collect(p) for p in Iterators.product(fill(-5:5, 2)...)]
        pts_v = sort(filter(v -> norm(v, Inf) >= 2 && v != [0, 2], pts_square)) # extra missing box
        pts_m = hcat(pts_v...)

        carrier = CubicalVRCarrier(pts_m)

        ind = fix_edge(carrier, [1, 2], [-1, 2])
        @test pts_v[ind[1]] == [1, 2]
        @test pts_v[ind[end]] == [-1, 2]
        @test pts_v[ind[2]] == [0, 3]


        boxes = sort(filter(p -> norm(p, Inf) == 2, pts_v), by=p -> atan(p[1], p[2])) # counter-clockwise order
        push!(boxes, boxes[1]) # close cycle

        c = induced_one_chain(carrier, boxes)
        # test that its a chain
        D = boundaryMatrix(carrier.cplx, 1)
        @test all(==(0), D * c)

        # test annotation
        @test annotate_chain(carrier, c) == [1] || annotate_chain(carrier, c) == [-1]
    end
end

# ==============================================================================
# CUBICAL COMPARISON SPACE TESTS
# ==============================================================================

@testset "CubicalComparisonSpace" begin
    # l_infty circle data
    circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(20, 1), dims=2)

    # test cubical_vr_comparison_space_via_cover
    boxsize = 0.5
    comp_space = cubical_vr_comparison_space_via_cover(circle_data, boxsize)
    @test all(==(2), norm.(eachcol(comp_space.carrier.pts), Inf))

    # test edge_boxes
    e_b = map(zip(eachcol(circle_data), Iterators.drop(eachcol(circle_data), 1))) do t
        edge_boxes(comp_space, t...)
    end

    for (bxs1, bxs2) in zip(e_b, Iterators.drop(e_b, 1))
        @test bxs1[end] == bxs2[1]
    end

    # test map_cycle
    n = size(circle_data, 2)
    simplices = map(i -> [i, mod1(i + 1, n)], 1:n)
    coeffs = ones(Int, length(simplices) + 1)

    cycle = map_cycle(comp_space, circle_data, simplices, coeffs)
    @test cycle == [1] || cycle == [-1]

    coeffs2 = ones(FF{7}, length(simplices) + 1)
    cycle2 = map_cycle(comp_space, circle_data, simplices, coeffs2)
    @test cycle2 == [FF{7}(1)] || cycle == [FF{7}(-1)]

    coeffs3 = -coeffs
    cycle3 = map_cycle(comp_space, circle_data, simplices, coeffs3)
    @test cycle3 == -cycle

    # test betti_1
    @test betti_1(comp_space) == 1
end

@testset "SBCubicalComparisonSpace" begin
    # l_infty circle data
    circle_data = mapslices(v -> normalize(v, Inf), circle_time_series(41, 1), dims=2)
    # discrete derivative
    circle_data_dd = circle_data[:, 2:end] - circle_data[:, 1:end-1]

    # data in utb
    utb_circle_data = [circle_data[:, 1:end-1]; circle_data_dd]

    boxsize = .5
    sb_radius = 2
    comp_space = sb_cubical_vr_comparison_space_via_cover(utb_circle_data, boxsize, sb_radius)

    #
    # test edge_boxes
    #
    e_b = map(zip(eachcol(utb_circle_data), Iterators.drop(eachcol(utb_circle_data), 1))) do t
        edge_boxes(comp_space, t...)
    end

    @test all(vcat(e_b...)) do box
        return norm(box[3:4], Inf) == sb_radius
    end

    @test sortslices(unique(hcat(vcat(e_b...)...), dims=2),dims=2) == comp_space.carrier.pts

    @test 0 == count(zip(e_b, Iterators.drop(e_b, 1))) do t
        bxs1, bxs2 = t
        return bxs1[end] != bxs2[1]
    end

    eb2 = edge_boxes(comp_space, utb_circle_data[:,1], utb_circle_data[:,5])
    @test all(==(4), length.(eb2))

    #
    # test map cycle
    #
    n = size(utb_circle_data, 2)
    simplices = map(i -> [i, mod1(i + 1, n)], 1:n)
    coeffs = ones(Int, length(simplices) + 1)

    cycle = map_cycle(comp_space, utb_circle_data, simplices, coeffs)
    @test cycle == [1] || cycle == [-1]

    coeffs2 = ones(FF{7}, length(simplices) + 1)
    cycle2 = map_cycle(comp_space, utb_circle_data, simplices, coeffs2)
    @test cycle2 == [FF{7}(1)] || cycle == [FF{7}(-1)]

    coeffs3 = -coeffs
    cycle3 = map_cycle(comp_space, utb_circle_data, simplices, coeffs3)
    @test cycle3 == -cycle
end
