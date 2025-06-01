using Test
using CyclingSignatures: edge_boxes, projection_initial_max_component, projection_max_switches
using CyclingSignatures: edge_boxes_sphere_bundle

@testset "edge_boxes" begin
    p1 = [0;0]
    p2 = [10;0]

    res = [[i;0] for i = 0:10]
    @test res == edge_boxes(p1,p2)

    p1 = [0;0.25]
    p2 = [5;5.25]
    res = [[div(i,2);div(i+1,2)] for i = 0:10]
    @test res == edge_boxes(p1,p2)

    p1 = [4.2;4.2]
    p2 = [4.2;4.2]
    @test length(edge_boxes(p1,p2)) == 1
end

@testset "edge_boxes_utb" begin
    p1 = [0,0]
    p2 = [0,0]
    v1 = [1,0]
    v2 = [1,1]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 1)
    @test eb == [[0, 0, 1, 0],[0, 0, 1, 1]]

    for k = 1:10
        p1 = [0,0]
        p2 = [0,0]
        v1 = [1,0]
        v2 = [0,1]
        eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, k)
        res = Vector{Int}[]
        for i=0:k
            push!(res, [0;0;k;i])
        end
        for i=k-1:-1:0
            push!(res, [0;0;i;k])
        end
        @test res == eb
    end

    p1 = [0,0]
    p2 = [0,0]
    v1 = [1,0]
    v2 = [0,1]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 2)
    @test eb == [[0, 0, 2, 0], [0, 0, 2, 1], [0, 0, 2, 2], [0, 0, 1, 2], [0, 0, 0, 2]]


    p1 = [0,0]
    p2 = [10,0]
    v1 = [1,0]
    v2 = [1,0]
    eb = edge_boxes_sphere_bundle(p1, p2, v1, v2, 1)
    res = [[i;0;1;0] for i = 0:10]
    @test res == eb
end

@testset "projection_initial_max_component" begin
    v1 = [0.8,0.9,0.9]
    v2 = [1.0,0.6,0.5]
    @test projection_initial_max_component(v1,v2) == 2

    v1 = [0.8, 0.9, 0.9]
    v2 = [1.0, -.6, -.5]
    @test projection_initial_max_component(v1,v2) == 3

    v1 = [0.8,0.9,-0.9]
    v2 = [1.0,0.6,-0.5]
    @test projection_initial_max_component(v1,v2) == 2

    v1 = [0.8,-0.9, 0.9]
    v2 = [1.0, .6, -.5]
    @test projection_initial_max_component(v1,v2) == 3
end

@testset "projection_max_switches" begin
    @testset "explicit cases" begin
        # line a over line b
        v1 = [2;1]
        v2 = [3;2]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls, [0.0, 1.0]; atol=1e-8)
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls, [0.0, 1.0]; atol=1e-8)
        ls = projection_max_switches(v1, -v2)
        @test isapprox(ls, [0.0, 0.375, 0.5, 1.0]; atol=1e-8)
        ls = projection_max_switches(-v1, v2)
        @test isapprox(ls, [0.0, 0.375, 0.5, 1.0]; atol=1e-8)

        # line a intersects line b
        v1 = [2;1]
        v2 = [1;2]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])
        ls = projection_max_switches(v1, -v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])

        # line a over line b, b has irrelevant zero intersection
        v1 = [2;1]
        v2 = [1;-.5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls,[0.0, 1.0])

        # line a positive, line b negative and takes max
        v1 = [2;0]
        v2 = [4;-6]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])

        # like previous but flipped sign at second component
        v1 = [2;0]
        v2 = [4;6]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1, -v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])

        # degenerate intersection on right endpoint
        v1 = [2;0]
        v2 = [-1;1]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])

        # start at same point, one in positive one in negative direction
        v1 = [1;1]
        v2 = [3;-5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])
        ls = projection_max_switches(-v1,-v2)
        @test isapprox(ls,[0.0, 0.5, 1.0])

        # two intersections
        v1 = [1;3]
        v2 = [1;-5]
        ls = projection_max_switches(v1, v2)
        @test isapprox(ls,[0.0, 0.25, 0.5, 1.0])
        ls = projection_max_switches(-v1,-v2)
        @test isapprox(ls,[0.0, 0.25, 0.5, 1.0])

        # intersection with constant zero function
        v1 = [1;0]
        v2 = [1;0]
        ls = projection_max_switches(v1,v2)
        @test isapprox(ls,[0.0, 1.0])
        ls = projection_max_switches(-v1,-v2)
        @test isapprox(ls,[0.0, 1.0])

        # multiple intersections
        # TODO!!!
    end
    @testset "random tests" begin
        for _ in 1:10
            v1 = rand(10)
            v2 = rand(10)
            ls = projection_max_switches(-v1,-v2)
            for l in ls[2:end-1]
                v = (1-l)*v1 + l*v2
                m = maximum(abs, v)
                m_i = findall(isapprox(m), v)
                @test length(m_i) >= 2
            end
        end
    end
end
