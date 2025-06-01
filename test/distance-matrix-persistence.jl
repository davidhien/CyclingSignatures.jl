using Test
using CyclingSignatures
using LinearAlgebra, Distances
import DataStructures: num_groups
import CyclingSignatures: curve_cycle, vertex_to_essential_bar, component_representatives
import CyclingSignatures: dm_components_first_pass_explicit, dm_components_second_pass_explicit
import CyclingSignatures: dm_components_first_pass_implicit, dm_components_second_pass_implicit
import CyclingSignatures: dm_components_explicit, dm_components_implicit
include("test-util.jl")

@testset "curve_cycle" begin
    field = FF{3}
    edges, coeffs = curve_cycle(1,5;F=field)

    @test edges == [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
    @test coeffs == field.([1;1;1;1;-1])
end

@testset "vertex_to_essential_bar" begin
    field = FF{3}
    trajPoints = [1, 2, 3, 4, 5]'[:,:]
    metric(x,y) = norm(x-y)
    node = (1,5)
    bar = vertex_to_essential_bar(node, trajPoints, metric; field=field)

    @test bar.birth == 4.0
    @test bar.death == Inf
    @test bar.simplex_list == [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
    @test bar.coeff_list == field.([1;1;1;1;-1])
end

@testset "dm components individual passes" begin
    points = [1 0 -1 0; 0 1 0 -1]
    metric = chebyshev
    fltThreshold = 1.1
    @testset "explicit" begin
        cc, cc_labels, d_mat = dm_components_first_pass_explicit(points, metric, fltThreshold)
        @test size(d_mat) == size(cc_labels)
        @test num_groups(cc) == 2

        min_vertex_dict = dm_components_second_pass_explicit(cc,d_mat,cc_labels)
        @test length(collect(values(min_vertex_dict))) == 2

        min_vertices = component_representatives(min_vertex_dict, cc_labels)
        @test min_vertices == [(1, 4)]
    end
    @testset "implicit" begin
        cc, cc_labels = dm_components_first_pass_implicit(points, metric, fltThreshold)
        @test size(cc_labels) == (4,4)
        @test num_groups(cc) == 2

        min_vertex_dict = dm_components_second_pass_implicit(cc,cc_labels, points, metric)
        @test length(collect(values(min_vertex_dict))) == 2

        min_vertices = component_representatives(min_vertex_dict, cc_labels)
        @test min_vertices == [(1, 4)]
    end
end

@testset "dm components" begin
    points = circle_time_series(9, 2)
    metric = euclidean
    flt_threshold = 0.8
    @testset "explicit" begin
        min_vertices_expl = dm_components_explicit(points, metric, flt_threshold)
        @test length(min_vertices_expl) == 2

        v = sort(map(t-> t[2] - t[1] + 1, min_vertices_expl)) # vertex diagonals
        @test v[1] >=  9
        @test v[1] <= 11
        @test v[2] >= 18
        @test v[2] <= 19
    end
    @testset "implicit" begin
        min_vertices_impl = dm_components_implicit(points, metric, flt_threshold)
        @test length(min_vertices_impl) == 2

        v = sort(map(t-> t[2] - t[1] + 1, min_vertices_impl))  # vertex diagonals
        @test v[1] >=  9
        @test v[1] <= 11
        @test v[2] >= 18
        @test v[2] <= 19
    end

    points = double_circle_time_series(10,[0,1,0])
    metric = euclidean
    flt_threshold = 0.8

    @testset "explicit" begin
        min_vertices_expl = dm_components_explicit(points, metric, flt_threshold)
        @test length(min_vertices_expl) == 5
    end
    @testset "implicit" begin
        min_vertices_impl = dm_components_implicit(points, metric, flt_threshold)
        @test length(min_vertices_impl) == 5
    end
end


@testset "trajectory_barcode" begin
    # TODO: add test for return value and type stability
end
