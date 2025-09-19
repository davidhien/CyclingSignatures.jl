using CyclingSignatures
using JLD2
using LinearAlgebra
using Distances
using BenchmarkTools

Y_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Y_lorenz")
Z_lorenz = load("test/CyclingSignaturesTest/lorenz-data.jld2", "Z_lorenz")




###
### benchmark edge_boxes_utb
###

function utb_edge_boxes_easy_1(boxsize, sb_radius, x1, x2)
    d = div(length(x1),2)
    # check if they are in same box and return the box if true
    p_dist = maximum(i->abs(round(Int, x1[i] / boxsize) - round(Int, x2[i] / boxsize)), 1:d)
    v_dist = maximum(i->abs(round(Int, sb_radius*x1[i]) - round(Int,sb_radius*x2[i])), 1+d:2*d)

    if p_dist == 0 && v_dist == 0
        return [round.(Int,[boxsize*x1[1:d];sb_radius*x1[d+1:2*d]]) ]
    elseif p_dist <= 1 && v_dist <= 1
        return [round.(Int, [boxsize*x1[1:d];sb_radius*x1[d+1:2*d]]), round.(Int,[boxsize*x2[1:d];sb_radius*x2[d+1:2*d]])]
    end

    return [round.(Int, [boxsize*x1[1:d];sb_radius*x1[d+1:2*d]]), round.(Int,[boxsize*x2[1:d];sb_radius*x2[d+1:2*d]])]
end

function test_utb_easy_1()
    M0 = -1 .+ 2* rand(4, 10000)
    M = cumsum(M0, dims=2)
    boxsize = 1
    sb_radius = 1
    b = @benchmarkable utb_edge_boxes_easy_1($boxsize, $sb_radius, x1,x2) setup=(i = rand(1:size(M, 2)-1); x1=@view(M[:, i]); x2=@view(M[:, i+1]))
    return run(b)
end

function utb_edge_boxes_easy_2(boxsize, sb_radius, x1, x2)
    d = div(length(x1),2)
    p1 = x1[1:d]
    p2 = x2[1:d]
    v1 = x1[d+1:2*d]
    v2 = x2[d+1:2*d]

    b1 = round.(Int, [boxsize*p1; sb_radius*v1])
    b2 = round.(Int, [boxsize*p2; sb_radius*v2])

    if b1 == b2
        return [b1]
    elseif abs(maximum(b1-b2)) <= 1
        return [b1, b2]
    end

    return [b1,b2]
end

function test_utb_easy_2()
    M0 = -1 .+ 2* rand(4, 10000)
    M = cumsum(M0, dims=2)
    boxsize = 1
    sb_radius = 1
    b = @benchmarkable utb_edge_boxes_easy_2($boxsize, $sb_radius, $M[:, i], $M[:, i+1]) setup=(i = rand(1:size(M, 2)-1))
    return run(b)
end

function utb_edge_boxes_easy_3(boxsize, sb_radius, x1, x2)
    d = div(length(x1),2)
    p1 = @view(x1[1:d])
    p2 = @view(x2[1:d])
    v1 = @view(x1[d+1:2*d])
    v2 = @view(x2[d+1:2*d])

    b1 = round.(Int, [boxsize*p1; sb_radius*v1])
    b2 = round.(Int, [boxsize*p2; sb_radius*v2])

    if b1 == b2
        return [b1]
    elseif abs(maximum(b1-b2)) <= 1
        return [b1, b2]
    end

    return [b1,b2]
end

function test_utb_easy_3()
    M0 = -1 .+ 2* rand(4, 10000)
    M = cumsum(M0, dims=2)
    boxsize = 1
    sb_radius = 1
    b = @benchmarkable utb_edge_boxes_easy_3($boxsize, $sb_radius, $M[:, i], $M[:, i+1]) setup=(i = rand(1:size(M, 2)-1))
    return run(b)
end

function utb_edge_boxes_easy_4(boxsize, sb_radius, x1, x2)
    d = div(length(x1), 2)
    b1 = Vector{Int}(undef, 2d)
    b2 = Vector{Int}(undef, 2d)
    @inbounds for i in 1:d
        b1[i] = round(Int, boxsize * x1[i])
        b2[i] = round(Int, boxsize * x2[i])
        b1[d+i] = round(Int, sb_radius * x1[d+i])
        b2[d+i] = round(Int, sb_radius * x2[d+i])
    end

    # Compute distances
    p_dist = maximum(abs.(b1[1:d] .- b2[1:d]))
    v_dist = maximum(abs.(b1[d+1:end] .- b2[d+1:end]))

    if p_dist == 0 && v_dist == 0
        return [b1]
    elseif p_dist <= 1 && v_dist <= 1
        return [b1, b2]
    end

    return [b1, b2]
end

function test_utb_easy_4()
    M0 = -1 .+ 2* rand(4, 10000)
    M = cumsum(M0, dims=2)
    boxsize = 1
    sb_radius = 1
    b = @benchmarkable utb_edge_boxes_easy_4($boxsize, $sb_radius, $M[:, i], $M[:, i+1]) setup=(i = rand(1:size(M, 2)-1))
    return run(b)
end

function utb_edge_boxes_easy_4(boxsize, sb_radius, x1, x2)
    d = div(length(x1), 2)
    b1 = Vector{Int}(undef, 2d)
    b2 = Vector{Int}(undef, 2d)
    @inbounds for i in 1:d
        b1[i] = round(Int, boxsize * x1[i])
        b2[i] = round(Int, boxsize * x2[i])
        b1[d+i] = round(Int, sb_radius * x1[d+i])
        b2[d+i] = round(Int, sb_radius * x2[d+i])
    end

    # Compute distances
    dist = maximum(i -> abs(b1[i] - b2[i]), 1:d)

    if dist == 0
        return [b1]
    elseif dist == 1
        return [b1, b2]
    end

    return [b1, b2]
end

function test_utb_easy_4()
    M0 = -1 .+ 2* rand(4, 10000)
    M = cumsum(M0, dims=2)
    boxsize = 1
    sb_radius = 1
    b = @benchmarkable utb_edge_boxes_easy_4($boxsize, $sb_radius, $M[:, i], $M[:, i+1]) setup=(i = rand(1:size(M, 2)-1))
    return run(b)
end

test_utb_easy_1()
test_utb_easy_2()
test_utb_easy_3()
test_utb_easy_4()
