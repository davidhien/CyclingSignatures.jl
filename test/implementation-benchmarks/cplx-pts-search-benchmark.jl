using BenchmarkTools
using StatsBase

# What is a good way to quantize in the unit tangent bundle
# - use eachcol -> Iterators.map to quantize -> unique
# - copy, do operations and then unique with dim=2

# First, test for normal space without utb
function quantize_sb_and_sort_simple(pts, boxsize, sb_radius)
    d = div(size(pts, 1), 2)

    quantized_pts = Matrix{Int}(undef, size(pts))
    quantized_pts[1:d,:] = round.(Int, pts[1:d,:] ./ boxsize)
    quantized_pts[d+1:2*d,:] = mapslices(v-> round.(Int, normalize(v,Inf) * sb_radius), pts[d+1:2*d,:], dims=1)
    comp_space_pts = sortslices(unique(quantized_pts, dims=2), dims=2)
    return comp_space_pts
end

function quantize_sb_and_sort_simple_2(pts, boxsize, sb_radius)
    d = div(size(pts, 1), 2)

    quantized_pts = Matrix{Int}(undef, size(pts))
    quantized_pts[1:d,:] = round.(Int, pts[1:d,:] ./ boxsize)

    for j in 1:size(pts, 2)
        quantized_pts[d+1:2*d,j] = round.(Int, normalize(pts[d+1:2*d,j],Inf) * sb_radius)
    end

    comp_space_pts = sortslices(unique(quantized_pts, dims=2), dims=2)
    return comp_space_pts
end

pts = 20*rand(4,100000)
boxsize = 0.25
sb_radius = 2

b1 = @benchmarkable quantize_sb_and_sort_simple($pts, $boxsize, $sb_radius)
res1 = run(b1)

b2 = @benchmarkable quantize_sb_and_sort_simple_2($pts, $boxsize, $sb_radius)
res2 = run(b2)

b3 = @benchmarkable quantize_and_sort_it($pts, $boxsize)
res3 = run(b3)


res1
res2
res3
