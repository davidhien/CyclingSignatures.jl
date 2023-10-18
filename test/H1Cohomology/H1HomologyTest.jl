using Test

using SmithNormalForm
using SparseArrays
using LinearAlgebra
using Distances

include("../../src/AT-Tools-include-file.jl")

function circle_without_center(n)
    z = zeros(Int,2,n*n)

    for i = 0:n-1
        for j = 0:n-1
            if i ==j && i== (n-1)/2
            else
                z[:,1 + j + n*i] = [i;j]
            end
        end
    end
    return z
end

function unitCube(k)
    v = [parse.(Int,split(bitstring(i)[end-k+1:end],"")) for i = 0:2^k-1]
    return foldl(hcat,v)'
end

function doubleCircle(n,k)
    x1 = [[i;0] for i=0:n]
    x2 = [[n;i] for i=1:k]
    x3 = [[i;k] for i=n-1:-1:0]
    x4 = [[0;i] for i=k-1:-1:1]
    x = [x1;x2;x3;x4]
    x = foldl(hcat,x)
    return [x [-1;1].*x]
end

function cubicalCoverPoints(dataq)
    k = size(dataq,1)
    v = [parse.(Int,split(string(i,base=3,pad=k)[end-k+1:end],""))-ones(Int,k) for i = 0:3^k-1]
    v = foldl(hcat,v)
    datac = [dataq[:,i] .+ v for i=1:size(dataq,2)]
    return foldl(hcat,datac)
end

function complexFromPoints(mat)
    return vrIncremental(mat, Distances.chebyshev, 1)
end

function testHomology1()
    v = [parse.(Int,split(bitstring(i)[end-2:end],"")) for i = 0:7]
    data = foldl(hcat,v)'
    k = 1
    data = [data .+ [0 0 0]; data .+ [k 0 0]; data .+ [2*k 0 0]; data .+ [0 k 0]; data .+ [0 2*k 0]]
    k = 4
    data = [data .+ [0 0 0]; data .+ [k 0 0]; data .+ [0 k 0]; data .+ [k k 0]]

    data = unique(data, dims=1)
    c = complexFromPoints(data')

    gen = cohomologyGenerators(c,1)

    return size(gen,2) == 1
end

function testHomology2()
    v = [parse.(Int,split(bitstring(i)[end-2:end],"")) for i = 0:7]
    data = foldl(hcat,v)'
    k = 1
    data = [data .+ [0 0 0]; data .+ [k 0 0]; data .+ [2*k 0 0]; data .+ [0 k 0]; data .+ [0 2*k 0]]
    k = 4
    data = [data .+ [0 0 0]; data .+ [k 0 0]; data .+ [0 k 0]; data .+ [k k 0]]
    k = 8
    data = [data .+ [0 0 0]; data .+ [k 0 0]; data .+ [0 k 0]; data .+ [k k 0]]

    data = unique(data, dims=1)
    c = complexFromPoints(data')

    gen = cohomologyGenerators(c,1)

    return size(gen,2) == 9
end

function testHomology3()
    c = complexFromPoints(unitCube(3)')
    gen = cohomologyGenerators(c,1)
    return size(gen,2) == 0
end

function testHomology4()
    c = complexFromPoints(unique(doubleCircle(7,7), dims=2))
    gen = cohomologyGenerators(c,1)
    return size(gen,2) == 2
end

function testHomology5()
    cubes = unique(doubleCircle(7,7), dims=2)
    cubes = unique(cubicalCoverPoints(cubes),dims=2)
    c = complexFromPoints(cubes)
    gen = cohomologyGenerators(c,1)
    return size(gen,2) == 2
end

@testset "Cohomology Computations" begin
    @test testHomology1()
    @test testHomology2()
    @test testHomology3()
    @test testHomology4()
    @test testHomology5()
end;
