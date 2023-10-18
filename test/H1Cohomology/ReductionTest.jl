using Test
using SparseArrays

include("../../src/H1Cohomology/Reduction.jl")

"""
test example 1
"""
function testFindMaxQinChain1()
    c = [1;0;1;0]
    o = [0;1;2;3]

    ind, weight = findMaxQinChain(c,o)

    return ind = 3 && weight = 2
end

"""
test example 2
"""
function testFindMaxQinChain2()
    c = [1;0;1;0]
    o = [3;2;1;0]

    ind, weight = findMaxQinChain(c,o)

    return ind = 1 && weight = 3
end

"""
test default behaviour
"""
function testFindMaxQinChain3()
    c = [1;0;1;0]
    o = [0;1;0;3]

    ind, weight = findMaxQinChain(sparse(c),o)

    return ind == 0 && weight == -1
end

function gammaTestExample1()
    D = spzeros(Int, 9, 15)
    D[7,15] = D[9,14] = D[8,13] = D[9,12] = D[8,11] = D[7,10] = D[4,9] = D[6,8] = D[5,7] = D[6,6] = D[5,5] = D[4,4] = D[1,3] = D[3,2] = D[2,1] = 1
    D[9,15] = D[8,14] = D[7,10] = D[6,12] = D[5,11] = D[4,10] = D[6,9] = D[5,8] = D[4,7] = D[3,6] = D[2,5] = D[1,4] = D[3,3] = D[2,2] = D[1,1] = -1

    order = [-1 1 2 3 4 5 6 7 8]
    w = Dict{Int,Int}()
    w[2] = 1
    w[3] = 2
    w[4] = 4
    w[5] = 5
    w[6] = 6
    w[7] = 10
    w[8] = 11
    w[9] = 12

    return w, order, D
end

function gammaTestExample2()

end

"""
test example 1
0 dim chain in cylinder
"""
function testGammaAlgorithm1()
    w, order, D = gammaTestExample1()
    c = zeros(Int, 9)
    c[5] = 1

    v = gammaAlgorithm(c, w, order, D)
    z = c + D*v

    return z[1] != 0 && all(z[2:end] .== 0)
end

"""
test example 1
1 dim chain in cylinder
"""
function testGammaAlgorithm2()
    error("TODO")
end

"""
test example 1
0 dim chain in torus
"""
function testGammaAlgorithm3()
    error("TODO")
end

"""
test example 1
1 dim chain in torus
"""
function testGammaAlgorithm4()
    error("TODO")
end


function getTorus()
    X = CubicalSet(2)
    for i = 1:3
        for j = 1:3
            c = makeCube2([i;j;i+1;j+1])
            addCube(X,c)
        end
    end

    return X
end
