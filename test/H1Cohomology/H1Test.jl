using Test
using LinearAlgebra
using SparseArrays
using Distances
using Random

#=
include("../../src/H1Cohomology/H1.jl")
include("../../src/H1Cohomology/ChainComplex.jl")
include("../../src/H1Cohomology/AbstractCell.jl")
include("../../src/H1Cohomology/Cubical.jl")
include("../../src/H1Cohomology/Simplicial.jl")
#include("../../src/H1Cohomology/CubicalConstructions.jl")
include("../../src/H1Cohomology/Complex.jl")
=#
include("../../test/H1Cohomology/ZMatrixAlgebra.jl")
include("../../test/H1Cohomology/H1LinealMeschulam.jl")
include("../../src/AT-Tools-include-file.jl")

# util methods
function findZeroCols(A)
	return findall(mapslices(v -> all(v .== 0),A,dims=1 )[:])
end

function nzpart(A)
    ind2 = mapslices(v -> count(v .!= 0), A, dims=1)[:]
    ind1 = mapslices(v -> count(v .!= 0), A, dims=2)[:]
    ind2 = findall(ind2 .!= 0)
    ind1 = findall(ind1 .!= 0)
    return A[ind1,ind2]
end

function testDenseRowReduce1()
	A = [8;4;2;8;16]'
	B = copy(A)
	piv, ko = reduceRow!(1,B)

	ker = Matrix(Int(1)*I,5,5)[:,findall(B[:] .== 0)]
	for op in ko
		restore!(ker,op)
	end

	return all(A*ker .== 0)
end

function testDenseRowReduce2()
	A = [10 12 -2 10; -20 6 14 16]
	ignoreCol = [true;false;false;false]
	B = copy(A)
	piv, ko = reduceRow!(1,B, ignoreCol=ignoreCol)
	ignoreCol[piv] = true
	piv, ko = reduceRow!(2,B, ko, ignoreCol=ignoreCol)

	ker = Matrix(Int(1)*I,4,4)[:,findZeroCols(B)]
	for op in ko
		restore!(ker,op)
	end

	return all(A*ker .== 0) && gcd(ker) == 1
end

function testGetMinEntryRowAndNumberOfZeroCols()
	A =[
		2  -2 20 0  2  -2;
		3  -3  0 0  3  -3;
	    1   1  1 0  1   1]
	colHasPivot = [false;false;true;false;false;false]
	rowHasPivot = [false;false;true]
	minRow, nzcols = getMinEntryRowAndNumberOfZeroCols(A,rowHasPivot,colHasPivot)

	return minRow ==1 && nzcols == 1
end

function testDenseKernel1()
	A = [1 -12 11 0; -12 -13 -7 16]
	B = copy(A)
	ker = d1Kernel!(copy(A))

	return all(A*ker .== 0)
end

function testDenseKernel2()
	A = [-14 -5 -3 -5 -8 -7 17 -12 7 -10 -14 -1 -19 13 -7 13 7 -8 14 -17;
		12 7 18 -18 -3 7 9 -18 2 9 -7 12 18 11 15 3 16 -19 -16 -4;
		6 -11 -16 -12 -16 -1 -18 14 3 -18 16 -13 -10 -11 -11 3 5 10 -19 -3;
		9 1 16 6 1 -15 16 -20 11 19 9 -19 -8 -1 14 -2 -13 -11 -2 -13;
		15 -5 -11 2 14 -15 -20 2 5 7 -12 18 -1 4 19 16 13 -3 -3 11;
		14 3 -11 7 8 -11 19 -1 19 20 -19 -15 14 16 -10 17 7 3 10 -7;
		10 -20 20 -19 14 6 18 13 20 -3 14 -12 19 5 -9 -9 -5 -16 -16 11;
		-14 9 -10 -4 -8 3 -17 7 11 -16 -14 9 -4 7 -6 0 11 19 9 5;
		18 2 10 19 11 19 3 -17 6 -7 -20 -18 11 13 10 -18 19 -6 1 -5;
		20 -11 17 -11 -7 16 -7 7 -19 9 -19 -19 -6 -2 4 20 -13 20 -15 8] # a 10x20 matrix
	B = copy(A)
	ker = d1Kernel!(copy(A))

	snf,_,_ = smithForm(ZMatrix(Matrix(ker)))
	return size(ker,2) == 10 && all(A*ker .== 0) && all(diag(snf.B) .== 1)
end

function testDenseKernel3()
	A = [0 0 0 0 0; 18 6 -8 2 0; 20 -15 8 7 0]
	B = copy(A)

	ker = d1Kernel!(B, [false;false;false;false;true])

	return all(A*ker .== 0)
end

function testRestoreBasicReduction()
	A = [2 -4 -4 3; 18 6 -8 2; 20 -15 8 7]

	A = [A [17;0;0]]
	B = copy(A)

	redOp = ReductionOperation(5,A[1,:])
	B[1,:] = zeros(Int,1,5)

	ker = d1Kernel!(B, [false;false;false;false;true])

	restore!(ker, redOp)
	snf,_,_ = smithForm(ZMatrix(Matrix(ker)))
	return all(A*ker .== 0) && all(diag(snf.B) .== 1)
end

function testSparseKernel1()
	A = [1 -12 11 0; -12 -13 -7 16]
	B = sparse(copy(A))
	ker = d1Kernel!(B)

	return all(A*ker .== 0)
end

function testSparseKernel2()
	A = [-14 -5 -3 -5 -8 -7 17 -12 7 -10 -14 -1 -19 13 -7 13 7 -8 14 -17;
		12 7 18 -18 -3 7 9 -18 2 9 -7 12 18 11 15 3 16 -19 -16 -4;
		6 -11 -16 -12 -16 -1 -18 14 3 -18 16 -13 -10 -11 -11 3 5 10 -19 -3;
		9 1 16 6 1 -15 16 -20 11 19 9 -19 -8 -1 14 -2 -13 -11 -2 -13;
		15 -5 -11 2 14 -15 -20 2 5 7 -12 18 -1 4 19 16 13 -3 -3 11;
		14 3 -11 7 8 -11 19 -1 19 20 -19 -15 14 16 -10 17 7 3 10 -7;
		10 -20 20 -19 14 6 18 13 20 -3 14 -12 19 5 -9 -9 -5 -16 -16 11;
		-14 9 -10 -4 -8 3 -17 7 11 -16 -14 9 -4 7 -6 0 11 19 9 5;
		18 2 10 19 11 19 3 -17 6 -7 -20 -18 11 13 10 -18 19 -6 1 -5;
		20 -11 17 -11 -7 16 -7 7 -19 9 -19 -19 -6 -2 4 20 -13 20 -15 8] # a 10x20 matrix
	B = sparse(copy(A))
	ker = d1Kernel!(B)

	snf,_,_ = smithForm(ZMatrix(Matrix(ker)))
	return size(ker,2) == 10 && all(A*ker .== 0) && all(diag(snf.B) .== 1)
end

function testRestoreBasicReductionSparse()
	A = [2 -4 -4 3; 18 6 -8 2; 20 -15 8 7]

	A = [A [17;0;0]]
	B = sparse(copy(A))

	redOp = ReductionOperation(5,B[1,:])
	B[1,:] = spzeros(Int,1,5)

	ker = d1Kernel!(B, [false;false;false;false;true])

	restore!(ker, redOp)
	snf,_,_ = smithForm(ZMatrix(Matrix(ker)))
	return all(A*ker .== 0) && all(diag(snf.B) .== 1)
end

@testset "Unit Tests" begin
	@testset "Dense Row Reduce" begin
		@test testDenseRowReduce1()
		@test testDenseRowReduce2()
	end
	@testset "Dense Kernel" begin
		@test testGetMinEntryRowAndNumberOfZeroCols()
		@test testDenseKernel1()
		@test testDenseKernel2()
		@test testDenseKernel3()
	end
	@testset "Restore Basic Reduction" begin
		@test testRestoreBasicReduction()
	end
	@testset "Sparse Kernel" begin
		@test testSparseKernel1()
		@test testSparseKernel2()
	end
	@testset "Sparse Basic Reduction" begin
		@test testRestoreBasicReductionSparse()
	end
end;

###
### Thick circle
###

function getThickCircle()
    Q1 = makeCube1(0,1,0,1,0,0)
    Q2 = makeCube1(0,1,0,1,1,1)
    Q3 = makeCube1(0,1,0,0,0,1)
    Q4 = makeCube1(0,1,1,1,0,1)
    cubSet = CubicalSet([Q1;Q2;Q3;Q4])
    complex = createCellComplex(cubSet)

    return sparse(coboundaryMatrix(complex,0)), sparse(coboundaryMatrix(complex,1))
end

function testThickCircle()
    D0, D1 = getThickCircle()
    gen = firstCohomology(D0,D1)

    return size(gen,2) == 1 && all(D1*gen .== 0)
end

###
### torus
###

function getTorus(k)
    Q = makeCube1(0,1,0,1)
    torus = primaryFaces(Q)
    for i = 2:k
        torus = torus*primaryFaces(Q)
    end

    return torus
end

function testTorus(k; verbose=false)
    torus = getTorus(k)
    toruscomplex = createCellComplex(torus)
    D0 = sparse(coboundaryMatrix(toruscomplex,0))
    D1 = sparse(coboundaryMatrix(toruscomplex,1))

    t = @elapsed gen = firstCohomology(D0,D1)

    return size(gen,2) == k && all(D1*Matrix(gen) .== 0)
end

function vrcircle(n,r)
    theta = 2*pi*collect(0:1:n)/n
    data = [cos.(theta)'; sin.(theta)']

    return vrIncremental(data, Distances.euclidean, r, maxdim=2)
end

function testVRCircle(n,r)
    if euclidean([1;0],[cos(1/n);sin(1/n)]) > r
        error()
    end
    complex = vrcircle(n,r)
    D0 = sparse(coboundaryMatrix(complex,0))
    D1 = sparse(coboundaryMatrix(complex,1))

    gen = firstCohomology(D0,D1)

    return size(gen,2) == 1 && all(D1*gen .== 0)
end

function testLinealMeschulamUsingSNF(n,k, runs,seed=2189)
    Random.seed!(2189)
    for i = 1:runs
        D0,D1 = lmsTestdata(n,k)

        # betti nr via snf
        snf, _, _ = smithForm(ZMatrix(Matrix(D1)))
        kerOfD1 = count(mapslices(v -> all(v .== 0), snf.B,dims=1))
        betti_1 = kerOfD1 - size(D0,2)+1

        h1 = firstCohomology(D0,D1)
        h1snf,_,_ = smithForm(ZMatrix(h1))

        correctBetti = betti_1 == size(h1,2)
        correctBasis = all(diag(h1snf.B)[1:betti_1] .== 1) && all(D1*h1 .== 0)
        if !(correctBetti && correctBasis)
            return false
        end
    end

    return true
end

function testLinealMeschulamUsingNullspace(n,k, runs,seed=2189)
    Random.seed!(2189)
    for i = 1:runs
        D0,D1 = lmsTestdata(n,k)

        # betti nr via snf
		kerOfD1 = size(nullspace(Matrix(D1)),2)
        betti_1 = kerOfD1 - size(D0,2)+1

        h1 = firstCohomology(D0,D1)

        correctBetti = betti_1 == size(h1,2)
        correctBasis = rank(Matrix(h1)) == betti_1
        if !(correctBetti && correctBasis)
            return false
        end
    end

    return true
end
#=
@testset "Lineal Meschulam With SNF" begin
	@testset "15 vertices" for i = 15:10:455
		@test testLinealMeschulamUsingSNF(15, i, 1)
	end
	@testset "25 vertices 1" for i = 20:20:100
		@test testLinealMeschulamUsingSNF(25, i, 1)
	end
end;=#
@testset "Module Tests" begin
	@testset "Thick Circle" begin
	    @test testThickCircle()
	end
	@testset "Torus Tests" for i = 1:5
	    @test testTorus(i)
	end
	@testset "VR Circle" begin
	    @test testVRCircle(15,1/2)
	    @test testVRCircle(20, 1/3)
	    @test testVRCircle(100, 1/3)
	    @test testVRCircle(300, 1/4)
	end
	@testset "Lineal Meschulam With Nullspace" begin
		@testset "15 vertices" for i = 15:10:455
			@test testLinealMeschulamUsingNullspace(15, i, 1)
		end
		@testset "25 vertices 1" for i = 20:20:260
			@test testLinealMeschulamUsingNullspace(25, i, 1)
		end
		@testset "25 vertices 2" for i = 300:200:2300
			@test testLinealMeschulamUsingNullspace(25, i, 1)
		end
		@testset "35 vertices 1" for i = 100:20:300
			@test testLinealMeschulamUsingNullspace(25, i, 1)
		end
		@testset "35 vertices 2" for i = 300:100:1200
			@test testLinealMeschulamUsingNullspace(25, i, 1)
		end
	end
end;


###
### Test Custom Sparse Matrix Data Structure
###
function randomCoredSparseMatrix()
	m,n = 100, 100
	testcsc = sprand(Int,m, n, .75)
	testcrs = CoredSparseMatrix(testcsc)

	initialNPerRow = nzPerRow(testcsc)

	nDeletedPerRow = zeros(Int, m)
	for _ = 1:7
		initnnz = nnz(testcsc)
		delvec = rand(Bool, initnnz)
	
		counter = 1
		for (i,j,v) in collect(zip(findnz(testcsc)...))
			if delvec[counter] && nDeletedPerRow[i] < initialNPerRow[i]
				testcsc[i,j] = 0
				testcrs[i,j] = 0
				nDeletedPerRow[i] += 1
			end
			counter += 1
		end

		if testcsc != testcrs
			return false
		end

		for i = 1:m
			for k = 1:nDeletedPerRow[i]
				if rand(Bool)
					j = rand(1:n)
					v = 1
					testcsc[i,j] = v
					testcrs[i,j] = v
					nDeletedPerRow[i] -= 1
				end
			end
		end
	end

	return testcsc == testcrs
end

function randomCoredSparseMatrixCoredOpTest()
	m,n = 100, 100
	testcsc = sprand(Int,m, n, .8)
	testcrs = CoredSparseMatrix(testcsc)

	# pick random entries and apply coredOp
	for _ = 1:100
		# construct random CoredOp and apply to both matrices
		i = rand(findall(nzPerRow(testcsc) .>= 2))
		currow = findall(testcrs[i,:] .!= 0)
		i1, i2 = rand(findall(currow .!= 0), 2)
		v1, v2 = currow[[i1;i2]]
		crOp = CoreductionOperation(i1,v1,i2,v2)

		apply!(testcsc, crOp)		
		apply!(testcrs, crOp)

		if testcsc != testcrs
			return false
		end
	end
	return true
end

@testset "custom sparse matrix test" begin
	@test randomCoredSparseMatrix()
	@test randomCoredSparseMatrixCoredOpTest()
end
