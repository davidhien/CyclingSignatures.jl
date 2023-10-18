abstract type KernelOperation end

# models restoring row from combination
struct ReductionOperation{T} <: KernelOperation where T <: Integer
	row::Int
	combination::AbstractVector{T}
end

# data for two reductions
struct CoreductionOperation{T} <: KernelOperation where T <: Integer
	keepColIndex::Int
	keepColValue::T
	discardColIndex::Int
	discardColValue::T
end

# models adding a times col j to col i
struct ColAddOperation{T} <: KernelOperation where T <: Integer
	i::Int
	j::Int
	a::T
end

function firstCohomology(D0::AbstractSparseMatrix{T0,Int},
                                D1::AbstractSparseMatrix{T1,Int}) where {T0 <: Integer, T1 <: Integer}
    D1 = copy(D1)
	@info "Basic Reductions"
    basicRed, colRemoved = basicReductions!(D1)
	@info "Removing Max Spanning Tree"
	removeMaxSpanningTree!(D0,D1,colRemoved)

	@info "kernel computation started..."
	return d1Kernel!(D1,basicRed,colRemoved)
end

"""
	function d1Kernel!(D1::AbstractSparseMatrix{T,Int},kerOp::Stack{KernelOperation},preprocess=colRemoved) where {T <: Integer}

Computes the kernel of a 1st boundary matrix.
"""
function d1Kernel!(A::AbstractSparseMatrix{T,Int},kerOp::Stack{KernelOperation}=Stack{KernelOperation}(),colRemoved::Vector{Bool}=zeros(Bool,size(A,2))) where {T <: Integer}
	m,n = size(A)
	@info "coreductions started"
	redCored!(A, kerOp, colRemoved) # start with reductions and coreductions
	@info string("coreductions removed ", count(colRemoved),"/" ,size(A,2)," columns.")
	while ((i,rowNNZ,colNNZ) = findSmallestWeightReducableRow!(A))[1] !== nothing
		nzCols = count(colNNZ .!= 0)
		nzRows = count(rowNNZ .!= 0)
		density = sum(rowNNZ)/(nzCols*nzRows) # TODO: add support for direct sums

		if (density > .1 || (nzCols < 500 && nzRows < 500 )) && rowNNZ[i] > 3
			# go to dense solver
			break
		end

		# remove row
		j,_ = reduceRow!(i, A, kerOp)
		colRemoved[j] = true
		A[:,j] = spzeros(T,m)

		# dropzeros rule
		if nnz(A) - sum(rowNNZ) > 1000
			dropzeros!(A)
		end
	end

	kerCols = findall((colNNZ .== 0) .& (.! colRemoved))

	ker = getUnitVectors(kerCols, n,T)

	if size(A,1) != 0 && any(!=(0),A)
		ind1 = findall(nzPerRow(A) .!= 0)
		ind2 = findall(nzPerCol(A) .!= 0)

		denseA = Matrix(A[ind1,ind2])
		kerDenseA = d1Kernel!(denseA)
		kerDenseBig = spzeros(T,n, size(kerDenseA,2))
		kerDenseBig[ind2,:] = kerDenseA
		ker = [ker sparse(kerDenseBig)]
	end
	ker = restore!(ker, kerOp)
	return ker
end

function findSmallestWeightReducableRow!(A::AbstractSparseMatrix{T,Int}) where {T <: Integer}
	# TODO: refactor
	m,n = size(A)
	rowNNZ = zeros(Int,m)
	rowGCD = zeros(T,m)
	colNNZ = zeros(Int,n)

	# initialize
	for (i,j,v) in zip(findnz(A)...)
		if v != 0
			rowNNZ[i] += 1
			colNNZ[j] += 1
			rowGCD[i] = gcd(rowGCD[i],v)
		end
	end

	# remove gcd
	for i = 1:m
		if abs(rowGCD[i]) > 1
			A[i,:] = div.(A[i,:],rowGCD[i])
			rowGCD[i] = 1
		end
	end

	minEntry = (0,0)
	minEntryWeight = 0
	rowWeights = zeros(Int,m)
	for (i,j,v) in zip(findnz(A)...)
		if abs(v) == 1
			if rowNNZ[i] == 1
				# this entry is free
				return i, rowNNZ, colNNZ
			end
			# calculate entry weight and update min entry weight
			curWeight = rowNNZ[i]*colNNZ[j] - rowNNZ[i] + 1
			if minEntryWeight == 0 || curWeight < minEntryWeight
				minEntry = (i,j)
				minEntryWeight = curWeight
			end
		end
		if v != 0
			# update row weight
			rowWeights[i] += colNNZ[j]
		end
	end

	rowWeights = rowWeights .* rowNNZ
	minRow,	minRowWeight = findminnz(rowWeights)
	if minRow === nothing
		return nothing, rowNNZ, colNNZ
	elseif minEntryWeight != 0 && minEntryWeight <= minRowWeight
		return minEntry[1], rowNNZ, colNNZ
	else
		return minRow, rowNNZ, colNNZ
	end
end

function findminnz(v)
	#TODO: refactor
	minIndex = nothing
	minValue = nothing

	for i = 1:length(v)
		if v[i] != 0
			if minIndex === nothing || v[i]< minValue
				minIndex = i
				minValue = v[i]
			end
		end
	end

	return minIndex, minValue
end

function d1Kernel!(A::AbstractMatrix{T},colRemoved::Vector{Bool}) where {T <: Integer}
	return d1Kernel!(A, Stack{KernelOperation}(), colRemoved)
end

"""
    function d1Kernel!(A::Array{2,T},kerOp::Stack{KernelOperation}=Stack{KernelOperation}(),colRemoved::Vector{Bool}=zeros(Bool,size(D0,1))) where {T <: Integer}

Computes the kernel using row reduction.
"""
function d1Kernel!(A::AbstractMatrix{T},kerOp::Stack{KernelOperation}=Stack{KernelOperation}(), colRemoved::Vector{Bool}=zeros(Bool,size(A,2))) where {T <: Integer}
	@info "dense kernel algorithm started"
	m,n = size(A)
	approxNullspace = nullspace(Float64.(A))

	kerDim = size(approxNullspace,2)

	rowHasPivot = zeros(Bool,m)
	colHasPivot = zeros(Bool,n)

	i, curKer = getMinEntryRowAndNumberOfZeroCols(A, rowHasPivot, colHasPivot)

	while i !== nothing && curKer < kerDim
		j,_ = reduceRow!(i,A,kerOp, ignoreCol = colHasPivot)
		rowHasPivot[i] = true
		colHasPivot[j] = true
		i, curKer = getMinEntryRowAndNumberOfZeroCols(A, rowHasPivot, colHasPivot)
	end

	colIsZero = mapslices(v -> all(v .== 0), A, dims=1)[:]
	kerCols = findall( .!colRemoved .& colIsZero)
	ker = getUnitVectors(kerCols, n,T)

	return restore!(ker,kerOp)
end

function getMinEntryRowAndNumberOfZeroCols(A, rowHasPivot, colHasPivot)
	m,n = size(A)
	isZeroCol = .! colHasPivot
	minRow = nothing
	minRowSum = 0

	for i = 1:m
		rowsum = 0
		if !rowHasPivot[i]
			for j = 1:n
				if !colHasPivot[j]
					rowsum += abs(A[i,j])
					isZeroCol[j] = isZeroCol[j] & (A[i,j] == 0)
				end
			end
		end

		if rowsum != 0 && (rowsum < minRowSum || minRowSum == 0)
			minRow = i
			minRowSum = rowsum
		end
	end
	return minRow, count(isZeroCol)
end

"""
	function getUnitVectors(indices, m, T=Int)

Returns a sparse (m,lenth(indices))-matrix of type T with i-th column being indices[i]-th standard basis vector.
"""
function getUnitVectors(indices, m, T=Int)
	uv = spzeros(T, m, length(indices))

	for i = 1:length(indices)
		uv[indices[i],i] = 1
	end

	return uv
end

"""
	function redCored!(D1::AbstractSparseMatrix{T,Int}, ) where {T <: Integer}

"""
function redCored!(A::AbstractSparseMatrix{T,Int}, kerOp::Stack{KernelOperation}=Stack{KernelOperation}(),colRemoved::Vector{Bool}=zeros(Bool,size(A,2))) where {T <: Integer}
	while size(A,1)>0
		if hasNZColWithLEQThanKEntries(A,2)
			basicReductions!(A, kerOp, colRemoved)
		elseif hasNZRowWithLEQThanKEntries(A,2)
			basicCoreductions!(A, kerOp, colRemoved)
		else
			break
		end
	end
end

function hasNZRowWithLEQThanKEntries(A::AbstractMatrix, k)
	v = nzPerRow(A)
	return findfirst(i -> i <= k && i > 0, v) != nothing
end

function hasNZColWithLEQThanKEntries(A::AbstractMatrix, k)
	v = nzPerCol(A)
	return findfirst(i -> i <= k && i > 0, v) != nothing
end

function nEntriesInMinimumNZCol(A::AbstractMatrix)
	return minimum(nonzeros(dropzeros(nzPerCol(A))))
end

function nzPerRow(A::AbstractMatrix)
	return mapslices(v -> count(v .!= 0), A, dims=2)[:]
end

function nzPerRow(A::AbstractSparseArray)
	count = zeros(Int, size(A,1))
	for (i,j,v) in zip(findnz(A)...)
		if v != 0
			count[i] += 1
		end
	end
	return count
end

function nzPerCol(A::AbstractMatrix)
	return mapslices(v -> count(v .!= 0), A, dims=1)[:]
end

function nzPerCol(A::AbstractSparseArray)
	count = zeros(Int, size(A,2))
	for (i,j,v) in zip(findnz(A)...)
		if v != 0
			count[j] += 1
		end
	end
	return count
end


"""
    function removeMaxSpanningTree!(D0::AbstractSparseMatrix{T0,Int},D1::AbstractSparseMatrix{T1,Int}; colRemoved::Vector{Bool}=zeros(Bool,size(D0,1))) where {T0 <: Integer, T1 <: Integer}

Finds a maximal spanning tree in the graph defined by D0 with the edges in colRemoved ignored.
It then removes the spanning tree by setting the corresponding columns in D1 to zero.
"""
function removeMaxSpanningTree!(D0::AbstractSparseMatrix{T0,Int},D1::AbstractSparseMatrix{T1,Int}, colRemoved::Vector{Bool}=zeros(Bool,size(D0,1))) where {T0 <: Integer, T1 <: Integer}
	st = findMaxSpanningTree(D0,D1,colRemoved)
	sort!(st)
	I,J,K = findnz(D1)
	for i in st
		colRemoved[i] = true
	end
	# set columns which are removed to zero in sparse matrix data structure
	ct = 1
	for l = 1:length(st)
		while ct <= length(J) && J[ct] <= st[l]
			if J[ct] == st[l]
				D1.nzval[ct] = 0
			end
			ct+=1
		end
	end

	dropzeros!(D1)

	return D1, colRemoved
end

"""
    function findMaxSpanningTree(D0,D1; colRemoved::Vector{Bool}=zeros(Bool,size(D0,1)))

Finds a maximal spanning tree in D0. Edge weight is number of entries in the corresponding column in D1.
"""
function findMaxSpanningTree(D0,D1, colRemoved::Vector{Bool}=zeros(Bool,size(D0,1)))
    colnnz = nzPerCol(D1)

    m,n = size(D0)

    # generate edge list
    edgelist = zeros(Int,m,2)
    for (i,j,v) in zip(findnz(D0)...)
        if v != 0
            if edgelist[i,1] == 0
                edgelist[i,1] = j
            else
                edgelist[i,2] = j
            end
        end
    end

    # generate 1-skeleton as weighted graph
    g = MetaGraph(n)
    for i = 1:m
		if !colRemoved[i]
	        MetaGraphs.add_edge!(g, edgelist[i,1], edgelist[i,2])
	        set_prop!(g, edgelist[i,1], edgelist[i,2], :weight, colnnz[i])
	        set_prop!(g, edgelist[i,1], edgelist[i,2], :index, i)
		end
    end

    mst = MetaGraphs.kruskal_mst(g, minimize=false)
    return map(e -> get_prop(g, e, :index), mst)
end


"""
	function basicReductions!(D1::AbstractSparseMatrix{T,Int}, kernelOperations::Stack{KernelOperation}=Stack{KernelOperation}, colRemoved::Vector{Bool}=zeros(Bool,size(D1,2))) where T <: Integer

Performs basic reductions, i.e. reductions of columns with up to two entries.
"""
function basicReductions!(D1::AbstractSparseMatrix{T,Int}, kernelOperations::Stack{KernelOperation}=Stack{KernelOperation}(), colRemoved::Vector{Bool}=zeros(Bool,size(D1,2))) where T <: Integer
    m, n = size(D1)

	colnnz = Vector(nzPerCol(D1))
	colsWithOneNNZ = findall(x -> x == 1, colnnz)
	colsWithTwoNNZ = findall(x -> x == 2, colnnz)

	while !isempty(colsWithOneNNZ) || !isempty(colsWithTwoNNZ)
		while !isempty(colsWithOneNNZ)
			j = pop!(colsWithOneNNZ) # col index
			if colnnz[j] == 1
				i = findnz(dropzeros(D1[:,j]))[1][1] # row index
				w = D1[i,:] # the row
				colRemoved[j] = true

				for (k,v) in zip(findnz(w)...)
					if v != 0
						D1[i,k] = 0
						colnnz[k] -= 1

						if colnnz[k] == 1
							push!(colsWithOneNNZ, k)
						elseif colnnz[k] == 2
							push!(colsWithTwoNNZ, k)
						end
					end
				end

				# update kernel operations
				push!(kernelOperations, ReductionOperation(j,w))
			end
		end

		while !isempty(colsWithTwoNNZ)
			j = pop!(colsWithTwoNNZ)
			if colnnz[j] == 2
				# do Morse-type reduction
				i, vals = findnz(dropzeros(D1[:,j]))

				d,a,b = gcdx(vals...)

				# we have for vals = [x;y] that d = ax+by
				# therefore the basic reduction vector is
				redVector = a*D1[i[1],:] + b*D1[i[2],:]

				# update matrix and heap
				alpha, beta = div(vals[2],d), div(vals[1],d)
				D1[i[1],:] *= alpha
				row1  = D1[i[1],:]
				row2  = D1[i[2],:]

				colRemoved[j] = true
				# update matrix and queues
				for (k,v) in zip(findnz(row2)...)
					if v != 0
						u = row1[k]
						newv = u - beta*v

						D1[i[1],k] = newv
						D1[i[2],k] = 0

						if (newZeros = count([u != 0 && v != 0; newv == 0] )) != 0
							colnnz[k] -= newZeros
							if colnnz[k] == 1
								push!(colsWithOneNNZ, k)
							elseif colnnz[k] == 2
								push!(colsWithTwoNNZ, k)
							end
						end
					end
				end

				# update basic reductions vector
				push!(kernelOperations, ReductionOperation(j,redVector))
			end
		end
	end
	dropzeros!(D1)

    return kernelOperations, colRemoved
end

"""
    function basicReductions!(D1::AbstractSparseMatrix{T,Int}; kernelOperations::Stack{KernelOperation}) where T <: Integer

Performs basic coreductions, i.e. reductions of rows with up to two entries.
"""
function basicCoreductions!(D1::AbstractSparseMatrix{T,Int}, kernelOperations::Stack{KernelOperation}=Stack{kernelOperations}(),
									colRemoved::Vector{Bool}=zeros(Bool,size(A,2))) where {T <: Integer}

	# Start
	A = CoredSparseMatrix(D1)
	#A = D1
	m, n = size(A)
	rownnz = Vector(nzPerRow(D1))
	rowsWithOneNNZ = findall(x -> x == 1, rownnz)
	rowsWithTwoNNZ = findall(x -> x == 2, rownnz)

	onept = length(rowsWithOneNNZ)
	twopt = length(rowsWithTwoNNZ)

	rowsWithOneNNZ = [rowsWithOneNNZ;zeros(Int,m-length(rowsWithOneNNZ))]
	rowsWithTwoNNZ = [rowsWithTwoNNZ;zeros(Int,m-length(rowsWithTwoNNZ))]

	while onept != 0 || twopt != 0
		while onept != 0
			# pop
			i = rowsWithOneNNZ[onept] # row index

			onept -= 1

			if rownnz[i] == 1
				w = dropzeros(A[i,:]) # current row
				j = findnz(w)[1][1] # relevant col index
				colRemoved[j] = true

				for (k,v) in zip(findnz(A[:,j])...)
					if v != 0
						# update matrix
						A[k,j] = 0
						rownnz[k] -= 1
						# update queues
						if rownnz[k] == 1
							onept += 1
							rowsWithOneNNZ[onept] = k
						elseif rownnz[k] == 2
							twopt += 1
							rowsWithTwoNNZ[twopt] = k
						end
					end
				end
			end
		end

		while twopt != 0
			i = rowsWithTwoNNZ[twopt] # row index
			twopt -= 1
			if rownnz[i] == 2 # it might have decresed but the array was not updated
				currow = dropzeros(A[i,:])
				j, vals = findnz(currow)

				# divide out gcd
				d = gcd(vals[1], vals[2])
				vals = div.(vals, d)
				A[i,j] = vals

				# update arrays
				col1 = dropzeros(A[:,j[1]])
				col2 = dropzeros(A[:,j[2]])
				for (k,v) in zip(findnz(col2)...)
					if v != 0
						w = col1[k]
						newv = vals[2]*w - vals[1]*v

						if (newZeros = count([w != 0; newv == 0])) != 0
							rownnz[k] -= newZeros
							if rownnz[k] == 1
								onept += 1
								rowsWithOneNNZ[onept] = k
							elseif rownnz[k] == 2
								twopt += 1
								rowsWithTwoNNZ[twopt] += 1
							end
						end
					end
				end

				# update basic reduction vector
				colRemoved[j[2]] = true # keep 1 remove 2
				crOp = CoreductionOperation(j[1], vals[1], j[2], vals[2])
				push!(kernelOperations, crOp)
				apply!(A, crOp)
			end
		end
	end
	D1 .= SparseMatrixCSC(A)
#	dropzeros!(A)

	return kernelOperations, colRemoved
end

function reduceRow!(i,A::AbstractSparseMatrix{T0,Int},kerStack::Stack{KernelOperation}=Stack{KernelOperation}()) where {T0 <: Integer}
	# TODO: test this better
	innz = 0     # nnz in row i
	minIndex = 0 # store index of entry with min abs value in row i
	minValue = 0 # abs value of entry with min abs value in row i

	# initialize
	iIndices, iVals = findnz(A[i,:]) # those will be kept independently
	colnnz = map(j -> count(A[:,j] .!= 0) , iIndices)
	for (k,v,jnnz) in zip(iIndices, iVals, colnnz)
		if v != 0
			innz += 1

			# initialize with arbitrary key-value pair
			minIndex = k
			minValue = (v,jnnz)
		end
	end

	# main loop
	while innz > 1
		# compute minIndex and minValue
		colnnz = map(j -> count(A[:,j] .!= 0) , iIndices)
		for (k,v,jnnz) in zip(iIndices, iVals, colnnz)
			if v != 0 && (abs(v),jnnz) < abs.(minValue)
				minIndex = k
				minValue = (v,jnnz)
			end
		end

		# reduce row minIndex
		for j = 1:length(iIndices)
			k = iIndices[j]
			v = iVals[j]
			redVal = minValue[1]

			if k != minIndex && v != 0
				colOp = ColAddOperation(k, minIndex, -div(v, redVal))
				apply!(A, colOp)
				push!(kerStack, colOp)
				iVals[j] -= div(v, redVal)*redVal

				if iVals[j] == 0
					innz -= 1
				end
			end
		end
	end

	return minIndex, kerStack
end

function reduceRow!(i,A::AbstractArray{T,2}, kerStack::Stack{KernelOperation}=Stack{KernelOperation}(); ignoreCol::Vector{Bool}=zeros(Bool,size(A,2))) where {T <: Integer}
	# TODO: something here does not work
	m,n = size(A)
	innz = count(A[i,.!ignoreCol] .!= 0)     # nnz in row i
	if innz == 1
		return findfirst(A[i,.!ignoreCol] .!= 0), kerStack
	end
	minIndex = 0 # store index of entry with min abs value in row i
	minValue = (0,0) # abs value of entry with min abs value in row i

	while innz > 1
		# find col to reduce
		for j = 1:n
			if !ignoreCol[j] && A[i,j] != 0
				v = sum(abs.(A[:,j]))
				if minValue[1] == 0 || (abs(A[i,j]),v) < minValue
					minIndex = j
					minValue = (abs(A[i,j]),v)
				end
			end
		end
		v = A[i,minIndex]

		for j = 1:n
			if !ignoreCol[j] && A[i,j] != 0 && j != minIndex
				a = -div(A[i,j], v)
				op = ColAddOperation(j, minIndex, a)

				apply!(A, op)
				push!(kerStack, op)

				if A[i,j] == 0
					innz -= 1
				end
			end
		end
	end
	return minIndex, kerStack
end

"""
    function apply!(A::Array{T,2}, addOp::ColAddOperation{T}) where T <: Integer

Applies a ColAddOperation.
"""
function apply!(A::AbstractArray{T,2}, addOp::ColAddOperation{T}) where T <: Integer
	i,j,a = opData(addOp)
	A[:,i] += a*A[:,j]
end

"""
    function apply!(A::AbstractArray{T,2}, crOp::CoreductionOperation{T}) where T <: Integer

Applies the coreduction operation. See restore for details.
"""
function apply!(A::AbstractArray{T,2}, crOp::CoreductionOperation{T}) where T <: Integer
	keepColIndex, keepColValue, discardColIndex, discardColValue = opData(crOp)
	# TODO: implement this at the level of the sparse matrix representation
	A[:,keepColIndex] = discardColValue*A[:,keepColIndex] - keepColValue*A[:,discardColIndex]
	A[:,discardColIndex] = spzeros(T, size(A,1))

	return A
end

function restore!(A::AbstractArray{T,2}, kerOp::Stack{KernelOperation}) where T <: Integer
	for op in kerOp
		restore!(A, op)
	end

	return A
end

"""
    function restore!(A::Array{T,2}, basicReduction::ReductionOperation{T}) where T <: Integer

Restores the kernel after a basic reduction operation (where a row gets removed).
"""
function restore!(A::AbstractArray{T,2}, basicReduction::ReductionOperation{T}) where T <: Integer
	i,comb = opData(basicReduction)

	a = comb[i]
	w = comb'*A

	if all(mod.(w,a).== 0)
    	A[i,:] = -div.(w,a)
	else
		j, kerStack = reduceRow!(1,w)
		val = w[j]
		for op in Iterators.reverse(kerStack)
			apply!(A, op)
		end

		g = gcd(a,val)
		A[:,j] *= div(a,g)
		A[i,j] = -div(val,g)
    end
end

"""
    function restore!(A::Array{T,2}, addOp::ColAddOperation{T}) where T <: Integer

Restores a ColAddOperation. Therefore, it is the row(!) operation modeled by
the same matrix.
"""
function restore!(A::AbstractArray{T,2}, addOp::ColAddOperation{T}) where T <: Integer
	i,j,a = opData(addOp)
	A[j,:] += a*A[i,:]
end

"""
    function restore!(A::Array{T,2}, addOp::ColAddOperation{T}) where T <: Integer

Restores a coreduction two-reduction. Start is
     i     j
[****|*****|****
 ****|*****|****
 --- x --- y ---
 ****|*****|****
 ****|*****|****]

The reduction had to be performed with
discardColValue*keepCol - keepColValue*DiscardCol
"""
function restore!(A::AbstractArray{T,2}, crOp::CoreductionOperation{T}) where T <: Integer
	keepColIndex, keepColValue, discardColIndex, discardColValue = opData(crOp)

	A[discardColIndex,:] = A[keepColIndex,:]*(-keepColValue)
	A[keepColIndex,:] *= discardColValue
end

function opData(addOp::ColAddOperation)
	return addOp.i, addOp.j, addOp.a
end

function opData(br::ReductionOperation)
	return br.row, br.combination
end

function opData(cr::CoreductionOperation)
	return cr.keepColIndex, cr.keepColValue, cr.discardColIndex, cr.discardColValue
end


"""
This matrix type assumes that the number of  entries in each row never increases.

Because of this, we have
- getindex for a row in log(nnz)
- getindex for a column in #(nnz in column)
TODO!
"""
struct CoredSparseMatrix <: AbstractMatrix{Int}
	m::Int # Number of rows
	n::Int # Number of columns

	# Contains 0 if column has no stored entry, index of first stored entry in column in collnk.
	colinit::Vector{Int} 

	# Contains the index of the next stored value in the column or 0 if there is no next index.
	collnk::Vector{Int}

	# Contains the column index of the cooresponding entry in nzval.
	# Note: colval[rowptr[j]:(rowptr[j]-1)] is not necessarily sorted!
	colval::Vector{Int}

	# Row j is in rowptr[j]:(rowptr[j]-1)
	rowptr::Vector{Int}

	# Contains the row index of the cooresponding entry in nzval
	rowval::Vector{Int}

	# Stored values
	nzval::Vector{Int}
end

function CoredSparseMatrix(A::SparseMatrixCSC{Int,Int})
	m,n = size(A)
	nnzA = nnz(A)
	
	colinit = zeros(Int,n)
	collnk = zeros(Int, nnzA)
	colval = zeros(Int,nnzA)
	rowptr = zeros(Int,m+1)
	rowval = zeros(Int,nnzA)
	nzval = zeros(Int,nnzA)

	rowcounter = zeros(Int,m)
	for j in findnz(A)[1]
		rowcounter[j] += 1
	end
	rowptr[1] = 1
	for i = 2:length(rowptr)
		rowptr[i] = rowptr[i-1] + rowcounter[i-1]
	end

	# decrease rowcounter when an entry is added to the row
	lastind = 0
	for (i,j,v) in zip(findnz(A)...)
		ind = rowptr[i+1] - rowcounter[i]
		rowcounter[i] -= 1

		if colinit[j] == 0
			colinit[j] = ind
			lastind = ind
		else
			collnk[lastind] = ind
			lastind = ind
		end
		rowval[ind] = i
		colval[ind] = j
		nzval[ind] = v
	end
	return CoredSparseMatrix(m,n,colinit,collnk,colval,rowptr,rowval,nzval)
end

getcolinit(A::CoredSparseMatrix) = getfield(A, :colinit)
getcollnk(A::CoredSparseMatrix) = getfield(A, :collnk)
getcolval(A::CoredSparseMatrix) = getfield(A, :colval)
getrowptr(A::CoredSparseMatrix) = getfield(A, :rowptr)
getrowval(A::CoredSparseMatrix) = getfield(A, :rowval)
getnzval(A::CoredSparseMatrix) = getfield(A, :nzval)

# methods for AbstractArray
Base.size(A::CoredSparseMatrix) = (getfield(A,:m), getfield(A,:n))

function Base.getindex(A::CoredSparseMatrix, i::Int, j::Int)
	@boundscheck checkbounds(A, i, j)
	ind = getcolinit(A)[j]
	if ind == 0
		return 0
	end
	collnk = getcollnk(A)
	rowval = getrowval(A)
	nzval = getnzval(A)

	currowval = rowval[ind]
	while 0 < ind && currowval <= i
		currowval = rowval[ind]
		if currowval == i
			return nzval[ind]
		end
		ind = collnk[ind]
	end
	return 0
end

Base.getindex(A::CoredSparseMatrix, I::Tuple{Int,Int}) = getindex(A, I[1], I[2])

function getindex_col(A::CoredSparseMatrix, j)
	@boundscheck checkbounds(A, :, j)
	ind = getcolinit(A)[j]
	if ind == 0
		return spzeros(Int, size(A,1))
	end
	ptrnew = Int[]
	nzvalnew = Int[]

	collnk = getcollnk(A)
	rowval = getrowval(A)
	nzval = getnzval(A)
	while 0 != ind
		push!(ptrnew, rowval[ind])
		push!(nzvalnew, nzval[ind])
		ind = collnk[ind]
	end
	return sparsevec(ptrnew, nzvalnew, size(A,1))
end

function getindex_row(A::CoredSparseMatrix, i)
	@boundscheck checkbounds(A, i, :)
	rowptr = getrowptr(A)
	rng = rowptr[i]:(rowptr[i+1]-1)

	colval = getcolval(A)
	nzval = getnzval(A)
	return sparsevec(colval[rng], nzval[rng], size(A,2))
end

Base.getindex(A::CoredSparseMatrix, i, ::Colon) = getindex_row(A, i)
Base.getindex(A::CoredSparseMatrix, ::Colon, j) = getindex_col(A, j)

Base.IndexStyle(::CoredSparseMatrix) = IndexCartesian()

function Base.setindex!(A::CoredSparseMatrix, v::Int, i::Int, j::Int)
	if !((1 <= i <= size(A, 1)) & (1 <= j <= size(A, 2)))
        throw(BoundsError(A, (i,j)))
    end
	colinit = getcolinit(A)

	collnk = getcollnk(A)
	rowval = getrowval(A)
	
	ind	= colinit[j]
	if ind != 0 && rowval[ind] > i
		ind = 0
	end
	if ind != 0
		while ((indnext = collnk[ind]) != 0) && rowval[indnext] <= i
			ind = collnk[ind]
		end
	end
	# now ind is either 0 or the largest index such that rowval[ind] <= i
	nzval = getnzval(A)
	if ind != 0 && rowval[ind] == i
		# there is a value at the position, simply overwrite.
		nzval[ind] = v
		return
	elseif v == 0
		# there is no value at (i,j) therefore we can return
		return
	end

	# now either ind = 0, rowval[ind] < i
	
	# find insert position
	rowrngi = rowrng(A::CoredSparseMatrix, i)
	if isempty(rowrngi)
		error("Trying to set value to empty row!")
	end
	firstZeroInd = findfirst(==(0), nzval[rowrngi])
	if firstZeroInd === nothing
		error("Row has too many entries already!")
	end
	insertind = first(rowrngi) + firstZeroInd - 1

	# modify data structure according to input
	# 1. remove reference to insertind from previous column if necessary
	colval = getcolval(A)
	oldcolval = colval[insertind]

	if oldcolval != 0
		# there is a reference from a column which needs to be removed
		oldcolind = colinit[oldcolval]
		if oldcolind == insertind
			colinit[oldcolval] = collnk[insertind]
		else
			while (nextind = collnk[oldcolind]) != 0 && nextind != insertind
				oldcolind = nextind
			end
			# now collnk[oldcolind] == insertind, now modify it to the next entry
			collnk[oldcolind] = collnk[insertind]
		end
	end

	# 2. make sure new column references new entry
	if ind == 0
		if colinit[j] != 0
			collnk[insertind] = colinit[j]
		else
			collnk[insertind] = 0
		end
		colinit[j] = insertind
	else
		collnk[insertind] = collnk[ind]
		collnk[ind] = insertind
	end


	# insert value
	getrowval(A)[insertind] = i
	getcolval(A)[insertind] = j
	nzval[insertind] = v

	return v
end

function rowrng(A::CoredSparseMatrix, i)
	if !(1 <= i <= size(A,1))
		throw(BoundsError(A, (i,:)))
	end
	rowptr = getrowptr(A)
	return rowptr[i]:(rowptr[i+1]-1)
end

function SparseMatrixCSC(A::CoredSparseMatrix)
	m,n = size(A)
	colinit = getcolinit(A)
	collnk = getcollnk(A)
	rowval = getrowval(A)
	nzval = getnzval(A)
	nnz = length(collnk)

	if (m == 0 || n == 0) || nnz == 0
		return spzeros(Int,m,n)
	end

	rowvalnew = Vector{Int}(undef, nnz)
	colvalnew = Vector{Int}(undef, nnz)
	nzvalnew = Vector{Int}(undef, nnz)

	counter = 1
	for j = 1:n
		ind = colinit[j]
		while ind != 0
			rowvalnew[counter] = rowval[ind]
			colvalnew[counter] = j
			nzvalnew[counter] = nzval[ind]

			ind = collnk[ind]
			counter += 1
		end
	end

	return sparse(rowvalnew, colvalnew, nzvalnew, m, n)
end

# code for Coreductions
function apply!(A::CoredSparseMatrix, crOp::CoreductionOperation{T}) where T <: Integer
	keepColIndex, keepColValue, discardColIndex, discardColValue = opData(crOp)
	# The following is implemented at the level of the sparse matrix representation
	#A[:,keepColIndex] = discardColValue*A[:,keepColIndex] - keepColValue*A[:,discardColIndex]
	#A[:,discardColIndex] = spzeros(T, size(A,1))

	colinit = getcolinit(A)
	collnk = getcollnk(A)
	colval = getcolval(A)
	rowval = getrowval(A)
	nzval = getnzval(A)

	indKeep, indDiscard = colinit[keepColIndex], colinit[discardColIndex]
	prevIndKeep, prevIndDiscard = 0, 0

	while !(indKeep == 0 && indDiscard == 0)
		if indDiscard == 0
			nzval[indKeep] *= discardColValue
			prevIndKeep = indKeep
			indKeep = collnk[indKeep]
		elseif indKeep == 0 || rowval[indDiscard] < rowval[indKeep]
			# move entry at indDiscard from colDiscard to colKeep
			# need to ensure the references
			# 1. prevIndKeep -a-> indDiscard -b-> indKeep (and c: prevIndKeep := indDiscard)
			# 2. prevIndDiscard -a-> collnk[indDiscard] (and b: indDiscard := collnk[indDiscard])

			# fix arrow 1a
			if prevIndKeep == 0
				colinit[keepColIndex] = indDiscard
			else
				collnk[prevIndKeep] = indDiscard
			end

			# fix arrow 2a
			if prevIndDiscard == 0
				colinit[discardColIndex] = collnk[indDiscard]
			else
				collnk[prevIndDiscard] = collnk[indDiscard]
			end

			# fix 1c
			prevIndKeep = indDiscard # fix 1c. Note: after the next line, indDiscard is in prevIndKeep
			# fix arrow 2b.
			indDiscard = collnk[indDiscard]
			# fix 1b: since variable indDiscard changed, need to use prevIndKeep
			collnk[prevIndKeep] = indKeep
			
			colval[prevIndKeep] = keepColIndex
			nzval[prevIndKeep] = -keepColValue*nzval[prevIndKeep]
		elseif rowval[indDiscard] == rowval[indKeep]
			nzval[indKeep] = discardColValue*nzval[indKeep] - keepColValue*nzval[indDiscard]
			nzval[indDiscard] = 0
			prevIndDiscard = indDiscard
			indDiscard = collnk[indDiscard]
			prevIndKeep = indKeep
			indKeep = collnk[indKeep]
		else
			@assert rowval[indDiscard] > rowval[indKeep]
			# modify
			nzval[indKeep] *= discardColValue
			prevIndKeep = indKeep
			indKeep = collnk[indKeep]
		end
#		show(colinit')
#		show(collnk')
#		show(rowval)
#		show(nzval)
	end

	return A
end