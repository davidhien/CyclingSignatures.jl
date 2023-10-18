LSQR_TOL = 1e-11

struct H1
    # Note: thus far it is implicitely assumed that the ordering of the basis
    # of H₀ is consistent with the order of the points in the set x.
    generatorMatrix::SparseMatrixCSC{Int64,Int64} # of first cohomology, rows are generators
    D0::SparseMatrixCSC{Int64,Int64} # first COboundary matrix
    D1::SparseMatrixCSC{Int64,Int64} # second COboundary matrix
end

function betti_1(h::H1)
    return size(h.generatorMatrix,2)
end

function coboundaryMatrix(h::H1, i)
    if i == 0
        return h.D0
    elseif i == 1
        return h.D1
    end
    error("Not defined")
end

function generatorMatrix(h::H1)
    return h.generatorMatrix
end

"""
    function norms(h::H1)
computes the l2 norms of the generators
"""
function norms(h::H1, A=I; LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    liftMatrix = liftOfCircleValuedFunction(h, A, LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    realcocycles = h.D0 * liftMatrix' + generatorMatrix(h)*A

    return mapslices(norm,realcocycles,dims=[1])
end

function harmonicCocycles(h::H1, A=I; LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    liftMatrix = liftOfCircleValuedFunction(h, A, LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    return h.D0 * liftMatrix' + generatorMatrix(h)*A
end

"""
    function circularCoordinates(h, A)

Computes a matrix A such that A[i,:] contains the values of a circle valued
function corresponding to the i-th generator of h on the basis of C^0 of h
which was used to compute h.D0.
"""
function circularCoordinates(h::H1, A=I; LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    coord = liftOfCircleValuedFunction(h, A, LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    return mod.(coord, 1)
end

"""
    function liftOfCircleValuedFunction(h::H1, A, LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)

Computes a matrix A such that A[i,:] contains the values of a lift of the circle valued
function corresponding to the i-th generator of h on the basis of C^0 of h
which was used to compute h.D0. This is useful for example to compute the norm
of a generator.
"""
function liftOfCircleValuedFunction(h::H1, A=I; LSQR_ATOL = LSQR_TOL, LSQR_RTOL = LSQR_TOL)
    # TODO: rename
    D0 = coboundaryMatrix(h, 0)

    B = generatorMatrix(h)*A # contains generators in rows
    β1 = size(B, 2)

    coordinates = zeros(β1, size(D0, 2)) #preallocate memory

    for i = 1:β1
        cocycle = convert(Vector,B[:,i])
        ls_sol = lsqr(D0, (-1)*cocycle, atol = LSQR_ATOL, btol = LSQR_RTOL) #solve least squares
        coordinates[i,:] = ls_sol # do not force coord values in [0,1)
    end

    return coordinates
end

function liftCoordinates(coordinates, t_indices)
    return mapslices(x -> liftCoordinate(x,t_indices), coordinates, dims=[2])
end

function liftCoordinate(c, t_indices)
    # TODO: refactor, make this function accept only a t->c vector, this makes more sense
    # the values of c are assumed to be in [0,1]
    n = length(t_indices)
    newc = zeros(length(t_indices))

    newc[1] = 0
    for i = 2:n
        newc[i] = newc[i-1] + signedS1Distance(c[t_indices[i-1]], c[t_indices[i]])
    end
    return newc
end