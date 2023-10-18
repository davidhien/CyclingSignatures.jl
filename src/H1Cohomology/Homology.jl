function cohomologyGenerators(c::AbstractChainComplex, k)
    # TODO: overload getindex for cubicalcomplex
    D0 = Matrix(coboundaryMatrix(c, k-1))
    D1 = Matrix(coboundaryMatrix(c, k))

    return cohomologyGenerators(D0,D1)
end

"""
    function cohomologyGenerators(D0, D1)
Compute generators for cohomology.
"""
function cohomologyGenerators(D0, D1)
    S0 = smith(D0)
    S1 = smith(D1)

    r0 = findlast(x-> x != 0, S0.SNF)  # rank of im(D0)
    r1 = findfirst(x-> x == 0, S1.SNF) # rank of ker(D1)

    if r1 === nothing
        r1 = size(D1,2) - length(S1.SNF)
    else
        r1 = size(D1,2) - r1 + 1
    end


    if r1 == 0
        return spzeros(Int, size(D0,2), 0)
    end

    kerD1 = S1.Tinv[:,end-r1+1:end]
    Ainv = S1.T[end-r1+1:end,:]

    if r0 === nothing
        return sparse(kerD1)
    end

    imD0 = S0.S[:,1:r0]
    im_in_ker_basis = Ainv*imD0

    S2 = smith(im_in_ker_basis)

    r2 = findlast(x-> x != 0, S2.SNF) # rank of im(D0)

    if r2 === nothing
        r2 = 0 # TODO: this should never happen, maybe make this an exception
    end

    return sparse(kerD1*S2.S[:,r2+1:end])
end
