using CyclingSignatures
using Test
import CyclingSignatures: compute_pivots

@testset "basic_reduction!" begin
    M = Rational[1 1 1]
    basic_reduction!(M)
    @test M == Rational[1 0 0]

    M = Float64[1 1 1]
    basic_reduction!(M)
    @test M == Float64[1 0 0]

    M = FF{5}[1 1 1]
    basic_reduction!(M)
    @test M == FF{5}[1 0 0]

    # TODO: add some larger examples with known result

    # test with random matrices
    for _ = 1:10
        M = FF{7}.(rand(1:10, 5, 5))
        N = copy(M)
        basic_reduction!(M)
        
        # test that nonzero pivots are unique
        pivots = compute_pivots(M)
        nz_pivots = filter(!iszero,pivots)
        @test length(nz_pivots) == length(unique(nz_pivots))

        # test that 
        M2 = [M N]
        basic_reduction!(M2)
        @test M2[:,1:5] == M
        @test all(==(0), M2[:,6:end]) 
    end
end

@testset "colspace_normal_form(M)" begin
    M = Rational[1 1; 1 0]
    N = colspace_normal_form(M)
    @test N == Rational[1 0; 0 1]

    M = FF{7}[1 1; 1 0]
    N = colspace_normal_form(M)
    @test N == Rational[1 0; 0 1]

    M = FF{7}[1 1 4; 2 0 0]
    N = colspace_normal_form(M)
    @test N == Rational[1 0; 0 1]

    M = FF{7}[1 1 4; 0 0 0; 2 0 0]
    N = colspace_normal_form(M)
    @test N == Rational[1 0; 0 0; 0 1]
end
