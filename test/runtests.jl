using Test, SafeTestsets
using CyclingSignatures


@testset "Test" begin
    @safetestset "ff" begin include("ff.jl") end
    @safetestset "distance-matrix-persistence" begin include("distance-matrix-persistence.jl") end
    @safetestset "lin-alg-util" begin include("lin-alg-util.jl") end
    @safetestset "dynamic-distance.jl" begin include("dynamic-distance.jl") end
end
