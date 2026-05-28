using Test, SafeTestsets
using CyclingSignatures


@testset "Test" begin
    @safetestset "ff" begin include("ff.jl") end
    @safetestset "distance-matrix-persistence" begin include("distance-matrix-persistence.jl") end
    @safetestset "lin-alg-util" begin include("lin-alg-util.jl") end
    @safetestset "dynamic-distance" begin include("dynamic-distance.jl") end
    @safetestset "comparison-space" begin include("comparison-space.jl") end
    @safetestset "trajectory" begin include("trajectory.jl") end
    @safetestset "cycling-signature" begin include("cycling-signature.jl") end
    @safetestset "ripserer-extension" begin include("ripserer-extension.jl") end
    @safetestset "sample-tools" begin include("sample-tools.jl") end
    @safetestset "subsegment-experiments" begin include("subsegment-experiments.jl") end
    @safetestset "plotting-interface" begin include("plotting-interface.jl") end
end
