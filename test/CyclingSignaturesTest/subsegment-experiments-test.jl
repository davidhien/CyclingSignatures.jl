@testset "Simple tests" begin
    segs = sampleSegments(SubsegmentSampleParameter(5,5),1:5)

    @test all(==(1:5),segs)
end
