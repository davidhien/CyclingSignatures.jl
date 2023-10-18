function testSBDistance1()
    dist = SBDistance(1,2)
    x1 = [4;0]
    y1 = [0;0]
    if !isapprox(dist(x1,y1),4)
        return false
    end
    return true
end

function testSBDistance2()
    dist = SBDistance(1,2)
    x2 = [4;0]
    y2 = [0;3]
    if !isapprox(dist(x2,y2),6)
        return false
    end
    return true
end

function testSBDistancePairwise()
    X = [7 1; 4 0]
    dist = SBDistance(1,2)
    
    return isapprox(pairwise(dist, X),[0 8; 8 0])
end

