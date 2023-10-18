function circleTimeSeries(nTurns=1)
    # TODO: only one helper method for this
    a = 0:.05:2*pi*nTurns
    return [cos.(a)'; sin.(a)']
end

function plotCircleTimeSeries()
    pts = circleTimeSeries(nTurns=1)
    pts_Q = quantize(pts, .1)
    plt = scatter(10*pts)
    scatter!(plt.axis, pts_Q)
    return plt
end