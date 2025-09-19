function circle_time_series(subdivision, turns)
    a = LinRange(0, 2*π*turns, subdivision*turns + 1)
    return [cos.(a)'; sin.(a)']
end

function double_circle_time_series(subdivision, bits)
    a = zeros(2,0)
    template = circle_time_series(subdivision, 1)[:,1:end-1]
    l = template .+ [-1,0]
    r =  [-1,1] .*template .+ [1,0]
    for b in bits
        if b == 0
            a = hcat(a, l)
        elseif b == 1
            a = hcat(a, r)
        end
    end
    #note: l[:,1] == r[:,1]
    return [a l[:,1]]
end
