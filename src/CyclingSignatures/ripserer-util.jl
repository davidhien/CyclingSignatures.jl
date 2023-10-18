function trajectoryBarcode(::Val{:Ripserer}, trajPoints, metric, fltThreshold, field=Mod{2})
    dm_traj = pairwise(metric, trajPoints)
    fltTraj = Rips(dm_traj, threshold=fltThreshold)
    trajDiag = ripserer(fltTraj, alg=:involuted, reps=true, field=field)[2]

    return ripsererDiagToTrajectoryBars(trajDiag)
end

function ripsererDiagToTrajectoryBars(trajDiag)
    return map(trajDiag) do int
        rep = representative(int)
        simplex_list, coeff_list = ripsererRepToEdgeData(rep)
        return TrajectoryBar(int.birth,int.death,simplex_list,coeff_list)
    end
end

function ripsererDiagToDiagWithRepLists(trajDiag)
    newIntervals = map(trajDiag) do int
        rep = representative(int)
        simplex_list, coeff_list = ripsererRepToEdgeData(rep)
        return PersistenceInterval(int,simplex_list=simplex_list, coeff_list=coeff_list)
    end
    return PersistenceDiagram(newIntervals)
end

function ripsererRepToEdgeData(rep::Ripserer.Chain{F,G}) where {F,G}
    rep_dict = ripsererRepToEdgeDict(rep)
    simplex_list = collect(keys(rep_dict))
    coeff_list = map(s-> rep_dict[s], simplex_list)
    return simplex_list, coeff_list
end

function ripsererRepToEdgeDict(rep::Ripserer.Chain{F,G}) where {F,G}
    dict = Dict{Tuple{Int,Int},F}()
    for (k,v) in rep
        d_key = Tuple(collect(Ripserer.vertices(k)))
        if !haskey(dict, d_key)
            dict[d_key] = v
        else
            dict[d_key] += v
        end
    end

    return dict
end