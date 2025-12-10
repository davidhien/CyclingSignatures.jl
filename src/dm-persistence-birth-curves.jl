struct CyclingBirthCurve{T, S <: AbstractVector{T}}
    birth_curve::Vector{Tuple{Int,Int,Float64}} # (i,j,dist)
    inc_rep::S # for one representative
end

function birth_curves(trajectory_space::TrajectorySpace, range, r_max=nothing; field=DEFAULT_FIELD)
    traj = get_trajectory(trajectory_space)
    traj_pts =  evaluate_interval(traj, first(range), last(range))

    metric = get_metric(trajectory_space)
    threshold = r_max === nothing ? trajectory_space.flt_max_heuristic : r_max

    dm_birth_curves = dm_component_birth_curves(traj_pts, metric, threshold; field=field)

    birth_curves = CyclingBirthCurve{field,Vector{field}}[]
    for bc in dm_birth_curves
        inc_rep = map_cycle(get_comparison_space(trajectory_space), traj_pts, bc.simplex_list, bc.coeff_list)

        if !isempty(inc_rep) && !all(==(0), inc_rep)
            inc_rep_normal_form = Vector(inc_rep) # colspace_normal_form(Vector(inc_rep)[:,:])
            push!(birth_curves, CyclingBirthCurve(bc.birth_curve, inc_rep_normal_form[:] ))
        end
    end
    # reindexing according to time series
    #t_vec_segment = t_vec_segment(traj, first(range), last(range))
    #t_vec_segment -= t_vec_segment[1]+1

    # combine birth curves with the same representative
    birth_curves = reduce_birth_curves(birth_curves)

    return birth_curves #reindex_birth_curves(bc_reduced, t_vec)
end

#function reindex_birth_curves(birth_curves, t_vec_segment)
#    t_vec_segment
#end

function reduce_birth_curves(birth_curves::Vector{CyclingBirthCurve{T,S}}) where {T, S <: AbstractVector{T}}
    combined_births = Dict{Vector{Int}, Vector{Tuple{Int,Int,Float64}}}()

    for bc in birth_curves
        in_rep_int = Int.(bc.inc_rep)
        if haskey(combined_births, in_rep_int)
            append!(combined_births[in_rep_int], bc.birth_curve)
        else
            combined_births[in_rep_int] = copy(bc.birth_curve)
        end
    end

    return map(collect(keys(combined_births))) do k
        red_births = reduce_births(combined_births[k])
        return CyclingBirthCurve{T,S}(red_births, T.(k))
    end
end

function reduce_births(vec::Vector{Tuple{Int,Int,Float64}})
    vec_new = Tuple{Int,Int,Float64}[]
    for t in vec
        i, j, r = t
        isnew = true
        for t_test in vec
            if t_test[1] >= i && t_test[2] <= j && t_test[3] <= r
                if t_test != t
                    isnew = false
                    break
                end
            end
        end
        isnew && push!(vec_new, t)
    end
    return vec_new
end


struct DMBirthCurve{T<:Integer,N}
    birth_curve::Vector{Tuple{Int,Int,Float64}} # (i,j,dist)
    simplex_list::Vector{NTuple{N,Int}} # for one representative
    coeff_list::Vector{T} # for one representative
end

function dm_component_birth_curves(points, metric, threshold; field=DEFAULT_FIELD)
    # can reuse first pass
    cc, cc_labels = dm_components_first_pass_implicit(points, metric, threshold)

    # second pass: get representatives
    birth_curves = dm_component_birth_curves(cc, cc_labels, points, metric)

    cycling_birth_curves = map(birth_curves) do bc
        i, j, _ = bc[1] # generate generator for one representative
        edges, coeffs = curve_cycle(i, j; F=field)
        return DMBirthCurve(bc, edges, coeffs)
    end

    return cycling_birth_curves
end

function dm_component_birth_curves(cc, cc_labels::Matrix{Int}, points, metric)
    birth_curves = Dict{Int,Vector{Tuple{Int,Int,Float64}}}()
    n = size(points, 2)
    @inbounds for j = 1:n
        for i = 1:j
            if cc_labels[i, j] > 0
                tmp_label = cc_labels[i, j]
                label = find_root!(cc, tmp_label)
                cur_metric = metric((@view points[:, i]), (@view points[:, j]))
                # now fix representant of component label
                if !haskey(birth_curves, label)
                    birth_curves[label] = [(i, j, cur_metric)]
                else
                    # check if we have new
                    comparisons = map(birth_curves[label]) do x
                        # if all three of these are met, x less than the current birth
                        # if none of these are met, it is greater than the current birth
                        return (x[1] > i) + (x[2] < j) + (x[3] < cur_metric)
                    end
                    # if any is 3, we do not add
                    # otherwise, we add but remove all that are dominated
                    have_new = !any(==(3), comparisons)

                    counter = 1
                    while counter <= length(comparisons)
                        if comparisons[counter] == 0
                            deleteat!(birth_curves[label], counter)
                            deleteat!(comparisons, counter)
                        else
                            counter += 1
                        end
                    end
                    if have_new
                        push!(birth_curves[label], (i, j, cur_metric))
                    end
                end
            end
        end
    end
    return values(birth_curves)
end
