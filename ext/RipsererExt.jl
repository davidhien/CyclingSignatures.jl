module RipsererExt

using CyclingSignatures
using Distances: pairwise
using Ripserer

function __init__()
    @warn "Ripserer currently has a bug which makes this method return incorrect representatives. Consider using Val(:RipsererNoThreshold) or Val(:RipsererManualReconstruct) instead."
end

_ripserer_field(::Type{CyclingSignatures.FF{M}}) where {M} = Ripserer.Mod{M}
_ripserer_field(::Type{Ripserer.Mod{M}}) where {M} = Ripserer.Mod{M}
_ripserer_field(field::Type) = field

_coefficient_field(::Type{CyclingSignatures.FF{M}}) where {M} = CyclingSignatures.FF{M}
_coefficient_field(::Type{Ripserer.Mod{M}}) where {M} = Ripserer.Mod{M}
_coefficient_field(field::Type) = field

function CyclingSignatures.trajectory_barcode(::Val{:Ripserer}, points, metric, flt_threshold, field)
    return ripserer_involuted_threshold_bars(points, metric, flt_threshold, field)
end

function CyclingSignatures.trajectory_barcode(
    ::Val{:RipsererNoThreshold},
    points,
    metric,
    flt_threshold,
    field,
)
    return ripserer_involuted_no_threshold_bars(points, metric, flt_threshold, field)
end

function CyclingSignatures.trajectory_barcode(
    ::Val{:RipsererManualReconstruct},
    points,
    metric,
    flt_threshold,
    field,
)
    return ripserer_manual_reconstruct_bars(points, metric, flt_threshold, field)
end

function ripserer_involuted_threshold_bars(points, metric, flt_threshold, field)
    d_mat = ripserer_distance_matrix(points, metric)
    filtration = Ripserer.Rips(d_mat; threshold=flt_threshold)
    diagram = Ripserer.ripserer(
        filtration;
        dim_max=1,
        alg=:involuted,
        reps=true,
        field=_ripserer_field(field),
    )[2]

    return ripserer_representative_bars(
        active_intervals(diagram, flt_threshold),
        _coefficient_field(field),
    )
end

function ripserer_involuted_no_threshold_bars(points, metric, flt_threshold, field)
    d_mat = ripserer_distance_matrix(points, metric)
    filtration = Ripserer.Rips(d_mat)
    diagram = Ripserer.ripserer(
        filtration;
        dim_max=1,
        alg=:involuted,
        reps=true,
        field=_ripserer_field(field),
    )[2]

    return ripserer_representative_bars(
        active_intervals(diagram, flt_threshold),
        _coefficient_field(field);
        force_infinite_death=true,
    )
end

function ripserer_manual_reconstruct_bars(points, metric, flt_threshold, field)
    d_mat = ripserer_distance_matrix(points, metric)
    filtration = Ripserer.Rips(d_mat; threshold=flt_threshold)
    diagram = Ripserer.ripserer(filtration; dim_max=1, reps=true, field=_ripserer_field(field))[2]

    return ripserer_reconstructed_bars(
        active_intervals(diagram, flt_threshold),
        filtration,
        flt_threshold,
        _coefficient_field(field),
    )
end

function ripserer_distance_matrix(points, metric)
    d_mat = pairwise(metric, points; dims=2)
    edge_zero = eps(float(eltype(d_mat)))
    for j in axes(d_mat, 2), i in (j + 1):lastindex(axes(d_mat, 1))
        if iszero(d_mat[i, j])
            d_mat[i, j] = edge_zero
            d_mat[j, i] = edge_zero
        end
    end
    return d_mat
end

function active_intervals(diagram, flt_threshold)
    return filter(
        bar -> Ripserer.birth(bar) <= flt_threshold && Ripserer.death(bar) >= flt_threshold,
        diagram,
    )
end

function ripserer_representative_bars(diagram, field; force_infinite_death=false)
    return map(diagram) do interval
        representative = Ripserer.representative(interval)
        simplex_list, coeff_list = ripserer_representative_to_edge_data(representative, field)
        death = force_infinite_death ? Inf : Ripserer.death(interval)
        return CyclingSignatures.TrajectoryBar(
            Ripserer.birth(interval),
            death,
            simplex_list,
            coeff_list,
        )
    end
end

function ripserer_reconstructed_bars(diagram, filtration, flt_threshold, field)
    return map(diagram) do interval
        cycle = Ripserer.reconstruct_cycle(filtration, interval, flt_threshold)
        simplex_list, coeff_list = ripserer_reconstructed_cycle_to_edge_data(cycle, field)
        return CyclingSignatures.TrajectoryBar(
            Ripserer.birth(interval),
            Ripserer.death(interval),
            simplex_list,
            coeff_list,
        )
    end
end

function ripserer_reconstructed_cycle_to_edge_data(cycle, field)
    simplex_list = [Tuple(collect(Ripserer.vertices(simplex))) for simplex in cycle]
    coeff_list = fill(one(field), length(simplex_list))
    if !isempty(coeff_list)
        coeff_list[1] = -one(field)
    end
    return simplex_list, coeff_list
end

function ripserer_representative_to_edge_data(representative, field)
    edge_coefficients = Dict{Tuple{Int,Int},field}()
    for (simplex, coefficient) in representative
        edge = Tuple(collect(Ripserer.vertices(simplex)))
        value = field(Int(coefficient))
        edge_coefficients[edge] = get(edge_coefficients, edge, zero(field)) + value
    end

    simplex_list = collect(keys(edge_coefficients))
    coeff_list = [edge_coefficients[simplex] for simplex in simplex_list]
    keep = findall(!iszero, coeff_list)

    return simplex_list[keep], coeff_list[keep]
end

end
