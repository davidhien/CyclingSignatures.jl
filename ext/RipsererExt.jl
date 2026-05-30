module RipsererExt

using CyclingSignatures
using Distances: pairwise
using Ripserer

_ripserer_field(::Type{CyclingSignatures.FF{M}}) where {M} = Ripserer.Mod{M}
_ripserer_field(::Type{Ripserer.Mod{M}}) where {M} = Ripserer.Mod{M}
_ripserer_field(field::Type) = field

_coefficient_field(::Type{CyclingSignatures.FF{M}}) where {M} = CyclingSignatures.FF{M}
_coefficient_field(::Type{Ripserer.Mod{M}}) where {M} = Ripserer.Mod{M}
_coefficient_field(field::Type) = field

function CyclingSignatures.trajectory_barcode(::Val{:Ripserer}, points, metric, flt_threshold, field)
    d_mat = ripserer_distance_matrix(points, metric)
    filtration = Ripserer.Rips(d_mat)
    diagram = Ripserer.ripserer(
        filtration;
        dim_max=1,
        alg=:involuted,
        reps=true,
        field=_ripserer_field(field),
    )[2]

    filtered_diagram = filter(bar -> bar.birth <= flt_threshold, diagram)
    filtered_diagram = filter(bar -> bar.death >= flt_threshold, filtered_diagram)

    return ripserer_diagram_to_trajectory_bars(filtered_diagram, _coefficient_field(field))
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

function ripserer_diagram_to_trajectory_bars(diagram, field)
    return map(diagram) do interval
        simplex_list, coeff_list = ripserer_representative_to_edge_data(
            Ripserer.representative(interval),
            field,
        )
        return CyclingSignatures.TrajectoryBar(
            Ripserer.birth(interval),
            Inf,
            simplex_list,
            coeff_list,
        )
    end
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
