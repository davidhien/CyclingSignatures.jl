"""
Sampling and refinement helpers for trajectories and solution objects.

These utilities refine sample points so consecutive points are close in either
metric distance or adjacent-box sense, optionally combined with a direction
constraint (e.g. unit tangent vectors).
"""

function _to_vector(x)
    return x isa AbstractVector ? Vector(x) : vec(collect(x))
end

function _ds_stepper(ds)
    mod = parentmodule(typeof(ds))
    for name in (:set_state!, :step!, :get_state)
        if !isdefined(mod, name)
            throw(ArgumentError("`$(nameof(mod)).$name` is not defined; provide a DS type whose defining module implements `set_state!`, `step!`, and `get_state`."))
        end
    end
    set_state_f = getfield(mod, :set_state!)
    step_f = getfield(mod, :step!)
    get_state_f = getfield(mod, :get_state)
    return (x, Δt) -> begin
        set_state_f(ds, x)
        step_f(ds, Δt, true)
        return copy(get_state_f(ds))
    end
end

function _boxes_adjacent(x, y, r)
    qx = round.(Int, x ./ r)
    qy = round.(Int, y ./ r)
    return norm(qx - qy, Inf) <= 1
end

function _safe_unit(v, norm_f)
    nv = norm_f(v)
    if nv == 0
        return v, false
    end
    return v / nv, true
end

function _make_consistency_check(;
    r,
    mode::Symbol=:distance,
    metric=euclidean,
    pp=identity,
    direction_f=nothing,
    direction_r=r,
    direction_mode::Symbol=mode,
    direction_metric=metric,
    direction_norm=norm,
    direction_pp=identity,
)
    function space_ok(x, y)
        if mode === :distance
            return metric(pp(x), pp(y)) <= r
        elseif mode === :boxes
            return _boxes_adjacent(pp(x), pp(y), r)
        else
            throw(ArgumentError("Unknown mode: $mode"))
        end
    end

    if isnothing(direction_f)
        return space_ok
    end

    function direction_ok(x, y)
        vx = direction_f(x)
        vy = direction_f(y)
        ux, okx = _safe_unit(vx, direction_norm)
        uy, oky = _safe_unit(vy, direction_norm)
        if !(okx && oky)
            return false
        end
        ux = direction_pp(ux)
        uy = direction_pp(uy)
        if direction_mode === :distance
            return direction_metric(ux, uy) <= direction_r
        elseif direction_mode === :boxes
            return _boxes_adjacent(ux, uy, direction_r)
        else
            throw(ArgumentError("Unknown direction_mode: $direction_mode"))
        end
    end

    return (x, y) -> space_ok(x, y) && direction_ok(x, y)
end

function _resample_segment_stepper(
    ic,
    ec,
    dt,
    is_close,
    stepper;
    max_depth=512,
    verbose=false,
)
    x0 = _to_vector(ic)
    points = Vector{typeof(x0)}(undef, 1)
    points[1] = x0

    t_cur = zero(dt)
    t_end = dt
    x_cur = points[end]
    point_ct = 0
    max_points = 100 * max_depth

    while !is_close(x_cur, ec)
        point_ct += 1
        if point_ct > max_points
            @warn "Counter exhausted!"
            break
        end
        dt_rem = t_end - t_cur
        dt_step = dt_rem
        x_next = stepper(x_cur, dt_step)
        depth = 1
        while !is_close(x_cur, x_next)
            depth += 1
            if depth > max_depth
                @warn "Went over max depth"
                break
            end
            dt_step /= 2
            x_next = stepper(x_cur, dt_step)
            verbose && print(depth)
        end
        push!(points, _to_vector(x_next))
        t_cur += dt_step
        x_cur = points[end]
    end
    verbose && println()
    return reduce(hcat, points)
end

function _resample_segment_time(
    t_start,
    t_end,
    state_at,
    is_close;
    max_depth=512,
    verbose=false,
)
    x_start = _to_vector(state_at(t_start))
    points = Vector{typeof(x_start)}(undef, 1)
    points[1] = x_start
    x_end = state_at(t_end)

    t_cur = t_start
    x_cur = x_start
    point_ct = 0
    max_points = 100 * max_depth

    while !is_close(x_cur, x_end)
        point_ct += 1
        if point_ct > max_points
            @warn "Counter exhausted!"
            break
        end
        dt_rem = t_end - t_cur
        dt_step = dt_rem
        x_next = state_at(t_cur + dt_step)
        depth = 1
        while !is_close(x_cur, x_next)
            depth += 1
            if depth > max_depth
                @warn "Went over max depth"
                break
            end
            dt_step /= 2
            x_next = state_at(t_cur + dt_step)
            verbose && print(depth)
        end
        push!(points, _to_vector(x_next))
        t_cur += dt_step
        x_cur = points[end]
    end
    verbose && println()
    return reduce(hcat, points)
end

function _concat_segments(segments)
    new_lengths = map(m -> size(m, 2), segments)
    t_vec = zeros(Int, length(new_lengths) + 1)
    t_vec[1] = 1
    for i = 2:length(t_vec)
        t_vec[i] = t_vec[i-1] + new_lengths[i-1]
    end
    return reduce(hcat, segments), t_vec
end

"""
    get_interpolation_function(Y)

Build cubic-spline interpolation functions for a column-oriented trajectory
matrix `Y`. Returns `(c, dc)`, where `c(t)` evaluates the interpolated point and
`dc(t)` evaluates its derivative.
"""
function get_interpolation_function(Y)
    t_grid = collect(1:size(Y, 2))
    curves = map(eachrow(Y)) do row
        return CubicSpline(collect(row), t_grid)
    end
    return t -> map(curve -> curve(t), curves),
           t -> map(curve -> DataInterpolations.derivative(curve, t), curves)
end

"""
    interpolate_to_distance(Y, r, metric=euclidean; kwargs...)

Refine a column-oriented trajectory matrix `Y` using cubic-spline interpolation
so consecutive interpolated points are within distance `r`. Returns the refined
trajectory, the interpolated derivative samples, and `t_vec` indexing refined
points by original sample step.
"""
function interpolate_to_distance(
    Y,
    r,
    metric=euclidean;
    norm_sb=norm,
    sb_r=nothing,
    verbose=false,
    max_depth=50000,
)
    c, dc = get_interpolation_function(Y)
    point_type = eltype(c(1))
    derivative_type = eltype(dc(1))
    segments = Vector{Matrix{point_type}}(undef, size(Y, 2))
    derivative_segments = Vector{Matrix{derivative_type}}(undef, size(Y, 2))

    for i in 1:(size(Y, 2)-1)
        segments[i], derivative_segments[i] = interpolate_to_distance(
            c,
            dc,
            i,
            i + 1,
            r,
            metric;
            norm_sb=norm_sb,
            sb_r=sb_r,
            verbose=verbose,
            max_depth=max_depth,
        )
    end

    segments[end] = convert(Matrix{point_type}, reshape(_to_vector(Y[:, end]), :, 1))
    derivative_segments[end] = convert(
        Matrix{derivative_type},
        reshape(_to_vector(dc(size(Y, 2))), :, 1),
    )
    Y_new, t_vec = _concat_segments(segments)
    derivatives_new = reduce(hcat, derivative_segments)
    return Y_new, derivatives_new, t_vec
end

function interpolate_to_distance(
    c,
    dc,
    a,
    b,
    r,
    metric=euclidean;
    norm_sb=norm,
    sb_r=nothing,
    verbose=false,
    max_depth=50000,
)
    function space_consistent(x, y)
        return metric(x, y) <= r
    end

    function direction_consistent(x, y)
        if isnothing(sb_r)
            return true
        end
        ux, okx = _safe_unit(x, norm_sb)
        uy, oky = _safe_unit(y, norm_sb)
        return okx && oky && norm_sb(ux - uy) <= sb_r
    end

    t_values = [float(a)]
    target = c(b)
    target_derivative = dc(b)
    counter = 1

    while !(
        space_consistent(c(t_values[end]), target) &&
        direction_consistent(dc(t_values[end]), target_derivative)
    ) && (counter += 1) < max_depth
        push!(t_values, float(b))
        depth = 2
        while !(
            space_consistent(c(t_values[end-1]), c(t_values[end])) &&
            direction_consistent(dc(t_values[end-1]), dc(t_values[end]))
        )
            t_values[end] = t_values[end-1] + (b - t_values[end-1]) / depth
            depth *= 2
            if depth > max_depth
                @warn "Went over max depth"
                break
            end
            verbose && print(depth)
        end
    end

    if counter >= max_depth
        @warn "Counter exhausted!"
    end
    verbose && println()

    points = reduce(hcat, (_to_vector(c(t)) for t in t_values))
    derivatives = reduce(hcat, (_to_vector(dc(t)) for t in t_values))
    return points, derivatives
end

getInterpolationFunction(Y) = get_interpolation_function(Y)

function interpolateToDistance(Y, r, metric=euclidean; kwargs...)
    return interpolate_to_distance(Y, r, metric; kwargs...)
end

function interpolateToDistance(c, dc, a, b, r, metric=euclidean; kwargs...)
    return interpolate_to_distance(c, dc, a, b, r, metric; kwargs...)
end

"""
    resample_to_distance(ds, dt, Y, r; metric=euclidean, kwargs...)

Refines the time series `Y` so consecutive points are close according to a
distance or adjacency criterion. Returns the refined trajectory and `t_vec`
indexing refined points by original step.
"""
function resample_to_distance(
    ds,
    dt::Real,
    Y,
    r,
    ;
    metric=euclidean,
    pp=identity,
    mode=:distance,
    direction_f=nothing,
    direction_r=r,
    direction_metric=metric,
    direction_norm=norm,
    direction_pp=identity,
    norm_sb=norm,
    sb_r=nothing,
    sb_fct=nothing,
    verbose=false,
    max_depth=512,
)
    if isnothing(direction_f) && !isnothing(sb_fct) && !isnothing(sb_r)
        direction_f = sb_fct
        direction_r = sb_r
        direction_norm = norm_sb
        direction_metric = (x, y) -> direction_norm(x - y)
    end
    is_close = _make_consistency_check(
        r=r,
        mode=mode,
        metric=metric,
        pp=pp,
        direction_f=direction_f,
        direction_r=direction_r,
        direction_mode=mode,
        direction_metric=direction_metric,
        direction_norm=direction_norm,
        direction_pp=direction_pp,
    )

    n = size(Y, 2)
    segments = Vector{Matrix{eltype(Y)}}(undef, n)
    @views for i in 1:(n-1)
        segments[i] = resample_to_distance(
            ds,
            dt,
            Y[:, i],
            Y[:, i+1],
            r;
            metric=metric,
            pp=pp,
            mode=mode,
            direction_f=direction_f,
            direction_r=direction_r,
            direction_metric=direction_metric,
            direction_norm=direction_norm,
            direction_pp=direction_pp,
            norm_sb=direction_norm,
            sb_r=direction_r,
            sb_fct=direction_f,
            verbose=verbose,
            max_depth=max_depth,
        )
    end
    @views segments[end] = reshape(_to_vector(Y[:, end]), :, 1)
    return _concat_segments(segments)
end

function resample_to_distance(
    ds,
    dt::Real,
    ic,
    ec,
    r,
    ;
    metric=euclidean,
    pp=identity,
    mode=:distance,
    direction_f=nothing,
    direction_r=r,
    direction_metric=metric,
    direction_norm=norm,
    direction_pp=identity,
    norm_sb=norm,
    sb_r=nothing,
    sb_fct=nothing,
    verbose=false,
    max_depth=512,
)
    if isnothing(direction_f) && !isnothing(sb_fct) && !isnothing(sb_r)
        direction_f = sb_fct
        direction_r = sb_r
        direction_norm = norm_sb
        direction_metric = (x, y) -> direction_norm(x - y)
    end
    is_close = _make_consistency_check(
        r=r,
        mode=mode,
        metric=metric,
        pp=pp,
        direction_f=direction_f,
        direction_r=direction_r,
        direction_mode=mode,
        direction_metric=direction_metric,
        direction_norm=direction_norm,
        direction_pp=direction_pp,
    )
    stepper = _ds_stepper(ds)
    return _resample_segment_stepper(
        _to_vector(ic),
        _to_vector(ec),
        dt,
        is_close,
        stepper;
        max_depth=max_depth,
        verbose=verbose,
    )
end

"""
    resample_to_distance(sol, dt, r; kwargs...)

Refines a solution object `sol` by evaluating it at refined times.
"""
function resample_to_distance(
    sol,
    dt::Real,
    r,
    ;
    metric=euclidean,
    t0=sol.t[1],
    t_max=sol.t[end],
    pp=identity,
    mode=:distance,
    direction_f=nothing,
    direction_r=r,
    direction_metric=metric,
    direction_norm=norm,
    direction_pp=identity,
    verbose=false,
    max_depth=512,
)
    is_close = _make_consistency_check(
        r=r,
        mode=mode,
        metric=metric,
        pp=pp,
        direction_f=direction_f,
        direction_r=direction_r,
        direction_mode=mode,
        direction_metric=direction_metric,
        direction_norm=direction_norm,
        direction_pp=direction_pp,
    )
    state_at = t -> _to_vector(sol(t))

    if t_max <= t0
        X = reshape(state_at(t0), :, 1)
        return X, [1, 2]
    end

    n_intervals = ceil(Int, (t_max - t0) / dt)
    segments = Vector{Matrix{eltype(state_at(t0))}}(undef, n_intervals + 1)
    for i in 1:n_intervals
        t_start = t0 + (i - 1) * dt
        t_end = min(t_start + dt, t_max)
        segments[i] = _resample_segment_time(
            t_start,
            t_end,
            state_at,
            is_close;
            max_depth=max_depth,
            verbose=verbose,
        )
    end
    segments[end] = reshape(state_at(t_max), :, 1)
    return _concat_segments(segments)
end

"""
    resample_to_consistent(ds, Y, r, dt; kwargs...)

Refines the time series `Y` so consecutive points are adjacent boxes in
quantized space, with optional sphere-bundle consistency.
"""
function resample_to_consistent(ds, Y, r, dt::Real; pp=identity, sb_radius=nothing, sb_fct=nothing, verbose=false, max_depth=512)
    segments = Vector{Matrix{eltype(Y)}}(undef, size(Y, 2))
    @views for i in 1:(size(Y, 2)-1)
        segments[i] = resample_to_consistent(
            ds,
            Y[:, i],
            Y[:, i+1],
            r,
            dt;
            pp=pp,
            sb_radius=sb_radius,
            sb_fct=sb_fct,
            verbose=verbose,
            max_depth=max_depth,
        )
    end
    @views segments[end] = reshape(_to_vector(Y[:, end]), :, 1)
    return _concat_segments(segments)
end

function resample_to_consistent(ds, ic, ec, r, dt::Real; pp=identity, sb_radius=nothing, sb_fct=nothing, verbose=false, max_depth=512)
    direction_f = sb_fct
    direction_r = 1
    direction_mode = :boxes
    direction_norm = v -> norm(v, Inf)
    direction_pp = identity
    if !isnothing(sb_radius)
        direction_pp = v -> sb_radius * v
    else
        direction_f = nothing
    end
    is_close = _make_consistency_check(
        r=r,
        mode=:boxes,
        metric=euclidean,
        pp=pp,
        direction_f=direction_f,
        direction_r=direction_r,
        direction_mode=direction_mode,
        direction_metric=euclidean,
        direction_norm=direction_norm,
        direction_pp=direction_pp,
    )
    stepper = _ds_stepper(ds)
    return _resample_segment_stepper(
        _to_vector(ic),
        _to_vector(ec),
        dt,
        is_close,
        stepper;
        max_depth=max_depth,
        verbose=verbose,
    )
end

function is_dyn_consistent(p1, p2, r)
    p1q = round.(Int, p1 / r)
    p2q = round.(Int, p2 / r)
    return norm(p2q - p1q, Inf) <= 1
end

function sample_desol_to_distance_sb(sol, vf, sb_distance, r, dt::Real; max_depth=1000, t0=sol.t[1], t_max=sol.t[end])
    state_at = t -> begin
        x = _to_vector(sol(t))
        v = _to_vector(vf(x))
        return vcat(x, v)
    end

    if t_max <= t0
        X = reshape(state_at(t0), :, 1)
        return X, [1, 2]
    end

    n_intervals = ceil(Int, (t_max - t0) / dt)
    segments = Vector{Matrix{eltype(state_at(t0))}}(undef, n_intervals + 1)
    for i in 1:n_intervals
        t_start = t0 + (i - 1) * dt
        t_end = min(t_start + dt, t_max)
        segments[i] = _resample_segment_time(
            t_start,
            t_end,
            state_at,
            (a, b) -> sb_distance(a, b) <= r;
            max_depth=max_depth,
            verbose=false,
        )
    end
    segments[end] = reshape(state_at(t_max), :, 1)
    return _concat_segments(segments)
end

function sample_desol_to_distance_sb(sol, vf, dt::Real, t0, dist, r; max_depth=10000)
    t_end = t0 + dt
    state_at = t -> begin
        x = _to_vector(sol(t))
        v = _to_vector(vf(x))
        return vcat(x, v)
    end
    return _resample_segment_time(
        t0,
        t_end,
        state_at,
        (a, b) -> dist(a, b) <= r;
        max_depth=max_depth,
        verbose=false,
    )
end

function count_dynamic_inconsistencies(boxIt, d=1)
    box_it = Iterators.Stateful(boxIt)
    a = peek(box_it)
    if a === nothing
        return 0
    end
    a = first(box_it)
    return count(box_it) do b
        is_inconsistent = chebyshev(a, b) > d
        a = b
        return is_inconsistent
    end
end

function find_inconsistent_boxes(X, t, d=1)
    @views Xt = X[:, t]
    @views diff = Xt[:, 2:end] - Xt[:, 1:end-1]
    is_inconsistent = falses(size(diff, 2))
    for (j, col) in enumerate(eachcol(diff))
        is_inconsistent[j] = norm(col, Inf) > d
    end
    return unique([t[[is_inconsistent; false]] t[[false; is_inconsistent]]], dims=2)
end

function quantize_sb(VF_t, r)
    VF_t_rescaled = similar(VF_t)
    for (j, col) in enumerate(eachcol(VF_t))
        @views VF_t_rescaled[:, j] = r .* col ./ norm(col, Inf)
    end
    return quantize(VF_t_rescaled, 1)
end

function remove_successive_duplicates(t)
    t_new = [t[1]]
    for i = 2:length(t)
        if t[i] != t[i-1]
            push!(t_new, t[i])
        end
    end
    return t_new
end
