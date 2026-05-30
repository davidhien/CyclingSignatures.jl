using Test
using LinearAlgebra
using Distances
using OrdinaryDiffEqTsit5

using CyclingSignatures

mutable struct TestStepper{F,T,P}
    f!::F
    u::Vector{T}
    p::P
    t::T
end

function set_state!(ds::TestStepper, u)
    ds.u .= u
    return ds
end

get_state(ds::TestStepper) = copy(ds.u)

function step!(ds::TestStepper, Δt, _)
    u = ds.u
    p = ds.p
    t = ds.t

    internal_dt = 1e-3
    n = max(1, ceil(Int, Δt / internal_dt))
    h = Δt / n

    du1 = similar(u)
    du2 = similar(u)
    du3 = similar(u)
    du4 = similar(u)
    tmp = similar(u)

    for _ in 1:n
        ds.f!(du1, u, p, t)
        @. tmp = u + (h / 2) * du1
        ds.f!(du2, tmp, p, t + h / 2)
        @. tmp = u + (h / 2) * du2
        ds.f!(du3, tmp, p, t + h / 2)
        @. tmp = u + h * du3
        ds.f!(du4, tmp, p, t + h)

        @. u = u + (h / 6) * (du1 + 2 * du2 + 2 * du3 + du4)
        t += h
    end
    ds.t = t
    return ds
end

function trajectory(ds::TestStepper, nsteps::Int, Δt)
    d = length(ds.u)
    X = Matrix{eltype(ds.u)}(undef, d, nsteps + 1)
    X[:, 1] = ds.u
    for i in 1:nsteps
        step!(ds, Δt, true)
        X[:, i + 1] = ds.u
    end
    return X
end

@testset "dynamical systems based refinement" begin

    function lorenz!(du, u, p, t)
        σ, ρ, β = p
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return du
    end

    function dadras!(du, u, p, t)
        a, b, c = p
        du[1] = a * u[1] - u[2] * u[3] + u[4]
        du[2] = u[1] * u[3] - b * u[2]
        du[3] = u[1] * u[2] - c * u[3] + u[1] * u[4]
        du[4] = -u[2]
        return du
    end

    function max_consecutive_distance(X, dist)
        return maximum(dist(X[:, i], X[:, i + 1]) for i in 1:(size(X, 2) - 1))
    end

    function max_consecutive_boxes(X, r)
        return maximum(
            norm(round.(Int, X[:, i] ./ r) .- round.(Int, X[:, i + 1] ./ r), Inf) for
            i in 1:(size(X, 2) - 1)
        )
    end

    for (f!, u0, p) in (
        (lorenz!, [1.0, 0.0, 0.0], (10.0, 28.0, 8 / 3)),
        (dadras!, [10.0, 1.0, 10.0, 1.0], (8.0, 40.0, 14.9)),
    )
        ds = TestStepper(f!, copy(u0), p, 0.0)
        X = trajectory(ds, 25, 0.05)

        r = max_consecutive_distance(X, euclidean) / 2
        Y, t_vec = CyclingSignatures.resample_to_distance(ds, 0.05, X, r; metric=euclidean, max_depth=256)
        @test t_vec[1] == 1
        @test t_vec[end] == size(Y, 2) + 1
        @test max_consecutive_distance(Y, euclidean) ≤ r + 1e-10

        r_box = 0.25
        Z, _ = CyclingSignatures.resample_to_distance(ds, 0.05, X, r_box; metric=euclidean, mode=:boxes, max_depth=256)
        @test max_consecutive_boxes(Z, r_box) ≤ 1

        W, _ = CyclingSignatures.resample_to_consistent(ds, X, r_box, 0.05; max_depth=256)
        @test all(
            CyclingSignatures.is_dyn_consistent(W[:, i], W[:, i + 1], r_box) for
            i in 1:(size(W, 2) - 1)
        )
    end
end

@testset "OrdinaryDiffEqTsit5 based refinement" begin
    function max_consecutive_distance(X, dist)
        return maximum(dist(X[:, i], X[:, i + 1]) for i in 1:(size(X, 2) - 1))
    end

    function rescale(x)
        return x ./ sqrt(norm(x))
    end

    function check_refinement(sol, vf, p; dt=0.1, t0=0.0, t_max=1.0)
        tgrid = t0:dt:t_max
        Xcoarse = reduce(hcat, (Vector(sol(t)) for t in tgrid))
        r = maximum(euclidean(Xcoarse[:, i], Xcoarse[:, i + 1]) for i in 1:(size(Xcoarse, 2) - 1)) / 2
        X, _ = CyclingSignatures.resample_to_distance(sol, dt, r; metric=euclidean, t0=t0, t_max=t_max)
        @test max_consecutive_distance(X, euclidean) ≤ r + 1e-10

        r_box = 1.0
        Y, _ = CyclingSignatures.resample_to_distance(sol, dt, r_box; metric=euclidean, t0=t0, t_max=t_max, mode=:boxes)
        @test all(
            CyclingSignatures.is_dyn_consistent(Y[:, i], Y[:, i + 1], r_box) for
            i in 1:(size(Y, 2) - 1)
        )

        Z, _ = CyclingSignatures.resample_to_distance(
            sol,
            dt,
            r,
            metric=euclidean,
            t0=t0,
            t_max=t_max,
            pp=rescale,
            direction_f=x -> vf(x, p),
            direction_r=1.5,
        )
        @test maximum(euclidean(rescale(Z[:, i]), rescale(Z[:, i + 1])) for i in 1:(size(Z, 2) - 1)) ≤ r + 1e-10
    end

    @testset "Lorenz" begin
        function lorenz!(du, u, p, t)
            σ, ρ, β = p
            du[1] = σ * (u[2] - u[1])
            du[2] = u[1] * (ρ - u[3]) - u[2]
            du[3] = u[1] * u[2] - β * u[3]
            return du
        end

        function lorenz_v(u, p)
            du = similar(u)
            lorenz!(du, u, p, 0.0)
            return du
        end

        u0 = [1.0, 0.0, 0.0]
        p = (10.0, 28.0, 8 / 3)
        prob = ODEProblem(lorenz!, u0, (0.0, 1.0), p)
        sol = solve(prob, Tsit5(); abstol=1e-7, reltol=1e-7)
        check_refinement(sol, lorenz_v, p)
    end

    @testset "Dadras" begin
        function dadras!(du, u, p, t)
            a, b, c = p
            du[1] = a * u[1] - u[2] * u[3] + u[4]
            du[2] = u[1] * u[3] - b * u[2]
            du[3] = u[1] * u[2] - c * u[3] + u[1] * u[4]
            du[4] = -u[2]
            return du
        end

        function dadras_v(u, p)
            du = similar(u)
            dadras!(du, u, p, 0.0)
            return du
        end

        u0 = [10.0, 1.0, 10.0, 1.0]
        p = (8.0, 40.0, 14.9)
        prob = ODEProblem(dadras!, u0, (0.0, 1.0), p)
        sol = solve(prob, Tsit5(); abstol=1e-7, reltol=1e-7)
        check_refinement(sol, dadras_v, p)
    end
end

@testset "spline interpolation based refinement" begin
    function max_consecutive_distance(X, dist)
        return maximum(dist(X[:, i], X[:, i + 1]) for i in 1:(size(X, 2) - 1))
    end

    t = range(0, 2π; length=9)
    Y = reduce(hcat, ([cos(τ), sin(τ)] for τ in t))
    r = 0.25

    Y_refined, derivatives, t_vec = CyclingSignatures.interpolate_to_distance(Y, r, euclidean)

    @test size(Y_refined, 1) == size(Y, 1)
    @test size(derivatives) == size(Y_refined)
    @test t_vec[1] == 1
    @test t_vec[end] == size(Y_refined, 2) + 1
    @test max_consecutive_distance(Y_refined, euclidean) ≤ r + 1e-10

    c, dc = CyclingSignatures.getInterpolationFunction(Y)
    @test c(1) ≈ CyclingSignatures.get_interpolation_function(Y)[1](1)
    @test dc(1) ≈ CyclingSignatures.get_interpolation_function(Y)[2](1)

    Y_compat, derivatives_compat, t_vec_compat =
        CyclingSignatures.interpolateToDistance(Y, r, euclidean)
    @test Y_compat ≈ Y_refined
    @test derivatives_compat ≈ derivatives
    @test t_vec_compat == t_vec
end
