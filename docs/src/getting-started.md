```@meta
CurrentModule = CyclingSignatures
```

# Getting Started

We start with a time series from the Lorenz equations.


```@example lorenz
using OrdinaryDiffEq
using Plots, StatsPlots

const σ = 10.0
const ρ = 28.0
const β = 8/3

function lorenz!(du, u, p, t)
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 40.0))
sol  = solve(prob)

plot(sol, idxs=(1, 2, 3), xlabel="x", ylabel="y", zlabel="z", legend=false)
```


Although we usually recommend it, for this example we do not enforce any condition on the sampling density through the integration.

## Creating a TrajectorySpace

A `TrajectorySpace` represents a time series and a cubical cover `` Y``.
The cycling signature of a segment is the subspace of ``H_1(Y)`` induced by including a thickening of the segment into the ``Y``.

We generate time series data by evaluating the integrated trajectory on an equidistant grid.
We then compute the direction (i.e. normalized tangent vector) at each point.

```@example lorenz
using LinearAlgebra

dt = 0.01
tgrid = 0:dt:40

# time series vector
X = Array(sol(tgrid))

# compute unit tangent vectors
TX = mapslices(X, dims=1) do x
    v = zeros(3)
    lorenz!(v, x, 0, 0)
    return normalize(v)
end
nothing # hide
```

We can construct a `TrajectorySpace` simply by specifying the time series and two parameters for the size of the boxes used in the cubical cover. 
More precisely, the cover uses boxes of size `boxsize` in space and `1/sb_radius` in the unit tangent direction (more precisely, the unit tangent vectors are scaled to `sb_radius` and then covered with boxes with unit side length).

```@example lorenz

using CyclingSignatures

x_boxsize = 8.0
sb_radius = 1

traj_space = utb_trajectory_space_from_trajectory(X, TX, x_boxsize, sb_radius)
```



## Subsegment Experiments

We analyze cycling properties for a collection of randomly sampled segments of different lengths.
In the following, we select 100 segments of lengths `10:10:500` from the given dataset.

```@example lorenz

segment_lengths = collect(10:10:500)
n_runs = 100
exp = RandomSubsegmentExperiment(traj_space, segment_lengths, n_runs, 42)

results = run_experiment(exp)
```

The resulting collection of cycling signatures can be analyzed to obtain a coarse descrption of the dynamics on the Lorenz attractor.


```@example lorenz
evaluation_radius = 3.0

plt_rank = plot_rank_distribution_at_r(results, evaluation_radius)
sig1, plt1 = plot_subspace_frequency_at_r(results, 1, evaluation_radius; n_subspaces=3)
sig2, plt2 = plot_subspace_frequency_at_r(results, 2, evaluation_radius)
pltinc = plot_cycspace_inclusion(sig1, sig2)

combined = plot(plt_rank, plt1, plt2, pltinc; size=(800, 600), layout = (2, 2))
```

The top left shows the distribution of cycling ranks. While very short segment have rank 0 (which makes sense because they are too short to wrap around anything), after a certain time span all segments have nontrivial cycling rank. Therefore cycling is typical in this time series.

we see three rank 1 signatures.
one much more frequent than



here:
coarse description of the Lorenz system, we see the three oscillations, two wings and transitions,
comment that the time series is uneven and that is reflected in the cycling signatures,
longer dataset reveals more fine grained information
