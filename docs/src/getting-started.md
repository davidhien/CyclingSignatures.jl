```@meta
CurrentModule = CyclingSignatures
```

# Getting Started

We start with a time series from the Lorenz equations REF


```julia
using DifferentialEquations

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
```


Although we usually recommend it, for this example we do not enforce any condition on the sampling density through the integration.
TODO: Link!

## Creating a TrajectorySpace


## Subsegment Experiments


## Analysis

