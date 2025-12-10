```@meta
CurrentModule = CyclingSignatures
```

# CyclingSignatures.jl

Documentation for [CyclingSignatures.jl](https://github.com/davidhien/CyclingSignatures.jl).

## About

Cycling signatures are a tool to identify and classify oscillatory dynamics in (multivariate) time series. 

Given a time seris ``x_0,x_1,\dots, x_m \in \mathbb{R}^n``, the basic pipeline is:
1. Generate a comparison space ``Y`` by covering the time series and compute its homology ``H_1(Y)``.
2. Compute the homology ``H_1(\gamma)`` of thickened segments ``\gamma = \{x_i, x_{i+1},\dots, x_j\}``.
3. Compute the image of the map ``H_1(\gamma) \rightarrow H_1(Y)``, this is the *Cycling Signature*.


Collections of cycling signatures can provide a coarse description of recurrent dynamics.
For instance, frequently appearing one dimensional cycling signatures indicate elementary cycling motions.
These can be related using higher dimensional cycling signatures.

## Installation

Install the package via

```julia
using Pkg
Pkg.add(url="https://github.com/davidhien/CyclingSignatures.jl")

```
and load it via
```julia
using CyclingSignatures
```

## References

U. Bauer, D. Hien, O. Junge, and K. Mischaikow. Cycling signatures: Identifying cycling motions in time series using algebraic topology. Journal of Computational Dynamics 12(4), 554–595 (2025). DOI: 10.3934/jcd.2025007.