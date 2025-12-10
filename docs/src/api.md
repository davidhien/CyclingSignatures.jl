# Comparison Spaces

The idea of a comparison space is to have a space ``Y`` and a map ``\Gamma \rightarrow Y`` such that the maps

``
    H_1( O_r(\gamma) \rightarrow Y_r)
``

yield a meaningful classifier for a cycle in the thickning of ``\gamma``.
Currently, only static spaces ``Y`` are supported.

## Comparison Space Interface

The comparison space interface is defined by:

```
abstract type AbstractComparisonSpace end
```

The following methods are required:

- `betti_1(cs)`: Returns the first Betti number of the comparison space.
- `map_cycle(cs, points, simplices, coeffs)`: Returns a vector of length `betti_1(cs)` representing the image of the cycle mapped to the comparison space.

## Cubical Comparison Space

Cubical comparison spaces are comparison spaces which are cubical complexes. 
We implemented cubical comparison spaces in space alone and in the unit tangent bundle.
Furthermore, we provide an intereface which can be extended, e.g. to support different domains or periodic boundaries.

### Interface

The interface for cubical comparison spaces is defined by:

```
abstract type AbstractCubicalComparisonSpace <: AbstractComparisonSpace end
```

It assumes that a cubical acyclic carrier (i.e. a subtype of `AbstractCubicalAcyclicCarrier`) is used to map the cycle.

The following methods are required:
- `edge_boxes(cs, p1, p2)`: Compute the sequence of boxes covering the edge between points `p1` and `p2`.
- `carrier(cs)`: Returns the cubical acyclic carrier.
- `betti_1(comparison_space)`: Returns the first Betti number of the comparison space.

The implementation of `map_cycles` for `AbstractCubicalComparisonSpace`s does the following:

1. For each edge:
    1. Compute the boxes which it intersects (using `edge_boxes`).
    2. Use the carrier from `carrier(cs)` to get a 1-chain in the comparison space.
2. Sum and return the 1-chains.

### CubicalComparisonSpace and SBCubicalComparisonSpace

We provide two implementations of cubical comparison spaces:
- `CubicalComparisonSpace`: The cubes are assumed to cover the data in space only.
- `SBCubicalComparisonSpace`: The cubes are assumed to cover the data in the unit tangent bundle


## Acyclic Carrier Interface

The logic of representing cubical complexes is separated from the logic of the inclusion map.
Currently, only a Vietoris--Rips type implementation is supported.

# Sampleable Trajectory

*Under Construction.*
