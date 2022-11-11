# Public Interface

## Inlined

- Index: Computes linearized index.

## Virtual 

### DensityFunction (All)

- PointerAtCell [**<code>i,j</code>**]: Returns a pointer to coefficients in a spatial cell. Requires spatial indices only.
- CellStride [**<code>dim</code>**]: Returns the stride in the global array between subsequent spatial cells, moving in the specified dimension.
- DOFPerCell [**<code>i,j</code>**]: Returns the number of degrees of freedom per mesh element. This includes spatial, temporal, and angular degrees of freedom (when applicable).

### OrdinateFlux

- PointerAt[Quadrant/Octant]Cell [**<code>i,j,[quad/oct]</code>**]: Returns a pointer to coefficients in a spatial cell for the given quadrant/octant. Requires spatial and quadrant/octant indices.
- PointerAtOrdinateCell [**<code>i,j,q</code>**],: Returns a pointer to coefficients in a spatial cell for the given angular ordinate. Requires spatial and ordinate indices.