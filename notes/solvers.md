# Sweep Cell Solvers

The SweepCellSolve class template hierarchy provides implementations of methods used for solving the coupled systems in each mesh element during a transport sweep. The components of this hierarchy are outlined below.

### SweepCellSolveBase

- Abstract class template.
- Derives from:
	- DomainDecomposition
	- Quadrule::OrdinateSet
- Defines virtual interface for solving systems:
	- `CalcA`
	- `CalcB`
	- `SolveSystem`
	- `StoreResult`
	- `GetDimOfA`
	- `GetDimOfB`
	- `GetDimOfWork`
	- `SetDt`
- Owns member variables:
	- dt
	- system_dim
- Constructor initializes:
	- dt
	
### SweepCellSolve

- Abstract class template.
- Derives from:
	- SweepCellSolveBase
- Specialized for each discretization type (RKDG, STDG, etc.).
- Owns member variables:
	- RKDG:
		- DG_degree
	- STDG:
		- DG_degree_x
		- DG_degree_t
- Specializations implement:
	- Indexing functions.
	- `StoreResult` (override).
- Constructors initialize:
	- **All**:
		- system_dim
	- RKDG:
		- DG_degree
	- STDG:
		- DG_degree_x
		- DG_degree_t
	
### SweepLinearSolve

- Abstract class template.
- Derives from:
	- SweepCellSolve
- Defines interface for:
	- `ComputeMatrices`
- Implements:
	- `GetDimOfA`
	- `GetDimOfB`
	- `GetDimOfWork`
- Specialized member functions implement:
	- `ComputeMatrices` to precompute angular component of matrices.
	- `CalcB` (override).
	- `CalcA` (override).
- Owns member variables:
	- Aq
	- sizeof_A
	- sizeof_B
	- dimof_A
	- dimof_B
- Constructors initialize:
	- Aq
	- sizeof_A
	- sizeof_B
	- dimof_A
	- dimof_B
	
### Sweep<Type\>Solve

- Concrete class template.
- Derives from:
	- SweepLinearSolve
- Implements:
	- `SolveSystem` (override).
	- `Print` (override).
- Specialized member functions implement:
	- `CalcA` (override) if needed (e.g., Hessenberg).
	- `ComputeMatrices` (override) to augment precomputation done by SweepLinearSolve initialization (e.g., Hessenberg reduction).
- Owns member variables:
	- tau (Hess).
	
	
# Sweep Pattern

The SweepPattern class template hierarchy provides implementations of methods used for traversing spatial meshes during a transport sweep. The components of this hierarchy are outlined below.

### SweepPattern

- Abstract class template.
- Defines virtual interface for performing sweeps:
	- `Sweep`
- Contains private helper routines:
	- `SetSigns`
- Owns member variables:
	- sweep_order
- Constructor initializes:
	- sweep_order
	
### SweepPatternUniCart

Encapsulates logic for performing transport sweeps on uniform Cartesian meshes.

- Abstract class template.
- Contains private helper routines:
	- ComputeDiagLength
- Prints:
	- Sweep order.
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	