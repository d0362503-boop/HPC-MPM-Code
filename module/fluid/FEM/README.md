# `fluid/FEM/` — Stabilized FEM Fluid Solver

## Overview

This directory contains the **stabilized finite-element (FEM) fluid solver** for the incompressible Navier–Stokes equations, combined with a **volume-of-fluid (VOF) / phase-field** interface tracking method.

The solver uses a **SUPG/PSPG stabilized formulation** with **generalized-α time integration**.  Two `CrsMat` sparse-matrix objects are owned by `StabilizedFEM`:

| Matrix | DOF | Purpose |
|--------|-----|---------|
| `NS_`  | 4 (u, v, w, p) | Momentum + continuity system |
| `PF_`  | 1 (φ)          | Phase-field / volume-fraction (prepared, currently unused in the shown source) |

## File Inventory

| File | Responsibility |
|------|----------------|
| `stabilized_fem.h` | Class declaration and inline default constructor |
| `solve_ns_fem.cpp` | Time integration, stabilization-coefficient calculation, system assembly, and linear solve |
| `phase_field.cpp`  | Phase-field initialization, volume-fraction projection (`Particle2NodePhi`), material-property blending (`SetPFDomain`), liquid-volume calculation (`CalLiquidVol`) |
| `fluid_fem_data_io.cpp` | Restart I/O (`RestartInput` / `RestartOutput`), per-step field output (`OutputMeshDataVTKHDF`), VTK nodal interpolation (`Cp2NodeVTK`) |

Default boundary-condition methods are implemented inline in
`stabilized_fem.h`; derived owners may extend these virtual hooks.

## Boundary Conditions

`StabilizedFEM` inherits from `MaterialPoint`. `CrsMat` calls its owner
for velocity/pressure constraints through two default hooks in
`stabilized_fem.h`:

- `BuildPetscBCList(CrsMat&)`: collects global IDs from `ubc/vbc/wbc/pbc`
  for PETSc row/column elimination.
- `BCResidualSet(std::vector<double>&)`: zeros constrained residual entries.

Derived owners can extend these hooks.

## Numerical Method

`SolveNS` calls the phase helpers below: `ComputeAdvectionVel`,
`MakeNSStabCoeff`, and `AssembleSystem`.

### 1. Generalized-α Predictor
```cpp
adv_vel = (1 + α_f) · nvel_old  –  α_f · nvel_older
```
This advection velocity is frozen during the nonlinear/linear solve of the current time step.

### 2. SUPG/PSPG Stabilization Coefficients
Per-element values `tau1` (SUPG/PSPG) and `tau2` (shock capturing) are computed from:
- Element size `he`
- Local advection velocity magnitude
- Kinematic viscosity `rnue = rmue / rhoe`

### 3. System Assembly (`AssembleSystem`)
For each particle inside each element, shape functions are evaluated and the following matrices are accumulated into `NS_.amat`:

| Block | Physics |
|-------|---------|
| Mass | `emd` (Galerkin + SUPG) |
| Advection | `akg`, `akpu/v/w` (SUPG + PSPG) |
| Diffusion | `su/v/w`, cross terms (viscous stress) |
| Pressure | `esgu/v/w` (SUPG), `elt` (PSPG) |
| Continuity | `Cowu/v/w` (divergence) |
| Shock Capturing | `sc` (artificial diffusion) |

The RHS `NS_.b_rhs` is built from the old velocity, old acceleration, and body-force terms consistent with the generalized-α scheme.

### 4. Linear Solve (`SolveNSSystem`)
```cpp
if (NS_.use_petsc)  →  PETSc KSP
else                 →  GPBiCGAR (custom iterative solver)
```
After convergence, the solution vector `NS_.x_lhs` is unpacked back into `nvel` (u, v, w) and `npres` (p).

`ConfigurePreconditioner()` selects `PCFIELDSPLIT` directly for `NS_`, applying a lower
Schur field split to the four-component stabilized fluid system: components 0--2 form the
velocity block and component 3 forms the pressure block. Both child blocks use
HYPRE/BoomerAMG with PETSc's unshifted SELFP pressure Pmat and default scalar
diagonal inverse approximation and no additional pressure shift.
The scalar `PF_` system delegates PC selection to the base BoomerAMG callback;
the four-component NS system uses the velocity/pressure split. `PF_` remains
reserved rather than advanced by the current module workflow.

FEM retains the full equation layout; MPM equation compression does not apply.
The NS constructor sets `rebuild_freq=20`. Between scheduled rebuilds,
matrix values update while the numerical PC can be reused; FieldSplit groups
are preserved by checking the actual PC type. See the
[shared solver guide](../../solver/README.md).

### 5. Time-History Update (`UpdateNodalVar`)
At the end of each time step, the velocity / pressure history levels are moved down
one level to prepare for the next step:
```
nvel_older ← nvel_old
nvel_old   ← nvel
npres_old  ← npres
```

## Phase-Field / Two-Fluid Treatment

`Particle2NodePhi` projects a sharp particle indicator (liquid vs. gas) onto the mesh via standard FEM interpolation, then clips to `[0, 1]`.  `SetPFDomain` blends element-wise density `rhoe` and viscosity `rmue` from the liquid (`rhol`, `rmul`) and gas (`rhog`, `rmug`) values using the nodal volume fraction `nphi`.

## Data I/O

| Routine | Format | Content |
|---------|--------|---------|
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`) | `naccel`, `nvel`, `nvel_old`, `nvel_older`, `npres`, `nphi` |
| `OutputMeshDataVTKHDF` | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Global `Points`, hexahedron cells, `Velocity`, `Pressure`, `Phi` written collectively over `MPI_COMM_WORLD` |
| `Cp2NodeVTK` | Internal interpolation | Copies cell-corner data to nodes for VTK output, with clipping of `nphi_vtk` near 0 or 1 |
