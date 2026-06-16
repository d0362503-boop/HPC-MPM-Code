# `fluid/FEM/` — Stabilized FEM Fluid Solver

## Overview

This directory contains the **stabilized finite-element (FEM) fluid solver** for the incompressible Navier–Stokes equations, combined with a **volume-of-fluid (VOF) / phase-field** interface tracking method.

The solver uses a **SUPG/PSPG stabilized formulation** with **generalized-α time integration**.  Two `CrsMat` sparse-matrix objects are owned by `StabilizedFEM`:

| Matrix | DOF | Purpose |
|--------|-----|---------|
| `NS_`  | 4 (u, v, w, p) | Momentum + continuity system |
| `PF_`  | 1 (φ)          | Phase-field / volume-fraction (prepared, currently unused in the shown source) |

---

## File Inventory

| File | Responsibility |
|------|----------------|
| `stabilized_fem.h` | Class declaration, inline constructor, and global singleton `wfem` |
| `solve_ns_fem.cpp` | Time integration, stabilization-coefficient calculation, system assembly, and linear solve |
| `phase_field.cpp`  | Phase-field initialization, volume-fraction projection (`Particle2NodePhi`), material-property blending (`SetPFDomain`), liquid-volume calculation (`CalLiquidVol`) |
| `fem_data_io.cpp`  | Restart I/O (`RestartInput` / `RestartOutput`), per-step field output (`OutputMeshDataVTKHDF`), VTK nodal interpolation (`Cp2NodeVTK`) |

> **Note:** Fluid-specific boundary-condition overrides (`BuildPetscBCList`, `BCResidualSet`) are **declared** in `stabilized_fem.h` but **implemented** in `src_fsi/bc_setting.cpp` alongside the solid BC logic.

---

## Class Architecture

```
MaterialPoint  (base)
    └── StabilizedFEM
            ├── NS_  : CrsMat  (ndof = 4, owner_ = this)
            └── PF_  : CrsMat  (ndof = 1, owner_ = this)
```

`StabilizedFEM` **declares** two virtual hooks that the generic `CrsMat` calls back via the `owner_` pointer.  Their **implementations** live in `src_fsi/bc_setting.cpp` (shared with the solid BC logic):

### `BuildPetscBCList(CrsMat& mat)`
Collects all Dirichlet/essential boundary-condition global IDs for the PETSc `MatZeroRows` pass:
- `ubc`, `vbc`, `wbc`, `pbc` — standard velocity/pressure BCs
- `fsi_intf` — FSI interface constraints (applied to u, v, w)

### `BCResidualSet(std::vector<double>& rr)`
Zeroes out constrained degrees of freedom in a residual / RHS vector after a matrix-vector product:
- `fsi_intf` on u, v, w
- `wbc` on u, v, w
- `pbc` on p

These overrides replace the former free-function `BCResidualSet` that contained a hard-coded `if (ndof == 4) … else …` branch, achieving full polymorphic decoupling between `CrsMat` and fluid-specific logic.

---

## Numerical Method

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

### 3. System Assembly (`AssembleNSSystem`)
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

### 5. Time-Advance (`UpdateNodeVar`)
Shifts the velocity / pressure history arrays for the next time step:
```
nvel_older ← nvel_old
nvel_old   ← nvel
npres_old  ← npres
```

---

## Phase-Field / Two-Fluid Treatment

`Particle2NodePhi` projects a sharp particle indicator (liquid vs. gas) onto the mesh via standard FEM interpolation, then clips to `[0, 1]`.  `SetPFDomain` blends element-wise density `rhoe` and viscosity `rmue` from the liquid (`rhol`, `rmul`) and gas (`rhog`, `rmug`) values using the nodal volume fraction `nphi`.

---

## Data I/O

| Routine | Format | Content |
|---------|--------|---------|
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`) | `naccel`, `nvel`, `nvel_old`, `nvel_older`, `npres`, `nphi` |
| `OutputMeshDataVTKHDF` | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Global `Points`, hexahedron cells, `Velocity`, `Pressure`, `Phi` written collectively over `MPI_COMM_WORLD` |
| `Cp2NodeVTK` | Internal interpolation | Copies cell-corner data to nodes for VTK output, with clipping of `nphi_vtk` near 0 or 1 |

---

## Integration with the Rest of the Codebase

- **`CrsMat`** (`module/solver/`) — generic sparse-matrix wrapper; fluid-specific BC logic is injected through the virtual overrides declared here but implemented in `src_fsi/bc_setting.cpp`, so `crsmat.cpp` does **not** include `stabilized_fem.h`.
- **`MaterialPoint`** (`module/material_point.h`) — base class providing the `BuildPetscBCList` and `BCResidualSet` virtual hooks.
- **Solvers** (`module/solver/solver.cpp`) — generic GPBiCG* routines call `mat.MatVecMult(xx)` and `mat.owner_->BCResidualSet(rr)` without knowing whether the matrix belongs to fluid or solid physics.

---

## Coding Style Notes

- Trailing underscore for class members (`owner_`, `NS_`, `PF_`).
- Global singleton `wfem` is declared inline in `stabilized_fem.h`.
- `using namespace stabilizedfem;` appears in `.cpp` files only.
