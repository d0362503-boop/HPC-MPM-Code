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
- `BCNRSet()` (virtual): writes Dirichlet values into the n+1 iterate
  `nvel_k` / `npres_k` before each assembly. Accelerations get no separate
  BC; they follow from the generalized-α velocity relation.

Derived owners can extend these hooks (e.g. `FSIFluid::BCNRSet` adds the
solid interface velocity).

## Numerical Method

`SolveNS` runs `SolveNewton()` over the `NonlinearProblem` hooks:

```
BeginNewtonStep        nvel_k ← nvel, npres_k ← npres   (same-velocity predictor)
loop:
  AssembleNewtonSystem BCNRSet → ComputeAdvectionVel → ComputeNodeAccelFromVel
                       → MakeNSStabCoeff → AssembleSystem
  CheckNewtonConvergence
  SolveNewtonCorrection
  UpdateNRIncrement    nvel_k += Δv, npres_k += Δp
EndNewtonStep
```

`SolveNS` does not commit; the owner calls `CommitNodalKinematics` once the
time step is accepted.

### 1. Generalized-α Time Levels
`nvel` / `npres` / `naccel` hold time level n; `nvel_k` / `npres_k` / `naccel_k` hold the
n+1 Newton iterate. The advection velocity is evaluated at n+α_f:
```cpp
adv_vel  = (1 - α_f) · nvel   +  α_f · nvel_k     // v^{n+α_f}
p_af     = (1 - α_f) · npres  +  α_f · npres_k    // p^{n+α_f}
accel_am = (1 - α_m) · naccel +  α_m · naccel_k   // a^{n+α_m}
```
Only ∂v/∂t is evaluated at n+α_m; convection, viscosity, pressure, and the
stabilization terms are evaluated at n+α_f (pressure follows Scheme-2 of
Liu et al. 2020).

### 2. SUPG/PSPG Stabilization Coefficients
Per-element values `tau_supg` (SUPG/PSPG, Tezduyar type) and `tau_lsic` (LSIC) are
recomputed from v^{n+α_f} at every nonlinear iteration, using:
- Element size `he`
- Local advection velocity magnitude
- Kinematic viscosity `rnue = rmue / rhoe`

Their derivatives are not included in the tangent (as in Bazilevs et al. 2008).

### 3. System Assembly (`AssembleSystem`)
For each particle inside each element, shape functions are evaluated and the following matrices are accumulated into `NS_.amat`:

| Block | Physics |
|-------|---------|
| Mass | `emd` (Galerkin + SUPG) |
| Advection | `akg`, `akpu/v/w` (SUPG + PSPG) |
| Diffusion | `su/v/w`, cross terms (viscous stress) |
| Pressure | `esgu/v/w` (SUPG), `elt` (PSPG) |
| Continuity | `Cowu/v/w` (divergence) |
| LSIC | `sc` (grad-div, `tau_lsic`) |

Tangent coefficients (unknowns Δv^{n+1}, Δp^{n+1}):

| Term | Coefficient |
|------|-------------|
| Mass `emd`, PSPG mass `egt*` | α_m / (γΔt) |
| Convection, diffusion, LSIC, continuity, PSPG convection | α_f |
| Pressure `esg*`, `elt` | α_f |

The convection tangent is Picard: only (v·∇)δv is kept; (δv·∇)v and the
velocity derivatives of the SUPG/PSPG weights are omitted.

The RHS `NS_.b_rhs` is the negative residual, evaluated with `accel_am`,
`adv_vel`, `p_af`, and the body force.

### 4. Linear Solve (`SolveNewtonCorrection`)
```cpp
if (NS_.use_petsc)  →  PETSc KSP
else                 →  GPBiCGAR (custom iterative solver)
```
`NS_.x_lhs` is the Newton correction; `UpdateNRIncrement` adds it to `nvel_k`
(u, v, w) and `npres_k` (p).

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

### 5. Time-History Update (`CommitNodalKinematics`)
Called once per accepted time step (after FSI block convergence, not inside `SolveNS`):
```
naccel ← (nvel_k − nvel)/(γΔt) − (1−γ)/γ · naccel
nvel   ← nvel_k
npres  ← npres_k
```

## Phase-Field / Two-Fluid Treatment

`Particle2NodePhi` projects a sharp particle indicator (liquid vs. gas) onto the mesh via standard FEM interpolation, then clips to `[0, 1]`.  `SetPFDomain` blends element-wise density `rhoe` and viscosity `rmue` from the liquid (`rhol`, `rmul`) and gas (`rhog`, `rmug`) values using the nodal volume fraction `nphi`.

## Data I/O

| Routine | Format | Content |
|---------|--------|---------|
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`) | `naccel`, `nvel`, `npres`, `nphi` |
| `OutputMeshDataVTKHDF` | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Global `Points`, hexahedron cells, `Velocity`, `Pressure`, `Phi` written collectively over `MPI_COMM_WORLD` |
| `Cp2NodeVTK` | Internal interpolation | Copies cell-corner data to nodes for VTK output, with clipping of `nphi_vtk` near 0 or 1 |
