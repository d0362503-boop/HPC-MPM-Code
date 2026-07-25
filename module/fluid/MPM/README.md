# `fluid/MPM/` — Stabilized MPM Fluid Solver

## Overview

This directory contains the **stabilized Material Point Method (MPM) fluid solver** for the incompressible Navier–Stokes equations.  Particles carry the primary fluid state and are mapped to a background grid via PIC / FLIP / TPIC / APIC.  The solver treats the fluid velocity as the time derivative of a particle displacement, so the NS system is solved for a displacement increment and pressure using a **Newmark-β / generalized-α time integrator** inside a **Newton–Raphson (NR) loop**.

## File Inventory

| File                       | Responsibility                                                                                              |
| -------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `stabilized_mpm.h`         | `StabilizedMPM` class declaration, inflow buffer `ifp`, NS system owner, stabilization selector             |
| `var_trans_fluid.cpp`      | `Particle2Node` (P2G) and `Node2Particle` (G2P), including APIC `Dmat`/`InvDmat` and particle-shifting call |
| `solve_ns_mpm.cpp`         | Newton–Raphson loop, VMS/PSPG stabilization coefficients, NS assembly/solve                                 |
| `fluid_point_inflow.cpp`   | Empty-mesh and filled-mesh inflow particle generation                                                       |
| `fluid_material_point.cpp` | `InitializePointData`, `MoveParticle` (migration + inflow), `MigrateParticleData`, `ApplyDLB`, `RebuildBoundaryConditions` |
| `fluid_mpm_data_io.cpp`    | Input/output of BC/point data and plain-text restart; captures global BC IDs when `do_dlb` is on                     |

Shared helpers live outside this directory:

- `module/particle_shifting.cpp` — `DeltaCorrectionParticleShifting()` and `PairwiseRepulsiveParticleShifting()`.
- `module/nonlinear_time_integration.cpp` — `PredictNewmarkBetaVelAndAccel`, `CommitNodalKinematics`, `CommitParticleKinematics`.
- `module/map_and_interpolate.cpp` — `VelP2G`, `AccelP2G`, `PICFamilyVelG2P`, `PICFamilyAccelG2P`.

## Class Architecture

```
MaterialPoint (base)
    └── StabilizedMPM
            ├── NS_            : CrsMat (ndof = 4, FEM_flag = false, owner_ = this)
            ├── ifp            : inflow particle buffer
            ├── tau1 / tau2    : per-particle stabilization coefficients
            ├── stab_coeff     : compile-time selector (VMS or PSPG)
            └── solswitch      : PIC / FLIP / TPIC / APIC (MapScheme enum)
```

`StabilizedMPM` differs from the FEM fluid path in a few key defaults:

| Setting                | MPM fluid default                        | FEM fluid default    |
| ---------------------- | ---------------------------------------- | -------------------- |
| `ode_order`            | `2` (displacement/velocity/acceleration) | `1` (velocity only)  |
| `NS_.FEM_flag`         | `false`                                  | `true`               |
| `NS_.use_petsc`        | `true` (PETSc Schur field-split)         | `true` (PETSc Schur field-split) |
| `NS_.use_schur_fieldsplit` | `true`                              | `true`               |
| Pressure Pmat          | PETSc `SELFP + epsilon I`                | PETSc `SELFP`        |
| `NS_.amg_rebuild_freq` | `1`                                      | `20`                 |

## Driver Flow

The standalone MPM fluid driver is `work/src_fluid/MPM/main.cpp`:

```
DataInput()
RestartInput()          (if rstflag == 1 or 3)
BuildMesh()
BuildControlPoint()
NS_.BuildCrsMat(16)
MakNodalVol()
for istep = ista .. iend
    MeshPointLinklist()
    Particle2Node()
    SolveNS()
    Node2Particle()
    if (do_dlb && istep % iout == 0)
        ApplyDLB()          (if nprocs > 1: MoveParticle + repartition + rebuilds)
    else
        MoveParticle()      (if nprocs > 1)
    output / restart    (if istep % iout == 0)
```

`ApplyDLB()` redistributes the background mesh across ranks from the particle
distribution and rebuilds the partition-dependent data (mesh, `nvol`, BCs, `NS_`).
See `module/DLB/README.md` for the full pipeline and its invariants.

## Time Integration

### 1. Newmark-β / Generalized-α predictor

At the start of each NR iteration the nodal displacement increment `ndispl` is used to predict velocity and acceleration:

```cpp
PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k);
```

Parameters (`alpha_f`, `alpha_m`, `gamma_nb`, `beta_nb`) are set from `spec_rad` by `GeneralizedAlphaParaSet` / `NewmarkBetaParaSet` in the base class.  Because `ode_order == 2`, the fluid MPM uses the second-order generalized-α formulas.

### 2. Newton–Raphson loop (`SolveNS`)

```cpp
for NR_it = 0 .. iter_max
    BCNRSet();                              // apply Dirichlet displacement/pressure values
    PredictNewmarkBetaVelAndAccel(...);     // predictor
    AssembleNSSystem(nvel_k, naccel_k);     // build NS_.amat / NS_.b_rhs
    iter = NS_.SolveSystem(NR_it);          // PETSc or native linear solve
    UpdateNRIncrement();                    // ndispl += x_lhs[0:3*nodec], npres += x_lhs[npc:]
    if CheckNRConvergence(...) break;
```

`CrsMat::SolveSystem` dispatches to native `GPBiCGAR` when `use_petsc == false` and to PETSc otherwise. PETSc divergence is reported by the KSP reason; it does not automatically switch solvers.

## Numerical Method

### P2G (`Particle2Node`)

For each particle inside each element:

- mass → `nmass`
- volume → `nvof`
- pressure (mass-weighted) → `npres_old`
- velocity (scheme-aware: PIC/FLIP standard, TPIC/APIC with gradients/B-matrix) → `nmome`
- acceleration (scheme-aware) → `nforce`

After overlap communication (`NodeVarComm`), nodal velocity and acceleration are obtained by `CutOffSmallNodalVar`, and velocity/acceleration BCs are applied.

### Stabilization coefficients (`MakNSStabCoeff`)

`tau1` and `tau2` are computed at each particle from the local velocity magnitude and element size.  The compile-time selector chooses between:

- `StabCoeff::VMS` (default) — VMS-style `τ`.
- `StabCoeff::PSPG` — PSPG-style `τ`.

### System assembly (`AssembleNSSystem`)

For each particle:

1. Evaluate shape functions and gradients at the particle position.
2. Apply `ImplicitDsfCorr` to correct gradients for the deformed configuration.
3. Interpolate velocity, acceleration, pressure, and pressure gradient to the α-level using `alpha_f` / `alpha_m`.
4. Compute the total stress tensor `stress_k` from the velocity gradient: the deviatoric part plus `-pres_k` on the diagonal.
5. Accumulate the 4×4 block matrix into `NS_.amat` and the RHS into `NS_.b_rhs`.

The RHS has a Galerkin part (`RHS_G`) and a stabilized part (`RHS_S`).  After the element loop, `NodeVarComm` synchronizes the RHS and `AddInertialForceToRHS` adds the inertial contribution.

### Linear solve (`SolveSystem`)

```cpp
if (NS_.use_petsc)  →  PETSc KSP
else                 →  GPBiCGAR (native iterative solver)
```

The MPM fluid default is PETSc with `use_schur_fieldsplit = true`. Its monolithic (u,p)
system is preconditioned by a lower Schur-complement field split: velocity and pressure
blocks each get an HYPRE/BoomerAMG hierarchy, with `SELFP` as the Schur preconditioner.
Its `A_00` inverse approximation is PETSc's 3x3 `(u,v,w)` `blockdiag`, rather than a
scalar diagonal.
Because its `FEM_flag` is false, the MPM pressure Pmat is the regularized `SELFP + epsilon I`, with
`epsilon = 1.0e-4 * abs(trace_active(SELFP)) / n_active`; inactive identity rows are excluded from
the scale. The code keeps PETSc's `PC_FIELDSPLIT_SCHUR_PRE_SELFP` selection and shifts its internally
built pressure Pmat in place after field-split setup; it does not switch to `SCHUR_PRE_USER`. The outer
FGMRES system remains unchanged.
The FEM fluid `NS_` system uses unshifted `SELFP` through `FEM_flag == true`; the MPM and FEM
rebuild frequencies differ.
Set `NS_.use_petsc = false` to fall back to the native solver.

### G2P and particle motion (`Node2Particle`)

1. Predict `nvel_k` / `naccel_k` and commit nodal kinematics (`CommitNodalKinematics`).
2. Reset particle velocity/acceleration and TPIC/APIC state.
3. Interpolate pressure, displacement, acceleration, and velocity back to particles (scheme-aware).
4. Recompute APIC `inv_Dmat` if needed.
5. Compute a particle-shifting correction (`PairwiseRepulsiveParticleShifting` is currently active; `DeltaCorrectionParticleShifting` is available but commented out in `Node2Particle`).
6. Commit particle kinematics (`CommitParticleKinematics`), which updates FLIP velocity, advects positions, and optionally applies the shifting correction while rejecting shifts that cross the domain boundary.

## Particle Shifting

Two particle-shifting strategies are implemented in `module/particle_shifting.cpp`:

- **`PairwiseRepulsiveParticleShifting()`** — currently active in `Node2Particle`.  Computes a pairwise repulsive correction between neighboring particles within a local support radius and returns a correction displacement.  It uses a uniform spatial-hash grid with cell size `dxy` and a 3×3×3 cell neighbor search to avoid the $O(N^2)$ all-pairs cost; per-particle `support = cbrt(vol[ip])` is still used for the actual distance filter and for deciding which boundary particles to send as MPI ghosts.
- **`DeltaCorrectionParticleShifting()`** — alternative.  Builds a volume-deficit field on control points and returns a correction displacement that pushes particles away from over-dense regions.  Currently commented out in `Node2Particle`.

## Inflow Particles

Inflow generation follows the same three-layer dispatch as before (`InflowMeshisFilled` → `GenerateInflowParticles` → `AssignUniqueInflowIds`).  The MPM-specific differences are:

- `InflowParticles()` allocates the inflow buffer including scheme-specific TPIC/APIC state.
- `GenerateInflowParticlesEmptyMesh` fills empty boundary cells with a Gaussian sub-grid.
- `GenerateInflowParticlesFilledMesh` clones particles that moved one cell inward; it copies the scheme-specific state (velocity/acceleration gradients or APIC B-matrices) by G2P interpolation.
- `MoveParticle` communicates the scheme-specific state for both regular migration and inflow appending.

> **Calling-order note:** `MeshPointLinklist()` is called **before** `SolveNS()`/`Node2Particle()`, while inflow runs inside `MoveParticle()` **after** particle advection.  The linked list is therefore intentionally one step behind; the filled-mesh generator relies on this.

## Data I/O

| Routine                          | Format                                               | Content                                                                           |
| -------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`)                    | `coord`, `id`, `matid`, `mass`, `vol`, `pres`, `vel`, `accel`, TPIC/APIC matrices |
| `OutputPointDataVTKHDF`          | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Particle coordinates, velocity, pressure, ID                                      |

## Integration with the Rest of the Codebase

- `module/solver/crsmat.h`: generic sparse-matrix wrapper; BC hooks injected through virtual overrides declared in `stabilized_mpm.h`.
- `module/map_and_interpolate.*`: PIC/FLIP/TPIC/APIC transfer logic.
- `module/material_point.h`: base class providing `PredictNewmarkBetaVelAndAccel`, `CommitNodalKinematics`, `CommitParticleKinematics`, `SolveSystem`, and the virtual BC/inflow hooks.
- `module/nonlinear_time_integration.cpp`: shared Newmark-β / generalized-α helpers.
- `module/particle_shifting.cpp`: shared particle-shifting implementations.
- `module/mpi_data.h`: particle migration and overlap-node communication.

## Coding Style Notes

- Class lives in namespace `stabilizedmpm`; `.cpp` files use `using namespace stabilizedmpm;`.
- Members use trailing underscore per `AGENTS.md` style (e.g. `NS_`), although the base `MaterialPoint` uses plain names for historical reasons.
