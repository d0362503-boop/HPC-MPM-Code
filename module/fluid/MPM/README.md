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
| `fluid_material_point.cpp` | `InitializePointData`, `MoveParticle` (migration + inflow), `MigrateParticleData`, `ApplyDLB`, `RebuildBC`  |
| `fluid_mpm_data_io.cpp`    | Input/output of BC/point data and plain-text restart; captures global BC IDs when `do_dlb` is on            |

Shared helpers live outside this directory:

- `module/particle_shifting.cpp` — `DeltaCorrectionPST()` and `PairwiseRepulsivePST()`.
- `module/implicit/generalized_alpha_integrator.cpp` — Newmark/generalized-α coefficients and displacement-to-velocity/acceleration conversion.
- `module/material_point.cpp` — closed-boundary position helpers.
- `module/implicit/implicit_mpm.cpp` — `ImplicitMPM` (owned as `implicit_`): endpoint kinematics, inertia terms, `CommitNodalKinematics` and `CommitParticleKinematics`.
- `module/map_and_interpolate.cpp` — `VelP2G`, `AccelP2G`, `PICFamilyVelG2P`, `PICFamilyAccelG2P`.

## Class and Defaults

`StabilizedMPM` inherits from `MaterialPoint` and owns the four-component
`NS_` system, inflow buffer `ifp`, and per-particle `tau1/tau2`.
The compile-time `stab_coeff` selects VMS or PSPG; `solswitch` selects
PIC, FLIP, TPIC, or APIC transfers.

| Setting | MPM fluid | FEM fluid |
|---------|-----------|-----------|
| `ode_order` | 2 (displacement, velocity, acceleration) | 1 (velocity) |
| Equation layout | Full system | Full system |
| `NS_.rebuild_freq` | 1 | 20 |

## Driver Flow

The module's time-step operations are called in this order:

```
DataInput()
RestartInput()          (if rstflag == 1 or 3)
BuildMesh()
BuildControlPoint()
NS_.BuildCrsMat(16)
ComputeNodalVol()
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
See the [DLB guide](../../DLB/README.md) for the pipeline and its invariants.

## Time Integration

### 1. Newmark-β / Generalized-α predictor

Time-integration parameters belong to `implicit_.integrator_`
(`GeneralizedAlphaIntegrator`). `implicit_` is the `ImplicitMPM` component
(`module/implicit/implicit_mpm.h`) owned by the implicit MPM solvers. It holds only implicit
state; its kinematics, inertia, particle-commit and FLIP-update functions take the
particle object (`*this`) whose arrays they read or write. Explicit solvers do not carry it.

At the start of each NR iteration the nodal displacement increment `ndispl` is used to predict velocity and acceleration:

```cpp
implicit_.ComputeNodeVelAccelFromDispl(*this, nvel_k, naccel_k);
```

Parameters (`alpha_f`, `alpha_m`, `gamma_nb`, `beta_nb`) are set from `spec_rad` by `implicit_.integrator_.GeneralizedAlphaParaSet()` / `NewmarkBetaParaSet(dt)` in `DataInput`.  Because `ode_order == 2`, the fluid MPM uses the second-order generalized-α formulas.

### 2. Newton–Raphson loop (`SolveNS`)

`StabilizedMPM` is a `NonlinearProblem` (`module/implicit/newton_solver.h`);
`SolveNS()` calls `SolveNewton(*this, max_NR_it)` (`max_NR_it = 100`), which runs:

```cpp
BeginNewtonStep();                          // freeze stabilization, zero ndispl/npres
for NR_it = 0 .. max_NR_it
    AssembleNewtonSystem();                 // BCNRSet, predictor, zero amat/b_rhs/x_lhs, AssembleSystem
    if CheckNewtonConvergence(NR_it, ksp_its) break; // NS_.CheckNRConvergence on the physical RHS
    ksp_its = SolveNewtonCorrection(NR_it); // NS_.SolveSystem: PETSc or native linear solve
    UpdateNRIncrement();                    // ndispl += x_lhs[0:3*nodec], npres += x_lhs[npc:]
EndNewtonStep();                            // release NS_ storage
```

`CrsMat::SolveSystem` dispatches to native `GPBiCGAR` when `use_petsc == false` and to PETSc otherwise. PETSc divergence is reported by the KSP reason; it does not automatically switch solvers.

The NR monitor checks the assembled `b_rhs` before solving. After an increment,
the next pass reassembles at the updated state; each pass assembles only once.
Iteration zero saves the reference norm. For normalization and threshold caveats, see the
[solver limitations](../../solver/README.md#current-convergence-limitations).

## Numerical Method

### P2G (`Particle2Node`)

For each particle inside each element:

- mass → `nmass`
- volume → `nvof`
- pressure (mass-weighted) → `npres_old`
- velocity (scheme-aware: PIC/FLIP standard, TPIC/APIC with gradients/B-matrix) → `nmome`
- acceleration (scheme-aware) → `nforce`

After overlap communication (`NodeVarComm`), nodal velocity and acceleration are obtained by `CutOffSmallNodalVar`, and velocity/acceleration BCs are applied.

### Stabilization coefficients (`MakeNSStabCoeff`)

`tau1` and `tau2` are computed at each particle from the local velocity magnitude and element size.  The compile-time selector chooses between:

- `StabCoeff::VMS` (default) — VMS-style `τ`.
- `StabCoeff::PSPG` — PSPG-style `τ`.

### System assembly (`AssembleSystem`)

For each particle:

1. Evaluate shape functions and gradients at the particle's step-start position.
2. Interpolate velocity, acceleration, pressure, and pressure gradient to the α-level using `alpha_f` / `alpha_m`.
3. Compute the total stress tensor `stress_k` from the velocity gradient: the deviatoric part plus `-pres_k` on the diagonal.
4. Accumulate the 4×4 block matrix into `NS_.amat` and the RHS into `NS_.b_rhs`.

This fluid assembly uses the gradients from `MakeSF` directly; it does not call
the solid's `ImplicitDsfCorr` during NR assembly.

The RHS has a Galerkin part (`RHS_G`) and a stabilized part (`RHS_S`).  After the element loop, `NodeVarComm` synchronizes the RHS and `implicit_.AddInertialForce` adds the inertial contribution.

### Linear solve (`SolveSystem`)

The MPM fluid default is PETSc/FGMRES with a lower velocity/pressure Schur split.
Both child PCs use PREONLY + HYPRE/BoomerAMG; PETSc SELFP supplies the pressure Pmat.
PETSc uses its default scalar diagonal approximation for `A_00`, with no
additional pressure shift. FEM also uses unshifted SELFP and its own PC cadence.

The constructor sets `NS_.ndof=4`; native assembly and vectors remain
`nodec * 4`. `ConfigurePreconditioner()` declares velocity `{0,1,2}` and
pressure `{3}`; `CrsMat` owns field index sets, KSP/PC refresh and solution
recovery. See the [shared solver guide](../../solver/README.md).

### G2P and particle motion (`Node2Particle`)

1. Predict `nvel_k` / `naccel_k` and commit nodal kinematics (`implicit_.CommitNodalKinematics`).
2. Reset particle velocity/acceleration and TPIC/APIC state.
3. Interpolate pressure, displacement, acceleration, and velocity back to particles (scheme-aware).
4. Recompute APIC `inv_Dmat` if needed.
5. Compute a particle-shifting correction (`PairwiseRepulsivePST` is currently active; `DeltaCorrectionPST` is available but commented out in `Node2Particle`).
6. Commit particle kinematics (`implicit_.CommitParticleKinematics`), which updates FLIP velocity and adds displacement plus the supplied shifting correction to positions. The correction array must contain one entry per particle; the current commit routine does not reject shifting-induced domain crossings.

The shared commit does not apply a closed-face coordinate clamp. Applying nodal velocity/displacement BCs does not by itself
guarantee that a finite PIC displacement stays inside the domain. Particle
migration removes globally out-of-domain points; open-face outflow must be
distinguished from closed-face crossings when reviewing results.

## Particle Shifting

Two particle-shifting strategies are implemented in `module/particle_shifting.cpp`:

- **`PairwiseRepulsivePST()`** — currently active in `Node2Particle`.  Computes a pairwise repulsive correction between neighboring particles within a local support radius and returns a correction displacement.  It uses a uniform spatial-hash grid with cell size `dxy` and a 3×3×3 cell neighbor search to avoid the $O(N^2)$ all-pairs cost; per-particle `support = cbrt(vol[ip])` is still used for the actual distance filter and for deciding which boundary particles to send as MPI ghosts.
- **`DeltaCorrectionPST()`** — alternative.  Builds a volume-deficit field on control points and returns a correction displacement that pushes particles away from over-dense regions.  Currently commented out in `Node2Particle`.

## Inflow Particles

Inflow generation follows a three-layer dispatch (`InflowMeshisFilled` → `GenerateInflowParticles` → `AssignUniqueInflowIds`).  The MPM-specific differences are:

- `InflowParticles()` allocates the inflow buffer including scheme-specific TPIC/APIC state.
- `GenerateInflowParticlesEmptyMesh` fills empty boundary cells with a uniform sub-cell grid (see `SubcellOffsets`).
- `GenerateInflowParticlesFilledMesh` clones particles that moved one cell inward; it copies the scheme-specific state (velocity/acceleration gradients or APIC B-matrices) by G2P interpolation.
- `MoveParticle` communicates the scheme-specific state for both regular migration and inflow appending.

> **Calling-order note:** `MeshPointLinklist()` is called **before** `SolveNS()`/`Node2Particle()`, while inflow runs inside `MoveParticle()` **after** particle advection.  The linked list is therefore intentionally one step behind; the filled-mesh generator relies on this.

## Data I/O

| Routine                          | Format                                               | Content                                                                           |
| -------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`)                    | `coord`, `id`, `matid`, `mass`, `vol`, `pres`, `vel`, `accel`, TPIC/APIC matrices |
| `OutputPointDataVTKHDF`          | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Particle coordinates, velocity, pressure, ID                                      |
