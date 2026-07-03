# `fluid/MPM/` — Stabilized MPM Fluid Solver

## Overview

This directory contains the **stabilized Material Point Method (MPM) fluid solver** for the incompressible Navier–Stokes equations. It shares the same `StabilizedMPM` class hierarchy and VMS/PSPG stabilization concepts as the FEM fluid path, but uses particles as the primary simulation entities and maps them to a background grid via PIC / FLIP / TPIC / APIC.

## File Inventory

| File | Responsibility |
|------|----------------|
| `stabilized_mpm.h` | `StabilizedMPM` class declaration, inflow particle buffer `ifp`, NS system owner |
| `var_trans_fluid.cpp` | `Particle2Node` (P2G) and `Node2Particle` (G2P), including APIC `Dmat`/`InvDmat` |
| `solve_ns_mpm.cpp` | Newton–Raphson loop for the stabilized NS system, stabilization coefficients |
| `fluid_point_inflow.cpp` | Empty-mesh and filled-mesh inflow particle generation |
| `fluid_material_point.cpp` | `InitializePointData`, `Moveparticle` (migration + appending inflow particles) |
| `fluid_mpm_data_io.cpp` | Input/output of BC/point data and plain-text restart |

## Class Architecture

```
MaterialPoint (base)
    └── StabilizedMPM
            ├── NS_       : CrsMat (ndof = 4, owner_ = this)
            ├── ifp       : inflow particle buffer
            ├── tau1/tau2 : per-particle stabilization coefficients
            └── solswitch : PIC / FLIP / TPIC / APIC
```

## Inflow Particles

Inflow generation is split into three layers. The dispatch and generation helpers are `private` overrides in `StabilizedMPM`; they are invoked polymorphically through the `MaterialPoint` base class.

1. **Detection** — `MaterialPoint::InflowMeshisFilled` checks whether any marked boundary cell has a fill ratio above `0.95`.
2. **Dispatch** — `MaterialPoint::GenerateInflowParticles` switches between:
   - `GenerateInflowParticlesEmptyMesh`: fills marked boundary cells layer-by-layer on a Gaussian sub-grid.
   - `GenerateInflowParticlesFilledMesh`: clones particles that moved one cell inward back into the boundary cell.
3. **ID assignment** — `MaterialPoint::AssignUniqueInflowIds` computes the global maximum particle ID and uses `MPI_Exscan` to assign contiguous unique IDs across ranks.

> **Important calling-order note:** `MeshPointLinklist()` is called **before** `Node2Particle()` in the driver, while inflow generation runs inside `Moveparticle()` **after** particle advection. Therefore `idepf`/`idp2p` reflect pre-move cell occupancy. `GenerateInflowParticlesFilledMesh` relies on this: it walks the boundary cell's pre-move linked list and selects particles whose current position satisfies `iexp == ie + sign`, i.e. particles that crossed into the neighboring cell during the step.

## Numerical Method

- **P2G** (`Particle2Node`): mass, volume, momentum, and force are mapped from particles to control points.
- **NS solve** (`SolveNS`): stabilized Navier–Stokes with VMS/PSPG coefficients computed per particle from local velocity and element size.
- **G2P** (`Node2Particle`): nodal velocity, acceleration, pressure, and displacement are interpolated back to particles; APIC/TPIC high-order matrices are recomputed.
- **Particle motion**: `CommitParticleKinematics` updates particle positions and velocities.

## Data I/O

| Routine | Format | Content |
|---------|--------|---------|
| `RestartOutput` / `RestartInput` | Plain text, per-rank (`*_re.txt`) | `coord`, `id`, `matid`, `mass`, `vol`, `pres`, `vel`, `accel`, TPIC/APIC matrices |
| `OutputPointDataVTKHDF` | VTK HDF5, single shared file per view (`*-w.vtkhdf`) | Particle coordinates, velocity, pressure, ID |

## Integration with the Rest of the Codebase

- `module/solver/crsmat.h`: generic sparse-matrix wrapper; fluid-specific BC hooks are injected through virtual overrides declared in `stabilized_mpm.h`.
- `module/map_and_interpolate.*`: PIC/FLIP/TPIC/APIC transfer logic.
- `module/material_point.h`: base class providing inflow dispatch hooks and shared particle data layout.
- `module/mpi_data.h`: particle migration and overlap-node communication.

## Coding Style Notes

- Class lives in namespace `stabilizedmpm`; `.cpp` files use `using namespace stabilizedmpm;`.
- Members use trailing underscore per `AGENTS.md` style (e.g. `NS_`), although the base `MaterialPoint` uses plain names for historical reasons.
