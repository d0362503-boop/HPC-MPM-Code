# `module/solid/explicit/` — Explicit MPM Solid Solver

`module/solid/explicit/` implements the explicit time-integration solid MPM solver as the `ExplicitSolidMPM` class.

## Overview

`ExplicitSolidMPM` inherits from `SolidMaterialPointBase` and provides the standard explicit MPM cycle:

1. **P2G** — project particle mass, volume, momentum, and internal/external forces to control points.
2. **Solve** — compute nodal acceleration from force/mass and apply boundary conditions.
3. **G2P** — update particle velocity, deformation gradient, volume, stress, and position.

The class supports the same interpolation family as the fluid MPM path (PIC, FLIP, TPIC, APIC) through the shared `MapScheme` switch and helpers in `map_and_interpolate.h`.

## File Inventory

| File | Purpose |
|------|---------|
| `explicit_mpm_solid.h` | `ExplicitSolidMPM` class declaration. |
| `var_trans_solid_explicit.cpp` | P2G (`Particle2Node`), acceleration solve (`SolveSolid`), and G2P (`Node2Particle`). |
| `fbar_projection.cpp` | F-bar volumetric locking correction (`ComputeDefGradBar`). |
| `CMakeLists.txt` | Module build rules; all `.cpp` files are compiled into `mpm_modules`. |

## Entry Points

`DataInput(grid, params)` reads the case. `Particle2Node(ts)`, `SolveSolid(ts)`,
and `Node2Particle(ts)` perform the cycle above. `DoMUSL(ts)` reprojects particle
momentum after the velocity update, before updating position and deformation.

## F-bar Projection

When `Fbar_flag` is enabled, `ComputeDefGradBar()` corrects the volumetric
part of the incremental deformation gradient to reduce locking.

The projection consists of two steps:

1. Project the uncorrected updated Jacobian `J_new = det(F_bar_old) * det(dF)` weighted by particle volume to the control points, then average via the nodal volume `nvof`.
2. Interpolate the corrected Jacobian `Jbar` back to each particle and scale the incremental deformation gradient `dF` by `cbrt(Jbar / J)` so that `det(F_bar_new) = Jbar`.

## Integration

The explicit force/mass update does not use `CrsMat` or PETSc KSP; implicit
equation-compression and preconditioner policies do not apply.

- `SolveSolid` intentionally does **not** update `nvel`; velocity is advanced inside `Node2Particle` (or preserved for FLIP) before MUSL is called.
- `UpdateDefGrad` and `UpdateVolume` live in `SolidMaterialPointBase` and are shared with the implicit solid path.
- `CutOffSmallNodalVar`, `ApplyVelocityBC`, and `ApplyAccelerationBC` are defined in `MaterialPoint` and reused by the fluid and implicit solid solvers.
