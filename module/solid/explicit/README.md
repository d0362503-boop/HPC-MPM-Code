# Explicit Solid MPM

`module/solid/explicit/` implements the explicit time-integration solid MPM solver as the `ExplicitSolidMPM` class.

## Overview

`ExplicitSolidMPM` inherits from `SolidMaterialPointBase` and provides the standard explicit MPM cycle:

1. **P2G** — project particle mass, volume, momentum, and internal/external forces to control points.
2. **Solve** — compute nodal acceleration from force/mass and apply boundary conditions.
3. **G2P** — update particle velocity, deformation gradient, volume, stress, and position.

The class supports the same interpolation family as the fluid MPM path (PIC, FLIP, TPIC, APIC) through the shared `MapScheme` switch and helpers in `map_and_interpolate.h`.

## Files

| File | Purpose |
|------|---------|
| `explicit_mpm_solid.h` | `ExplicitSolidMPM` class declaration. |
| `var_trans_solid_explicit.cpp` | P2G (`Particle2Node`), acceleration solve (`SolveSolid`), MUSL (`DoMUSL`), and G2P (`Node2Particle`). |
| `fbar_projection.cpp` | F-bar volumetric locking correction (`ComputeDefGradBar`). |
| `CMakeLists.txt` | Module build rules; all `.cpp` files are compiled into `mpm_modules`. |

## Class responsibilities

- `DataInput()` — read standalone explicit-solid input parameters and particle data.
- `Particle2Node()` — reset nodal arrays, loop over particles to deposit mass/volume/momentum/force, then synchronize and apply velocity BCs.
- `SolveSolid()` — compute `naccel = nforce / nmass` with small-mass cutoff and acceleration BCs.
- `Node2Particle()` — update particle kinematics from the grid, run MUSL, update the deformation gradient (with optional F-bar correction), and update particle position/volume/stress.
- `DoMUSL()` — re-project particle momentum to the grid after the G2P velocity update so the subsequent position/deformation-gradient update uses a consistent grid velocity.
- `ComputeDefGradBar()` — nodal F-bar projection that corrects the volumetric part of the incremental deformation gradient to avoid locking.

## F-bar projection

When `Fbar_flag` is enabled, the solver applies the nodal F-bar projection described in:

> "Circumventing volumetric locking in explicit material point methods:  
> A simple, efficient, and general approach"

The projection consists of two steps:

1. Project the uncorrected updated Jacobian `J_new = det(F_bar_old) * det(dF)` weighted by particle volume to the control points, then average via the nodal volume `nvof`.
2. Interpolate the corrected Jacobian `Jbar` back to each particle and scale the incremental deformation gradient `dF` by `cbrt(Jbar / J)` so that `det(F_bar_new) = Jbar`.

## Integration

This directory is always built as part of `mpm_modules` through `module/solid/explicit/CMakeLists.txt`. It is consumed by:

- `work/src_fsi/block_fsi.cpp` for the FSI coupling path.
- `work/src_solid/explicit/main.cpp` for the standalone explicit-solid solver (currently still uses the old global interface and needs to be wired to the new class).

## Notes

- `SolveSolid` intentionally does **not** update `nvel`; velocity is advanced inside `Node2Particle` (or preserved for FLIP) before MUSL is called.
- `UpdateDefGrad` and `UpdateVolume` live in `SolidMaterialPointBase` and are shared with the implicit solid path.
- `CutOffSmallNodalVar`, `ApplyVelocityBC`, and `ApplyAccelerationBC` are defined in `MaterialPoint` and reused by the fluid and implicit solid solvers.
