# Monolithic MPM-MPM FSI

`MPMMPMMonolithicFSI` assembles fluid, solid, and interface multipliers into one
PETSc system. The ten components are fluid displacement (3), fluid pressure (1),
solid displacement (3), and interface multipliers (3), stored in 37 scalar blocks.

## Source responsibilities

| File | Responsibility |
|------|----------------|
| `main.cpp` | Time-step order, transfers, migration, output, and restart |
| `datain_para.cpp` | Coupled case parameters and phase input |
| `monolithic_fsi.cpp` | Newton loop, fluid/solid assembly wrappers, increments, and convergence |
| `monolithic_interface.cpp` | Solid contour, wet-interface weights, multiplier assembly, and correction |
| `monolithic_solver.cpp` | Active components, boundary constraints, scaling, and PETSc field split |
| `monolithic_fsi.h` | Class declarations and function documentation |
| `../../../module/interface_sdf.h` / `.cpp` | Indexed finite-triangle distance queries |

The former `fluid_mat.cpp` and `solid_mat.cpp` assembly wrappers live in
`monolithic_fsi.cpp`. `InterfaceSDF` accepts oriented triangles and bin spacing;
fluid particle domains and wet-interface logic belong to the FSI coordinator.
Its header declares the interface, and its `.cpp` implements the geometry.

CMake collects this directory's `.cpp` files automatically. After adding or
removing one, run `cmake -S . -B build` from the repository root, then build the
existing `build/` directory.

## Newton iteration

`SolveFSISystem()` builds the lumped interface weights and active-component mask
once at the beginning of the time step. Each Newton iteration then performs:

1. Clear the coupled matrix and RHS.
2. Apply phase boundary values and compute trial nodal velocity/acceleration.
3. Assemble fluid, solid, and interface contributions once.
4. Check the current nonlinear RHS and interface displacement difference.
5. If unconverged, solve a zero-initialized correction and update the unknowns.

The next iteration evaluates the updated state through the same assembly path.
`ComputeNRResidualNorms()` reads the unscaled, overlap-synchronized `b_rhs` for
fluid momentum, continuity, and solid momentum. The fourth field is the direct
`solid_.ndispl - fluid_.ndispl` difference, without lumped interface-area weights.

Each field uses an RMS over globally owned, unconstrained active components.
For the interface, the denominator is the active multiplier-component count,
not all `nodec` nodes or all mesh components. Shared components are counted
once. A field with no active components has zero RMS.

| Field | Residual | Acceptance threshold |
|-------|----------|----------------------|
| Fluid momentum | Nonlinear momentum RHS [N] | `max(1e-8, 1e-4 * initial_norm[0])` |
| Continuity | Nonlinear continuity RHS [m^3/s] | `max(1e-10, 1e-4 * initial_norm[1])` |
| Solid momentum | Nonlinear momentum RHS [N] | `max(1e-8, 1e-4 * initial_norm[2])` |
| Interface | Displacement-increment difference [m] | `1e-10` |

The first three references are saved inside `CheckNRConvergence()` at
`NR_it == 0`. All four RMS values must be finite and meet their thresholds.
The loop allows 100 Newton corrections; a failed check at `NR_it == 100`
aborts before another solve. `Monolithic_NR` and `Monolithic_converge` log the
Newton iteration, the preceding Krylov iteration count, then these four RMS
values in table order.

These checks assess the equations on the active components. Particle overlap
and the geometric correction require their own output checks.

The solid trial stress is currently initialized once before the loop. The
incremental linear-elastic constitutive model adds stress increments in place,
so repeated assembly can accumulate trial stress across iterations. This remains
a limitation of that material path; hyperelastic models evaluate stress from
the trial deformation gradient.

## Interface activation

`BuildSolidInterface(solid_phi)` reconstructs the `phi = 0.5` surface by sampling
the supplied nodal phase field, subdividing cells into tetrahedra, and exchanging
outward-oriented triangles across MPI ranks. `LumpedLagrangeMultiplier()` uses
the step-start `solid_.nphi` and nearby volume-equivalent fluid boxes to integrate
wet-interface area into `nlm_lump_local`; `NodeVarComm` produces `nlm_lump`.
Local weights assemble matrix contributions, and synchronized weights enter
the RHS and activation decisions.

Fluid components require `fluid_.nmass[n] > mtol`; solid components require
`solid_.nmass[n] > mtol`. A multiplier component needs both phases active and
`nlm_lump[n] > 0`. Phase boundary conditions and solid rigid constraints are
then applied. A multiplier component is fixed when both corresponding phase
displacement components are fixed.

The coupled matrix assembly uses the per-component `fixed_dof` mask; it does
not call the generic `CrsMat::BuildActiveRowMask()`.

Reducing `mtol` changes the mass activation test. It does not create interface
weights at control points where `nlm_lump` is zero. The assembled multiplier
equations enforce nodal displacement-increment continuity; their convergence
alone does not establish particle nonpenetration.

## Particle penetration correction

The driver calls the following sequence after the coupled solve:

```cpp
fsi.solid_.Node2Particle();
fsi.fluid_.Node2Particle();
fsi.solid_.MoveParticle();
fsi.CorrectFluidPenetration();
fsi.fluid_.MoveParticle();
```

`Node2Particle()` commits coordinates and kinematics, including the existing
particle shifting. `MoveParticle()` determines particle ownership and migrates
data; the fluid version also handles inflow. Correcting between the two migration
calls makes the final solid particles available on their owning ranks and lets
fluid ownership be determined from corrected coordinates.

`CorrectFluidPenetration()` maps current solid particle volume with
`LocateLocalElement`, `MakeSF`, and `StandardVarP2G` into a temporary nodal field.
After `NodeVarComm`, it divides by `nvol` and clamps the volume fraction to
`[0, 1]`. It passes this field to the existing `BuildSolidInterface()` routine;
the step-start contour cannot represent the final solid configuration.

For each fluid centre, `InterfaceSDF::FindSurfacePoint()` queries the closest
finite triangle within a band of `max(dxy)`. A negative signed distance causes
one direct coordinate assignment to that surface point. There is no iterative
retreat or fluid particle-to-grid mapping in this query. Triangle-interior
projection follows the normal; edge and vertex cases use the nearest finite
surface point.

The correction changes fluid coordinates only, preserving velocity,
acceleration, TPIC/APIC state, pressure, mass, and volume. It is a geometric
remapping separate from PST and does not supply a contact reaction. The distance
sign follows the nearest triangle rather than a global inside/outside test,
and the reconstructed `phi = 0.5` surface differs from a particle material-domain
boundary. Validation should identify which boundary and overlap measure it uses.
