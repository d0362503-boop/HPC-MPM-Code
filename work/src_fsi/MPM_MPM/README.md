# Monolithic MPM-MPM FSI

`MPMMPMMonolithicFSI` assembles fluid, solid, and interface multipliers into one
PETSc system. The ten components are fluid displacement (3), fluid pressure (1),
solid displacement (3), and interface multipliers (3), stored in 37 scalar blocks.

## Current solver strategy

| Layer | Implementation |
|-------|----------------|
| Native matrix interface | One `CrsMat::ndof` member; full `nodec * ndof` vectors and original component numbering |
| MPM linear solve | `CrsMat` removes exact zero-increment identity equations, solves the compact system and restores the full owned/ghost solution |
| Newton correction | Zero initial guess; physical RHS preserved before scaling and restored after solving |
| Root preconditioner | Lower Schur split: fluid / [solid + lambda] |
| Coupled preconditioning | PCKSP inner FGMRES acts on PETSc's assembled `selfp` Pmat, retaining solid/lambda coupling |
| Fluid block | FGMRES with CGS refinement always, rtol 0.5, max 100, restart 100; velocity/pressure lower Schur; velocity and pressure PREONLY + ASM/local ILU(1), RCM ordering, nonzero shift 1e-3 |
| Preconditioner refresh | Constructor sets `amg_rebuild_freq=1`; rebuild at NR=0 each step and reuse for NR>0 when the retained layout is unchanged |
| Solid block | FGMRES with CGS refinement always, `restart=30`, `rtol=0.1`, `max_it=100`; ASM overlap 1, local PREONLY/ILU(2), RCM ordering, nonzero shift `1e-3` |
| Mass activation | Each phase uses `1e-4` times its mean effective nodal mass; wet-interface and physical constraints determine active multiplier components |

```text
Outer FGMRES: fluid / [solid + lambda], lower Schur split
  fluid: FGMRES with CGS refinement always, velocity / pressure Schur split
    velocity: PREONLY + ASM / ILU(1), RCM ordering
    pressure: PREONLY + ASM overlap 1 / local PREONLY + ILU(1), RCM ordering
  [solid + lambda]: PREONLY + PCKSP
    inner FGMRES acts on the assembled selfp Pmat
    solid / lambda: lower Schur split
      solid: FGMRES + ASM overlap 1 / local PREONLY + ILU(2), RCM ordering
      lambda: PREONLY + ASM overlap 1 / local PREONLY + ILU(1)
```

The root `selfp` matrix uses the fluid diagonal inverse to approximate
elimination. The inner split retains the solid operator and its multiplier
coupling inside that approximation. The outer FGMRES solves the full coupled
equations. PETSc constructs the Schur preconditioning matrices; no custom
`BuildInterfaceSchur` implementation or global direct solver is used.

Velocity and pressure use ASM/local ILU(1) with RCM ordering and
nonzero shift 1e-3. Pressure uses PREONLY, ASM overlap 1 and local PREONLY/ILU(1).
Solid uses FGMRES with rtol 0.1 and max 100, ASM overlap 1 and local
PREONLY/ILU(2), retaining RCM ordering and nonzero shift 1e-3. The multiplier Schur PC retains
ASM/local ILU(1). These choices affect the preconditioner only.

The compiled fluid FGMRES default now includes
-fsi_fieldsplit_fluid_ksp_gmres_cgs_refinement_type refine_always.
A second classical Gram-Schmidt pass prevents loss of Krylov-basis
orthogonality on difficult pressure directions. The assembled VMS equations,
field splits, tolerances and iteration limits are unchanged.
Velocity RCM alone did not resolve the step 2395 fluid failure; refinement is
the additional change.

## Step 1720 checkpoint comparison (2026-10-20)

Both ten-step tests use the same validated 60-partition TPIC checkpoint,
pressure and velocity ASM/local ILU(1), unchanged physical parameters and
tolerances, and sh run.sh 60.

| Solid PC (FGMRES rtol 0.1, max 100) | Steps | Solver time | Outer total (30 corrections) | Max per correction | Negative reasons |
|---|---|---:|---:|---:|---:|
| ASM overlap 1 / ILU(1) | 1721-1730 | 253.6049 s | 1895 | 72 | 0 |
| ASM overlap 1 / ILU(2) | 1721-1730 | 244.5542 s | 1898 | 71 | 0 |

ILU(2) reduced solver time by about 3.6%, with essentially unchanged outer
iteration counts. Single-step screening also tried ASM overlap 2, HYPRE
Euclid levels 1/2, and HYPRE ILU with GMRES-ILUk. None improved total time
over ASM overlap 1 / ILU(2). The formal solid default is ASM/local ILU(2); velocity and pressure retain ILU(1).
No rigid-body near-nullspace code is retained. This restart comparison does
not establish convergence over the preceding or remaining simulation.

The current history starts from initial particle data with ista=1 and
rstflag=2. Steps 1-120 used solid ILU(1), with no negative KSP reasons.
A complete step 120 TPIC checkpoint was verified against both phase outputs
at t=0.12 s, including all 120 restart files. The continuation uses ista=121,
iend=5000, rstflag=3 and iout=10, launched by sh run.sh 60.
The 60-partition initial input contains 4140 solid and 626175 fluid particles.
Current defaults are velocity/pressure ASM/local ILU(1), solid FGMRES +
ASM/local ILU(2), and NR0 rebuild / NR>0 reuse. Physics, tolerances and all
other PC settings remain unchanged.

## Step 1260 checkpoint comparison (2026-10-20)

All rows start from the same validated 60-partition TPIC checkpoint and use
`sh run.sh 60`. Solver time includes initialization and output. Outer totals
in the two completed rows cover 20 linear corrections (two per time step).

| Velocity PC | Pressure PC | Solid PC | Refresh policy | Result / solver time | Outer total | Negative reasons |
|---|---|---|---|---|---:|---:|
| ASM/ILU(1) | BoomerAMG | ASM/ILU(1) | Every Newton correction | 1261–1270 passed / 249.1756 s | 1489 | 0 |
| BoomerAMG | BoomerAMG | ASM/ILU(1) | Every Newton correction | 1261–1270 passed / 341.9262 s | 1476 | 0 |
| BoomerAMG | BoomerAMG | BoomerAMG | NR0 rebuild, NR>0 reuse | Failed at 1261, NR0; stopped | No completed solve | 207 solid inner DIVERGED_ITS |
| ASM/ILU(1) | BoomerAMG | ASM/ILU(1) | NR0 rebuild, NR>0 reuse | 1261–1270 passed / 234.1054 s | 1490 | 0 |

The solid-AMG row follows the requested within-step reuse policy. Its failure
occurs at NR0 with a fresh preconditioner, before reuse can occur. Every failing
solid inner solve reached the unchanged limit of 100 iterations. This ordinary
scalar BoomerAMG configuration has not passed the ten-step test. The solid-AMG
trial was stopped; that comparison restored velocity and solid ASM/ILU(1),
retaining pressure BoomerAMG and within-step reuse. The
first two timings use per-correction rebuilding. The latest pressure-only AMG
row confirms that NR0-only rebuilding also passes the same ten steps without
negative KSP reasons. Physical parameters and convergence tolerances were held
fixed throughout.

## Launch configuration and logging

The 1261–1270 restart test uses the first Slosh block in
build/data/generate/fsi/MPM_MPM/input.txt:

    <sp.solswitch, rstflag, nlstep>
    TPIC  3  1
    <ista, iend, iout>
    1261  1270  10
    <dt, mtol, wp.spec_rad, sp.spec_rad>
    1.0e-3  1.0e-11  0.0e0  0.0e0

rstflag=3 reads the complete step 1260 TPIC checkpoint and continues writing
restart files. A fresh run uses rstflag=2 and ista=1 instead. This checkpoint
has 60 partitions; the MPI count must match. Verify both phase HDF time,
particle IDs/coordinates and all 120 restart files before restarting.

Launch with the compiled defaults and failure-only diagnostic viewers:

    cd /home/ps/HPC-MPM-Code/build
    export PETSC_OPTIONS="-fsi_fieldsplit_fluid_ksp_converged_reason ::failed -fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_converged_reason ::failed -fsi_fieldsplit_coupled_ksp_ksp_converged_reason ::failed"
    sh run.sh 60

The viewers do not override numerical settings or print successful inner KSP
solves. Newton records remain Monolithic_NR and Monolithic_converge; negative
outer KSP and Newton failures remain visible. The driver also prints startup,
final summaries and OutputMessage every 100 steps.

Current Slosh gravity uses real_time=dt*istep and
s=real_time^2/(2*0.51) before 0.51 s, then s=real_time-0.255.
The angle is theta=(pi/45)*sin(2*pi*s/1.211), with
gx=-9.81*sin(theta) and gz=-9.81*cos(theta).
This user-selected gravity law and VMS stabilization are unchanged by the
orthogonalization update.

## Restart verification

The current fluid block uses FGMRES, CGS refinement always, restart 100,
rtol 0.5 and max 100. Velocity and pressure use ASM/local ILU(1) with RCM ordering
and nonzero shift 1e-3. Pressure uses PREONLY and ASM overlap 1.
Solid retains FGMRES rtol 0.1 / max 100 with ASM overlap 1, local ILU(2),
RCM ordering and nonzero shift 1e-3.
The FSI constructor sets amg_rebuild_freq=1. The first Newton correction in
each step rebuilds the compressed KSP and its nested preconditioners; subsequent
corrections reuse them. A changed retained layout still requires rebuilding.
Iteration-growth forced rebuilding is disabled for this FSI configuration. Physical
parameters, VMS equations and nonlinear acceptance criteria are unchanged.
The tables below describe earlier configurations; current settings are listed above.

The following full-FSI tests start from the same validated step 2580 TPIC
checkpoint and use sh run.sh 61 with no numerical PETSC_OPTIONS overrides:

| Pressure local ILU ordering | Completed test steps | Negative KSP reasons | Ten-step time |
|---|---|---:|---|
| Natural, fluid restart 100 | 2581–2587; stopped on the first failure in 2588 | 1 | Not completed |
| RCM, fluid restart 100 | 2581–2590, normal finish | 0 | 271.2 s from launch to complete two-phase output |

The accepted RCM test's output-to-output interval is 267.8 s. Both phase
HDF files are complete at time 2.59 s with 4140 solid and 627210 fluid
particles. All 122 TPIC restart files match their particle positions and
velocities exactly; ownership validation reports no misplaced particles.
This passes the ten-step no-negative-reason and at-most-300-second targets.

The same compiled solver then restarts from that verified step 2590 state
with ista 2591, iend 5000 and rstflag 3. The earlier pre-2580 history contains
fluid-inner failures, so the new branch must not be described as a clean
5000-step run from the initial state. The ten-step test is not proof that
all remaining steps will satisfy the same limits.

## Source responsibilities

| File | Responsibility |
|------|----------------|
| `main.cpp` | Time-step order, transfers, migration, output, and restart |
| `datain_para.cpp` | Coupled case parameters and phase input |
| `monolithic_fsi.cpp` | Coupled lifecycle, Newton loop, increments, and convergence |
| `monolithic_assembly.cpp` | Fluid, solid, and multiplier contributions; 37-block PETSc insertion |
| `monolithic_interface.cpp` | Solid contour, wet-interface weights, and particle correction |
| `monolithic_solver.cpp` | Active components, constraints, scaling, linear correction, and physical preconditioner configuration |
| `monolithic_fsi.h` | Class declarations and function documentation |
| `../../../module/interface_sdf.h` / `.cpp` | Indexed finite-triangle distance queries |

The class header exposes driver operations and the physics callbacks used by
`CrsMat`. Assembly, scaling, convergence and geometry helpers are private.
`fsi_sys.ndof` is the single coupled dimension (10); `kStoredBlocks` defines its 37 stored blocks;
the constructor configures it in `monolithic_fsi.cpp`. `InterfaceSDF` accepts
oriented triangles and bin spacing; wet-interface activation remains in FSI.

CMake collects this directory's `.cpp` files automatically. After adding or
removing one, run `cmake -S . -B build` from the repository root, then build the
existing `build/` directory.

## Linear-system boundary

The FSI coordinator assembles with the original control points and ten components.
Native RHS and correction vectors keep their full `nodec * ndof` length and
component-major layout. It delegates equation compression and full solution
recovery to `CrsMat`; no compact equation numbering belongs in the assembly.

`CrsMat` removes only exact unit-diagonal equations with a zero RHS after
boundary elimination. It derives each field's compact index set from original
global equation IDs through `SetFieldSplitComponents()`. FSI declares only
PC choices and physical groups in `ConfigurePreconditioner()`:

```text
fluid {0,1,2,3}
  0 {0,1,2}
  1 {3}
coupled {4,5,6,7,8,9}
  solid {4,5,6}
  lambda {7,8,9}
```

These groups are stored in `CrsMat::field_splits`. The ordinary
`CrsMat::ConfigureReducedPreconditioner()` member owns compact root setup,
child-group numbering, nested PETSc setup and PCKSP operator selection.
It also handles pure-fluid flat groups and pure-solid ASM without groups.
FSI has no separate reduced-setup method or function callback.

The outer FGMRES splits fluid components (0..3) from solid and multiplier
components (4..9). PETSc `selfp` eliminates the fluid diagonal when assembling
the preconditioning matrix for the coupled solid-multiplier block. A PCKSP
application runs FGMRES on that assembled approximation, retaining the full solid
operator and its multiplier coupling. It does not repeatedly solve the fluid
block inside this inner iteration.

The inner coupled preconditioner splits solid (4..6) and multiplier (7..9)
components with another PETSc `selfp` Schur split. Solid remains FGMRES with
rtol 0.1 and max 100, ASM overlap 1, and local PREONLY/ILU(2) with nonzero shift
1e-3. The new coupled iteration uses rtol 0.5 and max 100. These are
preconditioning solves; the outer and physical Newton tolerances control the
final correction. `BuildInterfaceSchur` remains removed.

`SolveSystem(NR_it)` preserves the unscaled physical RHS, zeros the native
correction, scales the matrix/RHS and solves through `CrsMat`. It then multiplies
the recovered correction by `column_scale` and restores the physical RHS.
The outer maximum is 500. When any rank has positive lumped interface weight,
outer `rtol=1e-6` and `atol=1e-12`; otherwise they are `1e-10` and `1e-15`.
The divergence tolerance is `1e6` on the first Newton correction and
`PETSC_UNLIMITED` later. A logged negative reason remains a failed linear solve;
current code does not automatically abort or retry it.

Current option prefixes are `fsi_fieldsplit_fluid_` for fluid,
`fsi_fieldsplit_coupled_ksp_` for the coupled inner iteration, and
`fsi_fieldsplit_coupled_ksp_fieldsplit_solid_` /
`fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_` for its phase splits.
Old `fsi_fieldsplit_fields_` option names no longer configure these groups.

Current launch settings use the source defaults:

| Setting | Value |
|---------|-------|
| Mass activity | Relative phase mean times 1e-4 |
| Outer restart | 200 |
| Fluid tolerance / maximum / restart | 0.5 / 100 / 100 |
| Fluid orthogonalization | Classical Gram-Schmidt with refinement always |
| Velocity / pressure local ILU ordering | RCM / RCM |
| Solid tolerance / maximum / restart | 0.1 / 100 / 30 |
| Solid orthogonalization | Classical Gram-Schmidt with refinement always |
| Coupled tolerance / maximum / restart | 0.5 / 100 / 100 |
| Coupled orthogonalization | Classical Gram-Schmidt with refinement always |

The inner tolerances control preconditioning work. Outer linear convergence
and the four physical Newton criteria below determine correction acceptance.

The shared compression path is already used by pure solid MPM, pure fluid MPM
and this FSI system. `CrsMat::SolveReducedSystem()` creates the compact solve and
recovers the full correction. `SetFieldSplitComponents()` supplies field indices
from the same compact layout. Pure fluid defines its velocity/pressure split
with that API; pure solid applies ASM/ILU directly to the retained equations.
Only the physical group declarations differ between these systems;
`CrsMat` owns the reusable setup and numbering machinery.

The physical matrix, material parameters and nonlinear tolerances are separate
from this preconditioning layer. The full native CSR graph and full PETSc
assembly still exist; reducing their inactive storage and assembly communication
is a later change inside `CrsMat`.

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
| Fluid momentum | Nonlinear momentum RHS [N] | `max(1e-10, 1e-4 * initial_norm[0])` |
| Continuity | Nonlinear continuity RHS [m^3/s] | `max(1e-12, 1e-4 * initial_norm[1])` |
| Solid momentum | Nonlinear momentum RHS [N] | `max(1e-10, 1e-4 * initial_norm[2])` |
| Interface | Displacement-increment difference [m] | `1e-10` |

The first three references are saved inside `CheckNRConvergence()` at
`NR_it == 0`. All four RMS comparisons must meet their thresholds.
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

`BuildSolidInterface(solid_phi)` interpolates the supplied control-point phase
field to the existing mesh vertices, using the scalar interpolation from
`Cp2NodeVTK()` without its display cutoffs. Each `nc` hexahedron is split into
six tetrahedra around its `0--6` body diagonal; no extra subcell vertices are
generated. The extracted `phi = 0.5` triangles are oriented outward and exchanged
across MPI ranks. Evaluating the interpolation as `0.5 + sum(N * (phi - 0.5))`
preserves the constant contour level at partition faces without snapping values.
Mesh vertices are used for geometry only; the coupled fields
and interface-area weights remain on control points. `LumpedLagrangeMultiplier()` uses
the step-start `solid_.nphi` and nearby volume-equivalent fluid boxes to integrate
wet-interface area into `nlm_lump_local`; `NodeVarComm` produces `nlm_lump`.
Local weights assemble matrix contributions, and synchronized weights enter
the RHS and activation decisions.

By default, each phase uses 1e-4 times its mean effective nodal mass as its
activity threshold. `BuildActiveDOFs()` accumulates four statistics:
fluid mass sum/count in entries 0/1 and solid mass sum/count in entries 2/3.
Only nodal masses strictly above `mtol` contribute, weighted by `dbc` to count
shared nodes once; `MPI_Allreduce` combines the phase statistics.

The cutoffs are fluid_cut=1e-4*mass_stats[0]/mass_stats[1]
and solid_cut=1e-4*mass_stats[2]/mass_stats[3].
mtol affects the averaging population; it is not itself the relative cutoff
or a lower bound on it. The current FSI code uses this relative rule directly;
there is no launch-time absolute-mass override.
A multiplier component requires both phases active and `nlm_lump[n] > 0`. Phase boundary
conditions and solid rigid constraints are then applied; a multiplier component
is fixed when both corresponding phase displacement components are fixed.

The coupled matrix assembly uses the per-component `fixed_dof` mask; it does
not call the generic `CrsMat::BuildActiveRowMask()`.

Changing `mtol` changes the statistics population used by the relative cutoff. It does not create interface
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
