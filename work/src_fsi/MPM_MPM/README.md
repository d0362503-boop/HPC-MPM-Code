# `work/src_fsi/MPM_MPM/` — Monolithic MPM-MPM FSI

`MPMMPMMonolithicFSI` assembles fluid, solid, and interface multipliers into one
PETSc system. The ten components are fluid displacement (3), fluid pressure (1),
solid displacement (3), and interface multipliers (3), stored in 37 scalar blocks.
The multipliers enforce continuity of the phase endpoint velocities. After the
solve, a separate geometric correction moves fluid centres toward the final
solid volume-fraction interface.

## Solver Strategy

The outer FGMRES solves the full coupled equations. Its preconditioner uses
two lower Schur splits, with a coupled inner iteration on PETSc's assembled
`selfp` preconditioning matrix:

```text
fluid / [solid + lambda]
  fluid: FGMRES, velocity / pressure Schur
    velocity, pressure: PREONLY + ASM / local ILU(1)
  [solid + lambda]: PREONLY + PCKSP
    inner FGMRES on selfp; solid / lambda Schur
      solid: FGMRES + ASM / local ILU(1)
      lambda: PREONLY + ASM / local ILU(1)
```

Fluid and solid FGMRES use classical Gram-Schmidt with refinement always.
Their local ILU factors use RCM ordering and nonzero shift `1e-3`.
The inner tolerances below control preconditioning work; physical acceptance
is determined separately by the four Newton residuals.

## File Inventory

| File | Responsibility |
|------|----------------|
| `main.cpp` | Time-step order, transfers, migration, output, and restart |
| `datain_para.cpp` | Coupled case parameters and phase input |
| `monolithic_fsi.cpp` | Coupled lifecycle, Newton loop, increments, and convergence |
| `monolithic_assembly.cpp` | Fluid, solid, and multiplier contributions; 37-block PETSc insertion |
| `monolithic_interface.cpp` | Lumped interface weights, final solid phase mapping, and fluid coordinate correction |
| `monolithic_solver.cpp` | Active components, constraints, scaling, linear correction, and physical preconditioner configuration |
| `monolithic_fsi.h` | Class declarations and function documentation |
| `../../../module/interface_sdf.h` / `.cpp` | Local solid phase interpolation and one normal correction toward the 0.5 level set |

The class header documents the driver operations and physics callbacks.
CMake collects this directory's `.cpp` files; see the
[root build guide](../../../README.md#build) for configuration and compilation.

## Linear-System Boundary

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

These groups are stored in `CrsMat::field_splits`. The shared
`ConfigureReducedPreconditioner()` owns compact numbering, nested setup,
and PCKSP operator selection; FSI supplies the physical grouping and PC choices.
See the [shared solver guide](../../../module/solver/README.md).

The outer FGMRES splits fluid components (0..3) from solid and multiplier
components (4..9). PETSc `selfp` eliminates the fluid diagonal when assembling
the preconditioning matrix for the coupled solid-multiplier block. A PCKSP
application runs FGMRES on that assembled approximation, retaining the full solid
operator and its multiplier coupling. It does not repeatedly solve the fluid
block inside this inner iteration.

The inner coupled preconditioner splits solid (4..6) and multiplier (7..9)
components with another PETSc `selfp` Schur split. Solid remains FGMRES with
rtol 0.1 and max 100, ASM overlap 1, and local PREONLY/ILU(1) with nonzero shift
1e-3. The new coupled iteration uses rtol 0.5 and max 100. These are
preconditioning solves; the outer and physical Newton tolerances control the
final correction.

`SolveSystem(NR_it, istep)` preserves the unscaled physical RHS, zeros the native
correction, scales the matrix/RHS and solves through `CrsMat`. It then multiplies
the recovered correction by `column_scale` and restores the physical RHS.
The outer KSP uses `rtol=1e-8` and `atol=1e-13` when any control point
across all ranks has both phase masses greater than `mtol` and a positive
`nlm_lump`. Otherwise it uses `rtol=1e-10` and `atol=1e-15`.
The outer iteration limit is 300.
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
| Mass activity | Absolute nodal mass strictly greater than `mtol` |
| Outer restart | 200 |
| Fluid tolerance / maximum / restart | 0.5 / 200 / 100 |
| Fluid orthogonalization | Classical Gram-Schmidt with refinement always |
| Velocity / pressure local ILU fill / ordering | 1 / RCM for both |
| Solid tolerance / maximum / restart | 0.1 / 100 / 30 |
| Solid orthogonalization | Classical Gram-Schmidt with refinement always |
| Coupled tolerance / maximum / restart | 0.5 / 100 / 100 |
| Coupled orthogonalization | Classical Gram-Schmidt; no explicit refinement option (PETSc default: never) |

The inner tolerances control preconditioning work. Outer linear convergence
and the four physical Newton criteria below determine correction acceptance.

The native CSR graph, full PETSc assembly, material parameters, and physical
Newton criteria remain independent of this preconditioning layer.

## Newton Iteration

`MPMMPMMonolithicFSI` is a `NonlinearProblem` (`module/implicit/newton_solver.h`), and
`SolveFSISystem(ts)` calls `SolveNewton(*this, ts, max_NR_it)`. Its `BeginNewtonStep()`
builds the lumped interface weights and active-component mask
once at the beginning of the time step. Fluid stabilization coefficients
are computed from the step-start nodal velocity and remain frozen during Newton.
Each Newton iteration (`AssembleNewtonSystem`, `CheckNewtonConvergence`,
`SolveNewtonCorrection`, `UpdateNRIncrement`) then performs:

1. Clear the coupled matrix and RHS.
2. Apply phase boundary values and compute trial nodal velocity/acceleration.
3. Assemble fluid, solid, and interface contributions once.
4. Check the current nonlinear RHS and interface endpoint-velocity difference.
5. If unconverged, solve a zero-initialized correction and update the unknowns.

The next iteration evaluates the updated state through the same assembly path.
`ComputeNRResidualNorms()` reads the unscaled, overlap-synchronized `b_rhs` for
fluid momentum, continuity, and solid momentum. The fourth field divides the
interface RHS by its positive active `nlm_lump` weight to measure the endpoint
velocity difference in m/s without interface-area weights.

The interface equation is `v_s^(n+1) = v_f^(n+1)`, using the endpoint velocities
already computed from the phase displacement iterates. It is not imposed at
`n+alpha_f`. The fluid and solid constraint-row Jacobians contain
`-nlm_lump_local * gamma_f/(beta_f*dt)` and
`+nlm_lump_local * gamma_s/(beta_s*dt)`, respectively. Row scaling uses these
velocity derivatives; multiplier-column scaling retains the phase force weights.

Each field uses an RMS over globally owned, unconstrained active components.
For the interface, the denominator is the active multiplier-component count,
not all `nodec` nodes or all mesh components. Shared components are counted
once. A field with no active components has zero RMS.

| Field | Residual | Acceptance threshold |
|-------|----------|----------------------|
| Fluid momentum | Nonlinear momentum RHS [N] | `max(1e-10, 1e-5 * initial_norm[0])` |
| Continuity | Nonlinear continuity RHS [m^3/s] | `max(1e-13, 1e-5 * initial_norm[1])` |
| Solid momentum | Nonlinear momentum RHS [N] | `max(1e-7, 1e-5 * initial_norm[2])` |
| Interface | Endpoint-velocity difference [m/s] | `1e-8` |

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

## Interface Activation

`LumpedLagrangeMultiplier()` interpolates the step-start solid and fluid volume
fractions and their gradients with `MakeSF()` at quadrature points. It integrates
`||phi_s * grad(phi_f) - phi_f * grad(phi_s)||` against the shape functions
to form `nlm_lump_local`. After copying those weights to `nlm_lump`,
`NodeVarComm(grid_.part, nlm_lump, 0)` sums contributions at shared control
points. The constraint-row Jacobian uses the local weight. The endpoint-velocity
RHS and activation use the synchronized weight.

`BuildActiveDOFs()` uses the absolute nodal-mass threshold directly:

```cpp
const bool fluid_active = fluid_.nmass[n] > mtol;
const bool solid_active = solid_.nmass[n] > mtol;
```

Masses equal to or below `mtol` are inactive. There is no phase-mean calculation
or relative `1e-4` cutoff in the current implementation.
A multiplier component requires both phases active and `nlm_lump[n] > 0`. Phase boundary
conditions and solid rigid constraints are then applied; a multiplier component
is fixed when both corresponding phase displacement components are fixed.

The coupled matrix assembly uses the per-component `fixed_dof` mask; it does
not call the generic `CrsMat::BuildActiveRowMask()`.

Changing `mtol` changes the phase activity tests. It does not create interface
weights at control points where `nlm_lump` is zero. The assembled multiplier
equations enforce nodal endpoint-velocity continuity; their convergence
alone does not establish particle nonpenetration.

## Particle Penetration Correction

The driver calls the following sequence after the coupled solve:

```cpp
fsi.solid_.Node2Particle(ts);
fsi.fluid_.Node2Particle(ts);
fsi.solid_.MoveParticle();
fsi.CorrectFluidPenetration();
fsi.fluid_.MoveParticle();
```

`Node2Particle(ts)` commits coordinates and kinematics, including the existing
particle shifting. `MoveParticle()` determines particle ownership and migrates
data; the fluid version also handles inflow. Correcting between the two migration
calls makes the final solid particles available on their owning ranks and lets
fluid ownership be determined from corrected coordinates.

After solid migration, `CorrectFluidPenetration()` maps the final solid particle
volumes to control points with `LocateLocalElement()`, `MakeSF()` and
`StandardVarP2G()`. It sums shared contributions with
`NodeVarComm(grid_.part, final_phi, 0)`, divides by `nvol`, and clamps the
solid phase coefficients to `[0, 1]`.
This field represents the solid after the particle update and supplies the
geometric interface at `phi_s = 0.5`.

`InterfaceSDF::Evaluate()` uses the same background shape functions to
interpolate `phi_s(x) = 0.5 + sum_i N_i(x) * (phi_s,i - 0.5)` and its
physical gradient `g = grad(phi_s)(x)`. Subtracting 0.5 before interpolation
preserves a constant 0.5 field. For a centre with `phi_s(x) > 0.5` and a
finite gradient with `||g|| > 1e-12 / min(dxy)`,
`InterfaceSDF::FindSurfacePoint()` applies one update along decreasing solid phase.

```text
x_corrected = clamp(x - (phi_s(x) - 0.5) * g / dot(g, g), xyminw, xymaxw)
```

The clamp bounds each corrected coordinate by the global mesh. The update
linearizes the phase field at the original centre. Curvature or clipping can
leave a residual at the corrected location. Centres on the fluid side, or
with an unusable gradient, keep their coordinates.

A rank evaluates centres inside its own elements locally. For a centre in the
global mesh that belongs to another rank, the current partition regions identify
the owner. An independent `ParticleCommunication` instance exchanges query
counts with `MPI_Alltoall` and sends coordinates with `PointVarSendrecv()`.
The owning rank evaluates each query and returns the resulting coordinates
in the same order. The originating rank writes them back, then the driver's
fluid `MoveParticle()` migrates full particle state using the resulting positions.

The correction changes fluid coordinates only. Velocity, acceleration,
TPIC/APIC state, pressure, mass and volume retain their transferred values.
It is a geometric remap separate from PST and the multiplier solve.
A particle's finite material domain can still overlap the solid.

Nodal fixed-DOF conditions constrain nodal kinematics. The geometric clamp
applies only to corrected centres. Fluid centres already outside the global
mesh are skipped by the correction and removed by `MoveParticle()`.
