# Partitioned MPM-FEM FSI

`MPMFEMBlockFSI` couples a two-phase stabilized FEM fluid (`StabilizedFEM`,
liquid/gas on the background mesh) with an implicit MPM solid
(`ImplicitSolidMPM`) on the same background grid. The coupling is partitioned
and strong: in every time step the two sub-solvers are iterated in a block
Gauss-Seidel loop until the fluid and solid velocities agree on the interface
nodes.

- Fluid side: the solid interface velocity is imposed as a Dirichlet condition.
- Solid side: the fluid stress enters as a nodal force, plus an added-mass term.
- The interface velocity is relaxed with depth-one Anderson acceleration.

Each sub-solver keeps its own linear system: fluid `NS_` (FEM, 4 components)
and solid `SM_` (3 components). There is no monolithic matrix; compare the
[MPM-MPM FSI](../MPM_MPM/README.md) for the monolithic alternative.

## Source responsibilities

| File | Responsibility |
|------|----------------|
| `main.cpp` | Time-step order, solid transfers and motion, output and restart |
| `datain_para.cpp` | Coupled case parameters, grid/BC input for both phases, solid particles |
| `block_fsi.h` | `MPMFEMBlockFSI` coordinator and the `FSIFluid` / `FSISolid` overrides |
| `block_fsi.cpp` | Interface detection, block iteration, FSI force and added mass, solid overrides, Turek drag/lift |
| `bc_setting.cpp` | `FSIFluid` interface velocity BC (native and PETSc paths) |
| `../../../module/relaxation.h` / `.cpp` | Anderson (depth 1 and 2) and Aitken interface relaxation |

The sub-solvers are documented in the [fluid FEM](../../../module/fluid/FEM/README.md)
and [implicit solid](../../../module/solid/implicit/README.md) guides.

## Time step

`main.cpp` advances each step in this order:

```text
fluid_.SetPFDomain()          // blend element density/viscosity from the liquid fraction nphi
solid_.MeshPointLinklist()
solid_.Particle2Node()        // solid P2G, including solid volume fraction nvof
SolveFSISystem()              // block iteration (below)
solid_.Node2Particle()
solid_.MoveParticle()
output / restart every iout steps
```

The fluid is Eulerian: its liquid fraction is projected once at start-up by
`Particle2NodePhi()`, and no fluid particles move. The phase-field system `PF_`
is reserved and not advanced (see the fluid FEM guide).

## Block iteration (`SolveFSISystem`)

```text
DetectFSIInterface()       // solid nphi = clamp(nvof / nvol, 0, 1); interface nodes: nphi > phi_cut
UpdateFSIInterfaceBC()     // initial interface velocity = solid step-start nodal velocity
force a PC rebuild of NS_ and SM_
for block_it = 0 .. max_block_iter
    fluid_.SolveNS()                        // interface velocity imposed by FSIFluid
    CalFSIForce()                           // added mass and fluid stress force on the solid
    solid_.SolveSolid()                     // Newton with the FSISolid overrides
    v_s = solid endpoint velocity from the displacement increment
    Anderson_relaxation_M1(...)             // update the imposed interface velocity
    residual of v_s - v_f on the interface; converged -> break
    block_it == max_block_iter -> MPI_Abort
fluid_.naccel = ComputeNodeAccelFromVel(); fluid_.UpdateNodalVar()
```

### Fluid side (`FSIFluid`)

`BCSet()` applies the regular fluid BCs, then sets `nvel` at the interface
nodes to the relaxed interface velocity `fsi_intf.fbc`. `BuildPetscBCList()` and
`BCResidualSet()` add the interface velocity components to the constrained
rows, so the PETSc and native paths treat them as Dirichlet values.

### FSI force and added mass (`CalFSIForce`)

- Added mass: `added_mass[n] = sum over solid particles of sf_n * rhol * vol_p`.
- Fluid stress at `n + alpha_f`: `sigma = mu (grad v + grad v^T) - p I`, with
  `mu` blended from `rmul`/`rmug` by `nphi`.
- The force is integrated over the fluid quadrature points with weight
  `cell volume / (npxye_x * npxye_y * npxye_z)`, multiplied by the solid `nphi` of
  the node, and set to zero where `nphi < phi_cut`.

### Solid side (`FSISolid`)

| Override | Change to the implicit solid |
|----------|------------------------------|
| `ComputeNRLumpedMassMat` | adds `da_du * sf * rhol * vol_p` (added-mass tangent) |
| `ComputeExternalForce` | body force on `mass_p + vol_p * rhol` instead of `mass_p` |
| `AddInertialForceToRHS` | inertia with `nmass + added_mass`; subtracts `nfsi_force` |

The solid Newton loop itself is the shared one (`SolveNewton`, see the
[implicit core](../../../module/implicit/README.md)).

### Relaxation and convergence

`Anderson_relaxation_M1` uses the interface residual `R = v_s - v_bar`, where
`v_bar` is the currently imposed interface velocity:

```text
lambda = (dR . R) / (dR . dR)        dR = R_k - R_k-1, dbc-weighted, global sums
lambda = 0.5 on the first iteration
v_bar  = v_s - lambda * (v_s - v_s,old)
```

`Aitken_relaxation` (coefficient clamped to 0.05..0.8) and
`Anderson_relaxation_M2` (depth two) are kept in `block_fsi.cpp` as commented
alternatives.

`CalFSIResidual` measures `v_s - v_f` on the interface nodes:
`rkr = sqrt(sum dr^2 dbc)` and `rms = sqrt(sum dr^2 dbc / n_intf)`, with the
relative value `rkr / r0r`, where `r0r` comes from `block_it == 0`.
The block loop stops when any of these holds:

- `rkr / r0r < tol_ref`;
- `r0r < tol_ref`;
- `rms < tol_abs`.

Rank 0 prints `Block_converge: block_it  rkr/r0r  r0r  rms`. Failure at
`max_block_iter` prints `Block iteration can not converge` and aborts.

| Parameter (`block_fsi.h`) | Value | Meaning |
|---------------------------|-------|---------|
| `max_block_iter` | 100 | last block iteration |
| `phi_cut` | 0.25 | interface node threshold on the solid volume fraction |
| `tol_ref` | 1e-3 | relative interface residual (also applied to `r0r`) |
| `tol_abs` | 1e-6 | RMS interface velocity difference [m/s] |
| `relax_omega` | 0.1 | Aitken only |

## Linear solvers

Both systems use the shared `CrsMat` path; see the
[shared solver guide](../../../module/solver/README.md).

- Fluid `NS_`: velocity/pressure Schur FieldSplit with BoomerAMG.
- Solid `SM_`: ASM with local ILU(1).

`SolveFSISystem()` forces a PC rebuild of both at the start of each time step.

## Input and run

- Parameters: `data/generate/fsi/MPM_FEM/input.txt`. The blocks are:
  1. `solswitch rstflag nlstep`
  2. `ista iend iout`
  3. `dt mtol wfem.spec_rad sp.spec_rad`
  4. `xymin`, `xymax`, `xyelem`, `idimc`, `npxye`
  5. `sp.rho rhol rmul rhog rmug fs_height`
  6. the materials
  7. `bb`
- Data: build `makinput_fsi_mpm_fem` and `makdivide_fsi_mpm_fem`.
  - The grid file carries the grid, the solid BCs and the fluid BCs.
  - The point file carries the solid particles only. The fluid mesh data comes
    from `InitializeMeshData()`.
- Run from `build/` with the driver's `file.dat` in `build/work/src_fsi/MPM_FEM/`:
  `sh run.sh NP work/src_fsi/MPM_FEM`.
- Output: solid points via `OutputPointDataVTKHDF`; fluid mesh fields
  (velocity, pressure, phase) via `OutputMeshDataVTKHDF`.
- Restart: both phases write restart data at output steps when `rstflag` is 2 or 3.

`CalDragLiftCoeffForTurekCFD()` integrates the cylinder surface traction for the
Turek benchmark; its call in `SolveFSISystem()` is commented out.
