# `module/solid/implicit/` — Implicit MPM Solid Solver

## Overview

This directory implements the **implicit Material Point Method (MPM) solid mechanics solver** with the currently available rigid, linear-elastic and hyperelastic material models. It uses a **Newmark-β / generalized-α time integrator** and **Newton–Raphson (NR) iteration**.

Each NR iteration assembles a sparse tangent system and solves it with PETSc KSP
or the native iterative solver. This is not a direct factorization of the full tangent matrix.

## File Inventory

| File                             | Responsibility                                                                                   |
| -------------------------------- | ------------------------------------------------------------------------------------------------ |
| `implicit_mpm_solid.h`         | Class declaration, constructor wiring (`SM_.owner_ = this`)                                    |
| `solve_solid_implicit.cpp`     | Newton steps, system assembly, tangent-modulus computation, diagonal preconditioner                   |
| `var_trans_solid_implicit.cpp` | P2G (`Particle2Node`), G2P (`Node2Particle`), deformation-gradient update, particle shifting |

## Class Architecture

```
SolidMaterialPointBase
    └── ImplicitSolidMPM
            └── SM_ : CrsMat  (ndof = 3, owner_ = this, FEM_flag = false)
```

`ImplicitSolidMPM` exposes its time-step entry points and assembly helpers publicly.

### Public Interface

- `DataInput()` — loads the standalone implicit-solid case input (orchestration file, parameter file, mesh/time scalars, BC data, particle data).
- `Particle2Node()` — P2G transfer for mass, volume, momentum, and force.
- `Node2Particle()` — G2P transfer and particle kinematic update.
- `SolveSolid()` — driver for one implicit time step.

### Assembly Helpers / Overrides

The following overrides are invoked polymorphically through the `MaterialPoint` base class or called internally by `SolveSolid`:

- `BuildPetscBCList(CrsMat& mat)` — collects 3-DOF BC global IDs (`ubc`, `vbc`, `wbc`, `rigid_bc`) for PETSc.
- `BCResidualSet(std::vector<double>& rr)` — zeroes constrained DOFs in a residual vector.
- `BCNRSet()` — sets Dirichlet values for the current NR iteration.
- `UpdateNRIncrement()` — applies the converged Newton correction `SM_.x_lhs` to `ndispl`.
- `AssembleSystem(...)` — assembles the tangent matrix `SM_.amat` and residual `SM_.b_rhs`.
- `ConfigurePreconditioner(CrsMat &mat, PC pc)` — selects standalone solid ASM/ILU(1).
- `ComputeTangentModulus(...)` — computes the 3×3 nodal stiffness block from the material tangent and stress state.

## Numerical Method

### Time Integration — Newmark-β / Generalized-α

`implicit_.integrator_` owns the generalized-alpha/Newmark coefficients; `implicit_`
is the `ImplicitMPM` component (`module/implicit/implicit_mpm.h`) owned by `ImplicitSolidMPM`.
`implicit_.ComputeNodeVelAccelFromDispl()` predicts nodal kinematics from the current
displacement iterate; `implicit_.CommitNodalKinematics()` commits the converged state.
`ComputeNRLumpedMassMat()` and `AddInertialForceToRHS()` stay virtual on
`ImplicitSolidMPM` so FSI solids can add their added-mass terms.

### Newton–Raphson Loop (`SolveSolid`)

`ImplicitSolidMPM` is a `NonlinearProblem` (`module/implicit/newton_solver.h`);
`SolveSolid(ts)` sets `max_NR_it` from `NR_flag` and calls
`SolveNewton(*this, ts, max_NR_it)`, which runs the steps in this order:

```cpp
BeginNewtonStep();                          // zero ndispl, snapshot stress_k_
for NR_it = 0 .. max_NR_it
    AssembleNewtonSystem(ts);               // BCNRSet, predictor, zero amat/b_rhs/x_lhs, AssembleSystem
    if CheckNewtonConvergence(NR_it, ksp_its) break; // NR_flag: SM_.CheckNRConvergence; else print only
    ksp_its = SolveNewtonCorrection(NR_it, ts); // SM_.SolveSystem: PETSc or native linear solve
    UpdateNRIncrement();                    // ndispl += x_lhs
EndNewtonStep();                            // release SM_ storage
```

- `NR_flag = true` → nonlinear elasticity (`NR_it = 0..100`, at most 101 passes)
- `NR_flag = false` → linear elasticity (0 NR iterations, single linear solve)

### Linear Solver and Reduced Equations

The standalone PETSc configuration is FGMRES, restart 200, with ASM overlap 1.
Local sub-solvers apply PREONLY + ILU(1), nonzero shift `1e-3`, through the
`-solid_field_sub_*` options.
Standalone outer defaults are `rtol=1e-10`, `atol=1e-15` and `max_it=1000`.

The constructor sets `SM_.ndof=3`. Assembly and native vectors retain the
full `nodec * 3` layout; `CrsMat` handles compact PETSc equations and solution
recovery after boundary elimination. Activity uses the absolute `nmass > mtol`
test. Each correction refreshes compact values and the numerical PC; layout
changes or scheduled rebuilding recreate the KSP. See the
[shared solver guide](../../solver/README.md) for compression and lifecycle details.

### System Assembly (`AssembleSystem`)

For each particle inside each element:

1. **Shape function & gradient** (`MakeSF`)
2. **Trial deformation-gradient update** (`UpdateDefGrad`) using `alpha_f`:
   `F_trial = (I + alpha_f * grad(Delta u)) * F_n`, with gradients in the step-reference configuration.
3. **Volume and constitutive update** from the trial determinant and deformation gradient.
4. **Implicit gradient correction** (`ImplicitDsfCorr`) for the deformed configuration.
5. **Stress evaluation**: `sts_af = stress_k[pid]` uses the constitutive result at
   the trial configuration directly; there is no separate linear blend of endpoint stresses.

Then, for each node pair `(ni, nj)`:

- **Inertial tangent** (lumped): `alpha_m * sf * mass / (beta_nb * dt^2)` in the base solid path.
- **Tangent stiffness** `K_mat = ∇N · C · ∇N^T` via `ComputeTangentModulus`
  - For nonlinear problems adds the **geometric stiffness** term `σ_af ⊗ I`
- **Internal / external force vectors**: `f_int = -∇N · σ_af · vol`, `f_ext = sf · (mass · g + trac)`

The assembled matrix is stored in `SM_.amat`, RHS in `SM_.b_rhs`.

The stiffness contribution is multiplied by `alpha_f * vol`. Inertia in the RHS
uses `alpha_m * a_trial + (1-alpha_m) * a_n`. `Node2Particle()` later evaluates
the endpoint deformation gradient with coefficient 1 and commits particle state.
Assembly writes the member determinant and volume even while deformation-gradient
and stress outputs are temporary; these members are not isolated trial-state storage.

For nonlinear solves, `CrsMat::CheckNRConvergence()` checks the assembled `b_rhs`
before solving. After an increment, the next pass reassembles at the updated
displacement; each pass assembles only once. Iteration zero saves the reference
norm. Linear solves bypass this check and perform one solve. For norm and
threshold caveats, see the
[solver limitations](../../solver/README.md#current-convergence-limitations).

### Tangent Modulus (`ComputeTangentModulus`)

Computes the 3×3 nodal stiffness block from the 4th-order material tangent `C` (stored in `cm_.stif_mat`) and the current stress state:

```
K_ij = Σ_k,l  dsf[k][ni] · C[i][k][j][l] · dsf[l][nj]
```

If `NR_flag` is true, the geometric stiffness `σ_af[j][l] · δ_ik` is added.

### Variable Transfer (`var_trans_solid_implicit.cpp`)

#### `Particle2Node`

- P2G: mass, volume, momentum, force → `nmass`, `nvof`, `nmome`, `nforce`
- Nodal velocity / acceleration: `nvel = nmome / nmass`, `naccel = nforce / nmass`
- Apply BC values (`BCSetVal`) and zero out constrained accelerations (`BCSetZero`)

#### `Node2Particle`

- Predict velocity / acceleration with Newmark-β
- Commit nodal kinematics (`implicit_.CommitNodalKinematics`)
- G2P: update particle displacement gradient (`UpdateDefGrad`), stress, velocity, acceleration
- Update particle position: `coord += disp + disp_corr`
- Particle shifting (`DeltaCorrectionPST`) for uniform distribution
