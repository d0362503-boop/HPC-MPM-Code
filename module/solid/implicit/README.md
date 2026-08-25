# `solid/implicit/` — Implicit MPM Solid Solver

## Overview

This directory implements the **implicit Material Point Method (MPM) solid mechanics solver** for large-deformation hyperelastic/elastoplastic problems.  It uses a **Newmark-β / generalized-α time integrator** and solves the resulting nonlinear system with **Newton–Raphson (NR) iteration**.

Unlike the explicit MPM counterpart, this solver assembles and factorizes a sparse tangent-stiffness matrix every NR iteration.

| Matrix  | DOF         | Purpose                        |
| ------- | ----------- | ------------------------------ |
| `SM_` | 3 (u, v, w) | Momentum / displacement system |

---

## File Inventory

| File                             | Responsibility                                                                                   |
| -------------------------------- | ------------------------------------------------------------------------------------------------ |
| `implicit_mpm_solid.h`         | Class declaration, constructor wiring (`SM_.owner_ = this`)                                    |
| `solve_solid_implicit.cpp`     | NR loop, system assembly, tangent-modulus computation, diagonal preconditioner                   |
| `var_trans_solid_implicit.cpp` | P2G (`Particle2Node`), G2P (`Node2Particle`), deformation-gradient update, particle shifting |

---

## Class Architecture

```
SolidMaterialPointBase
    └── ImplicitSolidMPM
            └── SM_ : CrsMat  (ndof = 3, owner_ = this, FEM_flag = false)
```

`ImplicitSolidMPM` exposes the standard solver driver hooks and keeps the Newton–Raphson internals private.

### Public interface

- `DataInput()` — loads the standalone implicit-solid case input (orchestration file, parameter file, mesh/time scalars, BC data, particle data).
- `Particle2Node()` — P2G transfer for mass, volume, momentum, and force.
- `Node2Particle()` — G2P transfer and particle kinematic update.
- `SolveSolid()` — driver for one implicit time step.

### Private helpers / overrides

The following overrides are invoked polymorphically through the `MaterialPoint` base class or called internally by `SolveSolid`:

- `BuildPetscBCList(CrsMat& mat)` — collects 3-DOF BC global IDs (`ubc`, `vbc`, `wbc`, `rigid_bc`) for PETSc.
- `BCResidualSet(std::vector<double>& rr)` — zeroes constrained DOFs in a residual vector.
- `BCNRSet()` — sets Dirichlet values for the current NR iteration.
- `UpdateNRIncrement()` — applies the converged Newton correction `SM_.x_lhs` to `ndispl`.
- `AssembleSystem(...)` — assembles the tangent matrix `SM_.amat` and residual `SM_.b_rhs`.
- `ComputeTangentModulus(...)` — computes the 3×3 nodal stiffness block from the material tangent and stress state.

---

## Numerical Method

### 1. Time Integration — Newmark-β / Generalized-α

The solver inherits `ComputeNodeVelAccelFromDispl` and `CommitNodalKinematics` from `MaterialPoint` (shared with the fluid solver).  At the beginning of each time step:

```cpp
ComputeNodeVelAccelFromDispl(nvel_k, naccel_k);  // predictor
...
CommitNodalKinematics(nvel_k, naccel_k);          // commit converged state
```

Parameters (`alpha_f`, `alpha_m`, `gamma_nb`, `beta_nb`) are set via `GeneralizedAlphaParaSet` / `NewmarkBetaParaSet` in the base class.

### 2. Newton–Raphson Loop (`SolveSolid`)

```cpp
for NR_it = 0 .. iter_max
    BCNRSet();                              // apply Dirichlet values
    ComputeNodeVelAccelFromDispl(...);     // predictor step
    AssembleSystem(naccel_k, nvel_k, stress_k);
    iter = SolveSystem(SM_, NR_it);         // PETSc or native linear solve
    UpdateNRIncrement();                    // ndispl += x_lhs
    if CheckNRConvergence() break;
```

- `NR_flag = true` → nonlinear elasticity (up to 1000 NR iterations)
- `NR_flag = false` → linear elasticity (0 NR iterations, single linear solve)

### 3. System Assembly (`AssembleSystem`)

For each particle inside each element:

1. **Shape function & gradient** (`MakSf`)
2. **Deformation-gradient update** (`UpdateDefGrad`) — `F^{n+1} = (I + ∆u) · F^n`
3. **Constitutive model update** (`UpdateConstitutiveModel`) — stress `S` and tangent modulus `C`
4. **Implicit gradient correction** (`ImplicitDsfCorr`) — background-grid correction for MPM
5. **Stress interpolation** to the `α_f` time level: `σ_af = α_f · σ_k + (1-α_f) · σ_n`

Then, for each node pair `(ni, nj)`:

- **Mass matrix** (lumped): `M = γ · sf · mass`
- **Tangent stiffness** `K_mat = ∇N · C · ∇N^T` via `ComputeTangentModulus`
  - For nonlinear problems adds the **geometric stiffness** term `σ_af ⊗ I`
- **Internal / external force vectors**: `f_int = -∇N · σ_af · vol`, `f_ext = sf · (mass · g + trac)`

The assembled matrix is stored in `SM_.amat`, RHS in `SM_.b_rhs`.

### 4. Tangent Modulus (`ComputeTangentModulus`)

Computes the 3×3 nodal stiffness block from the 4th-order material tangent `C` (stored in `cm_.stif_mat`) and the current stress state:

```
K_ij = Σ_k,l  dsf[k][ni] · C[i][k][j][l] · dsf[l][nj]
```

If `NR_flag` is true, the geometric stiffness `σ_af[j][l] · δ_ik` is added.

### 5. Variable Transfer (`var_trans_solid_implicit.cpp`)

#### `Particle2Node`

- P2G: mass, volume, momentum, force → `nmass`, `nvof`, `nmome`, `nforce`
- Nodal velocity / acceleration: `nvel = nmome / nmass`, `naccel = nforce / nmass`
- Apply BC values (`BCSetVal`) and zero out constrained accelerations (`BCSetZero`)

#### `Node2Particle`

- Predict velocity / acceleration with Newmark-β
- Commit nodal kinematics (`CommitNodalKinematics`)
- G2P: update particle displacement gradient (`UpdateDefGrad`), stress, velocity, acceleration
- Update particle position: `coord += disp + disp_corr`
- Particle shifting (`DeltaCorrectionParticleShifting`) for uniform distribution

---

## Integration with the Rest of the Codebase

- **`CrsMat`** (`module/solver/`) — generic sparse-matrix wrapper; BC logic injected through `owner_` virtuals. `SM_.ndof = 3` configures the generic solver for 3-DOF solid mechanics.
- **`MaterialPoint`** (`module/material_point.h`) — base class providing `ComputeNodeVelAccelFromDispl`, `CommitNodalKinematics`, and the virtual BC hooks.
- **`SolidMaterialPointBase`** (`module/solid/solid_material_point.h`) — intermediate base providing constitutive-model interface (`UpdateConstitutiveModel`).
- **Solvers** (`module/solver/solver.cpp`) — `GPBiCGSafe` called on `SM_`; residual callbacks dispatched through `SM_.owner_->BCResidualSet(rr)`.

---

## Coding Style Notes

- Trailing underscore for class members (`owner_`, `SM_`).
- `using namespace implicitmpm;` appears in `.cpp` files only.
- The constructor enables PETSc (`use_petsc = true`) and sets `FEM_flag = false` to indicate an MPM (not FEM) matrix structure.
