# `module/solver/` - Linear Solver and PETSc Wrapper

> This note records the current `CrsMat` design, the PETSc/HYPRE bridge, and the production solver configuration used by the 32-rank FSI runs.

---

## 1. What Lives Here

`module/solver/` contains two layers:

1. `CrsMat`
   A custom CSR-style sparse matrix wrapper that stores:
   - `amat`
   - `b_rhs`
   - `x_lhs`
   - `matrow`
   - `matcolid`

2. PETSc / HYPRE bridge
   When `use_petsc = true`, the same `CrsMat` object also owns:
   - PETSc `Mat`
   - PETSc `Vec`
   - PETSc `KSP`
   - local-to-global mappings
   - BC index lists
   - local scatter objects for solution recovery

Main files:

| File | Role |
|------|------|
| `crsmat.h` | `CrsMat` definition |
| `crsmat.cpp` | PETSc setup, assembly, solve, BC handling |
| `solver.cpp` | native fallback solvers (`GPBiCG`, `GPBiCGSafe`, `GPBiCGAR`) |

---

## 2. Indexing, Layout, and CRS Structure

### 2.1 Local layout

Local unknowns are stored in component-major layout:

```cpp
int local_idx = natural_node_id + var * nodec;
```

Examples:
- fluid `ndof = 4`: `u, v, w, p`
- solid `ndof = 3`: `u, v, w`

### 2.2 Global PETSc numbering

Global PETSc scalar IDs are node-major:

```cpp
int global_id = global_node * ndof + var;
```

This mixed layout is legacy, but it is the current working design.

### 2.3 Critical rule

All solver-side DOFs use `nodec`, never `node`.

`node` is only for visualization / output. Using it in PETSc ownership or local row counts will break matrix sizes and assembly.

### 2.4 CRS block storage

`BuildCrsMat(num_block)` builds a static CSR graph from the background grid stencil implied by `idimc` (B-spline order). Each control point gets a dense neighbor list in every direction:

```
row idn  →  columns [idn_min ... idn_max]  (inclusive)
```

`matrow` stores row offsets, `matcolid` stores column indices, and `amat` stores the actual values in block-major order:

```cpp
amat[j + block_id[row_var * ndof + col_var]]
```

where `num_block = ndof * ndof`. For fluid (`ndof=4`) this means 16 scalar blocks per CSR entry; for solid (`ndof=3`) it is 9. This layout allows the assembly loop to fetch an entire dense `ndof×ndof` block with a single base offset.

### 2.5 Active vs inactive nodes (MPM)

In MPM, a control point is called *inactive* when its assembled matrix row carries negligible physical content. This is detected at solve time by `BuildActiveRowMask()`, which computes the absolute sum of all entries in the node's CSR row. If the sum is below `mtol`, the node is marked inactive.

**Critical invariant for parallel PETSc solves:** the active/inactive decision must be synchronized across overlap control points before assembly skips rows or inserts owned-row identity blocks. A shared control point is considered active if *any* overlapping rank assembled nontrivial row content for that node. `BuildActiveRowMask()` implements this by computing a local integer indicator, calling `NodeVarComm(..., 0)` to accumulate indicators on shared nodes, and rebuilding `active_row_mask` from the synchronized result.

This prevents owner-rank false negatives: a rank that happens to own a shared node but has no local element contribution would otherwise mark the row inactive and identity-fill it, producing an artificial Dirichlet-like wall along partition boundaries.

During assembly:
- inactive nodes are skipped in the first pass
- locally-owned inactive nodes receive an identity block in a second pass to keep the matrix well-conditioned

This replaces the older `nmass < mtol` heuristic with a matrix-content-based check that is consistent with the actual assembled system. See Section 3.5 for details.

> **Do not** classify active rows from purely rank-local `amat` content and immediately apply identity fill on owned rows. In overlapping decompositions this misclassifies shared rows whose physical contributions are split across ranks.

---

## 3. `CrsMat` Workflow

The PETSc path is:

1. `BuildCrsMat(num_block)`
2. `BuildLGMAP(ndof)`
3. `BuildPetscMat(ndof)`
4. `BuildKSPSolver()`
5. per solve:
   - `AssemblePetscMat(ndof)`
   - `UpdatePetscRhs(ndof)`
   - `SolveWithPetsc(ndof, NR_it)`

With `use_petsc = true`, `BuildCrsMat` first calls `ResetPetscSolver()` to release the
previous PETSc objects (`Mat`, `Vec`, `KSP`, lgmaps, scatter). This makes a matrix
rebuild safe after a DLB repartition changes `nodec` and the ownership layout
(see `module/DLB/README.md`).

### 3.1 Ownership

Shared control points are owned by a tie-break based on `aelemmin`.

`interior_list` stores locally owned control points.

Ghost control points are still assembled, because remote element contributions may live there before PETSc redistributes them during `MatAssemblyEnd`.

**Resolved: LGMAP size and PETSc ownership consistency.** `BuildLGMAP()` now builds explicit owned-first PETSc-local maps through `BuildPetscLocalMaps()`. `petsc_local_to_natural` places owned control points first, followed by ghost control points, and `petsc_local_block_gids` maps these `nodec` PETSc-local indices to global node IDs. `BuildPetscMat()` creates the matrix with local size `local_node * ndof` and attaches an LGMAP of size `nodec` (owned + ghost). PETSc interprets the first `local_node * ndof` entries as owned rows and uses the remaining entries as ghost/assembly indices. `CheckOwnershipMetadata()` validates that `owned_natural_ids`, `ghost_natural_ids`, and the local-to-natural map sum to the expected global control-point count.

### 3.2 Matrix setup

`BuildPetscMat()` currently uses:

```cpp
MatSetType(mat, MATAIJ);
MatSetBlockSize(mat, ndof);
MatMPIAIJSetPreallocation(...);
MatSeqAIJSetPreallocation(...);
MatSetOption(mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
```

`MatSetBlockSize(ndof)` is mandatory. Without it, BoomerAMG quality drops sharply.

`MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE` is **required for MPM**. See Section 3.5.

### 3.3 Blocked assembly path

`AssemblePetscMat()` treats each natural node as one block row and each CSR neighbour as one dense `ndof×ndof` block. The packing order is row-major:

```
(0,0) (0,1) ... (0,ndof-1)
(1,0) (1,1) ... (1,ndof-1)
...
```

All reusable buffers (`petsc_cols_buf`, `petsc_block_vals_buf`) are pre-allocated in `BuildPetscMat()` to their maximum size (`max_row_nnz * ndof`). This avoids repeated heap allocation in the hot loop.

### 3.4 Fake structure assembly

The code still performs one zero-value fake assembly during setup to lock the exact sparsity pattern.

This is not the most elegant approach, but it is stable with the current mixed local/global indexing.

### 3.5 MPM-specific matrix behavior

#### Dynamic nonzeros from moving particles

`BuildCrsMat()` builds a **static** CSR graph based on the background control-point stencil (`±idimc` in each direction). This graph is fixed for the entire simulation.

However, in MPM the particles move. When a particle crosses element boundaries, it can create **new coupling** between control points that were not connected in the initial static stencil. These new couplings appear as **new nonzero entries** during `MatAssemblyEnd()`.

By default, PETSc treats new nonzeros as a hard error:

```
PETSC ERROR: New nonzero at (12540,510) caused a malloc.
```

The fix is:

```cpp
MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
```

This tells PETSc to dynamically `malloc` storage for unexpected nonzeros instead of aborting. There is a small performance penalty, but it is the only robust way to handle particle-driven sparsity changes without over-allocating the entire matrix.

> **Tradeoff:** A more aggressive fix would be to enlarge the initial stencil in `BuildCrsMat()` to cover the maximum possible particle influence radius, but this would waste memory for the majority of rows that never see distant particles. The current dynamic-allow approach is simpler and memory-efficient.

#### Inactive-node identity fill

At the start of each assembly, `BuildActiveRowMask()` scans every control point's CSR row and computes a local active indicator. The indicator is then synchronized across overlap nodes via `NodeVarComm` so that a shared node is active if any overlapping rank has nonzero row content. FEM nodes are always active. Nodes that are still inactive after synchronization are skipped in the first assembly pass to avoid adding zero rows. For locally-owned inactive nodes, a second pass inserts an identity block:

```cpp
MatSetValuesBlocked(..., 1, &block_row, 1, &block_row, identity_block.data(), ADD_VALUES);
```

This keeps the matrix non-singular and prevents the linear solver from diverging on empty-space degrees of freedom, while ensuring that partition-boundary nodes with real physical contributions are not pinned by mistake.

#### Residual consistency for inactive MPM nodes

For the PETSc path, the Newton residual monitor must follow the **same active/inactive rule**
as the matrix assembly.

Current implementation (`ComputePetscResidualStats`):

- `AssemblePetscMat()` calls `BuildActiveRowMask()` at the start of each assembly,
  then skips node block rows where `!FEM_flag && active_row_mask[i] == 0`
- locally-owned skipped rows receive an identity block in the second pass
- `ComputeRefResidual()` and `ComputeAbsResidual()` now use `ComputePetscResidualStats()` for
  the PETSc path, which:
  1. Computes the residual vector via `MatMult(petsc_mat, petsc_x, residual)` and
     `VecAYPX(residual, -1.0, petsc_b)` (i.e. `r = b - A·x`)
  2. Builds a PETSc `active_mask` vector that zeros out inactive DOFs
  3. Applies the mask with `VecPointwiseMult` and computes the L2 norm via `VecNorm`
  4. Returns the active DOF count via `VecSum`

This replaces the older hand-rolled `MatVecMult` + `dbc`-weighted loop. The new path uses the
same PETSc operator that KSP sees, eliminating operator inconsistency between the solver and
the convergence monitor. `ComputeAbsResidual` normalizes by `sqrt(active_dof)` (unweighted RMS).

> **Note:** The native (non-PETSc) path still uses the original `adiag`-based filter and
> `dbc`-weighted norm. The two paths are intentionally kept independent.

---

## 4. Boundary Conditions

BC rows are handled by:

```cpp
MatZeroRowsColumns(...)
```

`CrsMat` does not hardcode BC lists itself. It asks its owner object:

```cpp
owner_->BuildPetscBCList(*this);
```

Current owner responsibilities:

| Physics | BC contributors |
|---------|-----------------|
| Fluid | `ubc`, `vbc`, `wbc`, `pbc`, `fsi_intf` |
| Solid | `ubc`, `vbc`, `wbc`, `rigid_bc` |

Because `fsi_intf` changes every step, BC lists must be rebuilt before every solve.

---

## 5. Current PETSc Solver Configuration

Monolithic AMG is the default preconditioner stack for PETSc-based `CrsMat` objects with `use_schur_fieldsplit = false`. The stabilized Navier--Stokes systems in both `StabilizedMPM` and `StabilizedFEM` set it to `true` and use the Schur field split described in Section 5.1a. The FEM phase-field system and the implicit-solid system keep the default `false` value.

### 5.1 Unified AMG preconditioner

- outer Krylov: `KSPFGMRES` (overridable via `-ksp_type`)
- preconditioner: `PCHYPRE`
- AMG type: `BoomerAMG`

This applies to every PETSc-based `CrsMat` with `use_schur_fieldsplit = false`, regardless of `ndof` (1, 3, or 4).

The code explicitly sets the following PETSc choices:

```cpp
KSPFGMRES
PCHYPRE + BoomerAMG
rtol = 1.0e-12, atol = 1.0e-15, max_it = 1000
```

The Schur field split sets its child-PC defaults with `PetscOptionsSetValue(...)` before
`KSPSetFromOptions()`. See Section 5.1a for the resulting PETSc option names.

### 5.1a Schur field split for stabilized fluid systems

`StabilizedMPM::NS_` and `StabilizedFEM::NS_` enable `use_schur_fieldsplit` in their
constructors. `ConfigurePreconditioner()` then configures:

```cpp
PCFIELDSPLIT with block size ndof (= 4)
velocity field: block components {0, 1, 2}; pressure field: {3}
PC_COMPOSITE_SCHUR + PC_FIELDSPLIT_SCHUR_FACT_LOWER
Schur preconditioner: PC_FIELDSPLIT_SCHUR_PRE_SELFP
both sub-solves: preonly + HYPRE BoomerAMG
```

The current split layout is specific to the four-DOF stabilized fluid system: the
`"velocity"` split contains components 0--2 and the `"pressure"` split contains component 3.
`ConfigurePreconditioner()` enforces this at runtime: with any `ndof` other than 4 it
falls back to plain BoomerAMG even when `use_schur_fieldsplit` is set.
The split names determine the PETSc option middle component, for example:

```text
-fieldsplit_velocity_pc_type hypre
-fieldsplit_pressure_pc_type hypre
```

`-fieldsplit_` and `_pc_type` are PETSc syntax; only `velocity` and `pressure` come from the
names passed to `PCFieldSplitSetFields()`.

This split exists because the stabilized (u,p) rows differ strongly in scale, and a single
monolithic AMG over the coupled block was not robust for the dam case. The split gives
velocity and pressure their own AMG hierarchies while `SELFP` builds the Schur
approximation for the stabilized `K_pp` term. The field split owns its child PCs, so the
AMG rebuild path does not call `PCReset`: `KSPSetOperators()` supplies the updated matrix and
`PCSetReusePreconditioner(pc, PETSC_FALSE)` makes PETSc rebuild the child preconditioners during
the next setup while retaining the split configuration.

### 5.2 Solver independence

Fluid and solid keep their AMG hierarchies independently.

There is no longer any global "release other AMG" behavior. The old fieldsplit flags (`use_fieldsplit`, `pressure_pc_use_amg`) were removed; the current Schur split is controlled solely by `use_schur_fieldsplit`. Fluid and solid AMG are allowed to coexist.

---

## 6. AMG Rebuild Policy

Each `CrsMat` controls its own rebuild cadence through:

```cpp
amg_rebuild_freq
```

Current production values:

| System | Value |
|--------|-------|
| MPM fluid `NS_` | `1` (PETSc + Schur field split is the default) |
| FEM fluid `NS_` / `PF_` | `20` |
| Implicit solid `SM_` | `1` |

Runtime logic in `SolveWithPetsc()`:

- periodic rebuild when `istep % amg_rebuild_freq == 0`
- at most one rebuild per time step: `NR_it > 0` always reuses the current
  preconditioner, so the AMG hierarchy is built on the first Newton iteration
  and lagged for the rest of the step
- optional forced rebuild when iterations deteriorate strongly: a solve whose
  KSP iteration count more than doubles versus the previous solve requests a
  rebuild on the next solve (`force_rebuild_next_`), including later Newton
  iterations
- rebuild is local to the current solver only

Implementation detail for a monolithic AMG preconditioner:

```cpp
if (need_rebuild) {
    PCReset(pc);
    KSPSetOperators(...);
    PCSetReusePreconditioner(pc, PETSC_FALSE);
} else {
    PCSetReusePreconditioner(pc, PETSC_TRUE);
}
```

This avoids keeping "old AMG + new AMG" simultaneously inside one solver for too long, while still allowing fluid and solid preconditioners to coexist.

For `use_schur_fieldsplit = true`, the same rebuild decision updates the operator and disables
preconditioner reuse, but deliberately skips `PCReset` and reconfiguration so that the existing
FieldSplit child-PC structure remains intact.

The `NR_it > 0` suppression yields to a pending forced rebuild:

```cpp
if (NR_it > 0 && !force_rebuild) { need_rebuild = false; }
```

Lagging the preconditioner within a time step is safe because the matrix
structure is fixed and its values only drift with the Newton iterate; the 2x
iteration-growth rule catches the cases where the lagged hierarchy actually
deteriorates. (The older proactive rule — rebuild whenever the previous solve
needed more than a few iterations — effectively rebuilt every Newton iteration
for the MPM fluid and has been removed.)

---

## 7. Solution Recovery

Older PETSc code gathered the full global solution onto every rank with `VecScatterCreateToAll`.

That was too memory-heavy for the 32-rank FSI runs.

Current behavior:

- each rank builds a scatter for only its own `nodec * ndof` entries
- PETSc solution is copied into a sequential local buffer `seq_x`
- `x_lhs` is filled from that local sequential vector

This reduces per-rank memory pressure significantly and was one of the key changes that made the 32-rank 10-step run stable.

---

## 8. Solid Memory Lifecycle

One important practical fix is in the solid solve path:

- after each solid solve, large temporary vectors are released again

Specifically:

- `SM_.adiag`
- `SM_.amat`
- `SM_.b_rhs`
- `SM_.x_lhs`

Without this release, solid work arrays remained resident across time steps and pushed the next fluid solve into OOM territory.

---

## 9. Current 32-Rank Baseline

### 9.1 Historical monolithic AMG baseline

The earlier Turek baseline used a monolithic AMG preconditioner stack (`PCHYPRE` + `BoomerAMG`)
for all PETSc-based solves. For the 10-step case:

- `NS_iter ~= 8`

This is the reference used to judge whether future preconditioner changes are acceptable.

### 9.2 Current recommended production run

Current configured default:

- stabilized fluid `NS_`: `Schur field split + FGMRES`
- FEM phase field / implicit solid: `monolithic AMG + FGMRES`
- local PETSc solution scatter
- solid temporary arrays released every step

Successful 16-rank cavity-flow run after tightening tolerances:

- fluid `NS_iter`: `6~10`, stabilizing around `6`
- FSI block iteration converges in `0~2` iterations per step

### 9.3 Default Krylov Policy

Current default policy in code is simple:

- all PETSc KSP objects default to `KSPFGMRES`
- stabilized MPM and FEM fluid `NS_` systems use the Schur field split (Section 5.1a)
- the FEM phase-field and implicit-solid systems use monolithic `BoomerAMG`
- remaining distinctions are the rebuild frequency and the physical BCs injected by each physics object

So the remaining distinction between systems is mainly in the physics (DOF count, BCs, rebuild cadence, and whether the system is a stabilized fluid `NS_` solve), not in the outer Krylov family.

---

## 10. Known Good / Known Bad Directions

### 10.1 Good directions

- keep `MatSetBlockSize(ndof)`
- keep fluid and solid AMG independent
- keep solid temporary arrays released after solve
- use local solution scatter, not global `CreateToAll`
- if mixed u-p problems show false convergence (very few or zero iterations), tighten `rtol`/`abstol` or add symmetric diagonal scaling before the PETSc solve

### 10.2 Bad or abandoned directions

- restoring global cross-solver AMG release logic
- gathering the full global PETSc solution on every rank
- using a loosely set `rtol` (e.g. `1.0e-8`) for unscaled mixed u-p systems, which can accept a trivial initial guess as "converged"
- treating the current solid matrix as safely CG-compatible without further proof
- editing headers without reconfiguring CMake afterward

---

## 11. Practical Notes for Future Changes

1. If you change any PETSc-related header path or solver header, reconfigure and rebuild: `rm -rf build/CMakeCache.txt build/CMakeFiles && cmake -S . -B build && cmake --build build -j8`.
2. If you change `nxyr`, you must regenerate partitions before solver benchmarking.
3. If a new preconditioner looks faster but `NS_iter` jumps far above the unified AMG baseline, do not keep it.
4. If a run becomes fast but starts printing `PETSc KSP diverged`, that result is invalid even if the wall time looks good.

---

## 12. Recommended Production Settings

For the current 32-rank Turek benchmark:

| Setting | Value |
|---------|-------|
| MPI ranks | `32` |
| Runner | `build/run.sh` |
| Fluid Krylov | `KSPFGMRES` |
| Stabilized fluid `NS_` PC | lower Schur field split; child PCs use `PCHYPRE` + `BoomerAMG` |
| FEM phase-field PC | `PCHYPRE` + `BoomerAMG` |
| Solid Krylov | `KSPFGMRES` |
| Solid PC | `PCHYPRE` + `BoomerAMG` |
| FEM fluid rebuild freq | `20` |
| Implicit-solid rebuild freq | `1` |

If future work targets more speed, the next high-value direction is likely a cleaner blocked PETSc assembly path, but it must preserve the current `NS_iter` behavior and 32-rank stability.
