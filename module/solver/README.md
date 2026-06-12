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

During assembly:
- inactive nodes are skipped in the first pass
- locally-owned inactive nodes receive an identity block in a second pass to keep the matrix well-conditioned

This replaces the older `nmass < mtol` heuristic with a matrix-content-based check that is consistent with the actual assembled system. See Section 3.5 for details.

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
   - `SolveWithPetsc(ndof, nr_it)`

### 3.1 Ownership

Shared control points are owned by a tie-break based on `aelemmin`.

`interior_list` stores locally owned control points.

Ghost control points are still assembled, because remote element contributions may live there before PETSc redistributes them during `MatAssemblyEnd`.

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

`AssemblePetscMat()` delegates to `AssemblePetscMatBlocked()`, which treats each natural node as one block row and each CSR neighbour as one dense `ndof×ndof` block. The packing order is row-major:

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

At the start of each assembly, `BuildActiveRowMask()` scans every control point's CSR row and marks the node inactive if the row's absolute entry sum falls below `mtol`. FEM nodes are always active. These inactive nodes are skipped in the first assembly pass to avoid adding zero rows. For locally-owned inactive nodes, a second pass inserts an identity block:

```cpp
MatSetValuesBlocked(..., 1, &block_row, 1, &block_row, identity_block.data(), ADD_VALUES);
```

This keeps the matrix non-singular and prevents the linear solver from diverging on empty-space degrees of freedom.

#### Residual consistency for inactive MPM nodes

For the PETSc path, the Newton residual monitor must follow the **same active/inactive rule**
as the matrix assembly.

Current implementation (`ComputePetscResidualStats`):

- `AssemblePetscMatBlocked()` calls `BuildActiveRowMask()` at the start of each assembly,
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

This section describes the current production configuration used by the successful 32-rank Turek FSI tests.

### 5.1 Fluid solver

Fluid uses:

- outer Krylov: `KSPBCGS`
- preconditioner: `PCFIELDSPLIT`
- split type: multiplicative
- velocity block: BoomerAMG
- pressure block: scalar BoomerAMG

Relevant code path:
- `use_fieldsplit = true`
- `ndof = 4`

PETSc options set in code:

```cpp
-pc_fieldsplit_type multiplicative
-pc_fieldsplit_block_size 4
-pc_fieldsplit_0_fields 0,1,2
-pc_fieldsplit_1_fields 3
-fieldsplit_0_ksp_type preonly
-fieldsplit_0_pc_type hypre
-fieldsplit_0_pc_hypre_type boomeramg
-fieldsplit_1_ksp_type preonly
-fieldsplit_1_pc_type hypre
-fieldsplit_1_pc_hypre_type boomeramg
```

Velocity AMG is intentionally lighter-weight than default:

```cpp
-fieldsplit_0_pc_hypre_boomeramg_coarsen_type HMIS
-fieldsplit_0_pc_hypre_boomeramg_interp_type ext+i
-fieldsplit_0_pc_hypre_boomeramg_agg_nl 2
-fieldsplit_0_pc_hypre_boomeramg_P_max 2
-fieldsplit_0_pc_hypre_boomeramg_strong_threshold 0.8
-fieldsplit_0_pc_hypre_boomeramg_max_levels 8
-fieldsplit_0_pc_hypre_boomeramg_nodal_coarsen 6
```

### 5.2 Solid solver

Solid stays monolithic:

- outer Krylov: `KSPBCGS`
- preconditioner: `PCHYPRE`
- AMG type: `BoomerAMG`

Relevant code path:
- `FEM_flag = false`
- `ndof = 3`
- `use_fieldsplit = false`

Current solid AMG tuning:

```cpp
-pc_hypre_boomeramg_coarsen_type HMIS
-pc_hypre_boomeramg_interp_type ext+i
-pc_hypre_boomeramg_agg_nl 2
-pc_hypre_boomeramg_P_max 2
-pc_hypre_boomeramg_strong_threshold 0.8
-pc_hypre_boomeramg_max_levels 8
-pc_hypre_boomeramg_nodal_coarsen 6
```

This keeps solid on AMG, but avoids the older cross-solver AMG-release logic.

### 5.3 Why fluid and solid are different

Fluid and solid now keep their AMG hierarchies independently.

There is no longer any global "release other AMG" behavior. Fluid and solid AMG are allowed to coexist.

---

## 6. AMG Rebuild Policy

Each `CrsMat` controls its own rebuild cadence through:

```cpp
amg_rebuild_freq
```

Current production values:

| System | Value |
|--------|-------|
| Fluid `NS_` | `20` |
| Solid `SM_` | `20` |

Runtime logic in `SolveWithPetsc()`:

- periodic rebuild when `istep % amg_rebuild_freq == 0`
- optional forced rebuild when iterations deteriorate strongly
- rebuild is local to the current solver only
- for implicit solid MPM, later Newton iterations may rebuild again when the previous KSP
  iteration count is already moderate (`prev_ksp_its_ >= 5`)

Implementation detail:

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

For solid, there is one extra Newton-specific rule in `SolveWithPetsc()`:

```cpp
if (nr_it > 0 && !force_rebuild) {
    need_rebuild = false;
    if (!this->FEM_flag && this->prev_ksp_its_ >= 5) { need_rebuild = true; }
}
```

This is intentionally more aggressive than the older "reuse until things are obviously bad"
policy. In implicit solid MPM the tangent matrix can change between Newton iterations even
when the Krylov solve still converges, so allowing an earlier BoomerAMG rebuild is safer
than repeatedly reusing an aging hierarchy.

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

### 9.1 Pure fluid AMG baseline

When fluid is temporarily switched back to monolithic AMG, the Turek 10-step baseline gives:

- `NS_iter ~= 8`

This is the reference used to judge whether the FieldSplit preconditioner is acceptable.

### 9.2 Current recommended production run

Current recommended configuration:

- fluid: `FieldSplit + velocity AMG + pressure AMG + BCGS`
- solid: `monolithic AMG + BCGS`
- local PETSc solution scatter
- solid temporary arrays released every step

Successful 32-rank 10-step run:

- file: `build/stdout/direct_speed_try6.txt`
- total time: about `180.30 s`
- fluid `NS_iter`: `7~9`
- solid linear iterations reported in `NR_converge`: about `1~4`

This is currently the best stable configuration reached in this branch.

### 9.3 Default Krylov Policy

Current default policy in code is simple:

- all PETSc KSP objects default to `KSPBCGS`
- fluid still changes the **preconditioner structure** through `FieldSplit`
- solid still changes the **preconditioner type** through monolithic BoomerAMG

So the distinction between fluid and solid is now mainly in the preconditioner, not in the Krylov family.

---

## 10. Known Good / Known Bad Directions

### 10.1 Good directions

- keep `MatSetBlockSize(ndof)`
- keep fluid velocity/pressure split physical
- keep fluid and solid AMG independent
- keep solid temporary arrays released after solve
- use local solution scatter, not global `CreateToAll`

### 10.2 Bad or abandoned directions

- restoring global cross-solver AMG release logic
- gathering the full global PETSc solution on every rank
- using fluid pressure Jacobi for the current Turek case
- treating the current solid matrix as safely CG-compatible without further proof
- editing headers without `make clean` afterward

---

## 11. Practical Notes for Future Changes

1. If you change any PETSc-related header path or solver header, do `make clean && make -j8`.
2. If you change `nxyr`, you must regenerate partitions before solver benchmarking.
3. If a new preconditioner looks faster but `NS_iter` jumps far above the pure AMG baseline, do not keep it.
4. If a run becomes fast but starts printing `PETSc KSP diverged`, that result is invalid even if the wall time looks good.

---

## 12. Recommended Production Settings

For the current 32-rank Turek benchmark:

| Setting | Value |
|---------|-------|
| MPI ranks | `32` |
| Runner | `build/run.sh` |
| Fluid Krylov | `KSPBCGS` |
| Fluid PC | `PCFIELDSPLIT` |
| Fluid split | velocity AMG + pressure AMG |
| Solid Krylov | `KSPBCGS` |
| Solid PC | `BoomerAMG` |
| Fluid rebuild freq | `20` |
| Solid rebuild freq | `20` |

If future work targets more speed, the next high-value direction is likely a cleaner blocked PETSc assembly path, but it must preserve the current `NS_iter` behavior and 32-rank stability.
