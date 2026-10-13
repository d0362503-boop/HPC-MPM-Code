# `module/solver/` - Linear Solver and PETSc Wrapper

> Updated 2026-10-13: shared matrix interfaces, MPM equation compression,
> physical preconditioner callbacks and solution recovery.

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
| `crsmat_reduced.cpp` | Internal MPM equation compression, field numbering, and full-vector expansion |
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
amat[j + block_id[b]]  // component pair: block_row[b], block_col[b]
```

Fluid (`ndof=4`) stores 16 scalar blocks; solid (`ndof=3`) stores 9.
Monolithic FSI has `ndof=10` and stores 37 selected scalar blocks.
For matrices storing only selected component pairs, `block_row` and `block_col`
define the layout. Block indices cannot generally be inferred as
`row_var * ndof + col_var`.

### 2.5 Active vs inactive nodes (MPM)

In the generic MPM PETSc assembly path, `BuildActiveRowMask()` marks a control
point active when its owning material's nodal mass satisfies

```cpp
owner_->nmass[nid] > mtol
```

`mtol` is an absolute nodal-mass threshold in the case's mass units. A mass equal
to or below it produces an inactive local indicator. FEM control points remain
active without a mass test.

The function first evaluates the mass test on all `nodec` control points,
including ghosts. It then calls `NodeVarComm(local_row_active, 0)` to add the
integer indicators across overlaps and sets `active_row_mask[nid]` from whether
the synchronized sum is positive. A shared control point is therefore active
if any overlapping rank marks it active. This is an OR decision; the indicator
is not multiplied by `dbc` or averaged.

Assembly and elimination use this synchronized mask, preventing a shared row
from being assembled as active on one rank and eliminated as inactive on its
owner. Nodal mass is also synchronized by the phase's P2G transfer.

Inactive block rows are skipped during assembly. After RHS and initial-guess
assembly, `SolveWithPetsc()` eliminates inactive non-Dirichlet rows and columns
with a unit diagonal and zero RHS. Physical boundary constraints are handled
separately. See Section 3.4 for the elimination and residual-mask details.

This is the single-phase mass rule. Monolithic FSI overrides assembly and uses
its per-component `fixed_dof` mask, relative phase-mass thresholds and a
wet-interface test; see the [FSI guide](../../work/src_fsi/MPM_MPM/README.md).

### 2.6 One dimension member

Physics constructors set `CrsMat::ndof`: solid `SM_=3`, fluid `NS_=4`,
FEM phase field `PF_=1`, coupled `fsi_sys=10`. Helpers read this member; no
separate FSI `kDofCount` or repeated dimension argument is needed.
`num_block` remains distinct because FSI stores 37 rather than 100 blocks.

| Operation | Current signature |
|-----------|-------------------|
| Diagonal scaling | `ExtractDiagonal()`, `ComputeDiagonalInverseSqrt()`, `BuildDiagonalPreconditioner()` |
| Numbering | `BuildLGMAP()`, `BuildNaturalGlobalMaps(node_l2g)`, `BuildPetscLocalMaps()` |
| Construction | `BuildPetscMat()`, `InitPetscSolver()` |
| Assembly/RHS | `AssemblePetscMat()`, `UpdatePetscRhs()` |
| Solve | `SolveSystem(NR_it)`, `SolveWithPetsc(NR_it)` |

The physics override is `MaterialPoint::AssemblePetscMat(CrsMat &mat)`;
its component count comes from `mat.ndof`.

---

## 3. `CrsMat` Workflow

The physics constructor sets `owner_`, `ndof` and the block layout.
`BuildCrsMat(num_block)` resets old PETSc state and calls `InitPetscSolver()`
to create maps, the full matrix/vectors and the KSP.

For each assembled correction:

1. `SolveSystem(NR_it)` calls `AssemblePetscMat()`, delegated to the owner.
2. `SolveWithPetsc(NR_it)` collects current BC IDs and copies the full native RHS.
3. Eliminate inactive and physical boundary rows/columns with prescribed values.
4. For MPM, discard explicit zeros in a solve copy, compress, solve from zero
   and expand into the full PETSc vector. FEM solves the full matrix.
5. `RecoverSolution()` fills native `x_lhs` through the owned/ghost scatter.

`ResetPetscSolver()` releases compact maps/matrices as well as full PETSc
objects. `BuildCrsMat()` calls it before a DLB ownership rebuild; see the
[DLB guide](../DLB/README.md).

### 3.1 Ownership

Shared control points are owned by a tie-break based on `aelemmin`.

`owned_natural_ids` stores locally owned control points; `ghost_natural_ids` stores ghosts.

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

The full matrix retains `MatSetBlockSize(mat, this->ndof)` metadata.
The compact MPM operator may not contain complete nodal blocks; its field
groups use explicit index sets rather than assuming the original block stride.

`MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE` is **required for MPM**. See Section 3.4.

### 3.3 Blocked assembly path

`AssemblePetscMat()` treats each natural node as one block row and each CSR neighbour as one dense `ndof×ndof` block. The packing order is row-major:

```
(0,0) (0,1) ... (0,ndof-1)
(1,0) (1,1) ... (1,ndof-1)
...
```

Reusable column and block-value buffers are sized in `BuildPetscMat()` from
the maximum row stencil and `this->ndof`. FSI overrides insertion for its 37
selected scalar blocks.

### 3.4 MPM-specific matrix behavior

#### Dynamic nonzeros from moving particles

`BuildCrsMat()` builds a **static** CSR graph based on the background control-point stencil (`±idimc` in each direction). This graph is fixed for the entire simulation.

Particle motion changes which entries within this stencil receive contributions.
PETSc preallocation can still require additional storage when assembled entries
exceed its estimates; this is distinct from enlarging the local CSR stencil.

By default, PETSc treats new nonzeros as a hard error:

```
PETSC ERROR: New nonzero at (12540,510) caused a malloc.
```

The fix is:

```cpp
MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
```

This permits PETSc to allocate additional storage instead of aborting when an
insertion exceeds preallocation. It does not create missing entries in `matcolid`
or make an out-of-stencil `FindIndex()` lookup valid.

The CSR graph is rebuilt when DLB changes the mesh layout. A change to interpolation
support must likewise be reflected in the graph, independently of PETSc's allocation option.

#### Inactive-node elimination

At the start of each assembly, `BuildActiveRowMask()` computes the mass-based active indicator (Section 2.5) and synchronizes it across overlap nodes via `NodeVarComm` so that a shared node is active if any overlapping rank marks it active. FEM nodes are always active. Nodes that are still inactive after synchronization are skipped in the first assembly pass.

After the RHS and initial guess are assembled, PETSc applies

```cpp
MatZeroRowsColumns(..., inactive_gids, 1.0, petsc_x, petsc_b);
```

to every inactive non-Dirichlet DOF. Their initial guess is zero, so PETSc sets their RHS to zero and removes their row and column from every active equation. Physical Dirichlet DOFs are excluded from this list and are handled by the normal BC pass afterwards. This avoids weak empty-space DOFs feeding back into active MPM equations.

#### Residual consistency for inactive MPM nodes

`ComputeResidualNormSq()` reads the freshly assembled native `b_rhs` for both
PETSc and native solves; it does not evaluate the linear residual `b-A*x`.
It zeroes constrained RHS entries through `BCResidualSet()` and applies the same
callback to an all-ones vector to identify free DOFs. Fixed DOFs are excluded
from both the squared residual sum and the DOF count, independently of residual values.

The native layout is component-blocked: entry `n` belongs to node `n % nodec`.
MPM includes free entries whose nodal mass exceeds `mtol`; FEM includes all free entries.
Over these entries, MPI sums `rtr = sum(b_rhs[n]^2 * dbc[n])` and
`active_dof = sum(dbc[n])`. The outputs are `sqrt(rtr)` and
`sqrt(rtr / active_dof)`, or both zero if no DOFs are counted.
Despite their parameter names (`ref_tol`, `abs_tol`), these are measured norms,
not tolerances. The RMS denominator counts only free active DOFs.

PETSc matrix assembly separately uses `BuildActiveRowMask()`, which synchronizes
the mass-based indicator across shared nodes. The RHS monitor tests local nodal
mass directly and does not consume that mask; consistency relies on synchronized
nodal masses. It does not depend on a previous PETSc assembly or solve.

---

### 3.5 Transparent MPM equation compression

All PETSc MPM paths (fluid, implicit solid, and FSI) use the same internal
compression in `SolveWithPetsc()`. External assembly retains the existing CSR
blocks, `nodec * ndof` RHS/correction vectors, and component-major indexing.

After physical and inactive boundary elimination, `BuildReducedDofMap()` removes
only rows containing exactly the unit diagonal with an exactly zero RHS.
A zero RHS by itself does not mark an equation inactive. Nonzero prescribed
boundary values are retained, and their column contributions have already been
transferred to the RHS.

The retained index set refers to full global equations. Its local entries follow
PETSc row ownership; an `MPI_Allreduce` makes any layout change rebuild the map
on every rank. `SolveReducedSystem()` constructs the retained operator, solves
from zero, and expands the correction into the full PETSc vector.
`RecoverSolution()` then uses the existing owned/ghost scatter to fill the full
native `x_lhs`.

Physical preconditioners call `SetFieldSplitComponents()` with original
components, optionally within a parent component group. This creates explicit
compact index sets; compact indices must never be classified using
`index % ndof`. Physical owners declare `field_splits` as a tree of split names,
original components and optional children. `ConfigureReducedPreconditioner()`
sets up that hierarchy inside `CrsMat`.
FEM keeps its full-system solve path. `HasReducedLayout()` reports whether the
compact map exists; `RetainedDofIds()` exposes original global equation IDs.
The map, `reduced_active_is_` and `reduced_matrix_` remain private state.

An empty retained system uses `PCNONE` and skips physical field-split setup;
full-vector recovery still fills removed components with zero.

The compact operator is rebuilt from current matrix values each solve. A changed
retained layout also rebuilds the KSP/preconditioner; otherwise the configured
reuse policy applies. Mesh/DLB reset releases compact maps and matrices together
with the full PETSc solver.

This stage reduces the equations seen by Krylov and the preconditioner.
Native CSR allocation, full PETSc preallocation, and assembly communication have
not yet been compressed. Their later optimization belongs inside `CrsMat`,
while preserving the external assembly and vector contract.

### 3.6 Shared preconditioner configuration

`CrsMat::ConfigurePreconditioner(PC)` delegates PC choices and the physical
`field_splits` tree to the owner, then converts its root groups to compact
index sets. Each `FieldSplit` contains a PETSc name, original component IDs
and optional child groups.

`CrsMat::ConfigureReducedPreconditioner()` is a normal member used by every
PETSc MPM solve. It sets compact operators, reads options, handles empty
systems with `PCNONE` and creates the root sub-solvers. The private
`ConfigureNestedFieldSplits()` converts child groups to their enclosing
numbering and sets up nested solvers. For a PCKSP child it uses the enclosing
preconditioning matrix as both inner operators, retaining the assembled-Pmat
strategy without repeated fluid elimination.

FSI declares fluid/[solid+lambda], velocity/pressure and solid/lambda groups.
Pure fluid declares velocity/pressure; pure solid has no field groups and keeps
ASM/ILU. No function callback or extra virtual function is needed.
The existing owner virtual `ConfigurePreconditioner()` supplies physical choices;
generic traversal, PETSc setup, compression and recovery stay in `CrsMat`.

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

`BuildKSPSolver()` initially selects outer `KSPFGMRES` with `rtol=1e-10`,
`atol=1e-15`, `divtol=1e6` and `max_it=1000`, then reads options.
A physical solve may set outer tolerances later, as FSI does.
The owner supplies the PC; `use_schur_fieldsplit=false` does not imply AMG.

| Owner/system | Current PC |
|--------------|------------|
| Base `MaterialPoint`, e.g. FEM phase field | HYPRE/BoomerAMG |
| `StabilizedMPM::NS_` | Velocity/pressure lower Schur, PETSc SELFP; child PREONLY + BoomerAMG |
| `StabilizedFEM::NS_` | Velocity/pressure lower Schur, PETSc SELFP; child PREONLY + BoomerAMG |
| `ImplicitSolidMPM::SM_` | ASM overlap 1; local PREONLY + ILU(1), nonzero shift `1e-3` |
| Coupled FSI | Fluid / [solid + lambda], nested SELFP/PCKSP; see the [FSI guide](../../work/src_fsi/MPM_MPM/README.md) |

### 5.1 Stabilized fluid

Original components are velocity `{0,1,2}` and pressure `{3}`.
Pure MPM fluid uses `SetFieldSplitComponents()` for compact index sets and
FieldSplit block size 1. FEM keeps its full block size 4 and field strides.

Both use lower Schur factorization and PETSc SELFP. Current code leaves the
inverse approximation at PETSc's default scalar diagonal:
`S_p = A_11 - A_10 diag(A_00)^(-1) A_01`, without an additional pressure shift.
The full outer operator remains the assembled physical matrix.

Pure-fluid child options use `-fieldsplit_velocity_*` and
`-fieldsplit_pressure_*`. These owners set child defaults with
`PetscOptionsSetValue()`; FSI has its separate `-fsi_fieldsplit_*` hierarchy.

### 5.2 Implicit solid

The standalone KSP remains FGMRES and sets restart 200.
`ImplicitSolidMPM::ConfigurePreconditioner()` selects ASM overlap 1 and PC prefix
`solid_field_`. Local sub-solves use PREONLY + ILU(1), nonzero shift `1e-3`;
their option prefix is `-solid_field_sub_*`.

PREONLY applies the local PC once; it is distinct from the outer FGMRES.
Standalone solid uses the shared outer tolerances above. FSI's solid block has
its own inner `rtol=0.1` and `max_it=100`; these are different solves and counts.

### 5.3 Independent resources

Each `CrsMat` owns its KSP/PC and numbering. Rebuilding one system does not
release another's preconditioner.

## 6. Preconditioner Rebuild Policy

The legacy name `amg_rebuild_freq` controls PC cadence, including ASM/FieldSplit.

| System | Frequency / iteration-growth policy |
|--------|-------------------------------------|
| Pure MPM fluid `NS_` | Every step; force next rebuild after more than 2x iteration growth |
| Implicit solid `SM_` | Every step; force next rebuild after more than 2x iteration growth |
| FEM fluid `NS_` / `PF_` | Every 20 steps; iteration-growth policy enabled |
| Coupled FSI `fsi_sys` | Every step; iteration-growth rebuilding disabled |

`NeedsPreconditionerRebuild(NR_it)` consumes a pending forced rebuild. Without
one, later Newton iterations (`NR_it>0`) reuse the PC; the first correction
follows the step cadence.

In MPM, changed retained equations on any rank force all ranks to rebuild
the compact layout and KSP/PC. Compact matrix values are recreated every solve
even when the PC is reused. A KSP rebuild restores its type, prefix and tolerances,
calls `ConfigureReducedPreconditioner()` to configure the physical PC, read
options and initialize its declared field hierarchy.

FEM retains the full operator. A scheduled rebuild disables PC reuse;
`PCReset()` and reconfiguration are used only for `use_schur_fieldsplit=false`,
preserving FieldSplit structure otherwise. Mesh/DLB reset releases both full
and compact state.

---

## 7. Solution Recovery and Memory

After a compact solve, retained components are copied into the full PETSc
solution vector and removed components remain zero. Each rank scatters only
its owned/ghost `nodec * ndof` entries into `seq_x`, then fills native `x_lhs`.
The full global solution is not gathered onto every rank.

After standalone solid/fluid solves, temporary `adiag`, `amat`, `b_rhs` and
`x_lhs` arrays are released. PETSc solver resources stay local to their owner
and are reset when their mesh/numbering is rebuilt.

## 8. Implementation Invariants

### Required properties

- retain full-matrix block metadata and explicit compact field index sets
- keep each physical solver's PC resources independent
- keep solid temporary arrays released after solve
- use local solution scatter, not global `CreateToAll`
- for poorly scaled mixed u-p linear systems, check KSP tolerances and scaling; separately verify the nonlinear residual at the updated state

### Avoid

- releasing another physical solver's PC resources
- gathering the full global PETSc solution on every rank
- using a loosely set `rtol` (e.g. `1.0e-8`) for unscaled mixed u-p systems, which can accept a trivial initial guess as "converged"
- treating the current solid matrix as safely CG-compatible without further proof
- classifying compressed equations with `compact_index % ndof`

---

## 9. Build and Validation

1. Reconfigure with `cmake -S . -B build` after adding/removing sources collected by CMake. Rebuild in the existing `build/` directory; normal header edits are handled by build dependencies.
2. If you change `nxyr`, you must regenerate partitions before solver benchmarking.
3. Compare time, inner/outer iterations and physical residuals on the same state when changing the PC.
4. If a run becomes fast but starts printing `PETSc KSP diverged`, that result is invalid even if the wall time looks good.

---

## 10. Current Settings by Solver

| Setting | Pure MPM fluid | Implicit solid | Monolithic FSI |
|---------|----------------|----------------|----------------|
| Components | 4 | 3 | 10, with 37 stored blocks |
| Equation compression | Shared `CrsMat` | Shared `CrsMat` | Shared `CrsMat` |
| Native vector length | `nodec * 4` | `nodec * 3` | `nodec * 10` |
| PC | Velocity/pressure Schur + BoomerAMG | ASM overlap 1 / ILU(1), shift `1e-3` | Fluid / [solid + lambda], nested SELFP/PCKSP |
| Node activation | Absolute `nmass > mtol` | Absolute `nmass > mtol` | Relative phase mean times `1e-4` by default, plus wet multiplier test |
| PC cadence | 1 step | 1 step | 1 step |

The FSI launch configuration omits the absolute-mass override. Existing
processes retain the options with which they were started. See the
[FSI guide](../../work/src_fsi/MPM_MPM/README.md).

## Current Convergence Limitations

`CrsMat::CheckNRConvergence()` must follow physical RHS assembly at the current
iterate, before the next linear solve. Iteration zero stores the initial L2 norm
in `r0r`; subsequent iterations compare the newly assembled norm `rkr` to it.
The relative value is `rkr/r0r` when `r0r > 1e-30`, otherwise zero.
The current acceptance rule is:

```cpp
rtr_ref < 1e-8 || r0r < 1e-10 || rtr_abs < 1e-12
```

The initial-norm branch can accept iteration zero without a linear solve.
Absolute thresholds depend on units, mesh and model scale. Mixed fluid momentum and
continuity entries are combined without field scaling. These limitations matter
when interpreting convergence; a small norm is not by itself an accuracy study.
Failure to meet the criteria at the iteration limit calls `MPI_Abort`.

`SolveWithPetsc()` logs negative KSP convergence reasons on rank zero, then still
scatters the solution and returns the iteration count. It does not currently
return a failure sentinel or automatically retry. Treat divergence logs as failed
linear solves even if the surrounding time-step loop continues.
