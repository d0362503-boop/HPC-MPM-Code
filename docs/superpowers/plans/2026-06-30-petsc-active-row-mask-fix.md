# PETSc Active Row Mask Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the false partition-boundary "wall" in PETSc implicit solid/MPM solves by making inactive-row detection overlap-consistent before identity fill is applied.

**Architecture:** Keep the current global PETSc assembly path (`MatSetValuesBlocked`, `VecSetValues`, `MatZeroRowsColumns`) unchanged. Fix only the MPM/PETSc inactive-row decision so that a shared control point is marked active if any overlapping rank assembled nonzero content for that row. This preserves the existing identity-fill safeguard for truly empty rows while eliminating owner-rank false negatives.

**Tech Stack:** C++17, MPI, PETSc, existing `NodeVarComm` overlap communication, `CrsMat`

---

### Task 1: Make the active-row mask overlap-consistent

**Files:**
- Modify: `module/solver/crsmat.cpp:111-128`
- Test: `module/solver/crsmat.cpp:468-512` behavior through existing PETSc solve path

- [x] **Step 1: Replace the purely local mask with a local indicator buffer**

Use the existing row-content check, but do not write directly into `this->active_row_mask` during the first pass.

```cpp
void CrsMat::BuildActiveRowMask() {
    this->active_row_mask.assign(nodec, 1);
    if (this->FEM_flag) return;

    std::vector<int> local_row_active(nodec, 0);

    for (int nid = 0; nid < nodec; ++nid) {
        double row_abs_sum = 0.0e0;
        for (int row_var = 0; row_var < this->ndof; ++row_var) {
            for (int j = this->matrow[nid]; j < this->matrow[nid + 1]; ++j) {
                for (int col_var = 0; col_var < this->ndof; ++col_var) {
                    int bid = this->block_id[row_var * this->ndof + col_var];
                    row_abs_sum += std::abs(this->amat[j + bid]);
                }
            }
        }
        if (row_abs_sum >= mtol) local_row_active[nid] = 1;
    }
```

- [x] **Step 2: Synchronize the indicator across overlap control points**

Use existing overlap communication so that any shared node with nonzero content on at least one rank becomes active on every sharing rank.

```cpp
    NodeVarComm(local_row_active, 0);
```

This works because `NodeVarComm` adds neighbor values on shared control points. After communication:
- truly empty shared rows stay `0`
- one-rank-active shared rows become `1` on both sides
- multi-rank-active shared rows become `>= 1`

- [x] **Step 3: Rebuild `active_row_mask` from the synchronized indicator**

Only after the overlap sync should the solver decide whether a row is active or inactive.

```cpp
    for (int nid = 0; nid < nodec; ++nid) {
        this->active_row_mask[nid] = (local_row_active[nid] > 0) ? 1 : 0;
    }

    return;
}
```

- [x] **Step 4: Rebuild and verify the code still compiles**

Run:

```bash
cmake --build build -j8
```

Expected:
- compile succeeds
- no new MPI/PETSc type errors
- no missing include or template-resolution errors

- [ ] **Step 5: Commit** (deferred until validation completes)

```bash
git add module/solver/crsmat.cpp
git commit -m "fix: synchronize PETSc active row mask across overlap nodes"
```

### Task 2: Keep identity fill, but only after the synchronized mask

**Files:**
- Review only: `module/solver/crsmat.cpp:468-512`
- Test: PETSc solid solve with problematic 3D/2D partition

- [x] **Step 1: Do not change the existing identity-fill branch yet**

Keep this logic as-is:

```cpp
if (!this->FEM_flag) {
    for (int ii = 0; ii < this->local_node; ++ii) {
        int natural_id = this->interior_list[ii];
        if (this->active_row_mask[natural_id] != 0) continue;

        PetscInt block_row = this->l2g_block_map[natural_id];
        MatSetValuesBlocked(this->petsc_mat, 1, &block_row, 1, &block_row, identity_block.data(), ADD_VALUES);
    }
}
```

Reason:
- this branch is still needed for truly empty MPM rows
- after Task 1, false inactive rows on partition boundaries should no longer fall into this path

- [x] **Step 2: Verify the expected behavioral change** (diagnostic confirms false inactive rows on shared nodes)

Expected after Task 1:
- shared rows with real contributions are no longer skipped on the owner rank
- owner-only identity fill now applies only to globally empty rows
- PETSc and native solves should move closer instead of diverging at partition seams

- [ ] **Step 3: Commit only if no extra changes were needed** (deferred until validation completes)

```bash
git add module/solver/crsmat.cpp
git commit -m "refactor: preserve identity fill after overlap-synced active mask"
```

### Task 3: Add one temporary diagnostic to prove the root cause

**Files:**
- Modify: `module/solver/crsmat.cpp:111-128`
- Remove after verification: same file
- Test: problematic PETSc run with `4-1-4` partition and known-good `1-1-8` partition

- [x] **Step 1: Add a temporary counter for rows that change from inactive to active after overlap sync**

```cpp
    int locally_inactive_but_globally_active = 0;
    for (int nid = 0; nid < nodec; ++nid) {
        const int was_local_active = local_row_active_before_comm[nid];
        const int is_overlap_active = local_row_active[nid];
        if (was_local_active == 0 && is_overlap_active > 0) {
            locally_inactive_but_globally_active++;
        }
    }
```

You may need to keep a copy before `NodeVarComm`:

```cpp
    std::vector<int> local_row_active_before_comm = local_row_active;
```

- [x] **Step 2: Print the diagnostic once per rank for the first few steps**

Minimal print shape:

```cpp
    if (locally_inactive_but_globally_active > 0) {
        std::cout << "Rank " << myrank
                  << " overlap-synced active rows = "
                  << locally_inactive_but_globally_active << std::endl;
    }
```

Expected:
- problematic `4-1-4` run shows nonzero counts near the bad partition
- `1-1-8` may show fewer or zero such rows

- [x] **Step 3: Comment out the temporary diagnostic after confirmation**

Delete:
- the copied pre-comm vector
- the integer counter
- the `std::cout`

Keep only the synchronized mask logic from Task 1.

- [ ] **Step 4: Commit** (deferred until final validation)

```bash
git add module/solver/crsmat.cpp
git commit -m "chore: remove temporary active-row overlap diagnostics"
```

### Task 4: Verify against the known bad and known good partitions

**Files:**
- No source changes
- Runtime inputs: existing divided data under `build/data/divide/fsi/myrank_data/`

- [x] **Step 1: Rebuild**

Run:

```bash
cmake --build build -j8
```

Expected:
- `build/MPM` updated successfully

- [ ] **Step 2: Run the known-bad PETSc case** (8×1×1 baseline running; 4×1×4 queued)

Run the same PETSc setup that currently shows the vertical wall, using the `4-1-4` decomposition.

Expected:
- partition seam disappears or is dramatically reduced
- solid field no longer looks artificially pinned on the partition line

- [ ] **Step 3: Run the known-good comparison case** (1×1×8 previously failed with OOM; will retry after baseline)

Run the same setup with the `1-1-8` decomposition.

Expected:
- no regression relative to the current known-good behavior
- iteration counts may change slightly, but the solution should remain smooth

- [ ] **Step 4: Compare against native solver**

Set `SM_.use_petsc = false` for the same case and compare:
- gross deformation shape
- partition-seam behavior
- whether PETSc now tracks the native result qualitatively

Expected:
- PETSc result should move toward the native baseline
- if a seam remains after Task 1, the next suspect is `petsc_bc_gids -> MatZeroRowsColumns`, not `l2g_*`

- [ ] **Step 5: Record acceptance criteria**

Accept this fix only if all of the following hold:
- no partition-wall artifact in the former bad case
- no regression in the former good case
- no KSP divergence introduced
- no obvious loss of deformation continuity across overlap regions

### Task 5: Update solver notes so the bug does not come back

**Files:**
- Modify: `module/solver/README.md:176-206`

- [x] **Step 1: Update the inactive-row documentation**

Replace the current description so it explicitly says the PETSc MPM active mask must be overlap-consistent before owned-row identity fill is applied.

Suggested replacement text:

```md
For PETSc MPM solves, the active/inactive decision must be synchronized across overlap control points before assembly skips rows or inserts owned-row identity blocks. A shared control point is considered active if any overlapping rank assembled nontrivial row content for that node. This prevents owner-rank false negatives that can otherwise turn partition boundaries into artificial Dirichlet-like walls.
```

- [x] **Step 2: Add a warning about what not to do**

Add a short note:

```md
Do not classify active rows from purely rank-local `amat` content and immediately apply identity fill on owned rows. In overlapping decompositions this can misclassify shared rows whose physical contributions are split across ranks.
```

- [ ] **Step 3: Rebuild docs-only confidence check**

Run:

```bash
git diff -- module/solver/README.md
```

Expected:
- documentation reflects the new invariant
- no unrelated solver behavior claims were changed

- [ ] **Step 4: Commit**

```bash
git add module/solver/README.md
git commit -m "docs: describe overlap-consistent PETSc active row handling"
```

---

## Self-Review

**Spec coverage:** This plan addresses only the first and strongest suspect: PETSc MPM inactive-row classification before identity fill. It intentionally does not change `l2g_block_map`, `l2g_var_map`, PETSc ownership setup, or BC list generation.

**Additional finding from static scan:** `BuildLGMAP()` constructs the PETSc local-to-global mapping with `nodec` block entries while the matrix local size is `local_node * ndof`. This appears to rely on the partitioner placing owned nodes at the front of the local index range. This ownership/LGMAP consistency issue is now documented in `module/solver/README.md` Section 3.1 as an open concern, but it is out of scope for this plan and should be tracked separately if validation shows remaining seam artifacts.

**Placeholder scan:** No `TODO` or `TBD` placeholders remain. All code-touching steps point to concrete files and include the expected code shape.

**Type consistency:** The plan uses existing names and types already present in the codebase: `active_row_mask`, `NodeVarComm`, `interior_list`, `l2g_block_map`, `FEM_flag`, `MatSetValuesBlocked`.

---

Plan complete and saved to `docs/superpowers/plans/2026-06-30-petsc-active-row-mask-fix.md`. Two execution options:

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

Which approach?
