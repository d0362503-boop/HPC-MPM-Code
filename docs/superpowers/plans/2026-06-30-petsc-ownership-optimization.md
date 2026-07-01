# PETSc Ownership Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Separate natural-node storage, overlap ownership metadata, and PETSc-local numbering so PETSc assembly no longer depends on ambiguous `interior_list` / `l2g_*` semantics.

**Architecture:** Keep all physics arrays (`amat`, `b_rhs`, `x_lhs`, `nmass`, `nforce`, etc.) in the current natural local-node order. Add an explicit ownership metadata layer plus a PETSc-local numbering layer that is owned-first then ghost. Use the PETSc-local layer only inside the PETSc bridge, while keeping global BC IDs and natural-order scatter/recovery explicit and readable.

**Tech Stack:** C++17, MPI, PETSc, existing overlap partition metadata (`aelemmin`, `aelemmax`, `dbc`, `NodeVarComm`)

**File Map:**
- `module/solver/crsmat.h`: declare ownership metadata, PETSc-local numbering buffers, and helper methods
- `module/solver/crsmat.cpp`: build ownership metadata/maps, initialize PETSc mappings, and switch assembly/RHS paths to explicit local/global semantics
- `module/solver/solver.cpp`: update residual-mask insertion to use PETSc-local scalar indices instead of implicit natural-order assumptions
- `module/solver/README.md`: document the three numbering spaces and the new invariants

---

### Task 1: Introduce explicit ownership metadata in `CrsMat`

**Files:**
- Modify: `module/solver/crsmat.h:24-172`
- Modify: `module/solver/crsmat.cpp:267-364`

- [ ] **Step 1: Add explicit ownership and numbering members**

Add the following members to `CrsMat` in place of the ambiguous `l2g_*` semantics. Keep `local_node` and `ghost_node` counts for now to minimize churn in downstream code.

```cpp
    int local_node, ghost_node;

    std::vector<int> owned_natural_ids;
    std::vector<int> ghost_natural_ids;
    std::vector<char> natural_is_owned;
    std::vector<int> natural_to_owned_pos;
    std::vector<int> natural_to_ghost_pos;
    std::vector<int> natural_to_petsc_local;
    std::vector<int> petsc_local_to_natural;

    // Natural-order global IDs: index with natural local node id.
    std::vector<PetscInt> natural_block_gids;
    std::vector<PetscInt> natural_var_gids;

    // PETSc-local global IDs: index with petsc-local position (owned first, ghosts after).
    std::vector<PetscInt> petsc_local_block_gids;
    std::vector<PetscInt> petsc_local_var_gids;

    ISLocalToGlobalMapping petsc_var_lgmap_ = nullptr;
    ISLocalToGlobalMapping petsc_block_lgmap_ = nullptr;
```

- [ ] **Step 2: Add helper declarations that make the numbering model explicit**

Add these declarations to `crsmat.h` near `BuildLGMAP()`:

```cpp
    void BuildNaturalGlobalMaps(int ndof, const std::vector<int> &node_l2g);
    void BuildOwnershipMetadata(const std::vector<int> &node_l2g);
    void BuildPetscLocalMaps(int ndof);
    void CheckOwnershipMetadata() const;

    inline PetscInt NaturalNodeToPetscLocalBlock(int natural_id) const {
        return static_cast<PetscInt>(this->natural_to_petsc_local[natural_id]);
    }

    inline PetscInt NaturalNodeVarToPetscLocalScalar(int natural_id, int var) const {
        return static_cast<PetscInt>(this->natural_to_petsc_local[natural_id] + var * nodec);
    }
```

- [ ] **Step 3: Keep `interior_list` only as a temporary compatibility alias or remove it immediately**

Preferred path: remove `interior_list` and replace its uses with `owned_natural_ids`, because its current name obscures that it is an ownership list, not a geometric subset.

If removing it immediately, the declaration in `crsmat.h` should become:

```cpp
    std::vector<PetscInt> petsc_bc_gids;
```

and every old `interior_list` usage will be updated in later tasks.

- [ ] **Step 4: Rebuild to catch declaration mistakes before implementation**

Run:

```bash
cmake --build build -j8
```

Expected:
- compile fails only because helper definitions and downstream replacements are not implemented yet
- no unrelated files are affected

- [ ] **Step 5: Commit**

```bash
git add module/solver/crsmat.h
git commit -m "refactor: declare explicit PETSc ownership metadata"
```

### Task 2: Split `BuildLGMAP()` into natural, ownership, and PETSc-local phases

**Files:**
- Modify: `module/solver/crsmat.cpp:267-364`

- [ ] **Step 1: Move natural-order global ID construction into `BuildNaturalGlobalMaps()`**

Implement a helper that only answers: "for this natural local node, what is its global node / scalar ID?"

```cpp
void CrsMat::BuildNaturalGlobalMaps(int ndof, const std::vector<int> &node_l2g) {
    VectorAssign(nodec, this->natural_block_gids);
    VectorAssign(nodec * ndof, this->natural_var_gids);

    for (int n = 0; n < nodec; ++n) {
        const PetscInt global_node = static_cast<PetscInt>(node_l2g[n]);
        this->natural_block_gids[n] = global_node;
        for (int var = 0; var < ndof; ++var) {
            this->natural_var_gids[n + var * nodec] = global_node * ndof + var;
        }
    }
}
```

- [ ] **Step 2: Move owner/ghost classification into `BuildOwnershipMetadata()`**

Reuse the current `aelemmin` tie-break, but fill all ownership metadata explicitly.

```cpp
void CrsMat::BuildOwnershipMetadata(const std::vector<int> &node_l2g) {
    this->owned_natural_ids.clear();
    this->ghost_natural_ids.clear();
    this->natural_is_owned.assign(nodec, 0);
    this->natural_to_owned_pos.assign(nodec, -1);
    this->natural_to_ghost_pos.assign(nodec, -1);

    std::vector<int> all_aelemmin(nprocs * 3);
    std::vector<int> all_aelemmax(nprocs * 3);
    MPI_Allgather(aelemmin.data(), 3, MPI_INT, all_aelemmin.data(), 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(aelemmax.data(), 3, MPI_INT, all_aelemmax.data(), 3, MPI_INT, MPI_COMM_WORLD);

    for (int n = 0; n < nodec; ++n) {
        const int global_node = node_l2g[n];
        const int kg = global_node / (xynodecw[0] * xynodecw[1]);
        const int rem = global_node % (xynodecw[0] * xynodecw[1]);
        const int jg = rem / xynodecw[0];
        const int ig = rem % xynodecw[0];

        bool is_owned = true;
        for (int p = 0; p < nprocs; ++p) {
            if (p == myrank) continue;

            bool in_p = true;
            for (int d = 0; d < 3; ++d) {
                const int p_min = all_aelemmin[p * 3 + d];
                const int p_max = all_aelemmax[p * 3 + d] + idimc[d];
                const int coord = (d == 0) ? ig : (d == 1) ? jg : kg;
                if (coord < p_min || coord > p_max) {
                    in_p = false;
                    break;
                }
            }

            if (in_p) {
                for (int d = 0; d < 3; ++d) {
                    const int my_min = aelemmin[d];
                    const int p_min = all_aelemmin[p * 3 + d];
                    if (p_min < my_min) {
                        is_owned = false;
                        break;
                    } else if (p_min > my_min) {
                        break;
                    }
                }
                if (!is_owned) break;
            }
        }

        if (is_owned) {
            this->natural_is_owned[n] = 1;
            this->natural_to_owned_pos[n] = static_cast<int>(this->owned_natural_ids.size());
            this->owned_natural_ids.push_back(n);
        } else {
            this->natural_to_ghost_pos[n] = static_cast<int>(this->ghost_natural_ids.size());
            this->ghost_natural_ids.push_back(n);
        }
    }

    this->local_node = static_cast<int>(this->owned_natural_ids.size());
    this->ghost_node = static_cast<int>(this->ghost_natural_ids.size());
}
```

- [ ] **Step 3: Build PETSc-local numbering in owned-first order**

This helper defines the only numbering PETSc-local APIs should see.

```cpp
void CrsMat::BuildPetscLocalMaps(int ndof) {
    this->petsc_local_to_natural.clear();
    this->petsc_local_to_natural.reserve(nodec);

    for (int nid : this->owned_natural_ids) this->petsc_local_to_natural.push_back(nid);
    for (int nid : this->ghost_natural_ids) this->petsc_local_to_natural.push_back(nid);

    this->natural_to_petsc_local.assign(nodec, -1);
    for (int pos = 0; pos < nodec; ++pos) {
        this->natural_to_petsc_local[this->petsc_local_to_natural[pos]] = pos;
    }

    VectorAssign(nodec, this->petsc_local_block_gids);
    VectorAssign(nodec * ndof, this->petsc_local_var_gids);
    for (int pos = 0; pos < nodec; ++pos) {
        const int natural_id = this->petsc_local_to_natural[pos];
        const PetscInt global_node = this->natural_block_gids[natural_id];
        this->petsc_local_block_gids[pos] = global_node;
        for (int var = 0; var < ndof; ++var) {
            this->petsc_local_var_gids[pos + var * nodec] = global_node * ndof + var;
        }
    }
}
```

- [ ] **Step 4: Add one invariant check helper**

This is required before touching assembly code.

```cpp
void CrsMat::CheckOwnershipMetadata() const {
    if (static_cast<int>(this->owned_natural_ids.size()) != this->local_node) MPI_Abort(MPI_COMM_WORLD, 1);
    if (static_cast<int>(this->ghost_natural_ids.size()) != this->ghost_node) MPI_Abort(MPI_COMM_WORLD, 1);
    if (static_cast<int>(this->petsc_local_to_natural.size()) != nodec) MPI_Abort(MPI_COMM_WORLD, 1);

    for (int n = 0; n < nodec; ++n) {
        if (this->natural_to_petsc_local[n] < 0) MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int global_nodec = xynodecw[0] * xynodecw[1] * xynodecw[2];
    int total_local_node = 0;
    MPI_Allreduce(&this->local_node, &total_local_node, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (total_local_node != global_nodec) MPI_Abort(MPI_COMM_WORLD, 1);
}
```

- [ ] **Step 5: Rewrite `BuildLGMAP()` as an orchestrator**

After computing `node_l2g`, `BuildLGMAP()` should just call the new helpers:

```cpp
    BuildNaturalGlobalMaps(ndof, node_l2g);
    BuildOwnershipMetadata(node_l2g);
    BuildPetscLocalMaps(ndof);
    CheckOwnershipMetadata();
```

- [ ] **Step 6: Rebuild**

Run:

```bash
cmake --build build -j8
```

Expected:
- helper definitions compile
- downstream failures are limited to old uses of `l2g_*` / `interior_list`

- [ ] **Step 7: Commit**

```bash
git add module/solver/crsmat.cpp
git commit -m "refactor: split natural and PETSc ownership maps"
```

### Task 3: Make PETSc setup use explicit map objects instead of ambiguous arrays

**Files:**
- Modify: `module/solver/crsmat.cpp:370-490`
- Modify: `module/solver/crsmat.h:133-166`

- [ ] **Step 1: Replace the old scalar `lgmap` member with explicit scalar and block mappings**

Use the class members added in Task 1:

```cpp
    ISLocalToGlobalMapping petsc_var_lgmap_ = nullptr;
    ISLocalToGlobalMapping petsc_block_lgmap_ = nullptr;
```

Remove the old field:

```cpp
    ISLocalToGlobalMapping lgmap = nullptr;
```

- [ ] **Step 2: Build class-owned LGMAPs from the PETSc-local owned-first arrays**

In `BuildPetscMat()` replace the temporary map creation with:

```cpp
    ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, ndof, nodec, this->petsc_local_block_gids.data(),
                                 PETSC_COPY_VALUES, &this->petsc_block_lgmap_);
    ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, 1, nodec * ndof, this->petsc_local_var_gids.data(),
                                 PETSC_COPY_VALUES, &this->petsc_var_lgmap_);

    MatSetLocalToGlobalMapping(this->petsc_mat, this->petsc_block_lgmap_, this->petsc_block_lgmap_);
    VecSetLocalToGlobalMapping(this->petsc_b, this->petsc_var_lgmap_);
    VecSetLocalToGlobalMapping(this->petsc_x, this->petsc_var_lgmap_);
```

- [ ] **Step 3: Keep natural-order global IDs for scatter and BCs**

In `InitPetscSolver()`, keep `from_is` as a global fetch list in natural order:

```cpp
    ISCreateGeneral(PETSC_COMM_SELF, nodec * ndof, this->natural_var_gids.data(), PETSC_COPY_VALUES, &from_is);
```

This preserves:
- natural-order `seq_x`
- natural-order `x_lhs[n + var * nodec]`
- no coupling between scatter output and PETSc-local owned-first numbering

- [ ] **Step 4: Update `AddBCComponent()` to use the renamed natural-order global IDs**

In `crsmat.h`:

```cpp
    void AddBCComponent(const BoundaryCondition &bc, int offset) {
        if (bc.ibc == 0) return;
        for (int i = 0; i < bc.ibc; ++i) {
            const int nid = bc.nbc[i];
            this->petsc_bc_gids.push_back(this->natural_var_gids[nid + offset]);
        }
    }
```

- [ ] **Step 5: Rebuild**

Run:

```bash
cmake --build build -j8
```

Expected:
- PETSc objects still build successfully
- no use remains of the old ambiguous `lgmap`, `l2g_block_map`, or `l2g_var_map`

- [ ] **Step 6: Commit**

```bash
git add module/solver/crsmat.h module/solver/crsmat.cpp
git commit -m "refactor: bind PETSc setup to explicit ownership maps"
```

### Task 4: Convert PETSc assembly and vector insertion to PETSc-local APIs

**Files:**
- Modify: `module/solver/crsmat.cpp:493-655`
- Modify: `module/solver/solver.cpp:444-470`

- [ ] **Step 1: Switch matrix assembly to `MatSetValuesBlockedLocal()`**

Replace natural-order global insertion with PETSc-local block indices:

```cpp
        const PetscInt block_row = this->NaturalNodeToPetscLocalBlock(natural_row);

        size_t block_col_idx = 0;
        for (int j = row_start; j < row_end; ++j) {
            const int natural_col = this->matcolid[j];
            this->petsc_block_cols_buf[block_col_idx++] = this->NaturalNodeToPetscLocalBlock(natural_col);
        }

        MatSetValuesBlockedLocal(this->petsc_mat, 1, &block_row, ncols,
                                 this->petsc_block_cols_buf.data(),
                                 this->petsc_block_vals_buf.data(), ADD_VALUES);
```

- [ ] **Step 2: Switch identity fill to PETSc-local block rows**

The second pass should use the same local block numbering:

```cpp
            const PetscInt block_row = this->NaturalNodeToPetscLocalBlock(natural_id);
            MatSetValuesBlockedLocal(this->petsc_mat, 1, &block_row, 1, &block_row,
                                     identity_block.data(), ADD_VALUES);
```

- [ ] **Step 3: Switch RHS insertion to `VecSetValuesLocal()`**

Build scalar local indices from the owned natural nodes:

```cpp
    size_t idx = 0;
    for (int i = 0; i < this->local_node; ++i) {
        const int natural_id = this->owned_natural_ids[i];
        for (int var = 0; var < ndof; ++var) {
            this->petsc_indices_buf[idx] = this->NaturalNodeVarToPetscLocalScalar(natural_id, var);
            this->petsc_values_buf[idx] = this->b_rhs[natural_id + var * nodec];
            ++idx;
        }
    }

    VecSetValuesLocal(this->petsc_b, static_cast<PetscInt>(buf_size),
                      this->petsc_indices_buf.data(),
                      this->petsc_values_buf.data(), INSERT_VALUES);
```

- [ ] **Step 4: Switch initial-guess insertion to `VecSetValuesLocal()`**

Use the same owned natural loop for `petsc_x`:

```cpp
    size_t idx = 0;
    for (int i = 0; i < this->local_node; ++i) {
        const int natural_id = this->owned_natural_ids[i];
        for (int var = 0; var < ndof; ++var) {
            this->petsc_indices_buf[idx] = this->NaturalNodeVarToPetscLocalScalar(natural_id, var);
            this->petsc_values_buf[idx] = this->x_lhs[natural_id + var * nodec];
            ++idx;
        }
    }

    VecSetValuesLocal(this->petsc_x, static_cast<PetscInt>(buf_size),
                      this->petsc_indices_buf.data(),
                      this->petsc_values_buf.data(), INSERT_VALUES);
```

- [ ] **Step 5: Update PETSc residual masking in `solver.cpp`**

The active-mask assembly should stop assuming global IDs come from `l2g_var_map`:

```cpp
    size_t idx = 0;
    for (int i = 0; i < this->local_node; ++i) {
        const int natural_id = this->owned_natural_ids[i];
        const PetscScalar is_active = (this->FEM_flag || this->active_row_mask[natural_id] != 0) ? 1.0 : 0.0;
        for (int var = 0; var < this->ndof; ++var) {
            indices[idx] = this->NaturalNodeVarToPetscLocalScalar(natural_id, var);
            values[idx] = is_active;
            ++idx;
        }
    }

    VecSetValuesLocal(active_mask, local_size, indices.data(), values.data(), INSERT_VALUES);
```

- [ ] **Step 6: Keep `MatZeroRowsColumns()` on global BC IDs in this refactor**

Do not switch BC application to the local API yet. The current global path is clearer and BCs are conceptually global constraints:

```cpp
    MatZeroRowsColumns(this->petsc_mat, static_cast<PetscInt>(this->petsc_bc_gids.size()),
                       this->petsc_bc_gids.data(), 1.0e0, this->petsc_x, this->petsc_b);
```

- [ ] **Step 7: Rebuild**

Run:

```bash
cmake --build build -j8
```

Expected:
- successful compile
- no unresolved references to removed map names
- no PETSc API type mismatch on local/global variants

- [ ] **Step 8: Commit**

```bash
git add module/solver/crsmat.cpp module/solver/solver.cpp
git commit -m "refactor: use PETSc-local numbering for assembly and vectors"
```

### Task 5: Verify ownership behavior and update docs

**Files:**
- Modify: `module/solver/README.md:115-158`
- Runtime: existing PETSc solid/fluid/FSI cases

- [ ] **Step 1: Rewrite the ownership section in the solver README**

Replace the current open-concern paragraph with an explicit three-space model:

```md
There are now three distinct numbering spaces:

1. natural local node order
   Used by physics arrays such as `amat`, `b_rhs`, `x_lhs`, `nmass`, `nforce`
2. natural-order global IDs
   Used for global BC row/column IDs and local solution recovery scatter
3. PETSc-local numbering
   Owned nodes first, ghost nodes after; used only by PETSc local assembly/insertion APIs

`owned_natural_ids` and `ghost_natural_ids` define ownership explicitly. `natural_to_petsc_local` is the only valid translation from physics-local indexing into PETSc-local insertion indices.
```

- [ ] **Step 2: Add a warning against reintroducing the old ambiguity**

```md
Do not reuse natural local node ids as if they were PETSc-local insertion indices. Do not assume `local_node` implies that the first `local_node` natural ids are owned. PETSc-local numbering must be constructed explicitly from `owned_natural_ids` and `ghost_natural_ids`.
```

- [ ] **Step 3: Verify with a compile-only smoke test**

Run:

```bash
cmake --build build -j8
```

Expected:
- build succeeds cleanly

- [ ] **Step 4: Verify with the two ownership-sensitive partition layouts**

Run the known cases:

1. PETSc solid/FSI case with the problematic `4-1-4` or equivalent 2D/3D decomposition
2. PETSc solid/FSI case with the known-good `1-1-8` decomposition

Expected:
- no partition-wall artifact
- no ownership-dependent corruption of rows/BCs/scatter
- PETSc and native solutions are qualitatively closer than before

- [ ] **Step 5: Verify no regression in the fluid path**

Run one fluid PETSc case that already works, such as Turek CFD/FSI.

Expected:
- no regression in convergence
- no broken scatter output
- no malformed BC rows

- [ ] **Step 6: Commit**

```bash
git add module/solver/README.md
git commit -m "docs: describe explicit PETSc ownership numbering"
```

---

## Self-Review

**Spec coverage:** This plan covers the ownership optimization described in the conversation: explicit metadata, owned-first PETSc-local numbering, clear separation between natural/global/local maps, and a staged migration to local PETSc APIs. It intentionally does not change BC semantics or inactive-row logic beyond the numbering required to support them safely.

**Placeholder scan:** No `TODO`, `TBD`, or vague “fix later” placeholders remain. Each task names concrete files, code shapes, and verification commands.

**Type consistency:** The plan consistently distinguishes `natural_*`, `owned_*`, `ghost_*`, and `petsc_local_*`. PETSc local APIs are only paired with PETSc-local indices; global BC APIs remain on global IDs.

---

Plan complete and saved to `docs/superpowers/plans/2026-06-30-petsc-ownership-optimization.md`. Two execution options:

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

Which approach?
