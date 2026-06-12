# MPM-Code PETSc Quick Reference

> Full documentation: `module/solver/README.md` (489 lines).

## Current Production Configuration (32-rank Turek FSI)

| Setting | Value |
|---------|-------|
| MPI ranks | `32` |
| Runner | `cd build && sh run.sh` |
| Fluid Krylov | `KSPBCGS` |
| Fluid PC | `PCFIELDSPLIT` (velocity AMG + pressure AMG) |
| Solid Krylov | `KSPBCGS` |
| Solid PC | `BoomerAMG` (monolithic) |
| AMG rebuild freq | `20` steps for both fluid and solid |

## Critical Rules

- **Always use `nodec`, never `node`** for PETSc DOFs, ownership, and local row counts.
- **`MatSetBlockSize(ndof)` is mandatory.** Without it, BoomerAMG quality drops sharply.
- **`MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE` is required for MPM.** Particles move and create new couplings outside the static CSR stencil.
- **Inactive nodes** are detected at solve time by `BuildActiveRowMask()`: if a node's assembled CSR row has absolute sum below `mtol`, it is marked inactive. Skipped in assembly pass 1, then receive identity blocks in pass 2.
- **Residual monitors** (`ComputeRefResidual` / `ComputeAbsResidual`) now use `ComputePetscResidualStats` for the PETSc path, which computes `r = b - A·x` via `MatMult` + `VecAYPX`, zeros inactive DOFs with a PETSc mask vector, and returns the L2 norm via `VecNorm`. This ensures the residual monitor uses the **same PETSc operator** that KSP sees.

## AMG Rebuild Policy

- Rebuild when `istep % 20 == 0`
- Solid only: also rebuild during Newton if previous KSP iterations were moderate
- Fluid and solid AMG hierarchies are **independent** — no cross-solver release logic

## Known Good Directions

- keep `MatSetBlockSize(ndof)`
- keep fluid velocity/pressure split physical
- keep fluid and solid AMG independent
- keep solid temporary arrays (`adiag`, `amat`, `b_rhs`, `x_lhs`) released after solve
- use local solution scatter (per-rank `nodec * ndof`), not global `VecScatterCreateToAll`

## Known Bad / Abandoned Directions

- restoring global cross-solver AMG release logic
- gathering full global PETSc solution on every rank
- using fluid pressure Jacobi for current Turek case
- treating solid matrix as CG-compatible without proof
- editing PETSc headers without `make clean`

## When Things Go Wrong

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| `PETSC ERROR: New nonzero caused a malloc` | Missing `MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE` | Add the `MatSetOption` |
| `NS_iter` jumps far above 9 | Preconditioner degraded | Check AMG rebuild, compare against pure AMG baseline |
| `KSP diverged` after speed tweak | Convergence sacrificed for speed | Revert; diverged = invalid even if wall time looks good |
| OOM on 32-rank runs | Solid temp arrays not released | Ensure `SM_.adiag/amat/b_rhs/x_lhs` freed after solve |
| Matrix size mismatch crash | Used `node` instead of `nodec` | Switch to `nodec` for all solver-side sizing |
