# AGENTS.md — MPM-Code

> Zero-prior-knowledge reference for AI coding agents.

## 1. Project Overview

C++17 MPM/FEM solver for solid/fluid mechanics and FSI. MPI-parallelized with PETSc backend.
Three modes: FLUID (MPM/FEM), SOLID (explicit/implicit), FSI (partitioned coupling).
Key deps: MPI, PETSc 3.24.5 (at `/home/pan/petsc-3.24.5`), GCC ≥11.

## 2. Build & Run

**Current:** Root `Makefile` with boolean flags (`USE_SRC_FSI`, `USE_FEM`, etc.). Build: `make -j8` → `build/MPM.exe`.
**Future:** CMake prepared but not active. `data/` and `vtk/` have standalone Makefiles.
Run: `mpiexec -np N ./MPM`. Inputs: `griddata.txt`, `pointdata.txt`, `input.txt`.

`build/file.dat` is read from the `build/` working directory at runtime.
Paths inside `file.dat` should point to files under `build/` such as `./data/generate/...`, `./data/divide/...`, and `./res/...`.
Do not point `file.dat` back at source-tree paths under `../data/...`.

## 3. Code Architecture

- **Global inline variables** in `dataset.h`, `mesh.h`, `mpi_data.h`. Functions operate on global state. Moving variables into classes risks ODR violations.
- **Class hierarchy:** `MaterialPoint` → `SolidMaterialPointBase` → `SolidMaterialPoint` (ImplicitMPM).
- **Interpolation:** FLIP / PIC / TPIC / APIC in `map_and_interpolate.h`, selected by `solswitch` string.
- **Linear algebra:** `CrsMat` in `module/solver/crsmat.h`. Wraps PETSc. Known issue: parallel ownership bug with shared overlap nodes (see `petsc_parallel_debug_report_v2.md`).

## 4. Code Style

Follow `.clang-tidy` (Google style). When editing legacy files, match surrounding style.

| Entity | Style | Example |
|--------|-------|---------|
| Types | `CamelCase` | `MaterialPoint` |
| Functions | `CamelCase` | `BuildMesh()` |
| Local vars | `lower_case` | `node_count` |
| Members | `lower_case_` | `num_` |
| Globals | `g_` prefix | `g_time_step` |
| Constants | `k` prefix | `kDaysInAWeek` |
| Macros | `UPPER_CASE` | `DEBUG_MODE` |

- In class member functions, call other member functions with explicit `this->`.

## 5. Critical Parallel Rules

**`node` vs `nodec`:** Computations and PETSc must use `nodec` (control points), never `node` (visualization nodes). Mixing them causes `MatAssemblyEnd` crashes.

**`NodeVarComm`:** Ghost sync uses `+=` on received values, then multiply by `dbc` (pre-computed `1.0/dbn`) to average. Do not replace with naive `MPI_Allreduce`.

**Overlap ownership:** `dbc` encodes overlap weight (1/shared-rank-count), not boolean mask. Tie-break: rank with smallest `aelemmin` owns shared control points.

## 6. Common Pitfalls

1. grep `inline` in `*.h` to find global variables before refactoring.
2. Do not assume class encapsulation.
3. When adding `.cpp` files, update both `CMakeLists.txt` and `Makefile`.
4. `data/` and `vtk/` tools use old SoA→AoS data structures; check for transposed indices.
5. `data/divide_fsi/` is out of sync with global inline vars (`sp`, `wfem`).

## 7. Agent Discipline

- **Think first:** State assumptions explicitly. If uncertain, ask. Surface tradeoffs, don't pick silently.
- **Simplicity first:** Minimum code that solves the problem. No speculative abstractions. No unrequested flexibility.
- **Surgical changes:** Touch only what you must. Match existing style. Clean up only what your changes made unused.
- **Goal-driven:** Transform tasks into verifiable goals. State brief plan with verification steps for multi-step tasks.
