# AGENTS.md — MPM-Code

> Zero-prior-knowledge reference for AI coding agents.

## 1. Project Overview

C++17 MPM/FEM solver for solid/fluid mechanics and FSI. MPI-parallelized with PETSc backend.
Three modes: FLUID (MPM/FEM), SOLID (explicit/implicit), FSI (partitioned coupling).
Key deps: MPI, PETSc 3.24.5 (auto-bootstrapped from `Ext/petsc-3.24.5.tar.gz` into `external/petsc`), GCC ≥11.

## 2. Build & Run

**Current:** Root CMake workflow. Solver source selection is done by uncommenting the relevant `add_subdirectory(...)` line in `work/CMakeLists.txt` and matching the call in `MPM_main.cpp`. Build: `cmake -S . -B build && cmake --build build -j8` → `build/MPM`.

**Never delete `build/`:** The `build/` directory contains generated runtime files, partitioned input data, job outputs, and local working state. Deleting it can destroy an ongoing or future simulation setup. If a clean configure is needed, remove only `CMakeCache.txt` and `CMakeFiles/` inside `build/` — never `rm -rf build`, and never use an alternative build directory.
Run: `mpirun -np N ./build/MPM` (or use `build/run.sh`). Inputs: orchestration file (`file.dat`), parameter file (`input.txt`), grid data (`griddata*.txt`), point data (`pointdata*.txt` / `wpdata*.txt` / `spdata*.txt`).
A convenience script `build/run.sh` exists but is gitignored. It runs `MPM` in the background via `nohup`, supports hyper-threading (`--use-hwthread-cpus --bind-to hwthread`), and accepts a process count via argument or `NP` env var. Run it from `build/`: `sh run.sh [N]`.

### CMake target dependency rules

- Link `mpm_cxx_compat` (directly or via `mpm_modules`) for any target that needs C++17 / `std::filesystem`. GCC < 9 gets `stdc++fs` automatically from this single interface target.
- `mpm_modules` is an `OBJECT` library; its PUBLIC dependencies propagate include directories and compile definitions, but final executables must still link `mpm_modules` to pull in the object files.
- Data tools (`makinput_*`, `makdivide_*`) should be created with `mpm_add_data_tool()` / `mpm_add_divide_tool()` from `cmake/MPMTools.cmake` so link rules stay consistent.
- The configuration stage prints an `MPM Configuration Summary` showing compiler, build type, active solver, MPI, PETSc, HDF5, and ZLIB.

## 3. Code Architecture

- **Global inline variables** in `dataset.h`, `mesh.h`, `mpi_data.h`. Functions operate on global state. Moving variables into classes risks ODR violations.
- **Class hierarchy:** `MaterialPoint` → `SolidMaterialPointBase` → `ImplicitSolidMPM` (implicitmpm).
- **Interpolation:** FLIP / PIC / TPIC / APIC in `map_and_interpolate.h` / `map_and_interpolate.cpp`, selected by `solswitch` (`MapScheme` enum parsed from the input string).
- **Linear algebra:** `CrsMat` in `module/solver/crsmat.h`. Wraps PETSc. Parallel ownership and inactive-row handling are now explicit; see `module/solver/README.md`.
- **Dynamic load balancing:** `module/DLB/` repartitions background elements across ranks from particle samples (fluid MPM only, `do_dlb` flag in `input.txt`). Regions are broadcast state — read them via `mpm_dlb::CurrentRegions()`; see `module/DLB/README.md`.

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

- **Doxygen documentation comments (`/** ... */` or `///`) belong in header files only.** Implementation files should not duplicate interface documentation.

### Comment rules

- **Variable comments: no Doxygen; keep to 3–5 words.** Use plain inline or trailing comments to explain purpose briefly, e.g. `int nstep;  // current time step`.
- **Function comments: always Doxygen in headers.** Document every public/protected function with `/** ... */` or `///`.
  - If the function takes parameters, explain the **physical meaning** of each parameter, not just its type.
  - If the function returns a value, explain the **physical meaning** of the returned value.

### Commit messages

- **English only.** No Chinese or other non-English text in commit messages.

## 5. Critical Parallel Rules

**`node` vs `nodec`:** Computations and PETSc must use `nodec` (control points), never `node` (visualization nodes). Mixing them causes `MatAssemblyEnd` crashes.

**`NodeVarComm`:** Ghost sync uses `+=` on received values, then multiply by `dbc` (pre-computed `1.0/dbn`) to average. Do not replace with naive `MPI_Allreduce`.

**Overlap ownership:** `dbc` encodes overlap weight (1/shared-rank-count), not boolean mask. Tie-break: rank with smallest `aelemmin` owns shared control points.

**Partition consistency:** `nprocs` must equal the partition topology stored in `griddata` (`nxyr` product), and the global mesh parameters in `input.txt` (`xyelem`, `xymax`) must match what the partitioned data was generated with. Mismatches do not fail cleanly: wrong `nprocs` leads to MPI invalid-rank aborts or segfaults, and a stale mesh size silently produces wrong physics.

## 6. Common Pitfalls

1. grep `inline` in `*.h` to find global variables before refactoring.
2. Do not assume class encapsulation.
3. When adding `.cpp` files, update the relevant `CMakeLists.txt`.

## 7. Agent Discipline

- **Think first:** State assumptions explicitly. If uncertain, ask. Surface tradeoffs, don't pick silently.
- **Simplicity first:** Minimum code that solves the problem. No speculative abstractions. No unrequested flexibility.
- **Surgical changes:** Touch only what you must. Match existing style. Clean up only what your changes made unused.
- **Goal-driven:** Transform tasks into verifiable goals. State brief plan with verification steps for multi-step tasks.
