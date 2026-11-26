# AGENTS.md — MPM-Code

> Zero-prior-knowledge reference for AI coding agents.

## 1. Project Overview

C++17 MPM/FEM solver for solid/fluid mechanics and FSI. MPI-parallelized with PETSc backend.
Three modes: FLUID (MPM/FEM), SOLID (explicit/implicit), FSI (partitioned coupling).
Key deps: MPI, PETSc 3.24.5 (auto-bootstrapped from `Ext/petsc-3.24.5.tar.gz` into `external/petsc`), GCC ≥11.

## 2. Build & Run

**Current:** Root CMake workflow. Every driver under `work/` is built as its own executable via `mpm_add_driver()` (`cmake/MPMTools.cmake`): each executable is written to its driver's build directory, e.g. `build/work/src_fsi/MPM_MPM/MPM_MPM_FSI`, next to that driver's `file.dat`. Build: `cmake -S . -B build && cmake --build build -j8`. Parameter files are copied into `build/data/` only when missing (`mpm_copy_input_once`).

**Never delete `build/`:** The `build/` directory contains generated runtime files, partitioned input data, job outputs, and local working state. Deleting it can destroy an ongoing or future simulation setup. If a clean configure is needed, remove only `CMakeCache.txt` and `CMakeFiles/` inside `build/` — never `rm -rf build`, and never use an alternative build directory.
Run from the driver directory, e.g. `build/work/src_fsi/MPM_MPM`: `mpirun -np N ./MPM_MPM_FSI` (or from `build/`: `sh run.sh N work/src_fsi/MPM_MPM`). Inputs: orchestration file (`file.dat`), parameter file (`input.txt`), grid data (`griddata*.txt`), point data (`pointdata*.txt` / `wpdata*.txt` / `spdata*.txt`).
A convenience script `build/run.sh` exists but is gitignored. It runs the single `MPM_*` executable found in `RUN_DIR` in the background via `nohup`, supports hyper-threading (`--use-hwthread-cpus --bind-to hwthread`), and accepts a process count via argument or `NP` env var. Run it from `build/`: `sh run.sh N [CPU_LIST] [RUN_DIR]`, e.g. `sh run.sh 61 work/src_fsi/MPM_MPM`; logs go to `RUN_DIR/stdout` and `RUN_DIR/stderr`.

### CMake target dependency rules

- Link `mpm_cxx_compat` (directly or via `mpm_modules`) for any target that needs C++17 / `std::filesystem`. GCC < 9 gets `stdc++fs` automatically from this single interface target.
- `mpm_modules` is an `OBJECT` library; its PUBLIC dependencies propagate include directories and compile definitions, but final executables must still link `mpm_modules` to pull in the object files.
- Data tools (`makinput_*`, `makdivide_*`) should be created with `mpm_add_data_tool()` / `mpm_add_divide_tool()` from `cmake/MPMTools.cmake` so link rules stay consistent.
- The configuration stage prints an `MPM Configuration Summary` showing compiler, build type, active solver, MPI, PETSc, HDF5, and ZLIB.

Each driver reads `file.dat` from its working directory, normally its driver directory such as `build/work/src_fsi/MPM_MPM/file.dat`.
Paths inside `file.dat` are relative to that directory and point into `build/`, such as `../../../data/generate/...`, `../../../data/divide/...`, and `../../../res/...`.
Never point `file.dat` at source-tree paths in the repository `data/` directory.

**Checkpoint restarts:** Use the output prefix already set for the selected case in the driver's `file.dat`. Change the prefix or output directory only when the user explicitly requests a different destination.

**User approval:** Obtain explicit user approval before changing any content in `file.dat` or creating any directory.

## 3. Code Architecture

- **Grid and run state are passed explicitly.** `Grid` (`mesh.h`: `Mesh`, `ControlMesh`, `Partition`), `RunParams` and `TimeState` (`dataset.h`) are owned by the driver `main` (data tools own them outright). Physics classes hold `const Grid &grid_` / `const RunParams &params_` set in their constructors; `TimeState` is passed per call (`SolveNS(ts)`, `Particle2Node(ts)`, …); only `ApplyDLB(Grid &)` and `DataInput(Grid &, RunParams &)` write the grid. Free functions take what they use (`MakeSF(grid, …)`, `NodeVarComm(part, …)`).
- **Transitional globals, CrsMat only:** `g_grid` / `g_params` / `g_time` and the few aliases left in `mesh.h` / `dataset.h` (`nodec`, `nuc`…`npc`, `xynodec`, `xynodecw`, `idimc`, `aelemmin/max`, `dbc`, `mtol`, `istep`) plus the two-argument `NodeVarComm` exist only because `CrsMat` still reads them; drivers therefore bind `grid`/`params`/`ts` to these instances. Never use them in new or touched non-`CrsMat` code; remove them when `CrsMat` takes the grid. `myrank` / `nprocs` stay global.
- **Class hierarchy:** `MaterialPoint` → `SolidMaterialPointBase` → `ImplicitSolidMPM` (implicitmpm) / `ExplicitSolidMPM`; `MaterialPoint` → `StabilizedMPM`. Shared implicit pieces live in `module/implicit/` (see its README). Implicit time integration (generalized-α integrator, endpoint kinematics, inertia terms, kinematic commits) lives in the `ImplicitMPM` component (`module/implicit/implicit_mpm.h`), owned as `implicit_` by `ImplicitSolidMPM` and `StabilizedMPM`; `StabilizedFEM` owns its own `integrator_`. The implicit MPM Newton loop is `SolveNewton()` over the `NonlinearProblem` interface (`module/implicit/newton_solver.h`), implemented by `ImplicitSolidMPM`, `StabilizedMPM`, `StabilizedFEM` and `MPMMPMMonolithicFSI`; failure aborts stay in their convergence checks.
- **Interpolation:** FLIP / PIC / TPIC / APIC in `map_and_interpolate.h` / `map_and_interpolate.cpp`, selected by `solswitch` (`MapScheme` enum parsed from the input string).
- **Linear algebra:** `CrsMat` in `module/solver/crsmat.h`. Wraps PETSc. Parallel ownership and inactive-row handling are now explicit; see `module/solver/README.md`.
- **Dynamic load balancing:** `module/DLB/` repartitions background elements across ranks from particle samples (fluid MPM only, `do_dlb` flag in `input.txt`). Regions are broadcast state stored in `Partition::regions`; see `module/DLB/README.md`.

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

1. Do not add `inline` globals; extend `Grid` / `RunParams` / `TimeState` (or the owning class) and pass them explicitly.
2. Do not assume class encapsulation.
3. When adding `.cpp` files, update the relevant `CMakeLists.txt`.

## 7. Agent Discipline

- **Clarification:** Resolve uncertainty from the existing context and focused inspection first. For routine, reversible choices within scope, use a reasonable default and continue; state material assumptions briefly. Ask only when missing information would materially change the requested outcome, correctness, or a consequential tradeoff that cannot reasonably be inferred. While waiting, continue work that does not depend on the answer.
- **Simplicity first:** Minimum code that solves the problem. No speculative abstractions. No unrequested flexibility.
- **Surgical changes:** Touch only what you must. Match existing style. Clean up only what your changes made unused.
- **Goal-driven:** Transform tasks into verifiable goals. State brief plan with verification steps for multi-step tasks.
