# Data CMake OOP Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild `data/` as an object-oriented subsystem with one selectable CMake-built generator target at a time and outputs placed under `build/data/`.

**Architecture:** Introduce a shared `DataGenerator` base class and a shared `output_util` layer, then migrate `fluid`, `solid`, and `fsi` onto concrete generator classes. Only after the code layout is stable, wire the targets into `data/CMakeLists.txt` and root CMake selection logic.

**Tech Stack:** C++17, existing `module/` backend, CMake, optional MPI/HDF5 helpers.

---

### Task 1: Restore Abstract Skeleton

**Files:**
- Create: `data/data_generator.h`
- Create: `data/data_generator.cpp`
- Create: `data/fluid/fluid_generator.h`
- Create: `data/fluid/fluid_generator.cpp`
- Create: `data/solid/solid_generator.h`
- Create: `data/solid/solid_generator.cpp`
- Create: `data/fsi/fsi_generator.h`
- Create: `data/fsi/fsi_generator.cpp`

- [ ] Add `DataGenerator` with a shared `Run()` workflow and pure virtual hooks.
- [ ] Add empty-but-compilable concrete generator shells for `fluid`, `solid`, and `fsi`.
- [ ] Do not change `data/CMakeLists.txt` yet.

### Task 2: Rebuild Shared Output Utilities

**Files:**
- Create: `data/output_util.h`
- Create: `data/output_util.cpp`

- [ ] Move repeated mesh-header and VTK/HDF5 helper logic into shared functions.
- [ ] Use `CamelCase` names consistent with project style.

### Task 3: Migrate Fluid

**Files:**
- Modify: `data/fluid/makinput.cpp`
- Modify: `data/fluid/fluid_generator.cpp`

- [ ] Reduce `data/fluid/makinput.cpp` to a thin `main()`.
- [ ] Move the current fluid generation body into `FluidGenerator`.

### Task 4: Migrate Solid

**Files:**
- Modify: `data/solid/makinput.cpp`
- Modify: `data/solid/solid_generator.cpp`

- [ ] Reduce `data/solid/makinput.cpp` to a thin `main()`.
- [ ] Move the current solid generation body into `SolidGenerator`.

### Task 5: Migrate FSI

**Files:**
- Modify: `data/fsi/makinput.cpp`
- Modify: `data/fsi/fsi_generator.cpp`

- [ ] Reduce `data/fsi/makinput.cpp` to a thin `main()`.
- [ ] Move the current FSI generation body into `FsiGenerator`.

### Task 6: CMake Integration

**Files:**
- Modify: `data/CMakeLists.txt`
- Modify: `CMakeLists.txt`
- Modify: `cmake/options.cmake`
- Delete: `data/Makefile`

- [ ] Create `mpm_data_core`.
- [ ] Add conditional executable creation for exactly one active `USE_DATA_*` flag.
- [ ] Route runtime outputs to `build/data/`.
- [ ] Keep the style aligned with `work/CMakeLists.txt`.
