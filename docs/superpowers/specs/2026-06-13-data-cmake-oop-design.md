# Data CMake OOP Design

**Goal**

Refactor `data/` into a small object-oriented generation subsystem with one selectable data target at a time, built only through CMake and emitted under `build/data/`.

**Scope**

- Replace the old single-file `makinput.cpp` organization with:
  - one shared abstract base class
  - one shared output utility layer
  - one concrete generator per mode: `fluid`, `solid`, `fsi`
- Keep `module/` as the backend and do not attempt to class-wrap its global state.
- Remove `data/Makefile` from the supported flow.
- Use `cmake/options.cmake` to choose exactly one active data generator target, following the same style as `work/CMakeLists.txt`.

**Architecture**

`data/` will be split into:

- `data/data_generator.h/.cpp`
  - abstract orchestration layer
- `data/output_util.h/.cpp`
  - shared output helpers
- `data/fluid/fluid_generator.h/.cpp`
- `data/solid/solid_generator.h/.cpp`
- `data/fsi/fsi_generator.h/.cpp`
- existing thin entrypoints:
  - `data/fluid/makinput.cpp`
  - `data/solid/makinput.cpp`
  - `data/fsi/makinput.cpp`

Each `makinput.cpp` becomes a thin `main()` that constructs the corresponding generator and calls `Run()`.

The base class owns the workflow:

1. Initialize runtime state if needed
2. Load case input
3. Build particles / mesh / BC data
4. Write text outputs
5. Write visualization outputs
6. Finalize runtime state if needed

**Naming**

Shared helper APIs in `output_util.h/.cpp` must follow the existing project `.clang-tidy` style:

- types: `CamelCase`
- functions: `CamelCase`
- globals: avoid introducing new ones unless required by legacy backend interop

So helper names such as `WriteGlobalMeshHeader`, `WriteGlobalBcData`, `WriteVtkHdf5Mesh`, and `WriteVtkHdf5Points` are acceptable; lowercase snake-like helper names are not.

**CMake**

`data/CMakeLists.txt` will mirror `work/CMakeLists.txt`:

- define a shared library `mpm_data_core`
- conditionally add only one executable:
  - `makinput_fluid`
  - `makinput_solid`
  - `makinput_fsi`
- select with:
  - `USE_DATA_FLUID`
  - `USE_DATA_SOLID`
  - `USE_DATA_FSI`

Output location:

- `build/data/`

This will be done by setting target runtime output directories explicitly.

**Non-Goals**

- Do not refactor `module/` global data layout.
- Do not merge the three data modes into one runtime-switched executable.
- Do not preserve Makefile as a primary path.
