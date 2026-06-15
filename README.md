# HPC-MPM-Code

`HPC-MPM-Code` is a C++17 MPM/FEM solver for solid mechanics, fluid mechanics, and fluid-structure interaction (FSI). The code is MPI-parallel and uses PETSc for linear algebra. Visualization output is standardized on VTK HDF5 (`.vtkhdf`).

This repository currently supports three top-level solver modes:

- `FLUID`
- `SOLID`
- `FSI`

It also contains standalone data generators and partitioners for building the input files consumed by those cases.

## Features

- C++17 solver code with MPI parallelism
- PETSc-backed sparse linear algebra
- Fluid solvers in both `FEM` and `MPM` variants
- Solid solvers in both `EXPLICIT` and `IMPLICIT` variants
- Partitioned strong-coupling FSI workflow
- Built-in third-party dependency bootstrap through CMake
- Standalone data generators and partitioners under `data/`
- VTK HDF5 output for mesh and particle visualization
- Compiler optimization enabled by default (`Release` build with `-march=native`)

## Repository Layout

The top-level directories are:

- `module/`: shared core modules, data structures, IO, solver wrappers, fluid modules, solid modules
- `work/`: top-level solver source selections for `FLUID`, `SOLID`, and `FSI`
- `data/`: standalone input generators (`data/generate/`) and partitioners (`data/divide/`) for `fluid`, `solid`, and `fsi`
- `cmake/`: build options and third-party dependency bootstrap logic
- `Ext/`: bundled source tarballs for external dependencies
- `external/`: installed third-party dependencies after CMake bootstrap
- `docs/`: notes, plans, and developer-facing documents
- `res/`: resource files referenced at runtime (e.g., `res/Turek/`)
- `build/`: required CMake binary directory

Important root files:

- `CMakeLists.txt`: root CMake entry
- `cmake/options.cmake`: central build option definitions
- `AGENTS.md`: project-specific engineering rules and pitfalls
- `MPM_main.cpp`: root solver entry point that dispatches to the selected driver

## Dependencies

The project expects:

- GCC with C++17 support
- MPI
- CMake 3.20 or newer
- Zlib

The repository can bootstrap these bundled dependencies automatically:

- `PETSc 3.24.5` from `Ext/petsc-3.24.5.tar.gz`
- `HDF5 1.14.5` from `Ext/hdf5-hdf5_1.14.5.tar.gz`
- optional `Hypre 3.0.0` from `Ext/hypre-3.0.0.tar.gz` during PETSc build

By default:

- PETSc is installed into `external/petsc`
- HDF5 is installed into `external/hdf5`

If those installs already exist, the CMake bootstrap logic skips rebuilding them.

## Build System

The project uses the root CMake workflow. There is no root `Makefile`; the root `CMakeLists.txt` is the primary build path.

One important repository convention is that the CMake binary directory must be named exactly `build`. Configuring into any other directory is rejected by the root `CMakeLists.txt`.

Typical configure step:

```bash
cmake -S . -B build
```

Typical build step:

```bash
cmake --build build -j8
```

This produces the main solver executable:

```bash
build/MPM
```

The default build type is `Release`, and `-march=native` is added to the C++ flags, so the default compile flags are effectively `-O3 -DNDEBUG -march=native`.

## Build Options

Most build switches live in `cmake/options.cmake`.

### 1. Solver Source Selection

The actual solver source tree is selected by uncommenting exactly one `add_subdirectory(...)` line in `work/CMakeLists.txt`:

```cmake
# add_subdirectory(src_fluid)
# add_subdirectory(src_solid)
add_subdirectory(src_fsi)
```

The `USE_SRC_FSI`, `USE_SRC_FLUID`, and `USE_SRC_SOLID` options in `cmake/options.cmake` are legacy toggles and are **not** the active selection mechanism.

You must also keep the call in `MPM_main.cpp` consistent with the selected source tree. For example, when building `work/src_fsi`, `MPM_main.cpp` should call `MPMBlockFSI()`:

```cpp
MPMBlockFSI();
```

### 2. Data Generator / Partitioner Selection

Similarly, data generators are selected in `data/generate/CMakeLists.txt`:

```cmake
# add_subdirectory(fluid)
# add_subdirectory(solid)
add_subdirectory(fsi)
```

And partitioners are selected in `data/divide/CMakeLists.txt`:

```cmake
# add_subdirectory(divide_fluid fluid)
# add_subdirectory(divide_solid solid)
add_subdirectory(divide_fsi fsi)
```

The `USE_DATA_*` and `USE_DIVIDE_*` options in `cmake/options.cmake` are legacy toggles and are **not** the active selection mechanism.

### 3. Solver Method Selection

Fluid solver method:

- `FLUID_METHOD=FEM`
- `FLUID_METHOD=MPM`

Solid solver method:

- `SOLID_METHOD=EXPLICIT`
- `SOLID_METHOD=IMPLICIT`

When building from `work/src_solid`, keep the uncommented subdirectory in `work/src_solid/CMakeLists.txt` consistent with `SOLID_METHOD`.

### 4. Dependency and Feature Toggles

- `BUILD_PETSC`: build PETSc from bundled tarball
- `USE_HDF5`: enable VTK HDF5 output
- `HDF5_ENABLE_PARALLEL`: build parallel HDF5
- `USE_MPI`: enable MPI support

## Common Build Examples

### Build the default FSI solver configuration

```bash
cmake -S . -B build
cmake --build build -j8
```

### Build a fluid solver configuration

Edit `work/CMakeLists.txt` to uncomment `add_subdirectory(src_fluid)` and edit `MPM_main.cpp` to call `StabilizedMixedMPM()`. Then:

```bash
cmake -S . -B build -DFLUID_METHOD=FEM
cmake --build build -j8
```

### Build a solid implicit configuration

Edit `work/CMakeLists.txt` to uncomment `add_subdirectory(src_solid)` and edit `MPM_main.cpp` to call `Solid_implicit_ULMPM()`. Then:

```bash
cmake -S . -B build -DSOLID_METHOD=IMPLICIT
cmake --build build -j8
```

### Build the FSI data generator

```bash
cmake -S . -B build
cmake --build build --target makinput_fsi -j8
```

## Running the Solver

The main executable is intended to run under MPI:

```bash
mpiexec -np 4 ./build/MPM
```

At runtime the solver reads an orchestration file, conventionally named `file.dat`, containing four lines:

1. Parameter file path (e.g., `input.txt`)
2. Grid file prefix (e.g., `griddata`)
3. Point file prefix (e.g., `pointdata`)
4. Output file prefix

The parameter file (`input.txt`) contains:

- Mapping scheme (`PIC`/`FLIP`/`TPIC`/`APIC`), restart flag, nonlinear-step count
- Start/end/output-step indices
- Time step, tolerance, spectral radius / Newmark parameters
- Domain bounds, element counts, B-spline orders (`idimc`)
- Particles per element per direction (`npxye`)
- Material densities, constitutive model flags and properties
- Body-force vector (`bb`)

Per-rank grid and point data are read from files such as:

- `griddata<rank>.txt`
- `pointdata<rank>.txt` or `wpdata<rank>.txt` / `spdata<rank>.txt`

The exact required files depend on the selected case and workflow.

## Data Generator / Partitioner Workflow

The `data/generate/` directory contains tools for generating input files for the three main case families:

- `data/generate/fluid`
- `data/generate/solid`
- `data/generate/fsi`

The partitioners in `data/divide/` split generated inputs into per-rank files under `myrank_data/`.

### Generator workflow

1. Edit `data/generate/CMakeLists.txt` to uncomment the desired case.
2. Build the target:

```bash
cmake --build build --target makinput_<case> -j8
```

3. Run from the case build directory, for example:

```bash
cd build/data/generate/fsi
./makinput_fsi
```

### Partitioner workflow

1. Edit `data/divide/CMakeLists.txt` to uncomment the desired case.
2. Build the target:

```bash
cmake --build build --target makdivide_<case> -j8
```

3. Run from the case build directory; it reads generator outputs and writes rank-split files under `myrank_data/`.

For more detail, see `data/README.md`.

## Output and Visualization

The project uses VTK HDF5 for visualization output.

Typical generated files include:

- mesh output: `grid.vtkhdf`
- particle/material point output: `wp.vtkhdf` (fluid) or `sp.vtkhdf` (solid)

The `data/` generators also emit text input files such as:

- `griddata.txt`
- `wpdata.txt`
- `spdata.txt`
- `pointdata.txt`

Restart output is written as plain-text per-rank files (`*_re.txt`) containing nodal/particle state.

Legacy `.vtk`, `.vtu`, `.pvtu`, and `.pvd` outputs are not the standard path in the current `data/` workflow.

## Architecture Notes

Some project-specific design choices matter a lot when modifying the code:

- A large amount of shared state lives in inline global variables declared in headers such as `module/dataset.h`, `module/mesh.h`, and `module/mpi_data.h`
- Many functions operate directly on this shared global state
- Refactoring those globals into class members without a broader redesign can easily break behavior or introduce ODR issues

The solver mode source trees are:

- `work/src_fluid`
- `work/src_solid`
- `work/src_fsi`

And the lower-level module trees are:

- `module/fluid/FEM`
- `module/fluid/MPM`
- `module/solid/explicit`
- `module/solid/implicit`
- `module/solver`

FSI is a partitioned strong-coupling solver with block iterations, Anderson relaxation, and FSI interface detection based on solid volume fraction.

## Parallel and Numerical Caveats

Several project rules are important enough to call out explicitly:

- Use `nodec` for computation and PETSc ownership logic; `node` is for visualization and mesh layout, not solver ownership
- `NodeVarComm` behavior is not equivalent to a naive `MPI_Allreduce`
- `dbc` stores overlap weights, not a boolean ownership mask
- Shared control-point ownership uses the smallest `aelemmin` tie-break rule

If you are touching parallel assembly or control-point communication, read `AGENTS.md` first.

## Development Notes

- Match the repository style rules in `.clang-tidy`
- In class member functions, call other member functions with explicit `this->`
- When adding new `.cpp` files, keep CMake targets in sync
- Be careful with `data/divide_fsi/`, which is not fully synchronized with the current global inline variable layout

## Known Current State

At the repository level, the current direction is:

- root CMake is the primary and only build path
- bundled tarballs in `Ext/` are the primary dependency source
- `data/` generators and partitioners are built through root CMake
- visualization output is centered on VTK HDF5
- solver source selection is done by uncommenting the relevant `add_subdirectory(...)` lines in `work/CMakeLists.txt` and keeping `MPM_main.cpp` consistent

Some older assumptions may still appear in historical notes or legacy subdirectories, so prefer the current root CMake workflow over older ad hoc build paths when in doubt.

## License

This repository includes a `LICENSE` file at the root. Check it before redistribution or external reuse.
