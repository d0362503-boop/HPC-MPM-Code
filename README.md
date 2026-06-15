# HPC-MPM-Code

`HPC-MPM-Code` is a C++17 MPM/FEM solver for solid mechanics, fluid mechanics, and fluid-structure interaction (FSI). The code is MPI-parallel and uses PETSc for linear algebra. Visualization output is standardized on VTK HDF5 (`.vtkhdf`).

This repository currently supports three top-level solver modes:

- `FLUID`
- `SOLID`
- `FSI`

It also contains standalone data generators for building input files for those cases.

## Features

- C++17 solver code with MPI parallelism
- PETSc-backed sparse linear algebra
- Fluid solvers in both `FEM` and `MPM` variants
- Solid solvers in both `EXPLICIT` and `IMPLICIT` variants
- Partitioned FSI workflow
- Built-in third-party dependency bootstrap through CMake
- Standalone data generators under `data/`
- VTK HDF5 output for mesh and particle visualization

## Repository Layout

The top-level directories are:

- `module/`: shared core modules, data structures, IO, solver wrappers, fluid modules, solid modules
- `work/`: top-level solver source selections for `FLUID`, `SOLID`, and `FSI`
- `data/`: standalone input generators for `fluid`, `solid`, and `fsi`
- `cmake/`: build options and third-party dependency bootstrap logic
- `Ext/`: bundled source tarballs for external dependencies
- `external/`: installed third-party dependencies after CMake bootstrap
- `docs/`: notes, plans, and developer-facing documents
- `res/`: resource files

Important root files:

- `CMakeLists.txt`: root CMake entry
- `cmake/options.cmake`: central build switches
- `AGENTS.md`: project-specific engineering rules and pitfalls

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

The project now uses the root CMake workflow.

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

## Build Options

Most build switches live in `cmake/options.cmake`.

### 1. Solver Source Selection

Solver source trees under `work/` are selected directly in:

- `work/CMakeLists.txt`

by commenting or uncommenting the relevant `add_subdirectory(...)` lines.

### 2. Data Generator Selection

Standalone data generators and divide tools are always configured through
`data/CMakeLists.txt`. Select cases directly in:

- `data/generate/CMakeLists.txt`
- `data/divide/CMakeLists.txt`

by commenting or uncommenting the relevant `add_subdirectory(...)` lines.

### 3. Solver Method Selection

Fluid solver method:

- `FLUID_METHOD=FEM`
- `FLUID_METHOD=MPM`

Solid solver method:

- `SOLID_METHOD=EXPLICIT`
- `SOLID_METHOD=IMPLICIT`

### 4. Dependency and Feature Toggles

- `BUILD_PETSC`: build PETSc from bundled tarball
- `USE_HDF5`: enable VTK HDF5 output
- `HDF5_ENABLE_PARALLEL`: build parallel HDF5
- `USE_MPI`: enable MPI support

## Common Build Examples

### Build the default solver configuration

```bash
cmake -S . -B build
cmake --build build -j8
```

### Build a fluid solver configuration

```bash
cmake -S . -B build -DFLUID_METHOD=FEM
cmake --build build -j8
```

### Build a solid implicit configuration

```bash
cmake -S . -B build -DSOLID_METHOD=IMPLICIT
cmake --build build -j8
```

When building from `work/src_solid`, keep the uncommented subdirectory in
`work/src_solid/CMakeLists.txt` consistent with `SOLID_METHOD`.

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

The solver expects input files such as:

- `griddata.txt`
- `pointdata.txt` or case-specific point data file
- `input.txt`

The exact required files depend on the selected case and workflow.

## Data Generator Workflow

The `data/` directory contains standalone tools for generating input files for the three main case families:

- `data/fluid`
- `data/solid`
- `data/fsi`

The generators are built through the root CMake workflow and produce executables in:

```bash
build/data/
```

Typical targets are:

- `makinput_fluid`
- `makinput_solid`
- `makinput_fsi`

Typical run pattern:

```bash
cd build/data
./makinput_fsi
```

This works because CMake copies the matching case `input.txt` into the expected subdirectory under `build/data/`.

For more detail, see `data/README.md`.

## Output and Visualization

The project uses VTK HDF5 for visualization output.

Typical generated files include:

- mesh output: `grid.vtkhdf`
- particle/material point output: `wp.vtkhdf` or `sp.vtkhdf`

The `data/` generators also emit text input files such as:

- `griddata.txt`
- `wpdata.txt`
- `spdata.txt`
- `pointdata.txt`

Legacy `.vtk`, `.vtu`, `.pvtu`, and `.pvd` outputs are not the standard path in the current `data/` workflow.

## Architecture Notes

Some project-specific design choices matter a lot when modifying the code:

- A large amount of shared state lives in inline global variables declared in headers such as `dataset.h`, `mesh.h`, and `mpi_data.h`
- Many functions operate directly on this shared global state
- Refactoring those globals into class members without a broader redesign can easily break behavior or introduce ODR issues

There are also mode-specific source trees:

- `work/src_fluid`
- `work/src_solid`
- `work/src_fsi`

And lower-level module trees:

- `module/fluid/FEM`
- `module/fluid/MPM`
- `module/solid/explicit`
- `module/solid/implicit`
- `module/solver`

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

- root CMake is the primary build path
- bundled tarballs in `Ext/` are the primary dependency source
- `data/` generators are class-based and built through root CMake
- visualization output is centered on VTK HDF5

Some older assumptions may still appear in historical notes or legacy subdirectories, so prefer the current root CMake workflow over older ad hoc build paths when in doubt.

## License

This repository includes a `LICENSE` file at the root. Check it before redistribution or external reuse.
