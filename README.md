# HPC-MPM-Code

C++17 MPM/FEM solver for solid mechanics, fluid mechanics, and fluid-structure interaction (FSI). MPI-parallel with PETSc linear algebra and VTK HDF5 visualization output.

Supported solver modes:

- `FLUID` — MPM or FEM fluid solver
- `SOLID` — explicit or implicit MPM solid solver
- `FSI` — partitioned strong-coupling fluid–solid interaction

## Features

- C++17 with MPI parallelism
- PETSc-backed sparse linear algebra
- Fluid: FEM or MPM
- Solid: explicit or implicit MPM
- Partitioned FSI with block iterations
- CMake-based build with bundled PETSc/HDF5 bootstrap
- Standalone input generators and partitioners under `data/`
- VTK HDF5 output (`*.vtkhdf`)

## Quick Start

```bash
cmake -S . -B build
cmake --build build -j8
mpiexec -np 4 ./build/MPM
```

The CMake binary directory **must** be named `build`; anything else is rejected.

Default compile flags are `-O3 -DNDEBUG` with `-march=native` enabled by default. Use `-DMPM_ENABLE_NATIVE_ARCH=OFF` for cross-architecture or cluster builds.

## Repository Layout

| Directory | Purpose |
|-----------|---------|
| `module/` | Shared core modules: mesh, material points, I/O, MPI comms, PETSc wrapper, fluid/solid mechanics implementations |
| `work/` | Top-level solver drivers: `src_fluid/`, `src_solid/`, `src_fsi/` |
| `data/` | Input generators (`data/generate/`) and partitioners (`data/divide/`) |
| `cmake/` | Build options and bundled-dependency bootstrap |
| `Ext/` | Bundled source tarballs (PETSc, HDF5, optional Hypre) |
| `external/` | Installed dependencies after bootstrap |
| `docs/` | Developer notes and plans |
| `res/` | Runtime resource files |
| `build/` | Required CMake binary directory |

Key files:

- `CMakeLists.txt` — root CMake entry
- `cmake/options.cmake` — build options
- `AGENTS.md` — engineering rules and pitfalls
- `MPM_main.cpp` — solver entry point; dispatches to the selected driver

## Dependencies

Required:

- GCC with C++17 support
- MPI
- CMake ≥ 3.20
- Zlib

Auto-bootstrapped from `Ext/`:

- PETSc 3.24.5 → `external/petsc`
- HDF5 1.14.5 → `external/hdf5`
- Hypre 3.0.0 (optional, during PETSc build)

If `external/petsc` and `external/hdf5` already exist, CMake skips rebuilding them.

## Build Options

Active options in `cmake/options.cmake`:

| Option | Default | Meaning |
|--------|---------|---------|
| `FLUID_METHOD` | `FEM` | `FEM` or `MPM` |
| `SOLID_METHOD` | `IMPLICIT` | `EXPLICIT` or `IMPLICIT` |
| `MPM_ENABLE_NATIVE_ARCH` | `ON` | Enable `-march=native` |
| `BUILD_PETSC` | `ON` | Build PETSc from bundled tarball |
| `USE_HDF5` | `ON` | Enable VTK HDF5 output |
| `HDF5_ENABLE_PARALLEL` | `ON` | Build parallel HDF5 |
| `USE_MPI` | `ON` | Enable MPI |

### Selecting the solver source

The active driver is selected by uncommenting **exactly one** `add_subdirectory(...)` line in `work/CMakeLists.txt`:

```cmake
# add_subdirectory(src_fluid)
# add_subdirectory(src_solid)
add_subdirectory(src_fsi)
```

You must also keep `MPM_main.cpp` consistent. For example, with `src_fsi` selected:

```cpp
MPMBlockFSI();
```

`work/src_solid/CMakeLists.txt` also requires its uncommented subdirectory to match `SOLID_METHOD`.

### Selecting data tools

Generators are selected in `data/generate/CMakeLists.txt`:

```cmake
# add_subdirectory(fluid)
# add_subdirectory(solid)
add_subdirectory(fsi)
```

Partitioners are selected in `data/divide/CMakeLists.txt`:

```cmake
# add_subdirectory(divide_fluid fluid)
# add_subdirectory(divide_solid solid)
add_subdirectory(divide_fsi fsi)
```

## Build Examples

Default FSI build:

```bash
cmake -S . -B build
cmake --build build -j8
```

Fluid build:

```bash
cmake -S . -B build -DFLUID_METHOD=FEM
cmake --build build -j8
```

Solid implicit build:

```bash
cmake -S . -B build -DSOLID_METHOD=IMPLICIT
cmake --build build -j8
```

FSI data generator:

```bash
cmake --build build --target makinput_fsi -j8
```

## Running the Solver

```bash
mpiexec -np 4 ./build/MPM
```

At runtime the solver reads an orchestration file, conventionally `file.dat`, with four lines:

1. Parameter file path (e.g., `input.txt`)
2. Grid file prefix (e.g., `griddata`)
3. Point file prefix (e.g., `pointdata`)
4. Output file prefix

Per-rank input files follow the prefix with the rank number, e.g. `griddata0.txt`, `pointdata0.txt`.

A convenience script is available at `build/run.sh` (gitignored; copy it if you want to keep a version in the repo):

```bash
cd build
sh run.sh [NP]
```

It runs `./MPM` under `nohup` in the background, uses all logical CPUs by default (hyper-threading aware via `--use-hwthread-cpus --bind-to hwthread`), and redirects output to `stdout/out_<PID>_<NP>.log` and errors to `stderr/err_<PID>_<NP>.log`. The number of ranks can be passed as an argument or via the `NP` environment variable.

## Data Generator / Partitioner Workflow

Build and run a generator:

```bash
cmake --build build --target makinput_fsi -j8
cd build/data/generate/fsi
./makinput_fsi
```

Build and run a partitioner:

```bash
cmake --build build --target makdivide_fsi -j8
cd build/data/divide/fsi
./makdivide_fsi
```

The partitioner writes rank-split data under `myrank_data/`.

## Output

- Visualization: VTK HDF5 (`grid.vtkhdf`, `wp.vtkhdf`, `sp.vtkhdf`)
- Text inputs/outputs: `griddata*.txt`, `wpdata*.txt`, `spdata*.txt`, `pointdata*.txt`
- Restart: per-rank `*_re.txt` files

## Important Notes

- A lot of shared state lives in inline global variables in `module/dataset.h`, `module/mesh.h`, and `module/mpi_data.h`. Refactoring them into classes risks ODR issues and behavior changes.
- Use `nodec` (control points) for computation and PETSc; `node` is for visualization only.
- `NodeVarComm` ghost exchange is not equivalent to a naive `MPI_Allreduce`.
- `dbc` stores overlap weights (`1/shared-rank-count`), not a boolean mask.
- Shared control-point ownership tie-break: smallest `aelemmin`.

See `AGENTS.md` for the full list of engineering rules and pitfalls.

## License

See `LICENSE` at the repository root.
