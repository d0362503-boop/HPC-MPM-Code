# HPC-MPM-Code

C++17 MPM/FEM solver for solid mechanics, fluid mechanics, and fluid-structure interaction (FSI). MPI-parallel with PETSc linear algebra and VTK HDF5 visualization output.

Supported solver modes:

- `FLUID` — MPM or FEM fluid solver
- `SOLID` — explicit or implicit MPM solid solver
- `FSI` — partitioned MPM–FEM or monolithic MPM–MPM fluid–solid interaction

## Features

- C++17 with MPI parallelism
- Dynamic load balancing (DLB) for fluid MPM via particle-coordinate sampling (`module/DLB/`)
- PETSc-backed sparse linear algebra
- Generalized-alpha/Newmark parameters and kinematics in `GeneralizedAlphaIntegrator`, owned by each `MaterialPoint`
- Fluid: FEM or MPM
- Solid: explicit or implicit MPM
- Partitioned MPM–FEM FSI and monolithic MPM–MPM FSI with interface multipliers
- CMake-based build with bundled PETSc/HDF5 bootstrap
- Standalone input generators and partitioners under `data/`
- VTK HDF5 output (`*.vtkhdf`)

## Quick Start

```bash
cmake -S . -B build
cmake --build build -j8
cd build
mpirun -np 4 ./MPM
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
| `build/res/` | Runtime visualization output |
| `build/` | Required CMake binary directory |

Key files:

- `CMakeLists.txt` — root CMake entry
- `cmake/options.cmake` — build options
- `AGENTS.md` — engineering rules and pitfalls
- `work/**/main.cpp` — each solver's executable entry point; CMake selects one

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

Active build options (most are defined in `cmake/options.cmake`; `MPM_ENABLE_NATIVE_ARCH` is defined in the root `CMakeLists.txt`):

| Option | Default | Meaning |
|--------|---------|---------|
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

The selected solver's `main.cpp` supplies `main()` directly; there is no separate dispatch entry.
Initialization and finalization are shared through `InitializeSimulation()` and `FinalizeSimulation()` in `module/dataset.h`.

For FSI, `work/src_fsi/CMakeLists.txt` selects either `MPM_FEM` or `MPM_MPM`.
`work/src_solid/CMakeLists.txt` similarly selects `explicit` or `implicit`.

### Selecting data tools

`data/generate/CMakeLists.txt` and `data/divide/CMakeLists.txt` configure all
standalone fluid/solid tools and both FSI variants by default. Build only the
target needed for a case:

```cmake
cmake --build build --target makinput_fsi_mpm_fem makdivide_fsi_mpm_fem -j8
```

## Build Examples

Current MPM--MPM FSI build (selected in `work/src_fsi/CMakeLists.txt`):

```bash
cmake -S . -B build
cmake --build build -j8
```

Fluid build (select `src_fluid` in `work/CMakeLists.txt`):

```bash
cmake -S . -B build
cmake --build build -j8
```

Solid implicit build (select `src_solid` and uncomment `implicit` in `work/src_solid/CMakeLists.txt`):

```bash
cmake -S . -B build
cmake --build build -j8
```

Example fluid data generator:

```bash
cmake --build build --target makinput_fluid -j8
```

## Running the Solver

```bash
cd build
mpirun -np 4 ./MPM
```

At runtime the solver reads an orchestration file, conventionally `file.dat`, with four lines:

1. Parameter file path (e.g., `input.txt`)
2. Grid file prefix (e.g., `griddata`)
3. Point file prefix (e.g., `pointdata`)
4. Output file prefix

Per-rank input files follow the prefix with the rank number, e.g. `griddata0.txt`, `pointdata0.txt`.

`file.dat` is read from the `build/` working directory. Its paths should refer to
generated files under `build/`, such as `./data/divide/...` and `./res/...`, not
source-tree paths under `../data/`. The MPI process count must match the partition
topology, and the input mesh dimensions must match the generated partition data.

A convenience script is available at `build/run.sh`:

```bash
cd build
sh run.sh [NP]
```

`NP` resolves in this order:

1. Command-line argument: `sh run.sh 16`
2. Environment variable: `NP=8 sh run.sh`
3. Default: all logical CPUs (`nproc`)

It launches `./MPM` with `mpirun` under `nohup` in the background (hyper-threading aware via `--use-hwthread-cpus --bind-to hwthread`) and redirects output to:

- `stdout/out_<PID>_<NP>.log`
- `stderr/err_<PID>_<NP>.log`

`build/run.sh` is gitignored; copy it elsewhere if you want to version a customized version.

> **Do not delete `build/`**: it contains generated runtime files, partitioned input data, and job outputs. If a clean CMake configure is needed, remove only `build/CMakeCache.txt` and `build/CMakeFiles/`.

## Data Generator / Partitioner Workflow

Build and run a generator:

```bash
cmake --build build --target makinput_fluid -j8
cd build/data/generate/fluid
./makinput_fluid
```

Build and run a partitioner:

```bash
cmake --build build --target makdivide_fluid -j8
cd build/data/divide/fluid
./makdivide_fluid
```

The partitioner writes rank-split data under `myrank_data/`.

For solid cases, each particle section stores a count, then all coordinate
triplets, then one `id matid surf_point mass vol0` row per particle. Coordinates
and attributes are separate blocks, not interleaved particle records.

## Output

- Visualization: VTK HDF5 (`grid.vtkhdf`, `wp.vtkhdf`, `sp.vtkhdf`)
- Text inputs/outputs: `griddata*.txt`, `pointdata*.txt`
- Restart: per-rank `*_re.txt` files

## Important Notes

Before modifying the solver core, read `AGENTS.md` — in particular §5 (Critical Parallel Rules) and §6 (Common Pitfalls). It is the single source of truth for the project's engineering rules.

## License

See `LICENSE` at the repository root.
