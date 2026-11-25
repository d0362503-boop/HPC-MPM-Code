# HPC-MPM-Code

C++17 MPM/FEM solver for solid mechanics, fluid mechanics, and fluid-structure interaction (FSI), parallelized with MPI and PETSc.

- `FLUID`: MPM or FEM fluid solver, with dynamic load balancing for fluid MPM
- `SOLID`: explicit or implicit MPM solid solver
- `FSI`: partitioned MPM–FEM or monolithic MPM–MPM fluid–solid interaction

Shared modules provide generalized-alpha/Newmark integration and VTK HDF5 visualization output.

## Build

Requires a C++17-capable GCC, MPI, CMake ≥ 3.20, and Zlib. CMake bootstraps PETSc 3.24.5 and HDF5 1.14.5 from `Ext/` into `external/`, reusing existing installations there. Hypre 3.0.0 is optional during the PETSc build.

### Solver executables

Every driver is built as its own executable in its own directory under `build/work/`, next to that driver's `file.dat`; no CMake edit is needed to switch solvers:

| Driver source | Executable |
|---------------|------------|
| `work/src_fluid/MPM` | `build/work/src_fluid/MPM/MPM_FLUID` |
| `work/src_solid/explicit` | `build/work/src_solid/explicit/MPM_SOLID_EXPLICIT` |
| `work/src_solid/implicit` | `build/work/src_solid/implicit/MPM_SOLID_IMPLICIT` |
| `work/src_fsi/MPM_FEM` | `build/work/src_fsi/MPM_FEM/MPM_FEM_FSI` |
| `work/src_fsi/MPM_MPM` | `build/work/src_fsi/MPM_MPM/MPM_MPM_FSI` |

`build/MPM` is a symlink to the default driver for the older `build/file.dat` workflow.
Select it with `-DMPM_DEFAULT_DRIVER=<executable>` at configure time (default `MPM_MPM_FSI`).
To build one driver only, use `cmake --build build --target <executable>`.

### Configure and compile

Run from the repository root:

```bash
cmake -S . -B build
cmake --build build -j8
```

The binary directory must be named `build`; other names are rejected. The default is a Release build with `-O3 -DNDEBUG` and `-march=native`. For cross-architecture or cluster builds, configure with `-DMPM_ENABLE_NATIVE_ARCH=OFF`.

Build options are defined in [cmake/options.cmake](cmake/options.cmake), except `MPM_ENABLE_NATIVE_ARCH` in the root [CMakeLists.txt](CMakeLists.txt):

| Option | Default | Purpose |
|--------|---------|---------|
| `MPM_ENABLE_NATIVE_ARCH` | `ON` | Enable `-march=native` |
| `BUILD_PETSC` | `ON` | Build bundled PETSc |
| `USE_HDF5` | `ON` | Enable VTK HDF5 output |
| `HDF5_ENABLE_PARALLEL` | `ON` | Build parallel HDF5 |
| `USE_MPI` | `ON` | Enable MPI |

> **Do not delete `build/`.** It contains generated inputs, partition data, outputs, and local runtime state. For a clean configure, remove only `build/CMakeCache.txt` and `build/CMakeFiles/`.

CMake copies parameter files (`input.txt`, `para_input_data.txt`) from `data/` into `build/data/` only when the build copy does not exist yet, so runtime edits under `build/` survive re-configures. Delete a build copy and re-configure to get the source version again.

## Prepare Input Data

Build and run the generator and partitioner for the selected case. For example, from the repository root:

```bash
cmake --build build --target makinput_fluid makdivide_fluid -j8
(cd build/data/generate/fluid && ./makinput_fluid)
(cd build/data/divide/fluid && ./makdivide_fluid)
```

Generators and partitioners for fluid, solid, and both FSI variants are configured by default. Select the targets needed for your case in [data/generate/](data/generate/CMakeLists.txt) and [data/divide/](data/divide/CMakeLists.txt); for example, `makinput_fsi_mpm_fem` and `makdivide_fsi_mpm_fem`. Partitioners write rank-split data to the output prefix configured in `para_input_data.txt`.

See the [data tools guide](data/README.md) for parameter files, particle formats, and restart inputs.

## Run

Run from the driver's directory under `build/work/` after generating and partitioning the case data:

```bash
cd build/work/src_fsi/MPM_MPM
mpirun -np 4 ./MPM_MPM_FSI  # example: input partitioned for four ranks
```

Each driver directory has its own orchestration file `file.dat`, which conventionally contains four lines:

1. Parameter file path
2. Grid file prefix
3. Point file prefix
4. Output file prefix

Per-rank files append the rank number to the input prefix, for example `griddata0.txt` and `pointdata0.txt`. Paths in `file.dat` are relative to the driver directory and must refer to runtime files under `build/`, such as `../../../data/divide/...` and `../../../res/...`, never to source paths in the repository `data/`.

**The MPI process count must match the partition topology, and the input mesh dimensions must match the generated partition data.**

### Background launcher and logging

The host-local, gitignored `build/run.sh` launches a driver in the background. Pass the driver directory as `RUN_DIR`;
the script runs the single `MPM_*` executable it finds there (falling back to `build/MPM`) with that directory's `file.dat`:

```bash
# from build/, with matching four-rank input
sh run.sh 4 work/src_fsi/MPM_MPM
```

Use the process-count argument or `NP`. The local script falls back to `nproc` if neither is supplied; this is valid only when it matches the input partition count. Cluster copies may require `NP` and use different CPU-binding policies, so inspect the installed script.

The launcher prints the process count, PID, and log paths, and inherits `PETSC_OPTIONS`:

- `<RUN_DIR>/stdout/out_<PID>_<NP>.log`
- `<RUN_DIR>/stderr/err_<PID>_<NP>.log`

For event timings, launch with `PETSC_OPTIONS="-log_view" sh run.sh 4`.
PETSc prints the summary at normal finalization. Event times can overlap and should not be summed as independent costs; ordinary NR records contain residuals and iteration counts, not separate KSP setup/solve timings. Options and rebuilt code affect new processes only.

## Output

- Visualization: `*.vtkhdf` under the configured output prefix, usually in `build/res/`
- Text inputs/outputs: `griddata*.txt`, `pointdata*.txt`
- Restart data: per-rank `*_re.txt` files

## Repository and Developer Guides

| Path | Purpose |
|------|---------|
| `module/` | Shared mesh, material points, mechanics, I/O, MPI communication, and PETSc wrapper |
| `work/` | Fluid, solid, and FSI executable drivers |
| `data/` | Input generators and partitioners |
| `cmake/` | Build options and dependency bootstrap |
| `Ext/`, `external/` | Bundled sources and installed dependencies |
| `docs/` | Developer notes and plans |
| `build/` | Required build directory and runtime working state |

Implementation details are maintained in the relevant guides:

- [Shared linear solver](module/solver/README.md): equation layout, compression, preconditioners, and rebuild policy
- [Implicit core](module/implicit/README.md): generalized-α integrator, implicit MPM kinematics, and the Newton driver
- [Implicit solid](module/solid/implicit/README.md): nonlinear solve and solid preconditioning
- [MPM–MPM FSI](work/src_fsi/MPM_MPM/README.md): coupled Newton loop, interface handling, and solver strategy
- [Dynamic load balancing](module/DLB/README.md): fluid MPM partitioning
- [AGENTS.md](AGENTS.md): engineering conventions and critical parallel rules; read before modifying the solver core

## License

See [LICENSE](LICENSE).
