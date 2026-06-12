# MPM-Code Technical Reference

Detailed reference for MPM-Code architecture, build system, and parallel data structures.
Consult this when working on parallel bugs, PETSc integration, or build-system changes.

## Directory Structure

```
MPM-Code/
├── MPM_main.cpp              # Entry point; selects solver mode by uncommenting calls
├── CMakeLists.txt            # Root CMake configuration
├── Makefile                  # Legacy manual build (still maintained; see below)
├── .clang-tidy               # Google C++ Style naming enforcement
│
├── module/                   # Core shared modules (compiled into mpm_modules)
│   ├── solver/               # Sparse matrix (CrsMat) & iterative solvers (GPBiCG*)
│   ├── fluid/                # Fluid-specific MPM/FEM kernels
│   │   ├── FEM/              # Phase-field VOF, NS-FEM solver
│   │   └── MPM/              # Stabilized mixed MPM fluid solver
│   ├── solid/                # Solid mechanics kernels
│   │   ├── explicit/         # Explicit UL-MPM with F-bar projection
│   │   └── implicit/         # Implicit MPM solid (Generalized-α / Newmark-β)
│   ├── bc.h / bc.cpp         # Boundary condition containers
│   ├── mesh.h / mesh.cpp     # Background grid & control-point generation
│   ├── material_point.h/.cpp # Base MaterialPoint class (P2G/G2P, FLIP/PIC/TPIC/APIC)
│   ├── mpi_data.h/.cpp       # MPI datatype helpers & particle/node communication
│   ├── dataset.h             # Global inline simulation parameters & VectorAssign helpers
│   ├── data_io.h/.cpp        # Text-based I/O routines for checkpoints & VTK
│   ├── map_and_interpolate.h # FLIP / PIC / TPIC / APIC interpolation schemes
│   └── ...
│
├── src_fluid/                # Fluid-only driver (datain_para, main, residual)
├── src_solid/                # Solid-only driver (explicit / implicit variants)
├── src_fsi/                  # FSI driver & coupling logic
│
├── data/                     # Pre-processing tools: grid/particle input generators
│   ├── fluid/                # Fluid benchmark inputs (cavity, dam-break, sloshing)
│   ├── solid/                # Solid benchmark inputs (beam, cook, disk, hyperwall)
│   ├── fsi/                  # FSI benchmark inputs (Turek)
│   └── divide/               # MPI domain-decomposition utility
│
├── vtk/                      # Post-processing: raw text → VTK converter
├── res/                      # Runtime output & restart files (per benchmark)
└── build/                    # CMake / Makefile build artifacts
```

## Build System

### Makefile (Currently Used)

A hand-maintained `Makefile` at the project root is the **current** build system. It uses boolean flags at the top of the file to manually include/exclude source directories:

```makefile
USE_SRC_FLUID = false
USE_SRC_SOLID = false
USE_SRC_FSI   = true
USE_FEM       = true
USE_MPM       = false
USE_EXPLICIT_SOLID = false
USE_IMPLICIT_SOLID = true
```

Build with:
```bash
cd MPM-Code
make -j8        # produces build/MPM.exe
```

This Makefile hard-codes the same PETSc path (`/home/pan/petsc-3.24.5`) and links HYPRE (`-lHYPRE`).

> **Note:** The project is currently built with this Makefile. CMake support exists in the repository but is **not actively used** and will be adopted in future development.

### CMake (Future)

CMake is prepared as the future primary build system. It uses cache variables to select the solver composition at configure time:

```bash
cd MPM-Code/build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPM_SOLVER_MODE=FSI \
  -DFLUID_METHOD=FEM \
  -DSOLID_METHOD=IMPLICIT
make -j$(nproc)
```

| Cache Variable | Allowed Values | Default | Description |
|----------------|----------------|---------|-------------|
| `MPM_SOLVER_MODE` | `FLUID`, `SOLID`, `FSI` | `FSI` | Which high-level driver to compile |
| `FLUID_METHOD` | `FEM`, `MPM` | `FEM` | Fluid discretization inside `module/fluid/` |
| `SOLID_METHOD` | `EXPLICIT`, `IMPLICIT` | `IMPLICIT` | Solid time integrator inside `module/solid/` |
| `CMAKE_BUILD_TYPE` | `Debug`, `Release` | `Release` | Debug adds `-g -O0 -Wall -Wextra`; Release uses `-O3 -march=native` |

CMake automatically sets `DATA_SUBDIR` (`fluid` / `solid` / `fsi`) based on `MPM_SOLVER_MODE`.

> **Note:** The `data/` and `vtk/` subdirectories are **commented out** in the root `CMakeLists.txt` because they have separate compilation issues. Use their standalone Makefiles instead (see below).

### Auxiliary Tools

| Tool | Build Command | Output |
|------|---------------|--------|
| Input generator | `cd data && make` | `makinput.exe` |
| VTK converter | `cd vtk && make` | `makvtk.exe` |
| MPI partitioner | `cd data/divide && make` | `divide.exe` |

## How to Run

The main executable is `build/MPM.exe` (Makefile). CMake produces `MPM` but is not currently in use. It reads plain-text input files (grid + point data) that must be generated beforehand by the pre-processing tools.

Example run script (`build/run.sh`):
```bash
mpiexec -np 2 --bind-to core --map-by core ./MPM
```

### Input File Conventions

- **Grid data** (`griddata.txt`) — Node coordinates, element connectivity, BC node lists.
- **Point data** (`pointdata.txt`) — Particle coordinates, material IDs, initial velocities.
- **Parameter data** (`input.txt` or `para_input_data.txt`) — Time step, density, viscosity, solver switches (`APIC` / `FLIP` / `TPIC`), restart flags, load-ramp steps, gravity vector.

Global parameters are stored as `inline` variables in `module/dataset.h` and populated by `DatainPara()` in the respective `src_*/datain_para.cpp`.

## Code Architecture Details

### Global Inline Variables

The codebase makes heavy use of C++17 `inline` global variables declared in headers:

- `module/dataset.h` — Time-stepping parameters (`dt`, `istep`, `iend`), material properties, file names.
- `module/mesh.h` — Grid dimensions (`node`, `nelem`, `dxy`), node coordinate arrays (`xyn`, `xyc`), control-point offsets (`nuc`, `nvc`, `nwc`, `npc`).
- `module/mpi_data.h` — MPI rank info (`nprocs`, `myrank`), communication buffers (`dbc`, `dbl`), and a template `MPIDatatypeCheck<T>`.

This design means most functions are effectively operating on global state. Be extremely careful when refactoring: moving variables into classes or changing header inclusion order can cause ODR violations or silent breakage.

### Class Hierarchy

```
MaterialPoint
├── SolidMaterialPointBase
│   └── SolidMaterialPoint (ImplicitMPM namespace)
└── (fluid point is just MaterialPoint, accessed via global wp / ifp)
```

- `MaterialPoint` — Holds particle state (mass, volume, stress, deformation gradient), P2G/G2P kernels, and time-integration helpers (Generalized-α, Newmark-β).
- `SolidMaterialPointBase` — Adds constitutive models (rigid, linear elasticity, Neo-Hookean, St. Venant-Kirchhoff, Mooney-Rivlin), rigid-BC detection, and restart I/O.
- `SolidMaterialPoint` (in `ImplicitMPM` namespace) — Implicit solid driver with Newton-Raphson loops, tangent modulus assembly, and a custom `CrsMat` for the solid subsystem.

### Mapping & Interpolation Schemes

Supported P2G/G2P schemes are implemented as small template classes in `module/map_and_interpolate.h`:

- **FLIP** — Particle-in-cell with velocity increment transfer.
- **PIC** — Direct velocity interpolation.
- **TPIC** — Taylor-series PIC with velocity gradient reconstruction.
- **APIC** — Affine PIC with B-matrix and D-matrix inversion.

The active scheme is selected at runtime by the string `solswitch` (e.g., `"APIC"`, `"FLIP"`) read from input.

### Linear Algebra

`module/solver/crsmat.h` defines `CrsMat`, a hand-rolled compressed sparse row matrix that also wraps PETSc objects:

- `BuildCrsMat(num_block)` — Builds local CSR structure from mesh connectivity.
- `BuildLGMAP()` / `BuildPetscMat()` / `BuildKSPSolver()` — PETSc setup.
- `SolveWithPetsc(ndof)` — KSP solve with gathered BC handling.
- Custom iterative solvers (`GPBiCG`, `GPBiCGSafe`, `GPBiCGAR`) exist as fallbacks when `use_petsc = false`.

> **Known issue:** The PETSc parallel path has a documented ownership bug where shared overlap nodes (`dbc < 1`) are excluded from `interior_list`, causing inconsistent row ownership across ranks. See `petsc_parallel_debug_report_v2.md` in the repo root for a detailed analysis.

## Security & Safety Considerations

- **No network exposure:** This is an offline HPC/scientific computing executable.
- **Input validation is minimal:** `datain_para.cpp` reads plain-text files with `operator>>`. Malformed input can cause out-of-bounds vector access or silent NaN propagation because sizes are read from the file itself.
- **Memory safety:** The code uses raw `std::vector` indexing extensively with few bounds checks. `-fsanitize=address,undefined` can be added to `CXXFLAGS` for debug builds.
- **File I/O:** Restart and output routines write to paths constructed from input strings. Ensure the working directory is isolated because the code will overwrite existing files without prompting.

## Parallel Data Structures & Mesh Partitioning

This section documents the MPI domain-decomposition data structures read by `InputParaGriddata()` in `module/data_io.cpp` and used by `NodeVarComm()` in `module/mpi_data.h`. Understanding these is **essential** for any parallel bug fix or PETSc integration work.

### Shape-Function Order (`idimc`)

- `idimc[i]` (read from `parafile`) is the **shape-function order** in direction *i*:
  - `idimc = 1` → **Linear** shape functions (standard FEM).
  - `idimc = 2` → **Quadratic B-Spline** shape functions (used in higher-order MPM/FEM).
- The control-point count in each direction is:  
  `xynodec[i] = xyelem[i] + idimc[i]`.
- When `idimc = 1`, control points coincide with visualization nodes (except at partition boundaries where ghost layers differ). When `idimc = 2`, the control-point stencil is wider than the visualization-node stencil.

### Two Node Types: `node` vs `nodec`

| Symbol | Header | Meaning | Formula |
|--------|--------|---------|---------|
| `node` | `mesh.h` | **Visualization nodes** (post-processing / VTK) | `xynode[0] * xynode[1] * xynode[2]` |
| `nodec`| `mesh.h` | **Control points** (real DOFs, including ghost) | `xynodec[0] * xynodec[1] * xynodec[2]` |

**Critical rule:** All FEM/MPM **computations** and **PETSc vectors/matrices** must use `nodec`, never `node`. The `node` arrays are only for output. Mixing them (e.g. using `node` to compute PETSc local row counts) causes global-size mismatches and `MatAssemblyEnd` crashes.

### `InputParaGriddata()` File Format

Each rank reads its own `gridfile + myrank + ".txt"`, written by `divide.exe`. The file contains, in order:

```
myrank
nelem  xyelem[0] xyelem[1] xyelem[2]
node   xynode[0] xynode[1] xynode[2]
nodec  xynodec[0] xynodec[1] xynodec[2]
nxyr[0] nxyr[1] nxyr[2]
xymin[0] xymin[1] xymin[2]
xymax[0] xymax[1] xymax[2]
aelemmin[0] aelemmin[1] aelemmin[2]
aelemmax[0] aelemmax[1] aelemmax[2]
isb
```

If `isb != 0` (multi-process with overlap), the following are appended:

| Variable | Size | Description |
|----------|------|-------------|
| `naid[isb]` | `isb` | IDs of overlapping neighbor ranks |
| `isubc` | scalar | **Total** number of overlapping control points |
| `nsbc[isb]` | `isb` | Overlapping control-point count **per neighbor** |
| `nsubc[isubc]` | `isubc` | Local IDs of overlapping control points (grouped by neighbor) |
| `dbn[nodec]` | `nodec` | How many ranks share each control point |
| `dbc[nodec*4]` | `nodec*4` | Control-point weights: `dbc[i+d*nodec] = 1.0 / dbn[i]` for `d=0,1,2,3` (u,v,w,p) |
| `isubl` | scalar | Total number of overlapping visualization nodes |
| `nsbl[isb]` | `isb` | Overlapping visualization-node count per neighbor |
| `nsubl[isubl]` | `isubl` | Local IDs of overlapping visualization nodes |
| `dbl[node*4]` | `node*4` | Visualization-node weights (same pattern as `dbc`) |

If `isb == 0` (single-process or no overlap), `naid`, `nsbc`, `nsubc`, `nsbl`, `nsubl` are empty and `dbc`/`dbl` are filled with `1.0`.

### `aelemmin` / `aelemmax` — Global Index Anchors

- `aelemmin[i]` = global starting element index of this rank in direction *i*.
- `aelemmax[i]` = global ending element index.
- These are used to map local control-point indices to global indices and are the basis for the **tie-break** rule when multiple ranks share a control point (the rank with the smallest `aelemmin` owns it).

### `NodeVarComm` — Ghost Synchronization

```cpp
template <typename T>
void NodeVarComm(std::vector<T>& dat, int nn);
```

This is the core MPI communication primitive for control-point data:

1. **Pack phase:** For each neighbor `i`, extract `dat[nsubc[j] + nn]` into send buffer `sdn`.
2. **Exchange phase:** Non-blocking send/receive (`MPI_Isend`/`MPI_Irecv`) to `naid[i]` with count `nsbc[i]`.
3. **Wait phase:** `MPI_Wait` on all requests.
4. **Unpack phase:** `dat[nsubc[j] + nn] += rvn[j]` — **adds** the received ghost contribution.

Because the operation is `+=`, each rank must later multiply by `dbc` (or divide by `dbn`) to get the correct averaged value. The `dbc` array is pre-computed as `1.0/dbn` so a single multiplication yields the final value.

**Important:** `nn` is an offset into `dat`. For example, `NodeVarComm(nvof, 0)` syncs the first `nodec` entries, while `NodeVarComm(nvel, nuc)` syncs the velocity `u` component starting at index `nuc`.

### Quick Reference: Parallel Variable Locations

| Variable | Declared In | Populated By | Used By |
|----------|-------------|--------------|---------|
| `node`, `nodec`, `nelem`, `xyelem`, `xynode`, `xynodec` | `mesh.h` | `InputParaGriddata()` | Mesh build, matrix assembly, I/O |
| `aelemmin`, `aelemmax` | `mesh.h` | `InputParaGriddata()` | Global L2G mapping, tie-break |
| `nxyr` | `mpi_data.h` | `InputParaGriddata()` | Partition geometry checks |
| `isb`, `isubc`, `isubl` | `mpi_data.h` | `InputParaGriddata()` | `NodeVarComm`, `ParticleCommunication` |
| `naid`, `nsbc`, `nsubc` | `mpi_data.h` | `InputParaGriddata()` | `NodeVarComm` (control-point ghost sync) |
| `nsbl`, `nsubl` | `mpi_data.h` | `InputParaGriddata()` | `NodeVarCommLowerOrder` (visualization-node ghost sync) |
| `dbc` | `mpi_data.h` | `InputParaGriddata()` | FEM/MPM assembly (overlap weighting), `SolveWithPetsc` |
| `dbl` | `mpi_data.h` | `InputParaGriddata()` | Visualization output weighting |

## Common Pitfalls for Agents

1. **Do not assume class encapsulation.** Many critical variables are global `inline` variables in headers. grep for `inline` in `*.h` files to find them.
2. **MPI communication is custom.** Do not replace `NodeVarComm` or `ParticleCommunication` with naive `MPI_Allreduce` without understanding the overlap-node weighting (`dbc`).
3. **PETSc ownership is subtle.** The `dbc` vector encodes *overlap weight* (reciprocal of shared-rank count), not a simple boolean ownership mask. Parallel correctness depends on this.
4. **CMake vs Makefile drift.** When adding new `.cpp` files, ensure both `CMakeLists.txt` and the root `Makefile` know about them. The CMake uses `GLOB` / `GLOB_RECURSE`, but the Makefile requires manual `find` logic.
5. **The `data/` and `vtk/` directories are excluded from CMake.** If you modify them, you must build them with their own Makefiles.
6. **Chinese comments in Makefiles / CMake.** Some build-system files contain Chinese comments for compiler selection and flag explanations. The C++ source code itself is commented in English.
7. **Old data structures in `vtk/` and `data/` subdirectories.** The post-processing (`vtk/`) and pre-processing (`data/`) tools were left behind during a codebase-wide migration from `std::vector<std::vector<double>>` (SoA) to `std::vector<std::array<double, 3>>` (AoS). If you encounter compilation errors there, check for:
   - Transposed indices: `vel[0][nid]` should be `vel[nid][0]`.
   - Obsolete resize calls: `coord.resize(3, vector<double>(N))` should be `coord.resize(N)`.
   - Missing `VectorAssign`/`InputVector` overloads for `std::vector<std::array<T, N>>`.
   - The `data/divide/` FSI partitioner (`divide_fsi/`) is additionally out of sync with global inline variables (`sp`, `wfem`) and may require interface refactoring beyond simple data-structure fixes.
