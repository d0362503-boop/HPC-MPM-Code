# Data Tools

`data/` contains the standalone input generators (`makinput_*`) and partition
tools (`makdivide_*`) for the `fluid`, `solid`, and `fsi` cases.

## Build

Use the root CMake entry. `data/CMakeLists.txt` always configures both
`generate/` and `divide/`.

```bash
cmake -S . -B build
cmake --build build --target makinput_fsi_mpm_mpm makdivide_fsi_mpm_mpm -j4
```

Fluid, solid, and both FSI variants (`MPM_FEM` and `MPM_MPM`) are configured
by default. Build only the targets needed for a case. FSI variants are selected in
[generate/fsi/CMakeLists.txt](generate/fsi/CMakeLists.txt) and
[divide/divide_fsi/CMakeLists.txt](divide/divide_fsi/CMakeLists.txt).

| Case | Generator under `build/data/generate/` | Partitioner under `build/data/divide/` |
|------|---------------------------------------|----------------------------------------|
| Fluid | `fluid/makinput_fluid` | `fluid/makdivide_fluid` |
| Solid | `solid/makinput_solid` | `solid/makdivide_solid` |
| MPM–FEM FSI | `fsi/MPM_FEM/makinput_fsi_mpm_fem` | `fsi/MPM_FEM/makdivide_fsi_mpm_fem` |
| MPM–MPM FSI | `fsi/MPM_MPM/makinput_fsi_mpm_mpm` | `fsi/MPM_MPM/makdivide_fsi_mpm_mpm` |

## Structure

- `generate/data_generator.h`: `DataGenerator::Run()` drives grid geometry,
  mesh, control points, BCs, and particle creation. Virtual hooks supply
  case-specific BC, particle, and visualization output.
- `generate/output_util.h`: shared text headers and VTK HDF5 helpers.
- `generate/fluid/` and `generate/solid/`: per-physics bases and standalone
  generators. FSI coordinators hold one solid and one fluid component.
- `divide/data_partitioner.h`: shared mesh partitioning and point/BC
  renumbering. `FluidDivider` and `SolidDivider` own their point containers;
  `PartitionPoints()` precedes mesh windowing and `PartitionBCs()` follows it.
- FSI dividers share one mesh-windowing pass. MPM fluid inflow BCs
  (`uinfbc/vinfbc/winfbc`) are element-indexed; FEM fluid has no particles
  and partitions only its four mesh BCs.
- `mpm_generate_common` and `mpm_divide_common` OBJECT libraries compile
  the common implementations once for all tools.

## Output

Run each tool from inside its own build directory, for example:

```bash
cd build/data/generate/fsi/MPM_MPM
./makinput_fsi_mpm_mpm
```

Each generator writes:

- `griddata.txt` (mesh + boundary conditions) and `pointdata.txt` (particles)
- VTK HDF5 files for visualization: `grid.vtkhdf` (mesh) and
  `sp.vtkhdf` / `wp.vtkhdf` (solid / fluid particles where applicable)

Each solid particle section contains the particle count, all coordinate triplets,
then one `id matid surf_point mass vol0` row per particle. The coordinates and
attributes occupy separate blocks; readers must preserve this ordering.

For divide runs, execute the tool from its own per-variant build directory:

```bash
cd build/data/divide/fsi/MPM_MPM
./makdivide_fsi_mpm_mpm
```

Each divider reads the matching generator outputs through the per-variant
`para_input_data.txt` copied next to the binary (rank topology, input paths,
and output prefix such as `./myrank_data_column/`), and writes rank-split
`griddata*.txt` / `pointdata*.txt` files under that prefix.

## Runtime parameters and reconfiguration

CMake copies generator parameter files such as `input.txt` from the source
directories into the matching `build/data/generate/` directories using
`configure_file(... COPYONLY)`. Reconfiguring can overwrite edits made only to
those runtime copies. Check the first active parameter block, time step,
restart mode and start step after configuration and before launching a solver.

The solver reads paths from `build/file.dat`; keep them inside `build/`.
Fresh runs use initial point data, while restart modes read the existing
per-rank restart files. Restart files are outputs of the solver, not products
of rerunning a partitioner. Rebuilding executables alone does not regenerate
the grid or initial particles; run the data tools explicitly when those inputs
need to change.
