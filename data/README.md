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

The `fluid` and `solid` cases are always built. FSI has two variants,
`MPM_FEM` (MPM solid + FEM fluid) and `MPM_MPM` (MPM solid + MPM fluid);
both are enabled by default. To build only one, comment out the
corresponding `add_subdirectory(...)` in `data/generate/fsi/CMakeLists.txt`
and `data/divide/divide_fsi/CMakeLists.txt`.

Executables are emitted per case/variant:

- `build/data/generate/fluid/makinput_fluid`
- `build/data/generate/solid/makinput_solid`
- `build/data/generate/fsi/MPM_FEM/makinput_fsi_mpm_fem`
- `build/data/generate/fsi/MPM_MPM/makinput_fsi_mpm_mpm`
- `build/data/divide/fluid/makdivide_fluid`
- `build/data/divide/solid/makdivide_solid`
- `build/data/divide/fsi/MPM_FEM/makdivide_fsi_mpm_fem`
- `build/data/divide/fsi/MPM_MPM/makdivide_fsi_mpm_mpm`

## Layout

Generator side:

- `data/generate/data_generator.h`: `DataGenerator` base class. `Run()` drives
  the fixed workflow; `BuildData()` is a non-virtual common sequence
  (grid geometry -> mesh -> control points -> `CreateBCs()` ->
  `CreateParticles()`), and case behavior is customized through the virtual
  hooks (`CreateBCs`, `CreateParticles`, `WriteBCData`, `WritePointData`,
  `WriteVisualizationOutputs`).
- `data/generate/output_util.h`: shared text-header and VTK HDF5 helpers.
- `data/generate/fluid`, `data/generate/solid`: per-physics base classes
  `FluidGenerator` / `SolidGenerator` (both `DataGenerator + MaterialPoint`),
  each also serving as the standalone case generator.
- `data/generate/fsi/MPM_FEM`, `data/generate/fsi/MPM_MPM`: per-variant
  per-physics subclasses plus a coordinator (`DataGenerator` subclass) holding
  `solid_` / `fluid_` members; shared BC/point writers are reused from the
  per-physics bases via `using` declarations.

Divider side:

- `data/divide/data_partitioner.h`: `DataPartitioner` base class. `Run()`
  drives the workflow; the partition engine (`MeshPartition`,
  `PointRenumber`) and the stateless `BCRenumber` (public static) live here.
- `data/divide/divide_fluid`, `data/divide/divide_solid`: per-physics
  components `FluidDivider` / `SolidDivider`, each owning its source
  (`points_`) and partitioned (`partition_points_`) containers. The per-rank
  work is split into `PartitionPoints()` (before `MeshPartition`) and
  `PartitionBCs(cp_min, cp_max)` (after `MeshPartition`) so that a coupled
  coordinator can share one mesh windowing pass between physics.
- `data/divide/divide_fsi/MPM_MPM`: `MPMMPMFSISolidDivider` /
  `MPMMPMFSIFluidDivider` inherit the per-physics components and expose their
  hooks via `using`; coordinator `MPMMPMFSIDivider` holds both and delegates.
  The MPM fluid additionally carries three inflow BC blocks
  (`uinfbc/vinfbc/winfbc`), which are **element-indexed** (not control-point
  indexed) and renumbered with the element window.
- `data/divide/divide_fsi/MPM_FEM`: same structure; the FEM fluid has no
  particles, so `MPMFEMFSIFluidDivider` partitions only the four mesh BCs.

Both sides share their base-class implementation through the
`mpm_generate_common` / `mpm_divide_common` OBJECT libraries, compiled once
and linked into every tool.

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

Each solid particle record contains: coordinates, `id`, `matid`,
`surf_point`, `mass`, `vol0`.

No legacy `.vtk`, `.vtu`, `.pvtu`, or `.pvd` output is maintained in `data/`.

For divide runs, execute the tool from its own per-variant build directory:

```bash
cd build/data/divide/fsi/MPM_MPM
./makdivide_fsi_mpm_mpm
```

Each divider reads the matching generator outputs through the per-variant
`para_input_data.txt` copied next to the binary (rank topology, input paths,
and output prefix such as `./myrank_data_column/`), and writes rank-split
`griddata*.txt` / `pointdata*.txt` files under that prefix.
