# Data Generators

`data/` contains the standalone input generators for the `fluid`, `solid`, and `fsi` cases.

## Build

Use the root CMake entry and select exactly one generator in `../cmake/options.cmake`:

- `USE_DATA_FLUID`
- `USE_DATA_SOLID`
- `USE_DATA_FSI`

Example:

```bash
cmake -S . -B build -DUSE_DATA_FLUID=OFF -DUSE_DATA_SOLID=OFF -DUSE_DATA_FSI=ON
cmake --build build --target makinput_fsi -j4
```

Only the selected case subdirectory is configured and built.
The generated Makefile and executable live inside that case directory:

- `build/data/fluid/Makefile`
- `build/data/solid/Makefile`
- `build/data/fsi/Makefile`

The executable is emitted beside the copied `input.txt`:

- `build/data/fluid/makinput_fluid`
- `build/data/solid/makinput_solid`
- `build/data/fsi/makinput_fsi`

## Layout

- `data/data_generator.h`: shared base class for the generator workflow.
- `data/output_util.h`: shared text-header and VTK HDF5 helpers.
- `data/fluid`, `data/solid`, `data/fsi`: case-specific generator classes, `makinput.cpp` entrypoints, and per-case CMake targets.

The base class also owns shared file-opening helpers so all generators fail fast on missing input paths or unwritable output paths.
The common grid-geometry setup at the top of each `BuildData()` is also centralized in the base class and writes directly into the shared inline mesh/dataset state, including `aelemmin` and `aelemmax`.

## Output

Run the selected generator from inside its own build directory, for example:

```bash
cd build/data/fsi
./makinput_fsi
```

Each generator writes:

- text inputs such as `griddata.txt`, `wpdata.txt`, `spdata.txt`, or `pointdata.txt`
- VTK HDF5 files for visualization

Visualization output is standardized on `vtkhdf` only:

- mesh: `grid.vtkhdf`
- particles/material points: `wp.vtkhdf` or `sp.vtkhdf`

No legacy `.vtk`, `.vtu`, `.pvtu`, or `.pvd` output is maintained in `data/`.
