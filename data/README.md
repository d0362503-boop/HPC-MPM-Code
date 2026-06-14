# Data Tools

`data/` contains the standalone input generators and partition tools for the `fluid`, `solid`, and `fsi` cases.

## Build

Use the root CMake entry and select exactly one generator in `../cmake/options.cmake`:

- `USE_DATA_FLUID`
- `USE_DATA_SOLID`
- `USE_DATA_FSI`

Use the same pattern for the divide tools:

- `USE_DIVIDE_FLUID`
- `USE_DIVIDE_SOLID`
- `USE_DIVIDE_FSI`

Example:

```bash
cmake -S . -B build \
  -DUSE_DATA_FLUID=OFF -DUSE_DATA_SOLID=OFF -DUSE_DATA_FSI=ON \
  -DUSE_DIVIDE_FLUID=OFF -DUSE_DIVIDE_SOLID=OFF -DUSE_DIVIDE_FSI=ON
cmake --build build --target makinput_fsi makdivide_fsi -j4
```

Only the selected case subdirectory is configured and built.

Generator executables are emitted as:

- `build/data/generate/fluid/makinput_fluid`
- `build/data/generate/solid/makinput_solid`
- `build/data/generate/fsi/makinput_fsi`

Divide executables are emitted as:

- `build/data/divide/fluid/makdivide_fluid`
- `build/data/divide/solid/makdivide_solid`
- `build/data/divide/fsi/makdivide_fsi`

## Layout

- `data/generate/data_generator.h`: shared base class for the generator workflow.
- `data/generate/output_util.h`: shared text-header and VTK HDF5 helpers.
- `data/generate/fluid`, `data/generate/solid`, `data/generate/fsi`: case-specific generator classes, `makinput.cpp` entrypoints, and per-case CMake targets.
- `data/divide/data_partitioner.h`: shared base class for the divide workflow and common mesh partition helpers.
- `data/divide/divide_fluid`, `data/divide/divide_solid`, `data/divide/divide_fsi`: case-specific divider classes, `makdivide.cpp` entrypoints, and per-case CMake targets.

Both base classes own shared file-opening helpers so the tools fail fast on missing input paths or unwritable output paths.
The generator base class also centralizes the common grid-geometry setup at the top of each `BuildData()` and writes directly into the shared inline mesh/dataset state, including `aelemmin` and `aelemmax`.

## Output

Run the selected generator from inside its own build directory, for example:

```bash
cd build/data/generate/fsi
./makinput_fsi
```

Each generator writes:

- text inputs such as `griddata.txt`, `wpdata.txt`, `spdata.txt`, or `pointdata.txt`
- VTK HDF5 files for visualization

Visualization output is standardized on `vtkhdf` only:

- mesh: `grid.vtkhdf`
- particles/material points: `wp.vtkhdf` or `sp.vtkhdf`

No legacy `.vtk`, `.vtu`, `.pvtu`, or `.pvd` output is maintained in `data/`.

For divide runs, execute the selected tool from its own build directory, for example:

```bash
cd build/data/divide/fsi
./makdivide_fsi
```

Each divider reads the matching generator outputs through its copied per-case `para_input_data.txt` and writes rank-split text files under `build/data/divide/<case>/myrank_data/`.
