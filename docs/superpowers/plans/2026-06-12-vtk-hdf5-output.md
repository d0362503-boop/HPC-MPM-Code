# VTK HDF5 Output Implementation Plan

**Goal:** Add a native VTK HDF5 (`.vtkhdf`) output path to the MPM-Code solver, leveraging the HDF5 1.14.5 dependency already bundled in `Ext/` and built by CMake. The new writer should eventually replace the legacy binary `.vtk` post-processing toolchain for FSI runs while keeping the old path available behind a runtime switch.

**Tech Stack:** C++17, MPI, HDF5 C API (parallel enabled), existing global mesh/material-point data structures, ParaView VTK HDF5 reader.

---

## 1. Background & Format Choice

### Why VTK HDF5?
- VTK HDF5 stores the same VTK data model as `.vtu`/`.vtk` but inside an HDF5 container.
- Better compression, chunked storage, and parallel I/O via HDF5 MPI-IO.
- ParaView natively opens `.vtkhdf` files; parallel collections can be opened with `.pvtkhdf`.
- A single rank-local `.vtkhdf` file is self-contained (geometry + fields), avoiding the text-shard merge step.

### VTK HDF5 UnstructuredGrid Layout (target schema)

```
/VTKHDF/
  Version           int64[2]        = [2, 0]
  Type              string          = "UnstructuredGrid"
  NumberOfPoints    int64
  NumberOfCells     int64
  Points            float64[n_points, 3]
  Connectivity      int64[...]      // cell connectivity, VTK-ordered
  Offsets           int64[n_cells+1]
  Types             uint8[n_cells]
  PointData/
    <array_name>    <type>[n_points] or [n_points, 3]
  CellData/
    <array_name>    <type>[n_cells]   or [n_cells, 3]
```

For parallel runs each rank writes its own piece file:

```
res-w-<step>_r<rank>.vtkhdf   // per-rank fluid mesh piece
res-s-<step>_r<rank>.vtkhdf   // per-rank solid particle piece
```

Rank 0 writes a lightweight `.pvtkhdf` XML index referencing the pieces.

---

## 2. Non-Negotiable Rules

- Keep all solver-side field preparation unchanged; only replace serialization.
- Use `nodec` for computation; use `node`/`xyn` only for visualization geometry.
- Do not change PETSc ownership, overlap, or `NodeVarComm` behavior.
- Keep legacy `.vtk`/text output behind a switch until the new path is verified.
- HDF5 calls must be guarded by `#ifdef HAVE_HDF5` so the code still compiles when `USE_HDF5=false`.
- All HDF5 errors must be checked and reported via `H5Eprint` or a small wrapper that aborts with a readable message.

---

## 3. Recommended Execution Order

### Task 1: Create a Minimal VTK HDF5 Writer Core

**Files to create:**
- `module/vtk_hdf5.h`
- `module/vtk_hdf5.cpp`

**Sub-tasks:**
- [ ] Implement a small RAII wrapper around HDF5 ids (`hid_t`) to prevent leaks.
- [ ] Implement helper functions:
  - `WriteHdf5AttributeInt(group, name, value)`
  - `WriteHdf5AttributeString(group, name, value)`
  - `WriteHdf5Dataset1D(group, name, data, hdf5_type)`
  - `WriteHdf5Dataset2D(group, name, data, rows, cols, hdf5_type)`
- [ ] Implement `WriteVtkHdf5UnstructuredGrid(...)` that writes the full `/VTKHDF` group from raw arrays:
  - points
  - connectivity
  - offsets
  - cell types
  - point data arrays
  - cell data arrays
- [ ] Default to HDF5 serial file creation (`H5Fcreate` with default plist); later upgrade to parallel plist.

**Verification:**
- [ ] Unit smoke test: create a 2-cell hex mesh with one scalar and one vector array, write to `test.vtkhdf`, open in ParaView.

---

### Task 2: Adapt Existing Output Data to VTK HDF5

**Files to read:**
- `module/fluid/FEM/fluid_fem_data_io.cpp`
- `module/solid/solid_data_io.cpp`
- `module/vtk.h` / `module/vtk.cpp`

**Sub-tasks:**
- [ ] Identify the exact arrays currently written by `OutputMeshData()` (fluid) and `OutputPointData()` (solid):
  - Fluid mesh: coordinates (`xyn`), connectivity (`nc`), velocity, pressure, phi, etc.
  - Solid particles: coordinates (`point.coord`), matid, velocity, stress, etc.
- [ ] Add converter functions that copy these arrays into contiguous `std::vector` layouts matching the VTK HDF5 schema.
- [ ] For fluid mesh: build `connectivity` (8-node hex cells), `offsets` (`0, 8, 16, ...`), and `types` (all `12` = `VTK_HEXAHEDRON`).
- [ ] For solid particles: build `connectivity` (`0, 1, 2, ...`), `offsets` (`0, 1, 2, ...`), and `types` (all `1` = `VTK_VERTEX`).

**Verification:**
- [ ] Single-rank FSI run writes `.vtkhdf` files that ParaView renders identically to the old `.vtk`.

---

### Task 3: Integrate into Runtime Output Path

**Files to modify:**
- `module/fluid/FEM/fluid_fem_data_io.cpp`
- `module/solid/solid_data_io.cpp`
- `module/dataset.h`
- relevant input parsing files

**Sub-tasks:**
- [ ] Add an output mode switch, e.g. `output_format = legacy_vtk | vtu | vtkhdf`.
- [ ] Default to `legacy_vtk` until VTK HDF5 is verified.
- [ ] When `vtkhdf` is selected:
  - Each rank writes its own piece file.
  - File naming: `vtk/<case>/res-w-<istep>_r<myrank>.vtkhdf` and `vtk/<case>/res-s-<istep>_r<myrank>.vtkhdf`.
- [ ] Keep the old text-shard path untouched when `legacy_vtk` is selected.

---

### Task 4: Add Parallel Collection Index (`.pvtkhdf`)

**Files to modify:**
- `module/vtk_hdf5.cpp`
- `src_fsi/main.cpp` (rank-0 orchestration)

**Sub-tasks:**
- [ ] After all ranks finish writing piece files, rank 0 writes `vtk/<case>/res-w-<istep>.pvtkhdf` and `res-s-<istep>.pvtkhdf`.
- [ ] The index file is a small XML file listing each piece with its file name and bounds (optional).
- [ ] Optionally write a top-level `.pvd` time-series file so ParaView can scrub all timesteps.

**Verification:**
- [ ] ParaView opens `.pvtkhdf` for a multi-rank run and shows a single assembled dataset.

---

### Task 5: Enable Parallel HDF5 Writes (Optional but Recommended)

**Files to modify:**
- `module/vtk_hdf5.cpp`

**Sub-tasks:**
- [ ] When running with MPI and `HDF5_ENABLE_PARALLEL=ON`, create files with `H5Pset_fapl_mpio`.
- [ ] Use collective `H5Dwrite` calls if all ranks contribute to a single shared file.
- [ ] Start with the simpler "one file per rank" approach; switch to shared-file only if benchmarking shows a benefit.

---

### Task 6: Update Build Systems

**Files to modify:**
- `Makefile`
- `CMakeLists.txt`
- `module/CMakeLists.txt`

**Sub-tasks:**
- [ ] Add `module/vtk_hdf5.cpp` to both Makefile source lists and `module/CMakeLists.txt`.
- [ ] Confirm `HAVE_HDF5` macro is visible to the new source file.
- [ ] Re-run `make test-hdf5` and a full CMake build.

---

### Task 7: Verification Sequence

**Checks to perform:**
- [ ] Smoke test: 1 rank, small mesh, write/read back with HDF5 C API, assert array sizes and values.
- [ ] ParaView check: 1 rank, one timestep, `.vtkhdf` geometry and fields match legacy `.vtk`.
- [ ] MPI check: 2 ranks, verify no missing or duplicated blocks at partition boundaries.
- [ ] MPI check: 4+ ranks, `.pvtkhdf` opens as a single dataset.
- [ ] Performance check: compare wall time and file size against legacy `.vtk`.
- [ ] Switch check: confirm `USE_HDF5=false` still compiles and legacy output still works.

---

## 4. Acceptance Criteria

- ParaView can open a single `.vtkhdf` file written by one rank.
- ParaView can open a `.pvtkhdf` index file for a multi-rank run and display the assembled dataset.
- Old and new outputs agree on geometry extents and primary fields within floating-point tolerance.
- The solver compiles with `USE_HDF5=true` and `USE_HDF5=false`.
- The old `.vtk`/text output remains functional when the new format is disabled.

---

## 5. Practical Priority

1. **Task 1** (writer core) — highest; blocks everything else.
2. **Task 2** (data adapters) — highest; connects writer to solver data.
3. **Task 3** (runtime integration) — high; makes output usable in real runs.
4. **Task 4** (`.pvtkhdf` index) — medium; needed for comfortable multi-rank visualization.
5. **Task 5** (parallel HDF5) — low; one-file-per-rank is sufficient for first pass.
6. **Task 6** (build systems) — high; must keep both Makefile and CMake working.
7. **Task 7** (verification) — high; confirms correctness.
