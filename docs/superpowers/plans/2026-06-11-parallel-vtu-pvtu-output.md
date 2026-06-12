# Parallel PVD/PVTU/VTU Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current text-shard plus legacy `.vtk` post-processing workflow with direct MPI-parallel `.vtu + .pvtu` binary output plus a `.pvd` time-series index, starting from the active FSI path.

**Architecture:** Keep all solver-side field generation unchanged and only replace the serialization layer. First validate the XML/binary layout using the existing post-processing path, then switch runtime output to per-rank `.vtu` piece files plus a rank-0 `.pvtu` index file.

**Tech Stack:** C++17, MPI, std::ofstream, existing global mesh/material-point data structures, ParaView XML UnstructuredGrid format.

---

## Current Baseline

- Runtime FSI output currently writes per-rank text shards:
  - `module/fluid/FEM/fluid_fem_data_io.cpp`
  - `module/solid/solid_data_io.cpp`
- Runtime FSI triggers are in:
  - `src_fsi/main.cpp`
- Legacy binary VTK writer is in:
  - `module/vtk.h`
  - `module/vtk.cpp`
- Serial post-processing stitchers are in:
  - `vtk/vtk_main.cpp`
  - `vtk/fluid/makvtk_wmesh.cpp`
  - `vtk/fluid/makvtk_wp.cpp`
  - `vtk/solid/makvtk_sp.cpp`

## Non-Negotiable Rules

- Use `nodec` and control-point fields for computation only; never change solver logic to compute on `node`.
- Use `node` / `xyn` only as the visualization mesh when writing fluid mesh output.
- Do not change PETSc ownership, overlap, or `NodeVarComm` behavior as part of this task.
- Keep old text output available behind a switch until the new `.vtu/.pvtu` path is verified.

## Recommended Execution Order

### Task 1: Freeze the Output Contract

**Files:**
- Read: `src_fsi/main.cpp`
- Read: `module/fluid/FEM/fluid_fem_data_io.cpp`
- Read: `module/solid/solid_data_io.cpp`
- Read: `module/vtk.h`
- Read: `module/vtk.cpp`

- [ ] Record the exact datasets that must survive the migration.
- [ ] Treat fluid mesh output as `UnstructuredGrid` with hexahedron cells on visualization nodes.
- [ ] Treat solid output as `UnstructuredGrid` with vertex cells on material points.
- [ ] Mark fluid particle output as lower priority unless standalone FLUID-MPM mode is in scope for the first pass.
- [ ] Write down target file naming, for example `fluid/res-w-<step>_r<rank>.vtu`, `fluid/res-w-<step>.pvtu`, `solid/res-s-<step>_r<rank>.vtu`, `solid/res-s-<step>.pvtu`, plus top-level `fluid/results-w.pvd` and `solid/results-s.pvd`.

### Task 2: Decide the Migration Strategy

**Files:**
- Read: `vtk/vtk_main.cpp`
- Read: `vtk/fluid/makvtk_wmesh.cpp`
- Read: `vtk/solid/makvtk_sp.cpp`

- [ ] Use a two-stage migration.
- [ ] Stage A: add `.vtu/.pvtu` generation to the existing `vtk/` post-processing path.
- [ ] Stage B: remove the text intermediate from runtime and write `.vtu` pieces directly in MPI ranks.
- [ ] Keep Stage A as a regression oracle until Stage B is validated.

### Task 3: Add a New XML Writer Layer

**Files:**
- Create: `module/vtk_xml.h`
- Create: `module/vtk_xml.cpp`
- Modify: `Makefile`
- Modify: `CMakeLists.txt` if present in active maintenance

- [ ] Implement a small writer dedicated to VTK XML `UnstructuredGrid`.
- [ ] Support appended binary payloads rather than legacy VTK 3.0 records.
- [ ] Support scalar arrays, 3-component vector arrays, connectivity, offsets, and cell types.
- [ ] Support writing one rank-local piece at a time without global gather.
- [ ] Support writing the rank-0 `.pvtu` file that references all piece files for one timestep.

### Task 4: Validate XML Output via Existing Post-Processor First

**Files:**
- Modify: `vtk/fluid/makvtk_wmesh.cpp`
- Modify: `vtk/solid/makvtk_sp.cpp`
- Optional later: `vtk/fluid/makvtk_wp.cpp`
- Modify: `vtk/vtk_main.cpp`

- [ ] Reuse existing text shard readers to emit one `.vtu` per rank instead of one merged `.vtk`.
- [ ] Emit one `.pvtu` on the host side after all piece files are written.
- [ ] Emit one `.pvd` that references the timestep `.pvtu` files in monotonically increasing output order.
- [ ] Open the `.pvtu` in ParaView and verify that all ranks load as one dataset.
- [ ] Compare scalar ranges and visible geometry against the old merged `.vtk`.
- [ ] Do this for fluid mesh first, then solid particles.

### Task 5: Wire Direct Runtime Output for Fluid Mesh

**Files:**
- Modify: `module/fluid/FEM/fluid_fem_data_io.cpp`
- Possibly create helper calls in `module/vtk_xml.cpp`
- Read: `src_fsi/main.cpp`

- [ ] Keep `Cp2NodeVTK()` as the field-preparation step.
- [ ] Replace or supplement `OutputMeshData()` so each rank writes its own fluid mesh `.vtu`.
- [ ] Restrict each rank’s piece geometry to its local visualization extent derived from `aelemmin/aelemmax`.
- [ ] Make rank 0 write the matching `.pvtu` after a barrier or after all ranks finish piece output.
- [ ] Make rank 0 append or rewrite the fluid `.pvd` so ParaView can open the full time series from one file.
- [ ] Keep the old `-w.txt` path behind a runtime flag until verified.

### Task 6: Wire Direct Runtime Output for Solid Particles

**Files:**
- Modify: `module/solid/solid_data_io.cpp`
- Read: `src_fsi/main.cpp`

- [ ] Replace or supplement `OutputPointData()` so each rank writes its own solid particle `.vtu`.
- [ ] Map each particle to a `VTK_VERTEX` cell in the local piece.
- [ ] Write the corresponding rank-0 `.pvtu`.
- [ ] Append or rewrite the solid `.pvd` from rank 0.
- [ ] Keep old `-s.txt` output behind a flag during transition.

### Task 7: Add Runtime Output Controls

**Files:**
- Modify: `module/dataset.h`
- Modify: relevant input parsing files under `module/*data_io*.cpp`

- [ ] Add a simple output mode switch, for example `legacy_txt`, `legacy_vtk`, `xml_vtu`.
- [ ] Default to the legacy path until the new writer is stable.
- [ ] Add directory and filename prefix controls if needed, but avoid over-designing this.

### Task 8: Verification Sequence

**Files:**
- Run-time verification only

- [ ] Serial check: 1 rank, one timestep, fluid mesh `.vtu/.pvtu` matches old `.vtk`.
- [ ] Serial check: 1 rank, one timestep, solid point `.vtu/.pvtu` matches old `.vtk`.
- [ ] Time-series check: open the `.pvd` in ParaView and verify timestep navigation works.
- [ ] MPI check: 2 ranks, verify no missing or duplicated blocks at partition boundaries.
- [ ] MPI check: 4+ ranks, verify ParaView opens `.pvtu` directly and arrays are intact.
- [ ] Boundary check: confirm visualization uses `node` geometry while field values still come from `nodec` interpolation.
- [ ] Performance check: compare wall time and output size against text-shard plus post-process workflow.

### Task 9: Expand to Other Modes Only After FSI Is Stable

**Files:**
- Modify later: `module/fluid/MPM/fluid_mpm_data_io.cpp`
- Modify later: `src_fluid/MPM/main.cpp`
- Modify later: `src_solid/implicit/main.cpp`

- [ ] Extend the same writer to FLUID-MPM particle output if needed.
- [ ] Extend the same writer to standalone SOLID mode.
- [ ] Only remove the `vtk/` post-processing executable after all active modes have a direct runtime writer.

## Practical Priority

- First priority: FSI fluid mesh output.
- Second priority: FSI solid particle output.
- Third priority: optional FLUID-MPM particle output.
- Fourth priority: none; `.pvd` is part of the baseline target architecture.

## Acceptance Criteria

- ParaView can open one `.pvtu` per timestep directly with no manual merge.
- ParaView can also open one `.pvd` per field family and scrub the whole time sequence.
- Each MPI rank writes only its own piece file; no rank gathers full global arrays to write output.
- The solver still computes on `nodec` and only visualizes on `node`.
- Old and new outputs agree on geometry extents and primary fields within normal floating-point tolerance.
- The old text-to-`vtk` toolchain can be removed only after the new runtime path is stable.
