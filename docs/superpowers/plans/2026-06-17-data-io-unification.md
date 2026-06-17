# Data I/O Unification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract the shared data structures and text I/O rules used by `data/generate` and `data/divide` into a single reusable core without merging the two tool flows.

**Architecture:** Keep `makinput_*` and `makdivide_*` as separate frontends. Introduce a shared `data/common/` layer that owns file-format schema, mesh/point datasets, and serialization/deserialization helpers. Migrate `solid`, `fluid`, and `fsi` case code to use the shared core incrementally, preserving existing text file compatibility at every stage.

**Tech Stack:** C++17, existing `module/data_io.h`, existing `MaterialPoint`/`BoundaryCondition`, root CMake data tool targets, optional later extension to HDF5.

---

## Target File Structure

**Create**

- `data/common/data_format.h`
- `data/common/data_format.cpp`
- `data/common/mesh_dataset.h`
- `data/common/point_dataset.h`
- `data/common/text_data_io.h`
- `data/common/text_data_io.cpp`
- `data/common/README.md`

**Modify**

- `data/CMakeLists.txt`
- `data/generate/data_generator.h`
- `data/generate/data_generator.cpp`
- `data/generate/solid/solid_generator.cpp`
- `data/generate/fluid/fluid_generator.cpp`
- `data/generate/fsi/fsi_generator.cpp`
- `data/divide/data_partitioner.h`
- `data/divide/data_partitioner.cpp`
- `data/divide/divide_solid/solid_divider.cpp`
- `data/divide/divide_fluid/fluid_divider.cpp`
- `data/divide/divide_fsi/fsi_divider.cpp`
- `data/README.md`

**Keep unchanged in the first pass**

- solver-side runtime readers under `module/` and `work/src_*`
- output filenames and current text file compatibility
- `makinput_*` / `makdivide_*` entrypoints

---

### Task 1: Freeze the existing text format as an explicit schema

**Files:**

- Create: `data/common/data_format.h`
- Create: `data/common/data_format.cpp`
- Modify: `data/README.md`

- [ ] List every current text payload that crosses the `generate -> divide -> runtime` boundary.
- [ ] Document exact field order for:
  - `griddata.txt`
  - `griddata<rank>.txt`
  - `wpdata.txt`
  - `spdata.txt`
  - `pointdata.txt`
- [ ] Add a small schema layer with named concepts instead of relying on scattered field order knowledge.
- [ ] Introduce version constants, for example:
  - `kGridTextFormatVersion`
  - `kSolidPointTextFormatVersion`
  - `kFluidPointTextFormatVersion`
  - `kFsiPointTextFormatVersion`
- [ ] Decide the compatibility rule:
  - first pass keeps existing file layout unchanged
  - version values are internal first, not yet written into files if that would break compatibility
- [ ] Update `data/README.md` so the format becomes canonical documentation rather than tribal knowledge.

**Verification**

- Confirm the schema document matches the current generators and dividers line-by-line.
- Confirm no runtime reader behavior changes in this task.

---

### Task 2: Define shared in-memory datasets for mesh and points

**Files:**

- Create: `data/common/mesh_dataset.h`
- Create: `data/common/point_dataset.h`
- Modify: `data/common/README.md`

- [ ] Define a mesh dataset type for the common information currently duplicated between generator and divider.
- [ ] Include only transport-level data, not solver behavior.
- [ ] Recommended split:
  - `GlobalMeshTextData`
  - `PartitionedMeshTextData`
  - `SolidPointTextData`
  - `FluidPointTextData`
  - `FsiPointTextData`
- [ ] Keep these types as plain structs with obvious field names.
- [ ] Explicitly include fields that have already bitten the codebase before, such as:
  - `surf_point`
  - per-point `id`
  - `matid`
  - `mass`
  - `vol0`
  - fluid state vectors
  - solid and fluid BC payloads
- [ ] Make these data structs independent of `DataGenerator` and `DataPartitioner` lifecycles so both flows can reuse them.

**Verification**

- Confirm every field currently written by a generator has a corresponding named field in the shared structs.
- Confirm every field currently read by a divider can be populated from the shared structs.

---

### Task 3: Build one shared text I/O layer

**Files:**

- Create: `data/common/text_data_io.h`
- Create: `data/common/text_data_io.cpp`
- Modify: `data/CMakeLists.txt`

- [ ] Move text serialization rules into one place.
- [ ] Split helpers by payload type instead of by tool:
  - `ReadGlobalMeshTextData(...)`
  - `WriteGlobalMeshTextData(...)`
  - `ReadPartitionedMeshTextData(...)`
  - `WritePartitionedMeshTextData(...)`
  - `ReadSolidPointTextData(...)`
  - `WriteSolidPointTextData(...)`
  - `ReadFluidPointTextData(...)`
  - `WriteFluidPointTextData(...)`
  - `ReadFsiPointTextData(...)`
  - `WriteFsiPointTextData(...)`
- [ ] Keep `OpenInputFile` / `OpenOutputFile` reuse where appropriate instead of duplicating stream error handling.
- [ ] Put all line-order assumptions in these helpers only.
- [ ] Add format sanity checks that do not break current files:
  - count must be non-negative
  - vector sizes must match declared counts
  - BC sizes must match `ibc`
  - FSI point sections must remain consistent with solid-point count and fluid-point count expectations

**Verification**

- Re-run a generator and divider pair for one case and diff outputs before and after this task.
- Expected result: byte-identical or formatting-identical text outputs.

---

### Task 4: Make generators write through shared I/O instead of hand-coded field order

**Files:**

- Modify: `data/generate/data_generator.h`
- Modify: `data/generate/data_generator.cpp`
- Modify: `data/generate/solid/solid_generator.cpp`
- Modify: `data/generate/fluid/fluid_generator.cpp`
- Modify: `data/generate/fsi/fsi_generator.cpp`

- [ ] Keep geometry generation and initial-condition logic in each case generator.
- [ ] Remove direct knowledge of text field order from case generators.
- [ ] Add a conversion boundary:
  - case generator builds native point/mesh state
  - native state is converted into shared text dataset structs
  - shared text I/O layer writes files
- [ ] Preserve `WriteBcData()` only if it still expresses case-specific content that cannot be normalized yet.
- [ ] If possible, replace `WriteBcData()` with shared dataset population plus one serializer call.
- [ ] Do not touch visualization output in this task unless it blocks the text-I/O cleanup.

**Verification**

- Build:
  - `cmake --build build --target makinput_solid makinput_fluid makinput_fsi -j4`
- Run each generator in its own build directory.
- Check that produced text files remain compatible with existing dividers.

---

### Task 5: Make dividers read and write through shared I/O

**Files:**

- Modify: `data/divide/data_partitioner.h`
- Modify: `data/divide/data_partitioner.cpp`
- Modify: `data/divide/divide_solid/solid_divider.cpp`
- Modify: `data/divide/divide_fluid/fluid_divider.cpp`
- Modify: `data/divide/divide_fsi/fsi_divider.cpp`

- [ ] Keep partition logic in `DataPartitioner` and case divider classes.
- [ ] Remove direct stream parsing for point payloads from case divider implementations.
- [ ] Replace with:
  - shared I/O reads full text dataset
  - divider converts shared dataset into existing partition workflow
  - shared I/O writes partitioned outputs
- [ ] Centralize the logic for preserving field order in split `griddata*.txt` and `spdata/wpdata/pointdata*.txt`.
- [ ] Ensure `surf_point` and any future extra fields cannot be forgotten in one case but not another.

**Verification**

- Build:
  - `cmake --build build --target makdivide_solid makdivide_fluid makdivide_fsi -j4`
- Run each divider against data generated by Task 4 outputs.
- Spot-check `myrank_data/*.txt` structure against prior known-good outputs.

---

### Task 6: Add round-trip regression checks for file compatibility

**Files:**

- Modify: `data/README.md`
- Optionally create: `docs/data-regression-checklist.md`

- [ ] Pick one stable sample case each for `solid`, `fluid`, and `fsi`.
- [ ] Record a repeatable regression workflow:
  - run generator
  - run divider
  - compare file counts
  - compare key file headers
  - compare selected point records
- [ ] If you do not want to add a full unit-test framework in `data/`, at least add a documented shell-based regression procedure.
- [ ] Capture a checklist for the fragile fields:
  - `surf_point`
  - BC record counts
  - overlap-node metadata
  - rank-split point counts

**Verification**

- Run the regression checklist once for all three cases.
- Save representative diffs or “no diff” evidence in your local notes or commit message.

---

### Task 7: Clean up old duplication after the shared layer is proven

**Files:**

- Modify: generator/divider case sources only where duplicated write/read logic remains
- Modify: `data/common/README.md`

- [ ] Remove obsolete helper code only after shared I/O is already in use and verified.
- [ ] Delete dead formatting code from case implementations.
- [ ] Keep one source of truth for each file layout.
- [ ] Document extension rules for adding a new point field:
  - update shared dataset struct
  - update shared serializer/deserializer
  - update generator conversion
  - update divider conversion
  - update docs

**Verification**

- Run full rebuild of all data tools.
- Re-run one end-to-end case per physics mode.

---

## Recommended Sequencing

1. Task 1: schema freeze
2. Task 2: shared datasets
3. Task 3: shared text I/O
4. Task 4: migrate generators
5. Task 5: migrate dividers
6. Task 6: regression checks
7. Task 7: duplication cleanup

This order keeps the risk low because behavior is frozen before abstractions are introduced.

## Non-Goals for This Plan

- Do not merge `generate` and `divide` executables into one tool.
- Do not rewrite runtime solver input in the same branch.
- Do not switch text bulk data to HDF5 in this first pass.
- Do not refactor global inline solver state as part of this work.

## Key Design Rules

- Shared layer owns format knowledge.
- Case classes own physics-specific data preparation and partition rules.
- Runtime solver readers are downstream clients and should not be broken by this refactor.
- Every new field must be added once in the shared schema and then propagated outward.

## Risk Notes

- `fsi` is the highest-risk case because it mixes solid and fluid payload conventions.
- `surf_point` is a good canary field because it was recently added and is easy to drop accidentally.
- Mesh overlap metadata in partitioned `griddata*.txt` must remain exact or downstream PETSc/MPI behavior can silently break.
- If round-trip compatibility is not checked continuously, this refactor can “look clean” while changing file semantics.

## Done Criteria

- `generate` and `divide` no longer each maintain their own independent text field-order logic.
- All current text outputs remain compatible with the runtime readers.
- `solid`, `fluid`, and `fsi` all build and run through generator+divider successfully.
- The rules for adding a new point field are documented and only require touching the shared schema/I-O path once.
