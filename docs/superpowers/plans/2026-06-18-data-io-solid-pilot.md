# Data I/O Solid Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the first working slice of the shared `data/common/` layer by migrating the `solid` generator/divider text I/O path without changing runtime file compatibility.

**Architecture:** Keep the current `makinput_solid` and `makdivide_solid` entrypoints. Add shared schema/data-struct/text-I/O helpers under `data/common/`, then route only the `solid` text path through that layer as a pilot. Leave `fluid` and `fsi` on the old path until the `solid` pilot is stable.

**Tech Stack:** C++17, root CMake data targets, existing `MaterialPoint` and `BoundaryCondition`, existing `module/data_io.h`.

---

## Why This Is The Next Step

- `solid` is the smallest migration surface, so it is the cheapest place to validate the abstraction.
- It already contains `surf_point`, which is a good canary for “field added in one place but forgotten in another.”
- Once `solid` works, the shared layer API will be much easier to apply to `fluid` and `fsi`.

## Target Scope

**In scope**

- `griddata.txt` global mesh writer for `solid`
- `spdata.txt` global point writer for `solid`
- `griddata<rank>.txt` partitioned mesh writer for `solid`
- `spdata<rank>.txt` partitioned point writer for `solid`
- shared schema and shared text serializer/deserializer for the `solid` path

**Out of scope**

- `fluid` and `fsi`
- solver-side runtime readers
- HDF5 migration
- visualization output refactor

---

### Task 1: Audit and freeze the current solid text formats

**Files:**

- Create: `data/common/data_format.h`
- Modify: `data/README.md`
- Reference only: `data/generate/solid/solid_generator.cpp`
- Reference only: `data/divide/data_partitioner.cpp`
- Reference only: `data/divide/divide_solid/solid_divider.cpp`

- [ ] Write down the exact current field order for `spdata.txt`.
- [ ] Write down the exact current field order for partitioned `spdata<rank>.txt`.
- [ ] Write down the exact current field order for `griddata.txt`.
- [ ] Write down the exact current field order for partitioned `griddata<rank>.txt`.
- [ ] Add named internal format constants:
  - `kGridTextFormatVersion`
  - `kSolidPointTextFormatVersion`
- [ ] Document the compatibility rule in `data/README.md`: the pilot must preserve the current file layout byte-for-byte as closely as possible.

**Expected result**

- The text schema stops being implicit and becomes one explicit contract.

---

### Task 2: Introduce shared transport structs for solid mesh and solid point data

**Files:**

- Create: `data/common/mesh_dataset.h`
- Create: `data/common/point_dataset.h`

- [ ] Create `GlobalMeshTextData` for the data written by `WriteGridDataFile(...)`.
- [ ] Create `PartitionedMeshTextData` for the rank-split mesh payload written by `OutputMeshData(...)`.
- [ ] Create `SolidPointTextData` with at least:
  - `num`
  - `coord`
  - `id`
  - `matid`
  - `surf_point`
  - `mass`
  - `vol0`
  - `ubc`
  - `vbc`
  - `wbc`
- [ ] Keep these structs passive. No solver logic, no geometry generation logic, no partition logic.
- [ ] Add brief comments explaining which file each struct corresponds to.

**Expected result**

- There is now one named in-memory representation for each solid text payload.

---

### Task 3: Implement shared solid text serializer/deserializer helpers

**Files:**

- Create: `data/common/text_data_io.h`
- Create: `data/common/text_data_io.cpp`
- Modify: `data/CMakeLists.txt`

- [ ] Add `WriteGlobalMeshTextData(std::ofstream&, const GlobalMeshTextData&)`.
- [ ] Add `ReadGlobalMeshTextData(std::ifstream&, GlobalMeshTextData&)`.
- [ ] Add `WritePartitionedMeshTextData(std::ofstream&, const PartitionedMeshTextData&)`.
- [ ] Add `ReadPartitionedMeshTextData(std::ifstream&, PartitionedMeshTextData&)`.
- [ ] Add `WriteSolidPointTextData(std::ofstream&, const SolidPointTextData&)`.
- [ ] Add `ReadSolidPointTextData(std::ifstream&, SolidPointTextData&)`.
- [ ] Keep all field-order knowledge inside these functions.
- [ ] Reuse `OpenInputFile`, `OpenOutputFile`, `InputVector`, and `OutputVector` where helpful.
- [ ] Add low-risk validation:
  - `num >= 0`
  - vector sizes match `num`
  - BC arrays are internally consistent with `ibc`

**Expected result**

- One source of truth now owns how solid text files are laid out.

---

### Task 4: Migrate the solid generator to write through the shared layer

**Files:**

- Modify: `data/generate/data_generator.h`
- Modify: `data/generate/data_generator.cpp`
- Modify: `data/generate/solid/solid_generator.cpp`

- [ ] Leave `BuildData()` untouched except where data extraction is needed.
- [ ] Add a conversion step from current generator state into `GlobalMeshTextData`.
- [ ] Add a conversion step from current generator state into `SolidPointTextData`.
- [ ] Replace the hand-written `spdata.txt` field output with `WriteSolidPointTextData(...)`.
- [ ] Replace the direct mesh text emission path with `WriteGlobalMeshTextData(...)`.
- [ ] Keep VTK/HDF output unchanged in this task.

**Expected result**

- `makinput_solid` still produces the same files, but no longer owns the field ordering itself.

---

### Task 5: Migrate the solid divider to read and write through the shared layer

**Files:**

- Modify: `data/divide/data_partitioner.h`
- Modify: `data/divide/data_partitioner.cpp`
- Modify: `data/divide/divide_solid/solid_divider.cpp`

- [ ] Replace manual parsing of `spdata.txt` in `SolidDivider::LoadPointData(...)` with `ReadSolidPointTextData(...)`.
- [ ] Convert the shared `SolidPointTextData` into the divider’s existing working state.
- [ ] Replace manual writing of `spdata<rank>.txt` with `WriteSolidPointTextData(...)` or a partitioned-solid variant if needed.
- [ ] Replace direct mesh parsing/writing in `InputMeshData(...)` and `OutputMeshData(...)` with the shared mesh text helpers.
- [ ] Keep partition math and renumbering logic untouched.

**Expected result**

- `makdivide_solid` still partitions exactly the same way, but no longer hand-maintains file layout rules.

---

### Task 6: Verify that the solid pilot is behavior-preserving

**Files:**

- Modify: `data/README.md`
- Optionally create: `docs/data-regression-checklist.md`

- [ ] Build:
  - `cmake --build build --target makinput_solid makdivide_solid -j4`
- [ ] Run `makinput_solid` before the refactor and keep a reference copy of:
  - `griddata.txt`
  - `spdata.txt`
- [ ] Run `makinput_solid` after the refactor and compare against the reference.
- [ ] Run `makdivide_solid` before the refactor and keep a reference copy of:
  - `myrank_data/griddata0.txt`
  - `myrank_data/spdata0.txt`
- [ ] Run `makdivide_solid` after the refactor and compare against the reference.
- [ ] Specifically verify:
  - point count
  - `surf_point`
  - BC count and ordering
  - mass and `vol0`
  - overlap metadata in `griddata<rank>.txt`

**Expected result**

- The new shared layer is proven on one full end-to-end path.

---

### Task 7: Prepare the follow-on migration for fluid and FSI

**Files:**

- Modify: `data/common/README.md`
- Modify: `docs/superpowers/plans/2026-06-17-data-io-unification.md`

- [ ] Record what the `solid` pilot taught us:
  - which helper signatures worked
  - which fields should be generalized
  - which assumptions were too `solid`-specific
- [ ] Decide whether `fluid` can reuse `SolidPointTextData` patterns with a separate struct, or whether one generic point-record abstraction is actually worthwhile.
- [ ] Use the `solid` pilot result to define the next branch scope:
  - `fluid` next
  - `fsi` last

**Expected result**

- The next migrations are smaller because the API was validated on the simplest case first.

---

## Recommended Execution Order

1. Freeze `solid` text schema
2. Create shared structs
3. Add solid serializer/deserializer
4. Migrate solid generator
5. Migrate solid divider
6. Run before/after diff checks
7. Only then expand to `fluid`

## Concrete “First Coding Session” Checklist

- [ ] Create `data/common/`
- [ ] Add `data/common/data_format.h`
- [ ] Add `data/common/mesh_dataset.h`
- [ ] Add `data/common/point_dataset.h`
- [ ] Add `data/common/text_data_io.h`
- [ ] Wire these files into `data/CMakeLists.txt`
- [ ] Do not migrate any call sites yet in the first session
- [ ] End the first session with a compiling skeleton and documented field schema

## Done Criteria For This Pilot

- `solid` generator/divider no longer hard-code solid text field ordering in multiple places
- `spdata.txt` and `spdata<rank>.txt` still carry `surf_point` correctly
- `griddata` and `griddata<rank>` layouts remain runtime-compatible
- `fluid` and `fsi` can now be migrated against a proven shared API instead of guesswork
