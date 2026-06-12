# HDF5 Auto-Build CMake Design

**Date:** 2026-06-12

**Goal**

Add an optional CMake-only dependency bootstrap path that builds HDF5 from a source tarball bundled with the project, while keeping PETSc external and leaving the current Makefile workflow unchanged.

## Scope

This design covers only:

- project-local HDF5 source tarball management
- optional HDF5 build orchestration from CMake
- exposing built HDF5 headers and libraries to MPM-Code CMake targets
- preserving the existing PETSc path and Makefile build behavior

This design does not cover:

- embedding PETSc in the repository
- replacing the root Makefile
- adding VTKHDF output yet
- refactoring the existing solver build graph beyond what is needed for optional HDF5

## Why This Boundary

HDF5 is the new dependency required by a future VTKHDF path, and it is realistic to build from source inside CMake. PETSc is already a core solver dependency with MPI, BLAS/LAPACK, HYPRE, and configure-time complexity; folding PETSc into the same first step would expand both risk and build time too much.

This boundary gives us a low-risk path:

- keep current solver builds working
- add HDF5 only when requested
- prepare for future VTKHDF work without destabilizing PETSc integration

## Proposed User Experience

Default build behavior remains unchanged.

Examples:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake -j8
```

This should behave like today and not require HDF5.

To enable project-managed HDF5:

```bash
cmake -S . -B build-cmake -DMPM_ENABLE_HDF5=ON
cmake --build build-cmake -j8
```

Optional override for a non-default tarball:

```bash
cmake -S . -B build-cmake \
  -DMPM_ENABLE_HDF5=ON \
  -DMPM_HDF5_ARCHIVE=/abs/path/to/hdf5-hdf5_1.14.5.tar.gz
```

## Repository Layout

Add a new third-party source archive directory:

```text
Ext/
  hdf5-hdf5_1.14.5.tar.gz
```

Add a dedicated helper module for CMake dependency bootstrap logic:

```text
cmake/
  MPMBuildHDF5.cmake
```

No changes are required to the root `Makefile` for this phase.

## Build Architecture

### 1. Root CMake Option Layer

Add the following cache options near the top-level CMake configuration:

- `MPM_ENABLE_HDF5` default `OFF`
- `MPM_HDF5_ARCHIVE` default `${PROJECT_SOURCE_DIR}/Ext/hdf5-hdf5_1.14.5.tar.gz`
- `MPM_HDF5_ENABLE_PARALLEL` default `ON`
- `MPM_HDF5_BUILD_HL` default `OFF`
- `MPM_HDF5_BUILD_TOOLS` default `OFF`

Rationale:

- `OFF` by default keeps current users unblocked
- tarball path override avoids hard-coupling every environment to the bundled archive
- parallel HDF5 should be controllable because future VTKHDF may want MPI-enabled HDF5

### 2. HDF5 Bootstrap Helper

Create `cmake/MPMBuildHDF5.cmake` to encapsulate all HDF5 orchestration. The helper should:

- validate the tarball path early
- create a deterministic build area under the CMake binary tree
- unpack the tarball into that build area
- configure HDF5 in an out-of-source sub-build
- build and install HDF5 into a private prefix under the CMake binary tree
- export include path and library path variables back to the main project

Recommended install prefix pattern:

```text
<binary_dir>/third_party/hdf5/install
```

Recommended source/build pattern:

```text
<binary_dir>/third_party/hdf5/src
<binary_dir>/third_party/hdf5/build
```

### 3. Orchestration Strategy

Use an `add_custom_command()` + `add_custom_target()` pattern similar to Athene for the first implementation.

Why this over `FetchContent`:

- HDF5 is a heavy dependency with its own build options
- installing to a private prefix is easier to reason about with an explicit external configure/build/install sequence
- it matches the working reference style already proven in Athene

Why this over `ExternalProject_Add` for version one:

- fewer moving parts for a codebase that is not yet fully standardized on CMake
- easier to debug with explicit configure/build/install commands

## Main Project Integration

When `MPM_ENABLE_HDF5=OFF`:

- no HDF5 tarball checks
- no HDF5 build target
- no HDF5 include directories or libraries linked anywhere

When `MPM_ENABLE_HDF5=ON`:

- the main module target gains a dependency on the HDF5 build target
- include directories point to the private HDF5 install tree
- target link libraries include the produced HDF5 static library
- if parallel HDF5 is enabled, the HDF5 sub-build is configured consistently with the already discovered MPI toolchain

The existing PETSc logic in `module/CMakeLists.txt` should remain intact and independent.

## PETSc Policy

PETSc stays external.

Allowed sources:

- `PETSC_DIR` and `PETSC_ARCH` from environment
- existing hard-coded project-local path fallback

This design deliberately avoids:

- unpacking PETSc tarballs during CMake
- driving PETSc `configure` from inside the project
- changing current solver link assumptions

## Error Handling

Configuration must fail early and clearly when:

- `MPM_ENABLE_HDF5=ON` and the tarball does not exist
- MPI is required for parallel HDF5 but MPI was not found
- HDF5 configure/build/install fails

Suggested failure style:

- state which file/path is missing
- state which stage failed: unpack, configure, build, or install
- suggest the exact cache variable the user may override, namely `MPM_HDF5_ARCHIVE`

## Testing Strategy

### Configuration Tests

1. CMake configure with default options
   - expected: success
   - expected: no HDF5 build triggered

2. CMake configure with `MPM_ENABLE_HDF5=ON`
   - expected: tarball validation
   - expected: HDF5 bootstrap target appears

3. CMake configure with invalid `MPM_HDF5_ARCHIVE`
   - expected: immediate readable failure

### Build Tests

1. Build without HDF5
   - expected: current CMake path still links PETSc as before

2. Build with HDF5 enabled
   - expected: HDF5 builds under the binary tree
   - expected: HDF5 installs to the private prefix
   - expected: the main project still builds

### Future Readiness Check

After this work, adding a small HDF5 compile check target should be straightforward:

- include `hdf5.h`
- call a minimal symbol such as `H5open()`

That check is optional for this phase, but the design should leave room for it.

## Risks

### Build-Time Cost

HDF5 is a heavy dependency. First-time builds with `MPM_ENABLE_HDF5=ON` will be noticeably slower.

Mitigation:

- keep it optional
- install to a stable private prefix so rebuilds can reuse outputs

### CMake / Makefile Divergence

This adds capability only to the CMake path while the project still primarily uses the root Makefile.

Mitigation:

- treat this as an additive path
- avoid touching Makefile behavior in this phase
- document clearly that HDF5 bootstrap is CMake-only

### Parallel HDF5 Toolchain Drift

If HDF5 is built with a different MPI toolchain than the main project, runtime issues may appear later.

Mitigation:

- feed the already discovered MPI compiler settings into the HDF5 sub-build
- keep PETSc external and unchanged so only one new dependency is moving

## Recommendation

Implement the HDF5 bootstrap exactly as a CMake-only optional feature, using a bundled source tarball in `Ext/` and an Athene-style explicit external build sequence.

This gives MPM-Code a clean stepping stone toward future VTKHDF work without widening the change surface to PETSc or the legacy Makefile flow.
