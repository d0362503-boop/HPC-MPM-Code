# HDF5 Auto-Build CMake Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional CMake-only HDF5 bootstrap path that builds HDF5 from a bundled archive, while keeping PETSc external and leaving the current Makefile build unchanged.

**Architecture:** First create a minimal top-level CMake entry point because the repository currently has no root `CMakeLists.txt`. Then add a dedicated `cmake/MPMBuildHDF5.cmake` helper that unpacks, configures, builds, and installs HDF5 into the binary tree, and wire the resulting include/library paths into the existing `module/` object-library build only when `MPM_ENABLE_HDF5=ON`.

**Tech Stack:** CMake 3.20+, MPI, external PETSc 3.24.5, optional HDF5 1.14.5 source archive, existing object-library CMake fragments under `module/` and `src_*`.

---

## File Structure

**Files:**
- Create: `CMakeLists.txt`
- Create: `cmake/MPMBuildHDF5.cmake`
- Create: `Ext/README.md`
- Modify: `module/CMakeLists.txt`
- Modify: `src_fluid/CMakeLists.txt`
- Modify: `src_solid/CMakeLists.txt`
- Modify: `src_fsi/CMakeLists.txt`
- Modify: `docs/superpowers/specs/2026-06-12-hdf5-autobuild-design.md`
- Verify against: `Makefile`

**Responsibilities:**
- `CMakeLists.txt`
  - top-level project entry
  - cache options
  - MPI discovery
  - solver-mode selection
  - HDF5 helper inclusion
  - executable assembly from existing object libraries
- `cmake/MPMBuildHDF5.cmake`
  - validate `MPM_HDF5_ARCHIVE`
  - stage unpack/build/install directories
  - define `mpm_hdf5_built` target
  - export HDF5 include/library variables back to the caller
- `module/CMakeLists.txt`
  - keep PETSc logic
  - consume optional HDF5 include/library variables only when enabled
- `src_* / CMakeLists.txt`
  - stay object-library based
  - rely on top-level variables instead of being standalone fragments
- `Ext/README.md`
  - document expected archive names and why archives are not used by `Makefile`

## Preconditions

- There is currently **no** root [CMakeLists.txt](<\\wsl.localhost\Ubuntu-22.04\home\pan\MPM-Code\CMakeLists.txt>), so the first task must create one.
- Existing CMake fragments already assume imported targets and top-level variables such as `USE_MPI` and `SOLID_METHOD`.
- The root `Makefile` remains the production build path during this phase and must not be broken or replaced.

### Task 1: Create the top-level CMake entry point

**Files:**
- Create: `CMakeLists.txt`
- Read: `module/CMakeLists.txt`
- Read: `src_fluid/CMakeLists.txt`
- Read: `src_solid/CMakeLists.txt`
- Read: `src_fsi/CMakeLists.txt`
- Read: `MPM_main.cpp`

- [ ] **Step 1: Write the failing configure command**

Run:

```bash
cmake -S . -B build-cmake
```

Expected: FAIL immediately because the repository has no root `CMakeLists.txt`.

- [ ] **Step 2: Create the minimal root CMake file**

Create `CMakeLists.txt` with a minimal driver that:

```cmake
cmake_minimum_required(VERSION 3.20)

project(MPMCode LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(USE_MPI "Enable MPI-aware build wiring" ON)
option(USE_SRC_FLUID "Build fluid solver objects" OFF)
option(USE_SRC_SOLID "Build solid solver objects" OFF)
option(USE_SRC_FSI "Build FSI solver objects" ON)
set(SOLID_METHOD "IMPLICIT" CACHE STRING "EXPLICIT or IMPLICIT")
set_property(CACHE SOLID_METHOD PROPERTY STRINGS EXPLICIT IMPLICIT)

if(USE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_subdirectory(module)

if(USE_SRC_FLUID)
  add_subdirectory(src_fluid)
endif()

if(USE_SRC_SOLID)
  add_subdirectory(src_solid)
endif()

if(USE_SRC_FSI)
  add_subdirectory(src_fsi)
endif()

add_executable(MPM
  MPM_main.cpp
  $<TARGET_OBJECTS:mpm_modules>
)

target_link_libraries(MPM PRIVATE mpm_modules)
```

- [ ] **Step 3: Extend the executable to include the selected source object libraries**

Update the top-level executable section so the correct object libraries are linked:

```cmake
target_sources(MPM PRIVATE
  $<TARGET_OBJECTS:mpm_modules>
  $<$<BOOL:${USE_SRC_FLUID}>:$<TARGET_OBJECTS:mpm_src_fluid>>
  $<$<BOOL:${USE_SRC_SOLID}>:$<TARGET_OBJECTS:mpm_src_solid>>
  $<$<BOOL:${USE_SRC_FSI}>:$<TARGET_OBJECTS:mpm_src_fsi>>
)

target_link_libraries(MPM PRIVATE mpm_modules)

if(USE_SRC_FLUID)
  target_link_libraries(MPM PRIVATE mpm_src_fluid)
endif()

if(USE_SRC_SOLID)
  target_link_libraries(MPM PRIVATE mpm_src_solid)
endif()

if(USE_SRC_FSI)
  target_link_libraries(MPM PRIVATE mpm_src_fsi)
endif()
```

- [ ] **Step 4: Run configure again**

Run:

```bash
cmake -S . -B build-cmake
```

Expected: configure now reaches the existing PETSc checks in `module/CMakeLists.txt` instead of failing at the repository root.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: add top-level cmake entry point"
```

### Task 2: Add optional HDF5 cache options and archive documentation

**Files:**
- Modify: `CMakeLists.txt`
- Create: `Ext/README.md`
- Test: configure cache output

- [ ] **Step 1: Add HDF5 cache options to the root CMake file**

Insert:

```cmake
option(MPM_ENABLE_HDF5 "Build and link project-managed HDF5" OFF)
set(MPM_HDF5_ARCHIVE
    "${CMAKE_CURRENT_SOURCE_DIR}/Ext/hdf5-hdf5_1.14.5.tar.gz"
    CACHE FILEPATH "Path to the HDF5 source archive used for project-managed builds")
option(MPM_HDF5_ENABLE_PARALLEL "Build HDF5 with MPI support" ON)
option(MPM_HDF5_BUILD_HL "Build HDF5 high-level library" OFF)
option(MPM_HDF5_BUILD_TOOLS "Build HDF5 tools" OFF)
```

- [ ] **Step 2: Document the expected archive in `Ext/README.md`**

Create:

```markdown
# Third-Party Source Archives

This directory is used only by the CMake build path.

Required for `-DMPM_ENABLE_HDF5=ON`:

- `hdf5-hdf5_1.14.5.tar.gz`

The root `Makefile` does not use archives from this directory.
PETSc remains an external dependency and is not bootstrapped from here in this phase.
```

- [ ] **Step 3: Verify the options appear in configure output**

Run:

```bash
cmake -S . -B build-cmake -L | rg "MPM_ENABLE_HDF5|MPM_HDF5_ARCHIVE|MPM_HDF5_ENABLE_PARALLEL"
```

Expected: all three variables are visible in the cache listing.

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt Ext/README.md
git commit -m "build: add hdf5 cache options"
```

### Task 3: Implement the HDF5 bootstrap helper

**Files:**
- Create: `cmake/MPMBuildHDF5.cmake`
- Modify: `CMakeLists.txt`
- Test: configure with an invalid archive path

- [ ] **Step 1: Create the helper skeleton**

Create `cmake/MPMBuildHDF5.cmake` with a function entry point:

```cmake
function(mpm_configure_hdf5_bootstrap)
  if(NOT MPM_ENABLE_HDF5)
    return()
  endif()
endfunction()
```

- [ ] **Step 2: Add archive validation**

Extend the function:

```cmake
if(NOT EXISTS "${MPM_HDF5_ARCHIVE}")
  message(FATAL_ERROR
    "MPM_ENABLE_HDF5=ON but HDF5 archive was not found at: ${MPM_HDF5_ARCHIVE}\n"
    "Set -DMPM_HDF5_ARCHIVE=/abs/path/to/hdf5-hdf5_1.14.5.tar.gz")
endif()
```

- [ ] **Step 3: Add the staged directory layout and custom build target**

Continue the function with:

```cmake
set(MPM_HDF5_STAGE_DIR   "${CMAKE_BINARY_DIR}/third_party/hdf5")
set(MPM_HDF5_SOURCE_DIR  "${MPM_HDF5_STAGE_DIR}/src")
set(MPM_HDF5_BUILD_DIR   "${MPM_HDF5_STAGE_DIR}/build")
set(MPM_HDF5_INSTALL_DIR "${MPM_HDF5_STAGE_DIR}/install")
set(MPM_HDF5_SENTINEL    "${MPM_HDF5_STAGE_DIR}/hdf5_built.stamp")

file(MAKE_DIRECTORY "${MPM_HDF5_STAGE_DIR}")
```

Add the external build rule:

```cmake
add_custom_command(
  OUTPUT "${MPM_HDF5_SENTINEL}"
  COMMAND ${CMAKE_COMMAND} -E rm -rf "${MPM_HDF5_SOURCE_DIR}" "${MPM_HDF5_BUILD_DIR}" "${MPM_HDF5_INSTALL_DIR}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${MPM_HDF5_STAGE_DIR}"
  COMMAND ${CMAKE_COMMAND} -E tar xzf "${MPM_HDF5_ARCHIVE}"
  WORKING_DIRECTORY "${MPM_HDF5_STAGE_DIR}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${MPM_HDF5_BUILD_DIR}"
  COMMAND ${CMAKE_COMMAND}
          -S "${MPM_HDF5_STAGE_DIR}/hdf5-hdf5_1.14.5"
          -B "${MPM_HDF5_BUILD_DIR}"
          -DBUILD_SHARED_LIBS=OFF
          -DHDF5_BUILD_TESTING=OFF
          -DHDF5_BUILD_EXAMPLES=OFF
          -DHDF5_BUILD_UTILS=${MPM_HDF5_BUILD_TOOLS}
          -DHDF5_BUILD_TOOLS=${MPM_HDF5_BUILD_TOOLS}
          -DHDF5_BUILD_HL_LIB=${MPM_HDF5_BUILD_HL}
          -DHDF5_ENABLE_PARALLEL=${MPM_HDF5_ENABLE_PARALLEL}
          -DCMAKE_INSTALL_PREFIX="${MPM_HDF5_INSTALL_DIR}"
  COMMAND ${CMAKE_COMMAND} --build "${MPM_HDF5_BUILD_DIR}" -j4
  COMMAND ${CMAKE_COMMAND} --install "${MPM_HDF5_BUILD_DIR}"
  COMMAND ${CMAKE_COMMAND} -E touch "${MPM_HDF5_SENTINEL}"
  COMMENT "Building project-managed HDF5"
)

add_custom_target(mpm_hdf5_built DEPENDS "${MPM_HDF5_SENTINEL}")
```

- [ ] **Step 4: Export include and library variables**

Append:

```cmake
set(MPM_HDF5_INCLUDE_DIR "${MPM_HDF5_INSTALL_DIR}/include" PARENT_SCOPE)
set(MPM_HDF5_LIBRARY     "${MPM_HDF5_INSTALL_DIR}/lib/libhdf5.a" PARENT_SCOPE)
set(MPM_HDF5_BUILT_TARGET mpm_hdf5_built PARENT_SCOPE)
```

- [ ] **Step 5: Call the helper from the root CMake file**

Insert in `CMakeLists.txt` before `add_subdirectory(module)`:

```cmake
include(MPMBuildHDF5)
mpm_configure_hdf5_bootstrap()
```

- [ ] **Step 6: Verify the invalid-archive failure**

Run:

```bash
cmake -S . -B build-cmake \
  -DMPM_ENABLE_HDF5=ON \
  -DMPM_HDF5_ARCHIVE=/tmp/does-not-exist.tar.gz
```

Expected: FAIL with a clear message that names `MPM_HDF5_ARCHIVE`.

- [ ] **Step 7: Commit**

```bash
git add CMakeLists.txt cmake/MPMBuildHDF5.cmake
git commit -m "build: add hdf5 bootstrap helper"
```

### Task 4: Wire optional HDF5 into the module build

**Files:**
- Modify: `module/CMakeLists.txt`
- Test: `cmake --build`

- [ ] **Step 1: Add conditional HDF5 integration to `module/CMakeLists.txt`**

Insert after PETSc setup:

```cmake
if(MPM_ENABLE_HDF5)
  add_dependencies(mpm_modules ${MPM_HDF5_BUILT_TARGET})
  target_include_directories(mpm_modules SYSTEM PUBLIC ${MPM_HDF5_INCLUDE_DIR})
  target_link_libraries(mpm_modules PUBLIC ${MPM_HDF5_LIBRARY} dl)
endif()
```

- [ ] **Step 2: Keep PETSc logic untouched**

Do not change:

```cmake
if(DEFINED ENV{PETSC_DIR})
  ...
else()
  set(PETSC_DIR "/home/pan/petsc-3.24.5")
  ...
endif()
```

Expected: PETSc source and link behavior remain identical to the current CMake fragment.

- [ ] **Step 3: Run a normal configure/build without HDF5**

Run:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake -j8
```

Expected: success with no HDF5 bootstrap target invoked.

- [ ] **Step 4: Run a configure/build with HDF5 enabled**

Run:

```bash
cmake -S . -B build-cmake -DMPM_ENABLE_HDF5=ON
cmake --build build-cmake -j8
```

Expected:

- HDF5 is unpacked under `build-cmake/third_party/hdf5`
- HDF5 installs under `build-cmake/third_party/hdf5/install`
- the `MPM` executable still links successfully

- [ ] **Step 5: Commit**

```bash
git add module/CMakeLists.txt
git commit -m "build: wire optional hdf5 into module target"
```

### Task 5: Stabilize the source-fragment assumptions

**Files:**
- Modify: `src_fluid/CMakeLists.txt`
- Modify: `src_solid/CMakeLists.txt`
- Modify: `src_fsi/CMakeLists.txt`

- [ ] **Step 1: Ensure all source fragment CMake files rely on the top-level driver**

Confirm these files only assume:

- `USE_MPI`
- `SOLID_METHOD`
- `mpm_modules`
- imported MPI target if enabled

If they contain standalone-project assumptions, remove them minimally.

- [ ] **Step 2: Verify each supported mode can still configure**

Run:

```bash
cmake -S . -B build-cmake-fsi    -DUSE_SRC_FSI=ON   -DUSE_SRC_FLUID=OFF -DUSE_SRC_SOLID=OFF
cmake -S . -B build-cmake-fluid  -DUSE_SRC_FSI=OFF  -DUSE_SRC_FLUID=ON  -DUSE_SRC_SOLID=OFF
cmake -S . -B build-cmake-solid  -DUSE_SRC_FSI=OFF  -DUSE_SRC_FLUID=OFF -DUSE_SRC_SOLID=ON -DSOLID_METHOD=IMPLICIT
```

Expected: each configure path succeeds at least through target generation.

- [ ] **Step 3: Commit**

```bash
git add src_fluid/CMakeLists.txt src_solid/CMakeLists.txt src_fsi/CMakeLists.txt
git commit -m "build: stabilize cmake source fragments"
```

### Task 6: Smoke-test HDF5 availability for future VTKHDF work

**Files:**
- Create: `cmake/hdf5_smoke.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Add a minimal compile-only smoke source**

Create:

```cpp
#include <hdf5.h>

int main() {
  return H5open();
}
```

- [ ] **Step 2: Add an opt-in smoke target**

Append to `CMakeLists.txt`:

```cmake
if(MPM_ENABLE_HDF5)
  add_executable(mpm_hdf5_smoke cmake/hdf5_smoke.cpp)
  add_dependencies(mpm_hdf5_smoke ${MPM_HDF5_BUILT_TARGET})
  target_include_directories(mpm_hdf5_smoke PRIVATE ${MPM_HDF5_INCLUDE_DIR})
  target_link_libraries(mpm_hdf5_smoke PRIVATE ${MPM_HDF5_LIBRARY} dl)
endif()
```

- [ ] **Step 3: Build the smoke target**

Run:

```bash
cmake -S . -B build-cmake -DMPM_ENABLE_HDF5=ON
cmake --build build-cmake --target mpm_hdf5_smoke -j8
```

Expected: the smoke target compiles and links successfully.

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt cmake/hdf5_smoke.cpp
git commit -m "build: add hdf5 smoke target"
```

### Task 7: Final verification and documentation pass

**Files:**
- Modify: `docs/superpowers/specs/2026-06-12-hdf5-autobuild-design.md`
- Verify: `Makefile`

- [ ] **Step 1: Re-run the two primary verification paths**

Run:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake -j8
cmake -S . -B build-cmake-hdf5 -DMPM_ENABLE_HDF5=ON
cmake --build build-cmake-hdf5 -j8
```

Expected:

- plain CMake build still works
- HDF5-enabled CMake build still works

- [ ] **Step 2: Confirm the Makefile path is untouched**

Run:

```bash
git diff -- Makefile
```

Expected: no changes.

- [ ] **Step 3: Update the spec if any final variable names or paths drifted**

Expected changes, if needed:

- `MPM_HDF5_ARCHIVE`
- actual HDF5 staging directory
- actual smoke-target name

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/specs/2026-06-12-hdf5-autobuild-design.md
git commit -m "docs: finalize hdf5 autobuild spec"
```
