# ----------------------------------------------------------------------------
# MPM-Code build options
# ----------------------------------------------------------------------------

# Solver source selection.
# Set exactly one of these to true to choose which driver code is compiled.
set(USE_SRC_FSI true CACHE BOOL "Build work/src_fsi sources")
set(USE_SRC_FLUID false CACHE BOOL "Build work/src_fluid sources")
set(USE_SRC_SOLID false CACHE BOOL "Build work/src_solid sources")

# Data generator selection.
option(USE_DATA_FLUID "Build makinput_fluid data generator" OFF)
option(USE_DATA_SOLID "Build makinput_solid data generator" OFF)
option(USE_DATA_FSI "Build makinput_fsi data generator" OFF)

# Solver method selection
set(FLUID_METHOD "FEM" CACHE STRING "Fluid solver method: FEM or MPM")
set_property(CACHE FLUID_METHOD PROPERTY STRINGS "FEM" "MPM")

set(SOLID_METHOD "IMPLICIT" CACHE STRING "Solid solver method: EXPLICIT or IMPLICIT")
set_property(CACHE SOLID_METHOD PROPERTY STRINGS "EXPLICIT" "IMPLICIT")

# Module toggles
option(USE_SOLVER "Build module/solver" ON)
option(USE_FLUID "Build module/fluid" ON)
option(USE_SOLID "Build module/solid" ON)

# PETSc options
option(BUILD_PETSC "Build PETSc from bundled Ext tarball" ON)
set(PETSC_DIR "${CMAKE_SOURCE_DIR}/external/petsc" CACHE PATH "PETSc installation directory")
set(PETSC_ARCH "" CACHE STRING "PETSc architecture subdirectory (leave empty for prefix installs)")

# HDF5 options
option(USE_HDF5 "Enable HDF5 for VTK HDF5 output" ON)
option(HDF5_ENABLE_PARALLEL "Build parallel HDF5 (requires MPI)" ON)
option(HDF5_BUILD_TOOLS "Build HDF5 tools" OFF)
option(HDF5_BUILD_HL_LIB "Build HDF5 high-level library" OFF)

# MPI option
option(USE_MPI "Enable MPI parallelism" ON)
