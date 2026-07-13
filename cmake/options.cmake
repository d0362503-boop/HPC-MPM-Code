# ----------------------------------------------------------------------------
# MPM-Code build options
# ----------------------------------------------------------------------------

# Solver source selection.
# The active solver source tree is selected by uncommenting exactly one
# add_subdirectory(...) line in work/CMakeLists.txt. Keep that file in sync
# with the function call that is uncommented in MPM_main.cpp.

# Data generator / partitioner selection.
# The active generator and partitioner cases are selected by uncommenting the
# relevant add_subdirectory(...) lines in data/generate/CMakeLists.txt and
# data/divide/CMakeLists.txt, respectively.

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
