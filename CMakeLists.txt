# ----------------------------------------------------------------------------
#  Root CMake file for MPM-Code
#
#  Build instructions:
#    mkdir build && cd build
#    cmake ..
#    make -j8
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20.0)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW) # DOWNLOAD_EXTRACT_TIMESTAMP option
endif()

project(
    MPM-Code
    VERSION 0.1.0
    DESCRIPTION "MPM/FEM/FSI solver with VTK HDF5 output support"
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------------------------------------------------
# Build options
# ----------------------------------------------------------------------------
option(USE_MPI "Enable MPI parallelism" ON)
option(USE_HDF5 "Enable HDF5 for VTK HDF5 output" ON)
option(HDF5_ENABLE_PARALLEL "Build parallel HDF5 (requires MPI)" ON)
option(HDF5_BUILD_TOOLS "Build HDF5 tools" OFF)
option(HDF5_BUILD_HL_LIB "Build HDF5 high-level library" OFF)

set(USE_SRC_FLUID false CACHE BOOL "Use work/src_fluid sources")
set(USE_SRC_SOLID false CACHE BOOL "Use work/src_solid sources")
set(USE_SRC_FSI true CACHE BOOL "Use work/src_fsi sources")
set(USE_SOLVER true CACHE BOOL "Build module/solver")
set(USE_FLUID true CACHE BOOL "Build module/fluid")
set(USE_SOLID true CACHE BOOL "Build module/solid")
set(USE_FEM true CACHE BOOL "Use fluid FEM mode")
set(USE_MPM false CACHE BOOL "Use fluid MPM mode")
set(SOLID_METHOD "IMPLICIT" CACHE STRING "Solid solver method: EXPLICIT or IMPLICIT")
set(FLUID_METHOD "FEM" CACHE STRING "Fluid solver method: FEM or MPM")

# ----------------------------------------------------------------------------
# MPI
# ----------------------------------------------------------------------------
if(USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
endif()

# ----------------------------------------------------------------------------
# Third-party dependencies (built from Ext tarballs into external/)
# ----------------------------------------------------------------------------
set(SRC_DIR ${CMAKE_SOURCE_DIR})
set(BIN_DIR ${CMAKE_BINARY_DIR})

if(USE_HDF5)
    find_package(ZLIB REQUIRED)
    include(cmake/BuildHDF5.cmake)
    add_custom_target(hdf5_lib_built DEPENDS hdf5_is_built)
endif()

include(cmake/BuildPETSc.cmake)
add_custom_target(petsc_lib_built)

# ----------------------------------------------------------------------------
# Subdirectories
# ----------------------------------------------------------------------------
add_subdirectory(module)
add_subdirectory(work)

# ----------------------------------------------------------------------------
# Main executable
# ----------------------------------------------------------------------------
add_executable(MPM MPM_main.cpp)

target_include_directories(MPM PRIVATE ${CMAKE_SOURCE_DIR}/module)

if(USE_MPI)
    target_link_libraries(MPM PRIVATE MPI::MPI_CXX)
endif()

target_link_libraries(MPM PRIVATE mpm_modules)

if(USE_SRC_FSI)
    target_link_libraries(MPM PRIVATE mpm_src_fsi)
endif()
if(USE_SRC_FLUID)
    target_link_libraries(MPM PRIVATE mpm_src_fluid)
endif()
if(USE_SRC_SOLID)
    target_link_libraries(MPM PRIVATE mpm_src_solid)
endif()
