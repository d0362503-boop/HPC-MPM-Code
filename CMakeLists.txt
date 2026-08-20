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

# Project headers use paths relative to the source root.
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------------------------------------------------
# Compiler compatibility helpers
# ----------------------------------------------------------------------------
# All targets that need C++17 and/or std::filesystem should link this
# interface library. It centralizes compiler-specific workarounds so they are
# maintained in exactly one place.
add_library(mpm_cxx_compat INTERFACE)

target_compile_features(mpm_cxx_compat INTERFACE cxx_std_17)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(mpm_cxx_compat INTERFACE stdc++fs)
endif()

# ----------------------------------------------------------------------------
# Optimization flags
# ----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(MPM_ENABLE_NATIVE_ARCH "Enable -march=native (disable for cross-architecture builds)" ON)
if(MPM_ENABLE_NATIVE_ARCH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# ----------------------------------------------------------------------------
# Build options
# ----------------------------------------------------------------------------
include(cmake/options.cmake)

# Only the 'build' directory is allowed as the CMake binary directory.
get_filename_component(MPM_BUILD_DIR_NAME ${CMAKE_BINARY_DIR} NAME)
if(NOT MPM_BUILD_DIR_NAME STREQUAL "build")
    message(FATAL_ERROR "CMake binary directory must be named 'build'. Current: ${MPM_BUILD_DIR_NAME}")
endif()

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

if(BUILD_PETSC)
    include(cmake/BuildPETSc.cmake)
    add_custom_target(petsc_lib_built)
endif()

# ----------------------------------------------------------------------------
# Subdirectories
# ----------------------------------------------------------------------------
add_subdirectory(module)
add_subdirectory(data)
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

if(TARGET mpm_src_fsi)
    target_link_libraries(MPM PRIVATE mpm_src_fsi)
endif()
if(TARGET mpm_src_fluid)
    target_link_libraries(MPM PRIVATE mpm_src_fluid)
endif()
if(TARGET mpm_src_solid)
    target_link_libraries(MPM PRIVATE mpm_src_solid)
endif()

# ----------------------------------------------------------------------------
# Configuration summary
# ----------------------------------------------------------------------------
# Determine active solver from the targets created in work/.
set(MPM_ACTIVE_SOLVER "unknown")
if(TARGET mpm_src_fsi)
    set(MPM_ACTIVE_SOLVER "FSI")
elseif(TARGET mpm_src_fluid)
    set(MPM_ACTIVE_SOLVER "FLUID")
elseif(TARGET mpm_src_solid)
    set(MPM_ACTIVE_SOLVER "SOLID")
endif()

# TODO: These will become cache options once Phase 4 (dependency strategy) is
# implemented. Until then, bundled tarballs are the only supported path.
set(MPM_PETSC_USE_SYSTEM OFF)
set(MPM_HDF5_USE_SYSTEM OFF)

message(STATUS "")
message(STATUS "MPM Configuration Summary:")
message(STATUS "  C++ compiler      : ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Native arch       : ${MPM_ENABLE_NATIVE_ARCH}")
message(STATUS "  Solver            : ${MPM_ACTIVE_SOLVER}")
message(STATUS "  MPI               : ${MPI_CXX_COMPILER} (USE_MPI=${USE_MPI})")
message(STATUS "  PETSc             : ${PETSC_DIR} (system=${MPM_PETSC_USE_SYSTEM})")
message(STATUS "  HDF5              : ${HDF5_INCLUDE_DIRS} (system=${MPM_HDF5_USE_SYSTEM})")
message(STATUS "  ZLIB              : ${ZLIB_VERSION}")
message(STATUS "")
