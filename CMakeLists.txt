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

if(USE_SRC_FSI)
    target_link_libraries(MPM PRIVATE mpm_src_fsi)
endif()
if(USE_SRC_FLUID)
    target_link_libraries(MPM PRIVATE mpm_src_fluid)
endif()
if(USE_SRC_SOLID)
    target_link_libraries(MPM PRIVATE mpm_src_solid)
endif()
