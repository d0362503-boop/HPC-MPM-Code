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

set(USE_SRC_FLUID false CACHE BOOL "Use src_fluid sources")
set(USE_SRC_SOLID false CACHE BOOL "Use src_solid sources")
set(USE_SRC_FSI true CACHE BOOL "Use src_fsi sources")
set(USE_SOLVER true CACHE BOOL "Build module/solver")
set(USE_FLUID true CACHE BOOL "Build module/fluid")
set(USE_SOLID true CACHE BOOL "Build module/solid")
set(USE_FEM true CACHE BOOL "Use fluid FEM mode")
set(USE_MPM false CACHE BOOL "Use fluid MPM mode")
set(SOLID_METHOD "IMPLICIT" CACHE STRING "Solid solver method: EXPLICIT or IMPLICIT")

# ----------------------------------------------------------------------------
# MPI
# ----------------------------------------------------------------------------
if(USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
endif()

# ----------------------------------------------------------------------------
# HDF5 (built from Ext tarball, Athene-style)
# ----------------------------------------------------------------------------
if(USE_HDF5)
    set(HDF5_TARBALL "${CMAKE_SOURCE_DIR}/Ext/hdf5-hdf5_1.14.5.tar.gz")
    if(NOT EXISTS ${HDF5_TARBALL})
        message(FATAL_ERROR "HDF5 tarball not found: ${HDF5_TARBALL}")
    endif()

    if(HDF5_ENABLE_PARALLEL AND NOT USE_MPI)
        message(FATAL_ERROR "HDF5_ENABLE_PARALLEL requires USE_MPI=ON")
    endif()

    # HDF5 needs ZLIB. Prefer a system ZLIB; if not found, fall back to
    # building one from Ext (not implemented here).
    find_package(ZLIB REQUIRED)

    set(TEMP_HDF5_DIR "${CMAKE_BINARY_DIR}/hdf5")
    set(HDF5_INSTALL_DIR "${TEMP_HDF5_DIR}/hdf5")
    set(HDF5_INCLUDE_DIRS "${HDF5_INSTALL_DIR}/include"
        CACHE STRING "HDF5 include directories" FORCE)
    set(HDF5_LIBRARIES "${HDF5_INSTALL_DIR}/lib/libhdf5.a"
        CACHE STRING "HDF5 library" FORCE)

    add_custom_command(
        OUTPUT hdf5_is_built
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TEMP_HDF5_DIR}
        COMMAND tar zxvf ${HDF5_TARBALL} -C ${TEMP_HDF5_DIR}/..
        COMMAND ${CMAKE_COMMAND} -E make_directory ${TEMP_HDF5_DIR}/build
        COMMAND
            ${CMAKE_COMMAND} -S ${TEMP_HDF5_DIR} -B ${TEMP_HDF5_DIR}/build
            -DHDF5_ENABLE_SZIP_SUPPORT=OFF
            -DHDF5_BUILD_UTILS=OFF
            -DBUILD_SHARED_LIBS=OFF
            -DHDF5_BUILD_TESTING=OFF
            -DHDF5_BUILD_EXAMPLES=OFF
            -DHDF5_ENABLE_PARALLEL=${HDF5_ENABLE_PARALLEL}
            -DHDF5_BUILD_TOOLS=${HDF5_BUILD_TOOLS}
            -DHDF5_BUILD_HL_LIB=${HDF5_BUILD_HL_LIB}
            -DHDF5_BUILD_CPP_LIB=OFF
            -DH5EX_BUILD_EXAMPLES=OFF
            -DH5EX_BUILD_TESTING=OFF
            -DH5EX_BUILD_HL_LIB=OFF
            -DTEST_SHELL_SCRIPTS=OFF
            -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        COMMAND ${CMAKE_COMMAND} --build ${TEMP_HDF5_DIR}/build
        COMMAND ${CMAKE_COMMAND} --install ${TEMP_HDF5_DIR}/build
        COMMAND ${CMAKE_COMMAND} -E touch hdf5_is_built
        COMMENT "Building HDF5 from Ext/hdf5-hdf5_1.14.5.tar.gz ...")

    add_custom_target(hdf5_lib_built DEPENDS hdf5_is_built)

    message(STATUS "HDF5 will be built from: ${HDF5_TARBALL}")
    message(STATUS "HDF5 install directory: ${HDF5_INSTALL_DIR}")
endif()

# ----------------------------------------------------------------------------
# Subdirectories
# ----------------------------------------------------------------------------
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
