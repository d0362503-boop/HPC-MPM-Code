# -----------------------------------------------------------------------------
# Build PETSc from Ext/petsc-3.24.5.tar.gz into external/petsc
# This script is idempotent: if external/petsc already exists, it does nothing.
# -----------------------------------------------------------------------------

set(PETSC_TARBALL "${SRC_DIR}/Ext/petsc-3.24.5.tar.gz")
set(PETSC_INSTALL_DIR "${SRC_DIR}/external/petsc")
set(MARKER "${BIN_DIR}/petsc_is_built")

if(EXISTS "${PETSC_INSTALL_DIR}/include/petsc.h")
    message(STATUS "PETSc already installed at ${PETSC_INSTALL_DIR}; skipping build.")
    file(TOUCH "${MARKER}")
    return()
endif()

if(NOT EXISTS ${PETSC_TARBALL})
    message(FATAL_ERROR "PETSc tarball not found: ${PETSC_TARBALL}")
endif()

set(HYPRE_TARBALL "${SRC_DIR}/Ext/hypre-3.0.0.tar.gz")
set(PETSC_HYPRE_ARG "")
if(EXISTS ${HYPRE_TARBALL})
    set(PETSC_HYPRE_ARG "--download-hypre=${HYPRE_TARBALL}")
    message(STATUS "PETSc will build Hypre from: ${HYPRE_TARBALL}")
endif()

find_package(MPI REQUIRED)

set(PETSC_BUILD_DIR "${BIN_DIR}/petsc")
set(PETSC_SRC_DIR "${PETSC_BUILD_DIR}/src")

message(STATUS "Building PETSc from ${PETSC_TARBALL} ...")

file(MAKE_DIRECTORY ${PETSC_SRC_DIR})
file(REMOVE_RECURSE ${PETSC_SRC_DIR})
file(MAKE_DIRECTORY ${PETSC_SRC_DIR})

execute_process(
    COMMAND tar zxvf ${PETSC_TARBALL} -C ${PETSC_SRC_DIR} --strip-components=1
    RESULT_VARIABLE tar_result)
if(NOT tar_result EQUAL 0)
    message(FATAL_ERROR "Failed to extract PETSc tarball")
endif()

execute_process(
    COMMAND
        ${CMAKE_COMMAND} -E env PETSC_DIR=${PETSC_SRC_DIR}
        ./configure
        --prefix=${PETSC_INSTALL_DIR}
        --with-petsc-arch=arch-linux-c-release
        --with-cc=mpicc
        --with-cxx=mpicxx
        --with-fc=mpifort
        --download-f2cblaslapack=1
        ${PETSC_HYPRE_ARG}
        --with-debugging=0
        COPTFLAGS=-O3
        CXXOPTFLAGS=-O3
        FOPTFLAGS=-O3
    WORKING_DIRECTORY ${PETSC_SRC_DIR}
    RESULT_VARIABLE configure_result)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "PETSc configure failed")
endif()

execute_process(
    COMMAND
        ${CMAKE_COMMAND} -E env PETSC_DIR=${PETSC_SRC_DIR} PETSC_ARCH=arch-linux-c-release
        make all
    WORKING_DIRECTORY ${PETSC_SRC_DIR}
    RESULT_VARIABLE build_result)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "PETSc build failed")
endif()

execute_process(
    COMMAND
        ${CMAKE_COMMAND} -E env PETSC_DIR=${PETSC_SRC_DIR} PETSC_ARCH=arch-linux-c-release
        make install
    WORKING_DIRECTORY ${PETSC_SRC_DIR}
    RESULT_VARIABLE install_result)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "PETSc install failed")
endif()

file(TOUCH "${MARKER}")
message(STATUS "PETSc installed at ${PETSC_INSTALL_DIR}")
