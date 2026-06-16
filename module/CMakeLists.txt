# -----------------------------------------------------------------------------
# Module Library - Core MPM modules
# -----------------------------------------------------------------------------

# Automatically collect all .cpp files in current directory
file(GLOB MODULE_BASE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# Create module library (object library for combining with other objects)
add_library(mpm_modules OBJECT ${MODULE_BASE_SOURCES})

# -----------------------------------------------------------------------------
# Set include directories
# -----------------------------------------------------------------------------
target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link MPI if available
if(USE_MPI)
    target_link_libraries(mpm_modules PUBLIC MPI::MPI_CXX)
endif()

# -----------------------------------------------------------------------------
# Find PETSc
# -----------------------------------------------------------------------------
# PETSC_DIR is expected to be set by cmake/options.cmake (default: external/petsc).
# If it is empty, fall back to the PETSC_DIR environment variable.
if(NOT PETSC_DIR OR PETSC_DIR STREQUAL "")
    if(DEFINED ENV{PETSC_DIR})
        set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "PETSc installation directory" FORCE)
    else()
        message(FATAL_ERROR "PETSc not found. Set PETSC_DIR or enable BUILD_PETSC.")
    endif()
endif()

# Only adopt PETSC_ARCH from the environment when PETSc is not already
# installed in the flat prefix layout (${PETSC_DIR}/include, ${PETSC_DIR}/lib).
if(NOT PETSC_ARCH AND DEFINED ENV{PETSC_ARCH})
    if(NOT EXISTS "${PETSC_DIR}/include/petsc.h")
        set(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "PETSc architecture subdirectory" FORCE)
    endif()
endif()

# If PETSc is installed directly under the prefix, ignore any leftover
# PETSC_ARCH from a source-tree build to avoid noisy warnings.
if(EXISTS "${PETSC_DIR}/include/petsc.h")
    set(PETSC_ARCH "" CACHE STRING "PETSc architecture subdirectory" FORCE)
elseif(PETSC_ARCH AND NOT EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petsc.h")
    message(WARNING "PETSC_ARCH '${PETSC_ARCH}' not found under ${PETSC_DIR}; ignoring it.")
    set(PETSC_ARCH "" CACHE STRING "PETSc architecture subdirectory" FORCE)
endif()

if(NOT EXISTS "${PETSC_DIR}/include/petsc.h")
    message(FATAL_ERROR "PETSc header not found at ${PETSC_DIR}/include/petsc.h")
endif()

message(STATUS "Using PETSc: ${PETSC_DIR}")

# PETSc include directories and libraries
target_include_directories(mpm_modules SYSTEM PUBLIC ${PETSC_DIR}/include)
if(PETSC_ARCH)
    target_include_directories(mpm_modules SYSTEM PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/include)
    target_link_directories(mpm_modules PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/lib)
else()
    target_link_directories(mpm_modules PUBLIC ${PETSC_DIR}/lib)
endif()
target_link_libraries(mpm_modules PUBLIC petsc)

if(TARGET petsc_lib_built)
    add_dependencies(mpm_modules petsc_lib_built)
endif()

# -----------------------------------------------------------------------------
# HDF5
# -----------------------------------------------------------------------------
if(USE_HDF5)
    add_dependencies(mpm_modules hdf5_lib_built)
    target_include_directories(mpm_modules SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(mpm_modules PUBLIC ${HDF5_LIBRARIES} ZLIB::ZLIB dl)
    target_compile_definitions(mpm_modules PUBLIC HAVE_HDF5)
    message(STATUS "HDF5 include: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 library: ${HDF5_LIBRARIES}")
endif()

# -----------------------------------------------------------------------------
# Add Subdirectories
# -----------------------------------------------------------------------------
add_subdirectory(solver)
add_subdirectory(fluid)
add_subdirectory(solid)
