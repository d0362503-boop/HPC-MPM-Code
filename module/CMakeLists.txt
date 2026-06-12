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
# Prefer project-local PETSc; fall back to PETSC_DIR environment variable.
set(MPM_LOCAL_PETSC "${CMAKE_SOURCE_DIR}/external/petsc")
if(EXISTS "${MPM_LOCAL_PETSC}/include/petsc.h")
    set(PETSC_DIR "${MPM_LOCAL_PETSC}")
    message(STATUS "Using project-local PETSc: ${PETSC_DIR}")
elseif(DEFINED ENV{PETSC_DIR})
    set(PETSC_DIR "$ENV{PETSC_DIR}")
    set(PETSC_ARCH "$ENV{PETSC_ARCH}")
    message(STATUS "Using PETSc from environment: ${PETSC_DIR}")
else()
    message(FATAL_ERROR "PETSc not found. Please install it to ${MPM_LOCAL_PETSC} or set PETSC_DIR.")
endif()

# 设置 PETSc 包含目录和库
target_include_directories(mpm_modules SYSTEM PUBLIC ${PETSC_DIR}/include)
if(PETSC_ARCH)
    target_include_directories(mpm_modules SYSTEM PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/include)
    target_link_directories(mpm_modules PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/lib)
    # PETSc 可能通过 --prefix 安装到独立目录
    if(EXISTS "${PETSC_DIR}-install/lib")
        target_link_directories(mpm_modules PUBLIC ${PETSC_DIR}-install/lib)
    endif()
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
