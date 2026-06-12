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
# 优先使用环境变量 PETSC_DIR 指定的版本，否则使用项目本地 3.24.5 版本
if(DEFINED ENV{PETSC_DIR})
    set(PETSC_DIR "$ENV{PETSC_DIR}")
    set(PETSC_ARCH "$ENV{PETSC_ARCH}")
    message(STATUS "Using PETSc from environment: ${PETSC_DIR}")
else()
    set(PETSC_DIR "/home/pan/petsc-3.24.5")
    set(PETSC_ARCH "arch-linux-c-release")
    if(EXISTS "${PETSC_DIR}/include/petsc.h")
        message(STATUS "Using project-local PETSc: ${PETSC_DIR}")
    else()
        message(FATAL_ERROR "PETSc not found at ${PETSC_DIR}. Please set PETSC_DIR environment variable.")
    endif()
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

# -----------------------------------------------------------------------------
# Add Subdirectories
# -----------------------------------------------------------------------------
add_subdirectory(solver)
add_subdirectory(fluid)
add_subdirectory(solid)
