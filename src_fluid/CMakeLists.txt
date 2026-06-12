# -----------------------------------------------------------------------------
# Fluid Solver Source - Main solver implementations for fluid problems
# -----------------------------------------------------------------------------

# Automatically collect all .cpp files in this directory and subdirectories
file(GLOB_RECURSE SRC_FLUID_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_fluid OBJECT ${SRC_FLUID_SOURCES})

# Include directories - no relative paths needed
# All module headers can be included directly: #include "bc.h"
target_include_directories(mpm_src_fluid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../module
)

# Link MPI if available
if(USE_MPI)
    target_link_libraries(mpm_src_fluid PUBLIC MPI::MPI_CXX)
endif()

# Link mpm_modules to inherit PETSc and other dependencies
target_link_libraries(mpm_src_fluid PUBLIC mpm_modules)
