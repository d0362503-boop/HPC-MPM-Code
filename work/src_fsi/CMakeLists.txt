# -----------------------------------------------------------------------------
# FSI Solver Source - Fluid-Structure Interaction coupling solver
# -----------------------------------------------------------------------------

# Automatically collect all .cpp files in this directory
file(GLOB SRC_FSI_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_fsi OBJECT ${SRC_FSI_SOURCES})

# Include directories - no relative paths needed
# All module headers can be included directly: #include "bc.h"
target_include_directories(mpm_src_fsi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../module
)

# Link MPI if available
if(USE_MPI)
    target_link_libraries(mpm_src_fsi PUBLIC MPI::MPI_CXX)
endif()

# Link mpm_modules to inherit PETSc and other dependencies
target_link_libraries(mpm_src_fsi PUBLIC mpm_modules)
