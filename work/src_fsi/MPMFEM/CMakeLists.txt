# -----------------------------------------------------------------------------
# FSI Solver Source - Fluid-Structure Interaction coupling solver
# -----------------------------------------------------------------------------

# Automatically collect all .cpp files in this directory
file(GLOB SRC_FSI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_fsi OBJECT ${SRC_FSI_SOURCES})

# Own headers are public to consumers of this solver library.
target_include_directories(mpm_src_fsi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# mpm_modules provides the module include directories and all transitive
# dependencies (MPI, PETSc, HDF5, compiler compatibility).
target_link_libraries(mpm_src_fsi PUBLIC mpm_modules)
