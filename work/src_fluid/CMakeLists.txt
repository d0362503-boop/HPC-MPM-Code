# -----------------------------------------------------------------------------
# Fluid Solver Source - Main solver implementations for fluid problems
# -----------------------------------------------------------------------------

# Automatically collect all .cpp files in this directory and subdirectories
file(GLOB_RECURSE SRC_FLUID_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_fluid OBJECT ${SRC_FLUID_SOURCES})

# Own headers are public to consumers of this solver library.
target_include_directories(mpm_src_fluid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# mpm_modules provides the module include directories and all transitive
# dependencies (MPI, PETSc, HDF5, compiler compatibility).
target_link_libraries(mpm_src_fluid PUBLIC mpm_modules)
