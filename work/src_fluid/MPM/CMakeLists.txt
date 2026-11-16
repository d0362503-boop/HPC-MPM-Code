# -----------------------------------------------------------------------------
# Fluid MPM driver source
# -----------------------------------------------------------------------------

file(GLOB SRC_FLUID_MPM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_fluid OBJECT ${SRC_FLUID_MPM_SOURCES})

# Own headers are public to consumers of this solver library.
target_include_directories(mpm_src_fluid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# mpm_modules provides the module include directories and all transitive
# dependencies (MPI, PETSc, HDF5, compiler compatibility).
target_link_libraries(mpm_src_fluid PUBLIC mpm_modules)
