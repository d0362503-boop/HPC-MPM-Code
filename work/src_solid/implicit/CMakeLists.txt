# -----------------------------------------------------------------------------
# Solid Implicit Solver Source
# -----------------------------------------------------------------------------

file(GLOB SRC_SOLID_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_solid OBJECT ${SRC_SOLID_SOURCES})

# Own headers are public to consumers of this solver library.
target_include_directories(mpm_src_solid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# mpm_modules provides the module include directories and all transitive
# dependencies (MPI, PETSc, HDF5, compiler compatibility).
target_link_libraries(mpm_src_solid PUBLIC mpm_modules)
