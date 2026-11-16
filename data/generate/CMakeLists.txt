# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

include(${CMAKE_SOURCE_DIR}/cmake/MPMTools.cmake)

add_library(mpm_generate_common OBJECT
    data_generator.cpp
    output_util.cpp
)

# mpm_generate_common reuses module utilities (data_io, mesh, global dataset)
# and therefore inherits the module's public interface. The compatibility
# layer is linked explicitly so the dependency is visible even if the
# relationship to mpm_modules changes later.
target_link_libraries(mpm_generate_common
    PUBLIC mpm_modules mpm_cxx_compat
)

# Build the standalone fluid and solid generator modules independently.
add_subdirectory(fluid)
add_subdirectory(solid)

# The FSI generator remains independently selected by its own subdirectory.
add_subdirectory(fsi)
