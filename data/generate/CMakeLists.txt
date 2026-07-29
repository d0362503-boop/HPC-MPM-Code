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

# Keep all candidate subdirectories visible here.
# Uncomment exactly the cases you want to build.
# add_subdirectory(fluid)
# add_subdirectory(solid)
add_subdirectory(fsi)
