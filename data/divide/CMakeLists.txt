# -----------------------------------------------------------------------------
# Data partition tools
# -----------------------------------------------------------------------------

# MPMTools.cmake is already included by data/generate, but include it here as
# well so this directory can be built independently.
include(${CMAKE_SOURCE_DIR}/cmake/MPMTools.cmake)

add_library(mpm_divide_common OBJECT
    data_partitioner.cpp
)

# mpm_divide_common reuses module utilities and inherits the module's public
# interface. The compatibility layer is linked explicitly so the dependency is
# visible even if the relationship to mpm_modules changes later.
target_link_libraries(mpm_divide_common
    PUBLIC mpm_modules mpm_cxx_compat
)

# Keep all candidate subdirectories visible here.
# Uncomment exactly the cases you want to build.
# add_subdirectory(divide_fluid fluid)
# add_subdirectory(divide_solid solid)
add_subdirectory(divide_fsi fsi)
