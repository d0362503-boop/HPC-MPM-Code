# -----------------------------------------------------------------------------
# Data partition tools
# -----------------------------------------------------------------------------

add_library(mpm_divide_common OBJECT
    data_partitioner.cpp
)

target_link_libraries(mpm_divide_common PUBLIC mpm_modules)

# Keep all candidate subdirectories visible here.
# Uncomment exactly the cases you want to build.
# add_subdirectory(divide_fluid fluid)
# add_subdirectory(divide_solid solid)
add_subdirectory(divide_fsi fsi)
