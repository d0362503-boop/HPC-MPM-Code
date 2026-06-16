# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

add_library(mpm_generate_common OBJECT
    data_generator.cpp
    output_util.cpp
)

target_link_libraries(mpm_generate_common PUBLIC mpm_modules)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(mpm_generate_common INTERFACE stdc++fs)
endif()

# Keep all candidate subdirectories visible here.
# Uncomment exactly the cases you want to build.
# add_subdirectory(fluid)
# add_subdirectory(solid)
add_subdirectory(fsi)
