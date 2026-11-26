# -----------------------------------------------------------------------------
# Fluid Module - Fluid dynamics solvers
# -----------------------------------------------------------------------------

# Include all fluid method implementations; every driver executable in work/
# links them through mpm_modules.
add_subdirectory(FEM)
add_subdirectory(MPM)
message(STATUS "Including FEM and MPM fluid methods in module/fluid")

# Automatically collect all .cpp files in current directory
file(GLOB FLUID_BASE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# Add sources to parent module library
target_sources(mpm_modules PRIVATE ${FLUID_BASE_SOURCES})

# Add include directories to parent module
target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
