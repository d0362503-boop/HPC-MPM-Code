# -----------------------------------------------------------------------------
# Fluid Module - Fluid dynamics solvers
# -----------------------------------------------------------------------------

# Select subdirectory based on fluid method
if(FLUID_METHOD STREQUAL "MPM")
    set(FLUID_SUBDIR MPM)
    message(STATUS "Using MPM fluid method in module/fluid")
else()
    set(FLUID_SUBDIR FEM)
    message(STATUS "Using FEM fluid method in module/fluid")
endif()

# Add selected subdirectory
add_subdirectory(${FLUID_SUBDIR})

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
