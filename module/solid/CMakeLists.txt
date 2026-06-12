# -----------------------------------------------------------------------------
# Solid Module - Solid mechanics solvers
# -----------------------------------------------------------------------------

# Select subdirectory based on solid method
if(SOLID_METHOD STREQUAL "EXPLICIT")
    set(SOLID_SUBDIR explicit)
    message(STATUS "Using EXPLICIT solid method in module/solid")
else()
    set(SOLID_SUBDIR implicit)
    message(STATUS "Using IMPLICIT solid method in module/solid")
endif()

# Add selected subdirectory
add_subdirectory(${SOLID_SUBDIR})

# Automatically collect all .cpp files in current directory (excluding subdirectories)
file(GLOB SOLID_BASE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# Add sources to parent module library
target_sources(mpm_modules PRIVATE ${SOLID_BASE_SOURCES})

# Add include directories to parent module
target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
