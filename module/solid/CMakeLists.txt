# -----------------------------------------------------------------------------
# Solid Module - Solid mechanics solvers
# -----------------------------------------------------------------------------

# Include all solid method implementations. The active driver is selected in
# work/CMakeLists.txt together with the matching call in MPM_main.cpp.
add_subdirectory(explicit)
add_subdirectory(implicit)
message(STATUS "Including EXPLICIT and IMPLICIT solid methods in module/solid")

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
