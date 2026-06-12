# -----------------------------------------------------------------------------
# Solid Solver Source - Main solver implementations for solid problems
# -----------------------------------------------------------------------------

# Select directory based on solid method
if(SOLID_METHOD STREQUAL "EXPLICIT")
    set(SOLID_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/explicit)
    message(STATUS "Using EXPLICIT solid solver sources")
else()
    set(SOLID_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/implicit)
    message(STATUS "Using IMPLICIT solid solver sources")
endif()

# Automatically collect all .cpp files in the selected directory
file(GLOB SRC_SOLID_SOURCES 
    ${SOLID_SRC_DIR}/*.cpp
)

add_library(mpm_src_solid OBJECT ${SRC_SOLID_SOURCES})

# Include directories - no relative paths needed
target_include_directories(mpm_src_solid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../module
)

# Link MPI if available
if(USE_MPI)
    target_link_libraries(mpm_src_solid PUBLIC MPI::MPI_CXX)
endif()

# Link mpm_modules to inherit PETSc and other dependencies
target_link_libraries(mpm_src_solid PUBLIC mpm_modules)
