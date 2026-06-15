# -----------------------------------------------------------------------------
# Solid Explicit Solver Source
# -----------------------------------------------------------------------------

file(GLOB SRC_SOLID_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(mpm_src_solid OBJECT ${SRC_SOLID_SOURCES})

target_include_directories(mpm_src_solid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../../module
)

if(USE_MPI)
    target_link_libraries(mpm_src_solid PUBLIC MPI::MPI_CXX)
endif()

target_link_libraries(mpm_src_solid PUBLIC mpm_modules)
