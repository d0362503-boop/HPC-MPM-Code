# -----------------------------------------------------------------------------
# Solver Module - Linear system solvers and sparse matrix operations
# -----------------------------------------------------------------------------

file(GLOB SOLVER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

target_sources(mpm_modules PRIVATE ${SOLVER_SOURCES})

target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
