# -----------------------------------------------------------------------------
# Fluid FEM Module - Finite Element Method for fluid dynamics
# -----------------------------------------------------------------------------

file(GLOB FLUID_FEM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

target_sources(mpm_modules PRIVATE ${FLUID_FEM_SOURCES})

target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
