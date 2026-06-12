# -----------------------------------------------------------------------------
# Fluid MPM Module - Material Point Method for fluid dynamics
# -----------------------------------------------------------------------------

file(GLOB FLUID_MPM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

target_sources(mpm_modules PRIVATE ${FLUID_MPM_SOURCES})

target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
