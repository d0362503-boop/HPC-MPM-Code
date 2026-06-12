# -----------------------------------------------------------------------------
# Solid Implicit Module - Implicit time integration for solid mechanics
# -----------------------------------------------------------------------------

file(GLOB SOLID_IMPLICIT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

target_sources(mpm_modules PRIVATE ${SOLID_IMPLICIT_SOURCES})

target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
