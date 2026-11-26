# -----------------------------------------------------------------------------
# Implicit Module - Time integration, implicit MPM kinematics and Newton driver
# -----------------------------------------------------------------------------

file(GLOB IMPLICIT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

target_sources(mpm_modules PRIVATE ${IMPLICIT_SOURCES})

target_include_directories(mpm_modules PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
