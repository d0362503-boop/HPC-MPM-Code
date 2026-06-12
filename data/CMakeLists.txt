# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

add_executable(makinput_fluid fluid/makinput.cpp)
add_executable(makinput_solid solid/makinput.cpp)
add_executable(makinput_fsi fsi/makinput.cpp)

foreach(target makinput_fluid makinput_solid makinput_fsi)
    target_include_directories(${target} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/module
        ${CMAKE_SOURCE_DIR}/module/fluid/MPM
        ${CMAKE_SOURCE_DIR}/module/fluid/FEM
        ${CMAKE_SOURCE_DIR}/module/solid/implicit
    )

    target_link_libraries(${target} PRIVATE mpm_modules)
endforeach()
