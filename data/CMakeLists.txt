# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

set(DATA_COMMON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/data_main.cpp
)

function(add_makinput_target case_name)
    set(target_name "makinput_${case_name}")
    set(case_source "${CMAKE_CURRENT_SOURCE_DIR}/${case_name}/makinput.cpp")

    if(NOT EXISTS ${case_source})
        message(WARNING "makinput source not found: ${case_source}")
        return()
    endif()

    add_executable(${target_name}
        ${DATA_COMMON_SOURCES}
        ${case_source}
    )

    target_include_directories(${target_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/module
        ${CMAKE_SOURCE_DIR}/module/fluid/MPM
        ${CMAKE_SOURCE_DIR}/module/fluid/FEM
        ${CMAKE_SOURCE_DIR}/module/solid/implicit
    )

    target_link_libraries(${target_name} PRIVATE mpm_modules)
endfunction()

add_makinput_target(fluid)
add_makinput_target(solid)
add_makinput_target(fsi)
