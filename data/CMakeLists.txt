# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

add_library(mpm_data_core
    data_generator.cpp
    output_util.cpp
    fluid/fluid_generator.cpp
    solid/solid_generator.cpp
    fsi/fsi_generator.cpp
)

target_include_directories(mpm_data_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/module
    ${CMAKE_SOURCE_DIR}/module/fluid/FEM
    ${CMAKE_SOURCE_DIR}/module/solid/implicit
)

target_link_libraries(mpm_data_core PUBLIC mpm_modules)

function(configure_makinput_target target_name case_dir)
    add_executable(${target_name} ${case_dir}/makinput.cpp)

    target_link_libraries(${target_name} PRIVATE mpm_data_core)

    set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data
    )

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${case_dir}/input.txt
         DESTINATION ${CMAKE_BINARY_DIR}/data/${case_dir})
endfunction()

if(USE_DATA_FLUID)
    configure_makinput_target(makinput_fluid fluid)
endif()

if(USE_DATA_SOLID)
    configure_makinput_target(makinput_solid solid)
endif()

if(USE_DATA_FSI)
    configure_makinput_target(makinput_fsi fsi)
endif()
