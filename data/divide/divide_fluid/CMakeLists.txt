add_executable(makdivide_fluid
    makdivide.cpp
    fluid_divider.cpp
    $<TARGET_OBJECTS:mpm_divide_common>
)

target_link_libraries(makdivide_fluid PRIVATE mpm_modules)

set_target_properties(makdivide_fluid PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(para_input_data.txt para_input_data.txt COPYONLY)
