add_executable(makinput_solid
    makinput.cpp
    solid_generator.cpp
    $<TARGET_OBJECTS:mpm_data_common>
)

target_link_libraries(makinput_solid PRIVATE mpm_modules)

set_target_properties(makinput_solid PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(input.txt input.txt COPYONLY)
