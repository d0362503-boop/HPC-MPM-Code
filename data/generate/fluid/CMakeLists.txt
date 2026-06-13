add_executable(makinput_fluid
    makinput.cpp
    fluid_generator.cpp
    $<TARGET_OBJECTS:mpm_generate_common>
)

target_link_libraries(makinput_fluid PRIVATE mpm_modules)

set_target_properties(makinput_fluid PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(input.txt input.txt COPYONLY)
