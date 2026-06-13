add_executable(makinput_fsi
    makinput.cpp
    fsi_generator.cpp
    $<TARGET_OBJECTS:mpm_generate_common>
)

target_link_libraries(makinput_fsi PRIVATE mpm_modules)

set_target_properties(makinput_fsi PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(input.txt input.txt COPYONLY)
