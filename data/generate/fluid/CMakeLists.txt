mpm_add_data_tool(
    TARGET makinput_fluid
    SOURCES makinput.cpp fluid_generator.cpp
)

configure_file(input.txt input.txt COPYONLY)
