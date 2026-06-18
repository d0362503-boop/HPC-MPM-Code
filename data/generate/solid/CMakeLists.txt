mpm_add_data_tool(
    TARGET makinput_solid
    SOURCES makinput.cpp solid_generator.cpp
)

configure_file(input.txt input.txt COPYONLY)
