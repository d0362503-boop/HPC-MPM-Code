mpm_add_data_tool(
    TARGET makinput_solid
    SOURCES makinput.cpp solid_generator.cpp
)

mpm_copy_input_once(input.txt)
