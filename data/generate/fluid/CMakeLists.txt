mpm_add_data_tool(
    TARGET makinput_fluid
    SOURCES makinput.cpp fluid_generator.cpp
)

mpm_copy_input_once(input.txt)
