mpm_add_divide_tool(
    TARGET makdivide_fluid
    SOURCES makdivide.cpp fluid_divider.cpp
)

configure_file(para_input_data.txt para_input_data.txt COPYONLY)
