mpm_add_divide_tool(
    TARGET makdivide_solid
    SOURCES makdivide.cpp solid_divider.cpp
)

configure_file(para_input_data.txt para_input_data.txt COPYONLY)
