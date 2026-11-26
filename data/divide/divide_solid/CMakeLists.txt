mpm_add_divide_tool(
    TARGET makdivide_solid
    SOURCES makdivide.cpp solid_divider.cpp
)

mpm_copy_input_once(para_input_data.txt)
