mpm_add_divide_tool(
    TARGET makdivide_fluid
    SOURCES makdivide.cpp fluid_divider.cpp
)

mpm_copy_input_once(para_input_data.txt)
