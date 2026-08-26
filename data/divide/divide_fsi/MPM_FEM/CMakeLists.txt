mpm_add_divide_tool(
    TARGET makdivide_fsi
    SOURCES makdivide.cpp fsi_divider.cpp
)

configure_file(para_input_data.txt para_input_data.txt COPYONLY)
