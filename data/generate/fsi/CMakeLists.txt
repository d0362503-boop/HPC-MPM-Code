mpm_add_data_tool(
    TARGET makinput_fsi
    SOURCES makinput.cpp fsi_generator.cpp
)

configure_file(input.txt input.txt COPYONLY)
