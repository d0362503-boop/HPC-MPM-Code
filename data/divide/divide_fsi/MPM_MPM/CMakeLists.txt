mpm_add_divide_tool(
    TARGET makdivide_fsi_mpm_mpm
    SOURCES
        makdivide.cpp
        fsi_divider.cpp
        ${CMAKE_SOURCE_DIR}/data/divide/divide_fluid/fluid_divider.cpp
        ${CMAKE_SOURCE_DIR}/data/divide/divide_solid/solid_divider.cpp
)

configure_file(para_input_data.txt para_input_data.txt COPYONLY)
