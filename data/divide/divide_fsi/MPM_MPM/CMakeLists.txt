mpm_add_divide_tool(
    TARGET makdivide_fsi_mpm_mpm
    SOURCES
        makdivide.cpp
        fsi_divider.cpp
        ${CMAKE_SOURCE_DIR}/data/divide/divide_fluid/fluid_divider.cpp
        ${CMAKE_SOURCE_DIR}/data/divide/divide_solid/solid_divider.cpp
)

mpm_copy_input_once(para_input_data.txt)
