mpm_add_data_tool(
    TARGET makinput_fsi_mpm_mpm
    SOURCES
        makinput.cpp
        fsi_generator.cpp
        ${CMAKE_SOURCE_DIR}/data/generate/fluid/fluid_generator.cpp
        ${CMAKE_SOURCE_DIR}/data/generate/solid/solid_generator.cpp
)

configure_file(input.txt input.txt COPYONLY)
