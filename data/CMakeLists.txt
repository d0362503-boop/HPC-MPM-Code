# -----------------------------------------------------------------------------
# Data generators (makinput tools)
# -----------------------------------------------------------------------------

set(g_selected_data_case_count 0)
if(USE_DATA_FLUID)
    math(EXPR g_selected_data_case_count "${g_selected_data_case_count} + 1")
endif()
if(USE_DATA_SOLID)
    math(EXPR g_selected_data_case_count "${g_selected_data_case_count} + 1")
endif()
if(USE_DATA_FSI)
    math(EXPR g_selected_data_case_count "${g_selected_data_case_count} + 1")
endif()

if(NOT g_selected_data_case_count EQUAL 1)
    message(FATAL_ERROR
        "Exactly one data generator must be enabled. "
        "Set exactly one of USE_DATA_FLUID, USE_DATA_SOLID, USE_DATA_FSI to ON.")
endif()

add_library(mpm_data_common OBJECT
    data_generator.cpp
    output_util.cpp
)

target_link_libraries(mpm_data_common PUBLIC mpm_modules)

file(REMOVE
    ${CMAKE_BINARY_DIR}/data/libmpm_data_core.a
    ${CMAKE_BINARY_DIR}/data/makinput_fluid
    ${CMAKE_BINARY_DIR}/data/makinput_solid
    ${CMAKE_BINARY_DIR}/data/makinput_fsi
)

if(NOT USE_DATA_FLUID)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/fluid)
endif()
if(NOT USE_DATA_SOLID)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/solid)
endif()
if(NOT USE_DATA_FSI)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/fsi)
endif()

if(USE_DATA_FLUID)
    add_subdirectory(fluid)
endif()

if(USE_DATA_SOLID)
    add_subdirectory(solid)
endif()

if(USE_DATA_FSI)
    add_subdirectory(fsi)
endif()
